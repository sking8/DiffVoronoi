#include "VoronoiField.h"
#include "TopoCellType.h"

namespace Meso {
	//////////////////////////////////////////////////////////////////////////
	////initialization
	template<int d> void VoronoiField<d>::Initialize(const Grid<d> _grid, const Array<VectorD>& points, const int _nb_n, int _beta, real _alpha, real _c)
	{
		int p_num = (int)points.size();
		grid = _grid;
		nb_n = _nb_n;
		particles.Resize(p_num);
		particles.xRef() = points;
		beta = _beta;
		ArrayFunc::Fill(particles.DRef(), _alpha * MatrixD::Identity());
		c = _c;

		////grid attributes
		rho.Init(grid);
		rho.Fill((real)0);
		soft_max_sum.Init(grid);
		soft_max_sum.Fill((real)0);
		active.Init(grid); //Fan: need to consider this, initialize outside of the class
		active.Fill((int)TopoCellType::Active);

		nbs_c.Init(grid);
		drho_dx.Init(grid);
	}

	//////////////////////////////////////////////////////////////////////////
	////field updates

	template<int d> void VoronoiField<d>::Advance(DriverMetaData& metadata)
	{
		Update_A();
		Update_Neighbors();
		Update_Softmax_Sum();
		Update_Rho();
		//Numerical_Derivative_DRho_DX();
	}

	template<int d> void VoronoiField<d>::Update_A()
	{
		particles.Exec_Points(
			[&](const int idx) {
				particles.A(idx) = particles.D(idx) * particles.D(idx).transpose() + epsi_A * MatrixD::Identity();
				/*Info("p: {}", idx);
				std::cout << "A: " << particles.A(idx) << std::endl;
				std::cout << "D: " << particles.D(idx) << std::endl;*/
			}
		);
	}

	template<int d> void VoronoiField<d>::Update_Neighbors()
	{
		////nb particles of each cell
		particles.Update_Searcher();

		grid.Exec_Nodes(
			[&](const VectorDi cell) {
				nbs_c(cell).clear();
				VectorD pos = grid.Position(cell);
				if (active(cell) != (int)TopoCellType::Active) { return; }

				particles.nbs_searcher.Find_K_Nearest_Nbs(pos, nb_n, nbs_c(cell));
				int nb_n = nbs_c(cell).size();
				if (nb_n > 0) {
					drho_dx(cell).resize(nb_n);
					//drho_dD[idx].resize(nb_c_n);
				}
			}
		);
	}

	template<int d> void VoronoiField<d>::Update_Softmax_Sum()
	{
		grid.Exec_Nodes(
			[&](const VectorDi cell) {
				VectorD pos = grid.Position(cell);
				if (active(cell) != (int)TopoCellType::Active) {
					soft_max_sum(cell) = (real)0.; return;
				}

				real sm_sum = (real)0;
				int nb_c_n = nbs_c(cell).size();
				for (int j = 0; j < nb_c_n; j++) {
					int pid = nbs_c(cell)[j];
					real dis = Dist(pid, pos);
					sm_sum += exp(-dis);
				}
				soft_max_sum(cell) = sm_sum;
			}
		);
	}

	template<int d> void VoronoiField<d>::Update_Rho()
	{
		grid.Exec_Nodes(
			[&](const VectorDi cell) {
				VectorD pos = grid.Position(cell);
				if (active(cell) != (int)TopoCellType::Active) { return; }

				int nb_n = nbs_c(cell).size();
				real s = (real)0;
				for (int j = 0; j < nb_n; j++) {
					int nb_p = nbs_c(cell)[j];
					s += pow((real)1 - Softmax(nb_p, cell), (real)beta);
				}
				rho(cell) = (real)(nb_n - c) - s;
			}
		);
	}

	//////////////////////////////////////////////////////////////////////////
	////sensitivities

	template<int d> void VoronoiField<d>::Update_DRho_DX()
	{
		grid.Exec_Nodes(
			[&](const VectorDi cell) {
				VectorD pos = grid.Position(cell);
				if (active(cell) != (int)TopoCellType::Active) { return; }

				int nb_c_n = nbs_c(cell).size();
				real s = (real)0;
				for (int j = 0; j < nb_c_n; j++) {
					int nb_p = nbs_c(cell)[j];
					VectorD drhodx = VectorD::Zero();
					for (int k = 0; k < nb_c_n; k++) {
						int m = nbs_c(cell)[k];
						VectorD dsdx = dS_dX(m, nb_p, cell);
						real S = Softmax(m, cell);
						drhodx += beta * pow(1 - S, beta - 1) * dsdx;
					}
					drho_dx(cell)[j] = drhodx;
				}
			}
		);

		Info("Derivative drho_dx is updated");
	}

	//////////////////////////////////////////////////////////////////////////
	////numerical derivatives

	template<int d> void VoronoiField<d>::Numerical_Derivative_DRho_DX()
	{
		real delta_x = (real)1e-6;
		int p_size = particles.Size();

		std::cout << "Numerical DRho_DX" << std::endl;
		Update_DRho_DX();
		Field<real, d> rho_test = rho;
		Field<Array<VectorD>, d> numeric_derv(grid);

		//should not parallelize here
		grid.Iterate_Nodes(
			[&](const VectorDi cell) {
				int nb_c_n = nbs_c(cell).size();
				numeric_derv(cell).resize(nb_c_n);
				for (int j = 0; j < nb_c_n; j++) {
					int nb_p = nbs_c(cell)[j];
					VectorD old_pos = particles.x(nb_p);
					for (int k = 0; k < d; k++) { //iterate through dimensions
						particles.x(nb_p)[k] += delta_x;
						Update_Softmax_Sum();
						Update_Rho();
						numeric_derv(cell)[j][k] = (rho(cell) - rho_test(cell)) / delta_x;
						particles.x(nb_p)[k] = old_pos[k];
					}
					if (!Meso::MathFunc::All_Close(numeric_derv(cell)[j], drho_dx(cell)[j], (real)1e-2, (real)1e-3)) {
						Warn("cell:{}, nb:{}, analytical:{}, numerical:{}", cell, j, drho_dx(cell)[j], numeric_derv(cell)[j]);
					}
				}
			}
		);

		Pass("Finished test of numerical derivative");
	}

	template class VoronoiField<2>;
	template class VoronoiField<3>;
}