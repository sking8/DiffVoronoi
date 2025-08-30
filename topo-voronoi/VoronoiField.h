#pragma once
#include "Grid.h"
#include "Field.h"
#include "Points.h"
#include "NeighborSearcher.h"
#include "Random.h"
#include "Interpolation.h"
#include "Timer.h"
#include "VoronoiParticles.h"
#include "TopoCellType.h"
#include "MetaData.h"
#include "Simulator.h"
#include "IOHelper.h"

namespace Meso{
	template<int d> class VoronoiField : public Simulator
	{
		Typedef_VectorD(d); Typedef_MatrixD(d);
	public:
		Grid<d> grid;
		VoronoiParticles<d, NeighborKDTree<d>> particles;			////main particles

		////field attributes
		Field<real, d> rho;						////Voronoi density field
		Field<int, d> active;					////active flag for cells: if =0, inactive; if =1, active; if =2, (outer) boundary; if =3, fixed. Use (int)VoronoiCellType to specify the values. Initialized in driver
		Field<real, d> soft_max_sum;			////precomputed denominator in softmax

		////parameters
		real epsi_A = (real)0;					////avoid singular A, not used by default
		int beta;								////power index
		real c;

		////neighbor searching
		int nb_n;
		Field<Array<int>,d> nbs_c;									////nb particles index for each cell
		Array<Array<int> > nbs_p;									////nb cells for each particle, Fan: may not be neccesary here

		////derivatives
		Field<Array<VectorD>,d> drho_dx;							////first index for cell; second index for particle

		//////////////////////////////////////////////////////////////////////////
		////Initialization
		void Initialize(const Grid<d> _grid, const Array<VectorD>& points, const int nb_n, int beta, real _alpha, real _c);

		//////////////////////////////////////////////////////////////////////////
		////Field update functions
		virtual void Output(DriverMetaData& metadata) {
			std::string vts_name = fmt::format("vts{:04d}.vts", metadata.current_frame);
			fs::path vtk_path = metadata.base_path / fs::path(vts_name);
			VTKFunc::Write_VTS(rho, vtk_path.string());

			std::string vtu_name = fmt::format("points{:04d}.vtu", metadata.current_frame);
			fs::path vtu_path = metadata.base_path / fs::path(vtu_name);
			VTKFunc::Write_Points<d>(particles.xRef(), vtu_path.string());
		}

		virtual real CFL_Time(const real cfl) { return 1; }
		virtual void Advance(DriverMetaData& metadata);		////updates state variables including nbs_c, softmax, and rho. It does not update any sensitivities.
		void Update_Neighbors();							////update point nbs_c 
		void Update_A();									////update A according to D
		void Update_Softmax_Sum();							////update s as a precomputed field
		void Update_Rho();									////update rho on the field
		void Update_DRho_DX();

		void Numerical_Derivative_DRho_DX();

	protected:
		inline real Delta(const int p, const int q)
		{
			return p == q ? 1 : 0;
		}

		inline real Dist(const int p_idx, const VectorD& pos)
		{
			VectorD vec = particles.x(p_idx) - pos; return sqrt(vec.transpose() * particles.A(p_idx) * vec);
		}

		inline real Softmax(const int p_idx, const VectorDi& cell)	////use precomputed softmax
		{
			return exp(-Dist(p_idx, grid.Position(cell))) / soft_max_sum(cell);
		}

		inline VectorD X_helper(const int p_idx, const VectorD& pos)
		{
			return (particles.x(p_idx) - pos) / (Dist(p_idx, pos));
		}

		inline VectorD Y_helper(const int p_idx, const VectorD& pos)
		{
			return particles.A(p_idx) * X_helper(p_idx, pos);
		}

		inline VectorD dS_dX(const int m, const int n, const VectorDi& cell)
		{
			VectorD pos = grid.Position(cell); return -Softmax(m, cell) * (Delta(m, n) - Softmax(n, cell)) * Y_helper(n, pos);
		}

		inline MatrixD dS_dD(const int m, const int n, const VectorDi& cell)
		{
			VectorD pos = grid.Position(cell);
			VectorD X = X_helper(n, pos);
			return -Softmax(m, cell) * (Delta(m, n) - Softmax(n, cell)) * Dist(n, pos) * X * X.transpose() * particles.D(n);
		}
	};
}