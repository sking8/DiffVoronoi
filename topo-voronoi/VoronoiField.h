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

using namespace Meso;
template<int d> class VoronoiField : public Simulator
{
	Typedef_VectorD(d); Typedef_MatrixD(d);
public:
	Grid<d> grid;
	VoronoiParticles<d> particles;	////main particles

	////field attributes
	Field<real, d> rho;						////Voronoi density field
	Field<int, d> active;					////active flag for cells: if =0, inactive; if =1, active; if =2, (outer) boundary; if =3, fixed. Use (int)VoronoiCellType to specify the values. Initialized in driver
	Field<real, d> soft_max_sum;			////precomputed denominator in softmax

	////parameters
	real epsi_A = (real)0;					////avoid singular A, not used by default
	real epsi_S = (real)0;					////for free boundary, not used by default
	int beta = 50;							////power index

	////neighbor searching
	int nb_n;
	std::shared_ptr<NeighborSearcher<d> > nbs_searcher;         ////radius search
	Array<Array<int> > nbs;										////nb particles index for each cell
	Array<Array<int> > nbs_p;									////nb cells for each particle, Fan: may not be neccesary here

	////derivatives
	Array<Array<VectorD> > drho_dx;								////first index for cell; second index for particle
	Array<Array<MatrixD> > drho_dD;								////first index for cell; second index for particle

	//////////////////////////////////////////////////////////////////////////
	////Initialization
	void Initialize(const Grid<d> _grid, const Array<VectorD>& points, const int nb_n);

	//////////////////////////////////////////////////////////////////////////
	////Field update functions
	virtual void Output(DriverMetaData& metadata) {
		std::string vts_name = fmt::format("vts{:04d}.vts", metadata.current_frame);
		bf::path vtk_path = metadata.base_path / bf::path(vts_name);
		VTKFunc::Write_VTS(rho, vtk_path.string());
	}
	virtual real CFL_Time(const real cfl) { return 1; }
	virtual void Advance(DriverMetaData& metadata);		////updates state variables including nbs, softmax, and rho. It does not update any sensitivities.
	void Update_Neighbors();							////update point nbs 
	void Update_A();									////update A according to D
	void Update_Softmax_Sum();							////update s as a precomputed field
	void Update_Rho();									////update rho on the field
	void Update_DRho_DX();
	void Update_DRho_DD();

	void Numerical_Derivative_DRho_DX();
	void Numerical_Derivative_DRho_DD();

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
		return exp(-Dist(p_idx, grid.Cell_Center(cell))) / soft_max_sum(cell);
	}

	inline real Softmax_Epsi(const VectorDi& cell)	////softmax for epsilon, free boundary case
	{
		return epsi_S / soft_max_sum(cell);
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
		VectorD pos = grid.Cell_Center(cell); return -Softmax(m, cell) * (Delta(m, n) - Softmax(n, cell)) * Y_helper(n, pos);
	}

	inline VectorD dS_Epsi_dX(const int n, const VectorDi& cell) //free boundary case
	{
		VectorD pos = grid.Cell_Center(cell); return Softmax_Epsi(cell) * Softmax(n, cell) * Y_helper(n, pos);
	}

	inline MatrixD dS_dD(const int m, const int n, const VectorDi& cell)
	{
		VectorD pos = grid.Cell_Center(cell);
		VectorD X = X_helper(n, pos);
		return -Softmax(m, cell) * (Delta(m, n) - Softmax(n, cell)) * Dist(n, pos) * X * X.transpose() * particles.D(n);
	}

	inline MatrixD dS_Epsi_dD(const int n, const VectorDi& cell) //free boundary case
	{
		VectorD pos = grid.Cell_Center(cell);
		VectorD X = X_helper(n, pos);
		return Softmax_Epsi(cell) * Softmax(n, cell) * Dist(n, pos) * X * X.transpose() * particles.D(n);
	}
};