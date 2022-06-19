//////////////////////////////////////////////////////////////////////////
// Linear Grid FEM
// Copyright (c) (2022-), Bo Zhu, Fan Feng
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "Hashtable.h"
#include "LinearFEMFunc.h"
#include "BoundaryConditionGrid.h"
#include "Field.h"
#include "Simulator.h"
#include "IOHelper.h"
//#include "GmgPcgSolverCPU.h"
//#ifdef USE_CUDA
//#include "GmgPcgSolverGPU.h"
//#endif
using namespace Meso;
template<int d> class LinearFEMGrid: public Simulator
{
	Typedef_VectorD(d); Typedef_MatrixD(d);
public:
	Grid<d> grid;					//this grid is a corner grid, represending the corner nodes
	SparseMatrix<real> K;			//global stiffness matrix
	VectorX u;						//displacement on corner
	VectorX f;						//force on corner

	Array<MatrixX> Ke;				//stiffness matrix of different materials
	Field<short, d> material_id;	//assuming limited number of materials existing
	BoundaryConditionGrid<d> bc;

	////multigrid and parallelization
	/*bool use_multigrid_solver = true;
	MultiGrid::Params multigrid_params;
	GMGPCG_Solver_CPU<d> gmg_solver_cpu;
#ifdef USE_CUDA
	GMGPCG_Solver_GPU<real, d> gmg_solver_gpu;
#endif*/

	//Array<int> colored_cell_ptr;	//mark where each color starts and ends
	//Array<int> colored_cell_indices;	//indices of colors in each cell

	virtual void Output(DriverMetaData& metadata); //Fan: to visualize the displacement field and strain field
	virtual void Advance(DriverMetaData& metadata) { Update_K_And_f(); Solve(); }
	virtual real CFL_Time(const real dt) { return 0; }

	virtual void Initialize(const Grid<d> _grid, const BoundaryConditionGrid<d>& _bc, const Array<std::tuple<real, real>>& _materials, const Field<short, d>& _material_id);
	void Allocate_K();
	virtual void Update_K_And_f();
	virtual void Solve();

	void Add_Material(const std::tuple<real, real> material);
	void Clear_Materials() { Ke.clear(); }

	void Compute_Cell_Displacement(const VectorX& u, const VectorDi& cell, VectorX& cell_u) const;

	virtual void Compute_Elastic_Energy(Field<real, d>& energy) const;

	//Need numerical verification for derivatives
};