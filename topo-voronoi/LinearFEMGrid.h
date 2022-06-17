//////////////////////////////////////////////////////////////////////////
// Linear Grid FEM
// Copyright (c) (2022-), Bo Zhu, Fan Feng
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "Hashtable.h"
#include "LinearFEMFunc.h"
#include "BoundaryConditionGrid.h"
#include "Field.h"
//#include "GmgPcgSolverCPU.h"
//#ifdef USE_CUDA
//#include "GmgPcgSolverGPU.h"
//#endif
using namespace Meso;
template<int d> class LinearFemGrid
{
	Typedef_VectorD(d); Typedef_MatrixD(d);
public:
	Grid<d> grid;					//this grid is a corner grid, represending the corner nodes
	SparseMatrix<real> K;
	VectorX u;
	VectorX f;

	Array<MatrixX> Ke;				//stiffness matrix of different materials
	Field<short, d> material_id;	//assumes multiple materials existing
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

	virtual void Initialize(const Grid<d> _grid);
	void Allocate_K();
	virtual void Update_K_And_f();
	virtual void Solve();

	void Add_Material(real youngs, real poisson);
	void Clear_Materials() { Ke.clear(); }

	virtual void Set_Fixed(const VectorDi& node);
	void Set_Displacement(const VectorDi& node, const VectorD& dis);
	virtual void Set_Force(const VectorDi& node, const VectorD& force);

	void Compute_Cell_Displacement(const VectorX& u, const VectorDi& cell, VectorX& cell_u) const;

	virtual void Compute_Elastic_Energy(Field<real, d>& energy) const;
};