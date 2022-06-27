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
	Field<VectorD, d> u;			//displacement on corner
	Field<VectorD, d> f;			//force on corner

	Array<int> compact_indices;		//convert index in the padded indices to compact indices
	Array<MatrixX> Ke;				//stiffness matrix of different materials
	Field<short, d> material_id;	//assuming limited number of materials existing
	BoundaryConditionGrid<d> bc;

	virtual void Output(DriverMetaData& metadata);
	virtual void Advance(DriverMetaData& metadata) { Update_K(); Solve(); }
	virtual real CFL_Time(const real dt) { return 1; }

	virtual void Initialize(const Grid<d> _grid, const BoundaryConditionGrid<d>& _bc, const Array<std::tuple<real, real>>& _materials, const Field<short, d>& _material_id);
	void Allocate_K();
	virtual void Update_K();
	virtual void Solve();

	void Add_Material(const std::tuple<real, real> material);
	void Clear_Materials() { Ke.clear(); }

	void Compute_Cell_Displacement(const VectorDi& cell, VectorX& cell_u) const;

	virtual void Compute_Elastic_Energy(Field<real, d>& energy) const;

	int Compact_Idx(int idx);
	//Need numerical verification for derivatives
};