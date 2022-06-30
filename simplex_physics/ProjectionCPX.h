//////////////////////////////////////////////////////////////////////////
// Project a vector field to divergence free on a MAC grid
// Copyright (c) (2018-), Bo Zhu, Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
//// This projection solver currently only takes zero Neumann bc. 
// The velocities on the Neumann boundary need to be set to the correct values before the projection.
//// The calculated pressure value is the real pressure value divided by delta_x.
//////////////////////////////////////////////////////////////////////////
#pragma once
#ifdef USE_CPX
#include <functional>
#include "MacGrid.h"
#include "FaceField.h"
#include "Field.h"
#include "TypeFunc.h"
#include "BoundaryCondition.h"
#include "Poisson.h"
#include "Json.h"

// example enum type declaration
enum SurfaceTensionMode {
	NONE = 0,//no surface tension
	EXPLICIT,//explicit surface tension, implemented by a jump condition at the fluid-air interface
	IMPLICIT,//implicit surface tension
};

// map TaskState values to JSON as strings
NLOHMANN_JSON_SERIALIZE_ENUM(SurfaceTensionMode, {
	{NONE, "ST_NONE"},
	{EXPLICIT, "ST_EXPLICIT"},
	{IMPLICIT, "ST_IMPLICIT"},
	})

template<class T, int d> class ProjectionCPX
{
	Typedef_VectorDii(d);
public:
	//data
	const MacGrid<d>* mac_grid = nullptr;
	const Field<ushort, d>* type = nullptr;
	const BoundaryConditionMacGrid<d>* bc = nullptr;

	//settings
	bool verbose=true;
	bool update_A=true;

	//flags
	bool A_initialized = false, cpx_initialized = false;
	
	////cpx
	Poisson<d> cpx_poisson;
	Field<T, d> cell_b;

	int max_iter;
	real relative_tolerance;

public:
	////constructors
	ProjectionCPX(json &j, const MacGrid<d>* _mac_grid, const Field<ushort, d>* _type = nullptr, const BoundaryConditionMacGrid<d>* _bc = nullptr) {
		Info("ProjectionCPX initialize from json: {}", j.dump());
		verbose = Json::Value<bool>(j, "verbose", true);
		update_A = Json::Value<bool>(j, "update_A", true);
		max_iter = Json::Value<int>(j, "max_iter", -1);
		relative_tolerance = Json::Value<real>(j, "relative_tolerance", 1e-4);

		mac_grid = _mac_grid;
		type = _type;
		bc = _bc;
	}

	~ProjectionCPX() {
	}

	////projection functions
	virtual void Build_A(const T dt = 0);
	virtual void Build_b(const T dt, const FaceField<T, d>& velocity);				////calculate b as div velocity
	virtual void Build(const T dt, const FaceField<T, d>& velocity);					////call allocate, update_A, and update_b
	virtual void Correction(const T dt, FaceField<T, d>& velocity);
	virtual void Solve(const T dt = 0);
	virtual void Project(const T dt, FaceField<T, d>& velocity);					////call both build, solve, and correction
	
	////Physical interface functions that defines the problem
	//NOTE: if you want to design a derived class of this, theoretically you only need to implement these 5 functions
	virtual T Off_Diag_Term(const T dt, const VectorDi& fluid_cell, const int& nbidx)const;
	virtual T Diag_Face_Term(const T dt, const int& axis, const VectorDi& face)const;
	virtual T Velocity_Correction(const T dt, const int& axis, const VectorDi& face)const;
	virtual bool Is_Valid_Cell(const VectorDi& cell) const {
		return mac_grid->grid.Valid_Cell(cell);
	}
	virtual bool Is_Fluid_Cell(const VectorDi& cell) const {
		return Is_Valid_Cell(cell) && (!type || (*type)(cell) == (ushort)CellType::Fluid);
	}
};
#endif
