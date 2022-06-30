//////////////////////////////////////////////////////////////////////////
// Project a vector field with free boundary to divergence free on a MAC grid
// Copyright (c) (2018-), Fan Feng, Shuqi Yang, Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////

#ifndef __ProjectionTwoPhaseCPX_h__
#define __ProjectionTwoPhaseCPX_h__

#ifdef USE_CPX

#include "ProjectionIrregularCPX.h"

template<class T, int d> class ProjectionTwoPhaseCPX : public ProjectionIrregularCPX<T, d>
{Typedef_VectorDii(d);
public:
	using Base=ProjectionIrregularCPX<T, d>;
	using Base::mac_grid;
	using Base::bc;
	using Base::verbose;
	using Base::cell_b;
	using Base::cpx_poisson;
	using Base::levelset;
	using Base::surface_tension_mode;
	using Base::sigma;
	using Base::apply_jump_condition;

	////divergence control
	using Base::use_vol_control;
	using Base::calc_current_vol;						////calculate the current vol within Apply_Vol_Control_To_b or not
	using Base::target_vol;						////ALWAYS NEEDS TO BE SET EXTERNALLY!
	using Base::current_vol;						////needs to be set externally or calculated within Apply_Vol_Control_To_b when calc_current_vol=true			
	using Base::vol_control_ks;

	//default narrow band width
	using Base::narrow_band_cell_num;

	//implicit gravity
	using Base::use_jump_gravity;
	using Base::gravity_coefficient;

	//face rho
	const FaceField<T, d>* rho;
	T rhoA;
	T rhoL;

public:
	////constructors
	ProjectionTwoPhaseCPX(json& j, const MacGrid<d>* _mac_grid, const LevelSet<d>* _levelset, const FaceField<T, d>* _rho, const Field<ushort, d>* _type = nullptr, const BoundaryConditionMacGrid<d>* _bc = nullptr)
		:Base(j, _mac_grid, _levelset, _type, _bc)
	{
		rhoA = j.value<T>("rhoA", 1e-3);
		rhoL = j.value<T>("rhoL", 1.);
		rho = _rho;
	}
	~ProjectionTwoPhaseCPX() {

	}

	////projection functions
	virtual void Apply_Jump_Condition_To_b(const T dt);
	virtual void Apply_Vol_Control_To_b();

	//////////////////////////////////////////////////////////////////////////
	////default callback functions: all levelset-related

	inline virtual T Pressure_Jump(const T dt, const Vector<T, d>& pos) const {
		T p_jump = 0;
		if (surface_tension_mode == SurfaceTensionMode::EXPLICIT) {
			T curvature = (*levelset).Curvature(pos);
			p_jump += sigma * curvature;
		}
		if (use_jump_gravity) {
			p_jump -= (rhoL - rhoA) * gravity_coefficient.dot(pos);
		}
		return p_jump * dt;
	}

	////Physical interface functions that defines the problem
	virtual T Off_Diag_Term(const T dt, const VectorDi& fluid_cell, const int& nbidx) const;
	virtual T Diag_Face_Term(const T dt, const int& axis, const VectorDi& face) const;
	virtual T Velocity_Correction(const T dt, const int& axis, const VectorDi& face) const;
	//Is_Valid_Cell: same as Projection<d>
	virtual bool Is_Fluid_Cell(const VectorDi& cell) const {return mac_grid->grid.Valid_Cell(cell);}
	virtual bool Is_Liquid_Cell(const VectorDi& cell) const {return Is_Fluid_Cell(cell)&&((*levelset).phi(cell) <= (T)0);}
	virtual bool Is_Air_Cell(const VectorDi& cell) const {return Is_Fluid_Cell(cell)&&((*levelset).phi(cell) > (T)0);}
};

#endif
#endif
