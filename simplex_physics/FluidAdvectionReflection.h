//////////////////////////////////////////////////////////////////////////
// advection-reflection fluid
// Copyright (c) (2021-), Fan Feng
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include "FluidEulerOld.h"
#ifndef __FluidAdvectionReflection_h__
#define __FluidAdvectionReflection_h__
template<int d> class FluidAdvectionReflection : public FluidEulerOld<d>
{
	Typedef_VectorDii(d); using Base = FluidEulerOld<d>;
public:
	using Base::velocity;
	using Base::mac_grid;
	using Base::projection;
	using Base::use_maccormack;

	FaceField<real, d> velocity_temp;

	FluidAdvectionReflection(const SolverType& _solver_mode = SolverType::MULTIGRID_AUTO) :Base(_solver_mode) {

	}

	virtual void Advance(const real dt)
	{
		Base::Advection(dt/(real)2);
		Base::Apply_Body_Forces(dt / (real)2);
		velocity_temp = velocity; //u_tilde
		Base::Enforce_Incompressibility();

		//u_hat^(1/2)=2u^(1/2)-u_tilde(1/2)
		velocity *= 2; 
		velocity -= velocity_temp;

		Base::Advection(dt/(real)2);
		Base::Enforce_Incompressibility();
		//Divergence();
	}
};
#endif
