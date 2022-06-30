//////////////////////////////////////////////////////////////////////////
// Incompressible fluid solver on an Eulerian grid, supporting two phase flow and surface tension
// Copyright (c) (2018-), Bo Zhu, Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "ProjectionTwoPhaseCPX.h"
#include "Timer.h"
#include "FluidEulerFreeSurface.h"

template<int d> class FluidEulerTwoPhase : public FluidEulerFreeSurface<d>
{Typedef_VectorDii(d);using Base=FluidEulerFreeSurface<d>;
public:
	using Base::mac_grid;
	using Base::velocity;
	using Base::type;			////supporting type: Fluid, NB, Air
	using Base::bc;
	using Base::g;
	using Base::use_body_force;	////default false
	using Base::projection;
	using Base::Enforce_Boundary_Conditions;
	using Base::use_implicit_surface_tension;
	using Base::levelset;
	using Base::narrow_band_width;
	using Base::narrow_band_cell_num;
	using Base::dirac_band_width;
	using Base::dirac_band_cell_num;
	using Base::verbose;

	real rhoA;
	real rhoL;
	FaceField<real,d> rho;

public:
	FluidEulerTwoPhase(json& j) :
		Base(j)
	{
		projection = std::make_shared<ProjectionTwoPhaseCPX<real, d>>(
			j.at("projection"),
			&mac_grid, &levelset, &rho, &type, &bc
			);
		rhoA = j.value<real>("rhoA", 0.001);
		rhoL = j.value<real>("rhoL", 1.);
	}

	virtual void Initialize(const VectorDi& cell_counts,const real dx,const VectorD& domain_min=VectorD::Zero())
	{
		Base::Initialize(cell_counts,dx,domain_min);
		rho.Resize(mac_grid.grid.cell_counts,(real)0);
	}

	virtual void Advance(const real dt)
	{
		Timer timer;						if(verbose)timer.Reset();
		Advection(dt);						if(verbose)timer.Elapse_And_Output_And_Reset("\t\ttimestep: advection");
		Apply_Body_Forces(dt);				if(verbose)timer.Elapse_And_Output_And_Reset("\t\ttimestep: body force");
		auto proj_ptr = std::dynamic_pointer_cast<ProjectionTwoPhaseCPX<real, d>>(projection);
		if (use_implicit_surface_tension) implicit_surface_tension_solver->Solve(dt, velocity, proj_ptr->sigma, narrow_band_width, dirac_band_width);
		Enforce_Incompressibility(dt);		if(verbose)timer.Elapse_And_Output_And_Reset("\t\ttimestep: projection");
	}



	virtual void Advection(const real dt)
	{
		Timer timer;if(verbose)timer.Reset();

		////advect velocity
		MacGrid<d> ghost_grid=mac_grid;FaceField<real,d> ghost_velocity=velocity;
		Advection::Semi_Lagrangian(dt,ghost_velocity,ghost_grid,ghost_velocity,mac_grid,velocity,use_zero_extrapolation);
		if(verbose)timer.Elapse_And_Output_And_Reset("\tadv: vel adv");

		////advect interface
		Field<real,d> ghost_phi=levelset.phi;
		Advection::Semi_Lagrangian(dt,ghost_velocity,ghost_grid,ghost_phi,mac_grid,levelset.phi,false);
		if(verbose)timer.Elapse_And_Output_And_Reset("\tadv: ls adv");

		levelset.Fast_Marching(narrow_band_width);
		if(verbose)timer.Elapse_And_Output_And_Reset("\tadv: ls fmm");

		Update_Cell_Types();
		if(verbose)timer.Elapse_And_Output_And_Reset("\tadv: celltype");

		Update_Rho();
		if(verbose)timer.Elapse_And_Output_And_Reset("\tadv: rho");
	}

	virtual void Update_Cell_Types()
	{
		int cell_num=mac_grid.grid.Number_Of_Cells();
		#pragma omp parallel for
		for(int i=0;i<cell_num;i++){
			const VectorDi& cell=mac_grid.grid.Cell_Coord(i);
			type(cell)=(ushort)CellType::Fluid;
			if(bc.Is_Psi_D(cell))type(cell)=bc.Psi_D_Type(cell);
		}
	}

	virtual void Update_Rho()
	{
		for(int axis=0;axis<d;axis++){int face_num=mac_grid.Number_Of_Faces(axis);
			#pragma omp parallel for
			for(int j=0;j<face_num;j++){
				VectorDi face=mac_grid.Face_Coord(axis,j);
                VectorDi cell_lr[2]; 
                real phi_lr[2];
                real tp[2];
                for(int i=0;i<2;i++){
                	cell_lr[i] = MacGrid<d>::Face_Incident_Cell(axis, face, i);
                	if(!mac_grid.grid.Valid_Cell(cell_lr[i])){phi_lr[i]=(real)0.;}
                	else{phi_lr[i]=levelset.phi(cell_lr[i]);}
                	tp[i]=std::max(phi_lr[i],(real)0.);}
				real tp_theta=(tp[0]+tp[1])/(abs(phi_lr[0])+abs(phi_lr[1]));
				rho.face_fields[axis](face)=rhoA*tp_theta+rhoL*((real)1.-tp_theta);}}
	}

	virtual void Apply_Body_Forces(const real dt)
	{
		if(!use_body_force)return;
		for(int axis=0;axis<d;axis++){int face_num=mac_grid.Number_Of_Faces(axis);
			#pragma omp parallel for
			for(int i=0;i<face_num;i++){VectorDi face=mac_grid.Face_Coord(axis,i);
				velocity.face_fields[axis](face)+=g[axis]*dt;}}
	}
};