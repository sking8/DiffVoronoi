//////////////////////////////////////////////////////////////////////////
// Incompressible fluid solver on an Eulerian grid, supporting free interface and surface tension
// Copyright (c) (2018-), Bo Zhu, Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "LevelSet.h"
#include "ProjectionIrregular.h"
#include "ProjectionIrregularCPX.h"
#include "Timer.h"
#include "FluidEuler.h"
#include "ImplicitSurfaceTension.h"
#include "SpraySystem.h"

template<int d> class FluidEulerFreeSurface : public FluidEuler<d>
{Typedef_VectorDii(d);using Base=FluidEuler<d>;
public:
	using Base::mac_grid;
	using Base::velocity;
	using Base::type;			////supporting type: Fluid, NB, Air
	using Base::bc;
	using Base::g;
	using Base::use_body_force;	////default false
	using Base::projection;
	using Base::Enforce_Boundary_Conditions;

	std::shared_ptr<ImplicitSurfaceTension<real, d>> implicit_surface_tension_solver;
	LevelSet<d> levelset;
	real narrow_band_width = 0;
	real dirac_band_width = 0;

	std::shared_ptr<SpraySystem<d>> spray_system = nullptr;
protected:
	bool use_implicit_surface_tension;
	int narrow_band_cell_num = 3;
	int dirac_band_cell_num=3;
	bool verbose=false;
public:
	FluidEulerFreeSurface(json& j) :
		Base(j)
	{
		projection = std::make_shared<ProjectionIrregularCPX<real, d>>(
			j.at("projection"),
			&mac_grid, &levelset, &type, &bc
			);
		auto proj_ptr = std::dynamic_pointer_cast<ProjectionIrregularCPX<real, d>>(projection);
		use_implicit_surface_tension = (proj_ptr->surface_tension_mode == SurfaceTensionMode::IMPLICIT);
		verbose = j.value<bool>("verbose", false);
		narrow_band_cell_num = j.value<int>("narrow_band_cell_num", 3);
		dirac_band_cell_num = j.value<int>("dirac_band_cell_num", 3);

		if (j.contains("spray_system")) {
			spray_system = std::make_shared<SpraySystem<d>>(j.at("spray_system"));
		}
	}

	virtual void Initialize(const VectorDi& cell_counts,const real dx,const VectorD& domain_min=VectorD::Zero())
	{
		Base::Initialize(cell_counts,dx,domain_min);
		levelset.Initialize(mac_grid.grid);
		narrow_band_width = mac_grid.grid.dx * narrow_band_cell_num;
		dirac_band_width = mac_grid.grid.dx * dirac_band_cell_num;
		json j_surface = json::value_t::object;
		implicit_surface_tension_solver = std::make_shared<ImplicitSurfaceTension<real, d>>(j_surface, mac_grid, &bc, &levelset);
	}

	virtual void Advance(const real time, const real dt)
	{
		Timer timer;						if(verbose)timer.Reset();
		Advection(dt);						if(verbose)timer.Elapse_And_Output_And_Reset("\t\ttimestep: advection");
		Apply_Body_Forces(dt);				if (verbose)timer.Elapse_And_Output_And_Reset("\t\ttimestep: body force");
		auto proj_ptr = std::dynamic_pointer_cast<ProjectionIrregularCPX<real, d>>(projection);
		if (use_implicit_surface_tension) implicit_surface_tension_solver->Solve(dt, velocity, proj_ptr->sigma, narrow_band_width, dirac_band_width);
		Enforce_Incompressibility(dt);		if(verbose)timer.Elapse_And_Output_And_Reset("\t\ttimestep: projection");
		Extrapolation(velocity);			if(verbose)timer.Elapse_And_Output_And_Reset("\t\ttimestep: extrapolation");
		if (spray_system) {
			spray_system->Advance(time, dt, mac_grid, velocity, levelset);
		}
	}

	virtual void Advection(const real dt)
	{
		Timer timer;if(verbose)timer.Reset();

		////advect velocity
		MacGrid<d> ghost_grid=mac_grid;FaceField<real,d> ghost_velocity=velocity;
		std::function<bool(const int,const VectorDi&)> valid_face=[&](const int axis,const VectorDi& face)->bool{return Is_Fluid_Nb_Face(axis,face);};
		Advection::Semi_Lagrangian<d>(dt,ghost_grid,ghost_velocity,mac_grid,velocity,valid_face);
		if(verbose)timer.Elapse_And_Output_And_Reset("\tadv: vel adv");

		////advect interface
		Field<real,d> ghost_phi=levelset.phi;
		std::function<bool(const VectorDi&)> valid_cell=[&](const VectorDi& cell)->bool{return Is_Fluid_Nb_Cell(cell);};
		Advection::Semi_Lagrangian(dt,ghost_velocity,ghost_grid,ghost_phi,mac_grid,levelset.phi,false,valid_cell);
		if(verbose)timer.Elapse_And_Output_And_Reset("\tadv: ls adv");

		levelset.Fast_Marching(narrow_band_width);
		if(verbose)timer.Elapse_And_Output_And_Reset("\tadv: ls fmm");

		Update_Cell_Types();
		if(verbose)timer.Elapse_And_Output_And_Reset("\tadv: celltype");
	}

	virtual void Update_Cell_Types()
	{
		int cell_num=mac_grid.grid.Number_Of_Cells();
		#pragma omp parallel for
		for(int i=0;i<cell_num;i++){
			if(levelset.phi.array[i]<(real)0)type.array[i]=(ushort)CellType::Fluid;
			else if(levelset.phi.array[i]<narrow_band_width)type.array[i]=(ushort)CellType::NB;
			else type.array[i]=(ushort)CellType::Air;}

		////update boundary cell types
		#pragma omp parallel for
		for(int i=0;i<cell_num;i++){
			if(type.array[i]!=(ushort)CellType::Fluid)continue;
			const VectorDi& cell=mac_grid.grid.Cell_Coord(i);
			for(int j=0;j<mac_grid.grid.Number_Of_Nb_C();j++){
				const VectorDi& nb=mac_grid.grid.Nb_C(cell,j);
				if(!mac_grid.grid.Valid_Cell(nb)||(type(nb)==(ushort)CellType::NB||type(nb)==(ushort)CellType::Air)){
					type.array[i]=(ushort)CellType::BD;}}}
	}

	virtual void Apply_Body_Forces(const real dt)
	{
		if(!use_body_force)return;
		for(int axis=0;axis<d;axis++){int face_num=mac_grid.Number_Of_Faces(axis);
			#pragma omp parallel for
			for(int i=0;i<face_num;i++){VectorDi face=mac_grid.Face_Coord(axis,i);
				VectorD pos=mac_grid.Face_Center(axis,face);real phi=levelset.Phi(pos);
				if(phi<narrow_band_width){velocity.face_fields[axis](face)+=g[axis]*dt;}}}
	}

	virtual void Enforce_Incompressibility(const real dt)
	{
		Enforce_Boundary_Conditions();
		projection->Project(dt, velocity);
	}

	////extrapolation with levelset
	virtual void Extrapolation(FaceField<real,d>& field_q)
	{
		Interpolation<d> intp(mac_grid);FaceField<real,d> ghost_field_q=field_q;

		for (int axis = 0;axis < d;axis++) {
			int face_num = mac_grid.Number_Of_Faces(axis);
#pragma omp parallel for
			for (int i = 0;i < face_num;i++) {
				VectorDi face = mac_grid.Face_Coord(axis, i);
				if (bc.Is_Psi_N(axis, face))continue;

				VectorD pos = mac_grid.Face_Center(axis, face);
				real phi = levelset.Phi(pos);
				if (phi > (real)0) {
					if (phi < narrow_band_width) {
						VectorD interface_pos = levelset.Closest_Point(pos);
						field_q(axis, face) = intp.Interpolate_Faces_With_Phi(ghost_field_q, std::bind(&LevelSet<d>::Phi, &levelset, std::placeholders::_1), interface_pos, axis);
					}
					else field_q(axis, face) = (real)0;
				}
			}
		}
	}

	//////extrapolation with levelset
	//virtual void Extrapolation(FaceField<real,d>& field_q)
	//{
	//	Interpolation<d> intp(mac_grid);FaceField<real,d> ghost_field_q=field_q;

	//	////Fill ghost fields with fluid velocities
	//	for(int axis=0;axis<d;axis++){int face_num=mac_grid.Number_Of_Faces(axis);
	//		#pragma omp parallel for
	//		for(int i=0;i<face_num;i++){VectorDi face=mac_grid.Face_Coord(axis,i);
	//			if(Is_Fluid_Face(axis,face))continue;

	//			int n=0;real fluid_v=(real)0;
	//			for(int i=0;i<Grid<d>::Number_Of_Nb_C();i++){
	//				VectorDi nb_face=Grid<d>::Nb_C(face,i);
	//				if(!mac_grid.face_grids[axis].Valid_Node(nb_face))continue;
	//				if(Is_Fluid_Face(axis,nb_face)){fluid_v+=ghost_field_q(axis,nb_face);n++;}}
	//			if(n>0)ghost_field_q(axis,face)=fluid_v/(real)n;}}

	//	
	//	for(int axis=0;axis<d;axis++){int face_num=mac_grid.Number_Of_Faces(axis);
	//		#pragma omp parallel for
	//		for(int i=0;i<face_num;i++){VectorDi face=mac_grid.Face_Coord(axis,i);
	//			if(bc.Is_Psi_N(axis,face))continue;

	//			VectorD pos=mac_grid.Face_Center(axis,face);real phi=levelset.Phi(pos);
	//			if(phi>(real)0){
	//				if(phi<narrow_band_width){
	//					VectorD interface_pos=levelset.Closest_Point(pos);
	//					field_q(axis,face)=intp.Interpolate_Faces(ghost_field_q,interface_pos,axis);}
	//				else field_q(axis,face)=(real)0;}}}
	//}

protected:
	//////////////////////////////////////////////////////////////////////////
	////helper functions

	inline bool Is_Fluid_Nb_Cell(const VectorDi& cell) const
	{return mac_grid.grid.Valid_Cell(cell)&&((type(cell)==(ushort)CellType::Fluid||type(cell)==(ushort)CellType::BD)||type(cell)==(ushort)CellType::NB);}

	inline bool Is_Fluid_Nb_Cell_Index(const int cell_idx) const 
	{return (type.array[cell_idx]==(ushort)CellType::Fluid||type.array[cell_idx]==(ushort)CellType::BD)||type.array[cell_idx]==(ushort)CellType::NB;}

	////a face is fluid if it is incident to a fluid cell
	inline bool Is_Fluid_Face(const int axis,const VectorDi& face) const
	{{VectorDi cell=MacGrid<d>::Face_Incident_Cell(axis,face,0);if(projection->Is_Fluid_Cell(cell))return true;}
	{VectorDi cell=MacGrid<d>::Face_Incident_Cell(axis,face,1);if(projection->Is_Fluid_Cell(cell))return true;}return false;}

	inline bool Is_Fluid_Nb_Face(const int axis,const VectorDi& face) const
	{{VectorDi cell=MacGrid<d>::Face_Incident_Cell(axis,face,0);if(Is_Fluid_Nb_Cell(cell))return true;}
	{VectorDi cell=MacGrid<d>::Face_Incident_Cell(axis,face,1);if(Is_Fluid_Nb_Cell(cell))return true;}return false;}

	inline bool Is_Interface_Face(const int axis, const VectorDi& face) const
	{
		VectorD pos = mac_grid.Face_Center(axis, face);
		if (bc.Is_Psi_N(axis, face)) return false;
		real phi = levelset.Phi(pos);
		return (phi > -narrow_band_width && phi < narrow_band_width);
	}

	inline bool Is_Interface_Face_Index(const std::pair<int, int> face_idx) const
	{
		int axis = face_idx.first;
		VectorDi face = mac_grid.face_grids[axis].Node_Coord(face_idx.second);
		VectorD pos = mac_grid.Face_Center(axis, face);
		if (bc.Is_Psi_N(axis, face)) return false;
		real phi = levelset.Phi(pos);
		return (phi > -narrow_band_width && phi < narrow_band_width);
	}

	//inline real Dirac(const real phi) const
	//{
	//	if (phi < -dirac_band_width) return 0;
	//	else if (phi > dirac_band_width) return 0;
	//	else return 0.5 * (1.0 + cos(pi * phi / dirac_band_width)) / dirac_band_width;
	//}
};