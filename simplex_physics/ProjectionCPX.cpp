//////////////////////////////////////////////////////////////////////////
// Project a vector field to divergence free
// Copyright (c) (2018-),Bo Zhu, Mengdi Wang
// This file is part of SimpleX,whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifdef USE_CPX
#include "ProjectionCPX.h"
#include "Timer.h"

template<class T, int d> void ProjectionCPX<T, d>::Build_A(const T dt)
{
	if (A_initialized && !update_A)return;
	A_initialized = true;

	int n = mac_grid->grid.cell_counts[0];
	//if (n % 4) AuxFunc::Crash_With_Info("ProjectionCPX::Prepare_CPX_System: cell size must can be divided by 8");

	//weird here, because cpx can only be initialized once
	//TODO: fix this
	if (!cpx_initialized) {
		cpx_initialized = true;

		cpx_poisson.Init(*mac_grid, max_iter, relative_tolerance);
		cpx_poisson.cg.verbose = verbose;
	}
	
	cpx_poisson.Update_Unknown(std::bind(&ProjectionCPX<T, d>::Is_Fluid_Cell, this, std::placeholders::_1));
	cpx_poisson.Update_Vol(std::bind(&ProjectionCPX<T, d>::Diag_Face_Term, this, dt, std::placeholders::_1, std::placeholders::_2));
	cpx_poisson.Send_To_Device();
	//cpx_poisson.closed = false;		//NOTE: set to not closed here
}

template<class T, int d> void ProjectionCPX<T, d>::Build_b(const T dt, const FaceField<T,d> &velocity)
{
	//cpx_poisson.Update_RHS<T>(
	//	[&](const VectorDi& cell) {
	//		if (!Is_Fluid_Cell(cell)) return (T)0;
	//		T div = (T)0;
	//		for (int axis = 0; axis < d; axis++) {
	//			div += velocity(axis, cell + VectorDi::Unit(axis)) - velocity(axis, cell);
	//		}
	//		////Attention: use negative div here. We are solving -lap p=-div u
	//		return -div;
	//		//cell_b(cell) = -div;
	//	}
	//);

	cell_b.Resize(mac_grid->grid.cell_counts,(T)0);

	int cell_num = mac_grid->grid.Number_Of_Cells();
	#pragma omp parallel for
	for (int i = 0; i < cell_num; i++) {
		VectorDi cell = mac_grid->grid.Cell_Coord(i);
		if (!Is_Fluid_Cell(cell)) continue;
		T div = (T)0;
		for (int axis = 0; axis < d; axis++) {
			div += velocity(axis, cell + VectorDi::Unit(axis)) - velocity(axis, cell);
		}
		////Attention: use negative div here. We are solving -lap p=-div u
		cell_b(cell) = -div;
	}
}

template<class T, int d> void ProjectionCPX<T, d>::Build(const T dt, const FaceField<T, d>& velocity)
{
	Timer timer;timer.Reset();
	Build_A(dt);
	if(verbose)timer.Elapse_And_Output_And_Reset("Update A");
	Build_b(dt, velocity);
	if(verbose)timer.Elapse_And_Output_And_Reset("Update b");
}

template<class T, int d> void ProjectionCPX<T, d>::Solve(const T dt)
{
	Assert(A_initialized && cpx_initialized, "ProjectionCPX<T, d>::Solve: uninitialized");
	cpx_poisson.Update_b(cell_b);
	cpx_poisson.Solve();
}

template<class T, int d> void ProjectionCPX<T, d>::Correction(const T dt, FaceField<T, d>& velocity)
{
	for (int axis = 0; axis < d; axis++) {
		int face_num = mac_grid->face_grids[axis].node_counts.prod();
#pragma omp parallel for
		for (int i = 0; i < face_num; i++) {
			VectorDi face = mac_grid->face_grids[axis].Node_Coord(i);
			velocity(axis, face) += Velocity_Correction(dt, axis, face);
		}
	}
}

template<class T, int d> void ProjectionCPX<T, d>::Project(const T dt, FaceField<T, d>& velocity)
{
	Timer timer;					if (verbose)timer.Reset();
	Build(dt, velocity);			if (verbose)timer.Elapse_And_Output_And_Reset("build");
	Solve(dt);						if (verbose)timer.Elapse_And_Output_And_Reset("solve");
	Correction(dt,velocity);		if (verbose)timer.Elapse_And_Output_And_Reset("correction");
}

//////////////////////////////////////////////////////////////////////////
////Physical interface functions that defines the problem
//see: https://wmdcstdio.com/projection-matrix-terms/
template<class T, int d> T ProjectionCPX<T, d>::Off_Diag_Term(const T dt, const VectorDi& fluid_cell,const int& nb_idx) const
{
	VectorDi nb_cell=Grid<d>::Nb_C(fluid_cell,nb_idx);
	int axis;VectorDi face;MacGrid<d>::Cell_Incident_Face(fluid_cell,nb_idx,axis,face);
	if(bc&&bc->Is_Psi_N(axis,face)) return 0;
	if(Is_Fluid_Cell(nb_cell)) return -1;
	return 0;
}

template<class T, int d> T ProjectionCPX<T, d>::Diag_Face_Term(const T dt, const int& axis,const VectorDi& face) const
{
	VectorDi cell[2];for(int i=0;i<2;i++)cell[i]=MacGrid<d>::Face_Incident_Cell(axis,face,i);
	if(!Is_Fluid_Cell(cell[0])) std::swap(cell[0],cell[1]);
	if(!Is_Fluid_Cell(cell[0])) return 0;
	if(bc&&bc->Is_Psi_N(axis,face)) return 0;
	return 1;
}

template<class T, int d> T ProjectionCPX<T, d>::Velocity_Correction(const T dt, const int& axis, const VectorDi& face) const
{
	if (bc&&bc->Is_Psi_N(axis, face)) return 0;
	VectorDi cell[2]; for (int i = 0; i < 2; i++)cell[i] = MacGrid<d>::Face_Incident_Cell(axis, face, i);
	T cell_p[2]; for (int i = 0; i < 2; i++)cell_p[i] = Is_Fluid_Cell(cell[i]) ? cpx_poisson.x(cell[i]) : 0;
	return -(cell_p[1] - cell_p[0]);
}

template class ProjectionCPX<Scalar, 2>;
template class ProjectionCPX<Scalar, 3>;

#endif