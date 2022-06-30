//////////////////////////////////////////////////////////////////////////
// Project a vector field with free boundary to divergence free on a MAC grid
// Copyright (c) (2018-),Fan Feng, Shuqi Yang, Bo Zhu
// This file is part of SimpleX,whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifdef USE_CPX

#include "ProjectionTwoPhaseCPX.h"
#include "Timer.h"

//////////////////////////////////////////////////////////////////////////
////projection functions
template<class T, int d> void ProjectionTwoPhaseCPX<T, d>::Apply_Jump_Condition_To_b(const T dt)
{
	if (!apply_jump_condition) return;
	T one_over_dx = 1 / mac_grid->grid.dx;

	int cell_num=mac_grid->grid.cell_counts.prod();
	#pragma omp parallel for
	for (int i = 0; i < cell_num; i++) {
		VectorDi cell = mac_grid->grid.Cell_Coord(i);
		if (!Is_Fluid_Cell(cell))continue;
		real phi0 = (*levelset).phi(cell);
		for (int i = 0; i < Grid<d>::Number_Of_Nb_C(); i++) {
			int axis;VectorDi face;MacGrid<d>::Cell_Incident_Face(cell,i,axis,face);
			VectorDi nb_cell = Grid<d>::Nb_C(cell, i);
			if (!Is_Fluid_Cell(nb_cell))continue;
			T coef=!(*bc).Is_Psi_N(axis,face)?(real)1:(real)0;
			T phi1 = (*levelset).phi(nb_cell);
			if (LevelSet<d>::Interface(phi0, phi1)) {
				real theta = LevelSet<d>::Theta(phi0, phi1);
				Vector<T, d> intf_pos = (1.0 - theta) * mac_grid->grid.Center(cell) + theta * mac_grid->grid.Center(nb_cell);
				real p_sign = phi0 < 0 ? (real)1 : (real)-1;
				cell_b(cell) += coef * p_sign * Pressure_Jump(dt, intf_pos) * one_over_dx / (*rho)(axis, face);
			}
		}
	}
}

////need to specify target_vol and current_vol 
template<class T, int d> void ProjectionTwoPhaseCPX<T, d>::Apply_Vol_Control_To_b()
{
	if(!use_vol_control)return;
	Assert(target_vol > 0, "ProjectionTwoPhaseCPX: target_vol not set");

	if (calc_current_vol)current_vol = levelset->Total_Volume();
	T vol_correction=(target_vol-current_vol)/current_vol;
	if (verbose)Info("vol correction: {}", vol_correction);

	int cell_num=mac_grid->grid.cell_counts.prod();
	#pragma omp parallel for
	for (int i = 0; i < cell_num; i++) {
		VectorDi cell = mac_grid->grid.Cell_Coord(i);
		if (!Is_Liquid_Cell(cell))continue;
		T cell_div = vol_correction;
		cell_b(cell) += vol_control_ks * mac_grid->grid.dx * cell_div;
	}
}

//////////////////////////////////////////////////////////////////////////
////Physical interface functions
//see: https://wmdcstdio.com/2021/07/11/projection-matrix-terms/
template<class T, int d> T ProjectionTwoPhaseCPX<T, d>::Off_Diag_Term(const T dt, const VectorDi& fluid_cell, const int& nbidx) const
{
	VectorDi nb_cell = Grid<d>::Nb_C(fluid_cell, nbidx);
	int axis; VectorDi face; MacGrid<d>::Cell_Incident_Face(fluid_cell, nbidx, axis, face);
	if (bc->Is_Psi_N(axis, face))return 0;
	if (Is_Fluid_Cell(nb_cell)) return -1./(*rho)(axis, face);
	return 0;
}

template<class T, int d> T ProjectionTwoPhaseCPX<T, d>::Diag_Face_Term(const T dt, const int& axis, const VectorDi& face) const
{
	VectorDi cell[2]; for (int i = 0; i < 2; i++)cell[i] = MacGrid<d>::Face_Incident_Cell(axis, face, i);
	if (!Is_Fluid_Cell(cell[0])) std::swap(cell[0], cell[1]);
	if (!Is_Fluid_Cell(cell[0])) return 0;
	if (bc->Is_Psi_N(axis, face)) return 0;
	return 1./(*rho)(axis, face);
}

template<class T, int d> T ProjectionTwoPhaseCPX<T, d>::Velocity_Correction(const T dt, const int& axis, const VectorDi& face) const
{
	if ((*bc).Is_Psi_N(axis, face)) { return 0; }

	VectorDi cell[2]; for (int i = 0; i < 2; i++)cell[i] = MacGrid<d>::Face_Incident_Cell(axis, face, i);
	T cell_p[2]; for (int i = 0; i < 2; i++)cell_p[i] = Is_Fluid_Cell(cell[i]) ? cpx_poisson.x(cell[i]) : (T)0.;		////void cell has zero pressure

	T p_jump = 0.;
	if (apply_jump_condition) {
		if (Is_Fluid_Cell(cell[0]) && Is_Fluid_Cell(cell[1])) {
			T one_over_dx = 1.0 / mac_grid->grid.dx;
			T phi[2]; for (int i = 0; i < 2; i++)phi[i] = (*levelset).phi(cell[i]);
			if (LevelSet<d>::Interface(phi[0], phi[1])){
				real theta = LevelSet<d>::Theta(phi[0], phi[1]);
				VectorD intf_pos = (1 - theta) * mac_grid->grid.Center(cell[0]) + theta * mac_grid->grid.Center(cell[1]);
				T p_sign = phi[0] < 0 ? (T)1 : (T)-1;
				p_jump = p_sign * Pressure_Jump(dt, intf_pos) * one_over_dx;}}}
	
	return -(cell_p[1] + p_jump - cell_p[0]) / (*rho)(axis, face);
}

template class ProjectionTwoPhaseCPX<Scalar, 2>;
template class ProjectionTwoPhaseCPX<Scalar, 3>;

#endif