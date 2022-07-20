//////////////////////////////////////////////////////////////////////////
// Topology Optimization with Voronoi Structure
// Copyright (c) (2022-), Fan Feng
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "Grid.h"
#include "Field.h"
#include "Optimizer.h"
#include "AuxFunc.h"
#include "IOHelper.h"
#include <mma/MMASolver.h>
#include "SoftBodyLinearFemGrid.h"
#include <thrust/functional.h>
#include "stencil.h"
#include "VoronoiField.h"

template<int d> class TopoOptVoronoi : public Meso::Optimizer {
	Typedef_VectorD(d);
public:
	std::shared_ptr<MMASolver>			mma_solver;			//mathematical optimizer
	SoftBodyLinearFemGrid<d>			linear_fem_grid;	//simulator
	Meso::VoronoiField<d>				voronoi_field;
	//general optimizer variables
	Meso::Array<real>		var;							//variables
	Meso::Array<real>		intmed_var;						//intermediate variables
	Meso::Array<real>		grad;							//gradients
	Meso::Array<real>		var_low_bounds;					//lower bound
	Meso::Array<real>		var_up_bounds;					//upper bound
	real					constraint;						//constraint value, only one constraint in this case
	Meso::Array<real>		constraint_grads;				//constraint gradients
	real					obj;							//current objective
	Meso::Grid<d>			grid;							//center grid
	Meso::Field<real, d>	hat_rho;						//density field
	Meso::Field<real, d>	dobj_dhat_rho;						//intermediate field dobj/dhat_rho

	////SIMP parameters
	Meso::Field<real, d>	fem_variable_coef;				//hat_rho^p, the multiplier in front of Young's modulus
	Meso::Field<real, d>	energy_f;
	int power;												//penalizing power
	real target_frac;
	real mov_lim;											//move limit for each iteration, can be tuned
	real rho_min = 1e-2;									//minimum rho allowed
	real rho_max = 1;									//maximum rho allowed
	real pos_min;											//minimum pos component
	real pos_max;											//maximum pos component
	real beta = 8;											//how sharp the filter is
	real vol_frac;											//current volume fraction

	virtual void Optimize(Meso::OptimizerDriverMetaData& meta_data) {
		Sync_Var_To_Voronoi();
		Compute_Objective(meta_data);
		Compute_Gradient();
		//Numerical_Derivative_Dobj_Dx();
		Compute_Constraint();
		Compute_Constraint_Grad();
		//Numerical_Derivative_Dconstraint_Dx();
		Compute_Bounds();
		intmed_var = var;							//record the intermediate variables before update

		mma_solver->Update(Meso::ArrayFunc::Data(var), Meso::ArrayFunc::Data(grad),
			&constraint, Meso::ArrayFunc::Data(constraint_grads), Meso::ArrayFunc::Data(var_low_bounds), Meso::ArrayFunc::Data(var_up_bounds));
	}

	virtual void Output(Meso::OptimizerDriverMetaData& meta_data) {
		if (meta_data.iter_count == 0) {
			//Meso::VTKFunc::Write_Boundary_Condition(linear_fem_grid.bc,grid, meta_data.base_path);
			meta_data.data_output << "iter,obj,frac,\n";
		}
		else {
			meta_data.data_output << meta_data.iter_count << "," << obj << "," << vol_frac << ",\n";
		}

		std::string vts_name = fmt::format("vts{:04d}.vts", meta_data.iter_count);
		Meso::bf::path vtk_path = meta_data.base_path / Meso::bf::path(vts_name);
		Meso::VTKFunc::Write_VTS(hat_rho, vtk_path.string());

		Grid<d> spx_grid(grid.Counts(), grid.dx, grid.Domain_Min(Meso::CENTER));
		Meso::Field<VectorD, d> meso_u(grid);
		meso_u.Calc_Nodes(
			[&](const VectorDi node)->VectorD {
				int idx = spx_grid.Index(node, spx_grid.node_counts);
				if constexpr (d == 2) return Vector2(linear_fem_grid.u[idx * 2], linear_fem_grid.u[idx * 2 + 1]);
				else if constexpr (d == 3) return Vector3(linear_fem_grid.u[idx * 3], linear_fem_grid.u[idx * 3 + 1], linear_fem_grid.u[idx * 3 + 2]);
			}
		);
		vts_name = fmt::format("u_vts{:04d}.vts", meta_data.iter_count);
		vtk_path = meta_data.base_path / Meso::bf::path(vts_name);
		Meso::VTKFunc::Write_Vector_Field(meso_u, vtk_path.string());

		std::string vtu_name = fmt::format("points{:04d}.vtu", meta_data.iter_count);
		Meso::bf::path vtu_path = meta_data.base_path / Meso::bf::path(vtu_name);
		Meso::VTKFunc::Write_Points<d>(voronoi_field.particles.xRef(), vtu_path.string());
	}

	virtual bool Is_Converged(Meso::OptimizerDriverMetaData& meta_data) {
		Meso::Array<real> temp = var;
		Meso::ArrayFunc::Minus(temp, intmed_var);
		real change = Meso::ArrayFunc::Max_Abs<real>(temp);

		if (change < meta_data.tol && meta_data.iter_count!=1) { Meso::Pass("Converged!"); return true; }
		return false;
	}

	void Init(const SoftBodyLinearFemGrid<d>& _linear_fem_grid,const Meso::VoronoiField<d>& _voronoi_field, real _target_frac, real _mov_lim, int _power) {
		linear_fem_grid = _linear_fem_grid;
		voronoi_field = _voronoi_field;
		grid = _voronoi_field.grid;
		target_frac = _target_frac;
		mov_lim = _mov_lim;
		pos_min = (real)2*grid.Domain_Min(Meso::CENTER)[0] - grid.Domain_Max(Meso::CENTER)[0]; //one more domain space
		pos_max = (real)2*grid.Domain_Max(Meso::CENTER)[0] - grid.Domain_Min(Meso::CENTER)[0];
		power = _power;
		int var_num = voronoi_field.particles.Size()*d;
		var.resize(var_num);
		//initialize the vars
		voronoi_field.particles.Exec_Points(
			[&](const int idx) {
				for (int dim = 0; dim < d; dim++) {
					var[idx * d + dim] = voronoi_field.particles.x(idx)[dim];
				}
			}
		);
		intmed_var.resize(var_num, (real)0);
		grad.resize(var_num, (real)1);
		var_up_bounds.resize(var_num);
		var_low_bounds.resize(var_num);
		constraint_grads.resize(var_num);
		fem_variable_coef.Init(grid, (real)0);
		hat_rho.Init(grid);
		dobj_dhat_rho.Init(grid);
		energy_f.Init(grid, (real)0);
		Assert(std::is_same<real, double>::value, "only double data type is supported by the mma solver");
		mma_solver = std::make_shared<MMASolver>(var_num, 1);	//one constraint
		Compute_Objective(Meso::OptimizerDriverMetaData());
	}

	//================== Essential MMA Functions================================
	void Compute_Objective(Meso::OptimizerDriverMetaData& meta_data) { //temporary objective is the sum of all rho
		Update_Hat_Rho(meta_data);
		grid.Exec_Nodes(
			[&](const VectorDi cell) {
				fem_variable_coef(cell) = pow(hat_rho(cell), (real)power);
			}
		);
		linear_fem_grid.Update_K_And_f(fem_variable_coef);
		linear_fem_grid.Solve();
		linear_fem_grid.Compute_Elastic_Energy(energy_f);
		Meso::ArrayFunc::Multiply(fem_variable_coef.Data(), energy_f.Data());
		obj = Meso::ArrayFunc::Sum(fem_variable_coef.Data());
		Meso::Info("Current objective: {}", obj);
	}

	void Update_Hat_Rho(Meso::OptimizerDriverMetaData& meta_data) {
		//update Voronoi points
		voronoi_field.Update_A();
		voronoi_field.Update_Neighbors();
		voronoi_field.Update_Softmax_Sum();
		voronoi_field.Update_Rho();

		if (meta_data.iter_count != 0 && meta_data.iter_count % 50 == 0) {
			beta *= 2;
			Info("beta now is: ", beta);
		}

		grid.Exec_Nodes(
			[&](const VectorDi cell) {
				hat_rho(cell) = Hat(voronoi_field.rho(cell));
			}
		);
	}

	void Sync_Var_To_Voronoi() {
		voronoi_field.particles.Exec_Points(
			[&](const int idx) {
				for (int dim = 0; dim < d; dim++) {
					voronoi_field.particles.x(idx)[dim] = var[idx*d+dim];
				}
			}
		);
	}

	////var -> dobj_drho
	void Compute_Gradient() {
		grid.Exec_Nodes(
			[&](const VectorDi cell) {
				////elastic objective
				dobj_dhat_rho(cell) = -(real)power * (pow(hat_rho(cell), (real)power - (real)1)) * energy_f(cell);
			}
		);

		voronoi_field.Update_DRho_DX();

		Meso::ArrayFunc::Fill(grad, (real)0);
		//should not parallelize here
		grid.Iterate_Nodes(
			[&](const VectorDi cell) {
				int nb_num = voronoi_field.nbs_c(cell).size();
				for (int nb_idx = 0; nb_idx < nb_num; nb_idx++) { // jth particle neighbor of the cell
					int nb_p = voronoi_field.nbs_c(cell)[nb_idx];
					VectorD dobj_dx = Hat_Grad(voronoi_field.rho(cell)) * voronoi_field.drho_dx(cell)[nb_idx]*dobj_dhat_rho(cell);
					for (int dim = 0; dim < d; dim++) { grad[d * nb_p + dim] += dobj_dx[dim]; }
				}
			}
		);
	}

	//volume constraints
	void Compute_Constraint() {
		real sum = Meso::ArrayFunc::Sum(hat_rho.Data());
		vol_frac = sum / (real)grid.Counts().prod();
		Meso::Info("Current fraction: {}", vol_frac);
		constraint = vol_frac - target_frac;
	}

	void Compute_Constraint_Grad() {
		Meso::ArrayFunc::Fill(constraint_grads, (real)0);
		//should not parallelize here
		grid.Iterate_Nodes(
			[&](const VectorDi cell) {
				int nb_num = voronoi_field.nbs_c(cell).size();
				for (int nb_idx = 0; nb_idx < nb_num; nb_idx++) { // jth particle neighbor of the cell
					int nb_p = voronoi_field.nbs_c(cell)[nb_idx];
					VectorD dvol_dx = Hat_Grad(voronoi_field.rho(cell)) * voronoi_field.drho_dx(cell)[nb_idx];
					for (int dim = 0; dim < d; dim++) { constraint_grads[d * nb_p + dim] += dvol_dx[dim] / (real)grid.Counts().prod(); }
				}
			}
		);
	}

	void Compute_Bounds() {
		var_up_bounds = var;
		var_low_bounds = var;
		Meso::ArrayFunc::Add_Scalar(var_up_bounds, mov_lim);
		Meso::ArrayFunc::Add_Scalar(var_low_bounds, -mov_lim);
		Meso::ArrayFunc::Unary_Transform(var_up_bounds, [=](const real a) {return std::min(a, pos_max); }, var_up_bounds);
		Meso::ArrayFunc::Unary_Transform(var_low_bounds, [=](const real a) {return std::max(a, pos_min); }, var_low_bounds);
	}

	////relaxed heaviside projection
	real Hat(const real x) const {
		if (x < 0) { return rho_min; }
		else if (x > 1) { return rho_max; }
		real threshol_proj=(tanh(beta / (real)2) + tanh(beta * (x - (real)0.5))) / tanh(beta / (real)2) / (real)2;
		return rho_min+threshol_proj*(rho_max - rho_min);

	}

	real Hat_Grad(const real x) const {
		if (x < 0 || x > 1 ) { return 0; }
		real threshol_proj_grad = beta * ((real)1 - pow(tanh(beta * (x - (real)0.5)), 2)) / tanh(beta / (real)2) / (real)2;
		return threshol_proj_grad * (rho_max - rho_min);
	}

	void Numerical_Derivative_Dobj_Dx()
	{
		Info("Numerical Dobj_Dx");
		real delta_x = (real)1e-5;
		Array<VectorD> numerical_dobj_dx(voronoi_field.particles.Size(), VectorD::Zero());
		real old_obj = obj;

		////should not parallelize here
		voronoi_field.particles.Iterate_Points(
			[&](const int idx) {
				VectorD old_x = voronoi_field.particles.x(idx);
				for(int dim=0;dim<d;dim++){
					voronoi_field.particles.x(idx)[dim] += delta_x;
					Compute_Objective(Meso::OptimizerDriverMetaData());
					real new_obj = obj;
					numerical_dobj_dx[idx][dim] = (new_obj - old_obj) / delta_x;
					voronoi_field.particles.x(idx) = old_x;
				}
				VectorD anatytial_dobj_dx = Meso::MathFunc::V<d>(grad[idx * d], grad[idx * d + 1], grad[idx * d + 2]);
				if (!Meso::MathFunc::All_Close(numerical_dobj_dx[idx], anatytial_dobj_dx, (real)1e-2, (real)1e-3)) {
					Meso::Warn("point:{},{}, analytical:{}, numerical:{}", idx,voronoi_field.particles.x(idx), anatytial_dobj_dx, numerical_dobj_dx[idx]);
				}
			}
		);

		Meso::Pass("Finished test of dobj dx");
	}

	void Numerical_Derivative_Dconstraint_Dx() {
		Info("Numerical Dconstraint_Dx");
		real delta_x = (real)1e-5;
		Array<VectorD> numerical_dconstraint_dx(voronoi_field.particles.Size(), VectorD::Zero());
		real old_constraint = constraint;

		////should not parallelize here
		voronoi_field.particles.Iterate_Points(
			[&](const int idx) {
				VectorD old_x = voronoi_field.particles.x(idx);
				for (int dim = 0; dim < d; dim++) {
					voronoi_field.particles.x(idx)[dim] += delta_x;
					Update_Hat_Rho(Meso::OptimizerDriverMetaData());
					Compute_Constraint();
					real new_constraint = constraint;
					numerical_dconstraint_dx[idx][dim] = (new_constraint - old_constraint) / delta_x;
					voronoi_field.particles.x(idx) = old_x;
				}
				VectorD anatytial_dconstraint_dx = Meso::MathFunc::V<d>(constraint_grads[idx * d], constraint_grads[idx * d + 1], constraint_grads[idx * d + 2]);
				if (!Meso::MathFunc::All_Close(numerical_dconstraint_dx[idx], anatytial_dconstraint_dx, (real)1e-2, (real)1e-5)) {
					Meso::Warn("point:{},{}, analytical:{}, numerical:{}", idx, voronoi_field.particles.x(idx), anatytial_dconstraint_dx, numerical_dconstraint_dx[idx]);
				}
			}
		);

		Meso::Pass("Finished test of dconstraint dx");
	}
};