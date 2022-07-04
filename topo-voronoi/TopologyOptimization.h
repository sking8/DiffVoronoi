#pragma once
#include "Grid.h"
#include "Field.h"
#include "Optimizer.h"
#include "AuxFunc.h"
#include "IOHelper.h"
#include <mma/MMASolver.h>
#include "SoftBodyLinearFemGrid.h"
#include <thrust/functional.h>

template<int d> class TopologyOptimization : public Meso::Optimizer {
	Typedef_VectorD(d);
public:
	std::shared_ptr<MMASolver>			mma_solver;			//mathematical optimizer
	SoftBodyLinearFemGrid<d>			linear_fem_grid;	//simulator

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
	Meso::Field<real, d>	rho;							//density field

	////SIMP parameters
	Meso::Field<real, d>	fem_variable_coef;				//rho^p, the multiplier in front of Young's modulus
	Meso::Field<real, d>	energy_f;
	int power;												//penalizing power
	real target_frac;
	real mov_lim;											//move limit for each iteration, can be tuned 
	real rho_min = (real)1e-3;							//the rho_min corresponds to minimum young's modulus = young's modulus * rho_min, optimized rho ranges from [0,1]
	real rho_max = (real)1;
	real beta = 1;									//how sharp the filter is
	real vol_frac;									//current volume fraction

	virtual void Optimize(Meso::OptimizerDriverMetaData & meta_data) {
		Sync_Var_Opt_To_Fem();
		Compute_Objective();
		Compute_Gradient();
		//Numerical_Derivative_DObj_DRho();
		Compute_Constraint();
		Compute_Constraint_Grad();
		Compute_Bounds();
		intmed_var = var;							//record the intermediate variables before update

		Assert(std::is_same<real, double>::value, "only double data type is supported by the mma solver");
		mma_solver->Update(Meso::ArrayFunc::Data(var), Meso::ArrayFunc::Data(grad),
			&constraint, Meso::ArrayFunc::Data(constraint_grads), Meso::ArrayFunc::Data(var_low_bounds), Meso::ArrayFunc::Data(var_up_bounds));
	}

	virtual void Output(Meso::OptimizerDriverMetaData& meta_data) {
		if(meta_data.iter_count==0){
			//Meso::VTKFunc::Write_Boundary_Condition(linear_fem_grid.bc,grid, meta_data.base_path);
			meta_data.data_output << "iter,obj,frac,\n";
		}
		else {
			Info("obj:{}", obj);
			meta_data.data_output << meta_data.iter_count << "," << obj << "," << vol_frac << ",\n";
		}

		std::string vts_name = fmt::format("vts{:04d}.vts", meta_data.iter_count);
		Meso::bf::path vtk_path = meta_data.base_path / Meso::bf::path(vts_name);
		Meso::VTKFunc::Write_VTS(rho, vtk_path.string());
		
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
	}

	virtual bool Is_Converged(Meso::OptimizerDriverMetaData& meta_data) {
		Meso::Array<real> temp = var;
		Meso::ArrayFunc::Minus(temp, intmed_var);
		real change = Meso::ArrayFunc::Max_Abs<real>(temp);

		if (change < meta_data.tol) { Meso::Pass("Converged!"); return true; }
		return false;
	}

	void Init(const SoftBodyLinearFemGrid<d>& _linear_fem_grid, Meso::Grid<d> _grid, real _target_frac,real _mov_lim,int _power) {
		linear_fem_grid = _linear_fem_grid;
		grid = _grid;
		target_frac = _target_frac;
		mov_lim = _mov_lim;
		power = _power;
		int cell_num = grid.Counts().prod();
		var.resize(cell_num,(real)_target_frac);
		intmed_var.resize(cell_num,(real)0);
		grad.resize(cell_num, (real)1);
		var_up_bounds.resize(cell_num);
		var_low_bounds.resize(cell_num);
		constraint_grads.resize(cell_num);
		fem_variable_coef.Init(grid, (real)0);
		rho.Init(grid,(real)_target_frac);						//padding cells also have _target_frac
		energy_f.Init(grid,(real)0);
		mma_solver = std::make_shared<MMASolver>(cell_num, 1);	//one constraint
	}

	//================== Essential MMA Functions================================
	void Compute_Objective() { //temporary objective is the sum of all rho
		/*grid.Exec_Nodes(
			[&](const VectorDi node) {
				fem_variable_coef(node) = pow(rho(node), power);
			}
		);

		linear_fem_grid.Update_K_And_f(fem_variable_coef);
		linear_fem_grid.Solve();
		linear_fem_grid.Compute_Elastic_Energy(energy_f);
		Meso::ArrayFunc::Multiply(energy_f.Data(), fem_variable_coef.Data());
		obj = Meso::ArrayFunc::Sum(energy_f.Data());*/

		real obj = 0;
		grid.Iterate_Nodes(
			[&](const VectorDi node) {
				if (grid.Index(node) % 2 == 0) { obj += rho(node); }
				//else{ obj += ((real)1-rho(node)); }
			}
		);
		Info("obj:{}", obj);
	}

	void Sync_Var_Opt_To_Fem() {
		grid.Iterate_Nodes(
			[&](const VectorDi node) {
				int idx = grid.Index(node);
				rho(node)= var[idx];
			}
		);
	}

	////var -> dobj_drho
	void Compute_Gradient() {
		//linear_fem_grid.Compute_Elastic_Energy(energy_f); //may not calculate twice

		//grid.Iterate_Nodes(
		//	[&](const VectorDi node) {
		//		int idx = grid.Index(node);
		//		grad[idx] = -power * (pow(rho(node), power - (real)1)) * energy_f(node);
		//	}
		//);
		grid.Exec_Nodes(
			[&](const VectorDi node) {
				int idx = grid.Index(node);
				if (idx % 2 == 0) { grad[idx] = (real)1; }
				else { grad[idx]= (real) - 1; }
			}
		);
	}

	//volume constraints
	void Compute_Constraint() {
		real sum = Meso::ArrayFunc::Sum(rho.Data());
		Info("Sum:{}", sum);
		vol_frac = sum / (real)var.size();
		constraint = vol_frac - target_frac;
	}

	void Compute_Constraint_Grad() {
		Meso::ArrayFunc::Fill(constraint_grads, (real)1 / constraint_grads.size());
	}

	void Compute_Bounds() {
		var_up_bounds = var;
		var_low_bounds = var;
		Meso::ArrayFunc::Add_Scalar(var_up_bounds, mov_lim);
		Meso::ArrayFunc::Add_Scalar(var_low_bounds, -mov_lim);
		Meso::ArrayFunc::Unary_Transform(var_up_bounds, [=](const real a) {return std::min(a, rho_max); }, var_up_bounds);
		Meso::ArrayFunc::Unary_Transform(var_low_bounds, [=](const real a) {return std::max(a, rho_min); }, var_low_bounds);
	}

	////relaxed heaviside projection
	real Hat(const real x) const {
		return (tanh(beta / (real)2) + tanh(beta * (x - (real)0.5))) / tanh(beta / (real)2) / (real)2;
	}
	real Hat_Grad(const real x) const {
		return beta * ((real)1 - pow(tanh(beta * (x - (real)0.5)), 2)) / tanh(beta / (real)2) / (real)2;
	}

	void Numerical_Derivative_DObj_DRho()
	{
		real delta_rho = (real)1e-3;

		Info("Numerical DObj_DRho");
		Array<real> numerical_dobj_drho(grid.Counts().prod(),(real)0);

		////should not parallelize here
		int counter = 0;
		grid.Iterate_Nodes(
			[&](const VectorDi cell) {
				real old_rho = rho(cell);
				real old_obj = obj;
				rho(cell) += delta_rho;
				Compute_Objective();
				real new_obj = obj;
				numerical_dobj_drho[counter] = (new_obj - old_obj)/delta_rho;
				rho(cell) = old_rho;
				if (!Meso::MathFunc::Close(numerical_dobj_drho[counter], grad[counter], (real)1e-4, (real)1e-3)) {
					Meso::Warn("cell:{}, analytical:{}, numerical:{}", cell, grad[counter], numerical_dobj_drho[counter]);
				}
				counter++;
			}
		);

		Meso::Pass("Finished test of dobj drho");
	}
};