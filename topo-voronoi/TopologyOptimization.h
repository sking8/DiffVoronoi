#pragma once
#include "Grid.h"
#include "Field.h"
#include "Optimizer.h"
#include "AuxFunc.h"
#include "IOHelper.h"
#include <mma/MMASolver.h>
#include "LinearFEMGrid.h"
#include <thrust/functional.h>

using namespace Meso;
template<int d> class TopologyOptimization : public Optimizer {
	Typedef_VectorD(d); Typedef_MatrixD(d);
public:
	std::shared_ptr<MMASolver>	mma_solver;			//mathematical optimizer
	LinearFEMGrid<d>			linear_fem_grid;	//simulator

	//general optimizer variables
	Array<real>		var;							//variables
	Array<real>		intmed_var;						//intermediate variables
	Array<real>		grad;							//gradients
	Array<real>		var_low_bounds;					//lower bound
	Array<real>		var_up_bounds;					//upper bound
	real			constraint;						//constraint value, only one constraint in this case
	Array<real>		constraint_grads;				//constraint gradients
	real			obj;							//current objective
	Grid<d> grid;
	Field<real, d> rho;								//density field

	////SIMP parameters
	Field<real, d> fem_variable_coef;				//rho^p, the multiplier in front of Young's modulus
	int power;										//penalizing power
	real target_frac;
	real mov_lim;									//move limit for each iteration, can be tuned 
	real rho_min = (real)0;							//the rho_min corresponds to minimum young's modulus = young's modulus * rho_min, optimized rho ranges from [0,1]
	real rho_max = (real)1;
	real beta = 1;									//how sharp the filter is
	real vol_frac;									//current volume fraction

	virtual void Optimize(OptimizerDriverMetaData& meta_data) {
		Compute_Objective();
		Compute_Gradient();
		Compute_Constraint();
		Compute_Constraint_Grad();
		Compute_Bounds();
		intmed_var = var;							//record the intermediate variables before update

		Assert(std::is_same<real, double>::value, "only double data type is supported by the mma solver");
		mma_solver->Update(ArrayFunc::Data(var), ArrayFunc::Data(grad), &constraint, ArrayFunc::Data(constraint_grads), ArrayFunc::Data(var_low_bounds), ArrayFunc::Data(var_up_bounds));
	}

	virtual void Output(OptimizerDriverMetaData& meta_data) {
		if(meta_data.iter_count==0){
			VTKFunc::Write_Boundary_Condition(linear_fem_grid.bc,linear_fem_grid.grid, meta_data.base_path);
			meta_data.data_output << "iter,obj,frac,\n";
		}
		else {
			meta_data.data_output << meta_data.iter_count << "," << obj << "," << vol_frac << ",\n";
		}

		std::string vts_name = fmt::format("vts{:04d}.vts", meta_data.iter_count);
		bf::path vtk_path = meta_data.base_path / bf::path(vts_name);
		VTKFunc::Write_VTS(rho, vtk_path.string());
		
		vts_name = fmt::format("u_vts{:04d}.vts", meta_data.iter_count);
		vtk_path = meta_data.base_path / bf::path(vts_name);
		VTKFunc::Write_Vector_Field(linear_fem_grid.u, vtk_path.string());
	}

	virtual bool Is_Converged(OptimizerDriverMetaData& meta_data) {
		Array<real> temp = var;
		ArrayFunc::Minus(temp, intmed_var);
		real change = ArrayFunc::Max_Abs<real>(temp);

		if (change < meta_data.tol) { return true; }
		return false;
	}

	void Init(const LinearFEMGrid<d>& _linear_fem_grid, real _target_frac,real _mov_lim,int _power) {
		linear_fem_grid = _linear_fem_grid;
		grid = _linear_fem_grid.grid.Cell_Grid();
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
		mma_solver = std::make_shared<MMASolver>(cell_num, 1);	//one constraint
	}

	//================== Essential MMA Functions================================
	void Compute_Objective() { //temporary objective is the sum of all rho
		Sync_Var_Opt_To_Fem();

		grid.Exec_Nodes(
			[&](const VectorDi node) {
				fem_variable_coef(node) = pow(rho(node), power) * ((real)1 - (real)1e-9) + (real)1e-9;
			}
		);

		linear_fem_grid.Update_K(fem_variable_coef);
		linear_fem_grid.Solve();
		Field<real,d> energy_f;
		linear_fem_grid.Compute_Elastic_Energy(energy_f);
		ArrayFunc::Multiply(energy_f.Data(), fem_variable_coef.Data());
		obj = ArrayFunc::Sum(energy_f.Data());
	}

	void Sync_Fem_To_Var_Opt() {
		int counter = 0;
		grid.Iterate_Nodes(
			[&](const VectorDi node) {
				var[counter] = rho(node);
				counter++;
			}
		);
	}

	void Sync_Var_Opt_To_Fem() {
		int counter = 0;
		grid.Iterate_Nodes(
			[&](const VectorDi node) {
				rho(node)= var[counter];
				counter++;
			}
		);
	}

	////var -> dobj_drho
	void Compute_Gradient() {
		Field<real, d> energy_f;
		linear_fem_grid.Compute_Elastic_Energy(energy_f); //may not calculate twice

		int counter = 0;
		grid.Iterate_Nodes(
			[&](const VectorDi node) {
				grad[counter] = -power * (pow(rho(node), power - (real)1)) * ((real)1 - (real)1e-9) * energy_f(node);
				counter++;
			}
		);
	}

	//volume constraints
	void Compute_Constraint() {
		real sum = ArrayFunc::Sum(rho.Data());
		vol_frac = sum / (real)var.size();
		constraint = vol_frac - target_frac;
	}

	void Compute_Constraint_Grad() {
		ArrayFunc::Fill(constraint_grads, (real)1 / constraint_grads.size());
	}

	void Compute_Bounds() {
		var_up_bounds =var;
		var_low_bounds = var;
		ArrayFunc::Add_Scalar(var_up_bounds, mov_lim);
		ArrayFunc::Add_Scalar(var_low_bounds, -mov_lim);
		ArrayFunc::Unary_Transform(var_up_bounds, [=](const real a) {return std::min(a, rho_max); }, var_up_bounds);
		ArrayFunc::Unary_Transform(var_low_bounds, [=](const real a) {return std::max(a, rho_min); }, var_low_bounds);
	}

	////relaxed heaviside projection
	real Hat(const real x) const {
		return (tanh(beta / (real)2) + tanh(beta * (x - (real)0.5))) / tanh(beta / (real)2) / (real)2;
	}
	real Hat_Grad(const real x) const {
		return beta * ((real)1 - pow(tanh(beta * (x - (real)0.5)), 2)) / tanh(beta / (real)2) / (real)2;
	}
};