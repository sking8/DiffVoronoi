#pragma once
#include "Grid.h"
#include "Field.h"
#include "Optimizer.h"
#include "AuxFunc.h"
#include "IOHelper.h"
#include <mma/MMASolver.h>

using namespace Meso;
template<int d> class TopologyOptimization : public Optimizer {
	Typedef_VectorD(d); Typedef_MatrixD(d);
public:
	std::shared_ptr<MMASolver> mma_solver;			//mathematical optimizer
	
	//general optimizer variables
	Array<real>		var;							//variables
	Array<real>		intmed_var;						//intermediate variables
	Array<real>		grad;							//gradients
	Array<real>		var_low_bounds;					//lower bound
	Array<real>		var_up_bounds;					//upper bound
	real			constraint;						//constraint value, only one constraint in this case
	Array<real>		constraint_grads;				//constraint gradients

	Grid<d> grid;
	Field<real, d> rho;								//density field

	////SIMP parameters
	real mov_lim = (real)0.05;						//move limit for each iteration, can be tuned 
	real rho_min = (real)0;							//the rho_min corresponds to minimum young's modulus = young's modulus * rho_min, optimized rho ranges from [0,1]
	real rho_max = (real)1;

	void Optimize(OptimizerDriverMetaData& meta_data) {
		Compute_Objective();
		Compute_Gradient();
		//Compute_Constraint();
		//Compute_Constraint_Grad();
		Compute_Bounds();
		intmed_var = var;							//record the intermediate variables before update

		Assert(std::is_same<real, double>::value, "only double data type is supported by the mma solver");
		mma_solver->Update(ArrayFunc::Data<real, DataHolder::HOST>(var), ArrayFunc::Data<real, DataHolder::HOST>(grad), nullptr, nullptr, ArrayFunc::Data<real, DataHolder::HOST>(var_low_bounds), ArrayFunc::Data<real, DataHolder::HOST>(var_up_bounds));
	}

	void Output(const bf::path base_path, const int iter) {
		std::string vts_name = fmt::format("vts{:04d}.vts", iter);
		bf::path vtk_path = base_path / bf::path(vts_name);
		VTKFunc::Write_VTS(rho, vtk_path.string());
	}

	bool Is_Converged(OptimizerDriverMetaData& meta_data) {
		Array<real> temp = var;
		ArrayFunc::Minus(temp, intmed_var);
		real change = ArrayFunc::Max_Abs<real>(temp);

		if (change < meta_data.tol) { return true; }
		return false;
	}

	void Init(Grid<d>& _grid) {
		grid = _grid;
		int cell_num = grid.counts.prod();
		var.resize(cell_num,(real)1);
		intmed_var.resize(cell_num,(real)0); //Fan: be careful with the initialization of the intermediate variable value
		grad.resize(cell_num, (real)1);
		var_up_bounds.resize(cell_num);
		var_low_bounds.resize(cell_num);
		rho.Init(grid,(real)1);
		mma_solver = std::make_shared<MMASolver>(cell_num, 0);
	}

	//================== Essential MMA Functions================================

	real Compute_Objective() { //temporary objective is the sum of all rho
		Sync_Var_Opt_To_Fem();
		real obj = 0;
		for (int i = 0; i < var.size(); i++) { obj += var[i]; }
		return obj;
	}

	void Sync_Fem_To_Var_Opt() {
#pragma omp parallel for
		for (int i = 0; i < var.size(); i++) { var[i] = rho.Data()[i]; }
	}

	void Sync_Var_Opt_To_Fem() {
#pragma omp parallel for
		for (int i = 0; i < var.size(); i++) { rho.Data()[i] = var[i]; }
	}

	////var -> dobj_drho
	void Compute_Gradient() {
#pragma omp parallel for
		for (int i = 0; i < var.size(); i++) { grad[i] = (real)1; }
	}

	//volume constraints
	void Compute_Constraint() {

	}

	void Compute_Constraint_Grad() {
		
	}

	void Compute_Bounds() {
#pragma omp parallel for
		for (int i = 0; i < var.size(); i++) {
			var_up_bounds[i] = std::min(rho_max, var[i] + mov_lim);
			var_low_bounds[i] = std::max(rho_min, var[i] - mov_lim);
		}
	}
};