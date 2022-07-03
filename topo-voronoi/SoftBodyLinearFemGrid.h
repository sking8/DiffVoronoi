//////////////////////////////////////////////////////////////////////////
// Linear Grid FEM
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __SoftBodyLinearFemGrid_h__
#define __SoftBodyLinearFemGrid_h__
#include "SPX_Hashtable.h"
#include "SPX_LinearFemFunc.h"
#include "SPX_BoundaryCondition.h"
#include "SPX_Field.h"
#include "GmgPcgSolverCPU.h"
#include "Simulator.h"
#include "IOHelper.h"
#include "Field.h"
#ifdef USE_CUDA
#include "GmgPcgSolverGPU.h"
#endif

template<int d> class SoftBodyLinearFemGrid
{Typedef_VectorDii(d);
public:
	Grid<d> grid;
	SparseMatrixT K;
	VectorX u;
	VectorX f;

	Array<MatrixX> E;
	Array<MatrixX> Ke;
	Field<short,d> material_id;
	BoundaryConditionGrid<d> bc;
	
	////multigrid and parallelization
	bool use_multigrid_solver=true;
	MultiGrid::Params multigrid_params;
	GMGPCG_Solver_CPU<d> gmg_solver_cpu;
	#ifdef USE_CUDA
	GMGPCG_Solver_GPU<real,d> gmg_solver_gpu;
	#endif

	Array<int> colored_cell_ptr;
	Array<int> colored_cell_indices;
	
	//SoftBodyLinearFemGrid();

	virtual void Initialize(const Grid<d> _grid, const BoundaryConditionGrid<d>& _bc, const Array<std::tuple<real, real>>& _materials, const Field<short, d>& _material_id);
	void Allocate_K();
	virtual void Update_K_And_f();
	virtual void Update_K_And_f(const Meso::Field<real, d>& multiplier);
	virtual void Solve();

	virtual void Output(Meso::DriverMetaData& metadata) {
		std::string vts_name = fmt::format("u_vts{:04d}.vts", metadata.current_frame);
		Meso::Grid<d> meso_grid(grid.node_counts, grid.dx, grid.domain_min, Meso::GridType::CORNER);
		Meso::bf::path vtk_path = metadata.base_path / Meso::bf::path(vts_name);
		Meso::Field<VectorD, d> meso_u(meso_grid);
		meso_u.Calc_Nodes(
			[&](const VectorDi node)->VectorD {
				int idx = grid.Index(node, grid.node_counts);
				if constexpr (d == 2) return Vector2(u[idx * 2], u[idx * 2 + 1]);
				else if constexpr (d == 3) return Vector3(u[idx * 3], u[idx * 3 + 1], u[idx * 3 + 2]);
			}
		);
		Meso::VTKFunc::Write_Vector_Field(meso_u, vtk_path.string());
		/*std::ofstream myfile;
		myfile.open("64_32_cb.txt");
		myfile << u;
		myfile.close();*/
	}

	void Add_Material(const std::tuple<real, real> material);
	void Clear_Materials(){Ke.clear();}

	virtual void Set_Fixed(const VectorDi& node);
	void Set_Displacement(const VectorDi& node,const VectorD& dis);
	virtual void Set_Force(const VectorDi& node,const VectorD& force);

	void Compute_Cell_Displacement(const VectorX& u,const VectorDi& cell,VectorX& cell_u) const;
	void Compute_Cell_Displacement(const VectorX& u,const Array<int>& cell_node_matrix_indices,VectorX& cell_u) const;
    void Compute_Strain(Field<Meso::Matrix<real, d, d, Eigen::ColMajor>,d>& strains) const;
	void Compute_Stress(Field<Meso::Matrix<real, d, d, Eigen::ColMajor>,d>& stresses) const;
	void Compute_Von_Mises_Stress(Field<real,d>& von_mises,const Field<Meso::Matrix<real, d, d, Eigen::ColMajor>,d>* stresses=nullptr) const;
	void Compute_Strain(VectorX& strain,const VectorDi& cell) const;
	void Compute_Strain(Meso::Matrix<real, d, d, Eigen::ColMajor>& strain,const VectorDi& cell) const;
	void Compute_Stress(VectorX& stress,const VectorDi& cell,const MatrixX& E) const;
	void Compute_Stress(Meso::Matrix<real, d, d, Eigen::ColMajor>& stress,const VectorDi& cell,const MatrixX& E) const;

	virtual void Compute_Elastic_Compliance(Field<real, d>& energy);			////multiplying rho
	virtual void Compute_Elastic_Energy(Meso::Field<real, d>& energy) const;			////not multiplying rho

protected:
	int Node_Index_In_K(const VectorDi& node) const;
	void Compute_Cell_Tensor_Helper(const Array<MatrixX>& B,const VectorX& u,Field<Meso::Matrix<real, d, d, Eigen::ColMajor>,d>& tensors) const;
};



#endif
