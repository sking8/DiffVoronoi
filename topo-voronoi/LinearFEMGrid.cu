#include "LinearFEMGrid.h"
#include "ColorGrid.h"

#include <amgcl/backend/eigen.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/solver/cg.hpp>
using namespace Meso;

//=================================Helper function================================
template<int d>
Vector<int,d> Corner_Offset(const Vector<int, d>& center, int i) {
	Assert(i < pow(2, d), "Corner_Offset: index out of  range");
	if constexpr (d == 2) { return center + Vector2i(i & 0x1, (i >> 1) & 0x1); }
	else if constexpr (d == 3) { return center + Vector3i(i & 0x1, (i >> 1) & 0x1, (i >> 2) & 0x1); }
	else { Error("Corner_Offset: dimension not supported"); return Vector<int, d>(); }
}

template<int d> void Vector_To_Field(const VectorX& v, Field<Vector<real,d>,d>& f) {
	Typedef_VectorD(d);
	Assert(f.grid.Counts().prod() * d == v.size(), "Vector_To_Field: vector and field should have the same size");
	f.Exec_Nodes(
		[&](const VectorDi node) {
			int idx = f.grid.Index(node);
			f(node) = v.segment<d>(d * idx);
		}
	);
}

//================================Main body=======================================
template<int d> void LinearFEMGrid<d>::Output(DriverMetaData& metadata) {
	//output the displacement field
	Field<Vector<real, d>, d> u_field;
	u_field.Init(grid);
	Vector_To_Field<d>(u, u_field);
	std::string vts_name = fmt::format("u_vts{:04d}.vts", metadata.current_frame);
	bf::path vtk_path = metadata.base_path / bf::path(vts_name);
	VTKFunc::Write_Vector_Field(u_field, vtk_path.string());
}

template<int d> void LinearFEMGrid<d>::Initialize(const Grid<d> _grid, const BoundaryConditionGrid<d>& _bc, const Array<std::tuple<real,real>>& _materials, const Field<short, d>& _material_id) //this is a corner grid
{
	//colored_cell_ptr.resize(Pow(2, d) + 1);
	//for (int i = 0; i < colored_cell_ptr.size(); i++)colored_cell_ptr[i] = i;

	grid = _grid;
	bc = _bc;
	material_id = _material_id;
	for (int i = 0; i < _materials.size(); i++) { Add_Material(_materials[i]); }
	int n = grid.Counts().prod() * d;								// be careful about the difference between cell counts and node counts
	K.resize(n, n); u.resize(n); u.fill((real)0); f.resize(n); f.fill((real)0);
	Allocate_K();
}

template<int d> void LinearFEMGrid<d>::Add_Material(const std::tuple<real, real> material)
{
	auto [youngs, poisson] = material;
	MatrixX Ke0; //stiffness matrix of the material
	LinearFEMFunc::Cell_Stiffness_Matrix<d>(youngs, poisson, grid.dx, Ke0);
	Ke.push_back(Ke0);
}

template<int d> void LinearFEMGrid<d>::Allocate_K()
{
	std::vector<Triplet<real>> elements; //Fan: Can only use std array 
	
	//need to check which side that data resides
	grid.Iterate_Nodes(
		[&](const VectorDi node) {
		int r = grid.Index(node);		//index here is the index in memory
		for (int nb_r = 0; nb_r < grid.Neighbor_Node_Number(); nb_r++) {
			int c = grid.Index(grid.Neighbor_Ring_Node(node, nb_r));
			if (!grid.Valid(c)) { continue; }
			for (int rr = r * d; rr < (r + 1) * d; rr++)for (int cc = c * d; cc < (c + 1) * d; cc++) { elements.push_back(Triplet<real>(rr, cc, (real)0)); }
		}
		}
	);

	K.setFromTriplets(elements.begin(), elements.end());
	K.makeCompressed();

	//ColorGrid::Color<d>(grid.Counts(), colored_cell_ptr, colored_cell_indices);
}

template<int d> void LinearFEMGrid<d>::Update_K_And_f()
{
	////Update K
//	int color_n = colored_cell_ptr.size()-1;
//	for (int c = 0; c < color_n; c++) {	//Fan: this can be improved as in PoissonFunc.h
//#pragma omp parallel for
//		for (int i = colored_cell_ptr[c]; i < colored_cell_ptr[c + 1]; i++) {
//			VectorDi cell = grid.Cell_Grid().Coord(colored_cell_indices[i]); //Fan: this need to be changed, the colored cell indices are not indices in memory
//			int mat_id = material_id(cell);
//			Array<int> corners(pow(2, d), 0);
//			for (int j = 0; j < corners.size(); j++) {
//				corners[j] = grid.Index(Corner_Offset<d>(cell, j));
//			}
//			LinearFEMFunc::Add_Cell_Stiffness_Matrix<d>(K, Ke[mat_id], corners);
//		}
//	}

	//Fan: not parallelized version, to be changed
	grid.Cell_Grid().Iterate_Nodes(
		[&](const VectorDi node) {
			Array<int> corners(pow(2, d), 0);
			for (int j = 0; j < corners.size(); j++) {
				corners[j] = grid.Index(Corner_Offset<d>(node, j));
			}
			int mat_id = material_id(node);
			LinearFEMFunc::Add_Cell_Stiffness_Matrix<d>(K, Ke[mat_id], corners); //Fan: exception here, may be some indexing issue!
		}
	);

	////Update rhs
	f.fill((real)0);
	for (auto& b : bc.forces) {
		VectorDi node = b.first; VectorD force = b.second;
		for (int axis = 0; axis < d; axis++) { int idx = grid.Index(node) * d + axis; f[idx] += force[axis]; } //Fan: Is the index usage right here. Yes
	}

	////Update bc
	for (auto& b : bc.psi_D_values) {
		VectorDi node = b.first; VectorD dis = b.second;
		for (int axis = 0; axis < d; axis++) {
			int idx = grid.Index(node) * d + axis;			//Fan: is the index usage right here. Yes
			LinearFEMFunc::Set_Dirichlet_Boundary_Helper(K, f, idx, dis[axis]);
		}
	}
}

template<int d> void LinearFEMGrid<d>::Solve()
{
	u.fill((real)0);

//	multigrid_params.use_auto_calculated_levels = true;
//	multigrid_params.dof_on_cell = false;
//	multigrid_params.block_size = d;
//	multigrid_params.use_gpu = true;
//	multigrid_params.init_hier_on_gpu = true;	////calculate hier on CPU to avoid GPU memory crash
//
//#ifdef USE_CUDA
//	if (multigrid_params.use_gpu) {
//		//GeometricMultiGrid::GMGPCG_GPU<d>(K,u,f,grid.node_counts,multigrid_params);
//		gmg_solver_gpu.update_A_levels = true;
//
//		gmg_solver_gpu.Initialize(K, grid.node_counts, multigrid_params, &material_id);
//		gmg_solver_gpu.Solve(u, f);
//	}
//	else {
//		gmg_solver_cpu.update_A_levels = true;
//		gmg_solver_cpu.Initialize(K, grid.node_counts, multigrid_params, &material_id);
//		gmg_solver_cpu.Solve(u, f);
//	}
//#else
//	gmg_solver_cpu.Initialize(K, grid.node_counts, multigrid_params, &material_id);
//	gmg_solver_cpu.Solve(u, f);
//#endif

	// use amgcl solver for solving the system
	// Setup the solver:
	typedef amgcl::make_solver<
		amgcl::amg<
		amgcl::backend::eigen<real>,
		amgcl::coarsening::smoothed_aggregation,
		amgcl::relaxation::spai0
		>,
		amgcl::solver::cg<amgcl::backend::eigen<real> >
	> Solver;

	Solver solve(K);
	std::cout << solve << std::endl;

	// Solve the system for the given RHS:
	auto [a, error]= solve(f, u);
	Info("Solver finished within {} iters with error {}.", a, error);
}

template<int d> void LinearFEMGrid<d>::Compute_Cell_Displacement(const VectorX& u, const VectorDi& cell, VectorX& cell_u) const
{
	int number_of_cell_nodes = pow(2, d);
	cell_u.resize(number_of_cell_nodes * d);
	for (int i = 0; i < number_of_cell_nodes; i++) {
		VectorDi nb_node = Corner_Offset<d>(cell, i); int nb_node_mtx_idx = grid.Index(nb_node);
		for (int j = 0; j < d; j++)cell_u(i * d + j) = u(nb_node_mtx_idx * d + j);
	}
}

template<int d> void LinearFEMGrid<d>::Compute_Elastic_Energy(Field<real, d>& energy) const
{
	energy.Init(grid.Cell_Grid().Counts(), (real)0);

	Grid<d> cell_grid = grid.Cell_Grid();
	cell_grid.Exec_Nodes(
		[&](const VectorDi cell) {
			int mat_id = material_id(cell); if (mat_id == -1)return;
			VectorX cell_u; Compute_Cell_Displacement(u, cell, cell_u);
			const MatrixX& K0 = Ke[mat_id];
			energy(cell) = (real)0.5*cell_u.dot(K0 * cell_u);
		}
	);
}

template class LinearFEMGrid<2>;
template class LinearFEMGrid<3>;
