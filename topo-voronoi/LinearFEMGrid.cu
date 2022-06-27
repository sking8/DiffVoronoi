#include "LinearFEMGrid.h"
#include <amgcl/backend/eigen.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/solver/cg.hpp>
using namespace Meso;

template<int d> int LinearFEMGrid<d>::Compact_Idx(int idx) {
	return compact_indices[idx];
}

//================================Main body=======================================
template<int d> void LinearFEMGrid<d>::Output(DriverMetaData& metadata) {
	//output the displacement field
	std::string vts_name = fmt::format("u_vts{:04d}.vts", metadata.current_frame);
	bf::path vtk_path = metadata.base_path / bf::path(vts_name);
	VTKFunc::Write_Vector_Field(u, vtk_path.string());
	VTKFunc::Write_Boundary_Condition(bc, grid, metadata.base_path);
}

template<int d> void LinearFEMGrid<d>::Initialize(const Grid<d> _grid, const BoundaryConditionGrid<d>& _bc, const Array<std::tuple<real,real>>& _materials, const Field<short, d>& _material_id) //this is a corner grid
{
	grid = _grid;
	bc = _bc;
	material_id = _material_id;
	for (int i = 0; i < _materials.size(); i++) { Add_Material(_materials[i]); }

	u.Init(grid, VectorD::Zero());
	f.Init(grid, VectorD::Zero());
	int n = grid.Counts().prod() * d;	// be careful about the difference between cell counts and node counts
	compact_indices.resize(grid.Memory_Size(), (int)-1);
	
	K.resize(n, n);

	//create the mapping
	int counter = 0;
	grid.Iterate_Nodes(
		[&](const VectorDi node) {
			if (grid.Valid(node)) { 
				compact_indices[grid.Index(node)] = counter; 
				counter++;
			}
		}
	);
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
		int r = Compact_Idx(grid.Index(node));
		for (int nb_r = 0; nb_r < grid.Neighbor_Ring_Number(); nb_r++) {
			VectorDi ring_node = Grid<d>::Neighbor_Ring_Node(node, nb_r);
			int c = grid.Index(ring_node);
			if (!grid.Valid(ring_node)) { continue; }
			c = Compact_Idx(c);
			for (int rr = r * d; rr < (r + 1) * d; rr++) {
				for (int cc = c * d; cc < (c + 1) * d; cc++) { 
					elements.push_back(Triplet<real>(rr, cc, (real)0)); 
				}
			}
		}
		}
	);

	K.setFromTriplets(elements.begin(), elements.end());
	K.makeCompressed();
}

template<int d> void LinearFEMGrid<d>::Update_K()
{
	//Update K
	//Fan: not parallelized version, to be changed
	grid.Cell_Grid().Iterate_Nodes(
		[&](const VectorDi node) {
			Array<int> corners(pow(2, d), 0);
			for (int j = 0; j < corners.size(); j++) {
				corners[j] = Compact_Idx(grid.Index(LinearFEMFunc::Corner_Offset<d>(node, j)));
			}
			int mat_id = material_id(node);
			LinearFEMFunc::Add_Cell_Stiffness_Matrix<d>(K, Ke[mat_id], corners);
		}
	);
}

template<int d> void LinearFEMGrid<d>::Solve()
{
	//VectorX f_v, u_v;
	VectorX f_v(K.cols());
	VectorX u_v(K.cols());
	f_v.fill((real)0);
	u_v.fill((real)0);

	grid.Exec_Nodes(
		[&](const VectorDi node) {
			int idx=Compact_Idx(grid.Index(node));
			if (idx) { f_v.segment<d>(idx * d) = f(node); }
		}
	);

	////Update rhs
	for (auto& b : bc.forces) {
		VectorDi node = b.first; VectorD force = b.second;
		int idx = Compact_Idx(grid.Index(node));
		f_v.segment<d>(idx*d) = force;
	}

	//Update bc
	for (auto& b : bc.psi_D_values) {
		VectorDi node = b.first; VectorD dis = b.second;
		for (int axis = 0; axis < d; axis++) {
			int idx = Compact_Idx(grid.Index(node)) * d + axis;		//need to be careful here when the grid has padding
			LinearFEMFunc::Set_Dirichlet_Boundary_Helper(K, f_v, idx, dis[axis]);
		}
	}
	
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

	Solver solve(K); //K is singular in 3D
	std::cout << solve << std::endl;

	// Solve the system for the given RHS:
	auto [iter, error]= solve(f_v, u_v);
	Info("AMGCL solver with Eigen backend finished within {} iters with error {}.", iter, error);

	grid.Exec_Nodes(
		[&](const VectorDi node) {
			int idx = Compact_Idx(grid.Index(node));
			if (idx) { u(node) = u_v.segment<d>(idx * d); }
		}
	);
}

template<int d> void LinearFEMGrid<d>::Compute_Cell_Displacement(const VectorDi& cell, VectorX& cell_u) const
{
	int number_of_cell_nodes = pow(2, d);
	cell_u.resize(number_of_cell_nodes * d);
	for (int i = 0; i < number_of_cell_nodes; i++) {
		VectorDi nb_node = LinearFEMFunc::Corner_Offset<d>(cell, i);
		cell_u.segment<d>(i*d) = u(nb_node);
	}
}

template<int d> void LinearFEMGrid<d>::Compute_Elastic_Energy(Field<real, d>& energy) const
{
	energy.Init(grid.Cell_Grid().Counts(), (real)0);

	Grid<d> cell_grid = grid.Cell_Grid();
	cell_grid.Exec_Nodes(
		[&](const VectorDi cell) {
			int mat_id = material_id(cell); if (mat_id == -1)return;
			VectorX cell_u; Compute_Cell_Displacement(cell, cell_u);
			const MatrixX& K0 = Ke[mat_id];
			energy(cell) = (real)0.5*cell_u.dot(K0 * cell_u);
		}
	);
}

template class LinearFEMGrid<2>;
template class LinearFEMGrid<3>;
