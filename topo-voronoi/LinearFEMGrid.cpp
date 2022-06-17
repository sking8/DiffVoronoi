#include "LinearFEMGrid.h"
#include "ColorGrid.h"
using namespace Meso;

//=================================Helper function================================
template<int d>
Vector<int,d> Corner_Offset(const Vector<int, d>& center, int i) {
	Assert(i < pow(2, d), "Corner_Offset: index out of  range");
	if constexpr (d == 2) { return center + Vector2i(i & 0x1, i >> 1 & 0x1); }
	else if constexpr (d == 3) { return center + Vector3i(i & 0x1, i >> 1 & 0x1, i >> 2 & 0x1); }
	else { Error("Corner_Offset: dimension not supported"); return Vector<int, d>(); }
}

//================================Main body=======================================
template<int d> void LinearFemGrid<d>::Initialize(const Grid<d> _node_grid)
{
	colored_cell_ptr.resize(Pow(2, d) + 1);
	for (int i = 0; i < colored_cell_ptr.size(); i++)colored_cell_ptr[i] = i;

	node_grid = _node_grid;
	Add_Material((real)1, (real).3);									//Should it be here?
	material_id.Init(_node_grid.Cell_Grid(), 0);						// be careful about the difference between cell counts and node counts

	int n = node_grid.Counts().prod() * d;								// be careful about the difference between cell counts and node counts
	K.resize(n, n); u.resize(n); u.fill((real)0); f.resize(n); f.fill((real)0);
}

template<int d> void LinearFemGrid<d>::Add_Material(real youngs, real poisson)
{
	MatrixX Ke0; LinearFEMFunc::Cell_Stiffness_Matrix<d>(youngs, poisson, node_grid.dx, Ke0); Ke.push_back(Ke0);
}

template<int d> void LinearFemGrid<d>::Allocate_K()
{
	std::vector<Triplet<real>> elements; //Fan: Can I use Array here?
	//Fan: to change
	//iterate_node(iter, node_grid) {
	//	const VectorDi& node = iter.Coord(); int r = node_grid.Node_Index(node);
	//	for (int i = 0; i < Grid<d>::Number_Of_Nb_R(); i++) {
	//		VectorDi nb = Grid<d>::Nb_R(node, i);
	//		if (node_grid.Valid_Node(nb)) {
	//			int c = Node_Index_In_K(nb);
	//			for (int rr = r * d; rr < (r + 1) * d; rr++)for (int cc = c * d; cc < (c + 1) * d; cc++) { elements.push_back(Triplet<real>(rr, cc, (real)0)); }
	//		}
	//	}
	//}
	
	//need to check which side that data resides
	node_grid.Iterate_Nodes(
		[&](const VectorDi node) {
		int r = node_grid.Index(node);
		for (int nb_r = 0; nb_r < node_grid.Neighbor_Node_Number(); nb_r++) {
			int c = node_grid.Index(node_grid.Neighbor_Ring_Node(node, nb_r));
			if (!node_grid.Valid(c)) { continue; }
			for (int rr = r * d; rr < (r + 1) * d; rr++)for (int cc = c * d; cc < (c + 1) * d; cc++) { elements.push_back(Triplet<real>(rr, cc, (real)0)); }
		}
		}
	);

	K.setFromTriplets(elements.begin(), elements.end());
	K.makeCompressed();

	ColorGrid::Color<d>(node_grid.Counts(), colored_cell_ptr, colored_cell_indices);
}

template<int d> void LinearFemGrid<d>::Update_K_And_f()
{
	////Update K
	int color_n = colored_cell_ptr.size()-1;
	for (int c = 0; c < color_n; c++) {
#pragma omp parallel for
		for (int i = colored_cell_ptr[c]; i < colored_cell_ptr[c + 1]; i++) {
			VectorDi cell = node_grid.Cell_Grid().Coord(colored_cell_indices[i]); int mat_id = material_id(cell);
			Array<int> cell_node_indices(pow(2, d), 0);
			for (int j = 0; j < cell_node_indices.size(); j++) {
				cell_node_indices[j] = node_grid.Index(Corner_Offset<d>(cell, j));
			}
			LinearFEMFunc::Add_Cell_Stiffness_Matrix<d>(K, Ke[mat_id], cell_node_indices);
		}
	}

	////Update rhs
	f.fill((real)0);
	for (auto& b : bc.forces) {
		VectorDi node = b.first; VectorD force = b.second;
		for (int axis = 0; axis < d; axis++) { int idx = node_grid.Index(node) * d + axis; f[idx] += force[axis]; }
	}

	////Update bc
	for (auto& b : bc.psi_D_values) {
		VectorDi node = b.first; VectorD dis = b.second;
		for (int axis = 0; axis < d; axis++) {
			int idx = node_grid.Index(node) * d + axis;
			LinearFEMFunc::Set_Dirichlet_Boundary_Helper(K, f, idx, dis[axis]);
		}
	}
}

template<int d> void LinearFemGrid<d>::Solve()
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
//		//GeometricMultiGrid::GMGPCG_GPU<d>(K,u,f,node_grid.node_counts,multigrid_params);
//		gmg_solver_gpu.update_A_levels = true;
//
//		gmg_solver_gpu.Initialize(K, node_grid.node_counts, multigrid_params, &material_id);
//		gmg_solver_gpu.Solve(u, f);
//	}
//	else {
//		gmg_solver_cpu.update_A_levels = true;
//		gmg_solver_cpu.Initialize(K, node_grid.node_counts, multigrid_params, &material_id);
//		gmg_solver_cpu.Solve(u, f);
//	}
//#else
//	gmg_solver_cpu.Initialize(K, node_grid.node_counts, multigrid_params, &material_id);
//	gmg_solver_cpu.Solve(u, f);
//#endif

	//Fan: to use the amgcl solver
}

template<int d> void LinearFemGrid<d>::Set_Fixed(const VectorDi& node)
{
	bc.psi_D_values[node] = VectorD::Zero();
}

template<int d> void LinearFemGrid<d>::Set_Displacement(const VectorDi& node, const VectorD& dis)
{
	bc.psi_D_values[node] = dis;
}

template<int d> void LinearFemGrid<d>::Set_Force(const VectorDi& node, const VectorD& force)
{
	bc.forces[node] = force;
}

template<int d> void LinearFemGrid<d>::Compute_Cell_Displacement(const VectorX& u, const VectorDi& cell, VectorX& cell_u) const
{
	int number_of_cell_nodes = pow(2, d);
	cell_u.resize(number_of_cell_nodes * d);
	for (int i = 0; i < number_of_cell_nodes; i++) {
		VectorDi nb_node = Corner_Offset<d>(cell, i); int nb_node_mtx_idx = node_grid.Index(nb_node);
		for (int j = 0; j < d; j++)cell_u(i * d + j) = u(nb_node_mtx_idx * d + j);
	}
}

template<int d> void LinearFemGrid<d>::Compute_Elastic_Energy(Field<real, d>& energy) const
{
	energy.Init(node_grid.Cell_Grid().Counts(), (real)0);
//	int cell_n = node_grid.Number_Of_Cells();
//#pragma omp parallel for
//	for (int i = 0; i < cell_n; i++) {
//		VectorDi cell = node_grid.Cell_Coord(i);
//		int mat_id = material_id(cell); if (mat_id == -1)continue;
//		VectorX cell_u; Compute_Cell_Displacement(u, cell, cell_u);
//		const MatrixX& K0 = Ke[mat_id];
//		energy(cell) = cell_u.dot(K0 * cell_u);
//		if (variable_coef) { energy(cell) *= (*variable_coef)(cell); }
//	}	////multiplying rho

	//Fan: be careful about the cell_grid size here
	Grid<d> cell_grid = node_grid.Cell_Grid();
	cell_grid.Exec_Nodes(
		[&](const VectorDi cell) {
			int mat_id = material_id(cell); if (mat_id == -1)return;
			VectorX cell_u; Compute_Cell_Displacement(u, cell, cell_u);
			const MatrixX& K0 = Ke[mat_id];
			energy(cell) = (real)0.5*cell_u.dot(K0 * cell_u);
		}
	);
}

template class LinearFemGrid<2>;
template class LinearFemGrid<3>;
