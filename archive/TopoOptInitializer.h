#pragma once
#include "TopologyOptimization.h"

using namespace Meso;
template<int d> class TopoOptInitializer {
	Typedef_VectorD(d);
public:
	void Apply(json& j, TopologyOptimization<d>& optimizer)
	{
		int scale = Json::Value(j, "scale", 32);
		real strength = Json::Value(j, "strength", (real)1);
		real target_frac = Json::Value(j, "target_frac", (real)0.35);
		real mov_lim = Json::Value(j, "mov_lim", (real)0.05);
		int power = Json::Value(j, "power", 3);
		real dx = 1.0 / scale;
		VectorDi grid_size = scale * VectorDi::Ones();
		grid_size[1] /= 2;
		Grid<d> grid(grid_size, dx, VectorD::Zero(), CENTER);
		Info("Grid:{}", grid);
		Grid<d> corner_grid(grid_size + VectorDi::Ones(), dx, VectorD::Zero(), CORNER);
		Info("Corner Grid:{}", corner_grid);
		//only one case now, cantilever beam
		BoundaryConditionGrid<d> bc;
		corner_grid.Iterate_Nodes(
			[&](const VectorDi node) {
				if (node[0] == 0) { bc.Set_Psi_D(node, VectorD::Zero()); }
				else if (node[0] == corner_grid.Counts()[0] - 1 && node[1] == 0) {
					bc.Set_Force(node, -VectorD::Unit(1) * strength);
					if constexpr (d == 3) { bc.forces[node] /= corner_grid.Counts()[2]; }
				}
			}
		);

		Array<std::tuple<real, real>> materials;
		materials.push_back({ (real)1,(real)0.3 });
		Field<short, d> material_id;
		material_id.Init(grid, 0); //use cell grid here
		LinearFEMGrid<d> linear_fem_grid;
		linear_fem_grid.Initialize(corner_grid, bc, materials, material_id);
		optimizer.Init(linear_fem_grid, target_frac, mov_lim, power);
	}
};