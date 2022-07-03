#pragma once
#include "TopologyOptimization.h"

template<int d> class TopoOptInitializer {
	Typedef_VectorD(d);
public:
	void Apply(json& j, TopologyOptimization<d>& optimizer)
	{
		int scale = Meso::Json::Value(j, "scale", 32);
		real strength = Meso::Json::Value(j, "strength", (real)1);
		real target_frac = Meso::Json::Value(j, "target_frac", (real)0.35);
		real mov_lim = Meso::Json::Value(j, "mov_lim", (real)0.05);
		int power = Meso::Json::Value(j, "power", 3);
		real dx = 1.0 / scale;
		VectorDi grid_size = scale * VectorDi::Ones();
		grid_size[1] /= 2;
		Meso::Grid<d> grid(grid_size, dx, VectorD::Zero(), Meso::CENTER);
		Meso::Info("Grid:{}", grid);
		Meso::Grid<d> corner_grid(grid_size + VectorDi::Ones(), dx, VectorD::Zero(), Meso::CORNER);
		Meso::Info("Corner Grid:{}", corner_grid);
		//only one case now, cantilever beam
		BoundaryConditionGrid<d> spx_bc;
		corner_grid.Iterate_Nodes(
			[&](const VectorDi node) {
				if (node[0] == 0) { spx_bc.Set_Psi_D(node, VectorD::Zero()); }
				else if (node[0] == corner_grid.Counts()[0] - 1 && node[1] == 0) {
					spx_bc.Set_Force(node, -VectorD::Unit(1) * strength);
					if constexpr (d == 3) { spx_bc.forces[node] /= corner_grid.Counts()[2]; }
				}
			}
		);

		Grid<d> spx_grid(grid_size, dx, VectorD::Zero());

		Array<std::tuple<real, real>> materials;
		materials.push_back({ (real)1,(real)0.3 });
		Field<short, d> material_id(spx_grid.cell_counts, 0);
		SoftBodyLinearFemGrid<d> linear_fem_grid;
		linear_fem_grid.Initialize(spx_grid, spx_bc, materials, material_id);
		optimizer.Init(linear_fem_grid, grid, target_frac, mov_lim, power);
	}
};