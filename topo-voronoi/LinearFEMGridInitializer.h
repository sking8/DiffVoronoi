#pragma once
#include "LinearFEMGrid.h"

using namespace Meso;
template<int d> class LinearFEMGridInitializer {
	Typedef_VectorD(d);
public:
	void Apply(json& j, LinearFEMGrid<d>& simulator)
	{
		int scale = Json::Value(j, "scale", 32);
		real dx = 1.0 / scale;
		VectorDi grid_size = scale * VectorDi::Ones();
		Grid<d> grid(grid_size, dx, VectorD::Zero(), CORNER);

		//only one case now, cantilever beam
		BoundaryConditionGrid<d> bc;
		grid.Iterate_Nodes(
			[&](const VectorDi node) {
				if (node[0] == 0) { bc.Set_Psi_D(node, VectorD::Zero()); }
				else if (node[0] == grid.Counts()[0] - 1 && node[1] == 0) {
					bc.Set_Force(node, -VectorD::Unit(1));
					if constexpr (d == 3) { bc.forces[node] /= grid.Counts()[2]; }
				}
			}
		);

		Array<std::tuple<real, real>> materials;
		materials.push_back({ (real)1,(real)0.3 });
		Field<short, d> material_id;
		material_id.Init(grid.Cell_Grid(), 0); //use cell grid here
		simulator.Initialize(grid,bc,materials,material_id);
	}
};