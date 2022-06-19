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
		VectorDi grid_size = (scale+1) * VectorDi::Ones(); //plus one for corner grid
		Grid<d> grid(grid_size, dx, VectorD::Zero(), CORNER);

		//only one case now, cantilever beam
		BoundaryConditionGrid<d> bc;
		grid.Iterate_Nodes(
			[&](const VectorDi node) {
				if (node[0] == 0) { bc.Set_Psi_D(node, VectorD::Zero()); }
				else if (node[0] == grid.Counts()[0] - 1 && node[1] == grid.Counts()[1] - 1) {
					bc.Set_Force(node, VectorD::Ones());
					if constexpr (d == 3) { bc.forces[node] /= grid.Counts()[2]; }
				}
			}
		);

		Array<std::tuple<real, real>> materials;
		Field<short, d> material_id;
		material_id.Init(grid.Cell_Grid(), 0); //use cell grid here
		materials.push_back({ (real)1,(real)0.3 });
		simulator.Initialize(grid,bc,materials,material_id);
	}
};