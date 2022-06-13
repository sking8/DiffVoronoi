#pragma once
#include "TopologyOptimization.h"

using namespace Meso;
template<int d> class TopoOptInitializer {
	Typedef_VectorD(d);
public:
	void Apply(json& j, TopologyOptimization<d>& optimizer)
	{
		int scale = Json::Value(j, "scale", 32);
		real dx = 1.0 / scale;
		VectorDi grid_size = scale * VectorDi::Ones();
		Grid<d> grid(grid_size, dx, VectorD::Zero(), COLLOC);
		optimizer.Init(grid);
	}
};