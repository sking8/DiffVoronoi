#pragma once
#include "VoronoiField.h"
#include "Random.h"

using namespace Meso;
template<int d> class VoronoiFieldInitializer {
	Typedef_VectorD(d);
public:
	void Apply(json& j, VoronoiField<d>& voronoi_field)
	{
		int scale = Json::Value(j, "scale", 32);
		int point_num = Json::Value(j, "point_num", 5);
		real dx = 1.0 / scale;
		VectorDi grid_size = scale * VectorDi::Ones();
		Grid<d> grid(grid_size, dx, VectorD::Zero(), COLLOC);
		Array<VectorD> points (point_num);
		for (int i = 0; i < points.size(); i++) {
			points[i] = Random::Uniform_In_Box(grid.Domain_Min(),grid.Domain_Max());
		}
		voronoi_field.Initialize(grid, points, point_num);
	}
};