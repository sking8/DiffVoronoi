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
		int beta = Json::Value(j, "beta", 20);
		real alpha = Json::Value(j, "alpha", (real)50);
		real c = Json::Value(j, "c", (real)1);
		real dx = 1.0 / (real)scale;
		VectorDi grid_size = scale * VectorDi::Ones();
		Grid<d> grid(grid_size, dx, VectorD::Zero(), CENTER);
		Array<VectorD> points (point_num);
		Info("domain min: {}, domain max: {}",grid.Domain_Min(CENTER), grid.Domain_Max(CENTER));
		for (int i = 0; i < points.size(); i++) {
			points[i] = Random::Uniform_In_Box(grid.Domain_Min(CENTER),grid.Domain_Max(CENTER));
			Info("point {}: {}", i, points[i]);
		}
		voronoi_field.Initialize(grid, points, point_num, beta, alpha, c);
	}
};