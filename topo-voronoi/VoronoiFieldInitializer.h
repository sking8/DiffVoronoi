#pragma once
#include "VoronoiField.h"
#include "Random.h"

namespace Meso {
	template<int d> class VoronoiFieldInitializer {
		Typedef_VectorD(d);
	public:
		/*void Apply(json& j, VoronoiField<d>& voronoi_field)
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
		}*/

		void Apply(json& j, VoronoiField<d>& voronoi_field)
		{
			int scale = Json::Value(j, "scale", 32);
			int beta = Json::Value(j, "beta", 20);
			real alpha = Json::Value(j, "alpha", (real)50);
			int point_num = Json::Value(j, "point_num", 5); //point num per row
			real c = Json::Value(j, "c", (real)1);
			real dx = 1.0 / (real)scale;
			VectorDi grid_size = scale * VectorDi::Ones();
			Grid<d> grid(grid_size, dx, VectorD::Zero(), CENTER);
			Array<VectorD> points(Pow(point_num, d) + Pow(point_num + 1, d));
			real interval = 1.0 / (real) point_num;
			int a = Pow(point_num+1, d - 1);
			int b = Pow(point_num, d - 1);
			if constexpr (d == 2) {
				c = 1;
				for (int i = 0; i < points.size(); i++) {
					int y = i / (point_num + (point_num + 1));
					if ((i % (point_num + (point_num + 1))) < (point_num+1)) {
						int x = i % (point_num + (point_num + 1));
						points[i] = Vector2((real)x * interval, (real)y * interval);
					}
					else {
						int x = i % (point_num + (point_num + 1))-(point_num+1);
						points[i] = Vector2(((real)x+(real)0.5) * interval, ((real)y+(real)0.5) * interval);
					}
				}
			}
			else if constexpr (d == 3) {
				c = 2;
				for (int i = 0; i < points.size(); i++) {
					int y = i / (a + b);
					if ((i % (a + b)) < a) {
						int x = i % (a + b) % (point_num+1);
						int z = i % (a + b) / (point_num+1);
						points[i] = Vector3((real)x * interval, (real)y * interval, (real)z*interval);
					}
					else {
						int x = (i % (a + b)-a) % point_num;
						int z = (i % (a + b)-a) / point_num;
						points[i] = Vector3(((real)x + (real)0.5) * interval, ((real)y + (real)0.5) * interval, ((real)z + (real)0.5) * interval);
					}
				}
			}
			voronoi_field.Initialize(grid, points, point_num, beta, alpha, c);
		}
	};
}