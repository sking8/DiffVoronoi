#pragma once
#include "TopoOptVoronoi.h"

template<int d> class TopoOptVoronoiInitializer {
	Typedef_VectorD(d);
public:
	void Apply(json& j, TopoOptVoronoi<d>& optimizer)
	{
		std::string test = Meso::Json::Value(j, "test", std::string("cantilever_beam"));
		if (test == "cantilever_beam") Init_Cantilever_Beam(j, optimizer);
		else Assert(false, "test {} not exist", test);
	}

	void Init_Cantilever_Beam(json& j, TopoOptVoronoi<d>& optimizer) {
		int scale = Meso::Json::Value(j, "scale", 32);
		real strength = Meso::Json::Value(j, "strength", (real)1);
		real target_frac = Meso::Json::Value(j, "target_frac", (real)0.35);
		real mov_lim = Meso::Json::Value(j, "mov_lim", (real)0.05);
		int power = Meso::Json::Value(j, "power", 1);
		int point_num = Meso::Json::Value(j, "point_num", 20);
		int alpha = Meso::Json::Value(j, "alpha", 100);
		int beta = Meso::Json::Value(j, "beta", 50);
		int c = Meso::Json::Value(j, "c", 1);	//controls the dimension of the voronoi junction
		real dx = 1.0 / scale;
		VectorDi grid_size = scale * VectorDi::Ones();
		grid_size[1] /= 2;
		Meso::Grid<d> grid(grid_size, dx, VectorD::Zero(), Meso::CENTER);
		Meso::Info("Grid:{}", grid);
		Meso::Grid<d> corner_grid(grid_size + VectorDi::Ones(), dx, VectorD::Zero(), Meso::CORNER);
		Meso::Info("Corner Grid:{}", corner_grid);

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

		Meso::Array<VectorD> points(point_num);
		for (int i = 0; i < points.size(); i++) {
			points[i] = Meso::Random::Uniform_In_Box(grid.Domain_Min(Meso::CENTER), grid.Domain_Max(Meso::CENTER));
			Meso::Info("point {}: {}", i, points[i]);
		}

		Grid<d> spx_grid(grid_size, dx, VectorD::Zero());

		Array<std::tuple<real, real>> materials;
		materials.push_back({ (real)1,(real)0.3 });
		Field<short, d> material_id(spx_grid.cell_counts, 0);
		SoftBodyLinearFemGrid<d> linear_fem_grid;
		linear_fem_grid.Initialize(spx_grid, spx_bc, materials, material_id);
		Meso::VoronoiField<d> voronoi_field;
		voronoi_field.Initialize(grid, points, point_num, beta, alpha, c);
		optimizer.Init(linear_fem_grid, voronoi_field, target_frac, mov_lim, power);
	}
};