//////////////////////////////////////////////////////////////////////////
// Auxillary Spray particles
// see: Real-Time Eulerian Water Simulation Using a Restricted Tall Cell Grid
// Copyright (c) (2018-), Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "Particles.h"
#include "LevelSet.h"
#include "RandomNumber.h"
#include "Advection.h"
#include "Json.h"

template<int d>
class MarkerPoints :public Points<d> {
public:
	Typedef_VectorDii(d);
	using Base = Points<d>;
	using Base::X;
public:
	MarkerPoints(const json& j = json::value_t::object) {
		New_Attributes();
	}
	Declare_Attribute(VectorD, V, v);//velocity
	Declare_Attribute(real, R, r);//radius
	Declare_Attribute(real, T, t);//time of creation

	Declare_Attribute_Inherent_Func(v, r, t);
};


template<int d>
class SpraySystem {
	Typedef_VectorDii(d);
public:
	MarkerPoints<d> trial_points;
	MarkerPoints<d> spray_points;
	MarkerPoints<d> foam_points;
	real phi_gen;//threshold of sampling trial points
	real spray_bandwidth;//if a trial particle excceeds this value, it becomes a spray particle
	int cell_sample_num;
	VectorD gravity_acc;
	real turn_foam_prob;
	real drag_coef;
	real trial_die_time;
	real foam_die_time;
	int closest_iter_num;

	std::function<bool(const VectorD)> extra_del_func = nullptr;

	//Run-time variables
	Field<int, d> cell_particle_counts;

	SpraySystem() {}
	SpraySystem(json& j) {
		phi_gen = Json::Value<real>(j, "phi_gen", -0.01);
		spray_bandwidth = Json::Value<real>(j, "spray_bandwidth", 0.01);
		cell_sample_num = Json::Value<int>(j, "cell_sample_num", 64);
		gravity_acc = Json::Value<VectorD>(j, "gravity_acc", -PhysicalConstants::g * VectorD::Unit(1));
		turn_foam_prob = Json::Value<real>(j, "turn_foam_prob", 0.5);
		drag_coef = Json::Value<real>(j, "drag_coef", 0.5);
		trial_die_time = Json::Value<real>(j, "trial_die_time", 1.0);
		foam_die_time = Json::Value<real>(j, "foam_die_time", 1.0);
		closest_iter_num = Json::Value<int>(j, "closest_iter_num", 5);
	}
	void Register_Cell_Counts(const Grid<d> &grid, Field<int, d>& cell_counts, const MarkerPoints<d>& particles) {
		//static int max_count; max_count = 0;
		for (int i = 0; i < particles.Size(); i++) {
			VectorDi cell = grid.Cell_Coord(particles.X(i));
			if (!grid.Valid_Cell(cell)) continue;
			cell_counts(cell)++;
			//max_count = std::max(max_count, cell_counts(cell));
		}
		//Info("max_count: {}", max_count);
	}
	void Sample_Trial_Points(const real time, const MacGrid<d>& mac_grid, const FaceField<real,d> &velocity, const LevelSet<d>& levelset) {
		Info("sample trial points with phi gen {} and {} points each cell", phi_gen, cell_sample_num);
		cell_particle_counts.Resize(mac_grid.grid.cell_counts, 0);
		Register_Cell_Counts(mac_grid.grid, cell_particle_counts, trial_points);
		Register_Cell_Counts(mac_grid.grid, cell_particle_counts, spray_points);
		Register_Cell_Counts(mac_grid.grid, cell_particle_counts, foam_points);

		Interpolation<d> intp(mac_grid);
		real dx = mac_grid.grid.dx;
		static Array<int> added_cell_indicies;added_cell_indicies.clear();
		int cell_num = mac_grid.grid.Number_Of_Cells();
		mac_grid.grid.Exec_Each(
			[&](const VectorDi& cell) {
				real cell_phi = levelset.phi(cell);
				if (phi_gen < cell_phi && cell_phi < 0) {
					int offset;
					int sample_num = std::max(0, cell_sample_num - cell_particle_counts(cell));
					if (sample_num == 0) return;
					for (int k = 0; k < sample_num; k++) {
						VectorD sample_point;
						int idx;
#pragma omp critical
						{
							idx = trial_points.Add_Element();//push_back() is better than resize()
							sample_point = mac_grid.grid.Center(cell) + RandomFunc::Random_Vector_Cartesian<d>(-0.5 * dx, 0.5 * dx);
						}
						trial_points.X(idx) = sample_point;
						trial_points.V(idx) = intp.Interpolate_Face_Vectors(velocity, sample_point);
						trial_points.R(idx) = -cell_phi;
						trial_points.T(idx) = time;
					}
				}
			}
		);
	}
	void Advect_Points(const real dt, const MacGrid<d> &mac_grid, const FaceField<real, d>& velocity, MarkerPoints<d> &points) {
		Interpolation<d> intp(mac_grid);
		points.Exec_Each(
			[&](const int i) {
				VectorD new_pos = Advection::Runge_Kutta2<d>(dt, intp, velocity, points.X(i));
				points.X(i) = new_pos;
			}
		);
	}
	void Move_Spray_Points(const real dt, const VectorD gravity_acc) {
		spray_points.Exec_Each(
			[&](const int i) {
				real v_norm = spray_points.V(i).norm();
				VectorD drag_acc = -v_norm * drag_coef * spray_points.V(i);//proportional to the square of rate
				spray_points.V(i) += dt * (gravity_acc + drag_acc);
				spray_points.X(i) += dt * spray_points.V(i);
			}
		);
	}
	void Trial_To_Spray(const real time, const MacGrid<d>& mac_grid, const LevelSet<d>& levelset) {
		static Array<int> trial_to_delete;
		trial_to_delete.resize(trial_points.Size());
		trial_points.Calc_Each(
			[&](const int i)->int {
				return levelset.Phi(trial_points.X(i)) > spray_bandwidth;
			},
			trial_to_delete
				);
		trial_points.Exec_Each(
			[&](const int i) {
				if (trial_to_delete[i]) {
					int idx;
#pragma omp critical
					idx = spray_points.Add_Element();
					spray_points.Copy_Element_From(idx, trial_points, i);
					spray_points.T(idx) = time;
				}
			}
		);
		trial_points.Delete_Elements(trial_to_delete);
	}
	void Spray_To_Foam(const real time, const MacGrid<d>& mac_grid, const LevelSet<d>& levelset, const real turn_probability) {
		static Array<int> spray_to_delete;
		spray_to_delete.resize(spray_points.Size());
		spray_points.Calc_Each(
			[&](const int i)->int {
				if (levelset.Phi(spray_points.X(i)) <= 0) {
					real p;
#pragma omp critical
					p = RandomFunc::Random_Real(0, 1);
					if (p <= turn_probability) return 1;
					else return 2;
				}
				return false;
			},
			spray_to_delete
				);
		spray_points.Exec_Each(
			[&](const int i) {
				if (spray_to_delete[i] == 1) {
					int idx;
#pragma omp critical
					idx = foam_points.Add_Element();
					foam_points.Copy_Element_From(idx, spray_points, i);
					VectorD pos = foam_points.X(idx);
					foam_points.X(idx) = levelset.Closest_Point_With_Iterations(pos, closest_iter_num);
					foam_points.T(idx) = time;
				}
			}
		);
		spray_points.Delete_Elements(spray_to_delete);
	}
	void Delete_Points(const Grid<d> &grid, const LevelSet<d> &levelset, const real earliest_time, const Box<d>& box, MarkerPoints<d>& points) {
		static Array<int> to_delete;
		to_delete.resize(points.Size()); std::fill(to_delete.begin(), to_delete.end(), 0);
		cell_particle_counts.Resize(grid.cell_counts, 0);
		for (int i = 0; i < points.Size(); i++) {
			VectorDi cell = grid.Cell_Coord(points.X(i));
			if (grid.Valid_Cell(cell)) {
				cell_particle_counts(cell)++;
				if (cell_particle_counts(cell) > cell_sample_num) to_delete[i] = 1;
			}
		}
		points.Exec_Each(
			[&](const int i) {
				if (extra_del_func && extra_del_func(points.X(i))) to_delete[i] = 1;
				else if (levelset.Phi(points.X(i)) <= phi_gen) to_delete[i] = 1;
				else if (!box.Inside(points.X(i)) || points.T(i) < earliest_time) to_delete[i] = 1;
			}
		);
		points.Delete_Elements(to_delete);
	}
	void Delete_Foam_Points(const Grid<d>& grid, const LevelSet<d>& levelset, const real earliest_time, const Box<d>& box, MarkerPoints<d>& points) {
		static Array<int> to_delete;
		to_delete.resize(points.Size()); std::fill(to_delete.begin(), to_delete.end(), 0);
		cell_particle_counts.Resize(grid.cell_counts, 0);
		for (int i = 0; i < points.Size(); i++) {
			VectorDi cell = grid.Cell_Coord(points.X(i));
			if (grid.Valid_Cell(cell)) {
				cell_particle_counts(cell)++;
				if (cell_particle_counts(cell) > cell_sample_num) to_delete[i] = 1;
			}
		}
		points.Exec_Each(
			[&](const int i) {
				if (extra_del_func && extra_del_func(points.X(i))) to_delete[i] = 1;
				else if (levelset.Phi(points.X(i)) <= phi_gen || levelset.Phi(points.X(i)) > spray_bandwidth) to_delete[i] = 1;
				else if (!box.Inside(points.X(i)) || points.T(i) < earliest_time) to_delete[i] = 1;
			}
		);
		points.Delete_Elements(to_delete);
	}
	void Advance(const real time, const real dt, const MacGrid<d> &mac_grid, const FaceField<real,d> &velocity, const LevelSet<d> &levelset) {
		////Step 1: Dynamics advance
		//1.1 advect trial points
		Advect_Points(dt, mac_grid, velocity, trial_points);
		//1.2 integrate spray points
		Move_Spray_Points(dt, gravity_acc);
		//1.3 advect foam points
		Advect_Points(dt, mac_grid, velocity, foam_points);
		foam_points.Exec_Each(
			[&](const int i) {
				foam_points.X(i) = levelset.Closest_Point_With_Iterations(foam_points.X(i), closest_iter_num);
			}
		);

		////Step 2: transition
		//2.1 trial => spray
		Trial_To_Spray(time, mac_grid, levelset);
		//2.2 spray => foam
		Spray_To_Foam(time, mac_grid, levelset, turn_foam_prob);
		 
		////Step 3: sample points
		Sample_Trial_Points(time, mac_grid, velocity, levelset);

		////Step 4: delete points outside the boundary
		Box<d> box(mac_grid.grid.domain_min, mac_grid.grid.domain_max);
		Delete_Points(mac_grid.grid, levelset, time - trial_die_time, box, trial_points);
		VectorD spray_max = mac_grid.grid.domain_max;
		spray_max[1] = 1e8;
		Delete_Points(mac_grid.grid, levelset, -1, Box<d>(mac_grid.grid.domain_min, spray_max), spray_points);//don't check time for spray points
		Delete_Foam_Points(mac_grid.grid, levelset, time - foam_die_time, box, foam_points);
		Info("Spray system with {} trial points, {} spray points and {} foam points", trial_points.Size(), spray_points.Size(), foam_points.Size());
	}
};
