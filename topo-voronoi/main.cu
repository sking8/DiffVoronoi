#include "Json.h"
#include "omp.h"
#include "TopoOptInitializer.h"
#include "TopologyOptimization.h"
#include "Driver.h"
#include "OptimizerDriver.h"
#include "VoronoiField.h"
#include "VoronoiFieldInitializer.h"
//#include "LinearFEMGrid.h"
//#include "LinearFEMGridInitializer.h"

template<int d>
void Run_Topology_Optimization(json& j) {
	TopologyOptimization<d> optimizer;
	TopoOptInitializer<d> scene;
	Meso::OptimizerDriver driver;
	driver.Run(j, scene, optimizer);
}

template<int d>
void Run_Voronoi_Field(json& j) {
	Meso::VoronoiField<d> voronoi_field;
	Meso::VoronoiFieldInitializer<d> scene;
	Meso::Driver driver;
	driver.Run(j, scene, voronoi_field);
}

//template<int d>
//void Run_Linear_FEM_Grid(json& j) {
//	LinearFEMGrid<d> voronoi_field;
//	LinearFEMGridInitializer<d> scene;
//	Driver driver;
//	driver.Run(j, scene, voronoi_field);
//}

int main(int argc, char** argv) {
	try {
		json j;

		if (argc > 1) {
			std::ifstream json_input(argv[1]);
			json_input >> j;
			json_input.close();
		}

		int thread_num = Meso::Json::Value(j, "thread_num", omp_get_max_threads());
		omp_set_num_threads(thread_num);
		Meso::Info("Using {} threads, out of {} available threads", omp_get_num_threads(), omp_get_max_threads());

		int dim = Meso::Json::Value(j, "dimension", 2);
		std::string app = Meso::Json::Value(j, "app", std::string("topology_optimization"));

		if (app == "topology_optimization") {
			if (dim == 2) { Run_Topology_Optimization<2>(j); }
			else if (dim == 3) { Run_Topology_Optimization<3>(j); }
			else { Meso::Error("Dimension not supported!"); }
		}
		else if (app == "voronoi_field") {
			if (dim == 2) { Run_Voronoi_Field<2>(j); }
			else if (dim == 3) { Run_Voronoi_Field<3>(j); }
			else { Meso::Error("Dimension not supported!"); }
		}
		/*else if (app == "linear_fem_grid") {
			if (dim == 2) { Run_Linear_FEM_Grid<2>(j); }
			else if (dim == 3) { Run_Linear_FEM_Grid<3>(j); }
			else { Error("Dimension not supported!"); }
		}*/
		else {
			Meso::Error("main: Invalid app");
		}
	}
	catch (nlohmann::json::exception& e)
	{
		Meso::Info("json exception {}", e.what());
	}
	return 0;
}