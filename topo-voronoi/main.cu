#include "Json.h"
#include "omp.h"
#include "TopoOptInitializer.h"
#include "TopologyOptimization.h"
#include "OptimizerDriver.h"

using namespace Meso;

template<int d>
void Run(json& j) {
	TopologyOptimization<d> optimizer;
	TopoOptInitializer<d> scene;
	OptimizerDriver driver;
	driver.Run(j, scene, optimizer);
}

int main(int argc, char** argv) {
	try {
		json j;

		if (argc > 1) {
			std::ifstream json_input(argv[1]);
			json_input >> j;
			json_input.close();
		}

		int thread_num = Json::Value(j, "thread_num", omp_get_max_threads());
		omp_set_num_threads(thread_num);
		Info("Using {} threads, out of {} available threads", omp_get_num_threads(), omp_get_max_threads());

		int dim = Json::Value(j, "dimension",2);
		if (dim == 2) {
			Run<2>(j);
		}
		else if (dim==3) {
			Run<3>(j);
		}
		else {
			Error("Dimension not supported!");
		}
	}
	catch (nlohmann::json::exception& e)
	{
		Meso::Info("json exception {}", e.what());
	}
	return 0;
}