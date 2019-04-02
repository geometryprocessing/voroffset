////////////////////////////////////////////////////////////////////////////////
#include <vor3d/CompressedVolume.h>
#include <vor3d/VoronoiVorPower.h>
#include <vor3d/VoronoiBruteForce.h>
#include <vor3d/Dexelize.h>
#include <vor3d/Timer.h>
#ifdef USE_TBB
#include <tbb/task_scheduler_init.h>
#endif
#include <CLI/CLI.hpp>
#include <json.hpp>
#include <geogram/basic/logger.h>
#include <geogram/basic/stopwatch.h>
#include <fstream>
#include <iostream>
#include <thread>
#include <vector>
////////////////////////////////////////////////////////////////////////////////

int main(int argc, char * argv[]) {
	// Default arguments
	struct {
		std::string input;
		std::string output_mesh = "output.obj";
		std::string output_json = "";
		std::string method = "ours";
		std::string operation = "dilation";
		double radius = 8;
		double dexels_size = 1;
		int padding = 0;
		int num_dexels = 256;
		unsigned int num_thread = std::max(1u, std::thread::hardware_concurrency());
		bool force = false;
		bool radius_in_mm = false;
	} args;

	// Parse arguments
	CLI::App app("Offset3D");
	app.add_option("input,-i,--input", args.input, "Input model")->required()->check(CLI::ExistingFile);
	app.add_option("output,-o,--output", args.output_mesh, "Output model", true);
	app.add_option("-j,--json", args.output_json, "Output json file");
	app.add_option("-d,--dexels_size", args.dexels_size, "Size of a dexel (in mm)", true);
	app.add_option("-n,--num_dexels", args.num_dexels, "Number of dexels (-1 to use dexel size instead)", true);
	app.add_option("-p,--padding", args.padding, "Padding (in #dexels)");
	app.add_option("-t,--num_thread", args.num_thread, "Number of threads", true);
	app.add_option("-r,--radius", args.radius, "Dilation/erosion radius (in #dexels)", true);
	app.add_set("-m,--method", args.method, {"ours","brute_force"}, "The method to use");
	app.add_set("-x,--apply", args.operation, {"noop","dilation","erosion","closing","opening"},
		"Morphological operation to apply", true);
	app.add_flag("-f,--force", args.force, "Overwrite output file");
	app.add_flag("-u,--radius_in_mm", args.radius_in_mm, "Radius is given in mm instead");
	try {
		app.parse(argc, argv);
	} catch (const CLI::ParseError &e) {
		return app.exit(e);
	}

	// Declare stuff
	vor3d::CompressedVolume input;
	vor3d::CompressedVolume output;
	double time = 0;
	double time_1, time_2;
	nlohmann::json json_data;

	// Load input model and dexelize
	Timer t;
	{
		GEO::Stopwatch W("Dexelize");
		input = vor3d::create_dexels(args.input, args.dexels_size, args.padding, args.num_dexels);
	}

	// Compute radius in #dexels if needed, and display some stats
	if (args.radius_in_mm) {
		args.radius /= input.spacing();
	}
	GEO::Logger::out("Stats") << "Grid size: " << input.gridSize().transpose() << std::endl;
	GEO::Logger::out("Stats") << "Spacing (in mm): " << input.spacing() << std::endl;
	GEO::Logger::out("Stats") << "Origin (in mm): " << input.origin().transpose() << std::endl;
	GEO::Logger::out("Stats") << "Extent (in mm): " << input.extent().transpose() << std::endl;
	GEO::Logger::out("Stats") << "Radius (in #dexels): " << args.radius << std::endl;
	GEO::Logger::out("Stats") << "Number of threads: " << args.num_thread << std::endl;

	if (!args.output_json.empty()) {
		json_data = {
			{ "method", args.method },
			{ "num_threads", args.num_thread },
			{ "model_name", args.input },
			{ "voxel_size", input.spacing() },
			{ "padding", input.padding() },
			{ "num_dexels", args.num_dexels },
			{ "radius", args.radius },
			{ "grid_size", { input.gridSize()(0),input.gridSize()(1) } },
			{ "num_segments", input.numSegments() },
			{ "operation", args.operation }
		};
	}

	#ifdef USE_TBB
	tbb::task_scheduler_init init(args.num_thread);
	#endif

	// Create offset operator
	GEO::Logger::div("Offseting");
	std::unique_ptr<vor3d::VoronoiMorpho> op;
	if (args.method == "ours") {
		op = std::make_unique<vor3d::VoronoiMorphoVorPower>();
	} else if (args.method == "brute_force") {
		op = std::make_unique<vor3d::VoronoiMorphoBruteForce>();
	} else {
		GEO::Logger::err("Offset") << "Invalid method: " << args.method << std::endl;
		return 1;
	}
	assert(op);

	// Apply operation
	if (args.operation == "noop") {
		output = input;
	} else if (args.operation == "erosion") {
		GEO::Stopwatch W("Erosion");
		op->erosion(input, output, args.radius, time_1, time_2);
	} else if (args.operation == "dilation") {
		GEO::Stopwatch W("Dilation");
		op->dilation(input, output, args.radius, time_1, time_2);
	} else if (args.operation == "closing") {
		GEO::Stopwatch W("Closing");
		vor3d::CompressedVolume tmp;
		op->dilation(input, tmp, args.radius, time_1, time_2);
		op->erosion(tmp, output, args.radius, time_1, time_2);
	} else if (args.operation == "opening") {
		GEO::Stopwatch W("Opening");
		vor3d::CompressedVolume tmp;
		op->erosion(input, tmp, args.radius, time_1, time_2);
		op->dilation(tmp, output, args.radius, time_1, time_2);
	} else {
		throw std::invalid_argument("Operation");
	}

	GEO::Logger::div("Saving");
	if (!args.output_mesh.empty()) {
		if (std::ifstream(args.output_mesh)) {
			// Output file exists!
			if (args.force) {
				GEO::Logger::out("Save") << "Overwriting output file: " << args.output_mesh << std::endl;
				GEO::Stopwatch W("Save");
				vor3d::dexel_dump(args.output_mesh, output);
			} else {
				GEO::Logger::out("Save") << "Output mesh already exists. Please use -f to force overwriting." << std::endl;
			}
		} else {
			GEO::Stopwatch W("Save");
			std::ofstream out(args.output_mesh.c_str());
			vor3d::dexel_dump(args.output_mesh, output);
		}
	}

	time = t.get();
	if (!args.output_json.empty()) {
		json_data["time"] = time;
		json_data["time_first_pass"] = time_1;
		json_data["time_second_pass"] = time_2;
	}

	if (!args.output_json.empty()) {
		if (std::ifstream(args.output_json)) {
			// Output file exists!
			if (args.force) {
				GEO::Logger::out("Save") << "Overwriting output file: " << args.output_json << std::endl;
				std::ofstream o(args.output_json);
				o << std::setw(4) << json_data << std::endl;
			} else {
				GEO::Logger::out("Save") << "Output json already exists. Please use -f to force overwriting." << std::endl;
			}
		} else {
			std::ofstream o(args.output_json);
			o << std::setw(4) << json_data << std::endl;
		}
	}
	return 0;
}
