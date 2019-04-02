////////////////////////////////////////////////////////////////////////////////
#include <CLI/CLI.hpp>
#include <YImage.hpp>
#include <vector>
#include <vor2d/CompressedImage.h>
#include <vor2d/Dexelize.h>
#include <vor2d/DoubleCompressedImage.h>
////////////////////////////////////////////////////////////////////////////////

#define WHITE YImage::YPixel({255, 255, 255, 255})
#define BLACK YImage::YPixel({0, 0, 0, 255})

bool pixel_is_white(YImage::YPixel p) {
	return p.r > 128 || p.g > 128 || p.b > 128 || p.a > 128;
}

int main(int argc, char *argv[]) {
	// Default arguments
	struct {
		std::string input;
		std::string output = "out.png";
		double radius = 0;
		bool erode = false;
		bool force = false;
		bool transpose = false;
		bool negate = false;
	} args;

	// Parse arguments
	CLI::App app("Offset2D");
	app.add_option("input,-i,--input", args.input, "Input image.")
			->required()
			->check(CLI::ExistingFile);
	app.add_option("output,-o,--output", args.output, "Output image.");
	app.add_option("-r,--radius", args.radius, "Dilation/erosion radius.");
	app.add_flag("-e,--erode", args.erode, "Erode instead of dilate.");
	app.add_flag("-f,--force", args.force, "Overwrite output file.");
	app.add_flag("-t,--transpose", args.transpose, "Transpose input image.");
	app.add_flag("-n,--negate", args.negate, "Negate input image.");
	try {
		app.parse(argc, argv);
	} catch (const CLI::ParseError &e) {
		return app.exit(e);
	}

	// Load SVG
	vor::DoubleCompressedImage dexels = vor::create_dexels(args.input.c_str());

	// Pre-processing operations
	if (args.transpose) {
		dexels.transposeInPlace();
	}
	if (args.negate) {
		dexels.negate();
	}

	// Offset
	if (args.radius > 0) {
		std::cout << "-- Performing offset by radius r = " << args.radius
							<< std::endl;
		if (args.erode) {
			dexels.erode(args.radius);
		} else {
			dexels.dilate(args.radius);
		}
	}

	// Save output image
	if (std::ifstream(args.output)) {
		// Output file exists!
		if (args.force) {
			std::cout << "-- Overwriting output file: " << args.output << std::endl;
			vor::dexel_dump(args.output.c_str(), dexels);
		} else {
			std::cerr << "-- Output file already exists. Please use -f to force overwriting."
				<< std::endl;
		}
	} else {
		std::cout << "-- Saving" << std::endl;
		vor::dexel_dump(args.output.c_str(), dexels);
	}
	return 0;
}
