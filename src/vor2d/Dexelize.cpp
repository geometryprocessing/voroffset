#include "Dexelize.h"
#define NANOSVG_IMPLEMENTATION
#include <nanosvg.h>
#include <geogram/basic/file_system.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
#include <geogram/basic/logger.h>
#include <geogram/basic/progress.h>
#include <geogram/basic/stopwatch.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_AABB.h>
#include <geogram/numerics/predicates.h>
#include <iostream>

#ifndef  DPI
#define DPI 90
#endif // ! DPI


voroffset::DoubleCompressedImage voroffset::create_dexels(const std::string &file)
{
	std::vector<Curve> contours;
	std::cerr << "[loading] " << file << " ... ";
	double tol = 1.5;
	NSVGimage* image = nsvgParseFromFile(file.c_str(), "mm", DPI);
	// For all curves in the file
	for (NSVGshape *shape = image->shapes; shape != NULL; shape = shape->next)
	{
		// Flatten the curve in a series of segments
		Curve poly;
		for (NSVGpath *path = shape->paths; path != NULL; path = path->next)
		{
			for (int i = 0; i < path->npts - 1; i++)
				poly.push_back(PointF(path->pts[2 * i], path->pts[2 * i + 1]));
		}
		contours.push_back(poly);
	}
	std::cerr << "Read a SVG of size : " << image->width << " x " << image->height
		<< " (" << contours.size() << " polygones)" << std::endl;
	DoubleCompressedImage dexels(std::ceil(image->height*25.4 / DPI), std::ceil(image->width*25.4 / DPI));
	dexels.fromImage(contours);
	nsvgDelete(image);
	return dexels;
}

void voroffset::dexel_dump(const std::string &filename, const DoubleCompressedImage &dexels)
{
	// Initialize the Geogram library
	GEO::initialize();

	// Import standard command line arguments, and custom ones
	GEO::CmdLine::import_arg_group("standard");

	GEO::Mesh mesh;

	for (int x = 0; x < dexels.height(); ++x)
	{
		for (int i = 0; 2 * i < dexels.at(x).size(); ++i)
		{
			GEO::vec3 xyz_min, xyz_max;
			xyz_min[0] = x;
			xyz_min[1] = dexels.at(x)[2 * i + 0];
			xyz_min[2] = 1;
			xyz_max[0] = x + 1;
			xyz_max[1] = dexels.at(x)[2 * i + 1];
			xyz_max[2] = 1;
			GEO::vec3 diff[4] =
			{
				GEO::vec3(0,0,0), GEO::vec3(1,0,0), GEO::vec3(0,1,0), GEO::vec3(1,1,0),
			};
			int v = mesh.vertices.nb();
			for (int lv = 0; lv < 4; ++lv)
			{
				for (int d = 0; d < 3; ++d)
				{
					diff[lv][d] = xyz_min[d] + diff[lv][d] * (xyz_max[d] - xyz_min[d]);
				}
				mesh.vertices.create_vertex(&diff[lv][0]);
			}
			mesh.facets.create_quad(v, v + 1, v + 2, v + 3);
		}
	}
	mesh.facets.compute_borders();
	mesh.facets.connect();
	mesh.vertices.remove_isolated();
	GEO::mesh_save(mesh, filename);
}
