#include "vor3d/Dexelize.h"
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
#include <random>
#include <chrono>
#include <algorithm>
#include <array>
#include <iterator>
//#include "tbb/blocked_range.h"
//#include "tbb/parallel_for.h"
//#include "tbb/task_scheduler_init.h"

////////////////////////////////////////////////////////////////////////////////
#ifndef TEST_WID
#define TEST_WID 25
#endif
#ifndef TEST_LEN
#define TEST_LEN 30
#endif
#ifndef TEST_HEI
#define TEST_HEI 40
#endif

namespace
{

	template <typename Scalar, size_t Rows>
	inline std::ostream& operator<<(std::ostream &out, std::array<Scalar, Rows> v)
	{
		out << "{";
		if (!v.empty())
		{
			std::copy(v.begin(), v.end() - 1, std::ostream_iterator<Scalar>(out, "; "));
			out << v.back();
		}
		out << "}";
		return out;
	}

	////////////////////////////////////////////////////////////////////////////////
	// NOTE: Function `point_in_triangle_2d` comes from SDFGen by Christopher Batty.
	// https://github.com/christopherbatty/SDFGen/blob/master/makelevelset3.cpp
	////////////////////////////////////////////////////////////////////////////////

	// calculate twice signed area of triangle (0,0)-(x1,y1)-(x2,y2)
	// return an SOS-determined sign (-1, +1, or 0 only if it's a truly degenerate triangle)
	int orientation(
		double x1, double y1, double x2, double y2, double &twice_signed_area)
	{
		twice_signed_area = y1 * x2 - x1 * y2;
		if (twice_signed_area > 0) return 1;
		else if (twice_signed_area < 0) return -1;
		else if (y2 > y1) return 1;
		else if (y2 < y1) return -1;
		else if (x1 > x2) return 1;
		else if (x1 < x2) return -1;
		else return 0; // only true when x1==x2 and y1==y2
	}

	// -----------------------------------------------------------------------------

	// robust test of (x0,y0) in the triangle (x1,y1)-(x2,y2)-(x3,y3)
	// if true is returned, the barycentric coordinates are set in a,b,c.
	bool point_in_triangle_2d(
		double x0, double y0, double x1, double y1,
		double x2, double y2, double x3, double y3,
		double &a, double &b, double &c)
	{
		x1 -= x0; x2 -= x0; x3 -= x0;
		y1 -= y0; y2 -= y0; y3 -= y0;
		int signa = orientation(x2, y2, x3, y3, a);
		if (signa == 0) return false;
		int signb = orientation(x3, y3, x1, y1, b);
		if (signb != signa) return false;
		int signc = orientation(x1, y1, x2, y2, c);
		if (signc != signa) return false;
		double sum = a + b + c;
		geo_assert(sum != 0); // if the SOS signs match and are nonzero, there's no way all of a, b, and c are zero.
		a /= sum;
		b /= sum;
		c /= sum;
		return true;
	}

	// -----------------------------------------------------------------------------

	// \brief Computes the (approximate) orientation predicate in 2d.
	// \details Computes the sign of the (approximate) signed volume of
	//  the triangle p0, p1, p2
	// \param[in] p0 first vertex of the triangle
	// \param[in] p1 second vertex of the triangle
	// \param[in] p2 third vertex of the triangle
	// \retval POSITIVE if the triangle is oriented positively
	// \retval ZERO if the triangle is flat
	// \retval NEGATIVE if the triangle is oriented negatively
	// \todo check whether orientation is inverted as compared to
	//   Shewchuk's version.
	inline GEO::Sign orient_2d_inexact(GEO::vec2 p0, GEO::vec2 p1, GEO::vec2 p2)
	{
		double a11 = p1[0] - p0[0];
		double a12 = p1[1] - p0[1];

		double a21 = p2[0] - p0[0];
		double a22 = p2[1] - p0[1];

		double Delta = GEO::det2x2(
			a11, a12,
			a21, a22
		);

		return GEO::geo_sgn(Delta);
	}

	////////////////////////////////////////////////////////////////////////////////


	/**
	 * @brief      { Intersect a vertical ray with a triangle }
	 *
	 * @param[in]  M     { Mesh containing the triangle to intersect }
	 * @param[in]  f     { Index of the facet to intersect }
	 * @param[in]  q     { Query point (only XY coordinates are used) }
	 * @param[out] z     { Intersection }
	 *
	 * @return     { description_of_the_return_value }
	 */
	template<int X = 0, int Y = 1, int Z = 2>
	int intersect_ray_z(const GEO::Mesh &M, GEO::index_t f, const GEO::vec3 &q, double &z)
	{
		using namespace GEO;

		index_t c = M.facets.corners_begin(f);
		const vec3& p1 = Geom::mesh_vertex(M, M.facet_corners.vertex(c++));
		const vec3& p2 = Geom::mesh_vertex(M, M.facet_corners.vertex(c++));
		const vec3& p3 = Geom::mesh_vertex(M, M.facet_corners.vertex(c));

		double u, v, w;
		if (point_in_triangle_2d(
			q[X], q[Y], p1[X], p1[Y], p2[X], p2[Y], p3[X], p3[Y], u, v, w))
		{
			z = u * p1[Z] + v * p2[Z] + w * p3[Z];
			auto sign = orient_2d_inexact(vec2(p1[X], p1[Y]), vec2(p2[X], p2[Y]), vec2(p3[X], p3[Y]));
			switch (sign)
			{
			case GEO::POSITIVE: return 1;
			case GEO::NEGATIVE: return -1;
			case GEO::ZERO:
			default: return 0;
			}
		}

		return 0;
	}

	// -----------------------------------------------------------------------------

	void compute_sign(const GEO::Mesh &M,
		const GEO::MeshFacetsAABB &aabb_tree, vor3d::CompressedVolume &dexels)
	{
		auto size = dexels.gridSize();

		try
		{
			GEO::ProgressTask task("Ray marching", 100);

			GEO::vec3 min_corner, max_corner;
			GEO::get_bbox(M, &min_corner[0], &max_corner[0]);

			const GEO::vec3 origin(dexels.origin().data());
			const double origin_z = origin[2];
			const double spacing = dexels.spacing();
			//GEO::parallel_for([&](int y)
			for (int y = 0; y < size[1]; ++y)
			{
				//if (GEO::Thread::current()->id() == 0) {
				// task.progress((int) (100.0 * y / size[1] * GEO::Process::number_of_cores()));
				//}
				for (int x = 0; x < size[0]; ++x)
				{
					GEO::vec2 center_xy(dexels.dexelCenter(x, y).data());
					GEO::vec3 center(center_xy[0], center_xy[1], 0);

					GEO::Box box;
					box.xyz_min[0] = box.xyz_max[0] = center[0];
					box.xyz_min[1] = box.xyz_max[1] = center[1];
					box.xyz_min[2] = min_corner[2] - spacing;
					box.xyz_max[2] = max_corner[2] + spacing;

					std::vector<double> inter;
					auto action = [&M, &inter, &center, &spacing](GEO::index_t f)
					{
						double z;
						if (intersect_ray_z(M, f, center, z))
						{
							inter.push_back(z / spacing);
						}
					};
					aabb_tree.compute_bbox_facet_bbox_intersections(box, action);
					std::sort(inter.begin(), inter.end());

					dexels.at(x, y).resize(inter.size());
					//dexels.at(x, y).resize(inter.size() / 2);
					size_t n = dexels.at(x, y).size();
					std::copy_n(inter.begin(), inter.size(), dexels.at(x, y).begin());
					/*for (int i = 0; i < n; i++)
					{
					dexels.at(x, y)[i] = vor3d::m_Segment(inter[2 * i], inter[2 * i + 1], radius);
					}*/
				}
			}//, 0, size[1]);
		}
		catch (const GEO::TaskCanceled&)
		{
			// Do early cleanup
		}
	}

}

////////////////////////////////////////////////////////////////////////////////

voroffset3d::CompressedVolume voroffset3d::create_dexels(
	const std::string &filename, double &voxel_size, int padding, int num_voxels)
{
	// Initialize the Geogram library
	GEO::initialize();

	// Import standard command line arguments, and custom ones
	GEO::CmdLine::import_arg_group("standard");

	// Declare a mesh
	GEO::Mesh M;

	// Load the mesh and display timings
	GEO::Logger::div("Loading");
	{
		GEO::Stopwatch W("Load");
		if (!GEO::mesh_load(filename, M))
		{
			throw std::runtime_error("Invalid input mesh.");
		}
		geo_assert(M.vertices.dimension() == 3);
	}

	// Initialize voxel grid and AABB tree
	GEO::vec3 min_corner, max_corner;
	GEO::get_bbox(M, &min_corner[0], &max_corner[0]);
	GEO::vec3 extent = max_corner - min_corner;
	if (num_voxels > 0)
	{
		// Force number of voxels along longest axis
		double max_extent = std::max(extent[0], std::max(extent[1], extent[2]));
		voxel_size = max_extent / num_voxels;
	}
	GEO::MeshFacetsAABB aabb_tree(M);

	// Dexelize the input mesh
	GEO::Logger::div("Dexelizing");
	CompressedVolume dexels(
		Eigen::Vector3d(min_corner[0], min_corner[1], min_corner[2]),
		Eigen::Vector3d(extent[0], extent[1], extent[2]),
		voxel_size, padding);
	compute_sign(M, aabb_tree, dexels);
	return dexels;
}

////////////////////////////////////////////////////////////////////////////////
// NOTE: Function `dexel_dump` comes from SDFGen by Christopher Batty.
// https://github.com/jdumas/geotools/blob/master/voxmesh
////////////////////////////////////////////////////////////////////////////////

bool endswith(const std::string &str, const std::string &suffix) {
	if (str.length() >= suffix.length()) {
		return (0 == str.compare(str.length() - suffix.length(), suffix.length(), suffix));
	} else {
		return false;
	}
}

void points_dump(const std::string &filename, const voroffset3d::CompressedVolume &dexels)
{
	GEO::Mesh mesh;

	for (int y = 0; y < dexels.gridSize()[1]; ++y) {
		for (int x = 0; x < dexels.gridSize()[0]; ++x) {
			for (int i = 0; 2 * i < dexels.at(x, y).size(); ++i) {
				GEO::vec3 xyz_min, xyz_max;
				xyz_min[0] = dexels.origin()[0] + (x + 0.5) * dexels.spacing();
				xyz_min[1] = dexels.origin()[1] + (y + 0.5) * dexels.spacing();
				xyz_min[2] = dexels.at(x, y)[2 * i + 0] * dexels.spacing();
				xyz_max[0] = dexels.origin()[0] + (x + 0.5) * dexels.spacing();
				xyz_max[1] = dexels.origin()[1] + (y + 0.5) * dexels.spacing();
				xyz_max[2] = dexels.at(x, y)[2 * i + 1] * dexels.spacing();
				mesh.vertices.create_vertex(&xyz_min[0]);
				mesh.vertices.create_vertex(&xyz_max[0]);
			}
		}
	}

	GEO::mesh_save(mesh, filename);
}

void voroffset3d::dexel_dump(const std::string &filename, const CompressedVolume &dexels)
{
	if (endswith(filename, ".xyz")) {
		points_dump(filename, dexels);
		return;
	}
	GEO::Mesh mesh;

	for (int y = 0; y < dexels.gridSize()[1]; ++y) {
		for (int x = 0; x < dexels.gridSize()[0]; ++x) {
			for (int i = 0; 2 * i < dexels.at(x, y).size(); ++i) {
				GEO::vec3 xyz_min, xyz_max;
				xyz_min[0] = dexels.origin()[0] + x * dexels.spacing();
				xyz_min[1] = dexels.origin()[1] + y * dexels.spacing();
				xyz_min[2] = dexels.at(x, y)[2 * i + 0] * dexels.spacing();
				xyz_max[0] = dexels.origin()[0] + (x + 1) * dexels.spacing();
				xyz_max[1] = dexels.origin()[1] + (y + 1) * dexels.spacing();
				xyz_max[2] = dexels.at(x, y)[2 * i + 1] * dexels.spacing();
				GEO::vec3 diff[8] = {
					GEO::vec3(0,0,1), GEO::vec3(1,0,1), GEO::vec3(0,1,1), GEO::vec3(1,1,1),
					GEO::vec3(0,0,0), GEO::vec3(1,0,0), GEO::vec3(0,1,0), GEO::vec3(1,1,0),
				};
				int v = mesh.vertices.nb();
				for (int lv = 0; lv < 8; ++lv) {
					for (int d = 0; d < 3; ++d) {
						diff[lv][d] = xyz_min[d] + diff[lv][d] * (xyz_max[d] - xyz_min[d]);
					}
					mesh.vertices.create_vertex(&diff[lv][0]);
				}
				mesh.cells.create_hex(v, v + 1, v + 2, v + 3, v + 4, v + 5, v + 6, v + 7);
			}
		}
	}

	mesh.cells.compute_borders();
	mesh.cells.connect();
	mesh.vertices.remove_isolated();

	GEO::mesh_save(mesh, filename);
}

