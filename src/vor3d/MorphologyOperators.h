#pragma once
#include "vor3d/Common.h"
// Dealing with the operation on the same sweepline, such as union segments, negate ray, append segment to the line and so on

namespace voroffset3d 
{
	// Appends segment [a,b] to the given sorted line
	void appendSegment(std::vector<double> &nl, double a, double b);
	void appendSegment(std::vector<SegmentWithRadius> &n1, SegmentWithRadius seg);

	// add segment at the end of the sorted segment vector
	template<typename Container>
	void addSegmentAtTheEnd(Container & dst, double start, double end);
	template<typename Container, typename SegmentType>
	void addSegmentAtTheEnd(Container & dst, SegmentType segment);
	

	// get the next segment of the given position
	void getNextSegment(const std::vector<double> & segs, std::vector<double>::const_iterator & it,
		double & left, double & right, bool & full);

	// Computes the union of two sorted lists of segments.
	// Uses a parallel sweep of both lists, akin to merge-sort.
	void unionSegs(const std::vector<double> & a, const std::vector<double> & b, std::vector<double> & result);
	template<typename SegmentType>
	void unionSegs(const std::vector<SegmentType> & a, const std::vector<SegmentType> & b, std::vector<SegmentType> & result);

	void unionPoints(const std::vector<PointWithRadius> & a, const std::vector<PointWithRadius> & b, std::vector<PointWithRadius> & result);

	// negate ray which is occluded by [y_min, y_max]
	void negate_ray(std::vector<double> &result, double y_min, double y_max);
	void negate_ray(std::vector<SegmentWithRadius> input, std::vector<SegmentWithRadius> &result, double y_min, double y_max);

	// negate ray based on the range [y_min, y_max], i.e. delete the segments which don't overlap with [y_min, y_max]
	// as well as negate the segments which are occluded by [y_min, y_max]
	void negate_ray_range(std::vector<double> &result, double y_min, double y_max);

	// calculate the xor of two rays
	void calculate_ray_xor(std::vector<SegmentWithRadius> ray_1, std::vector<SegmentWithRadius> ray_2, 
		std::vector<SegmentWithRadius> &ray_xor, double y_min, double y_max);
	void calculate_ray_xor(std::vector<double> ray_1, std::vector<double> ray_2, std::vector<double> &ray_xor, double y_min, double y_max);

	// remove the segment when its size is too small to be seen, which can modify the viwer of visualization
	void removepoint(const std::vector<SegmentWithRadius> segs, std::vector<SegmentWithRadius> &res);
	void removepoint(const std::vector<double> segs, std::vector<double> &res);

	// coordinate convertion
	int index3dToIndex2d(int i, int j, int xsize, int ysize);
	Eigen::Vector2i index2dToIndex3d(int index2d, int xsize, int ysize);

}//namespace voroffset

#include "MorphologyOperators.hpp"
