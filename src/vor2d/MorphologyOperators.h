#pragma once
#include "vor2d/Common.h"
// Dealing with the operation on the same sweepline, such as union segments, negate ray, append segment to the line and so on

namespace voroffset 
{
	// Appends segment [a,b] to the given sorted line
	void appendSegment(std::vector<double> &nl, double a, double b);
	void appendSegment(std::vector<m_Segment> &n1, m_Segment seg);

	// add segment at the end of the sorted segment vector
	template<typename Container>
	void addSegmentAtTheEnd(Container & dst, double start, double end);
	template<typename Container>
	void addSegmentAtTheEnd(Container & dst, m_Segment segment);
	

	// get the next segment of the given position
	void getNextSegment(const std::vector<double> & segs, std::vector<double>::const_iterator & it,
		double & left, double & right, bool & full);
	void getNextSegment(const std::vector<m_Segment> & segs, std::vector<m_Segment>::const_iterator & it,
		m_Segment & next_seg, bool & full);

	// Computes the union of two sorted lists of segments.
	// Uses a parallel sweep of both lists, akin to merge-sort.
	void unionSegs(const std::vector<double> & a, const std::vector<double> & b, std::vector<double> & result);
	void unionSegs(const std::vector<m_Segment> & a, const std::vector<m_Segment> & b, std::vector<m_Segment> & result);

	// negate ray
	void negate_ray(std::vector<double> &result, double y_min, double y_max);
	void negate_ray(std::vector<m_Segment> input, std::vector<m_Segment> &result, double y_min, double y_max);

	// calculate the xor of two rays
	void calculate_ray_xor(std::vector<m_Segment> ray_1, std::vector<m_Segment> ray_2, std::vector<m_Segment> &ray_xor, double y_min, double y_max);
	void calculate_ray_xor(std::vector<double> ray_1, std::vector<double> ray_2, std::vector<double> &ray_xor, double y_min, double y_max);

	// remove the segment when its size is too small to be seen, which can modify the viwer of visualization
	void removepoint(const std::vector<m_Segment> segs, std::vector<m_Segment> &res);
	void removepoint(const std::vector<double> segs, std::vector<double> &res);

}//namespace voroffset