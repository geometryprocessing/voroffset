#include"vor2d/MorphologyOperators.h"

using namespace voroffset;

#ifndef MIN_SEG_SIZE
#define MIN_SEG_SIZE 1e-10
#endif

// Appends segment [a,b] to the given sorted line
void voroffset::appendSegment(std::vector<double> &nl, double a, double b)
{
	// Ensures nl[-2] < a
	while (!nl.empty() && nl.rbegin()[1] >= a)
	{
		b = std::max(b, nl.back());
		nl.pop_back();
		nl.pop_back();
	}
	if (nl.empty() || nl.back() < a)
	{
		nl.push_back(a);
		nl.push_back(b);
	}
	else if (nl.back() < b)
	{
		nl.back() = b;
	}
	vor_assert(nl.size() % 2 == 0);
}

void voroffset::appendSegment(std::vector<m_Segment> &n1, m_Segment seg)
{
	while (!n1.empty() && n1.back()._begin_pt >= seg._begin_pt)
	{
		seg._end_pt = std::max(seg._end_pt, n1.back()._end_pt);
		n1.pop_back();
	}
	if (n1.empty() || n1.back()._end_pt < seg._begin_pt)
		n1.push_back(seg);
	else if (n1.back()._end_pt < seg._end_pt)
		n1.back()._end_pt = seg._end_pt;
}
//////////////////////////////////////////////////////////////////////////////

// add segment at the end of the sorted segment vector
template<typename Container>
void voroffset::addSegmentAtTheEnd(Container & dst, double start, double end)
{
	if (dst.empty() || (dst.back() < start))
	{
		dst.push_back(start);
		dst.push_back(end);
	}
	else if (end > dst.back())
	{
		dst.back() = end;
	}
}

template<typename Container>
void voroffset::addSegmentAtTheEnd(Container & dst, m_Segment segment)
{
	if (dst.empty() || (dst.back()._end_pt < segment._begin_pt))
		dst.push_back(segment);
	else if (segment._end_pt > dst.back()._end_pt)
	{
		dst.back()._end_pt = segment._end_pt;
	}
}
////////////////////////////////////////////////////////////////////

// get the next segment of the given position
void voroffset::getNextSegment(const std::vector<double> & segs, std::vector<double>::const_iterator & it,
	double & left, double & right, bool & full)
{
	if (it == segs.end())
	{
		full = false;
	}
	else
	{
		left = *it++;
		right = *it++;
	}
}

void voroffset::getNextSegment(const std::vector<m_Segment> & segs, std::vector<m_Segment>::const_iterator & it,
	m_Segment &next_seg, bool & full)
{
	if (it == segs.end())
	{
		full = false;
	}
	else
	{
		next_seg = *it++;
	}
}
////////////////////////////////////////////////////////////////////////////////////

// Computes the union of two sorted lists of segments.
// Uses a parallel sweep of both lists, akin to merge-sort.
void voroffset::unionSegs(const std::vector<double> & a, const std::vector<double> & b, std::vector<double> & result)
{
	if (a.empty())
	{
		result.insert(result.end(), b.begin(), b.end());
		return;
	}
	if (b.empty())
	{
		result.insert(result.end(), a.begin(), a.end());
		return;
	}
	std::vector<double>::const_iterator ia(a.begin()), ib(b.begin());
	double al, ar, bl, br;
	al = *ia++; ar = *ia++;
	bl = *ib++; br = *ib++;
	bool bfull(true), afull(true);
	while (bfull || afull)
	{
		if (bfull && bl <= al)
		{
			addSegmentAtTheEnd(result, bl, br);
			getNextSegment(b, ib, bl, br, bfull);
		}
		else if (afull)
		{
			addSegmentAtTheEnd(result, al, ar);
			if (ia == a.end())
			{
				afull = false;
				al = std::numeric_limits<int>::max();
			}
			else
			{
				al = *ia++;
				ar = *ia++;
			}
		}
	}
}

void voroffset::unionSegs(const std::vector<m_Segment> & a, const std::vector<m_Segment> & b, std::vector<m_Segment> & result)
{
	if (a.empty())
	{
		result.insert(result.end(), b.begin(), b.end());
		return;
	}
	if (b.empty())
	{
		result.insert(result.end(), a.begin(), a.end());
		return;
	}
	std::vector<m_Segment>::const_iterator ia(a.begin()), ib(b.begin());
	m_Segment a_1, b_1;
	a_1 = *ia++;
	b_1 = *ib++;
	bool bfull(true), afull(true);
	while (bfull || afull)
	{
		if (bfull && b_1._begin_pt <= a_1._begin_pt)
		{
			addSegmentAtTheEnd(result, b_1);
			getNextSegment(b, ib, b_1, bfull);
		}
		else if (afull)
		{
			addSegmentAtTheEnd(result, a_1);
			if (ia == a.end())
			{
				afull = false;
				a_1._begin_pt = std::numeric_limits<int>::max();
			}
			else
				a_1 = *ia++;
		}
	}
}
//////////////////////////////////////////////////////////////////////////////

// Unoin map
/*

template<typename Iterator1, typename Iterator2>
void voroffset::unionMap(Iterator1 a, Iterator2 b, std::vector<SegVector> &res)
{
	for (size_t i = 0; i < res.size(); ++i)
	{
		SegVector temp;
		unionSegs(a[i], b[i], temp);
		res[i].swap(temp);
	}
}*/
//////////////////////////////////////////////////////////////////////////////////

// negate ray
void voroffset::negate_ray(std::vector<double> &result, double y_min, double y_max)
{
	int size_ = result.size();

	if (result.empty())
	{
		result = { y_min, y_max };
	}
	else
	{
		// First event
		if (result[0] == y_min)
		{
			result.erase(result.begin());
		}
		else
		{
			result.insert(result.begin(), y_min);
		}
		// Last event
		if (result.back() == y_max)
		{
			result.pop_back();
		}
		else
		{
			result.push_back(y_max);
		}
	}
}

void voroffset::negate_ray(std::vector<m_Segment> input, std::vector<m_Segment> &result, double y_min, double y_max)
{
	int size_ = input.size();
	result.clear();
	if (input.empty())
	{
		result = { m_Segment(y_min, y_max ,1) };
	}
	else
	{
		// First event
		if (input[0]._begin_pt > y_min)
			result.push_back(m_Segment(y_min, input[0]._begin_pt, 1));
		for (int i = 0; i < size_ - 1; i++)
			result.push_back(m_Segment(input[i]._end_pt, input[i + 1]._begin_pt, 1));
		// End event
		if (input[size_ - 1]._end_pt < y_max)
			result.push_back(m_Segment(input[size_ - 1]._end_pt, y_max, 1));
	}
}
////////////////////////////////////////////////////////////////////////////////////////

// calculate the xor between two vectors of segments
void voroffset::calculate_ray_xor(std::vector<m_Segment> ray_1, std::vector<m_Segment> ray_2, std::vector<m_Segment> &ray_xor, double y_min, double y_max)
{
	ray_xor.clear();
	std::vector<m_Segment> record_1, record_2, tmpt_union_1, tmpt_union_2, before_remove_pts;
	negate_ray(ray_1, record_1, y_min, y_max);
	negate_ray(ray_2, record_2, y_min, y_max);
	unionSegs(record_1, ray_2, tmpt_union_1);
	unionSegs(record_2, ray_1, tmpt_union_2);
	negate_ray(tmpt_union_1, ray_1, y_min, y_max);
	negate_ray(tmpt_union_2, ray_2, y_min, y_max);
	unionSegs(ray_1, ray_2, before_remove_pts);
	removepoint(before_remove_pts, ray_xor);
}

void voroffset::calculate_ray_xor(std::vector<double> ray_1, std::vector<double> ray_2, std::vector<double> &ray_xor, double y_min, double y_max)
{
	ray_xor.clear();
	std::vector<double> record_1, record_2, tmpt_union_1, tmpt_union_2, before_remove_pts;
	record_1.resize(ray_1.size());
	record_2.resize(ray_2.size());
	std::copy_n(ray_1.begin(), ray_1.size(), record_1.begin());
	std::copy_n(ray_2.begin(), ray_2.size(), record_2.begin());
	negate_ray(record_1, y_min, y_max);
	negate_ray(record_2, y_min, y_max);
	unionSegs(record_1, ray_2, tmpt_union_1);
	unionSegs(record_2, ray_1, tmpt_union_2);
	negate_ray(tmpt_union_1, y_min, y_max);
	negate_ray(tmpt_union_2, y_min, y_max);
	unionSegs(tmpt_union_1, tmpt_union_2, before_remove_pts);
	removepoint(before_remove_pts, ray_xor);
}
/////////////////////////////////////////////////////////////////////////////

// remove the segment when its size is too small to be seen
void voroffset::removepoint(const std::vector<m_Segment> segs, std::vector<m_Segment> &res)
{
	int size_ = segs.size();
	for (int i = 0; i < size_; i++)
	{
		if (segs[i]._end_pt - segs[i]._begin_pt < MIN_SEG_SIZE)
			continue;
		res.push_back(segs[i]);
	}
}
void voroffset::removepoint(const std::vector<double> segs, std::vector<double> &res)
{
	int size_ = segs.size();
	for (int i = 0; i < size_; i += 2)
	{
		if (segs[i + 1] - segs[i] < MIN_SEG_SIZE)
			continue;
		res.push_back(segs[i]);
		res.push_back(segs[i + 1]);
	}
}