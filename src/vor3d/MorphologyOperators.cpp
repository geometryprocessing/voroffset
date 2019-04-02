#include"vor3d/MorphologyOperators.h"
using namespace voroffset3d;

#ifndef MIN_SEG_SIZE
#define MIN_SEG_SIZE 1e-10
#endif



namespace voroffset3d
{

	// Appends segment [a,b] to the given sorted line
	// Assume that inserted segment is larger than all of the segments in the segments vector
	void appendSegment(std::vector<double> &nl, double a, double b)
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

	// Appends segment [a,b] to the given sorted line
	void appendSegment(std::vector<SegmentWithRadius> &n1, SegmentWithRadius seg)
	{
		std::vector<SegmentWithRadius> _segs_vec;
		while (!n1.empty() && n1.back().y1 >= seg.y1)
		{
			if (n1.back().y1 > seg.y2)
			{
				_segs_vec.push_back(n1.back());
			}
			else
			{
				if (n1.back().y2 >= seg.y2)
				{
					if (n1.back().r > seg.r)
					{
						seg.y2 = n1.back().y1;
						_segs_vec.push_back(n1.back());
					}
					else if (n1.back().r == seg.r)
					{
						seg.y2 = n1.back().y2;
					}
					else
					{
						n1.back().y1 = seg.y2;
						_segs_vec.push_back(n1.back());
					}
				}
				else
				{
					if (n1.back().r > seg.r)
					{
						_segs_vec.push_back(vor3d::SegmentWithRadius(n1.back().y2, seg.y2, seg.r));
						_segs_vec.push_back(n1.back());
						seg.y2 = n1.back().y1;
					}
				}
			}
			n1.pop_back();
		}
		if (n1.empty() || n1.back().y2 < seg.y1)
			n1.push_back(seg);
		else if (n1.back().y2 < seg.y2)
		{
			if (n1.back().r == seg.r)
				n1.back().y2 = seg.y2;
			else if (n1.back().r > seg.r)
			{
				seg.y1 = n1.back().y2;
				_segs_vec.push_back(seg);
			}
			else
			{
				n1.back().y2 = seg.y1;
				_segs_vec.push_back(seg);
			}
		}
		n1.insert(n1.end(), _segs_vec.rbegin(), _segs_vec.rend());
	}
	
	// get the next segment of the given position
	void getNextSegment(const std::vector<double> & segs, std::vector<double>::const_iterator & it,
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

	////////////////////////////////////////////////////////////////////////////////////

	// Computes the union of two sorted and nonoverlapping lists of segments.
	// output is also sorted and nonoverlapping segments
	// Uses a parallel sweep of both lists, akin to merge-sort.
	void unionSegs(const std::vector<double> & a, const std::vector<double> & b, std::vector<double> & result)
	{
		if (result.size())
			result.clear();
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
	

	void unionPoints(const std::vector<PointWithRadius> & a, const std::vector<PointWithRadius> & b
		, std::vector<PointWithRadius> & result)
	{
		if (result.size())
			result.clear();
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
		std::vector<PointWithRadius>::const_iterator ia(a.begin()), ib(b.begin());
		PointWithRadius a_1, b_1;
		a_1 = *ia++;
		b_1 = *ib++;
		bool bfull(true), afull(true);
		while (bfull || afull)
		{
			if (bfull && b_1._pt <= a_1._pt)
			{
				result.push_back(b_1);
				if (ib == b.end())
				{
					bfull = false;
					b_1._pt = std::numeric_limits<int>::max();
				}
				else
				{
					b_1 = *ib++;
				}
			}
			else if (afull)
			{
				result.push_back(a_1);
				if (ia == a.end())
				{
					afull = false;
					a_1._pt = std::numeric_limits<int>::max();
				}
				else
				{
					a_1 = *ia++;
				}
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
	void negate_ray(std::vector<double> &result, double y_min, double y_max)
	{
		size_t size_ = result.size();

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

	void negate_ray(std::vector<SegmentWithRadius> input, std::vector<SegmentWithRadius> &result, double y_min, double y_max)
	{
		size_t size_ = input.size();
		result.clear();
		if (input.empty())
		{
			result = { SegmentWithRadius(y_min, y_max , 0) };
		}
		else
		{
			// First event
			if (input[0].y1 > y_min)
				result.push_back(SegmentWithRadius(y_min, input[0].y1, 0));
			for (int i = 0; i < size_ - 1; i++)
				result.push_back(SegmentWithRadius(input[i].y2, input[i + 1].y1, 0));
			// End event
			if (input[size_ - 1].y2 < y_max)
				result.push_back(SegmentWithRadius(input[size_ - 1].y2, y_max, 0));
		}
	}

	void negate_ray_range(std::vector<double> &result, double y_min, double y_max)
	{
		size_t size_ = result.size();

		if (result.empty())
		{
			result = { y_min, y_max };
		}
		else
		{
			int count_first = 0;
			int count_last = 0;
			while (result[0] <= y_min)
			{
				result.erase(result.begin());
				count_first++;
				if(result.empty())
					break;
			}
			if (count_first % 2 == 0)
				result.insert(result.begin(), y_min);

			while (result.back() >= y_max)
			{
				result.pop_back();
				count_last++;
				if(result.empty())
					break;
			}
			if (count_last % 2 == 0)
				result.push_back(y_max);

		}
	}
	////////////////////////////////////////////////////////////////////////////////////////

	// calculate the xor between two vectors of segments
	void calculate_ray_xor(std::vector<SegmentWithRadius> ray_1, std::vector<SegmentWithRadius> ray_2
		, std::vector<SegmentWithRadius> &ray_xor, double y_min, double y_max)
	{
		ray_xor.clear();
		std::vector<SegmentWithRadius> record_1, record_2, tmp_union_1, tmp_union_2, before_remove_pts;
		negate_ray(ray_1, record_1, y_min, y_max);
		negate_ray(ray_2, record_2, y_min, y_max);
		unionSegs(record_1, ray_2, tmp_union_1);
		unionSegs(record_2, ray_1, tmp_union_2);
		negate_ray(tmp_union_1, ray_1, y_min, y_max);
		negate_ray(tmp_union_2, ray_2, y_min, y_max);
		unionSegs(ray_1, ray_2, before_remove_pts);
		removepoint(before_remove_pts, ray_xor);
	}

	void calculate_ray_xor(std::vector<double> ray_1, std::vector<double> ray_2, std::vector<double> &ray_xor, double y_min, double y_max)
	{
		ray_xor.clear();
		std::vector<double> record_1, record_2, tmp_union_1, tmp_union_2, before_remove_pts;
		record_1.resize(ray_1.size());
		record_2.resize(ray_2.size());
		std::copy_n(ray_1.begin(), ray_1.size(), record_1.begin());
		std::copy_n(ray_2.begin(), ray_2.size(), record_2.begin());
		negate_ray(record_1, y_min, y_max);
		negate_ray(record_2, y_min, y_max);
		unionSegs(record_1, ray_2, tmp_union_1);
		unionSegs(record_2, ray_1, tmp_union_2);
		negate_ray(tmp_union_1, y_min, y_max);
		negate_ray(tmp_union_2, y_min, y_max);
		unionSegs(tmp_union_1, tmp_union_2, before_remove_pts);
		removepoint(before_remove_pts, ray_xor);
	}
	/////////////////////////////////////////////////////////////////////////////

	// remove the segment when its size is too small to be seen
	void removepoint(const std::vector<SegmentWithRadius> segs, std::vector<SegmentWithRadius> &res)
	{
		size_t size_ = segs.size();
		for (int i = 0; i < size_; i++)
		{
			if (segs[i].y2 - segs[i].y1 < MIN_SEG_SIZE)
				continue;
			res.push_back(segs[i]);
		}
	}
	void removepoint(const std::vector<double> segs, std::vector<double> &res)
	{
		size_t size_ = segs.size();
		for (int i = 0; i < size_; i += 2)
		{
			if (segs[i + 1] - segs[i] < MIN_SEG_SIZE)
				continue;
			res.push_back(segs[i]);
			res.push_back(segs[i + 1]);
		}
	}


	int index3dToIndex2d(int i, int j, int xsize, int ysize)
	{
		vor_assert(i < xsize && j < ysize && i >= 0 && j >= 0);
		return i + j*xsize;
	}

	Eigen::Vector2i index2dToIndex3d(int index, int xsize, int ysize)
	{
		vor_assert(index < xsize*ysize && index >= 0);
		Eigen::Vector2i res(index%xsize, index / xsize);
		return res;
	}

}