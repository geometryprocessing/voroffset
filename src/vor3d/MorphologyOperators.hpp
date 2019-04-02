#include"vor3d/MorphologyOperators.h"


namespace voroffset3d
{
	// add segment at the end of the sorted segment vector
	template<typename Container>
	void addSegmentAtTheEnd(Container & dst, double start, double end)
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

	template<typename Container, typename SegmentType>
	void addSegmentAtTheEnd(Container & dst, SegmentType segment)
	{
		if (dst.empty())
		{
			dst.push_back(segment);
			return;
		}
		else if (dst.back().isOcclude(segment)) // Attention: we have dst.back().y1<=segment.y1
		{
			return;
		}
		else if (segment.isOcclude(dst.back())) // Attention: we have dst.back().y1<=segment.y1
		{
			dst.pop_back();
			dst.push_back(segment);
			return;
		}
		else if (dst.back().y2 < segment.y1)
		{
			dst.push_back(segment);
			return;
		}
		else if (segment.y2 > dst.back().y2)
		{
			if (dst.back().isSplit(segment))
			{
				dst.back().y2 = segment.y1;
				dst.push_back(segment);
				return;
			}
			else
			{
				segment.y1 = dst.back().y2;
				dst.push_back(segment);
				return;
			}
		}
		else // segment is occluded by dst.back(), for we have assumed that dst.back()._begin_pt <= segment._begin_pt
		{
			if (dst.back().isSplit(segment))
			{
				SegmentType new_segment(dst.back());
				new_segment.y1 = segment.y2;
				dst.back().y2 = segment.y1;
				dst.push_back(segment);
				dst.push_back(new_segment);
				return;
			}
		}
	}

	// Computes the union of two sorted and nonoverlapping lists of segments.
	// output is also sorted and nonoverlapping segments
	template<typename SegmentType>
	void unionSegs(const std::vector<SegmentType> & a, const std::vector<SegmentType> & b
		, std::vector<SegmentType> & result)
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
		auto ia = a.begin();
		auto ib = b.begin();
		SegmentType a_1, b_1;
		a_1 = *ia++;
		b_1 = *ib++;
		bool bfull(true), afull(true);
		while (bfull || afull)
		{
			if (bfull && b_1.y1 <= a_1.y1)
			{
				addSegmentAtTheEnd(result, b_1);
				if (ib == b.end())
				{
					bfull = false;
					b_1.y1 = std::numeric_limits<int>::max();
				}
				else
					b_1 = *ib++;
			}
			else if (afull)
			{
				addSegmentAtTheEnd(result, a_1);
				if (ia == a.end())
				{
					afull = false;
					a_1.y1 = std::numeric_limits<int>::max();
				}
				else
					a_1 = *ia++;
			}
		}
		// clear the capacity
		/*a.clear();
		a.shrink_to_fit();
		b.clear();
		b.shrink_to_fit();*/
	}
}
///////////////////////////////////////////////////////////////////////////////

