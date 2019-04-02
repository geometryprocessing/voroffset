////////////////////////////////////////////////////////////////////////////////
#include "vor3d/SeparatePower2D.h"
#include "vor3d/MorphologyOperators.h"
#include <set>
#include <vector>
#include <complex>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <set>
////////////////////////////////////////////////////////////////////////////////

using namespace voroffset3d;

////////////////////////////////////////////////////////////////////////////////

namespace
{

	// -----------------------------------------------------------------------------

	template<typename T>
	T sign(T a)
	{
		return a > 0 ? T(1) : T(-1);
	}


	// Auxiliary methods
	template<typename T>
	T det(T a, T b, T c, T d)
	{
		return a * d - b * c;
	}
	// Returns true iff a ¡Ü sqrt(x)
	template<typename T>
	bool infSqrt(T a, T x)
	{
		return (a <= 0 || a*a <= x);
	}

	// Returns true iff sqrt(x) ¡Ü b
	template<typename T>
	bool supSqrt(T b, T x)
	{
		return (b >= 0 && x <= b*b);
	}

	static inline bool ceilIfGreateq(double x, double xmin, double &xdst)
	{
		if (x >= xmin)
		{
			xdst = x;
			return true;
		}
		else
		{
			return false;
		}
	}

	// Typedefs
	typedef SegmentWithRadiusX Segment;
	typedef PointWithRadiusX Point;

	////////////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////////////////////

	//template<typename Scalar>
	struct Ray
	{
		typedef double Scalar;

		Scalar a, b, c;

		Ray(Scalar _a, Scalar _b, Scalar _c)
			: a(_a), b(_b), c(_c)
		{}

		// Perpendicular bisector of p and q, equation is a*x + b*y = c
		Ray(Point p, Point q)
			: a(2 * (p.x - q.x))
			, b(2 * (p.y - q.y))
			, c(p.x * p.x - q.x * q.x + p.y * p.y - q.y * q.y + q.r * q.r - p.r * p.r)
		{
			//reduce(a, b, c);
		}

		inline double inter(Scalar x) const
		{
			//return divideDouble(c - a*x, b);
			return (c - a*x) / b;
		}

		inline bool interBefore(int x, double yr) const
		{
			// Precondition: b==0
			return ((c - a*x) * sign(b) <= yr * b * sign(b));
		}

		inline bool interAfter(int x, double yl) const
		{
			// Precondition: b==0
			return (yl * b * sign(b) <= (c - a*x) * sign(b));
		}
	};

} // anonymous namespace

  ////////////////////////////////////////////////////////////////////////////////
  // VoronoiMorpho
  ////////////////////////////////////////////////////////////////////////////////

  /*
  * Assumes that y_p < y_q < y_r
  */
double SeparatePowerMorpho2D::rayIntersect(PointWithRadiusX p, PointWithRadiusX q, PointWithRadiusX r) const
{
	if (std::abs(p.x - q.x) < EPS && std::abs(q.x - r.x) < EPS)
	{
		return m_XMax*1.0f;
	}
	else
	{
		Ray r1(p, q);
		Ray r2(q, r);
		if (r1.b*r2.b < 0)
		{
			r1.a = -r1.a;
			r1.b = -r1.b;
			r1.c = -r1.c;
		}
		const double denom = det(r1.a, r2.a, r1.b, r2.b);
		//std::cout << r1.a << " " << r1.b << " " << r1.c << " , " << r2.a << " " << r2.b << " " << r2.c << std::endl;
		if (std::abs(denom) < EPS)
		{
			// Intersection at x >= m_XMax are discarded;
			return m_XMax;
		}
		else if (denom > 0)
		{
			return m_XMax;
		}
		else
		{
			const double numer = det(r1.c, r2.c, r1.b, r2.b);
			if (numer / denom < q.x + EPS || numer / denom >(q.x + q.r))
			{
				return m_XMax;
			}
			else
			{
				return ((double)numer / denom);
			}
		}
	}
}

// -----------------------------------------------------------------------------

// Assumes that it != m_S.end()
void SeparatePowerMorpho2D::exploreLeft(P_const_iter it, PointWithRadiusX p, int i)
{
	while (it != m_P.begin())
	{
		double xinter = rayIntersect(*std::prev(it), *it, p);
		if (xinter < i + EPS)
		{
			it = std::prev(m_P.erase(it));
		}
		else
		{
			int x_sup = (int)std::floor(xinter) + 1;
			if (xinter < m_XMax)
			{
				if (x_sup < m_XMax)
				{
					m_Q_P[x_sup].push_back(*it);
				}
			}
			return;
		}
	}
}

// Assumes that it != m_S.end()
void SeparatePowerMorpho2D::exploreRight(P_const_iter it, PointWithRadiusX p, int i)
{
	while (std::next(it) != m_P.end())
	{
		double xinter = rayIntersect(p, *it, *std::next(it));
		if (xinter < i + EPS)
		{
			it = m_P.erase(it);
		}
		else
		{
			int x_sup = (int)std::floor(xinter) + 1;
			if (xinter < m_XMax)
			{
				if (x_sup < m_XMax)
				{
					m_Q_P[x_sup].push_back(*it);
				}
			}
			return;
		}
	}
}

// -----------------------------------------------------------------------------

// Insert a new seed segment (i,j1) -- (i,j2)
void SeparatePowerMorpho2D::insertSegment(int i, double j1, double j2, double r)
{
	if (std::abs(j1 - j2) < EPS)
	{
		m_P_Tmp.push_back(PointWithRadiusX(i, j1, r));
		return;
	}
	SegmentWithRadiusX new_seg(i, j1, j2, r);
	m_S_Tmp.push_back(new_seg);
	m_P_Tmp.push_back(new_seg.left());
	m_P_Tmp.push_back(new_seg.right());
}

void SeparatePowerMorpho2D::insertPoint(PointWithRadiusX p)
{
	//PointWithRadiusX p(i, j1, r);
	P_const_iter it = m_P.lower_bound(p);
	while (it != m_P.end() && p.contains(*it))
	{
		it = m_P.erase(it);
	}
	while (it != m_P.begin() && p.contains(*std::prev(it)))
	{
		m_P.erase(std::prev(it));
	}
	P_const_iter jt = m_P.insert(it, p);
	int xsup = (int)std::floor(p.x + p.r) + 1;
	if (xsup < m_XMax) { m_Q_P[xsup].push_back(p); }
	it = std::next(jt);
	if (it != m_P.end())
	{
		exploreRight(it, p, (int)p.x);
	}
	if (jt != m_P.begin())
	{
		exploreLeft(std::prev(jt), p, (int)p.x);
	}

}

// -----------------------------------------------------------------------------

// Remove seeds that are not contributing anymore to the current sweep line (ith)
void SeparatePowerMorpho2D::removeInactiveSegments(int i)
{
	// remove inactive segments
	int k, j;
	k = 0; j = 0;
	size_t size = m_S.size();
	for (; k < size; k++)
	{
		if ((int)std::floor(m_S[k].x + m_S[k].r) + 1 > i)
		{
			m_S[j] = m_S[k];
			j++;
		}
	}
	m_S.resize(j);
}

void SeparatePowerMorpho2D::removeInactivePoints(int i)
{
	std::sort(m_Q_P[i].begin(), m_Q_P[i].end());
	while (!m_Q_P[i].empty())
	{
		P_const_iter it = m_P.find(m_Q_P[i].back());
		m_Q_P[i].pop_back();
		if (it != m_P.end())
		{
			it = m_P.erase(it);
			if (it != m_P.begin() && it != m_P.end())
			{
				exploreLeft(std::prev(it), *it, i);
				exploreRight(it, *std::prev(it), i);
			}
		}
	}
	m_Q_P[i].clear();
}

void SeparatePowerMorpho2D::getLine(bool is_set_radii, int posX, int posY, int deltaX, int deltaY, int i,
	std::function<void(int, int, double, double, double)> _appendSegment)
{
	int pos_x = posX + i*deltaX;
	int pos_y = posY + i*deltaY;
	ray_P.clear();
	ray_S.clear();
	ray_U.clear();
	for (SegmentWithRadiusX s : m_S)
	{
		vor_assert((i - s.x) <= s.r + EPS);
		double a, b;
		if (i - s.x > s.r)
		{
			a = s.y1;
			b = s.y2;
		}
		else
		{
			const double dy = std::sqrt(s.r * s.r - (i - s.x) * (i - s.x));
			a = s.y1 - dy;
			b = s.y2 + dy;
		}
		
		appendSegment(ray_S, a, b); //all the radii are same for the first phase
	}
	for (PointWithRadiusX p : m_P)
	{
		vor_assert((i - p.x) <= p.r + EPS);
		if (i - p.x > p.r)
		{
			continue;
		}
		else
		{
			const double dy = std::sqrt(p.r * p.r - (i - p.x) * (i - p.x));
			double a = p.y - dy;
			double b = p.y + dy;
			appendSegment(ray_P, a, b);
		}
	}
	vor3d::unionSegs(ray_P, ray_S, ray_U);
	for (int k = 0; k < ray_U.size(); k += 2)
	{
		_appendSegment(pos_x, pos_y, ray_U[k], ray_U[k + 1], 0.f);
	}
}

void SeparatePowerMorpho2D::resetData()
{
	m_S.clear();
	m_S_New.clear();
	m_S_Tmp.clear();
	m_P.clear();
	current_line = -1;
}
void SeparatePowerMorpho2D::unionSegs()
{
	if (m_S_New.size())
		m_S_New.clear();
	if (m_S.empty())
	{
		m_S_New.insert(m_S_New.end(), m_S_Tmp.begin(), m_S_Tmp.end());
		return;
	}
	if (m_S_Tmp.empty())
	{
		m_S_New.insert(m_S_New.end(), m_S.begin(), m_S.end());
		return;
	}
	//std::vector<SegmentWithRadiusX>::const_iterator ia(m_S.begin()), ib(m_S_Tmp.begin());
	size_t size_a, size_b;
	size_a = m_S.size();
	size_b = m_S_Tmp.size();
	int i, j;
	i = 0; j = 0;
	SegmentWithRadiusX a_1(m_S[0]), b_1(m_S_Tmp[0]);
	i++;
	j++;
	bool bfull(true), afull(true);
	while (bfull || afull)
	{
		if (afull && a_1.y1 <= b_1.y1)
		{
			addSegment(a_1);
			if (i == size_a)
			{
				afull = false;
				a_1.y1 = std::numeric_limits<int>::max();
			}
			else
				a_1 = SegmentWithRadiusX(m_S[i++]);
		}
		else
		{
			// appendsegment
			addSegment(b_1);
			if (j == size_b)
			{
				bfull = false;
				b_1.y1 = std::numeric_limits<int>::max();
			}
			else
				b_1 = SegmentWithRadiusX(m_S_Tmp[j++]);
		}
	}
}

void SeparatePowerMorpho2D::addSegment(SegmentWithRadiusX segment)
{
	if (segment.isPoint())
		return;
	if (m_S_New.empty())
	{
		m_S_New.push_back(segment);
	}
	else if (m_S_New.back().isOcclude(segment)) // Attention: we have dst.back().y1<=segment.y1
	{
		// do nothing
	}
	else if (m_S_New.back().y2 < segment.y1)
	{
		m_S_New.push_back(segment);
	}
	else if (segment.y2 > m_S_New.back().y2)
	{
		if (m_S_New.back().isSplit(segment))
		{
			if (std::abs(m_S_New.back().y1 - segment.y1) < EPS)
				m_S_New.pop_back();
			else
			{
				m_S_New.back().y2 = segment.y1;
			}
			m_S_New.push_back(segment);
		}
		else
		{
			segment.y1 = m_S_New.back().y2;
			m_S_New.push_back(segment);
		}
	}
	else // segment is contained by m_S_New.back(), for we have assumed that dst.back()._begin_pt <= segment._begin_pt
	{
		if (m_S_New.back().isSplit(segment))
		{

			SegmentWithRadiusX new_segment(m_S_New.back());
			new_segment.y1 = segment.y2;
			m_S_New.back().y2 = segment.y1;
			if (m_S_New.back().isPoint())
				m_S_New.pop_back();
			m_S_New.push_back(segment);
			if (!new_segment.isPoint())
			{
				m_S_New.push_back(new_segment);
			}
		}
	}
}

void SeparatePowerMorpho2D::flushLine(int i)
{
	// removing inactive points when a new segment comes
	/*for (int k = 0; k < m_S_Tmp.size(); k++)
	{
		auto it = m_P.lower_bound(m_S_Tmp[k].left());
		if (it != m_P.begin() && std::abs(std::prev(it)->y - m_S_Tmp[k].y1) < EPS)
		{
			it = std::prev(it);
		}
		while (it != m_P.end())
		{
			if (it->y > m_S_Tmp[k].y2)
				break;
			double delta_x = m_S_Tmp[k].x - it->x;
			double delta_y = std::min(it->y - m_S_Tmp[k].y1, m_S_Tmp[k].y2 - it->y);
			if (it->r < std::sqrt(delta_x*delta_x + delta_y*delta_y))
			{
				it = m_P.erase(it);
			}
			else
			{
				it = std::next(it);
			}
		}
	}*/
	// removing the points lying on the current sweepline
	for (int k = 0; k < m_P_Tmp.size(); k++)
	{
		int pos = locatePoint(m_P_Tmp[k]);
		if (pos != -1)
		{
			double r_m = std::min(std::min(m_S[pos].r - m_P_Tmp[k].x + m_S[pos].x, m_P_Tmp[k].y - m_S[pos].y1), m_S[pos].y2 - m_P_Tmp[k].y);
			if (m_P_Tmp[k].r <= r_m)
				continue;
			else
				insertPoint(m_P_Tmp[k]);
		}
		else
		{
			insertPoint(m_P_Tmp[k]);
		}
	}
	unionSegs();
	m_S.swap(m_S_New);
	m_S_New.clear();
	m_S_Tmp.clear();
	m_P_Tmp.clear();

}
int SeparatePowerMorpho2D::locatePoint(PointWithRadiusX p)
/*	using binary search to find the first segment in m_S which occuludes the point p
	if succeed, return the index of that segment,
	otherwise, return -1
*/
{
	int left = 0;
	int right = (int)m_S.size() - 1;
	if (right < 0)
		return -1;
	int middle;
	if (p.y > m_S[right].y2 || p.y < m_S[left].y1)
		return -1;
	while (left < right - 1)
	{
		middle = left + ((right - left) >> 1);

		if (m_S[middle].y1 > p.y)
		{
			right = middle;
		}
		else if (m_S[middle].y1 < p.y)
		{
			left = middle;
		}
		else
		{
			left = middle;
			break;
		}
	}
	if (m_S[left].y2 > p.y)
		return left;
	else
		return -1;
}


#ifndef WIN32
// static_assert(std::is_standard_layout<VoronoiMorpho>::value, "VoronoiMorpho is not standard layout!");
#endif