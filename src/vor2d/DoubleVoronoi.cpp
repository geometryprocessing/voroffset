////////////////////////////////////////////////////////////////////////////////
#include "vor2d/DoubleVoronoi.h"
#include <set>
#include <vector>
#include <complex>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <set>
////////////////////////////////////////////////////////////////////////////////

using namespace voroffset;

////////////////////////////////////////////////////////////////////////////////

namespace
{

	// -----------------------------------------------------------------------------

	template<typename T>
	T sign(T a) { return a > 0 ? T(1) : T(-1); }

	// Works only for non-negative numbers (≥ 0)
	/*template<typename IntType>
	IntType ceilDiv(IntType x, IntType y) {
		return (x / y + (x % y != 0));
	}*/

	// Auxiliary methods
	template<typename T>
	T det(T a, T b, T c, T d)
	{
		return a * d - b * c;
	}

	// Returns true iff a ≤ sqrt(x)
	template<typename T>
	bool infSqrt(T a, T x)
	{
		return (a <= 0 || a*a <= x);
	}

	// Returns true iff sqrt(x) ≤ b
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
	typedef VoronoiMorphoF::SegmentF Segment;
	//typedef double real;

	double X(PointF p) { return p.real(); }
	double Y(PointF p) { return p.imag(); }

	////////////////////////////////////////////////////////////////////////////////

	//template<typename Scalar>
	struct Ray
	{
		typedef double Scalar;

		Scalar a, b, c;

		Ray(Scalar _a, Scalar _b, Scalar _c)
			: a(_a), b(_b), c(_c)
		{ }

		// Perpendicular bisector of p and q, equation is a*x + b*y = c
		Ray(PointF p, PointF q)
			: a(2 * (X(p) - X(q)))
			, b(2 * (Y(p) - Y(q)))
			, c(X(p) * X(p) - X(q) * X(q) + Y(p) * Y(p) - Y(q) * Y(q))
		{}

		inline double inter(Scalar x) const
		{
			return (c - a*x)/b;
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

	////////////////////////////////////////////////////////////////////////////////

	/*
	* Assumes that y_p < y_q < y_r
	* Returns true if the intersection is between y1 and y2
	*/
	bool ray_intersect_range(PointF p, PointF q, PointF r, double y1, double y2, double &xinter)
	{
		if (X(p) == X(q) && X(q) == X(r))
		{
			return false;
		}
		else
		{
			const Ray r1(p, q);
			const Ray r2(q, r);
			const double denom = det(r1.a, r2.a, r1.b, r2.b);
			if (denom == 0)
			{
				return false;
			}
			else
			{
				const double ynumer = det(r1.a, r2.a, r1.c, r2.c);
				if (y1 * denom * sign(denom) <= ynumer * sign(denom)
					&& ynumer * sign(denom) <= y2 * denom * sign(denom))
				{
					const double xnumer = det(r1.c, r2.c, r1.b, r2.b);
					if (xnumer * sign(denom) > X(q) * denom * sign(denom))
					{
						xinter = ((double)xnumer / denom);
						return true;
					}
				}
				return false;
			}
		}
	}

	bool ray_intersect_after(PointF p, PointF q, PointF r, double y, double &xinter)
	{
		if (X(p) == X(q) && X(q) == X(r))
		{
			return false;
		}
		else
		{
			const Ray r1(p, q);
			const Ray r2(q, r);
			const double denom = det(r1.a, r2.a, r1.b, r2.b);
			if (denom == 0)
			{
				return false;
			}
			else
			{
				const double ynumer = det(r1.a, r2.a, r1.c, r2.c);
				if (y * denom * sign(denom) <= ynumer * sign(denom))
				{
					const double xnumer = det(r1.c, r2.c, r1.b, r2.b);
					if (xnumer * sign(denom) > X(q) * denom * sign(denom))
					{
						xinter = ((double)xnumer / denom);
						return true;
					}
				}
				return false;
			}
		}
	}

	bool ray_intersect_before(PointF p, PointF q, PointF r, double y, double &xinter)
	{
		if (X(p) == X(q) && X(q) == X(r))
		{
			return false;
		}
		else
		{
			const Ray r1(p, q);
			const Ray r2(q, r);
			const double denom = det(r1.a, r2.a, r1.b, r2.b);
			if (denom == 0)
			{
				return false;
			}
			else
			{
				const double ynumer = det(r1.a, r2.a, r1.c, r2.c);
				if (ynumer * sign(denom) <= y * denom * sign(denom))
				{
					const double xnumer = det(r1.c, r2.c, r1.b, r2.b);
					if (xnumer * sign(denom) > X(q) * denom * sign(denom))
					{
						xinter = ((double)xnumer / denom);
						return true;
					}
				}
				return false;
			}
		}
	}

	/*
	* Finds a point that is at equal distance of p, q and s
	* Assumes that s.x < p.x && s.x < q.x
	*/
	bool voronoi_vertex(Segment s, PointF p, PointF q, double &xinter)
	{
		/*
		* We seek the next voronoi point within the range [s.y1, s.y2].
		* Let (xm, ym) be its coordinates. Per the previous sentence we must have:
		*   (a) xm >= s.x
		*   (b) s.y1 ≤ ym ≤ s.y2
		*
		* Let r be the bisector of p and q. We suppose that the point we seek is on r.
		* As the squared distance from (xm,ym) to s is given by (xm - s.x)^2, we must have:
		*
		* (xm - s.x)^2 = (xm - xp)^2 + (ym - yp)^2  (1)
		*              = (xm - xq)^2 + (ym - yq)^2  (2)
		*
		* The following code solves (1) for ym, and returns true iff (a) and (b)
		*/
		const Ray r(p, q);
		const double u = 2.0 * (X(p) - s.x);
		const double w = s.x * s.x - X(p) * X(p);
		double A = r.a;
		double B = r.b * u - 2.0 * r.a * Y(p);
		double C = r.a * Y(p) * Y(p) - w * r.a - r.c * u;
		//reduce(A, B, C);
		const double delta = B * B - 4.0 * A * C;
		//const ll alpha = 2 * A * A * s.y1 + 2 * A * B;
		//const ll beta = 2 * A * A * s.y2 + 2 * A * B;
		const double alpha = 2.0 * std::abs(A) * s.y1 + sign(A) * B;
		const double beta = 2.0 * std::abs(A) * s.y2 + sign(A) * B;
		const bool pos = (A > 0);
		// We must have alpha ≤ sign(a) * ±sqrt(delta) ≤ beta
		if (A == 0)
		{
			if (2.0 * s.y1 <= Y(p) + Y(q) && Y(p) + Y(q) <= 2.0 * s.y2)
			{
				vor_assert(X(p) != s.x);
				const double xm = ((Y(p) - Y(q)) * (Y(p) - Y(q)) - 4.0 * w) / (4.0 * u);
				return ceilIfGreateq(xm, s.x, xinter);
			}
			else
			{
				return false;
			}
		}
		else if (delta < 0)
		{
			return false;
		}
		else if (delta == 0 && alpha <= 0 && 0 <= beta)
		{
			// root y0 = -B / 2A is within range
			double xm = (2.0 * A * r.c + r.b * B) / (2.0 * r.a * A);
			return ceilIfGreateq(xm, s.x, xinter);
		}
		else if (
			pos ? infSqrt(alpha, delta) && supSqrt(beta, delta)
			: supSqrt(-alpha, delta) && infSqrt(-beta, delta))
		{
			// root y1 = (-B+sqrt(delta)) / 2A is within range
			const double xm = (2.0 * A * r.c + r.b * B - r.b * std::sqrt(delta)) / (2.0 * r.a * A);
			return ceilIfGreateq(xm, s.x, xinter);
		}
		else if (
			!pos ? infSqrt(alpha, delta) && supSqrt(beta, delta)
			: supSqrt(-alpha, delta) && infSqrt(-beta, delta))
		{
			// root y2 = (-B-sqrt(delta)) / 2A is within range
			const double xm = (2.0 * A * r.c + r.b * B + r.b * std::sqrt(delta)) / (2.0 * r.a * A);
			return ceilIfGreateq(xm, s.x, xinter);
		}
		else
		{
			return false;
		}
	}

	// -----------------------------------------------------------------------------

} // anonymous namespace

  ////////////////////////////////////////////////////////////////////////////////
  // VoronoiMorpho
  ////////////////////////////////////////////////////////////////////////////////

  /*
  * Assumes that y_p < y_q < y_r
  */
double VoronoiMorphoF::rayIntersect(PointF p, PointF q, PointF r) const
{
	if (X(p) == X(q) && X(q) == X(r))
	{
		return m_XMax*1.0f;
	}
	else
	{
		const Ray r1(p, q);
		const Ray r2(q, r);
		const double denom = det(r1.a, r2.a, r1.b, r2.b);
		if (denom == 0)
		{
			// Intersection at x >= m_XMax are discarded;
			return m_XMax*1.0f;
		}
		else
		{
			const double numer = det(r1.c, r2.c, r1.b, r2.b);
			if (numer * sign(denom) <= X(q) * denom * sign(denom))
			{
				return m_XMax*1.0f;
			}
			else
			{
				return ((double)numer / denom);
			}
		}
	}
}

/*
* Assumes that y_lp < y_ab < y_qr
*/
double VoronoiMorphoF::treatSegments(SegmentF lp, SegmentF ab, SegmentF qr) const
{
	if (ab.x >= lp.x && ab.x >= qr.x)
	{
		// Nothing to do
		return m_XMax;
	}
	else if (ab.isPoint())
	{
		if (lp.isPoint() && qr.isPoint())
		{

			// Case: lp = point, ab = point, qr = point
			return rayIntersect(lp.any(), ab.any(), qr.any());
		}
		else if (lp.isPoint())
		{
			// Case: lp = point, ab = point, qr = segment
			double xinter = m_XMax*1.0f;
			if (ray_intersect_before(lp.any(), ab.any(), qr.left(), qr.y1, xinter))
			{
				return xinter;
			}
			else if (ray_intersect_after(lp.any(), ab.any(), qr.right(), qr.y2, xinter))
			{
				return xinter;
			}
			else if (voronoi_vertex(qr, lp.any(), ab.any(), xinter))
			{
				return xinter;
			}
			else
			{
				return xinter;
			}
		}
		else if (qr.isPoint())
		{
			// Case: lp = segment, ab = point, qr = point
			double xinter = m_XMax*1.0f;
			if (ray_intersect_before(lp.left(), ab.any(), qr.any(), lp.y1, xinter))
			{
				return xinter;
			}
			else if (ray_intersect_after(lp.right(), ab.any(), qr.any(), lp.y2, xinter))
			{
				return xinter;
			}
			else if (voronoi_vertex(lp, ab.any(), qr.any(), xinter))
			{
				return xinter;
			}
			else
			{
				return xinter;
			}
		}
		else
		{
			// Case: lp = segment, ab = point, qr = segment
			double xinter = m_XMax*1.0f;
			if (ray_intersect_before(lp.left(), ab.any(), qr.left(), lp.y1, xinter))
			{
				return xinter;
			}
			else if (ray_intersect_after(lp.right(), ab.any(), qr.right(), qr.y2, xinter))
			{
				return xinter;
			}
			else if (ray_intersect_range(lp.right(), ab.any(), qr.left(), lp.y2, qr.y1, xinter))
			{
				return xinter;
			}
			else if (voronoi_vertex(lp, ab.any(), qr.left(), xinter))
			{
				return xinter;
			}
			else if (voronoi_vertex(qr, lp.right(), ab.any(), xinter))
			{
				return xinter;
			}
			else
			{
				return xinter;
			}
		}
	}
	else
	{
		if (lp.isPoint() && qr.isPoint())
		{
			// Case: lp = point, ab = segment, qr = point
			double xinter = m_XMax*1.0f;
			if (ray_intersect_before(lp.any(), ab.left(), qr.any(), ab.y1, xinter))
			{
				return xinter;
			}
			else if (ray_intersect_after(lp.any(), ab.right(), qr.any(), ab.y2, xinter))
			{
				return xinter;
			}
			else if (voronoi_vertex(ab, lp.any(), qr.any(), xinter))
			{
				return xinter;
			}
			else
			{
				return xinter;
			}
		}
		else if (lp.isPoint())
		{
			// Case: lp = point, ab = segment, qr = segment
			double xinter = m_XMax*1.0f;
			if (ray_intersect_before(lp.any(), ab.left(), qr.left(), ab.y1, xinter))
			{
				return xinter;
			}
			else if (ray_intersect_after(lp.any(), ab.right(), qr.right(), qr.y2, xinter))
			{
				return xinter;
			}
			else if (ray_intersect_range(lp.any(), ab.right(), qr.left(), ab.y2, qr.y1, xinter))
			{
				return xinter;
			}
			else if (voronoi_vertex(ab, lp.any(), qr.left(), xinter))
			{
				return xinter;
			}
			else if (voronoi_vertex(qr, lp.any(), ab.right(), xinter))
			{
				return xinter;
			}
			else
			{
				return xinter;
			}
		}
		else if (qr.isPoint())
		{
			// Case: lp = segment, ab = segment, qr = point
			double xinter = m_XMax*1.0f;

			if (ray_intersect_before(lp.left(), ab.left(), qr.any(), lp.y1, xinter))
			{
				return xinter;
			}
			else if (ray_intersect_after(lp.right(), ab.right(), qr.any(), ab.y2, xinter))
			{
				return xinter;
			}
			else if (ray_intersect_range(lp.right(), ab.left(), qr.any(), lp.y2, ab.y1, xinter))
			{
				return xinter;
			}
			else if (voronoi_vertex(ab, lp.right(), qr.any(), xinter))
			{
				return xinter;
			}
			else if (voronoi_vertex(lp, ab.left(), qr.any(), xinter)) {

				return xinter;
			}
			else
			{
				return xinter;
			}
		}
		else
		{
			// Case: lp = segment, ab = segment, qr = segment
			double xinter = m_XMax*1.0f;
			if (ray_intersect_before(lp.left(), ab.left(), qr.left(), lp.y1, xinter))
			{
				return xinter;
			}
			else if (ray_intersect_after(lp.right(), ab.right(), qr.right(), qr.y2, xinter))
			{
				return xinter;
			}
			else if (ray_intersect_range(lp.right(), ab.left(), qr.left(), lp.y2, ab.y1, xinter))
			{
				return xinter;
			}
			else if (ray_intersect_range(lp.right(), ab.right(), qr.left(), ab.y2, qr.y1, xinter))
			{
				return xinter;
			}
			else if (voronoi_vertex(ab, lp.right(), qr.left(), xinter))
			{
				return xinter;
			}
			else if (voronoi_vertex(lp, ab.left(), qr.left(), xinter))
			{
				return xinter;
			}
			else if (voronoi_vertex(qr, lp.right(), ab.right(), xinter))
			{
				return xinter;
			}
			else
			{
				return xinter;
			}
		}
	}
}

// -----------------------------------------------------------------------------

// Assumes that it != m_S.end()
void VoronoiMorphoF::exploreLeft(S_const_iter it, SegmentF qr, int i)
{
	while (it != m_S.begin())
	{
		double xinter = treatSegments(*std::prev(it), *it, qr);
		if (xinter <= i)
		{
			it = std::prev(m_S.erase(it));
		}
		else
		{
			if (std::ceil(xinter) < m_XMax)
			{
				m_Q[std::ceil(xinter)].push_back(*it);
			}
			return;
		}
	}
}

// Assumes that it != m_S.end()
void VoronoiMorphoF::exploreRight(S_const_iter it, SegmentF lp, int i)
{
	while (std::next(it) != m_S.end())
	{
		double xinter = treatSegments(lp, *it, *std::next(it));
		if (xinter <= i)
		{
			it = m_S.erase(it);
		}
		else
		{
			if (std::ceil(xinter) < m_XMax)
			{
				m_Q[std::ceil(xinter)].push_back(*it);
			}
			return;
		}
	}
}

// -----------------------------------------------------------------------------

// Insert a new seed segment (i,j1) -- (i,j2)
void VoronoiMorphoF::insertSegment(int i, double j1, double j2)
{
	const SegmentF s(i, j1, j2);
	// 1st step: remove all segments occluded by s
	S_const_iter it = m_S.lower_bound(s);
	while (it != m_S.end() && s.contains(*it))
	{
		it = m_S.erase(it);
	}
	while (it != m_S.begin() && s.contains(*std::prev(it)))
	{
		m_S.erase(std::prev(it));
	}
	// 2nd step: See if we have to split the neighboring segments on the right
	if (it != m_S.end())
	{
		if (it->y1 <= s.y2)
		{
			int xsup = it->x + (int)std::floor(m_Radius) + 1;
			// Special case: we may have to split the overlapping Segment
			if (it->y1 < s.y1)
			{
				SegmentF lp(it->x, it->y1, s.y1);
				m_S.insert(it, lp);
				if (xsup < m_XMax) { m_Q[xsup].push_back(lp); }
			}
			SegmentF ab(it->x, s.y2, it->y2);
			it = m_S.erase(it);
			it = m_S.insert(it, ab);
			if (xsup < m_XMax) { m_Q[xsup].push_back(ab); }
		}
	}
	// And now look at the neighboring segments on the left
	if (it != m_S.begin())
	{
		S_const_iter jt = std::prev(it);
		if (jt->y2 >= s.y1)
		{
			int xsup = jt->x + (int)std::floor(m_Radius) + 1;
			// Special case: we may have to split the overlapping Segment
			if (jt->y2 > s.y2)
			{
				SegmentF lp(jt->x, s.y2 , jt->y2);
				it = m_S.insert(it, lp);
				if (xsup < m_XMax) { m_Q[xsup].push_back(lp); }
			}
			SegmentF ab(jt->x, jt->y1, s.y1);
			jt = m_S.erase(jt);
			m_S.insert(jt, ab);
			if (xsup < m_XMax) { m_Q[xsup].push_back(ab); }
		}
	}
	// 3rd step: Insert the new segment in the list of seeds
	S_const_iter jt = m_S.insert(it, s);
	int xsup = s.x + (int)std::floor(m_Radius) + 1;
	if (xsup < m_XMax) { m_Q[xsup].push_back(s); }
	// 4th step: Inspect what's on the right and on the left
	if (it != m_S.end())
	{
		exploreRight(it, s, i);
	}
	if (jt != m_S.begin())
	{
		exploreLeft(std::prev(jt), s, i);
	}
}

// -----------------------------------------------------------------------------

// Remove seeds that are not contributing anymore to the current sweep line (x==i)
void VoronoiMorphoF::removeInactiveSegments(int i)
{
	std::sort(m_Q[i].begin(), m_Q[i].end());
	while (!m_Q[i].empty())
	{
		S_const_iter it = m_S.find(m_Q[i].back());
		m_Q[i].pop_back();
		if (it != m_S.end())
		{
			it = m_S.erase(it);
			if (it != m_S.begin() && it != m_S.end())
			{
				exploreLeft(std::prev(it), *it, i);
				exploreRight(it, *std::prev(it), i);
			}
		}
	}
	m_Q[i].clear();
}

////////////////////////////////////////////////////////////////////////////////

namespace {

	// Appends segment [a,b[ to the given line
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

} // anonymous namespace

  // -----------------------------------------------------------------------------

  // Extract the result for the current line
void VoronoiMorphoF::getLine(int i, std::vector<double> &nl, int ysize)
{
	for (Segment s : m_S)
	{
		vor_assert((i - s.x) <= m_Radius);
		const double dy = (std::sqrt(m_Radius * m_Radius - (i - s.x) * (i - s.x)));
		const double a = s.y1 - dy > 0 ? s.y1-dy : 0; // Clamp at 0 no matter what
		const double b = s.y2 + dy < ysize ? s.y2 + dy : ysize; // ysize is the actual line size
		if (a > ysize || b < 0)
			continue;			// The segment is totally beyond the range
		appendSegment(nl, a, b);
	}
}

////////////////////////////////////////////////////////////////////////////////

#ifndef WIN32
// static_assert(std::is_standard_layout<VoronoiMorpho>::value, "VoronoiMorpho is not standard layout!");
#endif
