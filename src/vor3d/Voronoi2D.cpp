////////////////////////////////////////////////////////////////////////////////
#include "vor3d/Voronoi2D.h"
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
#ifndef TOLENRANCE
#define TOLENRANCE 1e-6
#endif
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

	////////////////////////////////////////////////////////////////////////////////

	//template<typename Scalar>
	struct Ray
	{
		Scalar a, b, c;

		Ray(Scalar _a, Scalar _b, Scalar _c)
			: a(_a), b(_b), c(_c)
		{ }

		// Perpendicular bisector of p and q, equation is a*x + b*y = c
		Ray(Point p, Point q)
			: a(2 * (p.x - q.x))
			, b(2 * (p.y - q.y))
			, c(p.x * p.x - q.x * q.x + p.y * p.y - q.y * q.y)
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

	////////////////////////////////////////////////////////////////////////////////

	/*
	* Assumes that y_p < y_q < y_r
	* Returns true if the intersection is between y1 and y2
	*/
	bool ray_intersect_range(Point p, Point q, Point r, double y1, double y2, double &xinter)
	{
		if (std::abs(p.x - q.x)<EPS && std::abs(q.x - r.x)<EPS)
		{
			return false;
		}
		else
		{
			const Ray r1(p, q);
			const Ray r2(q, r);
			const double denom = det(r1.a, r2.a, r1.b, r2.b);
			if (std::abs(denom) < EPS)
			{
				return false;
			}
			else
			{
				const double ynumer = det(r1.a, r2.a, r1.c, r2.c);
				if (y1 < ynumer / denom + EPS
					&& ynumer / denom < y2 + EPS)
				{
					const double xnumer = det(r1.c, r2.c, r1.b, r2.b);
					if (xnumer / denom > q.x)
					{
						xinter = ((double)xnumer / denom);
						return true;
					}
				}
				return false;
			}
		}
	}

	bool ray_intersect_after(Point p, Point q, Point r, double y, double &xinter)
	{
		if (std::abs(p.x - q.x)<EPS && std::abs(q.x - r.x)<EPS)
		{
			return false;
		}
		else
		{
			const Ray r1(p, q);
			const Ray r2(q, r);
			const double denom = det(r1.a, r2.a, r1.b, r2.b);
			if (std::abs(denom) < EPS)
			{
				return false;
			}
			else
			{
				const double ynumer = det(r1.a, r2.a, r1.c, r2.c);
				if (y < ynumer / denom + EPS)
				{
					const double xnumer = det(r1.c, r2.c, r1.b, r2.b);
					if (xnumer / denom > q.x)
					{
						xinter = ((double)xnumer / denom);
						return true;
					}
				}
				return false;
			}
		}
	}

	bool ray_intersect_before(Point p, Point q, Point r, double y, double &xinter)
	{
		if (std::abs(p.x - q.x)<EPS && std::abs(q.x - r.x)<EPS)
		{
			return false;
		}
		else
		{
			const Ray r1(p, q);
			const Ray r2(q, r);
			const double denom = det(r1.a, r2.a, r1.b, r2.b);
			if (std::abs(denom) < EPS)
			{
				return false;
			}
			else
			{
				const double ynumer = det(r1.a, r2.a, r1.c, r2.c);
				if (ynumer / denom < y + EPS)
				{
					const double xnumer = det(r1.c, r2.c, r1.b, r2.b);
					if (xnumer / denom > q.x)
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
	bool voronoi_vertex(Segment s, Point p, Point q, double &xinter)
	{
		/*
		* We seek the next voronoi point within the range [s.y1, s.y2].
		* Let (xm, ym) be its coordinates. Per the previous sentence we must have:
		*   (a) xm >= s.x
		*   (b) s.y1 ¡Ü ym ¡Ü s.y2
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
		const double u = 2.0 * (p.x - s.x);
		const double w = s.x * s.x - p.x * p.x;
		double A = r.a;
		double B = r.b * u - 2.0 * r.a * p.y;
		double C = r.a * p.y * p.y - w * r.a - r.c * u;
		//reduce(A, B, C);
		const double delta = B * B - 4.0 * A * C;
		//const ll alpha = 2 * A * A * s.y1 + 2 * A * B;
		//const ll beta = 2 * A * A * s.y2 + 2 * A * B;
		const double alpha = 2.0 * std::abs(A) * s.y1 + sign(A) * B;
		const double beta = 2.0 * std::abs(A) * s.y2 + sign(A) * B;
		const bool pos = (A > 0);
		// We must have alpha ¡Ü sign(a) * ¡Àsqrt(delta) ¡Ü beta
		if (std::abs(A)<EPS)
		{
			if (2.0 * s.y1 <= p.y + q.y && p.y + q.y <= 2.0 * s.y2)
			{
				//vor_assert(X(p) != s.x);
				const double xm = ((p.y - q.y) * (p.y - q.y) - 4.0 * w) / (4.0 * u);
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
		else if (std::abs(delta)<EPS && alpha < EPS && -EPS < beta)
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
double VoronoiMorpho2D::rayIntersect(Point p, Point q, Point r) const
{
	if (std::abs(p.x - q.x) < EPS && std::abs(q.x - r.x) < EPS)
	{
		return m_XMax*1.0f;
	}
	else
	{
		const Ray r1(p, q);
		const Ray r2(q, r);
		const double denom = det(r1.a, r2.a, r1.b, r2.b);
		if (std::abs(denom) < EPS)
		{
			// Intersection at x >= m_XMax are discarded;
			return m_XMax*1.0f;
		}
		else
		{
			const double numer = det(r1.c, r2.c, r1.b, r2.b);
			if (numer /denom < q.x + EPS)
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
double VoronoiMorpho2D::treatSegments(Segment lp, Segment ab, Segment qr) const
{
	if (ab.x + EPS > lp.x && ab.x + EPS > qr.x)
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
			else if (voronoi_vertex(lp, ab.left(), qr.any(), xinter))
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
void VoronoiMorpho2D::exploreLeft(S_const_iter it, Segment qr, int i)
{
	while (it != m_S.begin())
	{
		double xinter = treatSegments(*std::prev(it), *it, qr);
		if (xinter < i + EPS)
		{
			it = std::prev(m_S.erase(it));
		}
		else
		{
			int x_sup = (int)std::floor(xinter) + 1;
			if (xinter < m_XMax)
			{
				if (x_sup < m_XMax)
				{
					m_Q[x_sup].push_back(*it);
				}
			}
			return;
		}
	}
}

// Assumes that it != m_S.end()
void VoronoiMorpho2D::exploreRight(S_const_iter it, Segment lp, int i)
{
	while (std::next(it) != m_S.end())
	{
		double xinter = treatSegments(lp, *it, *std::next(it));
		if (xinter < i + EPS)
		{
			it = m_S.erase(it);
		}
		else
		{
			int x_sup = (int)std::floor(xinter) + 1;
			if (xinter < m_XMax)
			{
				if (x_sup < m_XMax)
				{
					m_Q[x_sup].push_back(*it);
				}
			}
			return;
		}
	}
}

// -----------------------------------------------------------------------------

// Insert a new seed segment (i,j1) -- (i,j2)
void VoronoiMorpho2D::insertSegment(int i, double j1, double j2, double r)
{
	m_new_segs.clear();
	const Segment s(i, j1, j2);
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
			int xsup = (int)std::floor(it->x + m_Radius)+1;
			// Special case: we may have to split the overlapping Segment
			if (it->y1 < s.y1)
			{
				Segment lp(it->x, it->y1, s.y1);
				//m_S.insert(it, lp);
				//if (xsup < m_XMax) { m_Q[xsup].push_back(lp); }
				m_new_segs.push_back(lp);
			}
			Segment ab(it->x, s.y2, it->y2);
			it = m_S.erase(it);
			//it = m_S.insert(it, ab);
			//if (xsup < m_XMax) { m_Q[xsup].push_back(ab); }
			m_new_segs.push_back(ab);
		}
	}
	// And now look at the neighboring segments on the left
	if (it != m_S.begin())
	{
		S_const_iter jt = std::prev(it);
		if (jt->y2 >= s.y1)
		{
			int xsup = (int)std::floor(jt->x + m_Radius)+1;
			// Special case: we may have to split the overlapping Segment
			if (jt->y2 > s.y2)
			{
				Segment lp(jt->x, s.y2, jt->y2);
				//it = m_S.insert(it, lp);
				//if (xsup < m_XMax) { m_Q[xsup].push_back(lp); }
				m_new_segs.push_back(lp);
			}
			Segment ab(jt->x, jt->y1, s.y1);
			jt = m_S.erase(jt);
			//m_S.insert(jt, ab);
			//if (xsup < m_XMax) { m_Q[xsup].push_back(ab); }
			m_new_segs.push_back(ab);
		}
	}
	for (int k = 0; k < m_new_segs.size(); k++)
	{
		m_S.insert(it, m_new_segs[k]);
		int xsup = (int)std::floor(m_new_segs[k].x + m_Radius) + 1;
		if (xsup < m_XMax) { m_Q[xsup].push_back(m_new_segs[k]); }
	}
	// 3rd step: Insert the new segment in the list of seeds
	//it = m_S.lower_bound(s);
	S_const_iter jt = m_S.insert(it, s);
	it = std::next(jt);
	int xsup = (int)std::floor(s.x + m_Radius)+1;
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

// Remove seeds that are not contributing anymore to the current sweep line (ith)
void VoronoiMorpho2D::removeInactiveSegments(int i)
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



// Extract the result for the current line
void VoronoiMorpho2D::getLine(bool is_set_radii, int posX, int posY, int deltaX, int deltaY, int i,
	std::function<void(int, int, double, double, double)> _appendSegment)
{
	if (is_set_radii)
	{
		for (Segment s : m_S)
		{
			double a, b;
			double r;
			double delta_x = i - s.x < m_Radius ? i - s.x : m_Radius;
			if (i == s.x)
			{
				r = m_Radius;
			}
			else
			{
				r = std::sqrt(m_Radius * m_Radius - delta_x*delta_x);
			}
			a = s.y1;
			b = s.y2;
			int pos_x = posX + i*deltaX;
			int pos_y = posY + i*deltaY;
			_appendSegment(pos_x, pos_y, a, b, r);
		}
	}
	else
	{
		for (Segment s : m_S)
		{
			vor_assert((i - s.x) <= m_Radius + TOLENRANCE*m_DexelSize);
			double delta_x = i - s.x < m_Radius ? i - s.x : m_Radius;
			const double dy = (std::sqrt(m_Radius * m_Radius - delta_x * delta_x));
			const double a = s.y1 - dy;
			const double b = s.y2 + dy;
			int pos_x = posX + i*deltaX;
			int pos_y = posY + i*deltaY;
			_appendSegment(pos_x, pos_y, a, b, 0.f);
		}

	}
}

void VoronoiMorpho2D::resetData()
{
	m_S.clear();
}


#ifndef WIN32
// static_assert(std::is_standard_layout<VoronoiMorpho>::value, "VoronoiMorpho is not standard layout!");
#endif
