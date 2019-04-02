////////////////////////////////////////////////////////////////////////////////
#include "vor2d/DistanceTransform.h"
#include <vector>
#include <cassert>
#include <limits>
////////////////////////////////////////////////////////////////////////////////

namespace 
{

	constexpr int DIST_MAX = std::numeric_limits<int>::max();

	// v++ avoiding overflows
	int safe_incr(int v) { return (v == DIST_MAX ? v : v + 1); }

} // anonymous namespace

////////////////////////////////////////////////////////////////////////////////

// Meijster, A., Roerdink, J. B., & Hesselink, W. H. (2002).
// A general algorithm for computing distance transforms in linear time.
void voroffset::computeDistanceTransformExact(const Image<int> &img, Image<int> &dt) 
{
	const int width = (int) img.width();
	const int height = (int) img.height();

	// Initialization
	dt.resize(width, height);

	/////////////
	// PHASE 1 //
	/////////////

	#pragma omp parallel for
	for (int x = 0; x < width; ++x) 
	{
		// Scan 1
		dt(x, 0) = (img(x, 0) == 0 ? 0 : DIST_MAX);
		for (int y = 1; y < height; ++y) 
		{
			dt(x, y) = (img(x, y) == 0 ? 0 : safe_incr(dt(x, y - 1)));
		}
		// Scan 2
		for (int y = height - 2; y >= 0; --y) 
		{
			if (dt(x, y + 1) < dt(x,y)) 
			{
				dt(x, y) = 1 + dt(x, y + 1);
			}
		}
	}

	/////////////
	// PHASE 2 //
	/////////////

	#pragma omp parallel for
	for (int y = 0; y < height; ++y) 
	{
		Eigen::VectorXi s(width);
		Eigen::VectorXi t(width);
		Eigen::VectorXi r(width); // Temporary row

		// See [Meijster et al. 2002] for explanations of 'f' and 'sep'
		auto f = [y, &dt] (int x, int i) 
		{
			return (dt(i, y) == DIST_MAX ? DIST_MAX : (x - i)*(x - i) + dt(i, y));
		};
		auto sep = [y, &dt] (int i, int u)
		{
			if (dt(u, y) == DIST_MAX || dt(i, y) == DIST_MAX) 
			{
				return DIST_MAX;
			} 
			else
			{
				return (u*u - i*i + dt(u, y) - dt(i, y)) / (2 * (u - i));
			}
		};

		// Scan 3
		int q = 0; s[0] = 0; t[0] = 0;
		for (int u = 1; u < width; ++u) 
		{
			vor_assert(q < width);
			while (q >= 0 && f(t[q], s[q]) > f(t[q], u)) 
			{
				--q;
			}
			if (q < 0) 
			{
				q = 0; s[0] = u;
			} 
			else 
			{
				int w = safe_incr(sep(s[q], u));
				if (w < width) 
				{
					q = q + 1; s[q] = u; t[q] = w;
				}
			}
		}

		// Scan 4
		for (int u = width - 1; u >= 0; --u) 
		{
			r[u] = f(u, s[q]);
			if (u == t[q]) { --q; }
		}
		dt.data().row(y) = r;
	}
}

////////////////////////////////////////////////////////////////////////////////

// Danielsson, P. E. (1980). Euclidean distance mapping.
void voroffset::computeDistanceTransformApprox(const Image<int> &img, Image<Vector2i> &dist) 
{
	const int xmax = (int) dist.width();
	const int ymax = (int) dist.height();

	// Initialization
	dist.resize(xmax, ymax);
	for (int y = 0; y < ymax; ++y) 
	{
		for (int x = 0; x < xmax; ++x) 
		{
			dist(x, y) = (img(x, y) == 0 ? Vector2i(0, 0) : Vector2i(xmax, ymax));
		}
	}

	////////////
	// SCAN 1 //
	////////////

	for (int y = 1; y < ymax; ++y) 
	{
		for (int x = 0; x < xmax; ++x) 
		{
			Vector2i vt = dist(x, y-1) + Vector2i(0, -1);
			if (vt.squaredNorm() < dist(x, y).squaredNorm()) { dist(x, y) = vt; }
		}
		for (int x = 1; x < xmax; ++x)
		{
			Vector2i vt = dist(x-1, y) + Vector2i(-1, 0);
			if (vt.squaredNorm() < dist(x, y).squaredNorm()) { dist(x, y) = vt; }
		}
		for (int x = (int) xmax - 2; x >= 0; --x) 
		{
			Vector2i vt = dist(x+1, y) + Vector2i(1, 0);
			if (vt.squaredNorm() < dist(x, y).squaredNorm()) { dist(x, y) = vt; }
		}
	}

	////////////
	// SCAN 2 //
	////////////

	for (int y = ymax - 2; y >= 0; --y)
	{
		for (int x = 0; x < xmax; ++x) 
		{
			Vector2i vt = dist(x, y+1) + Vector2i(0, 1);
			if (vt.squaredNorm() < dist(x, y).squaredNorm()) { dist(x, y) = vt; }
		}
		for (int x = 1; x < xmax; ++x) 
		{
			Vector2i vt = dist(x-1, y) + Vector2i(-1, 0);
			if (vt.squaredNorm() < dist(x, y).squaredNorm()) { dist(x, y) = vt; }
		}
		for (int x = xmax - 2; x >= 0; --x) 
		{
			Vector2i vt = dist(x+1, y) + Vector2i(1, 0);
			if (vt.squaredNorm() < dist(x, y).squaredNorm()) { dist(x, y) = vt; }
		}
	}
}
