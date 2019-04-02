#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <vector>
#include <array>
#include <Eigen/Dense>
////////////////////////////////////////////////////////////////////////////////
#ifndef EPS
#define EPS 1e-8
#endif

namespace voroffset3d 
{
	// typedef//
	typedef double Scalar;
	typedef std::complex<Scalar> PointF;
	typedef std::vector<Scalar>::iterator data_iterator;
	///////////////////////////////////////////////////////////////////////////////////////////////////
	//Data Struct
	//////////////////////////////////////////////////////////////////////////////////////////////////
	//Point struct
	//////////////////////////////////////////////////////////////////////////////////////////////////
	struct Point	// for voronoi diagram
	{
		double x, y;
		Point(double _x, double _y)
			:x(_x), y(_y)
		{}
		Point() = default;
		bool operator <(const Point &p) const
		{
			return y < p.y || (y == p.y && x < p.x);
		}
		inline bool contains(const Point &o) const
		{
			return std::sqrt((x - o.x)*(x - o.x) + (y - o.y)*(y - o.y)) < EPS;
		}
	};

	struct PointWithRadius	// for initial data CompressedVolumeWithRadii
	{
		Scalar _pt, _radius;
		PointWithRadius() = default;
		PointWithRadius(Scalar pt, Scalar radius)
			:_pt(pt), _radius(radius)
		{}
		bool operator<(const PointWithRadius &p) const
		{
			return _pt < p._pt || (_pt == p._pt && _radius < p._radius);
		}
		inline bool contains(const PointWithRadius &o) const
		{
			return _radius - o._radius > std::abs(_pt - o._pt) - EPS;
		}
	};

	struct PointWithRadiusX		// for power diagram 
	{
		double x, y, r;
		PointWithRadiusX(double _x, double _y, double _r)
			:x(_x), y(_y), r(_r)
		{}
		PointWithRadiusX() = default;
		bool operator <(const PointWithRadiusX &p) const
		{
			return y < p.y || (y == p.y && x < p.x) || (y == p.y && x == p.x && r < p.r);
		}

		inline bool contains(const PointWithRadiusX &o) const
		{
			//return (y1 < o.y1 || std::abs(y1 - o.y1) < EPS) && (o.y2 < y2 || std::abs(o.y2 - y2));
			return r - o.r > std::sqrt(std::pow(x - o.x, 2.0) + std::pow(y - o.y, 2.0)) - EPS;
		}
	};

	/////////////////////////////////////////////////////////////////////////////////////////////////////
	//Segment struct
	////////////////////////////////////////////////////////////////////////////////////////////////////
	struct Segment	// for voronoi diagram
	{
		double x, y1, y2;

		Segment(double _x, double _y1, double _y2)
			: x(_x), y1(_y1), y2(_y2)
		{ }
		Segment() = default;
		inline Point left()  const { return Point(x, y1); }
		inline Point right() const { return Point(x, y2); }
		inline Point any()   const { return Point(x, y1); }

		bool operator <(const Segment &s) const
		{
			return (y1 + y2 < s.y1 + s.y2) || (y1 + y2 == s.y1 + s.y2 && x < s.x);
		}

		inline bool contains(const Segment &o) const
		{
			return (y1 <= o.y1) && (o.y2 <= y2);
		}

		inline bool isPoint() const
		{
			return std::abs(y1 - y2) < EPS;
		}
		inline bool isSplit(const Segment &o) const
			/* return true if we need to split the current segment when we append o, which overlaps with current segment
			Notice: the current segment and o have this below postion relationship:
			z1 <= o.z1
			*/
		{
			return x < o.x && !isPoint();
		}
		inline bool isOcclude(const Segment &o) const
		{
			return (y1 < o.y1 + EPS) && (o.y2 < y2 + EPS);
		}
	};

	struct SegmentWithRadius	// for initial data CompressedVolumeWithRadii
	{
		Scalar y1, y2, r;
		SegmentWithRadius() = default;
		SegmentWithRadius(Scalar begin_pt, Scalar end_pt, Scalar radius)
			:y1(begin_pt), y2(end_pt), r(radius)
		{}
		inline PointWithRadius left()  const { return PointWithRadius(y1, r); }
		inline PointWithRadius right() const { return PointWithRadius(y2, r); }
		inline PointWithRadius any()   const { return PointWithRadius(y1, r); }
		bool operator <(const SegmentWithRadius &s) const
		{
			return (y1 + y2 < s.y1 + s.y2) || (y1 + y2 == s.y1 + s.y2 && r < s.r);
		}
		inline bool contains(const SegmentWithRadius &o) const
		{
			//return (y1 < o.y1 || std::abs(y1 - o.y1) < EPS) && (o.y2 < y2 || std::abs(o.y2 - y2));
			return (y1 < o.y1 + EPS) && (o.y2 < y2 + EPS);
		}

		inline bool isPoint() const
		{
			return std::abs(y1 - y2) < EPS;
		}
		inline bool isSplit(const SegmentWithRadius &o) const
			/* return true if we need to split the current segment when we append o, which overlaps with current segment
			Notice: the current segment and o have this below postion relationship:
					z1 <= o.z1
			*/
		{
			return r < o.r && !isPoint();
		}
		inline bool isOcclude(const SegmentWithRadius &o) const
		{
			return (y1 < o.y1 + EPS) && (o.y2 < y2 + EPS) && (r > o.r - EPS);
		}
	};

	
	struct SegmentWithRadiusX	// for power diagram
	{
		double x, y1, y2, r;

		SegmentWithRadiusX(double _x, double _y1, double _y2, double _r)
			: x(_x), y1(_y1), y2(_y2), r(_r)
		{ }

		SegmentWithRadiusX() = default;

		inline PointWithRadiusX left()  const { return PointWithRadiusX(x, y1, r); }
		inline PointWithRadiusX right() const { return PointWithRadiusX(x, y2, r); }
		inline PointWithRadiusX any()   const { return PointWithRadiusX(x, y1, r); }

		bool operator <(const SegmentWithRadiusX &s) const
		{
			return (y1 + y2 < s.y1 + s.y2) || (y1 + y2 == s.y1 + s.y2 && x < s.x) || (y1 + y2 == s.y1 + s.y2 && x == s.x && r < s.r);
		}

		inline bool contains(const SegmentWithRadiusX &o) const
		{
			//return (y1 < o.y1 || std::abs(y1 - o.y1) < EPS) && (o.y2 < y2 || std::abs(o.y2 - y2));
			return (y1 < o.y1 + EPS) && (o.y2 < y2 + EPS);
		}

		inline bool isPoint() const
		{
			return std::abs(y1 - y2) < EPS;
		}
		 
		inline bool isSplit(const SegmentWithRadiusX &o) const
			/* return true if we need to split the current segment when we append o, which overlaps with current segment
			Notice: the current segment and o have this below postion relationship:
						z1 <= o.z1				
			*/
		{
			return (r + x < o.r + o.x) && !isPoint();
		}
		inline bool isOcclude(const SegmentWithRadiusX &o) const
		{
			return (y1 < o.y1 + EPS) && (o.y2 < y2 + EPS) && (r + x > o.r + o.x - EPS);
		}
	};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Custom assert function
[[noreturn]] void assertion_failed(
	const std::string& condition_string,
	const std::string& file, int line);

} // namespace voroffset

// Shortcut alias
namespace vor3d = voroffset3d;

////////////////////////////////////////////////////////////////////////////////

#define vor_assert_msg(x, s)									\
do																\
{																\
    if (!(x))													\
	{															\
        voroffset3d::assertion_failed((s), __FILE__, __LINE__);	\
    }															\
} while(0)

// -----------------------------------------------------------------------------

#define vor_assert(x)											\
do																\
{																\
    if(!(x))													\
	{															\
        voroffset3d::assertion_failed(#x, __FILE__, __LINE__);	\
    }															\
} while (0)

#ifdef voroffset_NO_DEBUG

#define vor_debug(x)

#else

#define vor_debug(x)											\
do																\
{																\
    if(!(x))													\
	{															\
        voroffset3d::assertion_failed(#x, __FILE__, __LINE__);	\
    }															\
} while (0)

#endif

