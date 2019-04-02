#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <vector>
#include <array>
#include <Eigen/Dense>
////////////////////////////////////////////////////////////////////////////////

namespace voroffset
{

// 2D index coordinates
typedef Eigen::Vector2i Vector2i;
typedef double Scalar;
typedef std::complex<double> PointF;
typedef std::vector<PointF> Curve;
struct m_Segment 
{
	Scalar _begin_pt, _end_pt, _radius;
	m_Segment() {};
	m_Segment(Scalar begin_pt, Scalar end_pt, Scalar radius)
		:_begin_pt(begin_pt), _end_pt(end_pt), _radius(radius)
	{}
	bool operator <(const m_Segment &s) const
	{
		return (_begin_pt + _end_pt < s._begin_pt + s._end_pt) || (_begin_pt + _end_pt == s._begin_pt + s._end_pt && _radius < s._radius);
	}
};

// 2D array of pairs of integers
typedef Eigen::Matrix<Vector2i, Eigen::Dynamic, Eigen::Dynamic> MatrixXci;

// Custom assert function
[[noreturn]] void assertion_failed(
	const std::string& condition_string,
	const std::string& file, int line);

} // namespace voroffset

// Shortcut alias
namespace vor = voroffset;

////////////////////////////////////////////////////////////////////////////////

#define vor_assert_msg(x, s)									\
do																\
{																\
    if (!(x))													\
	{															\
        voroffset::assertion_failed((s), __FILE__, __LINE__);	\
    }															\
} while(0)

// -----------------------------------------------------------------------------

#define vor_assert(x)											\
do																\
{																\
    if(!(x))													\
	{															\
        voroffset::assertion_failed(#x, __FILE__, __LINE__);	\
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
        voroffset::assertion_failed(#x, __FILE__, __LINE__);	\
    }															\
} while (0)

#endif
