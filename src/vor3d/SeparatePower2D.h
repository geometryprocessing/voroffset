#pragma once

////////////////////////////////////////////////////////////////////////////////
#include "vor3d/Common.h"
#include "vor3d/Morpho2D.h"
#include <iostream>
#include <vector>
#include <array>
#include <set>
////////////////////////////////////////////////////////////////////////////////

namespace voroffset3d
{
	// the coordinate of our sweep line are range from dexelcenter(0,0) to dexelcenter(gridsize(0),0) if scanning from direction y
	// else it should be from dexelcenter(0,0) to dexelcenter(0, gridsize(1))
	struct SeparatePowerMorpho2D : public Morpho2D
	{
		////////////////
		//////////////////
		// Data members //
		//////////////////

		// Constant parameters for the current dilation
		int m_XMax;
		double m_YMax, m_YMin;
		double m_Radius;
		double m_DexelSize;
		int current_line;
		// Seeds above sorted by ascending y of their midpoint
		std::vector<SegmentWithRadiusX> m_S;
		std::vector<SegmentWithRadiusX> m_S_New;
		std::vector<SegmentWithRadiusX> m_S_Tmp;
		std::set<PointWithRadiusX> m_P;
		std::vector<PointWithRadiusX> m_P_Tmp;

		// Candidate Voronoi vertices sorted by ascending x
		std::vector<std::vector<PointWithRadiusX>> m_Q_P;

		// Typedefs
		typedef std::vector<SegmentWithRadiusX>::const_iterator S_const_iter;
		typedef std::set<PointWithRadiusX>::const_iterator P_const_iter;

		/////////////
		// Methods //
		/////////////

		// Assumes that y_p < y_q < y_r
		double rayIntersect(PointWithRadiusX p, PointWithRadiusX q, PointWithRadiusX r) const;

		// Assumes that it != m_S.end()
		void exploreLeft(P_const_iter it, PointWithRadiusX p, int i);

		// Assumes that it != m_S.end()
		void exploreRight(P_const_iter it, PointWithRadiusX p, int i);

		// Insert a new seed segment
		virtual void insertSegment(int i, double j1, double j2, double r) override;
		void insertPoint(PointWithRadiusX p);

		// Remove seeds that are not contributing anymore to the current sweep line
		virtual void removeInactiveSegments(int i) override;
		void removeInactivePoints(int i);

		virtual void resetData() override;

		void flushLine(int i);

		int locatePoint(PointWithRadiusX p);
		/*	using binary search to find the first segment in m_S which occuludes the point p
			if succeed, return the index of that segment,
			otherwise, return -1
		*/

		// Extract the result for the current line
		void getLine(bool is_set_radii, int posX, int posY, int deltaX, int deltaY, int i,
			std::function<void(int, int, double, double, double)> _appendSegment);

		SeparatePowerMorpho2D(int _xmax, double _ymin, double _ymax, double _dexel_size)
			: m_XMax(_xmax)
			, m_YMax(_ymax)
			, m_YMin(_ymin)
			, m_DexelSize(_dexel_size)
			, m_Q_P(m_XMax)
		{
			current_line = -1;
		}

		SeparatePowerMorpho2D() = default;

	private:
		std::vector<double> ray_S;
		std::vector<double> ray_P;
		std::vector<double> ray_U;

	private:
		void unionSegs();	// union segments in m_S and m_S_Tmp, at mean time, we upate m_S
		void addSegment(SegmentWithRadiusX segment); // append the segment to the end

	};


	////////////////////////////////////////////////////////////////////////////////

} // namespace voroffset3d
