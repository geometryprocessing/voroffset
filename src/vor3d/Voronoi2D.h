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
	struct VoronoiMorpho2D : public Morpho2D
	{
		////////////////
		// Data types //
		// A horizontal segment


		////////////////
		//////////////////
		// Data members //
		//////////////////

		//void normalDilation(std::vector<std::vector<double>> &results);
		//void midphaseDilation(std::vector<std::vector<double>> input, std::vector<std::vector<MySegment>> &result);

		// Constant parameters for the current dilation
		int m_XMax;
		double m_YMax, m_YMin;
		double m_Radius;
		double m_DexelSize;
		// Seeds above sorted by ascending y of their midpoint
		std::set<Segment> m_S;

		// Candidate Voronoi vertices sorted by ascending x
		std::vector<std::vector<Segment> > m_Q;

		// Typedefs
		typedef std::set<Segment>::const_iterator S_const_iter;

		/////////////
		// Methods //
		/////////////

		// Assumes that y_p < y_q < y_r
		double rayIntersect(Point p, Point q, Point r) const;

		// Assumes that y_lp < y_ab < y_qr
		double treatSegments(Segment lp, Segment ab, Segment qr) const;

		// Assumes that it != m_S.end()
		void exploreLeft(S_const_iter it, Segment qr, int i);

		// Assumes that it != m_S.end()
		void exploreRight(S_const_iter it, Segment lp, int i);

		// Insert a new seed segment [(i,j1), (i,j2)]
		virtual void insertSegment(int i, double j1, double j2, double r) override;

		// Remove seeds that are not contributing anymore to the current sweep line (y==i)
		virtual void removeInactiveSegments(int i) override;

		virtual void resetData() override;

		// Extract the result for the current line
		void getLine(bool is_set_radii, int posX, int posY, int deltaX, int deltaY, int i,
			std::function<void(int, int, double, double, double)> _appendSegment);


		VoronoiMorpho2D(int _xmax, double _ymin, double _ymax, double _radius, double _dexel_size)
			: m_XMax(_xmax)
			, m_YMax(_ymax)
			, m_YMin(_ymin)
			, m_Radius(_radius)
			, m_DexelSize(_dexel_size)
			, m_Q(m_XMax)
		{}

		VoronoiMorpho2D()
		{}

		private:
			std::vector<Segment> m_new_segs;
	};


	////////////////////////////////////////////////////////////////////////////////

} // namespace voroffset3d
