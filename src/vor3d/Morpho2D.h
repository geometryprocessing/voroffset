#pragma once

////////////////////////////////////////////////////////////////////////////////
#include "vor3d/Common.h"
#include <iostream>
#include <vector>
#include <array>
#include <set>
////////////////////////////////////////////////////////////////////////////////


namespace voroffset3d
{
	class Morpho2D
	{
	public:
		// Constant parameters for the current dilation
		int m_XMax;
		double m_YMax, m_YMin;
		double m_Radius;
		double m_DexelSize;

	public:
		virtual ~Morpho2D() = default;

		/////////////
		// Methods //
		/////////////
	public:
		// Insert a new seed segment
		virtual void insertSegment(int i, double j1, double j2, double r) = 0;

		// Remove seeds that are not contributing anymore to the current sweep line
		virtual void removeInactiveSegments(int i) = 0;

		// Remove seeds that are not contributing anymore to the current sweep line
		virtual void removeInactivePoints(int i) {}

		virtual void flushLine(int i) {}

		// Extract the result for the current line
		virtual void getLine(bool is_set_radii, int posX, int posY, int deltaX, int deltaY, int i,
			std::function<void(int, int, double, double, double)> _appendSegment)=0;

		// Reset data
		virtual void resetData() = 0;
	};


	////////////////////////////////////////////////////////////////////////////////

} // namespace voroffset3d
