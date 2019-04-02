#pragma once

#include"vor3d/Voronoi.h"

namespace voroffset3d
{
	class VoronoiMorphoVorPower : public VoronoiMorpho
	{
	public:
		virtual void dilation(CompressedVolume input, CompressedVolume &result, double radius, double &time_1, double &time_2) override;

	private:
		double m_zmin, m_zmax;
	};
}