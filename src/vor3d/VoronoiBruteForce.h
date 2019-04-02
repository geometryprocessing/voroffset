#pragma once

#include"vor3d/Voronoi.h"

namespace voroffset3d
{
    class VoronoiMorphoBruteForce : public VoronoiMorpho
    {
		typedef std::pair<double, int> MyPoint;
	public:
		virtual void dilation(CompressedVolume input, CompressedVolume &result, double radius, double &time_1, double &time_2) override;
    private:
    	void unionMap(std::vector<MyPoint> vec, std::vector<double> &result);
    };
}