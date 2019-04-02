#include "vor3d/HalfDilationOperator.h"
#include "vor3d/MorphologyOperators.h"
namespace voroffset3d
{
	// Forward sweep for dilation operation
	void halfDilate(
		Morpho2D &vor,
		bool is_apply_power_alg,
		CompressedVolumeBase &input,
		CompressedVolumeBase &output,
		int x0, int y0, int deltaX, int deltaY)
	{
		//output.clear();
		for (int i = x0, j = y0, k = 0; i < input.gridSize()(0) && i >= 0 && j < input.gridSize()(1) && j >= 0; i += deltaX, j += deltaY, ++k)
		{
			vor.removeInactiveSegments(k);
			vor.removeInactivePoints(k);
			input.iterate(i, j, [&vor, k](double z1, double z2, double r)
			{
				vor.insertSegment(k, z1, z2, r);
			});
			vor.flushLine(i);
			vor.getLine(is_apply_power_alg, x0, y0, deltaX, deltaY, k, [&output](int posX, int posY, double z1, double z2, double r)
			{
				output.appendSegment(posX, posY, z1, z2, r);
			});
		}
		
	}

}