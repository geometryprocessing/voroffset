#include "vor3d/HalfDilationOperator.h"

namespace voroffset3d
{
	template<typename CompressedVolumeType>
	void unionMap(CompressedVolumeType voxel_1, CompressedVolumeType voxel_2, CompressedVolumeType &result)
	{
		int x_size = voxel_1.gridSize()(0);
		int y_size = voxel_1.gridSize()(1);
		for (int i = 0; i < x_size; i++)
			for (int j = 0; j < y_size; j++)
			{
				unionSegs(voxel_1.at(i, j), voxel_2.at(i, j), result.at(i, j));
			}
	}
}