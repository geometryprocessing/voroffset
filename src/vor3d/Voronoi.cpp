#include"vor3d/Voronoi.h"
#include"vor3d/MorphologyOperators.h"
using namespace voroffset3d;

#ifndef MY_TYPE_EPS
#define MY_TYPE_EPS 1e-10
#endif
void VoronoiMorpho::erosion(CompressedVolume input, CompressedVolume &result, double radius, double &time_1, double &time_2)
{
	double z_min = input.origin()(2) / input.spacing() - 1;
	double z_max = input.origin()(2) / input.spacing() + 2 * input.padding() + input.extent()(2) / input.spacing() + 1;
	negate(input, z_min, z_max);
	//result.copy_volume_from(input);
	result.reshape(input.gridSize()(0), input.gridSize()(1));
	dilation(input,result,radius, time_1, time_2);
	negateInv(result, z_min + 1, z_max - 1);
}
void VoronoiMorpho::negate(CompressedVolume &result, double z_min, double z_max)
{
	int xsize = result.gridSize()(0);
	int ysize = result.gridSize()(1);
	// First step, add a border to the current grid
	result.resize(xsize + 2, ysize + 2);
	int read = xsize*ysize - 1;
	int write = index3dToIndex2d(xsize, ysize, xsize + 2, ysize + 2);
	int count = 0;
	for (; read >= 0; read--)
	{
		count++;
		Eigen::Vector2i pre_coor = index2dToIndex3d(read, xsize + 2, ysize + 2);
		Eigen::Vector2i cur_coor = index2dToIndex3d(write, xsize + 2, ysize + 2);
		result.at(cur_coor(0), cur_coor(1)).resize(result.at(pre_coor(0), pre_coor(1)).size());
		std::copy_n(result.at(pre_coor(0), pre_coor(1)).begin(), result.at(pre_coor(0), pre_coor(1)).size(),
			result.at(cur_coor(0), cur_coor(1)).begin());
		if (count == xsize)
		{
			count = 0;
			write -= 3;
		}
		else
			write--;
	}
	for (int i = 0; i < xsize + 2; i++)
		result.at(i, 0).clear();
	for (int j = 0; j < ysize + 2; j++)
		result.at(0, j).clear();
	for (int j = 0; j < ysize + 2; j++)
		result.at(xsize + 1, j).clear();
	// Second step, negate ray
	for (int x = 0; x < result.gridSize()(0); x++)
		for (int y = 0; y < result.gridSize()(1); y++)
		{
			negate_ray(result.at(x, y), z_min, z_max);
		}
}

void VoronoiMorpho::negateInv(CompressedVolume &result, double z_min, double z_max)
{
	int xsize = result.gridSize()(0);
	int ysize = result.gridSize()(1);
	// First step, add a border to the current grid
	int read = index3dToIndex2d(1, 1, xsize, ysize);
	int write = 0;
	int count = 0;
	for (; read <= xsize*ysize - xsize - 2; write++)
	{
		count++;
		Eigen::Vector2i pre_coor = index2dToIndex3d(read, xsize, ysize);
		Eigen::Vector2i cur_coor = index2dToIndex3d(write, xsize, ysize);
		result.at(cur_coor(0), cur_coor(1)).resize(result.at(pre_coor(0), pre_coor(1)).size());
		std::copy_n(result.at(pre_coor(0), pre_coor(1)).begin(), result.at(pre_coor(0), pre_coor(1)).size(),
			result.at(cur_coor(0), cur_coor(1)).begin());
		if (count == xsize - 2)
		{
			count = 0;
			read += 3;
		}
		else
			read++;
			
	}
	result.resize(xsize - 2, ysize - 2);
	// Second step, negate ray
	for (int x = 0; x < result.gridSize()(0); x++)
		for (int y = 0; y < result.gridSize()(1); y++)
		{
			negate_ray_range(result.at(x, y), z_min, z_max);
		}
}

double VoronoiMorpho::calculateXor(CompressedVolume voxel_1, CompressedVolume voxel_2, CompressedVolume &result)
/**
* @brief      calculate the xor between two voxels, with the assumption that these two voxels have the same gridesize
*
* @param[in]  voxel_1		{ One of voxel. }
* @param[in]  voxel_2		{ Another voxel. }
* @param[in]  result		{ the result of xor operation. }
*
* @return     { The volume of the xor voxel. }
*/
{
	int x_size = voxel_1.gridSize()(0);
	int y_size = voxel_1.gridSize()(1);
	double z_min = voxel_1.origin()(2) / voxel_1.spacing();
	double z_max = voxel_1.origin()(2) / voxel_1.spacing() + 2 * voxel_1.padding() + voxel_1.extent()(2) / voxel_1.spacing();
	result.reset(voxel_1.origin(),voxel_1.extent(),voxel_1.spacing(),voxel_1.padding(), x_size, y_size);
	for (int x = 0; x < x_size; x++)
		for (int y = 0; y < y_size; y++)
			calculate_ray_xor(voxel_1.at(x, y), voxel_2.at(x, y), result.at(x, y),z_min,z_max);
	return result.get_volume();
}

