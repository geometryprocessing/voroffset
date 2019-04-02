#pragma once
#include "vor3d/Common.h"
#include "vor3d/CompressedVolumeBase.h"
#include "vor3d/CompressedVolume.h"
#include "vor3d/CompressedVolumeWithRadii.h"
#include <iostream>
#include <vector>
#include <array>
#include <set>

namespace voroffset3d
{
	class VoronoiMorpho
	{

	public:
		virtual ~VoronoiMorpho() = default;
		virtual void dilation(CompressedVolume input, CompressedVolume &result, double radius, double &time_1, double &time_2) = 0;
		/**
		* @brief      dilate the input with given radius
		*
		* @param[in]  input			{ Input data. }
		* @param[in]  result		{ Output data. }
		* @param[in]  radius		{ Dilated radius. }
		* @param[in]  time_1		{ time cost by first pass}
		* @param[in]  time_2		{ time cost by second pass}
		* 
		**/
		virtual void erosion(CompressedVolume input, CompressedVolume &result, double radius, double &time_1, double &time_2);
		double calculateXor(CompressedVolume voxel_1, CompressedVolume voxel_2, CompressedVolume &result);
		/**
		* @brief      calculate the xor between two voxels, with the assumption that these two voxels have the same gridesize
		*
		* @param[in]  voxel_1		{ One of voxel. }
		* @param[in]  voxel_2		{ Another voxel. }
		* @param[in]  result		{ the result of xor operation. }
		*
		* @return     { The volume of the xor voxel. }
		*/

	protected:
		void negate(CompressedVolume &result, double z_min, double z_max);
		void negateInv(CompressedVolume &result, double z_min, double z_max);
	};
}
	