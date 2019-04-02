#pragma once
#include "vor3d/CompressedVolumeBase.h"
#include "vor3d/CompressedVolume.h"
#include "vor3d/CompressedVolumeWithRadii.h"
#include "vor3d/Morpho2D.h"

namespace voroffset3d
{
	

	/**
	* @brief      Forward sweep for dilation operation.
	*
	* @param[in]  vor					{ The algorithm for computing the dilation. }
	* @param[in]  is_apply_power_alg	{ Whether we apply power algorithm, which means that we need to change the dilation radius. }
	* @param[in]  input					{ Input data. }
	* @param[in]  output				{ Output data. }
	* @param[in]  (x0,y0,deltaX, deltaY)
										{ Arguments for locate the plane as well as decide the sweepline direction.  }
	**/
	void halfDilate(
		Morpho2D &vor,
		bool is_apply_power_alg,
		CompressedVolumeBase &input,
		CompressedVolumeBase &output,
		int x0, int y0, int deltaX, int deltaY);

	template<typename CompressedVolumeType>
	void unionMap(CompressedVolumeType voxel_1, CompressedVolumeType voxel_2, CompressedVolumeType &result);
}

#include "vor3d/HalfDilationOperator.hpp"