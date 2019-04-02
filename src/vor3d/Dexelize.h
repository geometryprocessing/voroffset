#pragma once

////////////////////////////////////////////////////////////////////////////////
#include "vor3d/Common.h"
#include "vor3d/CompressedVolume.h"
#include <set>
#include <limits>
#include <cassert>
#include <vector>
////////////////////////////////////////////////////////////////////////////////

namespace voroffset3d
{

/**
 * @brief         Creates dexels from a triangle mesh.
 *
 * @param[in]     filename    { Filename of the mesh to load. }
 * @param[in,out] voxel_size  { Voxel size for the 2D dexel grid. }
 * @param[in]     padding     { Additional padding on each side. }
 * @param[in]     num_voxels  { Explicitly set the 2D grid size (max length), before padding. }
 *
 * @return        { The dexelized volume. }
 */
CompressedVolume create_dexels(const std::string &filename,
	double &voxel_size, int padding = 0, int num_voxels = -1);

void dexel_dump(const std::string &filename, const CompressedVolume &voxels);

} // namespace voroffset3d
