#include "vor3d/CompressedVolumeBase.h"

using namespace voroffset3d;

Eigen::Vector2d CompressedVolumeBase::dexelCenter(int x, int y) const
{
	Eigen::Vector2d pos;
	pos[0] = (x + 0.5) * m_Spacing;
	pos[1] = (y + 0.5) * m_Spacing;
	return pos + m_Origin.head<2>();
}

void CompressedVolumeBase::reset(Eigen::Vector3d origin, Eigen::Vector3d extent, Scalar voxel_size, int padding, 
	int xsize, int ysize)
{
	m_Origin = origin;
	m_Extent = extent;
	m_Spacing = voxel_size;
	m_Padding = padding;
	reshape(xsize, ysize);
}
