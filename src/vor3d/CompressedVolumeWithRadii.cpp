////////////////////////////////////////////////////////////////////////////////
#include "vor3d/CompressedVolumeWithRadii.h"
#include "vor3d/MorphologyOperators.h"
////////////////////////////////////////////////////////////////////////////////

using namespace voroffset3d;

////////////////////////////////////////////////////////////////////////////////

CompressedVolumeWithRadii::CompressedVolumeWithRadii(
	Eigen::Vector3d origin, Eigen::Vector3d extent, double voxel_size, int padding)
{
	m_Origin = origin;
	m_Extent = extent;
	m_Spacing = voxel_size;
	m_Padding = padding;
	m_Origin.array() -= padding * voxel_size;
	m_GridSize[0] =  (int)(std::ceil(extent[0] / m_Spacing) + 2 * padding);
	m_GridSize[1] =  (int)(std::ceil(extent[1] / m_Spacing) + 2 * padding);
	m_Data.assign(m_GridSize[0] * m_GridSize[1], {});
	
}
// ----------------------------------------------------------------------------
void CompressedVolumeWithRadii::iterate(int i, int j, std::function<void(Scalar, Scalar, Scalar)> func)
{
	auto &m_ray = at(i, j);
	for (const SegmentWithRadius s : m_ray)
		func(s.y1, s.y2, s.r);

}

void CompressedVolumeWithRadii::iterate(int i, int j, std::function<void(Scalar, Scalar)> func)
{
	auto &m_ray = at(i, j);
	for (const SegmentWithRadius s : m_ray)
		func(s.y1, s.y2);

}

void CompressedVolumeWithRadii::copy_volume_from(const CompressedVolumeWithRadii voxel)
{
	m_Origin = voxel.origin();
	m_Extent = voxel.extent();
	m_GridSize = voxel.gridSize();
	m_Padding = voxel.padding();
	m_Spacing = voxel.spacing();
	m_Data.assign(m_GridSize[0] * m_GridSize[1], {});
	for(int x=0;x<m_GridSize[0];x++)
		for (int y = 0; y < m_GridSize[1]; y++)
		{
			m_Data[x + y*m_GridSize[0]].resize(voxel.at(x, y).size());
			std::copy_n(voxel.at(x, y).begin(), voxel.at(x, y).size(), m_Data[x + y*m_GridSize[0]].begin());
		}

}

double CompressedVolumeWithRadii::get_volume()
{
	double vol = 0.f;
	double dexel_size = m_Spacing;
	for (int x = 0; x<m_GridSize(0); x++)
		for (int y = 0; y < m_GridSize(1); y++)
		{
			iterate(x, y, [&vol, dexel_size](Scalar begin_pt, Scalar end_pt)
			{
				vol = vol + dexel_size*dexel_size*dexel_size*(end_pt - begin_pt);
			});
		}
	return vol;
}

int CompressedVolumeWithRadii::numSegments()
{
	int num_segs=0;
	for (int x = 0; x < m_GridSize(0); x++)
		for (int y = 0; y < m_GridSize(1); y++)
			iterate(x, y, [&num_segs](Scalar begin_pt, Scalar end_pt){++num_segs; });
	return num_segs;
}

void CompressedVolumeWithRadii::appendSegment(int i, int j, Scalar begin_pt, Scalar end_pt, Scalar radius)
{
	SegmentWithRadius m_seg(begin_pt, end_pt, radius);
	vor3d::appendSegment(at(i, j), m_seg);
}


void CompressedVolumeWithRadii::reshape(int xsize, int ysize)
{
	m_GridSize << xsize, ysize;
	m_Data.assign(xsize*ysize, {});
}

void CompressedVolumeWithRadii::resize(int xsize, int ysize)
{
	m_GridSize << xsize, ysize;
	m_Data.resize(xsize*ysize);
}

void CompressedVolumeWithRadii::clear()
{
	for (int x = 0; x < m_GridSize[0]; x++)
	{
		for (int y = 0; y < m_GridSize[1]; y++)
		{
			m_Data[x + m_GridSize[0] * y].clear();
			m_Data[x + m_GridSize[0] * y].shrink_to_fit();
		}
	}
	m_Data.clear();
	m_Data.shrink_to_fit();
}