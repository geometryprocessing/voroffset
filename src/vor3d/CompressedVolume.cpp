////////////////////////////////////////////////////////////////////////////////
#include "vor3d/CompressedVolume.h"
#include "vor3d/MorphologyOperators.h"
#include <iostream>
////////////////////////////////////////////////////////////////////////////////

using namespace voroffset3d;

////////////////////////////////////////////////////////////////////////////////

CompressedVolume::CompressedVolume(
	Eigen::Vector3d origin, Eigen::Vector3d extent, Scalar voxel_size, int padding)
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

void CompressedVolume::iterate(int i, int j, std::function<void(Scalar, Scalar)> func)
{
	const auto &ray = at(i, j);
	for (size_t k = 0; k+1 < ray.size(); k+=2)
	{
		func(ray[k], ray[k + 1]);
	}
}

void CompressedVolume::iterate(int i, int j, std::function<void(Scalar, Scalar, Scalar)> func)
{
	const auto &ray = at(i, j);
	for (size_t k = 0; k + 1 < ray.size(); k += 2)
	{
		func(ray[k], ray[k + 1], 0.f);
	}
}

void CompressedVolume::copy_volume_from(const CompressedVolume voxel)
{
	m_Origin = voxel.origin();
	m_Extent = voxel.extent();
	m_GridSize = voxel.gridSize();
	m_Padding = voxel.padding();
	m_Spacing = voxel.spacing();
	m_Data.assign(m_GridSize[0] * m_GridSize[1], {});
	//m_Data_Points.assign(m_GridSize[0] * m_GridSize[1], {});
	for(int x=0;x<m_GridSize[0];x++)
		for (int y = 0; y < m_GridSize[1]; y++)
		{
			m_Data[x + y*m_GridSize[0]].resize(voxel.at(x, y).size());
			std::copy_n(voxel.at(x, y).begin(), voxel.at(x, y).size(), m_Data[x + y*m_GridSize[0]].begin());
		}

}

double CompressedVolume::get_volume()
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

int CompressedVolume::numSegments()
{
	int num_segs=0;
	for (int x = 0; x < m_GridSize(0); x++)
		for (int y = 0; y < m_GridSize(1); y++)
			iterate(x, y, [&num_segs](Scalar begin_pt, Scalar end_pt) {++num_segs; });
	return num_segs;
}

void CompressedVolume::appendSegment(int i, int j, Scalar begin_pt, Scalar end_pt, Scalar radius)
{
	vor3d::appendSegment(at(i, j), begin_pt, end_pt);
}

void CompressedVolume::reshape(int xsize, int ysize)
{
	m_GridSize << xsize, ysize;
	m_Data.assign(xsize*ysize, {});
}

void CompressedVolume::resize(int xsize, int ysize)
{
	m_GridSize << xsize, ysize;
	m_Data.resize(xsize*ysize);
}

void CompressedVolume::clear()
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

void CompressedVolume::load(std::istream &in)
{
	in >> m_Origin(0) >> m_Origin(1) >> m_Origin(2);
	in >> m_Extent(0) >> m_Extent(1) >> m_Extent(2);
	in >> m_GridSize(0) >> m_GridSize(1);
	in >> m_Padding;
	in >> m_Spacing;
	m_Data.assign(m_GridSize(0)*m_GridSize(1), {});
	for (auto & row : m_Data)
	{
		size_t size;
		in >> size;
		row.resize(size);
		for (auto & val : row)
		{
			in >> val;
		}
	}
}

void CompressedVolume::save(std::ostream &out) const
{
	out << m_Origin(0) << ' ' << m_Origin(1) << ' ' << m_Origin(2) << "\n";
	out << m_Extent(0) << ' ' << m_Extent(1) << ' ' << m_Extent(2) << "\n";
	out << m_GridSize(0) << ' ' << m_GridSize(1) << "\n";
	out << m_Padding << "\n";
	out << m_Spacing << "\n";
	for (const auto & row : m_Data)
	{
		out << row.size();
		for (const auto & val : row)
		{
			out << ' ' << val;
		}
		out << "\n";
	}
}
