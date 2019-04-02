#pragma once

////////////////////////////////////////////////////////////////////////////////
#include "vor3d/Common.h"
#include <set>
#include <limits>
#include <cassert>
#include <vector>
////////////////////////////////////////////////////////////////////////////////

namespace voroffset3d
{

	class CompressedVolumeBase
	{
	protected:
		Eigen::Vector3d m_Origin;
		Eigen::Vector3d m_Extent;
		Scalar          m_Spacing; // voxel size (in mm)
		int				m_Padding;
		Eigen::Vector2i m_GridSize;

	public:
		// Interface
		virtual ~CompressedVolumeBase() = default;
		//CompressedVolumeBase() {}

		// Apply a function to each segment in the structure
		virtual void iterate(int i, int j, std::function<void(Scalar, Scalar, Scalar)> func) = 0;
		virtual void iterate(int i, int j, std::function<void(Scalar, Scalar)> func) = 0;

		// Append a segment to the given segments vector
		virtual void appendSegment(int i, int j, Scalar begin_pt, Scalar end_pt, Scalar radius) = 0;

		virtual void reshape(int xsize, int ysize) = 0;
		void reset(Eigen::Vector3d origin, Eigen::Vector3d extent, Scalar voxel_size, int padding,
						int xsize, int ysize);
		virtual void resize(int xsize, int ysize) = 0;
		virtual void clear() = 0;

		int numDexels() const { return m_GridSize[0] * m_GridSize[1]; }
		Eigen::Vector2d dexelCenter(int x, int y) const;
		Eigen::Vector2i gridSize() const { return m_GridSize; }
		Eigen::Vector3d origin() const { return m_Origin; }
		double spacing() const { return m_Spacing; }
		Eigen::Vector3d extent() const { return m_Extent; }
		int padding() const { return m_Padding; }

	};

} // namespace voroffset3d