#pragma once

////////////////////////////////////////////////////////////////////////////////
#include "vor3d/CompressedVolumeBase.h"
////////////////////////////////////////////////////////////////////////////////

namespace voroffset3d
{

	class CompressedVolumeWithRadii : public CompressedVolumeBase
	{
	private:
		std::vector<std::vector<SegmentWithRadius>> m_Data;
		std::vector<SegmentWithRadius> m_segs_vec;

	public:
		// Interface
		CompressedVolumeWithRadii(Eigen::Vector3d origin, Eigen::Vector3d extent, double voxel_size, int padding);
		CompressedVolumeWithRadii() = default;

		void copy_volume_from(const CompressedVolumeWithRadii voxel);
		double get_volume();
		int numSegments();
		const std::vector<SegmentWithRadius> & at(int x, int y) const { return m_Data[x + m_GridSize[0] * y]; }
		std::vector<SegmentWithRadius> & at(int x, int y) { return m_Data[x + m_GridSize[0] * y]; }

		// Apply a function to each segment in the structure
		//void iterate(std::function<void(int, int, Scalar, Scalar)> func) const;
		virtual void iterate(int i, int j, std::function<void(Scalar, Scalar)> func) override;
		virtual void iterate(int i, int j, std::function<void(Scalar, Scalar, Scalar)> func1) override;

		virtual void appendSegment(int i, int j, Scalar begin_pt, Scalar end_pt, Scalar radius) override;

		virtual void reshape(int xsize, int ysize) override;
		virtual void resize(int xsize, int ysize) override;
		virtual void clear() override;

	};

} // namespace voroffset3d
