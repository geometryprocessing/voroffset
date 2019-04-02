#pragma once

////////////////////////////////////////////////////////////////////////////////
#include "vor3d/CompressedVolumeBase.h"
////////////////////////////////////////////////////////////////////////////////

namespace voroffset3d
{

	class CompressedVolume : public CompressedVolumeBase
	{
	private:
		// Member data
		std::vector<std::vector<Scalar> > m_Data;

	public:
		// Interface
		CompressedVolume(Eigen::Vector3d origin, Eigen::Vector3d extent, Scalar voxel_size, int padding);
		CompressedVolume() = default;

		// Marshalling
		void save(std::ostream &out) const;
		void load(std::istream &in);

		void copy_volume_from(const CompressedVolume voxel);
		double get_volume();
		int numSegments();
		const std::vector<Scalar> & at(int x, int y) const { return m_Data[x + m_GridSize[0] * y]; }
		std::vector<Scalar> & at(int x, int y) { return m_Data[x + m_GridSize[0] * y]; }


		// Apply a function to each segment in the structure
		virtual void iterate(int i, int j, std::function<void(Scalar, Scalar)> func) override;
		virtual void iterate(int i, int j, std::function<void(Scalar, Scalar, Scalar)> func1) override;

		virtual void appendSegment(int i, int j, Scalar begin_pt, Scalar end_pt, Scalar radius) override;

		virtual void reshape(int xsize, int ysize) override;
		virtual void resize(int xsize, int ysize) override;
		virtual void clear() override;
	};

} // namespace voroffset3d
