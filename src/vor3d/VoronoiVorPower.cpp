////////////////////////////////////////////////////////////////////////////////
#include "vor3d/VoronoiVorPower.h"
#include "vor3d/Morpho2D.h"
#include "vor3d/SeparatePower2D.h"
#include "vor3d/Voronoi2D.h"
#include "vor3d/MorphologyOperators.h"
#include "vor3d/HalfDilationOperator.h"
#include "vor3d/Timer.h"
////////////////////////////////////////////////////////////////////////////////
#ifdef USE_TBB
#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"
#include "tbb/task_scheduler_init.h"
////////////////////////////////////////////////////////////////////////////////
#ifndef GRAIN_SIZE
#define GRAIN_SIZE 10
#endif
#endif

////////////////////////////////////////////////////////////////////////////////

using namespace voroffset3d;

void VoronoiMorphoVorPower::dilation(CompressedVolume input, CompressedVolume &result, double radius, double &time_1, double &time_2)
{
	int xsize = input.gridSize()(0);
	int ysize = input.gridSize()(1);
	m_zmin = input.origin()(2) / input.spacing();
	m_zmax = input.origin()(2) / input.spacing() + 2 * input.padding() + input.extent()(2) / input.spacing();
	CompressedVolumeWithRadii output1, output2, mid_output;
	CompressedVolume output3, output4;
	output1.reshape(xsize, ysize);
	output2.reshape(xsize, ysize);
	output3.reshape(xsize, ysize);
	output4.reshape(xsize, ysize);
	mid_output.reshape(xsize, ysize);
	result.reset(input.origin(), input.extent(), input.spacing(), input.padding(), xsize, ysize);

	// 1st pass
	Timer time_pass_1;
#ifdef USE_TBB
	auto firstpass = [&](const tbb::blocked_range<uint32_t> &range)
	{
		VoronoiMorpho2D op_x(ysize, m_zmin, m_zmax, radius, input.spacing());
		for (uint32_t phaseIdx = range.begin(); phaseIdx < range.end(); ++phaseIdx)
		{
			const uint32_t x = phaseIdx;
#else
			// x-direction
			VoronoiMorpho2D op_x(ysize, m_zmin, m_zmax, radius, input.spacing());
		for (int x = 0; x < xsize; x++) {
#endif
			halfDilate(op_x, true, input, output1, x, 0, 0, +1);
			op_x.resetData();
			halfDilate(op_x, true, input, output2, x, ysize - 1, 0, -1);
			op_x.resetData();
		}
#ifdef USE_TBB
	};

	tbb::blocked_range<uint32_t> rangex(0u, (uint32_t)xsize, GRAIN_SIZE);
	tbb::parallel_for(rangex, firstpass);
#endif

	unionMap(output1, output2, mid_output);
	time_1 = time_pass_1.get();

	// 2nd pass
	Timer time_pass_2;
#ifdef USE_TBB
	auto secondpass = [&](const tbb::blocked_range<uint32_t> &range)
	{
		SeparatePowerMorpho2D op_y(xsize, m_zmin, m_zmax, input.spacing());
		for (uint32_t phaseIdy = range.begin(); phaseIdy < range.end(); ++phaseIdy)
		{
			const uint32_t y = phaseIdy;
#else
		SeparatePowerMorpho2D op_y(xsize, m_zmin, m_zmax, input.spacing());
		for(int y=0;y<ysize;y++) {
#endif
			//std::cout << y << std::endl;
			halfDilate(op_y, false, mid_output, output3, 0, y, +1, 0);
			op_y.resetData();
			halfDilate(op_y, false, mid_output, output4, xsize - 1, y, -1, 0);
			op_y.resetData();
		}
#ifdef USE_TBB
	};

	tbb::blocked_range<uint32_t> rangey(0u, (uint32_t)ysize, GRAIN_SIZE);
	tbb::parallel_for(rangey, secondpass);
#endif

	unionMap(output3, output4, result);
	time_2 = time_pass_2.get();
}



