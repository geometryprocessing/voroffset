#include "vor3d/VoronoiBruteForce.h"
#include "vor3d/MorphologyOperators.h"
#ifdef USE_TBB
#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"
#include "tbb/task_scheduler_init.h"

#ifndef GRAIN_SIZE
#define GRAIN_SIZE 10
#endif
#endif

using namespace voroffset3d;


void VoronoiMorphoBruteForce::dilation(CompressedVolume input, CompressedVolume &result, double radius, double &time_1, double &time_2)
{
	int x_size = input.gridSize()(0);
	int y_size = input.gridSize()(1);
	result.reset(input.origin(), input.extent(), input.spacing(), input.padding(), x_size, y_size);
#ifdef USE_TBB
	auto operation = [&](const tbb::blocked_range<uint32_t> &range)
	{
		std::vector<MyPoint> tmp_ray;
		for (uint32_t phaseIdx = range.begin(); phaseIdx < range.end(); ++phaseIdx)
		{

			int x = phaseIdx % x_size;
			int y = phaseIdx / x_size;
#else
			std::vector<MyPoint> tmp_ray;
            for(int i = 0;i < x_size * y_size;i++)
            {
			int x = i % x_size;
			int y = i / x_size;
#endif

            tmp_ray.clear();
            for (int range_y = (int) std::ceil(y - radius);
                 range_y <= (int) std::floor(y + radius) && range_y <= y_size - 1; range_y++)
            {
                for (int range_x = (int) std::ceil(x - radius);
                     range_x <= (int) std::floor(x + radius) && range_x <= x_size - 1; range_x++)
                {
                    //overflow judgement
                    range_x = range_x > 0 ? range_x : 0;
                    range_y = range_y > 0 ? range_y : 0;
                    if (pow(range_x - x, 2.0) + pow(range_y - y, 2.0) <= radius * radius)
                    {
                        double dz = std::sqrt(radius * radius - pow(range_x - x, 2.0) - pow(range_y - y, 2.0));
                        for (int k = 0; k < input.at(range_x, range_y).size(); k += 2)
                        {
                            double down_z = input.at(range_x, range_y)[k] - dz;
                            double up_z = input.at(range_x, range_y)[k + 1] + dz;
                            tmp_ray.push_back(MyPoint(down_z, -1));
                            tmp_ray.push_back(MyPoint(up_z, 1));
                        }
                    }
                }
            }
			unionMap(tmp_ray, result.at(x, y));
		}
#ifdef USE_TBB
	};
	tbb::blocked_range<uint32_t> range(0u, (uint32_t)x_size * y_size, GRAIN_SIZE*x_size);
	tbb::parallel_for(range, operation);
#endif
}


void VoronoiMorphoBruteForce::unionMap(std::vector<MyPoint> vec, std::vector<double> &result)
{
	result.clear();
	size_t inter_num = vec.size();
	if(inter_num == 0)
		return;		// do nothing
	vor_assert(inter_num % 2 == 0);
	std::sort(vec.begin(), vec.end(), [](auto &left, auto &right)
	{
		return left.first < right.first;
	});
	int flag_ = vec[0].second;
	result.push_back(vec[0].first);
	for (int k = 1; k < inter_num; )
	{
		while (flag_ != 0)
		{
			flag_ += vec[k].second;
			k++;
		}
		result.push_back(vec[k - 1].first);
		if (k < inter_num)
		{
			result.push_back(vec[k].first);
			flag_ = vec[k].second;
			k++;
		}
	}

}

