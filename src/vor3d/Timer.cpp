#include "vor3d/Timer.h"

double Timer::show() 
{
	typedef std::chrono::duration<double, std::milli> double_milli;
	auto t = elapsed();
	auto t_s = std::chrono::duration_cast<std::chrono::seconds>(t);
	auto t_ms = std::chrono::duration_cast<double_milli>(t - t_s);
	auto total_ms = std::chrono::duration_cast<double_milli>(t);
	return std::chrono::duration_cast<double_milli>(t).count();
}



double Timer::get() 
{
	typedef std::chrono::duration<double, std::milli> double_milli;
	auto t = elapsed();
	auto t_ms = std::chrono::duration_cast<double_milli>(t);
	return t_ms.count();
}