#pragma once
#include <iostream>
#include <chrono>

struct Timer 
{
public:
	// Public time point typedef
	typedef std::chrono::steady_clock::time_point TimePoint;

	// Compute time difference
	template<typename clock_t>
	static auto elapsed(std::chrono::time_point<clock_t> &t1)
	{
		using namespace std::chrono;
		auto t2 = clock_t::now();
		auto time_span = t2 - t1;
		t1 = t2;
		return time_span;
	}

private:
	// Member variables
	TimePoint   m_CurrentTime;

public:
	// Get current time
	static TimePoint now() { return std::chrono::steady_clock::now(); }

	// Constructor
	Timer()
		: m_CurrentTime(Timer::now())
	{ }

	// Compute elapsed time
	auto elapsed() { return Timer::elapsed(m_CurrentTime); }

	// Print elapsed time, return elapsed time (in ms)
	double show();

	// Retrieve elapsed time (in ms)
	double get();

};