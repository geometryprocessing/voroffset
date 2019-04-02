////////////////////////////////////////////////////////////////////////////////
#include "vor2d/CompressedImage.h"
#include "vor2d/Voronoi.h"
#include <map>
#include <set>
#include <complex>
#include <algorithm>
#include <thread>
#include <cassert>
#include <type_traits>
#include <iostream>
#include <atomic>
////////////////////////////////////////////////////////////////////////////////

using namespace voroffset;

////////////////////////////////////////////////////////////////////////////////
// Input/Output
////////////////////////////////////////////////////////////////////////////////

// Import from 2D image
void voroffset::CompressedImage::fromImage(std::function<bool(int,int)> read_pixel_func)
{
	for (int j = 0; j < (int) m_Rays.size(); ++j) 
	{
		m_Rays[j].clear();
		bool prev = false;
		for (int i = 0; i < m_XSize; ++i)
		{
			bool curr = read_pixel_func(i, j);
			if (curr ^ prev) 
			{
				m_Rays[j].push_back(i);
			}
			prev = curr;
		}
		if (prev != false) 
		{
			m_Rays[j].push_back(m_XSize);
		}
	}
}

// Export to 2D image
void voroffset::CompressedImage::toImage(std::function<void(int,int,bool)> write_pixel_func) const 
{
	for (int j = 0; j < (int) m_Rays.size(); ++j) 
	{
		int i = 0;
		bool status = false;
		for (const auto & val : m_Rays[j]) 
		{
			for (; i < val; ++i) 
			{
				write_pixel_func(i, j, status);
			}
			status = !status;
		}
		for (; i < m_XSize; ++i) 
		{
			write_pixel_func(i, j, status);
		}
	}
}

#include <iomanip>

// Marshalling
void voroffset::CompressedImage::save(std::ostream &out) const 
{
	//out << std::setprecision(std::numeric_limits<double>::max_digits10);
	out << m_XSize << " "  << m_Rays.size() << "\n";
	for (const auto & row : m_Rays) 
	{
		out << row.size();
		for (const auto & val : row) 
			out << ' ' << val;
		{
		}
		out << "\n";
	}
}

void voroffset::CompressedImage::load(std::istream &in) 
{
	unsigned int num_cols;
	in >> m_XSize >> num_cols;
	m_Rays.assign(num_cols, {});
	for (auto & row : m_Rays) 
	{
		size_t size; in >> size;
		row.resize(size);
		for (auto & val : row) 
		{
			in >> val;
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
// Access
////////////////////////////////////////////////////////////////////////////////

// Resize
void voroffset::CompressedImage::resize(int w, int h) 
{
	m_XSize = w;
	m_Rays.assign(h, {});
}

// Sanity check
bool voroffset::CompressedImage::isValid() const
{
	for (const auto & row : m_Rays) 
	{
		if (row.size() % 2 != 0) 
		{
			std::cerr << "row size % 2 != 0" << std::endl;
			return false;
		}
		Scalar lastEvent = -1;
		for (const auto &val : row) 
		{
			if (val <= lastEvent) 
			{
				std::cerr << "val > lastEvent" << std::endl;
				return false;
			}
			lastEvent = val;
		}
		if (lastEvent > m_XSize) 
		{
			std::cerr << "lastEvent <= xsize" << std::endl;
			return false;
		}
	}
	return true;
}

// Empty?
bool voroffset::CompressedImage::empty() const
{
	for (const auto & row : m_Rays) 
	{
		if (!row.empty()) { return false; }
	}
	return true;
}

// Apply a function to each segment in the structure
void voroffset::CompressedImage::iterate(std::function<void(int, Scalar, Scalar)> func) const 
{
	for (int i = 0; i < (int) m_Rays.size(); ++i) 
	{
		if (m_Rays.size() % 2 != 0) 
		{
			throw std::runtime_error("Invalid number of events along a ray.");
		}
		for (size_t j = 0; 2*j < m_Rays[i].size(); ++j)
		{
			func(i, m_Rays[i][2*j], m_Rays[i][2*j+1]);
		}
	}
}

// Writes a pixel at a given location
// return true if the value was changed
bool voroffset::CompressedImage::write(int i, int j, bool newVal) 
{
	if (j < 0 || j >= (int) m_Rays.size()) 
	{
		return false;
	}
	std::vector<int> &row = m_Rays[j];
	if (row.empty()) 
	{
		// poke in an empty line
		if (!newVal) { return false; }
		row = {i, i+1};
		return true;
	}
	if (row[0] > i) 
	{
		// poke to the left of the segments
		if (!newVal) { return false; }
		if (row[0] > i+1) 
		{
			row.insert(row.begin(), 2, i+1);
		}
		row[0] = i;
		return true;
	}
	size_t a = 0;
	size_t b = row.size() - 1;
	if (row[b] <= i) 
	{
		// poke to the right of the segments
		if (!newVal) { return false; }
		if (row[b] == i) 
		{
			row[b] = i+1;
		} 
		else 
		{
			row.push_back(i);
			row.push_back(i+1);
		}
		return true;
	}
	// Assumes that row[a] ≤ i < row[b]
	while (a + 1 < b) 
	{
		size_t c = (a + b) / 2;
		if (row[c] <= i) 
		{
			a = c;
		} 
		else 
		{
			b = c;
		}
	}
	bool isInside = ((a % 2) == 0);
	if (isInside == newVal) { return false; }
	int len = row[a+1] - row[a];
	if (len == 1) 
	{ // erase all
		row.erase(row.begin()+a, row.begin()+a+2);
	}
	else if (i == row[a]) 
	{ // move left end to the right
		row[a] = i+1;
	} 
	else if (i+1 == row[a+1]) 
	{ // move right end to the left
		row[a+1] = i;
	} 
	else 
	{
		row.insert(row.begin()+a+1, 2, i+1);
		row[a+1] = i;
	}
	return true;
}

// Dichotomy search
bool voroffset::CompressedImage::at(int i, int j) const 
{
	if (j < 0 || j >= (int) m_Rays.size())
	{
		return false;
	}
	const std::vector<int> &row = m_Rays[j];
	size_t a = 0;
	size_t b = row.size() - 1;
	if (row.empty() || row[a] > i) { return false; }
	if (row[b] <= i) { return (b % 2) == 0; }
	// Assumes that row[a] ≤ i < row[b]
	while (a + 1 < b) 
	{
		size_t c = (a + b) / 2;
		if (row[c] <= i) 
		{
			a = c;
		}
		else 
		{
			b = c;
		}
	}
	return (a % 2) == 0;
}

// Returns the min k such that row[k, i] = 0 (or i+1 if not possible)
int voroffset::CompressedImage::lowerBound(int i, int j) const
{
	const int kmin = std::numeric_limits<int>::min();
	if (j < 0 || j >= (int) m_Rays.size()) 
	{
		return kmin;
	}
	const std::vector<int> &row = m_Rays[j];
	size_t a = 0;
	size_t b = row.size() - 1;
	if (row.empty() || row[a] > i) { return kmin; }
	if (row[b] <= i) { return (b % 2 == 0 ? i + 1 : row[b]); }
	// Assumes that row[a] ≤ i < row[b]
	while (a + 1 < b) 
	{
		size_t c = (a + b) / 2;
		if (row[c] <= i)
		{
			a = c;
		} 
		else 
		{
			b = c;
		}
	}
	return (a % 2 == 0 ? i + 1 : row[a]);
}

// Returns the max k such that row[i, k-1] = 0 (or i if not possible)
int voroffset::CompressedImage::upperBound(int i, int j) const 
{
	const int kmax = std::numeric_limits<int>::max();
	if (j < 0 || j >= (int) m_Rays.size()) 
	{
		return kmax;
	}
	const std::vector<int> &row = m_Rays[j];
	size_t a = 0;
	size_t b = row.size() - 1;
	if (row.empty() || row[b] <= i) { return kmax; }
	if (row[a] > i)  { return (a % 2 == 1 ? i : row[a]); }
	// Assumes that row[a] ≤ i < row[b]
	while (a + 1 < b) 
	{
		size_t c = (a + b) / 2;
		if (row[c] <= i) 
		{
			a = c;
		} 
		else 
		{
			b = c;
		}
	}
	return (b % 2 == 1 ? i : row[b]);
}

// Returns 1 if row j has a non-zero element in [u,v]
bool voroffset::CompressedImage::between(int u, int v, int j) const 
{
	if (j < 0 || j >= (int) m_Rays.size()) 
	{
		return false;
	}
	const std::vector<int> &row = m_Rays[j];
	size_t a = 0;
	size_t b = row.size() - 1;
	if (row.empty()) { return false; }
	if (row[a] > v) { return false; }
	if (row[b] <= v) { return (b % 2) == 0 || (row[b] > u); }
	// Assumes that row[a] ≤ b < row[b]
	while (a + 1 < b) 
	{
		size_t c = (a + b) / 2;
		if (row[c] <= v) 
		{
			a = c;
		} 
		else 
		{
			b = c;
		}
	}
	return (a % 2) == 0 || (row[a] > u);
}

////////////////////////////////////////////////////////////////////////////////
// Global operations
////////////////////////////////////////////////////////////////////////////////

// Negates the map
void voroffset::CompressedImage::negate()
{
	for (auto & row : m_Rays) 
	{
		if (row.empty()) 
		{
			row = {0, m_XSize};
		} 
		else 
		{
			// First event
			if (row[0] == 0) 
			{
				row.erase(row.begin());
			} 
			else 
			{
				row.insert(row.begin(), 0);
			}
			// Last event
			if (row.back() == m_XSize)
			{
				row.pop_back();
			} 
			else 
			{
				row.push_back(m_XSize);
			}
		}
	}
}

// Transpose operation to get a column-compressed map
voroffset::CompressedImage voroffset::CompressedImage::transposed() const
{
	CompressedImage tr(*this);
	tr.transposeInPlace();
	return tr;
}

void voroffset::CompressedImage::transposeInPlace()
{
	// Store all events and sort them transposed
	std::set<std::pair<int, int> > allEvents;
	for (int j = 0; j < (int) m_Rays.size(); ++j)
	{
		for (const auto & val : m_Rays[j]) 
		{
			allEvents.emplace(val, j);
		}
	}
	// Resize internal array and swap dims
	{
		int tmp = m_XSize;
		m_XSize = (int) m_Rays.size();
		m_Rays.assign(tmp, {});
	}

	// Helper function: toggle an event in a line
	auto toggle = [](std::set<int> &s, int x) 
	{
		if (s.count(x)) 
		{
			s.erase(x);
		} 
		else 
		{
			s.insert(x);
		}
	};

#if 0
	int prevLine = 0;
	std::set<int> currentLine;
	for (auto it = allEvents.begin(); it != allEvents.end(); ++it) 
	{
		int i = it->first;
		int j = it->second;

		// Register current line
		if (prevLine != i) 
		{
			for (int ii = prevLine; ii < i; ++ii) 
			{
				for (auto jt = currentLine.begin(); jt != currentLine.end(); ++jt)
				{
					m_Rays[ii].push_back(*jt);
				}
			}
		}

		// Toggle segment (i, j) -- (i, j + 1)
		toggle(currentLine, j);
		toggle(currentLine, j + 1);
		prevLine = i;
	}
	// Register last line
	m_Rays[prevLine].insert(m_Rays[prevLine].end(), currentLine.begin(), currentLine.end());
#else
	int startCol = -1;
	int endCol = -1;
	int prevLine = -1;
	std::set<int> currentLine;
	for (auto it = allEvents.begin(); it != allEvents.end(); ++it) 
	{
		int i = it->first;
		int j = it->second;

		// Register current line
		if (prevLine != i) 
		{
			// Toggle segment (i, startCol) -- (i, endCol)
			toggle(currentLine, startCol);
			toggle(currentLine, endCol);
			startCol = endCol = j;

			if (prevLine != -1) 
			{
				for (int ii = prevLine; ii < i; ++ii) 
				{
					for (auto jt = currentLine.begin(); jt != currentLine.end(); ++jt)
					{
						m_Rays[ii].push_back(*jt);
					}
				}
			}
		}

		if (j != endCol)
		{
			// Toggle segment (i, startCol) -- (i, endCol)
			toggle(currentLine, startCol);
			toggle(currentLine, endCol);
			startCol = j;
		}
		endCol = j + 1;
		prevLine = i;
	}

	// Toggle segment (i, startCol) -- (i, endCol)
	toggle(currentLine, startCol);
	toggle(currentLine, endCol);

	// Register last line
	m_Rays[prevLine].insert(m_Rays[prevLine].end(), currentLine.begin(), currentLine.end());
#endif
}

////////////////////////////////////////////////////////////////////////////////
// Morphology operators
////////////////////////////////////////////////////////////////////////////////

namespace 
{
	typedef std::vector<int> IntVector;

	template<typename Container>
	void addSegmentAtTheEnd(Container & dst, int start, int end) 
	{
		if (dst.empty() || (dst.back() < start)) 
		{
			dst.push_back(start);
			dst.push_back(end);
		} 
		else if (end > dst.back()) 
		{
			dst.back() = end;
		}
	}

	void getNextSegment(const IntVector & segs, IntVector::const_iterator & it,
		int & left, int & right, bool & full)
	{
		if (it == segs.end()) 
		{
			full = false;
		} 
		else 
		{
			left = *it++;
			right = *it++;
		}
	}

	// Computes the union of two sorted lists of segments.
	// Uses a parallel sweep of both lists, akin to merge-sort.
	void unionSegs(const IntVector & a, const IntVector & b, IntVector & result)
	{
		if (a.empty())  
		{
			result.insert(result.end(), b.begin(), b.end());
			return;
		}
		if (b.empty()) 
		{
			result.insert(result.end(), a.begin(), a.end());
			return;
		}
		IntVector::const_iterator ia(a.begin()), ib(b.begin());
		int al, ar, bl, br;
		al = *ia++; ar = *ia++;
		bl = *ib++; br = *ib++;
		bool bfull(true), afull(true);
		while (bfull || afull) 
		{
			if (bfull && bl <= al)
			{
				addSegmentAtTheEnd(result, bl, br);
				getNextSegment(b, ib, bl, br, bfull);
			} 
			else if (afull) 
			{
				addSegmentAtTheEnd(result, al, ar);
				if (ia == a.end())
				{
					afull = false;
					al = std::numeric_limits<int>::max();
				} 
				else 
				{
					al = *ia++;
					ar = *ia++;
				}
			}
		}
	}

	template<typename Iterator1, typename Iterator2>
	void unionMap(Iterator1 a, Iterator2 b, std::vector<std::vector<int> > &res) 
	{
		for (size_t i = 0; i < res.size(); ++i) 
		{
			IntVector temp;
			unionSegs(a[i], b[i], temp);
			res[i].swap(temp);
		}
	}

} // anonymous namespace

////////////////////////////////////////////////////////////////////////////////

void CompressedImage::dilate(double radius) 
{
	typedef std::reverse_iterator<decltype(m_Rays.begin())> rev_t;

	std::vector<std::vector<int> > r1, r2;
	voronoi_half_dilate(height(), width(), radius, m_Rays.begin(), r1);
	voronoi_half_dilate(height(), width(), radius, rev_t(m_Rays.end()), r2);
	unionMap(r1.cbegin(), r2.crbegin(), m_Rays);
	vor_assert(isValid() == true);
}

// -----------------------------------------------------------------------------

void CompressedImage::erode(double radius) 
{
	typedef std::reverse_iterator<decltype(m_Rays.begin())> rev_t;

	std::vector<std::vector<int> > r1, r2;
	voronoi_half_erode(height(), width(), radius, m_Rays.begin(), r1);
	voronoi_half_erode(height(), width(), radius, rev_t(m_Rays.end()), r2);
	unionMap(r1.cbegin(), r2.crbegin(), m_Rays);
	negate();
	vor_assert(isValid());
}

// -----------------------------------------------------------------------------

void CompressedImage::close(double r) 
{
	dilate(r);
	erode(r);
}

// -----------------------------------------------------------------------------

void CompressedImage::open(double r)
{
	erode(r);
	dilate(r);
}
