////////////////////////////////////////////////////////////////////////////////
#include "vor2d/DoubleCompressedImage.h"
#include "vor2d/DoubleVoronoi.h"
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
#ifndef MIN_SEG_SIZE
#define MIN_SEG_SIZE 1e-10
#endif

////////////////////////////////////////////////////////////////////////////////
// Input/Output
////////////////////////////////////////////////////////////////////////////////

// Import from 2D image
void voroffset::DoubleCompressedImage::fromImage(std::vector<Curve> input_curves)
{
	int size_ = input_curves.size();
	std::vector<Scalar> intersections_;
	for (int j = 0; j < (int)m_Rays.size(); ++j) 
	{
		m_Rays[j].clear();
		for (int k = 0; k < size_; k++)
		{
			scanLine(intersections_,j,input_curves[k]);
			unionIntersections(m_Rays[j], intersections_);
			intersections_.clear();
		}
	}
}


void voroffset::DoubleCompressedImage::scanLine(std::vector<Scalar> &intersections, int line_x, Curve curve)
// Scanline algorithm for the line[i]
{
	intersections.clear();
	int i, j, k;
	Scalar s;
	Scalar y_in_;
	bool status;
	int size_ = curve.size();
	for (i = 0; i<size_; i++)
		if (curve[i].real() < line_x)
		{
			if (curve[(i + 1) % size_].real() > line_x)
			{
				s = (curve[(i + 1) % size_].real() - line_x) / ( curve[(i + 1) % size_].real() - curve[i % size_].real() );
				y_in_ = s*curve[i].imag() + (1 - s)*curve[(i + 1) % size_].imag();
				intersections.push_back(y_in_);
			}
			else
			{
				j = 1;
				while (curve[(i + j) % size_].real() == line_x)
					j++;
				if (curve[(i + j) % size_].real() > line_x)
					intersections.push_back(curve[(i + j - 1) % size_].imag());
			}
		}
		else if (curve[i].real() > line_x)
		{
			if (curve[(i + 1) % size_].real() < line_x)
			{
				s = (curve[(i + 1) % size_].real() - line_x) / (curve[(i + 1) % size_].real() - curve[i].real());
				y_in_ = s*curve[i].imag() + (1 - s)*curve[(i + 1) % size_].imag();
				intersections.push_back(y_in_);
			}
			else
			{
				j = 1;
				while (curve[(i + j) % size_].real() == line_x)
					j++;
				if (curve[(i + j) % size_].real() < line_x)
					intersections.push_back(curve[(i + j - 1) % size_].imag());
			}
		}

	sort(intersections.begin(), intersections.end());

}
void voroffset::DoubleCompressedImage::unionIntersections(std::vector<Scalar> &intersections_1, std::vector<Scalar> intersections_2)
// Due to we assume two closed curves have no intersections, so the imageinary part of the this two intersections have no overlaps
{
	if (intersections_2.size() == 0)
		return;
	if (intersections_1.size()==0)
	{
		for (int i = 0; i < intersections_2.size(); i++)
			intersections_1.push_back(intersections_2[i]);
	}
	else if(intersections_1[intersections_1.size()-1]<=intersections_2[0]) // intersections_2 is totally on the right side of the intersections_1
	{
		for (int i = 0; i < intersections_2.size(); i++)
			intersections_1.push_back(intersections_2[i]);
	}
	else    // intersections_2 is totally on the left side of the intersections_1
	{
		std::vector<Scalar>::iterator it_ = intersections_1.begin();
		for (int i = intersections_2.size() - 1; i >= 0; i--)
			it_ = intersections_1.insert(it_, intersections_2[i]);
	}
}

void voroffset::DoubleCompressedImage::copyFrom(const voroffset::DoubleCompressedImage source_dexel)
{
	m_XSize = source_dexel.m_XSize;
	int size_ = source_dexel.m_Rays.size();
	m_Rays.assign(size_, std::vector<Scalar>());
	for (int i = 0; i < size_; i++)
	{
		int inter_num = source_dexel.m_Rays[i].size();
		vor_assert(inter_num % 2 == 0);
		m_Rays[i].resize(inter_num);
		std::copy_n(source_dexel.m_Rays[i].begin(), inter_num, m_Rays[i].begin());
	}
}

// Export to 2D image
/*void voroffset::CompressedImageFloat::toImage(std::function<void(int, int, bool)> write_pixel_func) const {
	for (int j = 0; j < (int)m_Rays.size(); ++j) {
		int i = 0;
		bool status = false;
		for (const auto & val : m_Rays[j]) {
			for (; i < val; ++i) {
				write_pixel_func(i, j, status);
			}
			status = !status;
		}
		for (; i < m_XSize; ++i) {
			write_pixel_func(i, j, status);
		}
	}
}
*/
// Marshalling
void voroffset::DoubleCompressedImage::save(std::ostream &out) const 
{
	out << m_XSize << " " << m_Rays.size() << "\n";
	size_t num_rows = m_Rays.size();
	for (int i = 0; i < num_rows; i++)
	{
		size_t size=m_Rays[i].size();
		out << size;
		vor_assert(size % 2 == 0);
		for (int k = 0; k < size; k+=2)
		{
			out << ' ' << m_Rays[i][k];
			out << ' ' << m_Rays[i][k + 1];
		}
		out << '\n';

	}
}

void voroffset::DoubleCompressedImage::load(std::istream &in) 
{
	unsigned int num_cols;
	in >> m_XSize >> num_cols;
	m_Rays.assign(num_cols, {});
	int flag;
	for (int i = 0; i < num_cols; i++)
	{
		size_t size;
		in >> size;
		vor_assert(size % 2 == 0);
		m_Rays[i].resize(size);
		for (int k = 0; k < size; k += 2)
		{
			in >> m_Rays[i][k]>>m_Rays[i][k+1];
		}


	}
}

////////////////////////////////////////////////////////////////////////////////
// Access
////////////////////////////////////////////////////////////////////////////////

// Resize
void voroffset::DoubleCompressedImage::resize(int w, int h) 
{
	m_XSize = w;
	m_Rays.assign(h, {});
}

// Sanity check
bool voroffset::DoubleCompressedImage::isValid() const 
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
			if (val < lastEvent) 
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
bool voroffset::DoubleCompressedImage::empty() const 
{
	for (const auto & row : m_Rays) 
	{
		if (!row.empty()) { return false; }
	}
	return true;
}

// Apply a function to each segment in the structure
void voroffset::DoubleCompressedImage::iterate(std::function<void(int, Scalar, Scalar)> func) const
{
	for (int i = 0; i < (int)m_Rays.size(); ++i) 
	{
		if (m_Rays[i].size() % 2 != 0) 
		{
			throw std::runtime_error("Invalid number of events along a ray.");
		}
		for (size_t j = 0; 2 * j < m_Rays[i].size(); ++j)
		{
			func(i, m_Rays[i][2 * j], m_Rays[i][2 * j + 1]);
		}
	}
}

// Writes a pixel at a given location
// return true if the value was changed
bool voroffset::DoubleCompressedImage::write(int i, int j, bool newVal) 
{
/*	if (j < 0 || j >= (int)m_Rays.size()) {
		return false;
	}
	std::vector<Scalar> &row = m_Rays[j];
	if (row.empty()) {
		// poke in an empty line
		if (!newVal) { return false; }
		row = { i, i + 1 };
		return true;
	}
	if (row[0] > i) {
		// poke to the left of the segments
		if (!newVal) { return false; }
		if (row[0] > i + 1) {
			row.insert(row.begin(), 2, i + 1);
		}
		row[0] = i;
		return true;
	}
	size_t a = 0;
	size_t b = row.size() - 1;
	if (row[b] <= i) {
		// poke to the right of the segments
		if (!newVal) { return false; }
		if (row[b] == i) {
			row[b] = i + 1;
		}
		else {
			row.push_back(i);
			row.push_back(i + 1);
		}
		return true;
	}
	// Assumes that row[a] �� i < row[b]
	while (a + 1 < b) {
		size_t c = (a + b) / 2;
		if (row[c] <= i) {
			a = c;
		}
		else {
			b = c;
		}
	}
	bool isInside = ((a % 2) == 0);
	if (isInside == newVal) { return false; }
	int len = row[a + 1] - row[a];
	if (len == 1) { // erase all
		row.erase(row.begin() + a, row.begin() + a + 2);
	}
	else if (i == row[a]) { // move left end to the right
		row[a] = i + 1;
	}
	else if (i + 1 == row[a + 1]) { // move right end to the left
		row[a + 1] = i;
	}
	else {
		row.insert(row.begin() + a + 1, 2, i + 1);
		row[a + 1] = i;
	}
	return true;*/
	return true;
}

// Dichotomy search
bool voroffset::DoubleCompressedImage::at(int i, int j) const 
{
	if (j < 0 || j >= (int)m_Rays.size()) 
	{
		return false;
	}
	const std::vector<Scalar> &row = m_Rays[j];
	size_t a = 0;
	size_t b = row.size() - 1;
	if (row.empty() || row[a] > i) { return false; }
	if (row[b] <= i) { return (b % 2) == 0; }
	// Assumes that row[a] �� i < row[b]
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
int voroffset::DoubleCompressedImage::lowerBound(int i, int j) const
{
	const int kmin = std::numeric_limits<int>::min();
	if (j < 0 || j >= (int)m_Rays.size()) 
	{
		return kmin;
	}
	const std::vector<Scalar> &row = m_Rays[j];
	size_t a = 0;
	size_t b = row.size() - 1;
	if (row.empty() || row[a] > i) { return kmin; }
	if (row[b] <= i) { return (b % 2 == 0 ? i + 1 : row[b]); }
	// Assumes that row[a] �� i < row[b]
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
int voroffset::DoubleCompressedImage::upperBound(int i, int j) const
{
	const int kmax = std::numeric_limits<int>::max();
	if (j < 0 || j >= (int)m_Rays.size()) 
	{
		return kmax;
	}
	const std::vector<Scalar> &row = m_Rays[j];
	size_t a = 0;
	size_t b = row.size() - 1;
	if (row.empty() || row[b] <= i) { return kmax; }
	if (row[a] > i) { return (a % 2 == 1 ? i : row[a]); }
	// Assumes that row[a] �� i < row[b]
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
bool voroffset::DoubleCompressedImage::between(int u, int v, int j) const
{
	if (j < 0 || j >= (int)m_Rays.size())
	{
		return false;
	}
	const std::vector<Scalar> &row = m_Rays[j];
	size_t a = 0;
	size_t b = row.size() - 1;
	if (row.empty()) { return false; }
	if (row[a] > v) { return false; }
	if (row[b] <= v) { return (b % 2) == 0 || (row[b] > u); }
	// Assumes that row[a] �� b < row[b]
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
void voroffset::DoubleCompressedImage::negate() 
{
	for (auto & row : m_Rays) 
	{
		if (row.empty()) 
		{
			row = { 0, m_XSize*1.0f };
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
voroffset::DoubleCompressedImage voroffset::DoubleCompressedImage::transposed() const 
{
	DoubleCompressedImage tr(*this);
	tr.transposeInPlace();
	return tr;
}

void voroffset::DoubleCompressedImage::transposeInPlace()
{
	// Store all events and sort them transposed
	std::set<std::pair<int, int> > allEvents;
	for (int j = 0; j < (int)m_Rays.size(); ++j)
	{
		for (const auto & val : m_Rays[j]) 
		{
			allEvents.emplace(val, j);
		}
	}
	// Resize internal array and swap dims
	{
		int tmp = m_XSize;
		m_XSize = (int)m_Rays.size();
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
	typedef std::vector<double> DoubleVector;

	template<typename Container>
	void addSegmentAtTheEnd(Container & dst, double start, double end)
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

	void getNextSegment(const DoubleVector & segs, DoubleVector::const_iterator & it,
		double & left, double & right, bool & full)
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
	void unionSegs(const DoubleVector & a, const DoubleVector & b, DoubleVector & result)
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
		DoubleVector::const_iterator ia(a.begin()), ib(b.begin());
		double al, ar, bl, br;
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
	void unionMap(Iterator1 a, Iterator2 b, std::vector<std::vector<double> > &res) 
	{
		for (size_t i = 0; i < res.size(); ++i) 
		{
			DoubleVector temp;
			unionSegs(a[i], b[i], temp);
			res[i].swap(temp);
		}
	}

} // anonymous namespace

  ////////////////////////////////////////////////////////////////////////////////

void DoubleCompressedImage::dilate(double radius) 
{
	typedef std::reverse_iterator<decltype(m_Rays.begin())> rev_t;

	std::vector<std::vector<double> > r1, r2;
	voronoiF_half_dilate(height(), width(), radius * m_Rays.size(), m_Rays.begin(), r1);
	voronoiF_half_dilate(height(), width(), radius * m_Rays.size(), rev_t(m_Rays.end()), r2);
	unionMap(r1.cbegin(), r2.crbegin(), m_Rays);
	vor_assert(isValid() == true);
}

// -----------------------------------------------------------------------------

void DoubleCompressedImage::erode(double radius)
{
	typedef std::reverse_iterator<decltype(m_Rays.begin())> rev_t;

	std::vector<std::vector<double> > r1, r2;
	voronoiF_half_erode(height(), width(), radius, m_Rays.begin(), r1);
	voronoiF_half_erode(height(), width(), radius, rev_t(m_Rays.end()), r2);
	unionMap(r1.cbegin(), r2.crbegin(), m_Rays);
	negate();
	vor_assert(isValid());
}

// -----------------------------------------------------------------------------

void DoubleCompressedImage::close(double r) 
{
	dilate(r);
	erode(r);
}

// -----------------------------------------------------------------------------

void DoubleCompressedImage::open(double r)
{
	erode(r);
	dilate(r);
}

void DoubleCompressedImage::negateRay(std::vector<double> &result, double y_min, double y_max)
{
	int size_ = result.size();

	if (result.empty())
	{
		result = { y_min, y_max };
	}
	else
	{
		// First event
		if (result[0] == y_min)
		{
			result.erase(result.begin());
		}
		else
		{
			result.insert(result.begin(), y_min);
		}
		// Last event
		if (result.back() == y_max)
		{
			result.pop_back();
		}
		else
		{
			result.push_back(y_max);
		}
	}
}

void DoubleCompressedImage::removePoint(const std::vector<double> segs, std::vector<double> &res)
{
	int size_ = segs.size();
	for (int i = 0; i < size_; i += 2)
	{
		if (segs[i + 1] - segs[i] < MIN_SEG_SIZE)
			continue;
		res.push_back(segs[i]);
		res.push_back(segs[i + 1]);

	}
}

double DoubleCompressedImage::area()
{
	double area = 0.f;
	for (int i = 0; i < m_Rays.size(); i++)
		for (int k = 0; k < m_Rays[i].size(); k+=2)
			area += 1.0*(m_Rays[i][k + 1] - m_Rays[i][k]);
	return area;
}