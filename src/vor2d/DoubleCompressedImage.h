#pragma once

////////////////////////////////////////////////////////////////////////////////
#include "vor2d/Common.h"
#include "vor2d/Image.h"
#include <set>
#include <limits>
#include <cassert>
#include <vector>
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
namespace voroffset
{

class DoubleCompressedImage
{


private:
	//////////////////
	// Data members //
	//////////////////
	int m_XSize;
public:
	std::vector<std::vector<Scalar>> m_Rays;
	std::vector<std::vector<m_Segment>> m_Segs;

public:
	//////////////////
	// Input/Output //
	//////////////////

	DoubleCompressedImage(int w = 0, int h = 0) : m_XSize(w), m_Rays(h) { }
	DoubleCompressedImage(std::istream &in) { this->load(in); }


	// Import/export
	void fromImage(std::vector<Curve>input_curves);
	// Marshalling
	void save(std::ostream &out) const;
	void load(std::istream &in);

private:
	void scanLine(std::vector<Scalar> &intersections, int i, Curve curve);	// Scanline algorithm for the line[i]
	void unionIntersections(std::vector<Scalar> &intersections_1, std::vector<Scalar> intersections_2);
																			// Union intersections, the results store in the first one
public:
	////////////
	// Access //
	////////////
	// Accessors
	int width() const { return m_XSize; }
	int height() const { return (int) m_Rays.size(); }

	// Resize
	void resize(int w, int h);

	// Sanity check
	bool isValid() const;

	// Empty?
	bool empty() const;

	// Apply a function to each segment in the structure
	void iterate(std::function<void(int, Scalar, Scalar)> func) const;

	const std::vector<Scalar> & at(int x) const { return m_Rays[x]; }

	// Writes a pixel at a given location
	// return true if the value was changed
	bool write(int i, int j, bool newVal);

	// Dichotomy search
	bool at(int i, int j) const;

	// Returns the min k such that x[k, i] = 0 (or i+1 if not possible)
	int lowerBound(int i, int j) const;

	// Returns the max k such that x[i, k-1] = 0 (or i if not possible)
	int upperBound(int i, int j) const;

	// Returns 1 if column j has a non-zero element between u and v (inclusive)
	bool between(int u, int v, int j) const;

public:
	///////////////////////
	// Global operations //
	///////////////////////

	// Negates the map
	void negate();

	// Copy from another dexel
	void copyFrom(const DoubleCompressedImage source_dexel);

	// Transpose operation to get a column-compressed map
	DoubleCompressedImage transposed() const;

	// Transpose operation to get a column-compressed map
	void transposeInPlace();

public:
	//////////////////////////
	// Morphology operators //
	//////////////////////////

	void dilate(double r);
	void erode(double r);
	void close(double r);
	void open(double r);
	void negateRay(std::vector<double> &result, double y_min, double y_max);
	double area();
	void removePoint(const std::vector<double> segs, std::vector<double> &res);
};

} // namespace voroffset
