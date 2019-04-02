#pragma once

////////////////////////////////////////////////////////////////////////////////
#include "vor2d/Common.h"
#include "vor2d/Image.h"
#include <set>
#include <limits>
#include <cassert>
#include <vector>
////////////////////////////////////////////////////////////////////////////////

namespace voroffset 
{

class CompressedImage 
{
	typedef int Scalar;

private:
	//////////////////
	// Data members //
	//////////////////

	int m_XSize;
	std::vector<std::vector<Scalar>> m_Rays;

public:
	//////////////////
	// Input/Output //
	//////////////////

	CompressedImage(int w = 0, int h = 0) : m_XSize(w), m_Rays(h) { }
	CompressedImage(std::istream &in) { this->load(in); }

	// Import/export
	void fromImage(std::function<bool(int,int)> read_pixel_func);
	void toImage(std::function<void(int,int,bool)> write_pixel_func) const;

	// Marshalling
	void save(std::ostream &out) const;
	void load(std::istream &in);

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

	// Transpose operation to get a column-compressed map
	CompressedImage transposed() const;

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
};

} // namespace voroffset
