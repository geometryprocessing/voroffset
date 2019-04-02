#pragma once

////////////////////////////////////////////////////////////////////////////////
#include "vor2d/Common.h"
#include <vector>
#include <array>
#include <Eigen/Dense>
////////////////////////////////////////////////////////////////////////////////

namespace voroffset
{

// Thin wrapper around Eigen::Matrix<Scalar> for representing images.
// In a matrix M, elements are indexed by [row,col], while we are used to index
// elements of an image by [x,y] coordinates. So, we will use a small proxy so
// that we don't have to swap the [x,y] coefficients when accessing an element.
//
// Also, we store images as row-major matrices, as we expect pixels to be layout
// by increasing y*w+x index.

template<typename Scalar>
class Image 
{
private:
	typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixXs;

	MatrixXs _data;

public:
	Image(int w = 0, int h = 0) : _data(h, w) { }

	void resize(int w, int h) { _data.resize(h, w); }

	auto width() const { return _data.cols(); }
	auto height() const { return _data.rows(); }

	auto operator() (int x, int y) const { return _data(y, x); }
	auto & operator() (int x, int y) { return _data(y, x); }

	// low-level data access
	const MatrixXs & data() const { return _data; }
	MatrixXs & data() { return _data; }
};

} // namespace voroffset
