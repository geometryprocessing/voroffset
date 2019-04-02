#pragma once

////////////////////////////////////////////////////////////////////////////////
#include "vor2d/Common.h"
#include "vor2d/Image.h"
#include <vector>
#include <array>
////////////////////////////////////////////////////////////////////////////////

namespace voroffset 
{

//
// @brief      Compute the distance transform of an image.
//
// @param[in]  img   { Input image, where 0 represents a feature. }
// @param[out] dt    { Distance to the nearest feature in the image. }
//
// The input image must contain 0 where the pixels are on the object,
// and high values (int max works) for pixels on the outside.
//
// Reference:
//
//   Meijster, A., Roerdink, J. and Hesselink, W. 2000.
//   A general algorithm for computing distance transforms in linear time.
//   In Proceedings of the 5th International Conference on Mathematical
//   Morphology and its Applications to Image and Signal Processing.
//
void computeDistanceTransformExact(const Image<int> &img, Image<int> &dt);

//
// @brief      Compute Euclidean distance map.
//
// @param[in]  img   { Input image, where 0 represents a feature. }
// @param[out] dist  { 2D array of vectors corresponding to the distance to the
//                   nearest feature. }
// Reference:
//
//   Danielsson, Per-Erik. Euclidean distance mapping.
//   Computer Graphics and image processing 14.3 (1980): 227-248.
//
void computeDistanceTransformApprox(const Image<int> &img, Image<Vector2i> &dist);

} // namespace voroffset
