#pragma once

////////////////////////////////////////////////////////////////////////////////
#include "vor2d/Common.h"
#include "vor2d/DoubleCompressedImage.h"
#include <set>
#include <limits>
#include <cassert>
#include <vector>
////////////////////////////////////////////////////////////////////////////////

namespace voroffset
{
	DoubleCompressedImage create_dexels(const std::string &filename);

	void dexel_dump(const std::string &filename, const DoubleCompressedImage &dexels);

} // namespace voroffset3d
