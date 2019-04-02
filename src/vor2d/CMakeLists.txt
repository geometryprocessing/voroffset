################################################################################

cmake_minimum_required(VERSION 3.3)
project(vor2d)

################################################################################

add_library(${PROJECT_NAME}
		Common.cpp
		Common.h
		CompressedImage.cpp
		CompressedImage.h
		DistanceTransform.cpp
		DistanceTransform.h
		Dexelize.cpp
		Dexelize.h
		DoubleCompressedImage.cpp
		DoubleCompressedImage.h
		DoubleVoronoi.cpp
		DoubleVoronoi.h
		Image.h
		MorphologyOperators.cpp
		MorphologyOperators.h
		Voronoi.cpp
		Voronoi.h
)

target_include_directories(${PROJECT_NAME} PUBLIC ..)

################################################################################

# Let's get a little bit paranoid
include(SetWarnings)
target_compile_options(${PROJECT_NAME} PRIVATE ${ALL_WARNINGS})

# Use C++14
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_14)

# Generate position independent code
set_target_properties(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)

# Sanitizers
if(VOROFFSET_WITH_SANITIZERS)
	add_sanitizers(${PROJECT_NAME})
endif()

# Dependencies
target_link_libraries(${PROJECT_NAME}
	PUBLIC
		Eigen3::Eigen
		geogram::geogram
		nanosvg::nanosvg
)
