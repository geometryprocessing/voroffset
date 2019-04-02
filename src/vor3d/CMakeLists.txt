################################################################################

cmake_minimum_required(VERSION 3.3)
project(vor3d)

################################################################################

add_library(${PROJECT_NAME}
		Common.cpp
		Common.h
		CompressedVolume.cpp
		CompressedVolume.h
		CompressedVolumeBase.cpp
		CompressedVolumeBase.h
		CompressedVolumeWithRadii.cpp
		CompressedVolumeWithRadii.h
		Dexelize.cpp
		Dexelize.h
		HalfDilationOperator.cpp
		HalfDilationOperator.h
		HalfDilationOperator.hpp
		Morpho2D.cpp
		Morpho2D.h
		MorphologyOperators.cpp
		MorphologyOperators.h
		MorphologyOperators.hpp
		SeparatePower2D.cpp
		SeparatePower2D.h
		Timer.cpp
		Timer.h
		Voronoi.cpp
		Voronoi.h
		Voronoi2D.cpp
		Voronoi2D.h
		VoronoiBruteForce.cpp
		VoronoiBruteForce.h
		VoronoiVorPower.cpp
		VoronoiVorPower.h
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
target_link_libraries(${PROJECT_NAME} PUBLIC Eigen3::Eigen geogram::geogram)

# TBB library
if(VOROFFSET_WITH_TBB)
	target_link_libraries(${PROJECT_NAME} PUBLIC tbb::tbb)
endif()
