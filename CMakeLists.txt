################################################################################
# General Information
################################################################################

cmake_minimum_required(VERSION 3.8)
project(voroffset)

################################################################################

set(VOROFFSET_EXTERNAL ${CMAKE_CURRENT_SOURCE_DIR}/external/)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

if(NOT CMAKE_BUILD_TYPE)
	message(STATUS "No build type selected, default to Release")
	set(CMAKE_BUILD_TYPE "Release")
endif()

################################################################################

# Use folder in Visual Studio
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Color output
include(UseColors)

# Export compile flags(used for autocompletion of the C++ code)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

# CMake plugin for vscode
include(CMakeToolsHelpers OPTIONAL)

# Generate position independent code
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Build shared or static?
# set(BUILD_SHARED_LIBRARIES OFF)

################################################################################

# CMake options
option(VOROFFSET_WITH_SANITIZERS "Enable sanitizers"  OFF)
option(VOROFFSET_WITH_TBB        "Enable TBB"         ON)

option(SANITIZE_ADDRESS          "Sanitize Address"   OFF)
option(SANITIZE_MEMORY           "Sanitize Memory"    OFF)
option(SANITIZE_THREAD           "Sanitize Thread"    OFF)
option(SANITIZE_UNDEFINED        "Sanitize Undefined" OFF)

# Override cached options
# set(SANITIZE_ADDRESS          OFF CACHE BOOL "" FORCE)
# set(VOROFFSET_WITH_TBB        ON  CACHE BOOL "" FORCE)

################################################################################
# Dependencies
################################################################################

# Sanitizers
if(VOROFFSET_WITH_SANITIZERS)
	list(APPEND CMAKE_MODULE_PATH ${VOROFFSET_EXTERNAL}/sanitizers-cmake/cmake)
endif()

include(VoroffsetDependencies)

################################################################################

# 2D version prototype
add_subdirectory(src/vor2d)

# 3D version prototype
add_subdirectory(src/vor3d)

# Binary executables
add_subdirectory(app)
