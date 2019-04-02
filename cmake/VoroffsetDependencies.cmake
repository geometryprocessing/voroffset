################################################################################
# CMake download helpers
################################################################################

# Download external dependencies
include(VoroffsetDownloadExternal)

# Sanitizers
if(VOROFFSET_WITH_SANITIZERS)
	voroffset_download_sanitizers()
	find_package(Sanitizers)
endif()

################################################################################
# Required libraries
################################################################################

# Eigen
if(NOT TARGET Eigen3::Eigen)
	add_library(voroffset_eigen INTERFACE)
	voroffset_download_eigen()
	target_include_directories(voroffset_eigen SYSTEM INTERFACE ${VOROFFSET_EXTERNAL}/eigen)
	add_library(Eigen3::Eigen ALIAS voroffset_eigen)
endif()

# CL11
if(NOT TARGET CLI11::CLI11)
	voroffset_download_cli11()
	add_subdirectory(${VOROFFSET_EXTERNAL}/cli11)
	target_compile_definitions(CLI11 INTERFACE -DCLI11_STD_OPTIONAL=0)
	target_compile_definitions(CLI11 INTERFACE -DCLI11_EXPERIMENTAL_OPTIONAL=0)
endif()

# Geogram
if(NOT TARGET geogram::geogram)
	voroffset_download_geogram()
	set(GEOGRAM_SEARCH_PATHS "${VOROFFSET_EXTERNAL}/geogram")
	include(geogram)
endif()

# json
if(NOT TARGET json::json)
	add_library(voroffset_json INTERFACE)
	voroffset_download_json()
	target_include_directories(voroffset_json SYSTEM INTERFACE ${VOROFFSET_EXTERNAL}/json)
	target_include_directories(voroffset_json SYSTEM INTERFACE ${VOROFFSET_EXTERNAL}/json/nlohmann)
	add_library(json::json ALIAS voroffset_json)
endif()

# Nanosvg
if(NOT TARGET nanosvg::nanosvg)
	voroffset_download_nanosvg()
	add_library(nanosvg_nanosvg INTERFACE)
	add_library(nanosvg::nanosvg ALIAS nanosvg_nanosvg)
	target_include_directories(nanosvg_nanosvg SYSTEM INTERFACE ${VOROFFSET_EXTERNAL}/nanosvg/src)
endif()

# TBB
if(VOROFFSET_WITH_TBB AND NOT TARGET tbb::tbb)
	set(TBB_BUILD_STATIC ON CACHE BOOL " " FORCE)
	set(TBB_BUILD_SHARED OFF CACHE BOOL " " FORCE)
	set(TBB_BUILD_TBBMALLOC OFF CACHE BOOL " " FORCE)
	set(TBB_BUILD_TBBMALLOC_PROXY OFF CACHE BOOL " " FORCE)
	set(TBB_BUILD_TESTS OFF CACHE BOOL " " FORCE)

	voroffset_download_tbb()
	add_subdirectory(${VOROFFSET_EXTERNAL}/tbb tbb)
	set_property(TARGET tbb_static tbb_def_files PROPERTY FOLDER "dependencies")
	set_target_properties(tbb_static PROPERTIES COMPILE_FLAGS "-Wno-implicit-fallthrough -Wno-missing-field-initializers -Wno-unused-parameter -Wno-keyword-macro")

	add_library(voroffset_tbb INTERFACE)
	target_include_directories(voroffset_tbb SYSTEM INTERFACE ${VOROFFSET_EXTERNAL}/tbb/include)
	target_link_libraries(voroffset_tbb INTERFACE tbb_static)
	add_library(tbb::tbb ALIAS voroffset_tbb)
endif()

# yimg
if(NOT TARGET ymg::ymg)
	voroffset_download_yimg()
	add_library(yimg_yimg STATIC ${VOROFFSET_EXTERNAL}/yimg/YImage.cpp)
	target_include_directories(yimg_yimg PUBLIC ${VOROFFSET_EXTERNAL}/yimg)
	add_library(yimg::yimg ALIAS yimg_yimg)
endif()
