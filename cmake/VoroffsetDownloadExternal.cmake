################################################################################
include(DownloadProject)

# With CMake 3.8 and above, we can hide warnings about git being in a
# detached head by passing an extra GIT_CONFIG option
if(NOT (${CMAKE_VERSION} VERSION_LESS "3.8.0"))
	set(VOROFFSET_EXTRA_OPTIONS "GIT_CONFIG advice.detachedHead=false")
else()
	set(VOROFFSET_EXTRA_OPTIONS "")
endif()

# Shortcut function
function(voroffset_download_project name)
	download_project(
		PROJ         ${name}
		SOURCE_DIR   ${VOROFFSET_EXTERNAL}/${name}
		DOWNLOAD_DIR ${VOROFFSET_EXTERNAL}/.cache/${name}
		QUIET
		${VOROFFSET_EXTRA_OPTIONS}
		${ARGN}
	)
endfunction()

################################################################################

## Eigen
function(voroffset_download_eigen)
	voroffset_download_project(eigen
		URL     http://bitbucket.org/eigen/eigen/get/3.3.7.tar.gz
		URL_MD5 f2a417d083fe8ca4b8ed2bc613d20f07
	)
endfunction()

## CppNumericalSolvers
function(voroffset_download_cppoptlib)
	voroffset_download_project(CppNumericalSolvers
		GIT_REPOSITORY https://github.com/PatWie/CppNumericalSolvers.git
		GIT_TAG        7eddf28fa5a8872a956d3c8666055cac2f5a535d
	)
endfunction()

## CLI11
function(voroffset_download_cli11)
	voroffset_download_project(cli11
		URL     https://github.com/CLIUtils/CLI11/archive/v1.7.1.tar.gz
	)
endfunction()

## nanosvg
function(voroffset_download_nanosvg)
	voroffset_download_project(nanosvg
		GIT_REPOSITORY https://github.com/memononen/nanosvg.git
		GIT_TAG        2b08deeb553c723d151f908d786c64136d26d576
	)
endfunction()

## Sanitizers
function(voroffset_download_sanitizers)
	voroffset_download_project(sanitizers-cmake
		GIT_REPOSITORY https://github.com/arsenm/sanitizers-cmake.git
		GIT_TAG        99e159ec9bc8dd362b08d18436bd40ff0648417b
	)
endfunction()

## Json
function(voroffset_download_json)
	voroffset_download_project(json
		URL      https://github.com/nlohmann/json/releases/download/v3.1.2/include.zip
		URL_HASH SHA256=495362ee1b9d03d9526ba9ccf1b4a9c37691abe3a642ddbced13e5778c16660c
	)
endfunction()

## geogram
function(voroffset_download_geogram)
	voroffset_download_project(geogram
		GIT_REPOSITORY https://github.com/alicevision/geogram.git
		GIT_TAG        v1.6.9
	)
endfunction()

# yimg
function(voroffset_download_yimg)
	voroffset_download_project(yimg
		GIT_REPOSITORY https://github.com/jdumas/yimg.git
		GIT_TAG        6d924fdb038c86242d272d69d5e4525086fe60cc
	)
endfunction()

## TBB
function(voroffset_download_tbb)
    voroffset_download_project(tbb
        GIT_REPOSITORY https://github.com/wjakob/tbb.git
        GIT_TAG        4c3ffe5a5f37addef0dd6283c74c4402a3b4ebc9
    )
endfunction()
