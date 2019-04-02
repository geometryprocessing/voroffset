################################################################################

cmake_minimum_required(VERSION 3.3)
project(offset3d)

################################################################################

# Add executable
add_executable(${PROJECT_NAME} offset3d.cpp)

# Let's get a little bit paranoid
include(SetWarnings)
target_compile_options(${PROJECT_NAME} PRIVATE ${ALL_WARNINGS})

# Sanitizers
if(VOROFFSET_WITH_SANITIZERS)
	add_sanitizers(${PROJECT_NAME})
endif()

# Dependencies
target_link_libraries(${PROJECT_NAME} PRIVATE vor3d CLI11::CLI11 json::json)

# Output directory for binaries
set_target_properties(${PROJECT_NAME}
	PROPERTIES
	RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
)
