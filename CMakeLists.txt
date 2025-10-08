cmake_minimum_required(VERSION 3.23.0)
find_package(IntelDPCPP REQUIRED)
set(CMAKE_CXX_COMPILER mpiicpx)
project(SPIC_1D)

# Set C++ version
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)


# Build type and compiler flags -O3 -qopenmp -qmkl -xHost -fast -ipo -Wall
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -qopenmp -qmkl -fast -xHost -ipo -Wall")

# Include directories
include_directories(include)

# Gather all source files
file(GLOB_RECURSE SRC_FILES src/*.cpp)

# Create the executable
add_executable(pic1d ${SRC_FILES})
