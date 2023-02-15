cmake_minimum_required(VERSION 3.5)
project(pybind_example)
find_package(OpenMP)
find_package(Boost 1.66 COMPONENTS container REQUIRED)
set(CMAKE_CXX_STANDARD 20)

add_compile_options(-Wall -march=native -ffast-math -fPIC -O3 -DPERFECT_MATCHING_DOUBLE)

# when installing directly in CMake, redirect its intall root to dist dir
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set (CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/dist" CACHE PATH "default install path" FORCE )
endif()

add_subdirectory(pybind11)
add_subdirectory(blossom)
add_subdirectory(src)
