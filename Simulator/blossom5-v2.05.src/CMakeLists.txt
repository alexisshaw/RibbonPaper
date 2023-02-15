cmake_minimum_required(VERSION 2.14)
project(blossom_example)
set(CMAKE_CXX_STANDARD 17)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release ... FORCE)
endif()

add_compile_options(-Wall -march=native)

add_subdirectory(blossom)

add_executable(${PROJECT_NAME} example.cpp)

target_link_libraries(${PROJECT_NAME} blossom-v)
