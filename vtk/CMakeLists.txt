cmake_minimum_required(VERSION 3.14)

project(spoc_viewer)

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 20)

include_directories(../../spoc)

add_executable(spoc_viewer spoc_viewer.cpp)
target_link_libraries(spoc_viewer gmp z)
