cmake_minimum_required(VERSION 3.14)

project(spoc_viewer)

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)

include(FetchContent)
FetchContent_Declare(dlib
    GIT_REPOSITORY https://github.com/davisking/dlib.git
    GIT_TAG        v19.18)
FetchContent_MakeAvailable(dlib)

include_directories(../../spoc/include)

add_executable(spoc_viewer spoc_viewer.cpp)

target_link_libraries(spoc_viewer dlib::dlib)
