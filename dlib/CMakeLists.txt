cmake_minimum_required(VERSION 3.14)

project(spoc_viewer)

include(FetchContent)
FetchContent_Declare(dlib
    GIT_REPOSITORY https://github.com/davisking/dlib.git
    GIT_TAG        v19.18)
FetchContent_MakeAvailable(dlib)

add_executable(spoc_viewer spoc_viewer.cpp)
target_link_libraries(spoc_viewer dlib::dlib)
