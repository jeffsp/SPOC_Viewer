project(spoc)

cmake_minimum_required(VERSION 3.0.0)

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 20)

macro(add_app name)
    add_executable(${name} ${name}.cpp)
    target_link_libraries(${name} z)
endmacro()

add_app(test_spoc_viewer)