cmake_minimum_required(VERSION 3.20)

project(linear_algebra)

add_subdirectory(measure_duration)

add_library(
    linear_algebra
    matrix.hpp
    vector.hpp
    linear_algebra.cpp
)
target_include_directories(linear_algebra PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
set_target_properties(linear_algebra PROPERTIES CXX_STANDARD 23)

add_executable(linear_algebra_test
    main.cpp)

set_target_properties(linear_algebra_test PROPERTIES CXX_STANDARD 23)
target_link_libraries(linear_algebra_test measure_duration linear_algebra)
