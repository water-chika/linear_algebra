cmake_minimum_required(VERSION 3.20)

project(linear_algebra)

add_library(
    linear_algebra
    matrix.hpp
    vector.hpp
    vector_space.hpp
    fixsized_vector.hpp
    fixsized_matrix.hpp
    combined_reference_vector.hpp
    combined_reference_matrix.hpp
    linear_algebra.cpp
)
target_include_directories(linear_algebra PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
set_target_properties(linear_algebra PROPERTIES CXX_STANDARD 23)