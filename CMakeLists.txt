cmake_minimum_required(VERSION 3.20)

project(linear_algebra)

add_library(
    linear_algebra
    matrix/matrix.hpp
    vector/vector.hpp
    vector_space.hpp
    vector/fixsized_vector.hpp
    matrix/fixsized_matrix.hpp
    vector/resizeable_vector.hpp
    vector/combined_reference_vector.hpp
    matrix/combined_reference_matrix.hpp
    linear_algebra.cpp
)
target_include_directories(linear_algebra PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_features(linear_algebra PUBLIC cxx_std_23)

add_subdirectory(subprojects/complex_number)
target_link_libraries(linear_algebra PUBLIC complex_number)
