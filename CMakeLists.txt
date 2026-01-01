cmake_minimum_required(VERSION 3.20)

project(linear_algebra LANGUAGES HIP)

add_library(linear_algebra INTERFACE
    matrix/matrix.hpp
    vector/vector.hpp
    vector_space.hpp
    vector/fixsized_vector.hpp
    matrix/fixsized_matrix.hpp
    vector/resizeable_vector.hpp
    vector/combined_reference_vector.hpp
    matrix/combined_reference_matrix.hpp
)
target_include_directories(linear_algebra INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_features(linear_algebra INTERFACE cxx_std_20)

target_link_libraries(linear_algebra INTERFACE complex_number)
