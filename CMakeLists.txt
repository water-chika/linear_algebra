cmake_minimum_required(VERSION 3.20)

project(linear_algebra)

add_subdirectory(measure_duration)

add_executable(test
    main.cpp
    vector.hpp
    matrix.hpp)

set_target_properties(test PROPERTIES CXX_STANDARD 23)
target_link_libraries(test measure_duration)
