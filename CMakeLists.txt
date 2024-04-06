cmake_minimum_required(VERSION 3.20)

project(linear_algebra)

add_executable(test
    main.cpp
    vector.hpp
    matrix.hpp)

set_target_properties(test PROPERTIES CXX_STANDARD 23)