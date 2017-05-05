cmake_minimum_required(VERSION 3.7)
project(2Dwavecouple_C0)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-O2")
set(SOURCE_FILES main.cpp spline.cpp spline.h)
add_executable(2Dwavecouple_C0 ${SOURCE_FILES})