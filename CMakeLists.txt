cmake_minimum_required(VERSION 3.15)
project(
  "SSE for 2d RTFIM"
  LANGUAGES CXX
)

add_compile_options(-Wall -Ofast -std=c++20 -march=native -mtune=native)
add_executable(2d 2d.cpp util.cpp)
#target_compile_options(2d PUBLIC -std=c++20 -Wall -Ofast)
