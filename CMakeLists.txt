cmake_minimum_required(VERSION 3.22.1)
project(GAP)

add_executable(GAP main.cpp)
# set(CMAKE_BUILD_TYPE "Debug")
# target_compile_options(GAP PRIVATE -O0 -W -Wall -std=c++20)
set(CMAKE_BUILD_TYPE "Release")
target_compile_options(GAP PRIVATE -O3 -W -Wall -std=c++20)

install(TARGETS GAP RUNTIME DESTINATION bin)
