cmake_minimum_required(VERSION 2.6)
project(GAP)

add_executable(GAP main.cpp)
set(CMAKE_BUILD_TYPE "Debug")
target_compile_options(GAP PRIVATE -O0 -W -Wall -std=c++17)
# set(CMAKE_BUILD_TYPE "Release")
# target_compile_options(GAP PRIVATE -O3 -W -Wall -std=c++17)
target_link_libraries(GAP PRIVATE stdc++fs -Wl,--whole-archive -lpthread -Wl,--no-whole-archive)

install(TARGETS GAP RUNTIME DESTINATION bin)
