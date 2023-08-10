cmake_minimum_required(VERSION 3.22.1)
project(GAP)

add_executable(GAP main.cpp)
target_compile_options(GAP PRIVATE -std=c++20)
target_link_libraries(GAP PRIVATE /usr/local/lib/libboost_program_options.so.1.82.0 /usr/local/lib/libboost_graph.so.1.82.0)

install(TARGETS GAP RUNTIME DESTINATION bin)
