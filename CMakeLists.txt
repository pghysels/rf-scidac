cmake_minimum_required(VERSION 3.17)
project(STRUMPACK-RF-SciDAC VERSION 0.1 LANGUAGES CXX C Fortran)

find_package(STRUMPACK REQUIRED)
add_executable(main main.cpp)
target_link_libraries(main PRIVATE STRUMPACK::strumpack)
