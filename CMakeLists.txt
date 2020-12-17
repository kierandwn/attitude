cmake_minimum_required(VERSION 3.10)

project(attitude)

set(CMAKE_CXX_STANDARD 11)

set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
add_subdirectory("lib/googletest")

include_directories("include")
file(GLOB_RECURSE SOURCES LIST_DIRECTORIES true *.h *.cpp)

add_subdirectory("test")
