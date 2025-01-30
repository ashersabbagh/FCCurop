# - Main CMake buildscript for costheta
# Comments in a CMake Script are lines begining with a '#'

# - Set CMake requirements then declare project
cmake_minimum_required(VERSION 3.2)

# - Library Sources
file(GLOB toolsSrcFiles ${CMAKE_CURRENT_SOURCE_DIR}/*.cc)

# - Build library
add_library(tools SHARED ${toolsSrcFiles})


