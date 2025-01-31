# - Main CMake buildscript for costheta
# Comments in a CMake Script are lines begining with a '#'

# - Set CMake requirements then declare project
cmake_minimum_required(VERSION 3.2)

# Declare project
project(anaTools VERSION 0.1)

# - When Makefiles are generated, output all command lines by default
#   Do this to begin with so it's easy to see what compiler command/flags
#   are used. This can also be done by removing the 'set' command and
#   running make as 'make VERBOSE=1'.
set(CMAKE_VERBOSE_MAKEFILE ON)

# - Don't allow C++ Compiler Vendor Extensions
set(CMAKE_CXX_EXTENSIONS OFF)

# - Set C++ Standard directly to c++17
set(CMAKE_CXX_STANDARD 17)

# - Use our standard set of flags
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Werror -Wfatal-errors -Wshadow -pedantic")

# - ROOT installation
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

set(Boost_USE_STATIC_LIBS OFF)
set(BUILD_SHARED_LIBS ON)


# - Find ROOT
find_package(ROOT COMPONENTS XMLIO MathCore MathMore Minuit)
find_package(PythonInterp 3)
find_package(PythonLibs 3)
find_package(Python COMPONENTS Development)
find_package(Boost REQUIRED COMPONENTS python)

# - Debug messages to check ROOT dependencies
message(STATUS "ROOT inc dir =${ROOT_INCLUDE_DIRS}")
message(STATUS "ROOT lib path=${ROOT_LIBRARIES}")
message(STATUS "Python executable=${Python_EXECUTABLE}")
message(STATUS "Python inc dir=${Python_INCLUDE_DIRS}")
message(STATUS "Boost libs=${Boost_LIBRARIES}")
message(STATUS "Boost python=${Boost_PYTHON37_LIBRARY}")
message(STATUS "Boost found=${Boost_FOUND}")
message(STATUS "Boost inc dir=${Boost_INCLUDE_DIRS}")


# Get all libraries to same place
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

# - Include dir
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${ROOT_INCLUDE_DIRS})
include_directories(${Python_INCLUDE_DIRS})

add_subdirectory(src)
add_subdirectory(apps)

