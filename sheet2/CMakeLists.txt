## CMake settings
cmake_minimum_required(VERSION 3.18)

## Project settings
project(pps C CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})

## Subdirectories
add_subdirectory(src)

