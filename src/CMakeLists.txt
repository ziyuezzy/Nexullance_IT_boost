cmake_minimum_required(VERSION 3.10)
project(MainProject)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Include directories
include_directories(/users/ziyzhang/topology-research/nexullance/IT_boost/src/boost)

# Compiler flags
# For debugging
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g -fexceptions")
# For release
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -ffast-math -Wall -pedantic")

# Link directories and libraries
link_directories(/users/ziyzhang/topology-research/nexullance/IT_boost/src/boost/stage/lib)
set(LIBS boost_graph)

# Find all source files in the current directory
file(GLOB SRCS *.cpp)

# Create the executable target
add_executable(main.exe ${SRCS})

# Link libraries to the target
target_link_libraries(main.exe ${LIBS})
