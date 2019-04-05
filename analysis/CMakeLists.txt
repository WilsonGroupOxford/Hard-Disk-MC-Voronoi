cmake_minimum_required(VERSION 3.7)
project(BNAHDMC_ANALYSIS)

set(CMAKE_CXX_STANDARD 11)

file(GLOB SOURCE_FILES
        "*.h"
        "*.cpp"
        "*.tpp"
)

add_executable(analysis.x ${SOURCE_FILES})
target_link_libraries(analysis.x /Users/David/Documents/Work/DPhil/Code/development/2d_colloid_monte_carlo/voro++/src/libvoro++.a)
