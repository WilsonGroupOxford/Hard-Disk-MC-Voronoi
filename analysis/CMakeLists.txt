cmake_minimum_required(VERSION 3.7)
project(BNAHDMC_ANALYSIS)

set(CMAKE_CXX_STANDARD 11)

file(GLOB SOURCE_FILES
        "*.h"
        "*.cpp"
        "*.tpp"
        "../voro++/src/voro++.hh"
)

add_executable(analysis.x ${SOURCE_FILES})
