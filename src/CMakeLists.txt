cmake_minimum_required(VERSION 3.7)
project(HDMC)

set(CMAKE_CXX_STANDARD 11)

file(GLOB SOURCE_FILES
        "*.h"
        "*.cpp"
        "*.tpp"
)

add_executable(hdmc.x ${SOURCE_FILES})
