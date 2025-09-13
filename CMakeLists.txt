cmake_minimum_required(VERSION 3.10)

project(ToCoSiAn CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(gen1 src/toy_event_generator2.cpp)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

#optional
add_compile_options(-Wall -Wextra -O2)