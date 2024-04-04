cmake_minimum_required(VERSION 3.20)
project(core_shell_mie)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS "-O3 ${CMAKE_CXX_FLAGS}")

add_compile_definitions(_USE_MATH_DEFINES)

add_executable(core_shell_mie main.cpp mie_calc.cpp)
