cmake_minimum_required(VERSION 3.16)
project(untitled2)

set(CMAKE_CXX_STANDARD 17)
add_compile_options(-Wall -Wextra -pedantic)

add_library(big BigInt/BigUInt.cpp BigInt/BigUInt.h Time/timer.cpp Time/timer.h BigInt/BigUIntBase.cpp BigInt/BigUIntBase.h BigInt/BigInt.cpp BigInt/BigInt.h)
add_library(bench Time/BenchMark.cpp Time/BenchMark.h Time/timer.cpp Time/timer.h)

link_libraries(bench big gmp gmpxx)

add_executable(untitled2 main.cpp)

add_subdirectory(Catch_tests)