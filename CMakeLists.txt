
cmake_minimum_required(VERSION 3.20)

set(PROJECT_NAME "tesetCPP")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)
project(${PROJECT_NAME})

find_package(Eigen3 REQUIRED)

add_executable(iLQR iLQR.cpp)

add_library(iLQR_Include INTERFACE)
target_include_directories(iLQR_Include INTERFACE E ${EIGEN3_INCLUDE_DIRS} Integration Cost Dynamics Controller Approximation OptimalControl Rollout IntrusiveList AugmentedLagrangian Constraint Penalties RiccatiEquations SearchStrategy Misc Initialization ${CMAKE_SOURCE_DIR})
target_link_libraries(iLQR PUBLIC iLQR_Include)

set(CMAKE_C_FLAGS "-std=gnu11 -Wall -Wpedantic -Wextra -ffunction-sections -fdata-sections")
set(CMAKE_CXX_FLAGS "-Wall -Wpedantic -Wextra -ffunction-sections -fdata-sections")

set(CMAKE_EXE_LINKER_FLAGS "-Wl,--gc-sections -Wl,--start-group -Wl,--end-group")
set(CMAKE_C_FLAGS_DEBUG "-Og -g3 -ggdb")
set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-Og -g3 -ggdb")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

# 启用测试
enable_testing()

# 添加测试子目录
add_subdirectory(Tests)