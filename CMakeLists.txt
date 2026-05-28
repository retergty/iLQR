
cmake_minimum_required(VERSION 3.20)

set(PROJECT_NAME "iLQR")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)
project(${PROJECT_NAME} LANGUAGES CXX)

add_library(iLQR INTERFACE)
add_library(iLQR::iLQR ALIAS iLQR)

target_include_directories(iLQR INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
)
target_compile_features(iLQR INTERFACE cxx_std_17)
add_compile_options(-Wall -Wpedantic -Wextra -ffunction-sections -fdata-sections)
add_link_options(-Wl,--gc-sections -Wl,--start-group -Wl,--end-group)

# 使用 CMake 标准 BUILD_TESTING 选项控制是否构建测试。
include(CTest)
if(BUILD_TESTING)
  add_subdirectory(Tests)
endif()