
cmake_minimum_required(VERSION 3.20)

set(PROJECT_NAME "iLQR")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)
project(${PROJECT_NAME} LANGUAGES CXX)

if(DEFINED PROJECT_IS_TOP_LEVEL)
  set(ILQR_IS_TOP_LEVEL ${PROJECT_IS_TOP_LEVEL})
else()
  if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
    set(ILQR_IS_TOP_LEVEL ON)
  else()
    set(ILQR_IS_TOP_LEVEL OFF)
  endif()
endif()

option(ILQR_BUILD_TESTING "Build iLQR tests" ${ILQR_IS_TOP_LEVEL})
option(ILQR_ENABLE_DEV_OPTIONS "Enable iLQR development compile/link options" ${ILQR_IS_TOP_LEVEL})
if(ILQR_IS_TOP_LEVEL)
  set(ILQR_FORCE_OPTIMIZATION_DEFAULT OFF)
else()
  set(ILQR_FORCE_OPTIMIZATION_DEFAULT ON)
endif()
option(ILQR_FORCE_OPTIMIZATION "Force optimized code generation for targets consuming iLQR" ${ILQR_FORCE_OPTIMIZATION_DEFAULT})

add_library(iLQR INTERFACE)
add_library(iLQR::iLQR ALIAS iLQR)

target_include_directories(iLQR INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
)
target_compile_features(iLQR INTERFACE cxx_std_17)

if(ILQR_FORCE_OPTIMIZATION)
  target_compile_options(iLQR INTERFACE
    $<$<CXX_COMPILER_ID:GNU,Clang>:-O3>
  )
  target_compile_definitions(iLQR INTERFACE NDEBUG)
endif()

add_library(iLQROptions INTERFACE)
add_library(iLQR::Options ALIAS iLQROptions)

if(ILQR_ENABLE_DEV_OPTIONS)
  target_compile_options(iLQROptions INTERFACE
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall>
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Wpedantic>
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Wextra>
    $<$<CXX_COMPILER_ID:GNU,Clang>:-ffunction-sections>
    $<$<CXX_COMPILER_ID:GNU,Clang>:-fdata-sections>
  )
  target_link_options(iLQROptions INTERFACE
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Wl,--gc-sections>
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Wl,--start-group>
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Wl,--end-group>
  )
endif()

if(ILQR_BUILD_TESTING)
  include(CTest)
  if(BUILD_TESTING)
    add_subdirectory(Tests)
  endif()
endif()