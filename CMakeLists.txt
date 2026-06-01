
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
option(ILQR_ENABLE_ASAN "Enable AddressSanitizer for targets using iLQR::Options" OFF)
set(ILQR_MATRIX_TARGET "" CACHE STRING "External Matrix target used by iLQR")
if(ILQR_IS_TOP_LEVEL)
  set(ILQR_FORCE_OPTIMIZATION_DEFAULT OFF)
else()
  set(ILQR_FORCE_OPTIMIZATION_DEFAULT ON)
endif()
option(ILQR_FORCE_OPTIMIZATION "Force optimized code generation for targets consuming iLQR" ${ILQR_FORCE_OPTIMIZATION_DEFAULT})

add_library(iLQR INTERFACE)
add_library(iLQR::iLQR ALIAS iLQR)

if(ILQR_MATRIX_TARGET)
  if(NOT TARGET "${ILQR_MATRIX_TARGET}")
    message(FATAL_ERROR "ILQR_MATRIX_TARGET='${ILQR_MATRIX_TARGET}' is not a CMake target.")
  endif()

  target_include_directories(iLQR BEFORE INTERFACE
    $<TARGET_PROPERTY:${ILQR_MATRIX_TARGET},INTERFACE_INCLUDE_DIRECTORIES>
  )
  target_link_libraries(iLQR INTERFACE "${ILQR_MATRIX_TARGET}")
else()
  add_library(iLQRMatrix INTERFACE)
  add_library(iLQR::Matrix ALIAS iLQRMatrix)

  target_include_directories(iLQRMatrix INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
  )
  target_link_libraries(iLQR INTERFACE iLQR::Matrix)
endif()

target_include_directories(iLQR INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
)
target_compile_features(iLQR INTERFACE cxx_std_17)

if(ILQR_FORCE_OPTIMIZATION)
  target_compile_options(iLQR INTERFACE
    -O3
  )
  target_compile_definitions(iLQR INTERFACE NDEBUG)
endif()

add_library(iLQROptions INTERFACE)
add_library(iLQR::Options ALIAS iLQROptions)

if(ILQR_ENABLE_DEV_OPTIONS)
  target_compile_options(iLQROptions INTERFACE
    -Wall
    -Wpedantic
    -Wextra
    -ffunction-sections
    -fdata-sections
  )
  target_link_options(iLQROptions INTERFACE
    -Wl,--gc-sections
    -Wl,--start-group
    -Wl,--end-group
  )
endif()

if(ILQR_ENABLE_ASAN)
  target_compile_options(iLQROptions INTERFACE
    -fsanitize=address
    -fno-omit-frame-pointer
  )
  target_link_options(iLQROptions INTERFACE
    -fsanitize=address
  )
endif()

if(ILQR_BUILD_TESTING)
  include(CTest)
  if(BUILD_TESTING)
    add_subdirectory(Tests)
  endif()
endif()