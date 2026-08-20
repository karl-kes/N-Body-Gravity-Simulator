cmake_minimum_required(VERSION 3.25)

project(nbody LANGUAGES CXX CUDA)

option(FP32 "Use single-precision floating point storage and kernels" OFF)

if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    message(FATAL_ERROR
        "nbody requires GCC 15.x (found ${CMAKE_CXX_COMPILER_ID} "
        "${CMAKE_CXX_COMPILER_VERSION}).")
endif()

if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 15
   OR CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 16)
    message(FATAL_ERROR
        "nbody requires GCC 15.x (found ${CMAKE_CXX_COMPILER_VERSION}).")
endif()

if(NOT CMAKE_CUDA_COMPILER_ID STREQUAL "NVIDIA")
    message(FATAL_ERROR
        "nbody requires NVIDIA CUDA 13.3.x (found ${CMAKE_CUDA_COMPILER_ID} "
        "${CMAKE_CUDA_COMPILER_VERSION}).")
endif()

if(CMAKE_CUDA_COMPILER_VERSION VERSION_LESS 13.3
   OR CMAKE_CUDA_COMPILER_VERSION VERSION_GREATER_EQUAL 13.4)
    message(FATAL_ERROR
        "nbody requires NVIDIA CUDA 13.3.x "
        "(found ${CMAKE_CUDA_COMPILER_VERSION}).")
endif()

if(NOT DEFINED CMAKE_CUDA_HOST_COMPILER)
    message(FATAL_ERROR
        "Set CMAKE_CUDA_HOST_COMPILER to the same GCC 15 compiler used for C++. "
        "For example: -DCMAKE_CXX_COMPILER=g++-15 "
        "-DCMAKE_CUDA_HOST_COMPILER=g++-15")
endif()

file(REAL_PATH "${CMAKE_CXX_COMPILER}" NBODY_CXX_COMPILER_REAL)
file(REAL_PATH "${CMAKE_CUDA_HOST_COMPILER}" NBODY_CUDA_HOST_COMPILER_REAL)

if(NOT NBODY_CXX_COMPILER_REAL STREQUAL NBODY_CUDA_HOST_COMPILER_REAL)
    message(FATAL_ERROR
        "C++ and CUDA must use the same host compiler. C++ uses "
        "${NBODY_CXX_COMPILER_REAL}, while CUDA uses "
        "${NBODY_CUDA_HOST_COMPILER_REAL}.")
endif()

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CUDA_STANDARD 23)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_EXTENSIONS OFF)

include(FetchContent)

set(XPU_ENABLE_CUDA ON CACHE BOOL "Build xpu with CUDA support" FORCE)
set(XPU_BUILD_TESTS OFF CACHE BOOL "Build xpu's tests" FORCE)

FetchContent_Declare(
    xpu
    GIT_REPOSITORY https://github.com/UWHPC/xpu.git
    GIT_TAG 598efa0403f33566c30be51bd2e0c0595c48b036
    GIT_SHALLOW FALSE
)

FetchContent_MakeAvailable(xpu)

find_package(OpenMP REQUIRED)

set(SOURCES
    src/main.cu
    src/Output/Output.cu
    src/Force/Force.cu
    src/Force/BarnesHut.cu
    src/Integrator/Integrator.cu
    src/Simulation/Simulation.cu
)

set(CORE_SOURCES
    src/Force/Force.cu
    src/Force/BarnesHut.cu
    src/Integrator/Integrator.cu
)

set(TEST_SOURCES
    tests/unit_tests/main.cu
    tests/unit_tests/integrator.cu
    tests/unit_tests/force.cu
    tests/unit_tests/kepler.cu
    tests/unit_tests/particle.cu
    tests/unit_tests/barnes_hut.cu
)

function(configure_target target_name)
    target_include_directories(${target_name} PRIVATE src)
    if(FP32)
        target_compile_definitions(${target_name} PRIVATE FP32)
    endif()
    target_compile_options(${target_name} PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-O3;-march=native;-Wall;-Wextra;-Wpedantic;-Werror>
        $<$<COMPILE_LANGUAGE:CUDA>:-O3;-std=c++23;-Xcompiler=-Wall,-Wextra,-Wpedantic,-Werror,-fopenmp>
    )
    target_link_libraries(${target_name} PRIVATE OpenMP::OpenMP_CXX xpu::xpu)
endfunction()

add_executable(main ${SOURCES})
configure_target(main)

add_executable(tests ${TEST_SOURCES} ${CORE_SOURCES})
configure_target(tests)

add_executable(benchmark tests/benchmark.cu ${CORE_SOURCES})
configure_target(benchmark)
