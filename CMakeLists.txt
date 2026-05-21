cmake_minimum_required(VERSION 3.25)
project(nbody LANGUAGES CXX)

find_package(OpenMP REQUIRED)

set(SOURCES
    src/main.cpp
    src/Output/Output.cpp
    src/Force/Force.cpp
    src/Force/BarnesHut.cpp
    src/Integrator/Integrator.cpp
    src/Simulation/Simulation.cpp
)

set(CORE_SOURCES
    src/Force/Force.cpp
    src/Force/BarnesHut.cpp
    src/Integrator/Integrator.cpp
)

set(TEST_SOURCES
    tests/unit_tests/main.cpp
    tests/unit_tests/integrator.cpp
    tests/unit_tests/force.cpp
    tests/unit_tests/kepler.cpp
    tests/unit_tests/particle.cpp
    tests/unit_tests/barnes_hut.cpp
)

function(configure_target target_name)
    target_include_directories(${target_name} PRIVATE src)
    target_compile_features(${target_name} PRIVATE cxx_std_23)
    set_target_properties(${target_name} PROPERTIES CXX_EXTENSIONS OFF)
    if(MSVC)
        target_compile_options(${target_name} PRIVATE /O2 /arch:AVX2 /W4 /WX)
    else()
        target_compile_options(${target_name} PRIVATE -O3 -march=native -Wall -Wextra -Wpedantic -Werror)
    endif()
    target_link_libraries(${target_name} PRIVATE OpenMP::OpenMP_CXX)
endfunction()

add_executable(main ${SOURCES})
configure_target(main)

add_executable(tests ${TEST_SOURCES} ${CORE_SOURCES})
configure_target(tests)

add_executable(benchmark tests/benchmark.cpp ${CORE_SOURCES})
configure_target(benchmark)