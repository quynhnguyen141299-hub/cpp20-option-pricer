cmake_minimum_required(VERSION 3.20)
project(FXPricer VERSION 3.0.0 LANGUAGES CXX
    DESCRIPTION "Production C++20 FX Option Pricer")

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    add_compile_options(-Wall -Wextra -Wpedantic -fcoroutines -O2 -march=native)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    add_compile_options(-Wall -Wextra -Wpedantic -fcoroutines-ts -O2 -march=native)
endif()

add_library(pricer INTERFACE)
target_include_directories(pricer INTERFACE
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>)

# Link pthreads for std::jthread
find_package(Threads REQUIRED)

add_executable(fx_pricer src/main.cpp)
target_link_libraries(fx_pricer PRIVATE pricer Threads::Threads)

add_executable(tests tests/test_pricer.cpp)
target_link_libraries(tests PRIVATE pricer Threads::Threads)
enable_testing()
add_test(NAME unit COMMAND tests)

# ── Optional pybind11 Python bindings ────────────────────────────
option(BUILD_PYTHON "Build pybind11 Python bindings" OFF)

if(BUILD_PYTHON)
    find_package(pybind11 REQUIRED)
    pybind11_add_module(fx_pricer_py python/bindings.cpp)
    target_link_libraries(fx_pricer_py PRIVATE pricer Threads::Threads)
    set_target_properties(fx_pricer_py PROPERTIES OUTPUT_NAME fx_pricer)
endif()
