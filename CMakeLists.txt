cmake_minimum_required(VERSION 3.20)
project(FXPricer VERSION 3.0.0 LANGUAGES CXX
    DESCRIPTION "Production C++20 FX Option Pricer")

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    add_compile_options(-Wall -Wextra -Wpedantic -fcoroutines -O2)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    add_compile_options(-Wall -Wextra -Wpedantic -fcoroutines-ts -O2)
endif()

add_library(pricer INTERFACE)
target_include_directories(pricer INTERFACE
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>)

# Link pthreads for std::jthread
find_package(Threads REQUIRED)

# ── Optional pybind11 Python bindings ────────────────────────────
option(BUILD_PYTHON "Build pybind11 Python bindings" OFF)

if(BUILD_PYTHON)
    find_package(pybind11 REQUIRED)
    pybind11_add_module(fx_pricer_py python/bindings.cpp)
    target_link_libraries(fx_pricer_py PRIVATE pricer Threads::Threads)
    set_target_properties(fx_pricer_py PROPERTIES OUTPUT_NAME fx_pricer)

    # Embed the compiler's libstdc++ directory into the module's RPATH
    # so it finds GLIBCXX_3.4.32 at import time (needed on Colab/Ubuntu 22.04)
    execute_process(
        COMMAND ${CMAKE_CXX_COMPILER} -print-file-name=libstdc++.so.6
        OUTPUT_VARIABLE _GCC_STDCPP_PATH
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(_GCC_STDCPP_PATH)
        cmake_path(GET _GCC_STDCPP_PATH PARENT_PATH _GCC_STDCPP_DIR)
        set_target_properties(fx_pricer_py PROPERTIES
            BUILD_RPATH "${_GCC_STDCPP_DIR}"
            INSTALL_RPATH "${_GCC_STDCPP_DIR}")
    endif()
else()
    # Build executable and tests only when not building Python module
    add_executable(fx_pricer src/main.cpp)
    target_link_libraries(fx_pricer PRIVATE pricer Threads::Threads)

    add_executable(tests tests/test_pricer.cpp)
    target_link_libraries(tests PRIVATE pricer Threads::Threads)
    enable_testing()
    add_test(NAME unit COMMAND tests)
endif()
