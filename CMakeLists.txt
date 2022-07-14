cmake_minimum_required(VERSION 3.20)

project(GeneticAlgorithm VERSION 0.1 LANGUAGES CXX)

include(GNUInstallDirs)

option(GA_ENABLE_TESTS "Enable the tests for the library." OFF) # TODO
option(GA_ENABLE_EXAMPLES "Enable building the examples for the library." OFF)
option(GA_ENABLE_LTO "Enable link time optimizations for the library." ON)
option(GA_ENABLE_MARCH_NATIVE "Add -march=native to release builds when using gcc or clang." ON)


# Set the output directory for the build files
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/build")

# Set the default build type to release
if(NOT CMAKE_BUILD_TYPE MATCHES "^(Debug|Release)$")
    message("Unrecognized build type: ${CMAKE_BUILD_TYPE}. Defaulting to Release.")
    set(CMAKE_BUILD_TYPE "Release")
endif()


# Compiler/linker flags
if(MSVC)
    set(CMAKE_CXX_FLAGS "-std:c++latest -permissive- -Zc:preprocessor -Zc:rvalueCast -W4 -WX -TP -MP -MD -fp:precise -EHsc")
    set(CMAKE_CXX_FLAGS_DEBUG "-Od -ZI -RTC1")
    set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -O2 -Oi -Gy -Qpar")
    if(GA_ENABLE_LTO)
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -GL -LTCG")
    endif()
elseif()
    set(CMAKE_CXX_FLAGS "-std=c++20 -Wall -Wextra -Wshadow -Wpedantic -pedantic-errors -Werror")
    set(CMAKE_CXX_FLAGS_DEBUG "-Og") # -ltbb_debug
    set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -O3") # link intel tbb for par stl algos -ltbb
    if (GA_ENABLE_MARCH_NATIVE)
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native")
    endif()
    if (GA_ENABLE_LTO)
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -flto -ffat-lto-objects -Wno-lto-type-mismatch")
    endif()
endif()


# Library
file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/src/*.cpp")
file(GLOB_RECURSE HEADERS CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/src/*.hpp")

add_library(GeneticAlgorithm STATIC ${SOURCES} ${HEADERS})

target_include_directories(GeneticAlgorithm SYSTEM INTERFACE "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>"
                                                             "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/GeneticAlgorithm>"
                                            PRIVATE "${CMAKE_SOURCE_DIR}/src")

# Install library
install(TARGETS GeneticAlgorithm EXPORT GeneticAlgorithmConfig
        ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}")

install(DIRECTORY "src/" DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/GeneticAlgorithm")
install(DIRECTORY "docs/" DESTINATION "${CMAKE_INSTALL_DOCDIR}/GeneticAlgorithm")

export(TARGETS GeneticAlgorithm NAMESPACE GeneticAlgorithm:: FILE "${CMAKE_BINARY_DIR}/GeneticAlgorithmConfig.cmake")

install(EXPORT GeneticAlgorithmConfig
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/GeneticAlgorithm"
        NAMESPACE GeneticAlgorithm::)


# Benchmarks / tests
if(GA_ENABLE_TESTS)
    add_subdirectory("test")
endif()

# Examples
if(GA_ENABLE_EXAMPLES)
    add_subdirectory("examples")
endif()