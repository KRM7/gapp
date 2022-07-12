cmake_minimum_required(VERSION 3.20)

include(GNUInstallDirs)

project(GeneticAlgorithm VERSION 0.1 LANGUAGES CXX)

# Set the output directory for the build files
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/build")

# Set the default build type to release
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE "Release")
endif()

# Compilers
if(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    # add option for lto (/GL compiler and /LTCG linker flags)
    set(CMAKE_CXX_FLAGS "-std:c++latest -permissive- -Zc:preprocessor -Zc:rvalueCast -W4 -WX -TP -MP -MD -fp:precise -EHsc")
    set(CMAKE_CXX_FLAGS_DEBUG "-Od -ZI -RTC1")
    set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -O2 -Oi -Gy -Qpar")
else()
    # clang and gcc, maybe intel
    # TODO (intel TBB for par stl, flags)
    # message(FATAL_ERROR "Unsupported compiler: ${CMAKE_CXX_COMPILER_ID}")
    set(CMAKE_CXX_FLAGS "-std=c++20 -Wall -Wextra")
endif()


# Library
file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/src/*.cpp")
file(GLOB_RECURSE HEADERS CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/src/*.hpp")

add_library(GeneticAlgorithm STATIC ${SOURCES} ${HEADERS})

target_include_directories(GeneticAlgorithm PUBLIC 
                            "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>"
                            "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/GeneticAlgorithm>")

# Install library
install(TARGETS GeneticAlgorithm EXPORT GeneticAlgorithmConfig
        ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}")

install(DIRECTORY "src/" DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/GeneticAlgorithm")
install(DIRECTORY "docs/" DESTINATION "${CMAKE_INSTALL_DOCDIR}/GeneticAlgorithm")

export(TARGETS GeneticAlgorithm
       NAMESPACE GeneticAlgorithm:: FILE "build/GeneticAlgorithmConfig.cmake")

install(EXPORT GeneticAlgorithmConfig
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/GeneticAlgorithm"
        NAMESPACE GeneticAlgorithm::)


# Benchmarks / tests
#add_executable(benchmark test/benchmark.cpp)
#target_link_libraries(benchmark PRIVATE GeneticAlgorithm)