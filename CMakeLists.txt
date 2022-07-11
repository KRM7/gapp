cmake_minimum_required(VERSION 3.20)

project(GeneticAlgorithm VERSION 0.1 LANGUAGES CXX)

# Set the output directory for the build files
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/build)

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
    message(FATAL_ERROR "Unsupported compiler: ${CMAKE_CXX_COMPILER_ID}")
endif()


# Library
file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/src/*.cpp")
file(GLOB_RECURSE HEADERS CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/src/*.hpp")

add_library(GeneticAlgorithm STATIC ${SOURCES} ${HEADERS})
target_include_directories(GeneticAlgorithm PUBLIC ${CMAKE_SOURCE_DIR}/src)
# need to check if this works with cygwin, try to use it in another test project
install(TARGETS GeneticAlgorithm ARCHIVE DESTINATION GeneticAlgorithm/lib)
install(DIRECTORY src/ DESTINATION GeneticAlgorithm/src)


# Benchmarks / tests
add_executable(benchmark test/benchmark.cpp)
target_link_libraries(benchmark PRIVATE GeneticAlgorithm)