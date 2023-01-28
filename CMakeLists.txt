cmake_minimum_required(VERSION 3.16)

project(GeneticAlgorithm VERSION 0.1 LANGUAGES CXX)

include(GNUInstallDirs)
include(CTest)
include(CheckIPOSupported)

option(GA_BUILD_TESTS "Enable the tests for the library." ON)
option(GA_BUILD_EXAMPLES "Enable building the examples for the library." OFF) # TODO
option(GA_USE_WERROR "Treat warnings as errors during the build." ON)
option(GA_USE_LTO "Enable link time optimizations for the library." ON)
option(GA_USE_MARCH_NATIVE "Add -march=native option to release builds when using gcc or clang." ON)
option(GA_LINK_TBB "Add -ltbb to the compiler flags (needed for parallel stl algorithms when using gcc/clang)." OFF)
#option(GA_ENABLE_EXCEPTIONS "Enable exceptions in the library." ON) # TODO


# Set the output directory for the build files
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/build")

# Set the default build type to release
if(NOT CMAKE_BUILD_TYPE MATCHES "^(Debug|Release|RelWithDebInfo)$")
    message("Unrecognized build type: ${CMAKE_BUILD_TYPE}. Defaulting to Release.")
    set(CMAKE_BUILD_TYPE "Release")
endif()

# Shared library specific settings
if(BUILD_SHARED_LIBS)
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()

# Always generate compile_commands.json
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)


# Compiler/linker setup
if(MSVC) # MSVC style compiler interface
    set(GA_WARN_FLAGS "-permissive- -W4 -w44263 -w44265 -w44388 -w44471")
    if(GA_USE_WERROR)
        set(GA_WARN_FLAGS "${GA_WARN_FLAGS} -WX")
    endif()
    set(CMAKE_CXX_FLAGS "${GA_CXX_FLAGS} -std:c++20 -Zc:preprocessor -Zc:rvalueCast ${GA_WARN_FLAGS} -TP -MP -EHsc")
    set(CMAKE_CXX_FLAGS_DEBUG "-D_DEBUG -Z7 -RTC1 -Od")
    set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -Zi -O2 -Oi -Gy -Qpar")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-Zi -O2 -Oi -Gy -Qpar")
else() # GNU style compiler interface
    set(GA_WARN_FLAGS "-Wall -Wextra -Wpedantic -Wshadow")
    if(GA_USE_WERROR)
        set(GA_WARN_FLAGS "${GA_WARN_FLAGS} -Werror -pedantic-errors")
    endif()
    set(CMAKE_CXX_FLAGS "${GA_CXX_FLAGS} -std=c++20 ${GA_WARN_FLAGS}")
    set(CMAKE_CXX_FLAGS_DEBUG "-Og")
    set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -O3 -ffast-math -fno-finite-math-only")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -O3 -ffast-math -fno-finite-math-only")
    if (GA_USE_MARCH_NATIVE)
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native")
    endif()
    # gcc specific options
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
	    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wlogical-op")
    endif()
    # clang specific options
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libstdc++")
    endif()
endif()

# LTO setup
check_ipo_supported(RESULT LTO_SUPPORTED)
if(GA_USE_LTO AND LTO_SUPPORTED AND CMAKE_BUILD_TYPE MATCHES "(Release|RelWithDebInfo)")
    if (BUILD_SHARED_LIBS AND MSVC)
        message("LTO is not available for shared libs using MSVC.")
    else()
        message("Using link time optimizations.")
        set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
    endif()
endif()


# Library
file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/src/*.cpp")
file(GLOB_RECURSE HEADERS CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/src/*.hpp")

add_library(GeneticAlgorithm ${SOURCES} ${HEADERS})

target_include_directories(GeneticAlgorithm SYSTEM INTERFACE "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>"
                                                             "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/GeneticAlgorithm>"
                                            PRIVATE "${CMAKE_SOURCE_DIR}/src")

target_compile_definitions(GeneticAlgorithm PUBLIC "$<$<NOT:$<BOOL:${BUILD_SHARED_LIBS}>>:GA_BUILD_STATIC>")

# Install library
install(TARGETS GeneticAlgorithm EXPORT GeneticAlgorithmConfig
        ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}")

install(DIRECTORY "src/" DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/GeneticAlgorithm")
#install(DIRECTORY "docs/" DESTINATION "${CMAKE_INSTALL_DOCDIR}/GeneticAlgorithm")

export(TARGETS GeneticAlgorithm NAMESPACE GeneticAlgorithm:: FILE "${CMAKE_BINARY_DIR}/GeneticAlgorithmConfig.cmake")

install(EXPORT GeneticAlgorithmConfig
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/GeneticAlgorithm"
        NAMESPACE GeneticAlgorithm::)


# Benchmarks / tests
if(GA_BUILD_TESTS AND BUILD_TESTING)
    enable_testing()
    add_subdirectory("${CMAKE_SOURCE_DIR}/test/unit")
    add_subdirectory("${CMAKE_SOURCE_DIR}/test/benchmarks")
endif()


# Examples
if(GA_BUILD_EXAMPLES)
    add_subdirectory("${CMAKE_SOURCE_DIR}/examples")
endif()