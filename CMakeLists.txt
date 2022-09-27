cmake_minimum_required(VERSION 3.16)

project(GeneticAlgorithm VERSION 0.1 LANGUAGES CXX)

include(GNUInstallDirs)

option(GA_ENABLE_TESTING "Enable the tests for the library." OFF) # TODO
option(GA_ENABLE_EXAMPLES "Enable building the examples for the library." OFF) # TODO
option(GA_ENABLE_WERROR "Treat warnings as errors during the build." ON)
option(GA_ENABLE_LTO "Enable link time optimizations for the library." ON)
option(GA_ENABLE_MARCH_NATIVE "Add -march=native option to release builds when using gcc or clang." ON)
option(GA_ENABLE_LINK_TBB "Add -ltbb option when using gcc (needed for parallel stl algorithms when using gcc)." ON)
#option(GA_ENABLE_EXCEPTIONS "Enable exceptions in the library." ON)


# Set the output directory for the build files
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/build")

# Set the default build type to release
if(NOT CMAKE_BUILD_TYPE MATCHES "^(Debug|Release)$")
    message("Unrecognized build type: ${CMAKE_BUILD_TYPE}. Defaulting to Release.")
    set(CMAKE_BUILD_TYPE "Release")
endif()


# Compiler/linker setup
if(MSVC) # MSVC style compiler interface
    set(GA_WARN_FLAGS "-permissive- -W4 -w44263 -w44265 -w44388 -w44471 -wd4505")
    if(GA_ENABLE_WERROR)
        set(GA_WARN_FLAGS "${GA_WARN_FLAGS} -WX")
    endif()
    set(CMAKE_CXX_FLAGS "${GA_DEFINES} -std:c++latest -Zc:preprocessor -Zc:rvalueCast ${GA_WARN_FLAGS} -TP -MP -MD -fp:fast -EHsc")
    set(CMAKE_CXX_FLAGS_DEBUG "-Od -ZI -RTC1")
    set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -O2 -Oi -Gy -Qpar")
    if(GA_ENABLE_LTO)
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -GL -LTCG")
    endif()
else() # GNU style compiler interface
    set(GA_WARN_FLAGS "-Wall -Wextra -Wpedantic -Wshadow -Wlogical-op -Wnull-dereference -Wno-unused-function")
    if(GA_ENABLE_WERROR)
        set(GA_WARN_FLAGS "${GA_WARN_FLAGS} -Werror -pedantic-errors")
    endif()
    set(CMAKE_CXX_FLAGS "${GA_DEFINES} -std=c++20 ${GA_WARN_FLAGS}")
    set(CMAKE_CXX_FLAGS_DEBUG "-Og")
    set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -O3 -ffast-math -fno-finite-math-only")
    if (GA_ENABLE_LTO)
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -flto -ffat-lto-objects -Wno-lto-type-mismatch")
    endif()
    if (GA_ENABLE_MARCH_NATIVE)
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native")
    endif()
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU") # gcc specific options
        if (GA_ENABLE_LINK_TBB)
            set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -ltbb")
        endif()
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
if(GA_ENABLE_TESTING)
    add_subdirectory("test")
endif()

# Examples
if(GA_ENABLE_EXAMPLES)
    add_subdirectory("examples")
endif()