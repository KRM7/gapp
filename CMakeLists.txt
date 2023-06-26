cmake_minimum_required(VERSION 3.21)

project(gapp VERSION 0.1 LANGUAGES CXX)

include(CheckIPOSupported)
include(CTest)
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

option(GAPP_BUILD_TESTS "Build the tests for the library when ON." ON)
option(GAPP_BUILD_BENCHMARKS "Build the benchmarks for the library when ON." OFF)
option(GAPP_BUILD_EXAMPLES "Build the examples for the library when ON." OFF)
option(GAPP_USE_WERROR "Treat all warnings as errors during the build." ON)
option(GAPP_USE_LTO "Use link time optimizations for the library." ON)
option(GAPP_USE_MARCH_NATIVE "Optimize for the host architecture in release builds (only used for gcc and clang)." OFF)
option(GAPP_DISABLE_EXCEPTIONS "Disable exception support when building the library." OFF)
option(GAPP_DISABLE_RTTI "Disable run-time type information when building the library." OFF)
option(GAPP_LINK_TBB "Link Intel's TBB library (needed for the parallel stl algorithms when using gcc/clang)." OFF)


set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build")

# Set the default build type to release
if(NOT CMAKE_BUILD_TYPE MATCHES "^(Debug|Release|RelWithDebInfo)$")
    message(WARNING "The specified build type [${CMAKE_BUILD_TYPE}] is not recognized. Defaulting to RelWithDebInfo.")
    set(CMAKE_BUILD_TYPE "RelWithDebInfo")
endif()

# Shared library specific settings
if(BUILD_SHARED_LIBS)
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()

# Prefix static libraries with "lib" on all platforms
set(CMAKE_STATIC_LIBRARY_PREFIX "lib")

# Always generate compile_commands.json
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)


# Compiler / linker setup
if(MSVC) # MSVC style compiler interface
    set(GAPP_WARN_FLAGS "-permissive- -W4 -w44263 -w44265 -w44388 -w44471 -w44820 -w44986 -w45024 -w45025 -w45038 -w45204 -w45233 -w45258 -w45262 -w45263 -w45264 -external:W0")
    if(GAPP_USE_WERROR)
        set(GAPP_WARN_FLAGS "${GAPP_WARN_FLAGS} -WX")
    endif()

    set(GAPP_OPT_FLAGS "-O2 -Oi -Gy")
    set(GAPP_OPT_LINK_FLAGS "-OPT:REF -OPT:ICF -OPT:LBR")

    if(GAPP_DISABLE_EXCEPTIONS)
        set(GAPP_CXX_FLAGS "${GAPP_CXX_FLAGS} -D_HAS_EXCEPTIONS=0")
    else()
        set(GAPP_CXX_FLAGS "${GAPP_CXX_FLAGS} -EHsc")
    endif()

    if(GAPP_DISABLE_RTTI)
        set(GAPP_CXX_FLAGS "${GAPP_CXX_FLAGS} -GR-")
    endif()

    # msvc specific options
    if(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
        set(GAPP_CXX_FLAGS "${GAPP_CXX_FLAGS} -Zc:throwingNew")
        set(GAPP_WARN_FLAGS "${GAPP_WARN_FLAGS} -external:anglebrackets")
        set(GAPP_OPT_FLAGS "${GAPP_OPT_FLAGS} -fp:contract")
    endif()
    # clang-cl specific options
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set(GAPP_WARN_FLAGS "${GAPP_WARN_FLAGS}")
    endif()

    set(CMAKE_CXX_FLAGS "${GAPP_CXX_FLAGS} -Z7 -diagnostics:caret ${GAPP_WARN_FLAGS}")
    set(CMAKE_CXX_FLAGS_DEBUG "-D_DEBUG -RTC1 -Od")
    set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG ${GAPP_OPT_FLAGS}")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${GAPP_OPT_FLAGS}")

    link_libraries("-DEBUG:FULL" "-pdbaltpath:%_PDB%" "-NATVIS:${CMAKE_CURRENT_SOURCE_DIR}/gapp.natvis")
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} ${GAPP_OPT_LINK_FLAGS}")
    set(CMAKE_SHARED_LINKER_FLAGS_RELEASE "${CMAKE_SHARED_LINKER_FLAGS_RELEASE} ${GAPP_OPT_LINK_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO "${CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO} ${GAPP_OPT_LINK_FLAGS}")
    set(CMAKE_SHARED_LINKER_FLAGS_RELWITHDEBINFO "${CMAKE_SHARED_LINKER_FLAGS_RELWITHDEBINFO} ${GAPP_OPT_LINK_FLAGS}")

else() # GNU style compiler interface
    set(GAPP_WARN_FLAGS "-Wall -Wextra -Wpedantic -Wshadow")
    if(GAPP_USE_WERROR)
        set(GAPP_WARN_FLAGS "${GAPP_WARN_FLAGS} -Werror -pedantic-errors")
    endif()

    set(GAPP_OPT_FLAGS "-O3 -fno-math-errno -fno-signed-zeros -fno-trapping-math -freciprocal-math -fno-rounding-math")
    if(GAPP_USE_MARCH_NATIVE)
        set(GAPP_OPT_FLAGS "${GAPP_OPT_FLAGS} -march=native")
    endif()

    if(GAPP_DISABLE_EXCEPTIONS)
        set(GAPP_CXX_FLAGS "${GAPP_CXX_FLAGS} -fno-exceptions")
    endif()

    if(GAPP_DISABLE_RTTI)
        set(GAPP_CXX_FLAGS "${GAPP_CXX_FLAGS} -fno-rtti")
    endif()

    # gcc specific options
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
	    set(GAPP_WARN_FLAGS "${GAPP_WARN_FLAGS} -Wlogical-op")
    endif()
    # clang specific options
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set(GAPP_OPT_FLAGS "${GAPP_OPT_FLAGS} -fapprox-func")
    endif()

    set(CMAKE_CXX_FLAGS "${GAPP_CXX_FLAGS} ${GAPP_WARN_FLAGS}")
    set(CMAKE_CXX_FLAGS_DEBUG "-Og")
    set(CMAKE_CXX_FLAGS_RELEASE "-g -DNDEBUG ${GAPP_OPT_FLAGS}")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g ${GAPP_OPT_FLAGS}")
endif()

# LTO setup
check_ipo_supported(RESULT LTO_SUPPORTED)
if(GAPP_USE_LTO AND LTO_SUPPORTED AND CMAKE_BUILD_TYPE MATCHES "(Release|RelWithDebInfo)")
    if (BUILD_SHARED_LIBS AND MSVC)
        message(WARNING "LTO is not available for shared libs using MSVC because of CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS.")
    else()
        message(STATUS "Using link time optimizations.")
        set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
    endif()
endif()


# Library
file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")
file(GLOB_RECURSE HEADERS CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/*.hpp")

add_library(gapp ${SOURCES} ${HEADERS})
add_library(gapp::gapp ALIAS gapp)

target_include_directories(gapp SYSTEM INTERFACE "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>"
                                                 "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
                                PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")

target_compile_features(gapp PUBLIC "cxx_std_20")
target_compile_definitions(gapp PUBLIC "$<$<NOT:$<BOOL:${BUILD_SHARED_LIBS}>>:GAPP_BUILD_STATIC>")
target_compile_options(gapp PUBLIC "$<$<CXX_COMPILER_ID:MSVC>:-Zc:preprocessor>")

find_package(TBB)
if(TBB_FOUND AND GAPP_LINK_TBB AND NOT MSVC)
    target_link_libraries(gapp PUBLIC "TBB::tbb")
endif()


# Install
install(TARGETS gapp EXPORT gapp-config
        ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}/$<CONFIG>"
        LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}/$<CONFIG>"
        RUNTIME DESTINATION "${CMAKE_INSTALL_LIBDIR}/$<CONFIG>")

install(DIRECTORY "src/" DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/gapp" FILES_MATCHING PATTERN "*.hpp" PATTERN "*.txt")
install(DIRECTORY "docs/" DESTINATION "${CMAKE_INSTALL_DOCDIR}" FILES_MATCHING PATTERN "*.md")

if(MSVC AND BUILD_SHARED_LIBS)
    install(FILES "$<TARGET_PDB_FILE:gapp>" DESTINATION "${CMAKE_INSTALL_LIBDIR}/$<CONFIG>")
endif()

export(TARGETS gapp NAMESPACE gapp::
       FILE "${CMAKE_CURRENT_BINARY_DIR}/gapp-config.cmake")

install(EXPORT gapp-config
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/gapp"
        NAMESPACE gapp::)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/gapp-config-version.cmake"
    VERSION "${PROJECT_VERSION}"
    COMPATIBILITY "SameMajorVersion"
)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/gapp-config-version.cmake"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/gapp")


# Tests / examples
if(GAPP_BUILD_TESTS AND BUILD_TESTING AND PROJECT_IS_TOP_LEVEL)
    enable_testing()
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/test/unit")
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/test/integration")
endif()

if(GAPP_BUILD_BENCHMARKS AND BUILD_TESTING AND PROJECT_IS_TOP_LEVEL)
    enable_testing()
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/test/benchmark")
endif()

if(GAPP_BUILD_EXAMPLES AND PROJECT_IS_TOP_LEVEL)
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/examples")
endif()
