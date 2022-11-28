find_package(Catch2 3 REQUIRED)
find_package(TBB)

if(GA_LINK_TBB)
    message("Using Intel TBB library.")
    set(TBB_LIB "TBB::tbb")
endif()

# Add -fpermissive for gcc as a workaround for: https://github.com/catchorg/Catch2/issues/2529
if(NOT MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -wd4388")
endif()

file(GLOB_RECURSE BENCHMARK_SOURCES CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/test/benchmark/*.cpp")
file(GLOB_RECURSE BENCHMARK_HEADERS CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/test/benchmark/*.hpp")

add_executable(benchmarks ${BENCHMARK_SOURCES} ${BENCHMARK_HEADERS})

target_link_libraries(benchmarks PRIVATE GeneticAlgorithm ${TBB_LIB})

#include(Catch)
#catch_discover_tests(benchmarks)