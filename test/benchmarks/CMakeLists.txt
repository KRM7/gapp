#find_package(Catch2 3 REQUIRED)

if(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -wd4388")
endif()

file(GLOB_RECURSE BENCHMARK_SOURCES CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/test/benchmarks/*.cpp")
file(GLOB_RECURSE BENCHMARK_HEADERS CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/test/benchmarks/*.hpp")

add_executable(benchmarks ${BENCHMARK_SOURCES} ${BENCHMARK_HEADERS})

target_link_libraries(benchmarks PRIVATE GeneticAlgorithm)

#include(Catch)
#catch_discover_tests(benchmarks)