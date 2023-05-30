if(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -wd4388")
endif()

file(GLOB_RECURSE BENCHMARK_SOURCES CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/test/integration/*.cpp")
file(GLOB_RECURSE BENCHMARK_HEADERS CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/test/integration/*.hpp")

add_executable(integration_tests ${BENCHMARK_SOURCES} ${BENCHMARK_HEADERS})

target_link_libraries(integration_tests PRIVATE gapp)
