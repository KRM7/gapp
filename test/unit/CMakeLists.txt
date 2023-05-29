find_package(Catch2 3 REQUIRED)

if(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -wd4388")
endif()

file(GLOB_RECURSE TEST_SOURCES CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/test/unit/*.cpp")
file(GLOB_RECURSE TEST_HEADERS CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/test/unit/*.hpp")

add_executable(unit_tests ${TEST_SOURCES} ${TEST_HEADERS})

target_link_libraries(unit_tests PRIVATE Catch2::Catch2WithMain gapp)

include(Catch)
catch_discover_tests(unit_tests)
