include_directories("${CMAKE_SOURCE_DIR}/src")

find_package(Catch2 3 REQUIRED)

# Add -fpermissive for gcc as a workaround for: https://github.com/catchorg/Catch2/issues/2529
if(NOT MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive")
endif()

file(GLOB_RECURSE TEST_SOURCES CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/test/unit/*.cpp")
file(GLOB_RECURSE TEST_HEADERS CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/test/unit/*.hpp")

add_executable(unit_tests ${TEST_SOURCES} ${TEST_HEADERS})

target_link_libraries(unit_tests PRIVATE Catch2::Catch2WithMain)

add_test(NAME unit COMMAND unit_tests)