find_package(Catch2 3 REQUIRED)
include(Catch)

if(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -wd4388")
endif()

file(GLOB SOURCES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

foreach(SOURCE ${SOURCES})
    string(REPLACE ".cpp" "_benchmark" BENCHMARK ${SOURCE})
    add_executable(${BENCHMARK} ${SOURCE})
    target_link_libraries(${BENCHMARK} PRIVATE Catch2::Catch2WithMain gapp::gapp)
    catch_discover_tests(${BENCHMARK})
endforeach()
