
file(GLOB SOURCES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

foreach(SOURCE ${SOURCES})
    string(REPLACE ".cpp" "" EXAMPLE ${SOURCE})
    add_executable(${EXAMPLE} ${SOURCE})
    target_link_libraries(${EXAMPLE} PRIVATE gapp::gapp)
endforeach()
