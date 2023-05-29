add_executable(example1 "${CMAKE_SOURCE_DIR}/examples/1_intro_minimal.cpp")
add_executable(example2 "${CMAKE_SOURCE_DIR}/examples/2_intro_detailed.cpp")

target_link_libraries(example1 PRIVATE gapp)
target_link_libraries(example2 PRIVATE gapp)
