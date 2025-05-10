add_executable(reproducibility_test "${CMAKE_CURRENT_SOURCE_DIR}/reproducibility.cpp")

target_link_libraries(reproducibility_test PRIVATE gapp)
