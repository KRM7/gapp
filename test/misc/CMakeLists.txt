
add_executable(rcu_test "${CMAKE_CURRENT_SOURCE_DIR}/rcu.cpp")

target_link_libraries(rcu_test PRIVATE gapp)
