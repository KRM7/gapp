find_package(Catch2 3 REQUIRED)
find_package(TBB)

if(GA_LINK_TBB)
    message("Using Intel TBB library.")
    set(TBB_LIB "TBB::tbb")
endif()

add_executable(example1 "${CMAKE_SOURCE_DIR}/examples/1_intro_minimal.cpp")
add_executable(example2 "${CMAKE_SOURCE_DIR}/examples/2_intro_detailed.cpp")

target_link_libraries(example1 PRIVATE GeneticAlgorithm ${TBB_LIB})
target_link_libraries(example2 PRIVATE GeneticAlgorithm ${TBB_LIB})