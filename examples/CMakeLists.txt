find_package(Catch2 3 REQUIRED)
find_package(TBB)

if(GA_LINK_TBB)
    message("Using Intel TBB library.")
    set(TBB_LIB "TBB::tbb")
endif()

add_executable(example1 "${CMAKE_SOURCE_DIR}/examples/1_rcga_basic.cpp")

target_link_libraries(example1 PRIVATE GeneticAlgorithm ${TBB_LIB})