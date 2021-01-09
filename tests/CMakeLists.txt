file(GLOB_RECURSE sources  *.cpp)
add_executable(run_tests ${sources})
target_link_libraries(run_tests celap gfortran scifor openblas gtest pthread)
