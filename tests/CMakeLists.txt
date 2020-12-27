file(GLOB_RECURSE sources  *.cpp)
add_executable(run_tests ${sources})
target_link_libraries(run_tests /home/jimmy/code/beginner-bem/Debug/src/libcelap.a gfortran /home/jimmy/opt/scifor/gnu/4.6.33-2-gcb22c67/lib/libscifor.a openblas gtest pthread)
