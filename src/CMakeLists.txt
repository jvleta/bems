file(GLOB_RECURSE sources  *.f90)
add_executable(bem ${sources})
target_link_libraries(bem openblas lapack /home/jimmy/opt/scifor/gnu/4.6.33-2-gcb22c67/lib/libscifor.a)
