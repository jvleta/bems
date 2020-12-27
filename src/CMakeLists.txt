file(GLOB_RECURSE sources  *.f90)
add_library(${LIBNAME} ${sources})
target_link_libraries(${LIBNAME} openblas lapack /home/jimmy/opt/scifor/gnu/4.6.33-2-gcb22c67/lib/libscifor.a)
