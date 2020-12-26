file(GLOB_RECURSE sources  *.f90)
add_executable(bem ${sources})
target_link_libraries(bem blas lapack)
