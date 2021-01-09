cmake_minimum_required(VERSION 3.15)

project(beginnerbem)

enable_language(Fortran)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(dialect "-ffree-form -std=f2008 -fimplicit-none")
  set(bounds "-fbounds-check")
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(dialect "-stand f08 -free -implicitnone")
  set(bounds "-check bounds")
endif()

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")
set(CMAKE_Fortran_MODULE_DIRECTORY "/home/jimmy/opt/scifor/gnu/4.6.33-2-gcb22c67/include/")

set(LIBNAME celap)

add_subdirectory(external/googletest)
enable_testing()
include_directories(${CMAKE_SOURCE_DIR}/src ${CMAKE_SOURCE_DIR}/external/SciFortran ${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
add_subdirectory(external/SciFortran)

add_subdirectory(src)
add_subdirectory(tests)
