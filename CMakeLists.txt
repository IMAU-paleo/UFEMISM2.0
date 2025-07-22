cmake_minimum_required(VERSION 3.10)
project(UFEMISM LANGUAGES Fortran)

set(CMAKE_Fortran_STANDARD 2018)

file(GLOB_RECURSE UFEMISM_SOURCES "src/*.f90")

add_executable(UFEMISM_program ${UFEMISM_SOURCES})

# Enable preprocessing for all Fortran files
target_compile_options(UFEMISM_program PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-cpp>)

# Optionally enable assertions
option(DO_ASSERTIONS "Compile UFEMISM with assertions" ON)
if(DO_ASSERTIONS)
    target_compile_definitions(UFEMISM_program PRIVATE DO_ASSERTIONS)
endif()

# Optionally enable resource tracking
option(DO_RESOURCE_TRACKING "Compile UFEMISM with resource tracking" ON)
if(DO_RESOURCE_TRACKING)
    target_compile_definitions(UFEMISM_program PRIVATE DO_RESOURCE_TRACKING)
endif()

# Optionally provide extra compiler flags
set(EXTRA_Fortran_FLAGS "" CACHE STRING "Extra gfortran compiler flags")
if(EXTRA_Fortran_FLAGS)
    target_compile_options(UFEMISM_program PRIVATE ${EXTRA_Fortran_FLAGS})
endif()

# =============
# == OpenMPI ==
# =============

# Find MPI package
find_package(MPI REQUIRED Fortran)

# Add include directories and link libraries
target_link_libraries(UFEMISM_program PRIVATE MPI::MPI_Fortran)

# ===========
# == PETSc ==
# ===========

find_package(PkgConfig REQUIRED)
pkg_check_modules(PETSC REQUIRED PETSc)

include_directories(${PETSC_INCLUDE_DIRS})
link_directories(${PETSC_LIBRARY_DIRS})
add_definitions(${PETSC_CFLAGS_OTHER})

target_link_libraries(UFEMISM_program PRIVATE ${PETSC_LIBRARIES})

# ============
# == NetCDF ==
# ============

find_package(PkgConfig REQUIRED)
pkg_check_modules(NETCDF REQUIRED netcdf-fortran)

include_directories(${NETCDF_INCLUDE_DIRS})
link_directories(${NETCDF_LIBRARY_DIRS})
add_definitions(${NETCDF_CFLAGS_OTHER})

target_link_libraries(UFEMISM_program PRIVATE ${NETCDF_LIBRARIES})
