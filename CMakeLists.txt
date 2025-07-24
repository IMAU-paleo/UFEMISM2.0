cmake_minimum_required(VERSION 3.10)
project(UFEMISM LANGUAGES Fortran)

set(CMAKE_Fortran_STANDARD 2018)

# Collect UPSY/LADDIE/UFEMISM source files
file(GLOB_RECURSE UPSY_SOURCES    "external/UPSY/src/*.f90")
file(GLOB_RECURSE LADDIE_SOURCES  "external/LADDIE/src/*.f90")
file(GLOB_RECURSE UFEMISM_SOURCES "src/*.f90")

# Remove all *_configuration.f90 from LADDIE
foreach(src ${LADDIE_SOURCES})
    if(src MATCHES ".*/[^/]*_configuration[.]f90$")
        list(REMOVE_ITEM LADDIE_SOURCES "${src}")
    endif()
endforeach()

# Combine source files
set(ALL_SOURCES ${UPSY_SOURCES} ${LADDIE_SOURCES} ${UFEMISM_SOURCES})

# Remove all *_program.f90 except for UFEMISM_program.f90
foreach(src ${ALL_SOURCES})
    if(src MATCHES ".*/[^/]*_program[.]f90$" AND NOT src MATCHES ".*/UFEMISM_program[.]f90$")
        list(REMOVE_ITEM ALL_SOURCES "${src}")
    endif()
endforeach()

# Create executable with all sources
add_executable(UFEMISM_program ${ALL_SOURCES})

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

# Detect platform
if(APPLE)
    set(IS_MACOS TRUE)
elseif(UNIX)
    set(IS_LINUX TRUE)
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

if(IS_LINUX)
    target_link_libraries(UFEMISM_program PRIVATE ${PETSC_LIBRARIES})
elseif(IS_MACOS)
    target_link_libraries(UFEMISM_program PRIVATE ${PETSC_LIBRARY_DIRS}/libpetsc.dylib)
endif()

# ============
# == NetCDF ==
# ============

find_package(PkgConfig REQUIRED)
pkg_check_modules(NETCDF REQUIRED netcdf-fortran)

include_directories(${NETCDF_INCLUDE_DIRS})
link_directories(${NETCDF_LIBRARY_DIRS})
add_definitions(${NETCDF_CFLAGS_OTHER})

if(IS_LINUX)
    target_link_libraries(UFEMISM_program PRIVATE ${NETCDF_LIBRARIES})
elseif(IS_MACOS)
    target_link_libraries(UFEMISM_program PRIVATE ${NETCDF_LIBRARY_DIRS}/libnetcdff.dylib)
endif()
