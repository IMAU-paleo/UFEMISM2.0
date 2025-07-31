cmake_minimum_required(VERSION 3.10)
project( UPSY_models LANGUAGES Fortran)

set( CMAKE_Fortran_STANDARD 2018)

include( cmake/link_to_dependencies.cmake)
include( cmake/add_compiler_flags.cmake)

# Collect UPSY/LADDIE/UFEMISM source files
file( GLOB_RECURSE UPSY_sources_raw    "src/UPSY/*.f90")
file( GLOB_RECURSE LADDIE_sources_raw  "src/LADDIE/*.f90")
list(REMOVE_ITEM LADDIE_sources_raw "${CMAKE_SOURCE_DIR}/src/LADDIE/main/LADDIE_program.f90")

file( GLOB_RECURSE UFEMISM_sources_raw "src/UFEMISM/*.f90")
list(REMOVE_ITEM UFEMISM_sources_raw "${CMAKE_SOURCE_DIR}/src/UFEMISM/main/UFEMISM_program.f90")

# Create static library for UPSY
add_library( UPSY_static_library STATIC ${UPSY_sources_raw})
link_to_dependencies( UPSY_static_library)
add_compiler_flags( UPSY_static_library)

# Initialise an empty list of executables
set( all_programs "")

# UPSY programs
add_executable( UPSY_unit_test_program src/UPSY/UPSY_unit_test_program.f90)
list( APPEND all_programs UPSY_unit_test_program)

# UFEMISM programs
add_executable( UFEMISM_program ${LADDIE_sources_raw} ${UFEMISM_sources_raw} src/UFEMISM/main/UFEMISM_program.f90)
list( APPEND all_programs UFEMISM_program)

# Link all executables to the UPSY static library and the external dependencies,
# and set the compiler flags
foreach(prog IN LISTS all_programs)
    target_link_libraries(${prog} PRIVATE UPSY_static_library)
    link_to_dependencies(${prog})
    add_compiler_flags(${prog})
endforeach()
