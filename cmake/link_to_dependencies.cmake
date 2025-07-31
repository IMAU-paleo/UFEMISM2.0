function(link_to_dependencies target)

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
  target_link_libraries(${target} PRIVATE MPI::MPI_Fortran)

  # ===========
  # == PETSc ==
  # ===========

  find_package(PkgConfig REQUIRED)
  pkg_check_modules(PETSC REQUIRED PETSc)

  include_directories(${PETSC_INCLUDE_DIRS})
  link_directories(${PETSC_LIBRARY_DIRS})
  add_definitions(${PETSC_CFLAGS_OTHER})

  if(IS_LINUX)
      target_link_libraries(${target} PRIVATE ${PETSC_LIBRARIES})
  elseif(IS_MACOS)
      target_link_libraries(${target} PRIVATE ${PETSC_LIBRARY_DIRS}/libpetsc.dylib)
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
      target_link_libraries(${target} PRIVATE ${NETCDF_LIBRARIES})
  elseif(IS_MACOS)
      target_link_libraries(${target} PRIVATE ${NETCDF_LIBRARY_DIRS}/libnetcdff.dylib)
  endif()

endfunction()