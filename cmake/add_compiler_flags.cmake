function(add_compiler_flags target)

  # Enable preprocessing for all Fortran files
  target_compile_options(${target} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-cpp>)

  # Optionally enable assertions
  option(DO_ASSERTIONS "Compile with assertions" ON)
  if(DO_ASSERTIONS)
      target_compile_definitions(${target} PRIVATE DO_ASSERTIONS)
  endif()

  # Optionally enable resource tracking
  option(DO_RESOURCE_TRACKING "Compile with resource tracking" ON)
  if(DO_RESOURCE_TRACKING)
      target_compile_definitions(${target} PRIVATE DO_RESOURCE_TRACKING)
  endif()

  # Optionally provide extra compiler flags
  set(EXTRA_Fortran_FLAGS "" CACHE STRING "Extra gfortran compiler flags")
  if(EXTRA_Fortran_FLAGS)
      target_compile_options(${target} PRIVATE ${EXTRA_Fortran_FLAGS})
  endif()

endfunction()