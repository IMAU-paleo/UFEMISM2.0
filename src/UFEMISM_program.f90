PROGRAM UFEMISM_program
  !
  ! ===============================================================================
  ! = The main program of the Utrecht FinitE voluMe Ice Sheet Model (UFEMISM)     =
  ! =                                                                             =
  ! = Main developers:                                                            =
  ! =                                                                             =
  ! = dr. C. J. (Tijn) Berends                                                    =
  ! =   Affiliation: Institute for Marine and Atmospheric Research Utrecht (IMAU) =
  ! =   E-mail     : c.j.berends@uu.nl                                            =
  ! =                                                                             =
  ! = dr. J. A. (Jorge) Bernales                                                  =
  ! =   Affiliation: Institute for Marine and Atmospheric Research Utrecht (IMAU) =
  ! =   E-mail     : j.a.bernalesconcha@uu.nl                                     =
  ! =                                                                             =
  ! ===============================================================================
  !
  ! NOTE: the executable should be run using mpiexec to specify the number of
  !       cores n, and with the path to the config file as the only argument, e.g.:
  !
  !       mpi_exec  -n 2  UFEMISM_program  config-files/config_test

! ===== Preamble =====
! ====================

#include <petsc/finclude/petscksp.h>
  USE petscksp
  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync, initialise_parallelisation, &
                                                                     finalise_parallelisation
  USE petsc_basic                                            , ONLY: perr
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string, do_colour_strings, &
                                                                     initialise_control_and_resource_tracker, print_UFEMISM_start, print_UFEMISM_end
  USE model_configuration                                    , ONLY: C, initialise_model_configuration
  USE main_validation                                        , ONLY: run_all_unit_tests
  USE UFEMISM_main_model                                     , ONLY: type_model_region, run_model

  IMPLICIT NONE

! ===== Main variables =====
! ==========================

  ! The four model regions
  TYPE(type_model_region)                :: NAM, EAS, GRL, ANT

  ! Coupling
  REAL(dp)                               :: t_coupling, t_end_models

  ! Computation time tracking
  REAL(dp)                               :: tstart, tstop, tcomp

! ===== START =====
! =================

  ! Initialise MPI parallelisation
  CALL initialise_parallelisation
  CALL PetscInitialize( PETSC_NULL_CHARACTER, perr)

  ! Start the clock
  tstart = MPI_WTIME()

  ! Print the UFEMISM start message to the terminal
  CALL print_UFEMISM_start

  ! Initialise the control and resource tracker
  CALL initialise_control_and_resource_tracker

  ! Initialise the main model configuration
  CALL initialise_model_configuration




! ===== FINISH =====
! ==================

  ! Stop the clock
  tstop = MPI_WTIME()
  tcomp = tstop - tstart

  ! Print the UFEMISM end message to the terminal
  CALL print_UFEMISM_end( tcomp)

  ! Finalise PETSc and MPI parallelisation
  CALL sync
  CALL PetscFinalize( perr)
  CALL sync
  CALL finalise_parallelisation

END PROGRAM UFEMISM_program