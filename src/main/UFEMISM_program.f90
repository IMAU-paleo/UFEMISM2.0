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
                                                                     initialise_control_and_resource_tracker, reset_resource_tracker, &
                                                                     print_UFEMISM_start, print_UFEMISM_end
  USE model_configuration                                    , ONLY: C, initialise_model_configuration
  USE netcdf_resource_tracking                               , ONLY: create_resource_tracking_file, write_to_resource_tracking_file
  USE main_validation                                        , ONLY: run_all_unit_tests, run_all_benchmarks
  USE region_types                                           , ONLY: type_model_region
  USE UFEMISM_main_model                                     , ONLY: initialise_model_region, run_model_region

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

  ! Create the resource tracking output file
  CALL create_resource_tracking_file

  ! == Unit testing
  ! ===============

  IF (C%do_unit_tests) THEN

    ! Run all unit tests
    CALL run_all_unit_tests

    ! Write to resource tracking file
    CALL write_to_resource_tracking_file( 0._dp)

  ! == Benchmark testing
  ! ====================

  ELSEIF (c%do_benchmarks) THEN

    CALL run_all_benchmarks

    CALL write_to_resource_tracking_file( 0._dp)

  ! == Custom simulation
  ! ====================

  ELSE

    ! == Initialise the model regions
    ! ===============================

    IF (C%do_NAM) CALL initialise_model_region( NAM, 'NAM')
    IF (C%do_EAS) CALL initialise_model_region( EAS, 'EAS')
    IF (C%do_GRL) CALL initialise_model_region( GRL, 'GRL')
    IF (C%do_ANT) CALL initialise_model_region( ANT, 'ANT')

    ! == The coupling time loop
    ! =========================

    t_coupling = C%start_time_of_run

    DO WHILE (t_coupling < C%end_time_of_run)

      ! Run all model regions forward in time for one coupling interval
      t_end_models = MIN( C%end_time_of_run, t_coupling + C%dt_coupling)

      IF (C%do_NAM) CALL run_model_region( NAM, t_end_models)
      IF (C%do_EAS) CALL run_model_region( EAS, t_end_models)
      IF (C%do_GRL) CALL run_model_region( GRL, t_end_models)
      IF (C%do_ANT) CALL run_model_region( ANT, t_end_models)

      ! Advance coupling time
      t_coupling = t_end_models

      ! Write to resource tracking file
      CALL write_to_resource_tracking_file( t_coupling)
      CALL reset_resource_tracker

    END DO ! DO WHILE (t_coupling < C%end_time_of_run)

  END IF ! do_unit_test/do_benchmark/run

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
