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

  USE petscksp
  USE precisions                                             , ONLY: dp
  use mpi_basic, only: par, initialise_parallelisation
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string, do_colour_strings, &
                                                                     initialise_control_and_resource_tracker, reset_resource_tracker, &
                                                                     print_UFEMISM_start, print_UFEMISM_end
  USE model_configuration                                    , ONLY: C, initialise_model_configuration, initialise_model_configuration_unit_tests
  use netcdf_io_main
  USE region_types                                           , ONLY: type_model_region
  USE global_forcing_types                                   , ONLY: type_global_forcing
  USE UFEMISM_main_model                                     , ONLY: initialise_model_region, run_model_region
  use global_forcings_main                                   , ONLY: initialise_global_forcings, update_global_forcings
  use inversion_utilities, only: MISMIPplus_adapt_flow_factor
  use unit_tests, only: run_all_unit_tests
  use component_tests, only: run_all_component_tests
  use unit_tests_multinode, only: run_all_multinode_unit_tests
  use checksum_mod, only: create_checksum_logfile

  IMPLICIT NONE

! ===== Main variables =====
! ==========================

  ! The four model regions
  TYPE(type_model_region)                :: NAM, EAS, GRL, ANT

  ! The global forcings
  TYPE(type_global_forcing)              :: forcing

  ! Coupling
  REAL(dp)                               :: t_coupling, t_end_models

  ! Computation time tracking
  REAL(dp)                               :: tstart, tstop, tcomp

  ! Surface elevations for the automated flow factor tuning in MISMIP+
  REAL(dp)                               :: Hs_prev, Hs_cur

  ! Input argument
  character(len=1024)                    :: input_argument

  integer :: ierr, perr

! ===== START =====
! =================

  ! Get the input argument (either the path to the config file,
  ! or an instruction to run unit/component tests)
  if (iargc() == 1) then
    call getarg( 1, input_argument)
  else
    stop 'UFEMISM requires a single argument, being the path to the config file, e.g. "mpi_exec  -n 2  UFEMISM_program  config-files/config_test"'
  end if

  ! Initialise MPI parallelisation and PETSc
  call initialise_parallelisation( input_argument)
  CALL PetscInitialize( PETSC_NULL_CHARACTER, perr)

  ! Only the primary process "sees" the input argument; all the others are
  ! initialised by MPI without it. Broadcast it so they know what to do.
  call MPI_BCAST( input_argument, len(input_argument), MPI_CHAR, 0, MPI_COMM_WORLD, ierr)

  ! Start the clock
  tstart = MPI_WTIME()

  ! Print the UFEMISM start message to the terminal
  CALL print_UFEMISM_start

  ! Initialise the control and resource tracker
  CALL initialise_control_and_resource_tracker

  ! Special cases
  if (input_argument == 'unit_tests') then
    call initialise_model_configuration_unit_tests
    call run_all_unit_tests
  elseif (input_argument == 'unit_tests_multinode') then
    call initialise_model_configuration_unit_tests
    call run_all_multinode_unit_tests
  elseif (input_argument == 'component_tests') then
    call initialise_model_configuration_unit_tests
    call run_all_component_tests
  else ! An actual model simulation

    ! Initialise the main model configuration
    CALL initialise_model_configuration

    ! Create the resource tracking output file
    CALL create_resource_tracking_file( C%output_dir)
    call create_checksum_logfile( C%output_dir)

    ! Initialise surface elevations for the automated flow factor tuning in MISMIP+
    Hs_cur = 1._dp

    ! ===== Global forcings =====
    ! ===========================
    CALL initialise_global_forcings(forcing)

    ! == Initialise the model regions
    ! ===============================

    IF (C%do_NAM) CALL initialise_model_region( NAM, 'NAM', forcing, C%start_time_of_run)
    IF (C%do_EAS) CALL initialise_model_region( EAS, 'EAS', forcing, C%start_time_of_run)
    IF (C%do_GRL) CALL initialise_model_region( GRL, 'GRL', forcing, C%start_time_of_run)
    IF (C%do_ANT) CALL initialise_model_region( ANT, 'ANT', forcing, C%start_time_of_run)

    ! == The coupling time loop
    ! =========================

    t_coupling = C%start_time_of_run

    DO WHILE (t_coupling < C%end_time_of_run)

      ! Run all model regions forward in time for one coupling interval
      t_end_models = MIN( C%end_time_of_run, t_coupling + C%dt_coupling)

      CALL update_global_forcings(forcing, t_coupling)

      IF (C%do_NAM) CALL run_model_region( NAM, t_end_models, forcing)
      IF (C%do_EAS) CALL run_model_region( EAS, t_end_models, forcing)
      IF (C%do_GRL) CALL run_model_region( GRL, t_end_models, forcing)
      IF (C%do_ANT) CALL run_model_region( ANT, t_end_models, forcing)

      ! Advance coupling time
      t_coupling = t_end_models

      ! MISMIP+ flow factor tuning for GL position
      IF (C%refgeo_idealised_MISMIPplus_tune_A) THEN
        Hs_prev = Hs_cur
        Hs_cur  = MAXVAL( ANT%ice%Hs)
        CALL MPI_ALLREDUCE( MPI_IN_PLACE, Hs_cur, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
        IF (ABS( 1._dp - Hs_cur / Hs_prev) < 5.0E-3_dp) THEN
          ! The model has converged to a steady state; adapt the flow factor
          CALL MISMIPplus_adapt_flow_factor( ANT%mesh, ANT%ice)
        END IF
      END IF

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
  call PetscFinalize( perr)
  call MPI_FINALIZE( ierr)


END PROGRAM UFEMISM_program
