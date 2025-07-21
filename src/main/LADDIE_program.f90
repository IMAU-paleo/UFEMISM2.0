PROGRAM LADDIE_program
  !
  ! ===============================================================================
  ! = The main program of the one Layer Antarctic model for Dynamical Downscaling =
  ! = of Ice-ocean Exchanges (LADDIE)                                             =
  ! =                                                                             =
  ! = Main developers:                                                            =
  ! =                                                                             =
  ! = dr. E. (Erwin) Lambert                                                      =
  ! =   Affiliation: Royal Dutch Meteorological Institute (KNMI)                  =
  ! =   E-mail     : erwin dot lambert at knmi dot nl                             =
  ! =                                                                             =
  ! = dr. C. J. (Tijn) Berends                                                    =
  ! =   Affiliation: Institute for Marine and Atmospheric Research Utrecht (IMAU) =
  ! =   E-mail     : c.j.berends@uu.nl                                            =
  ! =                                                                             =
  ! ===============================================================================
  !
  ! NOTE: the executable should be run using mpiexec to specify the number of
  !       cores n, and with the path to the config file as the only argument, e.g.:
  !
  !       mpi_exec  -n 2  LADDIE_program  config-files/config_test

! ===== Preamble =====
! ====================

  USE petscksp
  USE precisions                                             , ONLY: dp
  use mpi_basic, only: par, initialise_parallelisation
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string, do_colour_strings, &
                                                                     initialise_control_and_resource_tracker, reset_resource_tracker, &
                                                                     print_LADDIE_start, print_LADDIE_end
  USE model_configuration                                    , ONLY: C, initialise_model_configuration, initialise_model_configuration_unit_tests
  use netcdf_io_main
  USE region_types                                           , ONLY: type_model_region
  USE mesh_types                                             , ONLY: type_mesh
  USE laddie_model_types                                     , ONLY: type_laddie_model
  USE reference_geometry_types                               , ONLY: type_reference_geometry
  USE global_forcing_types                                   , ONLY: type_global_forcing
  USE LADDIE_main_model                                      , ONLY: initialise_model_region, run_model_region
  use unit_tests                                             , only: run_laddie_unit_tests

  IMPLICIT NONE

! ===== Main variables =====
! ==========================

  ! The model region
  TYPE(type_model_region)                :: ANT

  ! The global forcings
  TYPE(type_global_forcing)              :: forcing

  type(type_laddie_model)                :: laddie
  type(type_reference_geometry)          :: refgeo
  type(type_mesh)                        :: mesh

  ! Coupling
  REAL(dp)                               :: t_coupling, t_end_models

  ! Computation time tracking
  REAL(dp)                               :: tstart, tstop, tcomp

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
    stop 'LADDIE requires a single argument, being the path to the config file, e.g. "mpi_exec  -n 2  LADDIE_program  config-files/config_test"'
  end if

  ! Initialise MPI parallelisation and PETSc
  call initialise_parallelisation( input_argument)
  CALL PetscInitialize( PETSC_NULL_CHARACTER, perr)

  ! Only the primary process "sees" the input argument; all the others are
  ! initialised by MPI without it. Broadcast it so they know what to do.
  call MPI_BCAST( input_argument, len(input_argument), MPI_CHAR, 0, MPI_COMM_WORLD, ierr)

  ! Start the clock
  tstart = MPI_WTIME()

  ! Print the LADDIE start message to the terminal
  CALL print_LADDIE_start

  ! Initialise the control and resource tracker
  CALL initialise_control_and_resource_tracker

  ! Special cases
  if (input_argument == 'unit_tests') then
    call initialise_model_configuration_unit_tests
    call run_laddie_unit_tests
  else ! An actual model simulation

    ! Initialise the main model configuration
    CALL initialise_model_configuration

    ! Create the resource tracking output file
    CALL create_resource_tracking_file( C%output_dir)

    ! == Initialise the model regions
    ! ===============================

    call initialise_model_region( ANT, 'ANT', laddie, refgeo, mesh, forcing, C%start_time_of_run)

    ! == The coupling time loop
    ! =========================

    t_coupling = C%start_time_of_run

    DO WHILE (t_coupling < C%end_time_of_run)

      ! Run all model regions forward in time for one coupling interval
      t_end_models = MIN( C%end_time_of_run, t_coupling + C%dt_coupling)

      CALL run_model_region( ANT, t_end_models, forcing)

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

  ! Print the LADDIE end message to the terminal
  CALL print_LADDIE_end( tcomp)

  ! Finalise PETSc and MPI parallelisation
  call PetscFinalize( perr)
  call MPI_FINALIZE( ierr)


END PROGRAM LADDIE_program
