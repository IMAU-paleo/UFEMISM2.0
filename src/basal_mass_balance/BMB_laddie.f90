MODULE BMB_laddie

  ! LADDIE model

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE control_resources_and_error_messaging                  , ONLY: crash, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  USE BMB_model_types                                        , ONLY: type_BMB_model
  USE netcdf_input                                           , ONLY: read_field_from_file_2D

  IMPLICIT NONE

CONTAINS

! ===== Main routines =====
! =========================

  SUBROUTINE run_BMB_model_laddie( mesh, BMB, time)
    ! Calculate the basal mass balance
    !
    ! Call the external LADDIE model

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_BMB_model),                   INTENT(INOUT) :: BMB
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'run_BMB_model_laddie'
    CHARACTER(LEN=256)                                    :: filename_BMB_laddie_output
    CHARACTER(LEN=256)                                    :: filename_BMB_laddie_configname
    CHARACTER(LEN=256)                                    :: laddieready
    LOGICAL                                               :: found_laddie_file

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (time > C%start_time_of_run) THEN

      ! Define local variables, paths to resulting melt pattern from laddie (BMB_latest.nc), and to file indicating whether laddie is ready (laddieready)
      filename_BMB_laddie_output      = TRIM(C%fixed_output_dir) // '/laddie_output/BMB_latest.nc'
      laddieready                     = TRIM(C%fixed_output_dir) // '/laddie_output/laddieready'
      filename_BMB_laddie_configname  = TRIM(C%fixed_output_dir) // '/' // TRIM(C%filename_BMB_laddie_configfile)

      ! Run LADDIE
      IF (par%master) THEN
        ! Different commands are needed to run laddie on different systems. local_mac or slurm_HPC (the latter is used for snellius)
        IF (C%choice_BMB_laddie_system == 'local_mac') THEN
          CALL system('cd ' // TRIM(C%dir_BMB_laddie_model) // '; conda activate laddie ; python3 src/runladdie.py ' // TRIM(filename_BMB_laddie_configname) // '; conda deactivate')
        ELSEIF (C%choice_BMB_laddie_system == 'slurm_HPC') THEN
          CALL system('srun --ntasks=1 --exact --overlap --cpu-bind=cores python3 src/runladdie.py ' // TRIM(filename_BMB_laddie_configname))
        ELSE
          CALL crash('C%BMB_laddie_system not recognized, should be "local_mac" or "slurm_HPC".')
        END IF
      END IF 

      ! Other cores wait for master core to finish
      CALL sync

      ! Check whether new LADDIE file is available, if not sleep
      found_laddie_file = .FALSE.
  
      DO WHILE (.NOT. found_laddie_file) ! Start sleep loop
        IF (par%master) THEN
          print*, 'Looking for LADDIE file...'
        END IF 
        CALL sync

        INQUIRE( EXIST = found_laddie_file, FILE = laddieready)
  
        CALL SLEEP(1)
  
      END DO ! End sleep loop if laddieready is found

      ! If laddieready is found, read in BMB data from LADDIE
      IF (found_laddie_file) THEN
        CALL read_field_from_file_2D( filename_BMB_laddie_output, 'BMBext', mesh, BMB%BMB_shelf)
      END IF

      IF (par%master) THEN
        CALL system('rm ' // TRIM(laddieready))
      END IF
      CALL sync

      ! Convert to m.i.e./yr
      ! EL GET VALUES FROM MODEL PARAMS
      BMB%BMB_shelf = 31557600._dp * BMB%BMB_shelf / 918._dp

    END IF

    ! EL INCLUDE CHECK FOR NAN

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE

  SUBROUTINE initialise_BMB_model_laddie( mesh, BMB)
    ! Initialise the BMB model
    !
    ! Call the external LADDIE model

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_BMB_model),                   INTENT(INOUT) :: BMB

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'initialise_BMB_model_laddie'
    CHARACTER(LEN=256)                                    :: filename_BMB_laddie_output
    CHARACTER(LEN=256)                                    :: filename_BMB_laddie_configname
    LOGICAL                                               :: found_laddie_file
    CHARACTER(LEN=256)                                    :: laddieready

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Define local variables, paths to resulting melt pattern from laddie (BMB_latest.nc), and to file indicating whether laddie is ready (laddieready)

    ! EL SET INTENT(OUT) SO ONLY NEEDED ONCE DURING INITIALISATION
    filename_BMB_laddie_output      = TRIM(C%fixed_output_dir) // '/laddie_output/BMB_latest.nc'
    laddieready                     = TRIM(C%fixed_output_dir) // '/laddie_output/laddieready'
    filename_BMB_laddie_configname  = TRIM(C%fixed_output_dir) // '/' // TRIM(C%filename_BMB_laddie_configfile)

    CALL sync

    ! Prepare config file for spinup
    IF (par%master) THEN
      ! Copy template to new file
      CALL system('cp ' // TRIM(C%filename_BMB_laddie_configtemplate) // ' ' // TRIM(filename_BMB_laddie_configname))
      ! Overwrite input files; 
      ! EL PROBABLY NEED TO ADD '' for local_mac
      CALL system ('sed -i s#@RUN_DAYS_laddie#' // TRIM(C%choice_BMB_laddie_run_days_spinup) // '#g ' // TRIM(filename_BMB_laddie_configname))
      CALL system ('sed -i s#@UFEMISM_DIRECTORY#' // TRIM(C%fixed_output_dir) // '#g ' // TRIM(filename_BMB_laddie_configname))
      CALL system ('sed -i s#@GEOMFILE#' // TRIM(C%filename_BMB_laddie_initgeometry) // '#g ' // TRIM(filename_BMB_laddie_configname))
      CALL system ('sed -i s#@T1#' // TRIM(C%choice_BMB_laddie_T1) // '#g ' // TRIM(filename_BMB_laddie_configname))
      CALL system ('sed -i s#@FROMRESTART#false#g ' // TRIM(filename_BMB_laddie_configname))
    END IF
    CALL sync

    ! Run LADDIE
    IF (par%master) THEN
      ! Different commands are needed to run laddie on different systems. local_mac or slurm_HPC (the latter is used for snellius)
      IF (C%choice_BMB_laddie_system == 'local_mac') THEN
        CALL system('cd ' // TRIM(C%dir_BMB_laddie_model) // '; conda activate laddie ; python3 runladdie.py ' // TRIM(filename_BMB_laddie_configname) // '; conda deactivate')
      ELSEIF (C%choice_BMB_laddie_system == 'slurm_HPC') THEN
        CALL system('srun --ntasks=1 --cpu-bind=cores --overlap --exact python3 src/runladdie.py ' // TRIM(filename_BMB_laddie_configname))
      ELSE
        CALL crash('C%BMB_laddie_system not recognized, should be "local_mac" or "slurm_HPC".')
      END IF
    END IF 

    ! Other cores wait for master core to finish
    CALL sync

    ! Check whether new LADDIE file is available, if not sleep
    found_laddie_file = .FALSE.

    DO WHILE (.NOT. found_laddie_file) ! Start sleep loop
      IF (par%master) THEN
        print*, 'Looking for LADDIE file...'
      END IF 
      CALL sync

      INQUIRE( EXIST = found_laddie_file, FILE = laddieready)

      CALL SLEEP(1)

    END DO ! End sleep loop

    IF (found_laddie_file) THEN
      CALL read_field_from_file_2D( filename_BMB_laddie_output, 'BMBext', mesh, BMB%BMB_shelf)
    END IF

    ! Remove laddie file
    IF (par%master) THEN
      CALL system('rm ' // TRIM(laddieready))
    END IF 
    CALL sync

    ! Convert to m.i.e./yr
    ! EL GET VALUES FROM MODEL PARAMS
    BMB%BMB_shelf = 31557600._dp * BMB%BMB_shelf / 918._dp

    ! EL INCLUDE CHECK FOR NAN

    ! Prepare config file for coupled run
    IF (par%master) THEN
      ! Copy template to new file
      CALL system('cp ' // TRIM(C%filename_BMB_laddie_configtemplate) // ' ' // TRIM(filename_BMB_laddie_configname))
      ! Overwrite input files
      CALL system ('sed -i s#@RUN_DAYS_laddie#' // TRIM(C%choice_BMB_laddie_run_days_coupled) // '#g ' // TRIM(filename_BMB_laddie_configname))
      CALL system ('sed -i s#@UFEMISM_DIRECTORY#' // TRIM(C%fixed_output_dir) // '#g ' // TRIM(filename_BMB_laddie_configname))
      CALL system ('sed -i s#@GEOMFILE#' // TRIM(C%fixed_output_dir) // '/main_output_ANT_grid.nc#g ' // TRIM(filename_BMB_laddie_configname))
      CALL system ('sed -i s#@T1#' // TRIM(C%choice_BMB_laddie_T1) // '#g ' // TRIM(filename_BMB_laddie_configname))
      CALL system ('sed -i s#@FROMRESTART#true#g ' // TRIM(filename_BMB_laddie_configname))
    END IF
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE


END MODULE BMB_laddie
