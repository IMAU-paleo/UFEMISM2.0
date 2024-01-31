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
    CHARACTER(LEN=256)                                    :: filename_BMB_laddie_runtime
    LOGICAL                                               :: found_laddie_file

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (time > C%start_time_of_run) THEN

      ! Run LADDIE
      IF (par%master) THEN
        !CALL system('echo hello '// output_dir_IMAUICE)
        CALL system('./run_laddie_runtime.sh')
      END IF 
      CALL sync

      ! Check whether new LADDIE file is available, if not sleep
      found_laddie_file = .FALSE.
  
      DO WHILE (.NOT. found_laddie_file) ! Start sleep loop
        IF (par%master) THEN
          print*, 'Looking for LADDIE file...'
        END IF 
        CALL sync

        INQUIRE( EXIST = found_laddie_file, FILE = 'output/MISMIPplus_5km_laddie/laddieready')
  
        CALL SLEEP(1)
  
      END DO ! End sleep loop

      IF (found_laddie_file) THEN
        CALL read_field_from_file_2D( C%filename_BMB_laddie_runtime, 'BMB', mesh, BMB%BMB)
      END IF

      ! Remove laddie file
      IF (par%master) THEN
        CALL system('rm output/MISMIPplus_5km_laddie/laddieready')
      END IF 

      ! Convert to m.i.e./yr
      BMB%BMB = 31557600._dp * BMB%BMB / 918._dp

    END IF

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
    CHARACTER(LEN=256)                                    :: filename_BMB_laddie_initial
    LOGICAL                                               :: found_laddie_file

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Run LADDIE
    IF (par%master) THEN
      CALL system('./run_laddie_spinup.sh')
    END IF 
    CALL sync

    ! Check whether new LADDIE file is available, if not sleep
    found_laddie_file = .FALSE.

    DO WHILE (.NOT. found_laddie_file) ! Start sleep loop
      IF (par%master) THEN
        print*, 'Looking for LADDIE file...'
      END IF 
      CALL sync

      INQUIRE( EXIST = found_laddie_file, FILE = 'output/MISMIPplus_5km_laddie/laddieready')

      CALL SLEEP(1)

    END DO ! End sleep loop

    IF (found_laddie_file) THEN
      CALL read_field_from_file_2D( C%filename_BMB_laddie_runtime, 'BMB', mesh, BMB%BMB)
    END IF

    ! Remove laddie file
    IF (par%master) THEN
      CALL system('rm output/MISMIPplus_5km_laddie/laddieready')
    END IF 
    CALL sync

    ! Read BMB from file
    CALL read_field_from_file_2D( C%filename_BMB_laddie_initial, 'BMB', mesh, BMB%BMB)

    ! Convert to m.i.e./yr
    BMB%BMB = 31557600._dp * BMB%BMB / 918._dp

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE


END MODULE BMB_laddie