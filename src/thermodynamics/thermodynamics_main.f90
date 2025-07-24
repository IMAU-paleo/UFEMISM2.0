MODULE thermodynamics_main

  ! The main thermodynamical model module.

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE region_types                                           , ONLY: type_model_region
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  USE climate_model_types                                    , ONLY: type_climate_model
  USE SMB_model_types                                        , ONLY: type_SMB_model
  USE BMB_model_types                                        , ONLY: type_BMB_model
  use netcdf_io_main
  USE thermodynamics_3D_heat_equation                        , ONLY: solve_3D_heat_equation, create_restart_file_thermo_3D_heat_equation, &
                                                                     write_to_restart_file_thermo_3D_heat_equation
  USE thermodynamics_utilities                               , ONLY: calc_pressure_melting_point, replace_Ti_with_robin_solution, calc_homologous_temperature

  IMPLICIT NONE

CONTAINS

! ===== Main routines =====
! =========================

  SUBROUTINE run_thermodynamics_model( region)
    ! Calculate ice temperature at the desired time, and update
    ! predicted temperature if necessary

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),                INTENT(INOUT) :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'run_thermodynamics_model'
    REAL(dp)                                              :: wt_prev, wt_next
    INTEGER                                               :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If the desired time is beyond the time of the next modelled ice temperature,
    ! run the thermodynamics model to calculate a new next modelled ice temperature.
    ! ==============================================================================

    IF (region%time == region%ice%t_Ti_next) THEN
      ! Need to calculate new predicted ice temperature

      ! Store previous modelled ice temperature
      region%ice%Ti_prev = region%ice%Ti_next
      region%ice%t_Ti_prev = region%ice%t_Ti_next
      region%ice%t_Ti_next = region%ice%t_Ti_prev + C%dt_thermodynamics

      ! Run the thermodynamics model to calculate a new next modelled ice temperature.
      IF     (C%choice_thermo_model == 'none') THEN
        ! No need to do anything
      ELSEIF (C%choice_thermo_model == '3D_heat_equation') THEN
        CALL solve_3D_heat_equation( region%mesh, region%ice, region%climate, region%SMB, C%dt_thermodynamics)
      ELSE
        CALL crash('unknown choice_thermo_model "' // TRIM( C%choice_thermo_model) // '"!')
      END IF

    ELSEIF (region%time > region%ice%t_Ti_next) THEN
      ! This should not be possible
      CALL crash('overshot the thermodynamics time step')
    ELSE
      ! We're within the current thermodynamics prediction window
    END IF ! IF (region%time == region%ice%t_next) THEN

    ! Interpolate between previous and next modelled ice temperature
    ! to find the temperature at the desired time
    ! ===========================================

    ! Calculate time interpolation weights
    wt_prev = (region%ice%t_Ti_next - region%time) / (region%ice%t_Ti_next - region%ice%t_Ti_prev)
    wt_next = 1._dp - wt_prev

    ! Interpolate modelled ice temperature to desired time
    DO vi = region%mesh%vi1, region%mesh%vi2
      region%ice%Ti( vi,:) = wt_prev * region%ice%Ti_prev( vi,:) + wt_next * region%ice%Ti_next( vi,:)
    END DO

    ! Calculate Ti_pmp
    CALL calc_pressure_melting_point( region%mesh, region%ice)

    ! Calculate Ti_hom
    CALL calc_homologous_temperature( region%mesh, region%ice)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_thermodynamics_model

  SUBROUTINE initialise_thermodynamics_model( region)
    ! Initialise the thermodynamics model

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),                INTENT(INOUT) :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'initialise_thermodynamics_model'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Initialise ice temperatures
    CALL initialise_ice_temperature( region%mesh, region%ice, region%climate, region%SMB, region%name)

    ! Model states for thermodynamics model
    region%ice%t_Ti_prev = C%start_time_of_run
    region%ice%t_Ti_next = C%start_time_of_run
    region%ice%Ti_prev   = region%ice%Ti
    region%ice%Ti_next   = region%ice%Ti

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_thermodynamics_model

  SUBROUTINE write_to_restart_file_thermo( mesh, ice, time)
    ! Write to the restart NetCDF file for the thermodynamics

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_model),                INTENT(IN)              :: ice
    REAL(dp),                            INTENT(IN)              :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'write_to_restart_file_thermo'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Write to the restart file for the chosen thermodynamics model
    IF     (C%choice_thermo_model == 'none') THEN
      ! No need to do anything
    ELSEIF (C%choice_thermo_model == '3D_heat_equation') THEN
      CALL write_to_restart_file_thermo_3D_heat_equation( mesh, ice, time)
    ELSE
      CALL crash('unknown choice_thermo_model "' // TRIM( C%choice_thermo_model) // '"')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_restart_file_thermo

  SUBROUTINE create_restart_file_thermo( mesh, ice)
    ! Create a restart NetCDF file for the thermodynamics
    ! Includes generation of the procedural filename (e.g. "restart_thermo_00001.nc")

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_model),                INTENT(INOUT)           :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'create_restart_file_thermo'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Write to the restart file for the chosen thermodynamics model
    IF     (C%choice_thermo_model == 'none') THEN
      ! No need to do anything
    ELSEIF (C%choice_thermo_model == '3D_heat_equation') THEN
      CALL create_restart_file_thermo_3D_heat_equation( mesh, ice)
    ELSE
      CALL crash('unknown choice_thermo_model "' // TRIM( C%choice_thermo_model) // '"')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_restart_file_thermo

! ===== Initialise englacial temperatures =====
! =============================================

  SUBROUTINE initialise_ice_temperature( mesh, ice, climate, SMB, region_name)
    ! Initialise the englacial temperature at the start of a simulation

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                      INTENT(IN)    :: mesh
    TYPE(type_ice_model),                 INTENT(INOUT) :: ice
    TYPE(type_climate_model),             INTENT(IN)    :: climate
    TYPE(type_SMB_model),                 INTENT(IN)    :: SMB
    CHARACTER(LEN=3),                     INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                       :: routine_name = 'initialise_ice_temperature'
    CHARACTER(LEN=256)                                  :: choice_initial_ice_temperature

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine choice of initial ice temperatures for this model region
    IF     (region_name == 'NAM') THEN
      choice_initial_ice_temperature = C%choice_initial_ice_temperature_NAM
    ELSEIF (region_name == 'EAS') THEN
      choice_initial_ice_temperature = C%choice_initial_ice_temperature_EAS
    ELSEIF (region_name == 'GRL') THEN
      choice_initial_ice_temperature = C%choice_initial_ice_temperature_GRL
    ELSEIF (region_name == 'ANT') THEN
      choice_initial_ice_temperature = C%choice_initial_ice_temperature_ANT
    ELSE
      CALL crash('unknown region_name "' // region_name // '"')
    END IF

    IF     (choice_initial_ice_temperature == 'uniform') THEN
      ! Simple uniform temperature
      CALL initialise_ice_temperature_uniform( mesh, ice, region_name)
    ELSEIF (choice_initial_ice_temperature == 'linear') THEN
      ! Simple linear temperature profile
      CALL initialise_ice_temperature_linear( mesh, ice, climate)
    ELSEIF (choice_initial_ice_temperature == 'Robin') THEN
      ! Initialise with the Robin solution
      CALL initialise_ice_temperature_Robin( mesh, ice, climate, SMB)
    ELSEIF (choice_initial_ice_temperature == 'read_from_file') THEN
      ! Initialise with the temperature field from a provided NetCDF file
      CALL initialise_ice_temperature_from_file( mesh, ice, region_name)
    ELSE
      CALL crash('unknown choice_initial_ice_temperature "' // TRIM( choice_initial_ice_temperature) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_ice_temperature

  SUBROUTINE initialise_ice_temperature_uniform( mesh, ice, region_name)
    ! Initialise the englacial temperature at the start of a simulation
    !
    ! Simple uniform temperature

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                      INTENT(IN)    :: mesh
    TYPE(type_ice_model),                 INTENT(INOUT) :: ice
    CHARACTER(LEN=3),                     INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                       :: routine_name = 'initialise_ice_temperature_uniform'
    REAL(dp)                                            :: uniform_initial_ice_temperature
    INTEGER                                             :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to terminal
    IF (par%primary) WRITE (0,*) '  Initialising ice temperatures with a uniform value...'

    ! Determine choice of initial uniform ice temperatures for this model region
    IF     (region_name == 'NAM') THEN
      uniform_initial_ice_temperature = C%uniform_initial_ice_temperature_NAM
    ELSEIF (region_name == 'EAS') THEN
      uniform_initial_ice_temperature = C%uniform_initial_ice_temperature_EAS
    ELSEIF (region_name == 'GRL') THEN
      uniform_initial_ice_temperature = C%uniform_initial_ice_temperature_GRL
    ELSEIF (region_name == 'ANT') THEN
      uniform_initial_ice_temperature = C%uniform_initial_ice_temperature_ANT
    ELSE
      CALL crash('unknown region_name "' // region_name // '"')
    END IF

    DO vi = mesh%vi1, mesh%vi2
      ice%Ti( vi,:) = uniform_initial_ice_temperature
    END DO

    ! Calculate Ti_pmp
    CALL calc_pressure_melting_point( mesh, ice)

    ! Calculate Ti_hom
    CALL calc_homologous_temperature( mesh, ice)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_ice_temperature_uniform

  SUBROUTINE initialise_ice_temperature_linear( mesh, ice, climate)
    ! Initialise the englacial temperature at the start of a simulation
    !
    ! Simple linear temperature profile

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                      INTENT(IN)    :: mesh
    TYPE(type_ice_model),                 INTENT(INOUT) :: ice
    TYPE(type_climate_model),             INTENT(IN)    :: climate

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                       :: routine_name = 'initialise_ice_temperature_linear'
    INTEGER                                             :: vi,k
    REAL(dp)                                            :: T_surf_annual, T_PMP_base

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to terminal
    IF (par%primary) WRITE (0,*) '  Initialising ice temperatures with a linear profile...'

    DO vi = mesh%vi1, mesh%vi2

      T_surf_annual = MIN( T0, SUM( climate%T2m( vi,:)) / REAL( SIZE( climate%T2m( vi,:),1),dp))
      T_PMP_base    = T0 - Clausius_Clapeyron_gradient * ice%Hi_eff( vi)

      IF (ice%Hi( vi) > 0._dp) THEN
        DO k = 1, mesh%nz
          ice%Ti( vi,k) = ((1._dp - mesh%zeta( k)) * T_surf_annual) + (mesh%zeta( k) * T_PMP_base)
        END DO
      ELSE
        ice%Ti( vi,:) = T_surf_annual
      END IF

    END DO

    ! Calculate Ti_pmp
    CALL calc_pressure_melting_point( mesh, ice)

    ! Calculate Ti_hom
    CALL calc_homologous_temperature( mesh, ice)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_ice_temperature_linear

  SUBROUTINE initialise_ice_temperature_Robin( mesh, ice, climate, SMB)
    ! Initialise the englacial temperature at the start of a simulation
    !
    ! Initialise with the Robin solution

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                      INTENT(IN)    :: mesh
    TYPE(type_ice_model),                 INTENT(INOUT) :: ice
    TYPE(type_climate_model),             INTENT(IN)    :: climate
    TYPE(type_SMB_model),                 INTENT(IN)    :: SMB

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                       :: routine_name = 'initialise_ice_temperature_Robin'
    INTEGER                                             :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to terminal
    IF (par%primary) WRITE (0,*) '  Initialising ice temperatures with the Robin solution...'

    ! Calculate Ti_pmp
    CALL calc_pressure_melting_point( mesh, ice)

    ! Initialise with the Robin solution
    DO vi = mesh%vi1, mesh%vi2
      CALL replace_Ti_with_robin_solution( mesh, ice, climate, SMB, ice%Ti, vi)
    END DO

    ! Calculate Ti_hom
    CALL calc_homologous_temperature( mesh, ice)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_ice_temperature_Robin

  SUBROUTINE initialise_ice_temperature_from_file( mesh, ice, region_name)
    ! Initialise the englacial temperature at the start of a simulation
    !
    ! Initialise from an external NetCDF file

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                      INTENT(IN)    :: mesh
    TYPE(type_ice_model),                 INTENT(INOUT) :: ice
    CHARACTER(LEN=3),                     INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                       :: routine_name = 'initialise_ice_temperature_from_file'
    CHARACTER(LEN=256)                                  :: filename_initial_ice_temperature
    REAL(dp)                                            :: timeframe_initial_ice_temperature
    INTEGER                                             :: nzeta_read
    REAL(dp), DIMENSION(:    ), ALLOCATABLE             :: zeta_read
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE             :: Ti_read
    INTEGER                                             :: k

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine choice of initial uniform ice temperatures for this model region
    IF     (region_name == 'NAM') THEN
      filename_initial_ice_temperature  = C%filename_initial_ice_temperature_NAM
      timeframe_initial_ice_temperature = C%timeframe_initial_ice_temperature_NAM
    ELSEIF (region_name == 'EAS') THEN
      filename_initial_ice_temperature  = C%filename_initial_ice_temperature_EAS
      timeframe_initial_ice_temperature = C%timeframe_initial_ice_temperature_EAS
    ELSEIF (region_name == 'GRL') THEN
      filename_initial_ice_temperature  = C%filename_initial_ice_temperature_GRL
      timeframe_initial_ice_temperature = C%timeframe_initial_ice_temperature_GRL
    ELSEIF (region_name == 'ANT') THEN
      filename_initial_ice_temperature  = C%filename_initial_ice_temperature_ANT
      timeframe_initial_ice_temperature = C%timeframe_initial_ice_temperature_ANT
    ELSE
      CALL crash('unknown region_name "' // region_name // '"')
    END IF

    ! Exception for when we want to flexible read the last output file of a previous UFEMISM simulation
    if (index( filename_initial_ice_temperature,'_LAST.nc') > 1) then
      call find_last_output_file( filename_initial_ice_temperature)
      call find_last_timeframe(   filename_initial_ice_temperature, timeframe_initial_ice_temperature)
    end if

    ! Print to terminal
    IF (par%primary) WRITE (0,*) '  Initialising ice temperatures from file "' // &
      colour_string( TRIM( filename_initial_ice_temperature), 'light blue') // '"...'

    ! Allocate
    ALLOCATE( zeta_read (mesh%nz), source=0.0_dp)
    ALLOCATE( Ti_read (mesh%vi1:mesh%vi2, mesh%nz), source=0.0_dp)

    ! Read data
    IF (timeframe_initial_ice_temperature == 1E9_dp) THEN
      ! Assume the file has no time dimension
      CALL read_field_from_file_3D( filename_initial_ice_temperature, field_name_options_Ti, mesh, C%output_dir, Ti_read, &
        nzeta = nzeta_read, zeta = zeta_read)
    ELSE
      ! Read the specified timeframe from the file
      CALL read_field_from_file_3D( filename_initial_ice_temperature, field_name_options_Ti, mesh, C%output_dir, Ti_read, &
        nzeta = nzeta_read, zeta = zeta_read, time_to_read = timeframe_initial_ice_temperature)
    END IF

    ! Zeta remapping is not implemented; if the zeta dimension of the file
    ! does not match that of the model, throw an error.

    IF (nzeta_read /= mesh%nz) THEN
      CALL crash('zeta dimension in initial ice temperature file does not match that of UFEMISM!')
    END IF
    DO k = 1, mesh%nz
      IF (ABS( 1._dp - MAX( 0.001_dp, zeta_read( k)) / MAX( 0.001_dp, mesh%zeta( k))) > 1E-5_dp) THEN
        CALL crash('zeta dimension in initial ice temperature file does not match that of UFEMISM!')
      END IF
    END DO

    ! The zetas match; hurray!
    ice%Ti = Ti_read

    ! Calculate Ti_pmp
    CALL calc_pressure_melting_point( mesh, ice)

    ! Calculate Ti_hom
    CALL calc_homologous_temperature( mesh, ice)

    ! Clean up after yourself
    DEALLOCATE( zeta_read)
    DEALLOCATE( Ti_read)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_ice_temperature_from_file

END MODULE thermodynamics_main
