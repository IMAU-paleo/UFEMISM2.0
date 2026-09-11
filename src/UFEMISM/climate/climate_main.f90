MODULE climate_main

  ! The main climate model module.

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  use UPSY_main, only: UPSY
  USE mpi_basic                                              , ONLY: par, sync
  USE call_stack_and_comp_time_tracking                  , ONLY: crash, init_routine, finalise_routine
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  use ice_model_data, only: atype_ice_model_data
  use ice_geometry_model_data, only: atype_ice_geometry_model_data
  use SMB_model, only: atype_SMB_model
  use SMB_snapshot_plus_anomalies, only: type_SMB_model_snp_p_anml
  use reference_geometry_types, only: type_reference_geometry
  use grid_types                                             , ONLY: type_grid
  USE climate_model_types                                    , ONLY: type_climate_model
  USE global_forcing_types                                   , ONLY: type_global_forcing
  USE climate_idealised                                      , ONLY: initialise_climate_model_idealised, run_climate_model_idealised
  USE climate_realistic                                      , ONLY: initialise_climate_model_realistic, run_climate_model_realistic, remap_climate_realistic
  USE climate_snapshot_plus_uniform_deltaT                   , ONLY: initialise_climate_model_snapshot_plus_uniform_deltaT, run_climate_model_snapshot_plus_uniform_deltaT, remap_climate_snapshot_plus_uniform_deltaT
  USE climate_snapshot_plus_transient_deltaT                 , ONLY: initialise_climate_model_snapshot_plus_transient_deltaT, run_climate_model_snapshot_plus_transient_deltaT, remap_climate_snapshot_plus_transient_deltaT
  USE climate_snapshot_plus_anomalies                        , ONLY: initialise_climate_model_snp_p_anml, run_climate_model_snp_p_anml, remap_climate_snp_p_anml
  USE reallocate_mod                                         , ONLY: reallocate_bounds
  use netcdf_io_main
  use climate_matrix                                         , only: run_climate_model_matrix, initialise_climate_matrix, remap_climate_matrix_model
  use checksum_mod, only: checksum

  IMPLICIT NONE

CONTAINS

! ===== Main routines =====
! =========================

  SUBROUTINE run_climate_model( mesh, grid, ice, geom, climate, forcing, region_name, time, SMB)
    ! Calculate the climate

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_grid),             target,    INTENT(IN)    :: grid
    class(atype_ice_model_data), target,    INTENT(IN)    :: ice
    class(atype_ice_geometry_model_data),   intent(in   ) :: geom
    TYPE(type_climate_model),    target,    INTENT(INOUT) :: climate
    TYPE(type_global_forcing),              INTENT(IN)    :: forcing
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name
    REAL(dp),                               INTENT(IN)    :: time
    class(atype_SMB_model),                 INTENT(INOUT) :: SMB

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'run_climate_model'
    CHARACTER(LEN=256)                                    :: choice_climate_model

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if we need to calculate a new climate
    IF (C%do_asynchronous_climate) THEN
      ! Asynchronous coupling: do not calculate a new climate in
      ! every model loop, but only at its own separate time step

      ! Check if this is the next climate time step
      IF (time == climate%t_next) THEN
        ! Go on to calculate a new climate
        climate%t_next = time + C%dt_climate
      ELSEIF (time > climate%t_next) THEN
        ! This should not be possible
        CALL crash('overshot the climate time step')
      ELSE
        ! It is not yet time to calculate a new climate
        CALL finalise_routine( routine_name)
        RETURN
      END IF

    ELSE ! IF (C%do_asynchronous_climate) THEN
      ! Synchronous coupling: calculate a new climate in every model loop
      climate%t_next = time + C%dt_climate
    END IF

    ! Determine which climate model to run for this region
    IF     (region_name == 'NAM') THEN
      choice_climate_model = C%choice_climate_model_NAM
    ELSEIF (region_name == 'EAS') THEN
      choice_climate_model = C%choice_climate_model_EAS
    ELSEIF (region_name == 'GRL') THEN
      choice_climate_model = C%choice_climate_model_GRL
    ELSEIF (region_name == 'ANT') THEN
      choice_climate_model = C%choice_climate_model_ANT
    ELSE
      CALL crash('unknown region_name "' // region_name // '"')
    END IF

    ! Run the chosen climate model
    SELECT CASE (choice_climate_model)
    CASE ('none')
      ! No need to do anything
    CASE ('idealised')
      CALL run_climate_model_idealised( mesh, geom, climate, time)
    CASE ('realistic')
      CALL run_climate_model_realistic( mesh, geom, climate, forcing, time)
    CASE ('snapshot_plus_uniform_deltaT')
      CALL run_climate_model_snapshot_plus_uniform_deltaT( mesh, geom, climate, time)
    CASE ('snapshot_plus_transient_deltaT')
      CALL run_climate_model_snapshot_plus_transient_deltaT( mesh, geom, climate, time)
    CASE ('snapshot_plus_anomalies')
      CALL run_climate_model_snp_p_anml( mesh, geom, climate, time)
    CASE ('matrix')
      call run_climate_model_matrix( mesh, grid, geom, SMB, climate, region_name, time, forcing)
    case ('SMB_snapshot_plus_anomalies')
      select type (snapshot_plus_anomalies => SMB)
      class default
        call crash('choice_climate_model = SMB_snapshot_plus_anomalies only works when choice_SMB_model = snapshot_plus_anomalies')
      class is (type_SMB_model_snp_p_anml)
        call snapshot_plus_anomalies%run( time, ice, geom, climate, grid)
      end select
    case ('ISMIP7')
      call climate%ISMIP7%run( geom, time)
      climate%T2m   ( mesh%vi1:mesh%vi2,:) = climate%ISMIP7%T2m   ( mesh%vi1:mesh%vi2,:)
      climate%Precip( mesh%vi1:mesh%vi2,:) = climate%ISMIP7%Precip( mesh%vi1:mesh%vi2,:)
    case ('snapshot_lapse')
      call climate%snapshot_lapse%run( geom, time)
      climate%T2m   ( mesh%vi1:mesh%vi2,:) = climate%snapshot_lapse%T2m   ( mesh%vi1:mesh%vi2,:)
      climate%Precip( mesh%vi1:mesh%vi2,:) = climate%snapshot_lapse%Precip( mesh%vi1:mesh%vi2,:)
    CASE DEFAULT
      CALL crash('unknown choice_climate_model "' // TRIM( choice_climate_model) // '"')
    END SELECT

    call checksum( mesh%pai_V, climate%T2m   , 'climate%T2m')
    call checksum( mesh%pai_V, climate%Precip, 'climate%Precip')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_climate_model

  SUBROUTINE initialise_climate_model( mesh, grid, ice, geom, climate, forcing, refgeo_PD, refgeo_init, region_name)
    ! Initialise the climate model

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    type(type_grid),                        intent(in)    :: grid
    class(atype_ice_model_data),            INTENT(IN)    :: ice
    class(atype_ice_geometry_model_data),   intent(in   ) :: geom
    TYPE(type_climate_model),               INTENT(OUT)   :: climate
    TYPE(type_global_forcing),              INTENT(IN)    :: forcing
    type(type_reference_geometry),          intent(in)    :: refgeo_PD
    type(type_reference_geometry),          intent(in)    :: refgeo_init
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'initialise_climate_model'
    CHARACTER(LEN=256)                                    :: choice_climate_model

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to terminal
    IF (par%primary)  WRITE(*,"(A)") '   Initialising climate model...'

    ! Determine which climate model to initialise for this region
    IF     (region_name == 'NAM') THEN
      choice_climate_model = C%choice_climate_model_NAM
    ELSEIF (region_name == 'EAS') THEN
      choice_climate_model = C%choice_climate_model_EAS
    ELSEIF (region_name == 'GRL') THEN
      choice_climate_model = C%choice_climate_model_GRL
    ELSEIF (region_name == 'ANT') THEN
      choice_climate_model = C%choice_climate_model_ANT
    ELSE
      CALL crash('unknown region_name "' // region_name // '"')
    END IF

    ! Allocate memory for main variables
    ALLOCATE( climate%snapshot%Hs(    mesh%vi1:mesh%vi2))
    ALLOCATE( climate%Hs(    mesh%vi1:mesh%vi2))
    ALLOCATE( climate%T2m(    mesh%vi1:mesh%vi2,12))
    ALLOCATE( climate%Precip( mesh%vi1:mesh%vi2,12))
    allocate( climate%Wind_LR( mesh%vi1:mesh%vi2, 12))
    allocate( climate%Wind_DU( mesh%vi1:mesh%vi2, 12))
    allocate( climate%Q_TOA( mesh%vi1:mesh%vi2, 12))
    climate%Hs     = 0._dp
    climate%snapshot%Hs     = 0._dp
    climate%T2m    = 0._dp
    climate%Precip = 0._dp
    climate%Wind_LR = 0._dp
    climate%Wind_DU = 0._dp
    climate%Q_TOA = 0._dp

    ! Set time of next calculation to start time
    climate%t_next = C%start_time_of_run

    ! Determine which climate model to initialise
    ! Print to terminal
    IF (par%primary)  WRITE(*,"(A)") '     Initialising climate model "' // &
      UPSY%stru%colour_string( TRIM( choice_climate_model),'light blue') // '"...'
    select case (choice_climate_model)
    case default
      call crash('unknown choice_climate_model "' // trim( choice_climate_model) // '"')
    case ('none')
      ! No need to do anything
    case ('idealised')
      call initialise_climate_model_idealised( mesh, climate)
    case ('realistic')
      call initialise_climate_model_realistic( mesh, geom, climate, forcing, region_name)
    case ('snapshot_plus_uniform_deltaT')
      call initialise_climate_model_snapshot_plus_uniform_deltaT( mesh, geom, climate, region_name)
    case ('snapshot_plus_transient_deltaT')
      call initialise_climate_model_snapshot_plus_transient_deltaT( mesh, geom, climate, region_name, C%start_time_of_run)
    CASE ('snapshot_plus_anomalies')
      CALL initialise_climate_model_snp_p_anml( mesh, geom, climate, region_name)
    case ('matrix')
      if (par%primary)  write(*,"(A)") '   Initialising climate matrix model...'
      call initialise_climate_matrix( mesh, grid, ice, geom, climate, region_name, forcing)
    case ('SMB_snapshot_plus_anomalies')
      ! No need to do anything (initialisation is handled by the SMB model)
    case ('ISMIP7')
      call climate%ISMIP7%allocate( region_name, mesh)
      call climate%ISMIP7%initialise( geom, refgeo_PD, refgeo_init, region_name)
    case ('snapshot_lapse')
      call climate%snapshot_lapse%allocate( region_name, mesh)
      call climate%snapshot_lapse%initialise( geom, refgeo_PD, refgeo_init, region_name)
      climate%Q_TOA   ( mesh%vi1:mesh%vi2,:) = climate%snapshot_lapse%Q_TOA   ( mesh%vi1:mesh%vi2,:)
      climate%Wind_LR ( mesh%vi1:mesh%vi2,:) = climate%snapshot_lapse%Wind_LR ( mesh%vi1:mesh%vi2,:)
      climate%Wind_DU ( mesh%vi1:mesh%vi2,:) = climate%snapshot_lapse%Wind_DU ( mesh%vi1:mesh%vi2,:)
    end select

    call checksum( mesh%pai_V, climate%T2m   , 'climate%T2m')
    call checksum( mesh%pai_V, climate%Precip, 'climate%Precip')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_climate_model

  SUBROUTINE write_to_restart_file_climate_model( mesh, climate, region_name, time)
    ! Write to the restart file for the climate model

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_climate_model),               INTENT(IN)    :: climate
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'write_to_restart_file_climate_model'
    CHARACTER(LEN=256)                                    :: choice_climate_model

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine which climate model is used for this region
    select case (region_name)
    case default
      call crash('unknown region_name "' // region_name // '"')
    case ('NAM')
      choice_climate_model = C%choice_climate_model_NAM
    case ('EAS')
      choice_climate_model = C%choice_climate_model_EAS
    case ('GRL')
      choice_climate_model = C%choice_climate_model_GRL
    case ('ANT')
      choice_climate_model = C%choice_climate_model_ANT
    end select

    ! Write to the restart file of the chosen climate model
    select case (choice_climate_model)
    case default
      call crash('unknown choice_climate_model "' // trim( choice_climate_model) // '"')
    case ('none', &
          'idealised', &
          'SMB_snapshot_plus_anomalies', &
          'ISMIP7', &
          'snapshot_lapse')
      ! No need to do anything
    case ('realistic', &
          'snapshot_plus_uniform_deltaT', &
          'snapshot_plus_transient_deltaT', &
          'snapshot_plus_anomalies', &
          'matrix')
      call write_to_restart_file_climate_model_region( mesh, climate, region_name, time)
    end select

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_restart_file_climate_model

  SUBROUTINE write_to_restart_file_climate_model_region( mesh, climate, region_name, time)
    ! Write to the restart NetCDF file for the climate model

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),          INTENT(IN) :: mesh
    TYPE(type_climate_model), INTENT(IN) :: climate
    CHARACTER(LEN=3),         INTENT(IN) :: region_name
    REAL(dp),                 INTENT(IN) :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER        :: routine_name = 'write_to_restart_file_climate_model_region'
    INTEGER                              :: ncid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If no NetCDF output should be created, do nothing
    IF (.NOT. C%do_create_netcdf_output) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Print to terminal
    IF (par%primary) WRITE(0,'(A)') '   Writing to climate restart file "' // &
      UPSY%stru%colour_string( TRIM( climate%restart_filename), 'light blue') // '"...'

    ! Open the NetCDF file
    CALL open_existing_netcdf_file_for_writing( climate%restart_filename, ncid)

    ! Write the time to the file
    CALL write_time_to_file( climate%restart_filename, ncid, time)

    ! ! Write the velocity fields to the file
    CALL write_to_field_multopt_mesh_dp_2D_monthly( mesh, climate%restart_filename, ncid, 'T2m',    climate%T2m)
    CALL write_to_field_multopt_mesh_dp_2D_monthly( mesh, climate%restart_filename, ncid, 'Precip', climate%Precip)

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_restart_file_climate_model_region

  SUBROUTINE create_restart_file_climate_model( mesh, climate, region_name)
    ! Create the restart file for the climate model

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_climate_model),               INTENT(INOUT) :: climate
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'create_restart_file_climate_model'
    CHARACTER(LEN=256)                                    :: choice_climate_model

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine which climate model to initialise for this region
    select case (region_name)
    case default
      call crash('unknown region_name "' // region_name // '"')
    case ('NAM')
      choice_climate_model = C%choice_climate_model_NAM
    case ('EAS')
      choice_climate_model = C%choice_climate_model_EAS
    case ('GRL')
      choice_climate_model = C%choice_climate_model_GRL
    case ('ANT')
      choice_climate_model = C%choice_climate_model_ANT
    end select

    ! Create the restart file of the chosen climate model
    select case (choice_climate_model)
    case default
      call crash('unknown choice_climate_model "' // trim( choice_climate_model) // '"')
    case ('none', &
          'idealised', &
          'SMB_snapshot_plus_anomalies', &
          'ISMIP7', &
          'snapshot_lapse')
      ! No need to do anything
    case ('realistic', &
          'snapshot_plus_uniform_deltaT', &
          'snapshot_plus_transient_deltaT', &
          'snapshot_plus_anomalies', &
          'matrix')
      call create_restart_file_climate_model_region( mesh, climate, region_name)
    end select

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_restart_file_climate_model

  SUBROUTINE create_restart_file_climate_model_region( mesh, climate, region_name)
    ! Create a restart NetCDF file for the climate submodel
    ! Includes generation of the procedural filename (e.g. "restart_climate_00001.nc")

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),          INTENT(IN)    :: mesh
    TYPE(type_climate_model), INTENT(INOUT) :: climate
    CHARACTER(LEN=3),         INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER           :: routine_name = 'create_restart_file_climate_model_region'
    CHARACTER(LEN=256)                      :: filename_base
    INTEGER                                 :: ncid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If no NetCDF output should be created, do nothing
    IF (.NOT. C%do_create_netcdf_output) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Set the filename
    filename_base = TRIM( C%output_dir) // 'restart_climate_' // region_name
    CALL generate_filename_XXXXXdotnc( filename_base, climate%restart_filename)

    ! Print to terminal
    IF (par%primary) WRITE(0,'(A)') '   Creating climate model restart file "' // &
      UPSY%stru%colour_string( TRIM( climate%restart_filename), 'light blue') // '"...'

    ! Create the NetCDF file
    CALL create_new_netcdf_file_for_writing( climate%restart_filename, ncid)

    ! Set up the mesh in the file
    CALL setup_mesh_in_netcdf_file( climate%restart_filename, ncid, mesh)

    ! Add a time dimension to the file
    CALL add_time_dimension_to_file( climate%restart_filename, ncid)

    ! Add a depth dimension to the file
    CALL add_month_dimension_to_file( climate%restart_filename, ncid)

    ! Add the data fields to the file
    CALL add_field_mesh_dp_2D_monthly( climate%restart_filename, ncid, 'T2m',    long_name = 'Near-surface air temperatures', units = 'degrees C')
    CALL add_field_mesh_dp_2D_monthly( climate%restart_filename, ncid, 'Precip', long_name = 'Precipitation rates', units = 'm/yr')

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_restart_file_climate_model_region

  SUBROUTINE remap_climate_model( mesh_old, mesh_new, climate, region_name, time, &
    refgeo_PD, refgeo_init, grid, ice, geom, forcing)
    ! Remap the climate model

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                        INTENT(IN)    :: mesh_new
    TYPE(type_climate_model),               INTENT(INOUT) :: climate
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name
    REAL(dp),                               INTENT(IN)    :: time
    type(type_reference_geometry),          intent(in)    :: refgeo_PD
    type(type_reference_geometry),          intent(in)    :: refgeo_init
    type(type_grid),                        intent(in)    :: grid
    class(atype_ice_model_data),            intent(in)    :: ice
    class(atype_ice_geometry_model_data),   intent(in   ) :: geom
    type(type_global_forcing),              intent(in   ) :: forcing

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'remap_climate_model'
    CHARACTER(LEN=256)                                    :: choice_climate_model

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to terminal
    IF (par%primary)  WRITE(*,"(A)") '    Remapping climate model data to the new mesh...'

    ! Determine which climate model to initialise for this region
    IF     (region_name == 'NAM') THEN
      choice_climate_model = C%choice_climate_model_NAM
    ELSEIF (region_name == 'EAS') THEN
      choice_climate_model = C%choice_climate_model_EAS
    ELSEIF (region_name == 'GRL') THEN
      choice_climate_model = C%choice_climate_model_GRL
    ELSEIF (region_name == 'ANT') THEN
      choice_climate_model = C%choice_climate_model_ANT
    ELSE
      CALL crash('unknown region_name "' // region_name // '"')
    END IF

    ! Reallocate memory for main variables
    CALL reallocate_bounds( climate%T2m,    mesh_new%vi1, mesh_new%vi2,12)
    CALL reallocate_bounds( climate%Precip, mesh_new%vi1, mesh_new%vi2,12)

    ! Determine which climate model to remap
    IF     (choice_climate_model == 'none') THEN
      ! No need to do anything
    ELSEIF (choice_climate_model == 'idealised') THEN
      ! No need to remap anything here
    ELSEIF (choice_climate_model == 'realistic') THEN
      call remap_climate_realistic(mesh_old, mesh_new, climate, region_name)
    ELSEIF (choice_climate_model == 'snapshot_plus_uniform_deltaT')  THEN
      call remap_climate_snapshot_plus_uniform_deltaT(mesh_old, mesh_new, geom, climate, region_name)
    ELSEIF (choice_climate_model == 'snapshot_plus_transient_deltaT')  THEN
      call remap_climate_snapshot_plus_transient_deltaT(mesh_old, mesh_new, geom, climate, region_name, time)
    ELSEIF (choice_climate_model == 'snapshot_plus_anomalies')  THEN
      call remap_climate_snp_p_anml(mesh_old, mesh_new, geom, climate, region_name, time)
    ELSEIF (choice_climate_model == 'matrix') THEN
      call remap_climate_matrix_model( mesh_new, climate, region_name, grid, ice, geom, forcing)
    elseif (choice_climate_model == 'ISMIP7') then
      call climate%ISMIP7%remap( mesh_new)
    elseif (choice_climate_model == 'snapshot_lapse') then
      call climate%snapshot_lapse%remap( mesh_new)
    ELSE
      CALL crash('unknown choice_climate_model "' // TRIM( choice_climate_model) // '"')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_climate_model

END MODULE climate_main
