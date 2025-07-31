MODULE GIA_ELRA

  ! ELRA module.

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE control_resources_and_error_messaging                  , ONLY: crash, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  USE GIA_model_types                                        , ONLY: type_GIA_model, type_ELRA_model
  USE region_types                                           , ONLY: type_model_region
  USE grid_basic                                             , ONLY: setup_square_grid
  USE reallocate_mod                                         , ONLY: reallocate_bounds
  use mpi_distributed_memory_grid, only: gather_gridded_data_to_primary, distribute_gridded_data_from_primary
  USE grid_types                                             , ONLY: type_grid
  USE reference_geometry_types                               , ONLY: type_reference_geometry
  use ice_geometry_basics, only: is_floating
  use remapping_main, only: map_from_mesh_vertices_to_xy_grid_2D, map_from_xy_grid_to_mesh_2D
  use kelvin_function, only: kelvin

  implicit none

  private

  public :: run_ELRA_model
  public :: calculate_ELRA_bedrock_deformation_rate
  public :: initialise_ELRA_model
  public :: remap_ELRA_model

contains

  ! The ELRA GIA model
  SUBROUTINE run_ELRA_model( region)
    ! Use the ELRA model to update bedrock elevation. Once every (dt_bedrock_ELRA) years,
    ! update deformation rates. In all other time steps, just incrementally add deformation.

    ! In/output variables:
    TYPE(type_model_region),         INTENT(INOUT)     :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_ELRA_model'
    INTEGER                                            :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate the bedrock deformation rate
    CALL calculate_ELRA_bedrock_deformation_rate( region%mesh, region%GIA%grid, region%ice, region%GIA, region%ELRA)

    ! Update bedrock with last calculated deformation rate
    DO vi = region%mesh%vi1, region%mesh%vi2
      region%ice%dHb_dt( vi) = (region%GIA%dHb_next( vi) - region%GIA%dHb_prev( vi)) / C%dt_GIA
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_ELRA_model
  SUBROUTINE calculate_ELRA_bedrock_deformation_rate( mesh, grid, ice, GIA, ELRA)
    ! Use the ELRA model to update bedrock deformation rates.

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_GIA_model),                INTENT(INOUT) :: GIA
    TYPE(type_ELRA_model),               INTENT(INOUT) :: ELRA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calculate_ELRA_bedrock_deformation_rate'
    INTEGER                                            :: vi,i,j,n,k,l,ii,jj
    REAL(dp)                                           :: Lr

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Influence radius of the lithospheric rigidity
    Lr = (C%ELRA_lithosphere_flex_rigidity / (C%ELRA_mantle_density * grav))**0.25_dp

    ! Calculate the absolute and relative surface loads on the mesh

    DO vi = mesh%vi1, mesh%vi2

      ! Absolute surface load
      IF (is_floating( ice%Hi( vi), ice%Hb( vi), ice%SL( vi))) THEN
        ELRA%surface_load_mesh( vi) = (ice%SL( vi) - ice%Hb( vi)) * grid%dx**2 * seawater_density
      ELSEIF (ice%Hi( vi) > 0._dp) THEN
        ELRA%surface_load_mesh( vi) =  ice%Hi( vi) * grid%dx**2 * ice_density
      ELSE
        ELRA%surface_load_mesh( vi) = 0._dp
      END IF

      ! Relative surface load
      GIA%relative_surface_load_mesh( vi) = ELRA%surface_load_mesh( vi) - ELRA%surface_load_GIAeq( vi)

    END DO

    ! Map relative surface load to the GIA grid
    CALL map_from_mesh_vertices_to_xy_grid_2D( mesh, grid, C%output_dir, GIA%relative_surface_load_mesh, GIA%relative_surface_load_grid)

    !! Gather data to primary
    call gather_gridded_data_to_primary( grid, GIA%relative_surface_load_grid, ELRA%relative_surface_load_grid_tot)

    n = ELRA%flex_prof_rad

    ! Let the primary do the actual work
    if (par%primary) then

      do i = 1, grid%nx
      do j = 1, grid%ny
        ELRA%dHb_eq_grid( i, j) = 0._dp
        do k = -n, n
        do l = -n, n
          ii = max( 1, min( grid%nx, i+k ))
          jj = max( 1, min( grid%ny, j+l ))
          ELRA%dHb_eq_grid( i, j) = ELRA%dHb_eq_grid( i, j) + &
            (0.5_dp * grav * Lr**2 /(pi * C%ELRA_lithosphere_flex_rigidity) * ELRA%relative_surface_load_grid_tot( ii, jj) * ELRA%flex_prof_grid( k+n+1,l+n+1))
        end do
        end do
      end do
      end do

    end if

    ! Map the actual bedrock deformation from the mesh to the grid
    call map_from_mesh_vertices_to_xy_grid_2D( mesh, grid, C%output_dir, ice%dHb, ELRA%dHb_grid_partial)

    ! gather data from all processors to primary, from partial grid vec to total 2D grid
    call gather_gridded_data_to_primary( grid, ELRA%dHb_grid_partial, ELRA%dHb_grid_tot)

	  ! Let the primary do the actual work
    if (par%primary) then

      ! Calculate the bedrock deformation rate from the difference between the current and the equilibrium deformation
      DO i = 1, grid%nx
      DO j = 1, grid%ny
        ELRA%dHb_dt_grid( i,j) = (ELRA%dHb_eq_grid( i,j) - ELRA%dHb_grid_tot( i,j)) / C%ELRA_bedrock_relaxation_time
      END DO
      END DO

    end if

    ! distribute from 2D grid data on primary to vector grid data on all processors
    call distribute_gridded_data_from_primary( grid, ELRA%dHb_dt_grid, ELRA%dHb_dt_grid_partial)

    ! remap from partial grid vec data to mesh model
    call map_from_xy_grid_to_mesh_2D( grid, mesh, C%output_dir, ELRA%dHb_dt_grid_partial, ELRA%dHb_dt_mesh)

    ! multiply the GIA time-step to calculate the bedrock deformation
    GIA%dHb_next = GIA%dHb_prev + ELRA%dHb_dt_mesh * C%dt_GIA

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calculate_ELRA_bedrock_deformation_rate

  SUBROUTINE initialise_ELRA_model( mesh, grid, ELRA, refgeo_GIAeq)
    ! Allocate and initialise the ELRA GIA model

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ELRA_model),               INTENT(INOUT) :: ELRA
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo_GIAeq

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_ELRA_model'
    INTEGER                                            :: i,j,n,k,l
    REAL(dp)                                           :: Lr, r

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (par%primary) WRITE (0,*) '    Initialising ELRA GIA model...'

    ! Allocate memory
    ALLOCATE( ELRA%surface_load_GIAeq( mesh%vi1:mesh%vi2))
    ALLOCATE( ELRA%relative_surface_load_grid_tot( grid%nx, grid%ny))
    ALLOCATE( ELRA%surface_load_mesh( mesh%vi1:mesh%vi2))
    ALLOCATE( ELRA%dHb_eq_grid( grid%nx, grid%ny))
    ALLOCATE( ELRA%dHb_grid_partial( grid%n1:grid%n2))
    ALLOCATE( ELRA%dHb_grid_tot( grid%nx, grid%ny))
    ALLOCATE( ELRA%dHb_dt_grid( grid%nx, grid%ny))
    ALLOCATE( ELRA%dHb_dt_grid_partial( grid%n1:grid%n2))
    ALLOCATE( ELRA%dHb_dt_mesh( mesh%vi1:mesh%vi2))

    ! Fill in the 2D flexural profile (= Kelvin function), with which
    ! a surface load is convoluted to find the surface deformation
    ! ============================================================

    ! Influence radius of the lithospheric rigidity
    Lr = (C%ELRA_lithosphere_flex_rigidity / (C%ELRA_mantle_density * grav))**0.25_dp

    ! Calculate radius (in number of grid cells) of the flexural profile

    IF (par%primary) THEN
      ELRA%flex_prof_rad = MIN( CEILING(grid%dx/2._dp), MAX(1, INT(6._dp * Lr / grid%dx) - 1))
      n = 2 * ELRA%flex_prof_rad + 1
      ALLOCATE( ELRA%flex_prof_grid( n, n))

      ! Calculate flexural profile
      DO i = -ELRA%flex_prof_rad, ELRA%flex_prof_rad
      DO j = -ELRA%flex_prof_rad, ELRA%flex_prof_rad
        l = i+ELRA%flex_prof_rad+1
        k = j+ELRA%flex_prof_rad+1
        r = grid%dx * SQRT( (REAL(i,dp))**2 + (REAL(j,dp))**2)
        ELRA%flex_prof_grid( l,k) = kelvin(r / Lr)
      END DO
      END DO
    END IF

    ! Calculate the reference load
    ! ===============================

    CALL initialise_ELRA_reference_load( mesh, grid, ELRA, refgeo_GIAeq)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_ELRA_model
  SUBROUTINE initialise_ELRA_reference_load( mesh, grid, ELRA, refgeo_GIAeq)

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ELRA_model),                INTENT(INOUT) :: ELRA
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo_GIAeq

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_ELRA_reference_load'
    INTEGER                                            :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate PD reference load on the mesh
    DO vi = mesh%vi1, mesh%vi2
      IF (is_floating( refgeo_GIAeq%Hi( vi), refgeo_GIAeq%Hb( vi), 0._dp)) THEN
        ELRA%surface_load_GIAeq( vi) = -refgeo_GIAeq%Hb( vi) * grid%dx**2 * seawater_density
      ELSEIF (refgeo_GIAeq%Hi( vi) > 0._dp) THEN
        ELRA%surface_load_GIAeq( vi) =  refgeo_GIAeq%Hi( vi) * grid%dx**2 * ice_density
      ELSE
        ELRA%surface_load_GIAeq( vi) = 0._dp
      END IF
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_ELRA_reference_load
  SUBROUTINE remap_ELRA_model( mesh_old, mesh_new, ELRA, refgeo_GIAeq, grid)
    ! Remap or reallocate all the data fields

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_new
    TYPE(type_ELRA_model),               INTENT(INOUT) :: ELRA
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo_GIAeq
    TYPE(type_grid),                     INTENT(IN)    :: grid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'remap_ELRA_model'
    INTEGER                                            :: int_dummy

    ! Add routine to path
    CALL init_routine( routine_name)
	CALL reallocate_bounds( ELRA%surface_load_GIAeq, mesh_new%vi1, mesh_new%vi2)
	CALL reallocate_bounds( ELRA%surface_load_mesh, mesh_new%vi1, mesh_new%vi2)

    ! Recalculate the reference load on the GIA grid
    CALL initialise_ELRA_reference_load( mesh_new, grid, ELRA, refgeo_GIAeq)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_ELRA_model

END MODULE GIA_ELRA
