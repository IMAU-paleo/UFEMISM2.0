MODULE mesh_creation

  ! Routines used to create a mesh.

! ===== Preamble =====
! ====================

  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE mesh_memory                                            , ONLY: allocate_mesh_primary
  USE mesh_utilities                                         , ONLY: update_triangle_circumcenter, calc_mesh_contour_as_line, calc_mesh_mask_as_polygons
  USE mesh_refinement                                        , ONLY: refine_mesh_uniform, refine_mesh_line, Lloyds_algorithm_single_iteration, &
                                                                     refine_mesh_polygon, refine_mesh_line_ROI, refine_mesh_polygon_ROI, &
                                                                     calc_polygon_Pine_Island_Glacier, calc_polygon_Thwaites_Glacier, &
                                                                     calc_polygon_Amery_ice_shelf, calc_polygon_Riiser_Larsen_ice_shelf, &
                                                                     calc_polygon_Siple_Coast, calc_polygon_Patagonia, calc_polygon_Larsen_ice_shelf, &
                                                                     calc_polygon_Transantarctic_Mountains, calc_polygon_Narsarsuaq, &
                                                                     calc_polygon_Tijn_test_ISMIP_HOM_A, &
                                                                     enforce_contiguous_process_domains
  USe mesh_parallel_creation                                 , ONLY: broadcast_mesh
  USE mesh_secondary                                         , ONLY: calc_all_secondary_mesh_data
  USE mesh_operators                                         , ONLY: calc_all_matrix_operators_mesh
  USE grid_basic                                             , ONLY: type_grid, gather_gridded_data_to_master_dp_2D, calc_grid_contour_as_line, &
                                                                     calc_grid_mask_as_polygons
  USE math_utilities                                         , ONLY: thickness_above_floatation
  USE mpi_distributed_memory                                 , ONLY: gather_to_master_dp_1D

  IMPLICIT NONE

CONTAINS

! ===== Subroutines =====
! =======================

! == The two main mesh creation routines that should be called
! ============================================================

  ! Create a mesh from ice geometry on a grid
  SUBROUTINE create_mesh_from_gridded_geometry( region_name, name, grid, Hi, Hb, Hs, SL, xmin, xmax, ymin, ymax, lambda_M, phi_M, beta_stereo, mesh)
    ! Create a mesh from ice geometry on a grid

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=3),           INTENT(IN)        :: region_name
    CHARACTER(LEN=256),         INTENT(IN)        :: name
    TYPE(type_grid),            INTENT(IN)        :: grid
    REAL(dp), DIMENSION(:    ), INTENT(IN)        :: Hi, Hb, Hs, SL
    REAL(dp),                   INTENT(IN)        :: xmin, xmax, ymin, ymax
    REAL(dp),                   INTENT(IN)        :: lambda_M, phi_M, beta_stereo
    TYPE(type_mesh),            INTENT(OUT)       :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'create_mesh_from_gridded_geometry'
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE       :: poly_mult_sheet
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE       :: poly_mult_shelf
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE       :: p_line_grounding_line
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE       :: p_line_calving_front
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE       :: p_line_ice_front
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE       :: p_line_coastline

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Reduce the gridded ice geometry to lines and polygons
    CALL reduce_gridded_ice_geometry( grid, Hi, Hb, Hs, SL, poly_mult_sheet, poly_mult_shelf, &
      p_line_grounding_line, p_line_calving_front, p_line_ice_front, p_line_coastline)

    ! Create a mesh from the reduced ice geometry
    CALL create_mesh_from_reduced_geometry( region_name, name, poly_mult_sheet, poly_mult_shelf, &
      p_line_grounding_line, p_line_calving_front, p_line_ice_front, p_line_coastline, &
      xmin, xmax, ymin, ymax, lambda_M, phi_M, beta_stereo, mesh)

    ! Clean up after yourself
    DEALLOCATE( poly_mult_sheet)
    DEALLOCATE( poly_mult_shelf)
    DEALLOCATE( p_line_grounding_line)
    DEALLOCATE( p_line_calving_front)
    DEALLOCATE( p_line_ice_front)
    DEALLOCATE( p_line_coastline)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_mesh_from_gridded_geometry

  ! Create a mesh from ice geometry on a mesh
  SUBROUTINE create_mesh_from_meshed_geometry( region_name, name, mesh_src, Hi, Hb, Hs, SL, xmin, xmax, ymin, ymax, lambda_M, phi_M, beta_stereo, mesh)
    ! Create a mesh from ice geometry on a mesh

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=3),           INTENT(IN)        :: region_name
    CHARACTER(LEN=256),         INTENT(IN)        :: name
    TYPE(type_mesh),            INTENT(IN)        :: mesh_src
    REAL(dp), DIMENSION(:    ), INTENT(IN)        :: Hi, Hb, Hs, SL
    REAL(dp),                   INTENT(IN)        :: xmin, xmax, ymin, ymax
    REAL(dp),                   INTENT(IN)        :: lambda_M, phi_M, beta_stereo
    TYPE(type_mesh),            INTENT(OUT)       :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'create_mesh_from_meshed_geometry'
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE       :: poly_mult_sheet
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE       :: poly_mult_shelf
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE       :: p_line_grounding_line
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE       :: p_line_calving_front
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE       :: p_line_ice_front
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE       :: p_line_coastline

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Reduce the meshed ice geometry to lines and polygons
    CALL reduce_meshed_ice_geometry( mesh_src, Hi, Hb, Hs, SL, poly_mult_sheet, poly_mult_shelf, &
      p_line_grounding_line, p_line_calving_front, p_line_ice_front, p_line_coastline)

    ! Create a mesh from the reduced ice geometry
    CALL create_mesh_from_reduced_geometry( region_name, name, poly_mult_sheet, poly_mult_shelf, &
      p_line_grounding_line, p_line_calving_front, p_line_ice_front, p_line_coastline, &
      xmin, xmax, ymin, ymax, lambda_M, phi_M, beta_stereo, mesh)

    ! Clean up after yourself
    DEALLOCATE( poly_mult_sheet)
    DEALLOCATE( poly_mult_shelf)
    DEALLOCATE( p_line_grounding_line)
    DEALLOCATE( p_line_calving_front)
    DEALLOCATE( p_line_ice_front)
    DEALLOCATE( p_line_coastline)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_mesh_from_meshed_geometry

! == Reduce a gridded/meshed ice geometry to grid-less polygons and lines
! =======================================================================

  SUBROUTINE reduce_gridded_ice_geometry( grid, Hi, Hb, Hs, SL, poly_mult_sheet, poly_mult_shelf, &
    p_line_grounding_line, p_line_calving_front, p_line_ice_front, p_line_coastline)

    ! "Reduce" the gridded ice-sheet geometry to a set of polygons
    ! (describing regions covered by grounded/floating ice) and lines
    ! (grounding line, calving front, etc.)
    !
    ! NOTE: gridded geometry fields Hi, Hb, Hs, SL are provided in distributed vector form
    !
    ! NOTE: polygons and lines are returned in full to all processes

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                         INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:    ),              INTENT(IN)    :: Hi, Hb, Hs, SL
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE, INTENT(OUT)   :: poly_mult_sheet
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE, INTENT(OUT)   :: poly_mult_shelf
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE, INTENT(OUT)   :: p_line_grounding_line
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE, INTENT(OUT)   :: p_line_calving_front
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE, INTENT(OUT)   :: p_line_ice_front
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE, INTENT(OUT)   :: p_line_coastline

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'reduce_gridded_ice_geometry'
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE       :: Hi_grid
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE       :: Hb_grid
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE       :: Hs_grid
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE       :: SL_grid
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE       :: TAF_grid
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE       :: Hb_minus_SL_grid
    INTEGER                                       :: i,j,ii,jj
    LOGICAL,  DIMENSION(:,:  ), ALLOCATABLE       :: mask_sheet
    LOGICAL,  DIMENSION(:,:  ), ALLOCATABLE       :: mask_shelf
    LOGICAL,  DIMENSION(:,:  ), ALLOCATABLE       :: mask_calc_grounding_line
    LOGICAL,  DIMENSION(:,:  ), ALLOCATABLE       :: mask_calc_calving_front
    LOGICAL,  DIMENSION(:,:  ), ALLOCATABLE       :: mask_calc_ice_front
    LOGICAL,  DIMENSION(:,:  ), ALLOCATABLE       :: mask_calc_coastline
    INTEGER                                       :: n_poly_mult_sheet
    INTEGER                                       :: n_poly_mult_shelf
    INTEGER                                       :: n_line_grounding_line
    INTEGER                                       :: n_line_calving_front
    INTEGER                                       :: n_line_ice_front
    INTEGER                                       :: n_line_coastline

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate memory for gathered ice geometry data
    IF (par%master) THEN
      ALLOCATE( Hi_grid( grid%nx, grid%ny))
      ALLOCATE( Hb_grid( grid%nx, grid%ny))
      ALLOCATE( Hs_grid( grid%nx, grid%ny))
      ALLOCATE( SL_grid( grid%nx, grid%ny))
    else
      ALLOCATE( Hi_grid( 0, 0))
      ALLOCATE( Hb_grid( 0, 0))
      ALLOCATE( Hs_grid( 0, 0))
      ALLOCATE( SL_grid( 0, 0))
    END IF

    ! Gather ice geometry data in grid form to the master
    CALL gather_gridded_data_to_master_dp_2D( grid, Hi, Hi_grid)
    CALL gather_gridded_data_to_master_dp_2D( grid, Hb, Hb_grid)
    CALL gather_gridded_data_to_master_dp_2D( grid, Hs, Hs_grid)
    CALL gather_gridded_data_to_master_dp_2D( grid, SL, SL_grid)

    ! Let the master calculate the polygons and lines

    IF (par%master) THEN

      ! Allocate memory for thickness above floatation and Hb-SL
      ALLOCATE( TAF_grid(         grid%nx, grid%ny))
      ALLOCATE( Hb_minus_SL_grid( grid%nx, grid%ny))

      ! Calculate thickness above floatation and Hb-SL
      DO i = 1, grid%nx
      DO j = 1, grid%ny
        TAF_grid(         i,j) = thickness_above_floatation( Hi_grid( i,j), Hb_grid( i,j), SL_grid( i,j))
        Hb_minus_SL_grid( i,j) = Hb_grid( i,j) - SL_grid( i,j)
      END DO
      END DO

      ! Fill in masks for floating/grounded ice
      ALLOCATE( mask_sheet( grid%nx, grid%ny), source = .FALSE.)
      ALLOCATE( mask_shelf( grid%nx, grid%ny), source = .FALSE.)

      DO i = 1, grid%nx
      DO j = 1, grid%ny
        IF (Hi_grid( i,j) > 0.1_dp) THEN
          IF (TAF_grid( i,j) > 0._dp) THEN
            mask_sheet( i,j) = .TRUE.
          ELSE
            mask_shelf( i,j) = .TRUE.
          END IF
        END IF
      END DO
      END DO

      ! Fill in masks for where to calculate lines
      ALLOCATE( mask_calc_grounding_line( grid%nx, grid%ny), source = .FALSE.)
      ALLOCATE( mask_calc_calving_front(  grid%nx, grid%ny), source = .FALSE.)
      ALLOCATE( mask_calc_ice_front(      grid%nx, grid%ny), source = .FALSE.)
      ALLOCATE( mask_calc_coastline(      grid%nx, grid%ny), source = .FALSE.)

      DO i = 1, grid%nx
      DO j = 1, grid%ny

        ! Grounding line should only be calculated where there's ice
        IF (Hi_grid( i,j) > 0.1_dp) THEN
          mask_calc_grounding_line( i,j) = .TRUE.
        END IF

        ! Calving front should only be calculated for ice (both floating and grounded) next to ocean
        IF (Hi_grid( i,j) > 0.1_dp) THEN
          ! This grid cell has ice
          DO ii = MAX( 1, i-1), MIN( grid%nx, i+1)
          DO jj = MAX( 1, j-1), MIN( grid%ny, j+1)
            IF (Hi_grid( ii,jj) < 0.1_dp .AND. Hb_grid( ii,jj) < SL_grid( ii,jj)) THEN
              ! This neighbour is open ocean
              mask_calc_calving_front( i,j) = .TRUE.
            END IF
          END DO
          END DO
        END IF

        ! Ice front is simply ice next to non-ice
        IF (Hi_grid( i,j) > 0.1_dp) THEN
          ! This grid cell has ice
          DO ii = MAX( 1, i-1), MIN( grid%nx, i+1)
          DO jj = MAX( 1, j-1), MIN( grid%ny, j+1)
            IF (Hi_grid( ii,jj) < 0.1_dp) THEN
              ! This neighbour is ice-free
              mask_calc_ice_front( i,j) = .TRUE.
            END IF
          END DO
          END DO
        END IF

        ! Coastline is ice-free land next to open ocean
        IF (Hi_grid( i,j) < 0.1_dp .AND. Hb_grid( i,j) > SL_grid( i,j)) THEN
          ! This grid cell is ice-free land
          DO ii = MAX( 1, i-1), MIN( grid%nx, i+1)
          DO jj = MAX( 1, j-1), MIN( grid%ny, j+1)
            IF (Hi_grid( ii,jj) < 0.1_dp .AND. Hb_grid( ii,jj) < SL_grid( ii,jj)) THEN
              ! This neighbour is open ocean
              mask_calc_coastline( i,j) = .TRUE.
            END IF
          END DO
          END DO
        END IF

      END DO
      END DO

      ! Calculate polygons enveloping sheet and shelf
      CALL calc_grid_mask_as_polygons( grid, mask_sheet, poly_mult_sheet)
      CALL calc_grid_mask_as_polygons( grid, mask_shelf, poly_mult_shelf)

      ! Get polygon sizes
      n_poly_mult_sheet = SIZE( poly_mult_sheet,1)
      n_poly_mult_shelf = SIZE( poly_mult_shelf,1)

      ! Calculate lines in line segment format
      CALL calc_grid_contour_as_line( grid, TAF_grid        , 0.0_dp, p_line_grounding_line, mask_calc_grounding_line)
      CALL calc_grid_contour_as_line( grid, Hi_grid         , 0.1_dp, p_line_calving_front , mask_calc_calving_front )
      CALL calc_grid_contour_as_line( grid, Hi_grid         , 0.1_dp, p_line_ice_front     , mask_calc_ice_front     )
      CALL calc_grid_contour_as_line( grid, Hb_minus_SL_grid, 0.0_dp, p_line_coastline     , mask_calc_coastline     )

      ! Get line sizes
      n_line_grounding_line = SIZE( p_line_grounding_line,1)
      n_line_calving_front  = SIZE( p_line_calving_front ,1)
      n_line_ice_front      = SIZE( p_line_ice_front     ,1)
      n_line_coastline      = SIZE( p_line_coastline     ,1)

    END IF ! IF (par%master) THEN
    CALL sync

    ! Broadcast polygon sizes to all processes
    CALL MPI_BCAST( n_poly_mult_sheet, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( n_poly_mult_shelf, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Allocate memory for the polygons on the other processes
    IF (.NOT. par%master) THEN
      ALLOCATE( poly_mult_sheet( n_poly_mult_sheet,2))
      ALLOCATE( poly_mult_shelf( n_poly_mult_shelf,2))
    END IF

    ! Broadcast polygons to all processes
    CALL MPI_BCAST( poly_mult_sheet, n_poly_mult_sheet*2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( poly_mult_shelf, n_poly_mult_shelf*2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! Broadcast line sizes to all processes
    CALL MPI_BCAST( n_line_grounding_line, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( n_line_calving_front , 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( n_line_ice_front     , 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( n_line_coastline     , 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Allocate memory for the lines on the other processes
    IF (.NOT. par%master) THEN
      ALLOCATE( p_line_grounding_line( n_line_grounding_line,4))
      ALLOCATE( p_line_calving_front(  n_line_calving_front ,4))
      ALLOCATE( p_line_ice_front(      n_line_ice_front     ,4))
      ALLOCATE( p_line_coastline(      n_line_coastline     ,4))
    END IF

    ! Broadcast lines to all processes
    CALL MPI_BCAST( p_line_grounding_line, n_line_grounding_line*4, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( p_line_calving_front , n_line_calving_front *4, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( p_line_ice_front     , n_line_ice_front     *4, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( p_line_coastline     , n_line_coastline     *4, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! Clean up after yourself
    DEALLOCATE( Hi_grid)
    DEALLOCATE( Hb_grid)
    DEALLOCATE( Hs_grid)
    DEALLOCATE( SL_grid)
    IF (par%master) THEN
      DEALLOCATE( TAF_grid)
      DEALLOCATE( Hb_minus_SL_grid)
      DEALLOCATE( mask_sheet)
      DEALLOCATE( mask_shelf)
      DEALLOCATE( mask_calc_grounding_line)
      DEALLOCATE( mask_calc_calving_front)
      DEALLOCATE( mask_calc_ice_front)
      DEALLOCATE( mask_calc_coastline)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE reduce_gridded_ice_geometry

  SUBROUTINE reduce_meshed_ice_geometry( mesh, Hi, Hb, Hs, SL, poly_mult_sheet, poly_mult_shelf, &
    p_line_grounding_line, p_line_calving_front, p_line_ice_front, p_line_coastline)

    ! "Reduce" the meshed ice-sheet geometry to a set of polygons
    ! (describing regions covered by grounded/floating ice) and lines
    ! (grounding line, calving front, etc.)
    !
    ! NOTE: meshed geometry fields Hi, Hb, Hs, SL are provided in distributed vector form
    !
    ! NOTE: polygons and lines are returned in full to all processes

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                         INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),              INTENT(IN)    :: Hi, Hb, Hs, SL
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE, INTENT(OUT)   :: poly_mult_sheet
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE, INTENT(OUT)   :: poly_mult_shelf
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE, INTENT(OUT)   :: p_line_grounding_line
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE, INTENT(OUT)   :: p_line_calving_front
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE, INTENT(OUT)   :: p_line_ice_front
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE, INTENT(OUT)   :: p_line_coastline

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'reduce_meshed_ice_geometry'
    REAL(dp), DIMENSION(:    ), ALLOCATABLE       :: Hi_tot
    REAL(dp), DIMENSION(:    ), ALLOCATABLE       :: Hb_tot
    REAL(dp), DIMENSION(:    ), ALLOCATABLE       :: Hs_tot
    REAL(dp), DIMENSION(:    ), ALLOCATABLE       :: SL_tot
    REAL(dp), DIMENSION(:    ), ALLOCATABLE       :: TAF_tot
    REAL(dp), DIMENSION(:    ), ALLOCATABLE       :: Hb_minus_SL_tot
    INTEGER                                       :: vi
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE       :: mask_sheet
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE       :: mask_shelf
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE       :: mask_calc_grounding_line
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE       :: mask_calc_calving_front
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE       :: mask_calc_coastline
    INTEGER                                       :: n_poly_mult_sheet
    INTEGER                                       :: n_poly_mult_shelf
    INTEGER                                       :: n_line_grounding_line
    INTEGER                                       :: n_line_calving_front
    INTEGER                                       :: n_line_ice_front
    INTEGER                                       :: n_line_coastline

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate memory for gathered ice geometry data
    IF (par%master) THEN
      ALLOCATE( Hi_tot( mesh%nV))
      ALLOCATE( Hb_tot( mesh%nV))
      ALLOCATE( Hs_tot( mesh%nV))
      ALLOCATE( SL_tot( mesh%nV))
    END IF

    ! Gather ice geometry data in grid form to the master
    CALL gather_to_master_dp_1D( Hi, Hi_tot)
    CALL gather_to_master_dp_1D( Hb, Hb_tot)
    CALL gather_to_master_dp_1D( Hs, Hs_tot)
    CALL gather_to_master_dp_1D( SL, SL_tot)

    ! Let the master calculate the polygons and lines

    IF (par%master) THEN

      ! Allocate memory for thickness above floatation and Hb-SL
      ALLOCATE( TAF_tot(         mesh%nV))
      ALLOCATE( Hb_minus_SL_tot( mesh%nV))

      ! Calculate thickness above floatation and Hb-SL
      DO vi = 1, mesh%nV
        TAF_tot(         vi) = thickness_above_floatation( Hi_tot( vi), Hb_tot( vi), SL_tot( vi))
        Hb_minus_SL_tot( vi) = Hb_tot( vi) - SL_tot( vi)
      END DO

      ! Fill in masks for floating/grounded ice
      ALLOCATE( mask_sheet( mesh%nV), source = .FALSE.)
      ALLOCATE( mask_shelf( mesh%nV), source = .FALSE.)

      DO vi = 1, mesh%nV
        IF (Hi_tot( vi) > 0.1_dp) THEN
          IF (TAF_tot( vi) > 0._dp) THEN
            mask_sheet( vi) = .TRUE.
          ELSE
            mask_shelf( vi) = .TRUE.
          END IF
        END IF
      END DO

      ! Fill in masks for where to calculate lines
      ALLOCATE( mask_calc_grounding_line( mesh%nV), source = .FALSE.)
      ALLOCATE( mask_calc_calving_front ( mesh%nV), source = .FALSE.)
      ALLOCATE( mask_calc_coastline     ( mesh%nV), source = .FALSE.)

      DO vi = 1, mesh%nV

        ! Grounding line should only be calculated where there's ice
        IF (Hi_tot( vi) > 0.1_dp) THEN
          mask_calc_grounding_line( vi) = .TRUE.
        END IF

        ! Calving front should only be calculated on floating ice
        IF (TAF_tot( vi) < 0.1_dp .AND. Hi_tot( vi) > 0.1_dp) THEN
           mask_calc_calving_front( vi) = .TRUE.
        END IF

        ! Coastline should only be calculated where there's no ice
        IF (Hi_tot( vi) < 0.1_dp ) THEN
           mask_calc_coastline( vi) = .TRUE.
        END IF

      END DO ! DO vi = 1, mesh%nV

      ! Calculate polygons enveloping sheet and shelf
      CALL calc_mesh_mask_as_polygons( mesh, mask_sheet, poly_mult_sheet)
      CALL calc_mesh_mask_as_polygons( mesh, mask_shelf, poly_mult_shelf)

      ! Get polygon sizes
      n_poly_mult_sheet = SIZE( poly_mult_sheet,1)
      n_poly_mult_shelf = SIZE( poly_mult_shelf,1)

      ! Calculate lines in line segment format
      CALL calc_mesh_contour_as_line( mesh, TAF_tot        , 0.0_dp, p_line_grounding_line, mask_calc_grounding_line)
      CALL calc_mesh_contour_as_line( mesh, Hi_tot         , 1.0_dp, p_line_calving_front , mask_calc_calving_front )
      CALL calc_mesh_contour_as_line( mesh, Hi_tot         , 1.0_dp, p_line_ice_front                               )
      CALL calc_mesh_contour_as_line( mesh, Hb_minus_SL_tot, 0.0_dp, p_line_coastline     , mask_calc_coastline     )

      ! Get line sizes
      n_line_grounding_line = SIZE( p_line_grounding_line,1)
      n_line_calving_front  = SIZE( p_line_calving_front ,1)
      n_line_ice_front      = SIZE( p_line_ice_front     ,1)
      n_line_coastline      = SIZE( p_line_coastline     ,1)

    END IF ! IF (par%master) THEN
    CALL sync

    ! Broadcast polygon sizes to all processes
    CALL MPI_BCAST( n_poly_mult_sheet, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( n_poly_mult_shelf, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Allocate memory for the polygons on the other processes
    IF (.NOT. par%master) THEN
      ALLOCATE( poly_mult_sheet( n_poly_mult_sheet,2))
      ALLOCATE( poly_mult_shelf( n_poly_mult_shelf,2))
    END IF

    ! Broadcast polygons to all processes
    CALL MPI_BCAST( poly_mult_sheet, n_poly_mult_sheet*2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( poly_mult_shelf, n_poly_mult_shelf*2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! Broadcast line sizes to all processes
    CALL MPI_BCAST( n_line_grounding_line, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( n_line_calving_front , 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( n_line_ice_front     , 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( n_line_coastline     , 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Allocate memory for the lines on the other processes
    IF (.NOT. par%master) THEN
      ALLOCATE( p_line_grounding_line( n_line_grounding_line,4))
      ALLOCATE( p_line_calving_front(  n_line_calving_front ,4))
      ALLOCATE( p_line_ice_front(      n_line_ice_front     ,4))
      ALLOCATE( p_line_coastline(      n_line_coastline     ,4))
    END IF

    ! Broadcast lines to all processes
    CALL MPI_BCAST( p_line_grounding_line, n_line_grounding_line*4, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( p_line_calving_front , n_line_calving_front *4, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( p_line_ice_front     , n_line_ice_front     *4, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( p_line_coastline     , n_line_coastline     *4, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! Clean up after yourself
    IF (par%master) THEN
      DEALLOCATE( Hi_tot)
      DEALLOCATE( Hb_tot)
      DEALLOCATE( Hs_tot)
      DEALLOCATE( SL_tot)
      DEALLOCATE( TAF_tot)
      DEALLOCATE( Hb_minus_SL_tot)
      DEALLOCATE( mask_sheet)
      DEALLOCATE( mask_shelf)
      DEALLOCATE( mask_calc_grounding_line)
      DEALLOCATE( mask_calc_calving_front)
      DEALLOCATE( mask_calc_coastline)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE reduce_meshed_ice_geometry

! == Create a mesh from reduced ice geometry
! ==========================================

  SUBROUTINE create_mesh_from_reduced_geometry( region_name, name, poly_mult_sheet, poly_mult_shelf, &
    p_line_grounding_line, p_line_calving_front, p_line_ice_front, p_line_coastline, &
    xmin, xmax, ymin, ymax, lambda_M, phi_M, beta_stereo, mesh)
    ! Create mesh from the reduced ice geometry

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=3),           INTENT(IN)        :: region_name
    CHARACTER(LEN=256),         INTENT(IN)        :: name
    REAL(dp), DIMENSION(:,:  ), INTENT(IN)        :: poly_mult_sheet
    REAL(dp), DIMENSION(:,:  ), INTENT(IN)        :: poly_mult_shelf
    REAL(dp), DIMENSION(:,:  ), INTENT(IN)        :: p_line_grounding_line
    REAL(dp), DIMENSION(:,:  ), INTENT(IN)        :: p_line_calving_front
    REAL(dp), DIMENSION(:,:  ), INTENT(IN)        :: p_line_ice_front
    REAL(dp), DIMENSION(:,:  ), INTENT(IN)        :: p_line_coastline
    REAL(dp),                   INTENT(IN)        :: xmin, xmax, ymin, ymax
    REAL(dp),                   INTENT(IN)        :: lambda_M, phi_M, beta_stereo
    TYPE(type_mesh),            INTENT(OUT)       :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'create_mesh_from_reduced_geometry'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Choose single-core or parallelised version
    IF (C%do_singlecore_mesh_creation) THEN
      CALL create_mesh_from_reduced_geometry_singlecore( region_name, name, poly_mult_sheet, poly_mult_shelf, &
        p_line_grounding_line, p_line_calving_front, p_line_ice_front, p_line_coastline, &
        xmin, xmax, ymin, ymax, lambda_M, phi_M, beta_stereo, mesh)
    ELSE
      CALL create_mesh_from_reduced_geometry_parallelised( region_name, name, poly_mult_sheet, poly_mult_shelf, &
        p_line_grounding_line, p_line_calving_front, p_line_ice_front, p_line_coastline, &
        xmin, xmax, ymin, ymax, lambda_M, phi_M, beta_stereo, mesh)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_mesh_from_reduced_geometry

  SUBROUTINE create_mesh_from_reduced_geometry_singlecore( region_name, name, poly_mult_sheet, poly_mult_shelf, &
    p_line_grounding_line, p_line_calving_front, p_line_ice_front, p_line_coastline, &
    xmin, xmax, ymin, ymax, lambda_M, phi_M, beta_stereo, mesh)
    ! Create mesh from the ice geometry lines
    !
    ! Single-core version; all processes generate the same mesh independently

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=3),           INTENT(IN)        :: region_name
    CHARACTER(LEN=256),         INTENT(IN)        :: name
    REAL(dp), DIMENSION(:,:  ), INTENT(IN)        :: poly_mult_sheet
    REAL(dp), DIMENSION(:,:  ), INTENT(IN)        :: poly_mult_shelf
    REAL(dp), DIMENSION(:,:  ), INTENT(IN)        :: p_line_grounding_line
    REAL(dp), DIMENSION(:,:  ), INTENT(IN)        :: p_line_calving_front
    REAL(dp), DIMENSION(:,:  ), INTENT(IN)        :: p_line_ice_front
    REAL(dp), DIMENSION(:,:  ), INTENT(IN)        :: p_line_coastline
    REAL(dp),                   INTENT(IN)        :: xmin, xmax, ymin, ymax
    REAL(dp),                   INTENT(IN)        :: lambda_M, phi_M, beta_stereo
    TYPE(type_mesh),            INTENT(OUT)       :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'create_mesh_from_reduced_geometry_singlecore'
    REAL(dp)                                      :: res_max_uniform_applied
    INTEGER                                       :: n1,nn,n2
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE       :: poly
    INTEGER                                       :: i

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Single-core mesh generation: let the master do this,
    ! and then broadcast its result to all the other processes.
    IF (par%master) THEN

      ! Allocate mesh memory
      CALL allocate_mesh_primary( mesh, name, 1000, 2000, C%nC_mem)

      ! Initialise the dummy mesh
      CALL initialise_dummy_mesh( mesh, xmin, xmax, ymin, ymax)

    ! == Refine to a uniform resolution; iteratively reduce this,
    ! == and smooth the mesh in between to get a nice, high-quality mesh
    ! ==================================================================

      res_max_uniform_applied = MAX( xmax-xmin, ymax-ymin)

      DO WHILE (.TRUE.)

        ! Reduce the applied uniform resolution
        res_max_uniform_applied = MAX( res_max_uniform_applied / 2._dp, C%maximum_resolution_uniform)

        ! Refine the mesh to the applied uniform resolution
        CALL refine_mesh_uniform( mesh, res_max_uniform_applied, C%alpha_min)

        ! Stop refining once we've reached the desired resolution
        IF (res_max_uniform_applied == C%maximum_resolution_uniform) EXIT

      END DO ! DO WHILE (.TRUE.)

    ! == Refine along the ice geometry lines (grounding line, calving front, etc.)
    ! ============================================================================

      ! Refine the mesh along the ice geometry lines
      CALL refine_mesh_line( mesh, p_line_grounding_line, C%maximum_resolution_grounding_line, C%grounding_line_width, C%alpha_min)
      CALL refine_mesh_line( mesh, p_line_calving_front , C%maximum_resolution_calving_front , C%calving_front_width , C%alpha_min)
      CALL refine_mesh_line( mesh, p_line_ice_front     , C%maximum_resolution_ice_front     , C%ice_front_width     , C%alpha_min)
      CALL refine_mesh_line( mesh, p_line_coastline     , C%maximum_resolution_coastline     , C%coastline_width     , C%alpha_min)

    ! == Refine along the ice geometry areas (sheet, shelf, etc.)
    ! ===========================================================

        ! Ice sheet
        ! =========

        n1 = 1
        n2 = 0

        DO WHILE (n2 < SIZE( poly_mult_sheet,1))

          ! Copy a single polygon from poly_mult
          nn = NINT( poly_mult_sheet( n1,1))
          n2 = n1 + nn
          ALLOCATE( poly( nn,2))
          poly = poly_mult_sheet( n1+1:n2,:)
          n1 = n2+1

          ! Refine mesh over this single polygon
          CALL refine_mesh_polygon( mesh, poly, C%maximum_resolution_grounded_ice, C%alpha_min)

          ! Clean up after yourself
          DEALLOCATE( poly)

        END DO ! DO WHILE (n2 < SIZE( poly_mult_sheet,1))

        ! Ice shelf
        ! =========

        n1 = 1
        n2 = 0

        DO WHILE (n2 < SIZE( poly_mult_shelf,1))

          ! Copy a single polygon from poly_mult
          nn = NINT( poly_mult_shelf( n1,1))
          n2 = n1 + nn
          ALLOCATE( poly( nn,2))
          poly = poly_mult_shelf( n1+1:n2,:)
          n1 = n2+1

          ! Refine mesh over this single polygon
          CALL refine_mesh_polygon( mesh, poly, C%maximum_resolution_floating_ice, C%alpha_min)

          ! Clean up after yourself
          DEALLOCATE( poly)

        END DO ! DO WHILE (n2 < SIZE( poly_mult_sheet,1))

    ! == Refine in regions of interest
    ! ================================

      CALL refine_mesh_in_regions_of_interest( region_name, poly_mult_sheet, poly_mult_shelf, &
        p_line_grounding_line, p_line_calving_front, p_line_ice_front, p_line_coastline, mesh)

    ! == Smooth the mesh
    ! ==================

      DO i = 1, C%nit_Lloyds_algorithm
        CALL Lloyds_algorithm_single_iteration( mesh, C%alpha_min)
      END DO

    ! == Enforce contiguous process domains
    ! =====================================

      CALL enforce_contiguous_process_domains( mesh)

    END IF ! IF (par%master) THEN

    ! Broadcast the Master's mesh
    CALL broadcast_mesh( mesh)

  ! == Calculate secondary geometry data
  ! ====================================

    ! Calculate all secondary geometry data
    CALL calc_all_secondary_mesh_data( mesh, lambda_M, phi_M, beta_stereo)

    ! Calculate all matrix operators
    CALL calc_all_matrix_operators_mesh( mesh)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_mesh_from_reduced_geometry_singlecore

  SUBROUTINE create_mesh_from_reduced_geometry_parallelised( region_name, name, poly_mult_sheet, poly_mult_shelf, &
    p_line_grounding_line, p_line_calving_front, p_line_ice_front, p_line_coastline, &
    xmin, xmax, ymin, ymax, lambda_M, phi_M, beta_stereo, mesh)
    ! Create mesh from the ice geometry lines
    !
    ! Parallelised version

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=3),           INTENT(IN)        :: region_name
    CHARACTER(LEN=256),         INTENT(IN)        :: name
    REAL(dp), DIMENSION(:,:  ), INTENT(IN)        :: poly_mult_sheet
    REAL(dp), DIMENSION(:,:  ), INTENT(IN)        :: poly_mult_shelf
    REAL(dp), DIMENSION(:,:  ), INTENT(IN)        :: p_line_grounding_line
    REAL(dp), DIMENSION(:,:  ), INTENT(IN)        :: p_line_calving_front
    REAL(dp), DIMENSION(:,:  ), INTENT(IN)        :: p_line_ice_front
    REAL(dp), DIMENSION(:,:  ), INTENT(IN)        :: p_line_coastline
    REAL(dp),                   INTENT(IN)        :: xmin, xmax, ymin, ymax
    REAL(dp),                   INTENT(IN)        :: lambda_M, phi_M, beta_stereo
    TYPE(type_mesh),            INTENT(OUT)       :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'create_mesh_from_reduced_geometry_parallelised'
    REAL(dp) :: dummy1
    CHARACTER :: dummy2

    ! Add routine to path
    CALL init_routine( routine_name)

    ! To prevent compiler warnings
    dummy2 = region_name( 1:1)
    dummy2 = name( 1:1)
    dummy1 = poly_mult_sheet( 1,1)
    dummy1 = poly_mult_shelf( 1,1)
    dummy1 = p_line_grounding_line( 1,1)
    dummy1 = p_line_calving_front( 1,1)
    dummy1 = p_line_ice_front( 1,1)
    dummy1 = p_line_coastline( 1,1)
    dummy1 = xmin
    dummy1 = xmax
    dummy1 = ymin
    dummy1 = ymax
    dummy1 = lambda_M
    dummy1 = phi_M
    dummy1 = beta_stereo

    ! DENK DROM
    CALL crash('whoopsiedaisy!')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_mesh_from_reduced_geometry_parallelised

! == Region of interest mesh refinement
! =====================================

  SUBROUTINE refine_mesh_in_regions_of_interest( region_name, poly_mult_sheet, poly_mult_shelf, &
    p_line_grounding_line, p_line_calving_front, p_line_ice_front, p_line_coastline, mesh)
    ! Refine the mesh in the specified regions of interest

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=3),           INTENT(IN)        :: region_name
    REAL(dp), DIMENSION(:,:  ), INTENT(IN)        :: poly_mult_sheet
    REAL(dp), DIMENSION(:,:  ), INTENT(IN)        :: poly_mult_shelf
    REAL(dp), DIMENSION(:,:  ), INTENT(IN)        :: p_line_grounding_line
    REAL(dp), DIMENSION(:,:  ), INTENT(IN)        :: p_line_calving_front
    REAL(dp), DIMENSION(:,:  ), INTENT(IN)        :: p_line_ice_front
    REAL(dp), DIMENSION(:,:  ), INTENT(IN)        :: p_line_coastline
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'refine_mesh_in_regions_of_interest'
    INTEGER                                       :: i
    CHARACTER(LEN=256)                            :: all_names_ROI, name_ROI
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE       :: poly_ROI
    INTEGER                                       :: n1,n2,nn
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE       :: poly

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If no regions of interest are specified, do nothing
    IF (C%choice_regions_of_interest == '') THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    all_names_ROI = C%choice_regions_of_interest

    DO WHILE (.TRUE.)

      ! == Parse list of input ROIs
      ! ===========================

      ! Get the first region of interest from the list
      i = INDEX( all_names_ROI, '||')
      IF (i == 0) THEN
        ! There is only one left in the list
        name_ROI = TRIM( all_names_ROI)
        all_names_ROI = ''
      ELSE
        ! Get the first first one from the list and remove it
        name_ROI = all_names_ROI( 1:i-1)
        all_names_ROI = all_names_ROI( i+2:LEN_TRIM( all_names_ROI))
      END IF

      ! == Check validity of requested ROIs
      ! ===================================

      ! Check if current region is indeed defined in the model
      SELECT CASE (name_ROI)
        CASE ('')
         ! No region requested: don't need to do anything
         EXIT
        CASE ('PineIsland','Thwaites','Amery','RiiserLarsen','SipleCoast', 'LarsenC','TransMounts', & ! Antarctica
              'Narsarsuaq', &                                                                         ! Greenland
              'Patagonia', &                                                                          ! Patagonia
              'Tijn_test_ISMIP_HOM_A')                                                                ! Idealised
          ! List of known regions of interest: these pass the test
        CASE DEFAULT
          ! Region not found
          CALL crash('unknown region of interest "' // TRIM( name_ROI) // '"!')
      END SELECT

      ! == Calculate ROIs
      ! =================

      ! Calculate the polygon describing the specified region of interest
      SELECT CASE (region_name)
        CASE ('NAM')
          ! North america

          SELECT CASE (name_ROI)
            CASE DEFAULT
              ! Requested area not in this model domain; skip
              CYCLE
          END SELECT

        CASE ('EAS')
          ! Eurasia

          SELECT CASE (name_ROI)
            CASE DEFAULT
              ! Requested area not in this model domain; skip
              CYCLE
          END SELECT

        CASE ('GRL')
          ! Greenland

          SELECT CASE (name_ROI)
            CASE ('Narsarsuaq')
              CALL calc_polygon_Narsarsuaq( poly_ROI)
            CASE DEFAULT
              ! Requested area not in this model domain; skip
              CYCLE
          END SELECT

        CASE ('ANT')

          SELECT CASE (name_ROI)
            CASE ('PineIsland')
              CALL calc_polygon_Pine_Island_Glacier( poly_ROI)
            CASE ('Thwaites')
              CALL calc_polygon_Thwaites_Glacier( poly_ROI)
            CASE ('Amery')
              CALL calc_polygon_Amery_ice_shelf( poly_ROI)
            CASE ('RiiserLarsen')
              CALL calc_polygon_Riiser_Larsen_ice_shelf( poly_ROI)
            CASE ('SipleCoast')
              CALL calc_polygon_Siple_Coast( poly_ROI)
            CASE ('LarsenC')
              CALL calc_polygon_Larsen_ice_shelf( poly_ROI)
            CASE ('TransMounts')
              CALL calc_polygon_Transantarctic_Mountains( poly_ROI)
            CASE ('Patagonia')
              CALL calc_polygon_Patagonia( poly_ROI)
            CASE ('Tijn_test_ISMIP_HOM_A')
              CALL calc_polygon_Tijn_test_ISMIP_HOM_A( poly_ROI)
            CASE DEFAULT
              ! Requested area not in this model domain; skip
              CYCLE
          END SELECT

        CASE DEFAULT
          CALL crash('unknown region name "' // region_name // '"!')
      END SELECT

      ! Refine the mesh in the specified region of interest
      ! ===================================================

      ! Uniform
      CALL refine_mesh_polygon( mesh, poly_ROI, C%ROI_maximum_resolution_uniform, C%alpha_min)

      ! Polygons: ice sheet, ice shelf

      ! Ice sheet
      ! =========

      n1 = 1
      n2 = 0

      DO WHILE (n2 < SIZE( poly_mult_sheet,1))

        ! Copy a single polygon from poly_mult
        nn = NINT( poly_mult_sheet( n1,1))
        n2 = n1 + nn
        ALLOCATE( poly( nn,2))
        poly = poly_mult_sheet( n1+1:n2,:)
        n1 = n2+1

        ! Refine mesh over this single polygon
        CALL refine_mesh_polygon_ROI( mesh, poly, C%ROI_maximum_resolution_grounded_ice, C%alpha_min, poly_ROI)

        ! Clean up after yourself
        DEALLOCATE( poly)

      END DO ! DO WHILE (n2 < SIZE( poly_mult_sheet,1))

      ! Ice shelf
      ! =========

      n1 = 1
      n2 = 0

      DO WHILE (n2 < SIZE( poly_mult_shelf,1))

        ! Copy a single polygon from poly_mult
        nn = NINT( poly_mult_shelf( n1,1))
        n2 = n1 + nn
        ALLOCATE( poly( nn,2))
        poly = poly_mult_shelf( n1+1:n2,:)
        n1 = n2+1

        ! Refine mesh over this single polygon
        CALL refine_mesh_polygon_ROI( mesh, poly, C%ROI_maximum_resolution_floating_ice, C%alpha_min, poly_ROI)

        ! Clean up after yourself
        DEALLOCATE( poly)

      END DO ! DO WHILE (n2 < SIZE( poly_mult_sheet,1))

      ! Lines: grounding line, calving front, ice front, coastline
      CALL refine_mesh_line_ROI( mesh, p_line_grounding_line, C%ROI_maximum_resolution_grounding_line, C%ROI_grounding_line_width, C%alpha_min, poly_ROI)
      CALL refine_mesh_line_ROI( mesh, p_line_calving_front , C%ROI_maximum_resolution_calving_front , C%ROI_calving_front_width , C%alpha_min, poly_ROI)
      CALL refine_mesh_line_ROI( mesh, p_line_ice_front     , C%ROI_maximum_resolution_ice_front     , C%ROI_ice_front_width     , C%alpha_min, poly_ROI)
      CALL refine_mesh_line_ROI( mesh, p_line_coastline     , C%ROI_maximum_resolution_coastline     , C%ROI_coastline_width     , C%alpha_min, poly_ROI)

      ! Clean up after yourself
      DEALLOCATE( poly_ROI)

      ! If no names are left, we are finished
      IF (all_names_ROI == '') EXIT

    END DO ! DO WHILE (.TRUE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE refine_mesh_in_regions_of_interest

! == Some useful tools
! ====================

  ! Mesh creation success message
  SUBROUTINE write_mesh_success( mesh)
    ! Write the mesh creation success message to the terminal

    USE control_resources_and_error_messaging, ONLY: insert_val_into_string_int, insert_val_into_string_dp

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),            INTENT(IN)        :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'write_mesh_success'
    CHARACTER(LEN=256)                            :: str

    ! Add routine to path
    CALL init_routine( routine_name)

    str = '     Set up ' // colour_string( TRIM( mesh%name),'light blue') // ' with {int_01} vertices and {int_02} triangles' // &
      ', with a resolution of {dp_01} to {dp_02} m'
    CALL insert_val_into_string_int( str, '{int_01}', mesh%nV)
    CALL insert_val_into_string_int( str, '{int_02}', mesh%nTri)
    CALL insert_val_into_string_dp(  str, '{dp_01}', MINVAL( mesh%R))
    CALL insert_val_into_string_dp(  str, '{dp_02}', MAXVAL( mesh%R))

    IF (par%master) WRITE(0,'(A)') TRIM( str)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_mesh_success

  ! Initialise the five-vertex dummy mesh
  SUBROUTINE initialise_dummy_mesh( mesh, xmin, xmax, ymin, ymax)
    ! Initialises a 5-vertex, 4-triangle "dummy"  mesh:
    !
    !   v4 - - - - - - - - v3   V          nC     C             niTri   iTri          edge_index
    !   | \              / |    -1 -1      3      2  5  4        2      1  4            6
    !   |  \            /  |     1 -1      3      3  5  1        2      2  1            4
    !   |   \    t3    /   |     1  1      3      4  5  2        2      3  2            2
    !   |    \        /    |    -1  1      3      1  5  3        2      4  3            8
    !   |     \      /     |     0  0      4      1  2  3  4     4      1  2  3  4      0
    !   |      \    /      |
    !   |       \  /       |    Tri           TriC
    !   |  t4    v5    t2  |    1  2  5      2  4  0
    !   |       /  \       |    2  3  5      3  1  0
    !   |      /    \      |    3  4  5      4  2  0
    !   |     /      \     |    4  1  5      1  3  0
    !   |    /        \    |
    !   |   /    t1    \   |
    !   |  /            \  |
    !   | /              \ |
    !   v1 - - - - - - - - v2
    !
    ! NOTE: memory must already be allocated for the mesh before calling this routine

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    REAL(dp),                   INTENT(IN)        :: xmin, xmax, ymin, ymax

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'initialise_dummy_mesh'
    REAL(dp), PARAMETER                           :: tol = 1E-9_dp

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Meta properties
    mesh%xmin         = xmin    ! Boundaries of the square domain.
    mesh%xmax         = xmax
    mesh%ymin         = ymin
    mesh%ymax         = ymax

    ! Horizontal distance of tolerance. Used for several small routines - points that lie within
    ! this distance of each other (vertices, triangle circumcenters, etc.) are treated as identical.
    mesh%tol_dist     = ((mesh%xmax - mesh%xmin) + (mesh%ymax - mesh%ymin)) * tol / 2._dp

    mesh%nV           = 5

    mesh%V            = 0._dp
    mesh%V( 1,:)      = [xmin, ymin]
    mesh%V( 2,:)      = [xmax, ymin]
    mesh%V( 3,:)      = [xmax, ymax]
    mesh%V( 4,:)      = [xmin, ymax]
    ! Make sure the central vertex is slightly off-centre, to prevent trivial Delaunay criteria
    ! in the early stages of mesh refinement (i.e. 4 or more vertices being cocircular), which
    ! can lead to different meshes depending on processor/compiler/phase of the moon.
    mesh%V( 5,:)      = [(xmin+xmax)/2._dp + (xmax-xmin)*pi*1e-3_dp, (ymin+ymax)/2._dp + (ymax-ymin)*pi*2.1e-3_dp]

    mesh%VBI          = 0
    mesh%VBI(1:5)     = [6, 4, 2, 8, 0]

    mesh%nC           = 0
    mesh%nC( 1:5)     = [3, 3, 3, 3, 4]

    mesh%C            = 0
    mesh%C( 1,1:4)    = [2, 5, 4, 0]
    mesh%C( 2,1:4)    = [3, 5, 1, 0]
    mesh%C( 3,1:4)    = [4, 5, 2, 0]
    mesh%C( 4,1:4)    = [1, 5, 3, 0]
    mesh%C( 5,1:4)    = [1, 2, 3, 4]

    mesh%niTri        = 0
    mesh%niTri( 1:5)  = [2, 2, 2, 2, 4]

    mesh%iTri         = 0
    mesh%iTri( 1,1:4) = [1, 4, 0, 0]
    mesh%iTri( 2,1:4) = [2, 1, 0, 0]
    mesh%iTri( 3,1:4) = [3, 2, 0, 0]
    mesh%iTri( 4,1:4) = [4, 3, 0, 0]
    mesh%iTri( 5,1:4) = [1, 2, 3, 4]

    mesh%nTri         = 4

    mesh%Tri          = 0
    mesh%Tri( 1,:)    = [1, 2, 5]
    mesh%Tri( 2,:)    = [2, 3, 5]
    mesh%Tri( 3,:)    = [3, 4, 5]
    mesh%Tri( 4,:)    = [4, 1, 5]

    mesh%TriC         = 0
    mesh%TriC( 1,:)   = [2, 4, 0]
    mesh%TriC( 2,:)   = [3, 1, 0]
    mesh%TriC( 3,:)   = [4, 2, 0]
    mesh%TriC( 4,:)   = [1, 3, 0]

    mesh%TriCC = 0._dp
    CALL update_triangle_circumcenter( mesh, 1)
    CALL update_triangle_circumcenter( mesh, 2)
    CALL update_triangle_circumcenter( mesh, 3)
    CALL update_triangle_circumcenter( mesh, 4)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_dummy_mesh

END MODULE mesh_creation
