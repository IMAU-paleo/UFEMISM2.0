module mesh_creation_main

  ! Routines used to create a mesh.

  use mpi
  use precisions, only: dp
  use mpi_basic, only: par, cerr, ierr, recv_status, sync
  use control_resources_and_error_messaging, only: warning, crash, happy, init_routine, finalise_routine, colour_string
  use model_configuration, only: C
  use parameters
  use mesh_types, only: type_mesh
  use mesh_memory, only: allocate_mesh_primary
  use mesh_utilities, only: update_triangle_circumcenter, calc_mesh_contour_as_line, calc_mesh_mask_as_polygons
  use mesh_refinement_basic, only: refine_mesh_uniform, refine_mesh_point, refine_mesh_line, refine_mesh_polygon
  use mesh_refinement_basic_ROI, only: refine_mesh_line_ROI, refine_mesh_polygon_ROI
  use mesh_contiguous_domains, only: enforce_contiguous_process_domains
  use mesh_ROI_polygons
  use mesh_parallel_creation, only: broadcast_mesh
  use mesh_secondary, only: calc_all_secondary_mesh_data
  use mesh_operators, only: calc_all_matrix_operators_mesh
  use grid_basic, only: type_grid, gather_gridded_data_to_master_dp_2D, calc_grid_contour_as_line, &
    calc_grid_mask_as_polygons
  use math_utilities, only: thickness_above_floatation
  use mpi_distributed_memory, only: gather_to_master_dp_1D
  use mesh_refinement_fun, only: refine_CalvMIP_shelf_donut
  use mesh_Lloyds_algorithm, only: Lloyds_algorithm_single_iteration
  use mesh_dummy_meshes, only: initialise_dummy_mesh_5

  implicit none

  private

  public :: create_mesh_from_gridded_geometry, create_mesh_from_meshed_geometry, write_mesh_success

contains

! == The two main mesh creation routines that should be called
! ============================================================

  ! Create a mesh from ice geometry on a grid
  subroutine create_mesh_from_gridded_geometry( region_name, name, &
    grid, Hi, Hb, Hs, SL, xmin, xmax, ymin, ymax, lambda_M, phi_M, beta_stereo, mesh)
    ! Create a mesh from ice geometry on a grid

    ! In/output variables:
    character(len=3),           intent(in   ) :: region_name
    character(len=256),         intent(in   ) :: name
    type(type_grid),            intent(in   ) :: grid
    real(dp), dimension(:    ), intent(in   ) :: Hi, Hb, Hs, SL
    real(dp),                   intent(in   ) :: xmin, xmax, ymin, ymax
    real(dp),                   intent(in   ) :: lambda_M, phi_M, beta_stereo
    type(type_mesh),            intent(  out) :: mesh

    ! Local variables:
    character(len=256), parameter           :: routine_name = 'create_mesh_from_gridded_geometry'
    real(dp), dimension(:,:  ), allocatable :: poly_mult_sheet
    real(dp), dimension(:,:  ), allocatable :: poly_mult_shelf
    real(dp), dimension(:,:  ), allocatable :: p_line_grounding_line
    real(dp), dimension(:,:  ), allocatable :: p_line_calving_front
    real(dp), dimension(:,:  ), allocatable :: p_line_ice_front
    real(dp), dimension(:,:  ), allocatable :: p_line_coastline

    ! Add routine to path
    call init_routine( routine_name)

    ! Reduce the gridded ice geometry to lines and polygons
    call reduce_gridded_ice_geometry( grid, Hi, Hb, Hs, SL, poly_mult_sheet, poly_mult_shelf, &
      p_line_grounding_line, p_line_calving_front, p_line_ice_front, p_line_coastline)

    ! Create a mesh from the reduced ice geometry
    call create_mesh_from_reduced_geometry( region_name, name, poly_mult_sheet, poly_mult_shelf, &
      p_line_grounding_line, p_line_calving_front, p_line_ice_front, p_line_coastline, &
      xmin, xmax, ymin, ymax, lambda_M, phi_M, beta_stereo, mesh)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_mesh_from_gridded_geometry

  ! Create a mesh from ice geometry on a mesh
  subroutine create_mesh_from_meshed_geometry( region_name, name, mesh_src, Hi, Hb, Hs, SL, xmin, xmax, ymin, ymax, lambda_M, phi_M, beta_stereo, mesh)
    ! Create a mesh from ice geometry on a mesh

    ! In/output variables:
    character(len=3),           intent(in)        :: region_name
    character(len=256),         intent(in)        :: name
    type(type_mesh),            intent(in)        :: mesh_src
    real(dp), dimension(:    ), intent(in)        :: Hi, Hb, Hs, SL
    real(dp),                   intent(in)        :: xmin, xmax, ymin, ymax
    real(dp),                   intent(in)        :: lambda_M, phi_M, beta_stereo
    type(type_mesh),            intent(out)       :: mesh

    ! Local variables:
    character(len=256), parameter                 :: routine_name = 'create_mesh_from_meshed_geometry'
    real(dp), dimension(:,:  ), allocatable       :: poly_mult_sheet
    real(dp), dimension(:,:  ), allocatable       :: poly_mult_shelf
    real(dp), dimension(:,:  ), allocatable       :: p_line_grounding_line
    real(dp), dimension(:,:  ), allocatable       :: p_line_calving_front
    real(dp), dimension(:,:  ), allocatable       :: p_line_ice_front
    real(dp), dimension(:,:  ), allocatable       :: p_line_coastline

    ! Add routine to path
    call init_routine( routine_name)

    ! Reduce the meshed ice geometry to lines and polygons
    call reduce_meshed_ice_geometry( mesh_src, Hi, Hb, Hs, SL, poly_mult_sheet, poly_mult_shelf, &
      p_line_grounding_line, p_line_calving_front, p_line_ice_front, p_line_coastline)

    ! Create a mesh from the reduced ice geometry
    call create_mesh_from_reduced_geometry( region_name, name, poly_mult_sheet, poly_mult_shelf, &
      p_line_grounding_line, p_line_calving_front, p_line_ice_front, p_line_coastline, &
      xmin, xmax, ymin, ymax, lambda_M, phi_M, beta_stereo, mesh)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_mesh_from_meshed_geometry

! == Reduce a gridded/meshed ice geometry to grid-less polygons and lines
! =======================================================================

  subroutine reduce_gridded_ice_geometry( grid, Hi, Hb, Hs, SL, poly_mult_sheet, poly_mult_shelf, &
    p_line_grounding_line, p_line_calving_front, p_line_ice_front, p_line_coastline)

    ! "Reduce" the gridded ice-sheet geometry to a set of polygons
    ! (describing regions covered by grounded/floating ice) and lines
    ! (grounding line, calving front, etc.)
    !
    ! NOTE: gridded geometry fields Hi, Hb, Hs, SL are provided in distributed vector form
    !
    ! NOTE: polygons and lines are returned in full to all processes

    ! In/output variables:
    type(type_grid),                         intent(in   ) :: grid
    real(dp), dimension(:    ),              intent(in   ) :: Hi, Hb, Hs, SL
    real(dp), dimension(:,:  ), allocatable, intent(  out) :: poly_mult_sheet
    real(dp), dimension(:,:  ), allocatable, intent(  out) :: poly_mult_shelf
    real(dp), dimension(:,:  ), allocatable, intent(  out) :: p_line_grounding_line
    real(dp), dimension(:,:  ), allocatable, intent(  out) :: p_line_calving_front
    real(dp), dimension(:,:  ), allocatable, intent(  out) :: p_line_ice_front
    real(dp), dimension(:,:  ), allocatable, intent(  out) :: p_line_coastline

    ! Local variables:
    character(len=256), parameter           :: routine_name = 'reduce_gridded_ice_geometry'
    real(dp), dimension(:,:  ), allocatable :: Hi_grid
    real(dp), dimension(:,:  ), allocatable :: Hb_grid
    real(dp), dimension(:,:  ), allocatable :: Hs_grid
    real(dp), dimension(:,:  ), allocatable :: SL_grid
    real(dp), dimension(:,:  ), allocatable :: TAF_grid
    real(dp), dimension(:,:  ), allocatable :: Hb_minus_SL_grid
    integer                                 :: i,j,ii,jj
    logical,  dimension(:,:  ), allocatable :: mask_sheet
    logical,  dimension(:,:  ), allocatable :: mask_shelf
    logical,  dimension(:,:  ), allocatable :: mask_calc_grounding_line
    logical,  dimension(:,:  ), allocatable :: mask_calc_calving_front
    logical,  dimension(:,:  ), allocatable :: mask_calc_ice_front
    logical,  dimension(:,:  ), allocatable :: mask_calc_coastline
    integer                                 :: n_poly_mult_sheet
    integer                                 :: n_poly_mult_shelf
    integer                                 :: n_line_grounding_line
    integer                                 :: n_line_calving_front
    integer                                 :: n_line_ice_front
    integer                                 :: n_line_coastline

    ! Add routine to path
    call init_routine( routine_name)

    ! allocate memory for gathered ice geometry data
    if (par%master) then
      allocate( Hi_grid( grid%nx, grid%ny))
      allocate( Hb_grid( grid%nx, grid%ny))
      allocate( Hs_grid( grid%nx, grid%ny))
      allocate( SL_grid( grid%nx, grid%ny))
    else
      allocate( Hi_grid( 0, 0))
      allocate( Hb_grid( 0, 0))
      allocate( Hs_grid( 0, 0))
      allocate( SL_grid( 0, 0))
    end if

    ! Gather ice geometry data in grid form to the master
    call gather_gridded_data_to_master_dp_2D( grid, Hi, Hi_grid)
    call gather_gridded_data_to_master_dp_2D( grid, Hb, Hb_grid)
    call gather_gridded_data_to_master_dp_2D( grid, Hs, Hs_grid)
    call gather_gridded_data_to_master_dp_2D( grid, SL, SL_grid)

    ! Let the master calculate the polygons and lines

    if (par%master) then

      ! allocate memory for thickness above floatation and Hb-SL
      allocate( TAF_grid(         grid%nx, grid%ny))
      allocate( Hb_minus_SL_grid( grid%nx, grid%ny))

      ! Calculate thickness above floatation and Hb-SL
      do i = 1, grid%nx
      do j = 1, grid%ny
        TAF_grid(         i,j) = thickness_above_floatation( Hi_grid( i,j), Hb_grid( i,j), SL_grid( i,j))
        Hb_minus_SL_grid( i,j) = Hb_grid( i,j) - SL_grid( i,j)
      end do
      end do

      ! Fill in masks for floating/grounded ice
      allocate( mask_sheet( grid%nx, grid%ny), source = .false)
      allocate( mask_shelf( grid%nx, grid%ny), source = .false)

      do i = 1, grid%nx
      do j = 1, grid%ny
        if (Hi_grid( i,j) > 0.1_dp) then
          if (TAF_grid( i,j) > 0._dp) then
            mask_sheet( i,j) = .true.
          else
            mask_shelf( i,j) = .true.
          end if
        end if
      end do
      end do

      ! Fill in masks for where to calculate lines
      allocate( mask_calc_grounding_line( grid%nx, grid%ny), source = .false)
      allocate( mask_calc_calving_front(  grid%nx, grid%ny), source = .false)
      allocate( mask_calc_ice_front(      grid%nx, grid%ny), source = .false)
      allocate( mask_calc_coastline(      grid%nx, grid%ny), source = .false)

      do i = 1, grid%nx
      do j = 1, grid%ny

        ! Grounding line should only be calculated where there's ice
        if (Hi_grid( i,j) > 0.1_dp) then
          mask_calc_grounding_line( i,j) = .true.
        end if

        ! Calving front should only be calculated for ice (both floating and grounded) next to ocean
        if (Hi_grid( i,j) > 0.1_dp) then
          ! This grid cell has ice
          do ii = max( 1, i-1), min( grid%nx, i+1)
          do jj = max( 1, j-1), min( grid%ny, j+1)
            if (Hi_grid( ii,jj) < 0.1_dp .AND. Hb_grid( ii,jj) < SL_grid( ii,jj)) then
              ! This neighbour is open ocean
              mask_calc_calving_front( i,j) = .true.
            end if
          end do
          end do
        end if

        ! Ice front is simply ice next to non-ice
        if (Hi_grid( i,j) > 0.1_dp) then
          ! This grid cell has ice
          do ii = max( 1, i-1), min( grid%nx, i+1)
          do jj = max( 1, j-1), min( grid%ny, j+1)
            if (Hi_grid( ii,jj) < 0.1_dp) then
              ! This neighbour is ice-free
              mask_calc_ice_front( i,j) = .true.
            end if
          end do
          end do
        end if

        ! Coastline is ice-free land next to open ocean
        if (Hi_grid( i,j) < 0.1_dp .AND. Hb_grid( i,j) > SL_grid( i,j)) then
          ! This grid cell is ice-free land
          do ii = max( 1, i-1), min( grid%nx, i+1)
          do jj = max( 1, j-1), min( grid%ny, j+1)
            if (Hi_grid( ii,jj) < 0.1_dp .AND. Hb_grid( ii,jj) < SL_grid( ii,jj)) then
              ! This neighbour is open ocean
              mask_calc_coastline( i,j) = .true.
            end if
          end do
          end do
        end if

      end do
      end do

      ! Calculate polygons enveloping sheet and shelf
      call calc_grid_mask_as_polygons( grid, mask_sheet, poly_mult_sheet)
      call calc_grid_mask_as_polygons( grid, mask_shelf, poly_mult_shelf)

      ! Get polygon sizes
      n_poly_mult_sheet = size( poly_mult_sheet,1)
      n_poly_mult_shelf = size( poly_mult_shelf,1)

      ! Calculate lines in line segment format
      call calc_grid_contour_as_line( grid, TAF_grid        , 0.0_dp, p_line_grounding_line, mask_calc_grounding_line)
      call calc_grid_contour_as_line( grid, Hi_grid         , 0.1_dp, p_line_calving_front , mask_calc_calving_front )
      call calc_grid_contour_as_line( grid, Hi_grid         , 0.1_dp, p_line_ice_front     , mask_calc_ice_front     )
      call calc_grid_contour_as_line( grid, Hb_minus_SL_grid, 0.0_dp, p_line_coastline     , mask_calc_coastline     )

      ! Get line sizes
      n_line_grounding_line = size( p_line_grounding_line,1)
      n_line_calving_front  = size( p_line_calving_front ,1)
      n_line_ice_front      = size( p_line_ice_front     ,1)
      n_line_coastline      = size( p_line_coastline     ,1)

    end if ! if (par%master) then
    call sync

    ! Broadcast polygon sizes to all processes
    call MPI_BCAST( n_poly_mult_sheet, 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST( n_poly_mult_shelf, 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)

    ! allocate memory for the polygons on the other processes
    if (.not. par%master) then
      allocate( poly_mult_sheet( n_poly_mult_sheet,2))
      allocate( poly_mult_shelf( n_poly_mult_shelf,2))
    end if

    ! Broadcast polygons to all processes
    call MPI_BCAST( poly_mult_sheet, n_poly_mult_sheet*2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST( poly_mult_shelf, n_poly_mult_shelf*2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! Broadcast line sizes to all processes
    call MPI_BCAST( n_line_grounding_line, 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST( n_line_calving_front , 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST( n_line_ice_front     , 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST( n_line_coastline     , 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)

    ! allocate memory for the lines on the other processes
    if (.not. par%master) then
      allocate( p_line_grounding_line( n_line_grounding_line,4))
      allocate( p_line_calving_front(  n_line_calving_front ,4))
      allocate( p_line_ice_front(      n_line_ice_front     ,4))
      allocate( p_line_coastline(      n_line_coastline     ,4))
    end if

    ! Broadcast lines to all processes
    call MPI_BCAST( p_line_grounding_line, n_line_grounding_line*4, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST( p_line_calving_front , n_line_calving_front *4, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST( p_line_ice_front     , n_line_ice_front     *4, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST( p_line_coastline     , n_line_coastline     *4, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine reduce_gridded_ice_geometry

  subroutine reduce_meshed_ice_geometry( mesh, Hi, Hb, Hs, SL, poly_mult_sheet, poly_mult_shelf, &
    p_line_grounding_line, p_line_calving_front, p_line_ice_front, p_line_coastline)

    ! "Reduce" the meshed ice-sheet geometry to a set of polygons
    ! (describing regions covered by grounded/floating ice) and lines
    ! (grounding line, calving front, etc.)
    !
    ! NOTE: meshed geometry fields Hi, Hb, Hs, SL are provided in distributed vector form
    !
    ! NOTE: polygons and lines are returned in full to all processes

    ! In/output variables:
    type(type_mesh),                         intent(in)    :: mesh
    real(dp), dimension(:    ),              intent(in)    :: Hi, Hb, Hs, SL
    real(dp), dimension(:,:  ), allocatable, intent(out)   :: poly_mult_sheet
    real(dp), dimension(:,:  ), allocatable, intent(out)   :: poly_mult_shelf
    real(dp), dimension(:,:  ), allocatable, intent(out)   :: p_line_grounding_line
    real(dp), dimension(:,:  ), allocatable, intent(out)   :: p_line_calving_front
    real(dp), dimension(:,:  ), allocatable, intent(out)   :: p_line_ice_front
    real(dp), dimension(:,:  ), allocatable, intent(out)   :: p_line_coastline

    ! Local variables:
    character(len=256), parameter                 :: routine_name = 'reduce_meshed_ice_geometry'
    real(dp), dimension(:    ), allocatable       :: Hi_tot
    real(dp), dimension(:    ), allocatable       :: Hb_tot
    real(dp), dimension(:    ), allocatable       :: Hs_tot
    real(dp), dimension(:    ), allocatable       :: SL_tot
    real(dp), dimension(:    ), allocatable       :: TAF_tot
    real(dp), dimension(:    ), allocatable       :: Hb_minus_SL_tot
    integer                                       :: vi
    logical,  dimension(:    ), allocatable       :: mask_sheet
    logical,  dimension(:    ), allocatable       :: mask_shelf
    logical,  dimension(:    ), allocatable       :: mask_calc_grounding_line
    logical,  dimension(:    ), allocatable       :: mask_calc_calving_front
    logical,  dimension(:    ), allocatable       :: mask_calc_coastline
    integer                                       :: n_poly_mult_sheet
    integer                                       :: n_poly_mult_shelf
    integer                                       :: n_line_grounding_line
    integer                                       :: n_line_calving_front
    integer                                       :: n_line_ice_front
    integer                                       :: n_line_coastline

    ! Add routine to path
    call init_routine( routine_name)

    ! allocate memory for gathered ice geometry data
    if (par%master) then
      allocate( Hi_tot( mesh%nV))
      allocate( Hb_tot( mesh%nV))
      allocate( Hs_tot( mesh%nV))
      allocate( SL_tot( mesh%nV))
    end if

    ! Gather ice geometry data in grid form to the master
    call gather_to_master_dp_1D( Hi, Hi_tot)
    call gather_to_master_dp_1D( Hb, Hb_tot)
    call gather_to_master_dp_1D( Hs, Hs_tot)
    call gather_to_master_dp_1D( SL, SL_tot)

    ! Let the master calculate the polygons and lines

    if (par%master) then

      ! allocate memory for thickness above floatation and Hb-SL
      allocate( TAF_tot(         mesh%nV))
      allocate( Hb_minus_SL_tot( mesh%nV))

      ! Calculate thickness above floatation and Hb-SL
      do vi = 1, mesh%nV
        TAF_tot(         vi) = thickness_above_floatation( Hi_tot( vi), Hb_tot( vi), SL_tot( vi))
        Hb_minus_SL_tot( vi) = Hb_tot( vi) - SL_tot( vi)
      end do

      ! Fill in masks for floating/grounded ice
      allocate( mask_sheet( mesh%nV), source = .false)
      allocate( mask_shelf( mesh%nV), source = .false)

      do vi = 1, mesh%nV
        if (Hi_tot( vi) > 0.1_dp) then
          if (TAF_tot( vi) > 0._dp) then
            mask_sheet( vi) = .true.
          else
            mask_shelf( vi) = .true.
          end if
        end if
      end do

      ! Fill in masks for where to calculate lines
      allocate( mask_calc_grounding_line( mesh%nV), source = .false)
      allocate( mask_calc_calving_front ( mesh%nV), source = .false)
      allocate( mask_calc_coastline     ( mesh%nV), source = .false)

      do vi = 1, mesh%nV

        ! Grounding line should only be calculated where there's ice
        if (Hi_tot( vi) > 0.1_dp) then
          mask_calc_grounding_line( vi) = .true.
        end if

        ! Calving front should only be calculated on floating ice
        if (TAF_tot( vi) < 0.1_dp .AND. Hi_tot( vi) > 0.1_dp) then
           mask_calc_calving_front( vi) = .true.
        end if

        ! Coastline should only be calculated where there's no ice
        if (Hi_tot( vi) < 0.1_dp ) then
           mask_calc_coastline( vi) = .true.
        end if

      end do ! do vi = 1, mesh%nV

      ! Calculate polygons enveloping sheet and shelf
      call calc_mesh_mask_as_polygons( mesh, mask_sheet, poly_mult_sheet)
      call calc_mesh_mask_as_polygons( mesh, mask_shelf, poly_mult_shelf)

      ! Get polygon sizes
      n_poly_mult_sheet = size( poly_mult_sheet,1)
      n_poly_mult_shelf = size( poly_mult_shelf,1)

      ! Calculate lines in line segment format
      call calc_mesh_contour_as_line( mesh, TAF_tot        , 0.0_dp, p_line_grounding_line, mask_calc_grounding_line)
      call calc_mesh_contour_as_line( mesh, Hi_tot         , 1.0_dp, p_line_calving_front , mask_calc_calving_front )
      call calc_mesh_contour_as_line( mesh, Hi_tot         , 1.0_dp, p_line_ice_front                               )
      call calc_mesh_contour_as_line( mesh, Hb_minus_SL_tot, 0.0_dp, p_line_coastline     , mask_calc_coastline     )

      ! Get line sizes
      n_line_grounding_line = size( p_line_grounding_line,1)
      n_line_calving_front  = size( p_line_calving_front ,1)
      n_line_ice_front      = size( p_line_ice_front     ,1)
      n_line_coastline      = size( p_line_coastline     ,1)

    end if ! if (par%master) then
    call sync

    ! Broadcast polygon sizes to all processes
    call MPI_BCAST( n_poly_mult_sheet, 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST( n_poly_mult_shelf, 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)

    ! allocate memory for the polygons on the other processes
    if (.not. par%master) then
      allocate( poly_mult_sheet( n_poly_mult_sheet,2))
      allocate( poly_mult_shelf( n_poly_mult_shelf,2))
    end if

    ! Broadcast polygons to all processes
    call MPI_BCAST( poly_mult_sheet, n_poly_mult_sheet*2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST( poly_mult_shelf, n_poly_mult_shelf*2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! Broadcast line sizes to all processes
    call MPI_BCAST( n_line_grounding_line, 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST( n_line_calving_front , 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST( n_line_ice_front     , 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST( n_line_coastline     , 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)

    ! allocate memory for the lines on the other processes
    if (.not. par%master) then
      allocate( p_line_grounding_line( n_line_grounding_line,4))
      allocate( p_line_calving_front(  n_line_calving_front ,4))
      allocate( p_line_ice_front(      n_line_ice_front     ,4))
      allocate( p_line_coastline(      n_line_coastline     ,4))
    end if

    ! Broadcast lines to all processes
    call MPI_BCAST( p_line_grounding_line, n_line_grounding_line*4, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST( p_line_calving_front , n_line_calving_front *4, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST( p_line_ice_front     , n_line_ice_front     *4, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST( p_line_coastline     , n_line_coastline     *4, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine reduce_meshed_ice_geometry

! == Create a mesh from reduced ice geometry
! ==========================================

  subroutine create_mesh_from_reduced_geometry( region_name, name, poly_mult_sheet, poly_mult_shelf, &
    p_line_grounding_line, p_line_calving_front, p_line_ice_front, p_line_coastline, &
    xmin, xmax, ymin, ymax, lambda_M, phi_M, beta_stereo, mesh)
    ! Create mesh from the reduced ice geometry

    ! In/output variables:
    character(len=3),           intent(in)        :: region_name
    character(len=256),         intent(in)        :: name
    real(dp), dimension(:,:  ), intent(in)        :: poly_mult_sheet
    real(dp), dimension(:,:  ), intent(in)        :: poly_mult_shelf
    real(dp), dimension(:,:  ), intent(in)        :: p_line_grounding_line
    real(dp), dimension(:,:  ), intent(in)        :: p_line_calving_front
    real(dp), dimension(:,:  ), intent(in)        :: p_line_ice_front
    real(dp), dimension(:,:  ), intent(in)        :: p_line_coastline
    real(dp),                   intent(in)        :: xmin, xmax, ymin, ymax
    real(dp),                   intent(in)        :: lambda_M, phi_M, beta_stereo
    type(type_mesh),            intent(out)       :: mesh

    ! Local variables:
    character(len=256), parameter                 :: routine_name = 'create_mesh_from_reduced_geometry'

    ! Add routine to path
    call init_routine( routine_name)

    ! Choose single-core or parallelised version
    if (C%do_singlecore_mesh_creation) then
      call create_mesh_from_reduced_geometry_singlecore( region_name, name, poly_mult_sheet, poly_mult_shelf, &
        p_line_grounding_line, p_line_calving_front, p_line_ice_front, p_line_coastline, &
        xmin, xmax, ymin, ymax, lambda_M, phi_M, beta_stereo, mesh)
    else
      call create_mesh_from_reduced_geometry_parallelised( region_name, name, poly_mult_sheet, poly_mult_shelf, &
        p_line_grounding_line, p_line_calving_front, p_line_ice_front, p_line_coastline, &
        xmin, xmax, ymin, ymax, lambda_M, phi_M, beta_stereo, mesh)
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_mesh_from_reduced_geometry

  subroutine create_mesh_from_reduced_geometry_singlecore( region_name, name, poly_mult_sheet, poly_mult_shelf, &
    p_line_grounding_line, p_line_calving_front, p_line_ice_front, p_line_coastline, &
    xmin, xmax, ymin, ymax, lambda_M, phi_M, beta_stereo, mesh)
    ! Create mesh from the ice geometry lines
    !
    ! Single-core version; all processes generate the same mesh independently

    ! In/output variables:
    character(len=3),           intent(in)        :: region_name
    character(len=256),         intent(in)        :: name
    real(dp), dimension(:,:  ), intent(in)        :: poly_mult_sheet
    real(dp), dimension(:,:  ), intent(in)        :: poly_mult_shelf
    real(dp), dimension(:,:  ), intent(in)        :: p_line_grounding_line
    real(dp), dimension(:,:  ), intent(in)        :: p_line_calving_front
    real(dp), dimension(:,:  ), intent(in)        :: p_line_ice_front
    real(dp), dimension(:,:  ), intent(in)        :: p_line_coastline
    real(dp),                   intent(in)        :: xmin, xmax, ymin, ymax
    real(dp),                   intent(in)        :: lambda_M, phi_M, beta_stereo
    type(type_mesh),            intent(out)       :: mesh

    ! Local variables:
    character(len=256), parameter                 :: routine_name = 'create_mesh_from_reduced_geometry_singlecore'
    real(dp)                                      :: res_max_uniform_applied
    integer                                       :: n1,nn,n2
    real(dp), dimension(:,:  ), allocatable       :: poly
    integer                                       :: i

    ! Add routine to path
    call init_routine( routine_name)

    ! Single-core mesh generation: let the master do this,
    ! and then broadcast its result to all the other processes.
    if (par%master) then

      ! allocate mesh memory
      call allocate_mesh_primary( mesh, name, 1000, 2000, C%nC_mem)

      ! Initialise the dummy mesh
      call initialise_dummy_mesh_5( mesh, xmin, xmax, ymin, ymax)

      ! == Refine to a uniform resolution; iteratively reduce this,
      ! == and smooth the mesh in between to get a nice, high-quality mesh
      ! ==================================================================

      res_max_uniform_applied = max( xmax-xmin, ymax-ymin)

      do while (.true.)

        ! Reduce the applied uniform resolution
        res_max_uniform_applied = max( res_max_uniform_applied / 2._dp, C%maximum_resolution_uniform)

        ! Refine the mesh to the applied uniform resolution
        call refine_mesh_uniform( mesh, res_max_uniform_applied, C%alpha_min)

        ! Stop refining once we've reached the desired resolution
        if (res_max_uniform_applied == C%maximum_resolution_uniform) exit

      end do ! do while (.true.)

      ! == Refine along the ice geometry lines (grounding line, calving front, etc.)
      ! ============================================================================

      ! Refine the mesh along the ice geometry lines
      call refine_mesh_line( mesh, p_line_grounding_line, C%maximum_resolution_grounding_line, C%grounding_line_width, C%alpha_min)
      call refine_mesh_line( mesh, p_line_calving_front , C%maximum_resolution_calving_front , C%calving_front_width , C%alpha_min)
      call refine_mesh_line( mesh, p_line_ice_front     , C%maximum_resolution_ice_front     , C%ice_front_width     , C%alpha_min)
      call refine_mesh_line( mesh, p_line_coastline     , C%maximum_resolution_coastline     , C%coastline_width     , C%alpha_min)

      ! == Refine along the ice geometry areas (sheet, shelf, etc.)
      ! ===========================================================

        ! Ice sheet
        ! =========

        n1 = 1
        n2 = 0

        do while (n2 < size( poly_mult_sheet,1))

          ! Copy a single polygon from poly_mult
          nn = nint( poly_mult_sheet( n1,1))
          n2 = n1 + nn
          allocate( poly( nn,2))
          poly = poly_mult_sheet( n1+1:n2,:)
          n1 = n2+1

          ! Refine mesh over this single polygon
          call refine_mesh_polygon( mesh, poly, C%maximum_resolution_grounded_ice, C%alpha_min)

          ! Clean up after yourself
          deallocate( poly)

        end do ! do while (n2 < size( poly_mult_sheet,1))

        ! Ice shelf
        ! =========

        n1 = 1
        n2 = 0

        do while (n2 < size( poly_mult_shelf,1))

          ! Copy a single polygon from poly_mult
          nn = nint( poly_mult_shelf( n1,1))
          n2 = n1 + nn
          allocate( poly( nn,2))
          poly = poly_mult_shelf( n1+1:n2,:)
          n1 = n2+1

          ! Refine mesh over this single polygon. Use the ice sheet polygon set as a
          ! "no-refinement" zone to avoid extreme cases where the ice shelf polygon
          ! encompases the ice sheet one (e.g. in the circular domain of CalvMIP)
          if (C%choice_refgeo_PD_ANT == 'idealised') then
            call refine_mesh_polygon( mesh, poly, C%maximum_resolution_floating_ice, C%alpha_min, poly_mult_sheet)
          else
            call refine_mesh_polygon( mesh, poly, C%maximum_resolution_floating_ice, C%alpha_min)
          end if

          ! Clean up after yourself
          deallocate( poly)

        end do ! do while (n2 < size( poly_mult_sheet,1))

      ! == Refine in regions of interest
      ! ================================

      call refine_mesh_in_regions_of_interest( region_name, poly_mult_sheet, poly_mult_shelf, &
        p_line_grounding_line, p_line_calving_front, p_line_ice_front, p_line_coastline, mesh)

      ! == Special cases
      ! ================

      ! DENK DROM : Not very elegant; remove this later and generalise it
      if (C%do_ANT .AND. C%choice_refgeo_PD_ANT == 'idealised' .AND. C%choice_refgeo_init_idealised == 'calvmip_circular') then
        call refine_CalvMIP_shelf_donut( mesh, C%maximum_resolution_grounding_line*2._dp, 70000._dp)
      elseif (C%do_ANT .AND. C%choice_refgeo_PD_ANT == 'idealised' .AND. C%choice_refgeo_init_idealised == 'calvmip_Thule') then
        call refine_CalvMIP_shelf_donut( mesh, C%maximum_resolution_grounding_line*2._dp, 50000._dp)
      end if

      ! == Smooth the mesh
      ! ==================

      do i = 1, C%nit_Lloyds_algorithm
        call Lloyds_algorithm_single_iteration( mesh, C%alpha_min)
      end do

      ! == Enforce contiguous process domains
      ! =====================================

      call enforce_contiguous_process_domains( mesh)

    end if ! if (par%master) then

    ! Broadcast the Master's mesh
    call broadcast_mesh( mesh)

    ! == Calculate secondary geometry data
    ! ====================================

    ! Calculate all secondary geometry data
    call calc_all_secondary_mesh_data( mesh, lambda_M, phi_M, beta_stereo)

    ! Calculate all matrix operators
    call calc_all_matrix_operators_mesh( mesh)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_mesh_from_reduced_geometry_singlecore

  subroutine create_mesh_from_reduced_geometry_parallelised( region_name, name, poly_mult_sheet, poly_mult_shelf, &
    p_line_grounding_line, p_line_calving_front, p_line_ice_front, p_line_coastline, &
    xmin, xmax, ymin, ymax, lambda_M, phi_M, beta_stereo, mesh)
    ! Create mesh from the ice geometry lines
    !
    ! Parallelised version

    ! In/output variables:
    character(len=3),           intent(in)        :: region_name
    character(len=256),         intent(in)        :: name
    real(dp), dimension(:,:  ), intent(in)        :: poly_mult_sheet
    real(dp), dimension(:,:  ), intent(in)        :: poly_mult_shelf
    real(dp), dimension(:,:  ), intent(in)        :: p_line_grounding_line
    real(dp), dimension(:,:  ), intent(in)        :: p_line_calving_front
    real(dp), dimension(:,:  ), intent(in)        :: p_line_ice_front
    real(dp), dimension(:,:  ), intent(in)        :: p_line_coastline
    real(dp),                   intent(in)        :: xmin, xmax, ymin, ymax
    real(dp),                   intent(in)        :: lambda_M, phi_M, beta_stereo
    type(type_mesh),            intent(out)       :: mesh

    ! Local variables:
    character(len=256), parameter                 :: routine_name = 'create_mesh_from_reduced_geometry_parallelised'
    real(dp) :: dummy1
    character :: dummy2

    ! Add routine to path
    call init_routine( routine_name)

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
    call crash('whoopsiedaisy!')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_mesh_from_reduced_geometry_parallelised

! == Region of interest mesh refinement
! =====================================

  subroutine refine_mesh_in_regions_of_interest( region_name, poly_mult_sheet, poly_mult_shelf, &
    p_line_grounding_line, p_line_calving_front, p_line_ice_front, p_line_coastline, mesh)
    ! Refine the mesh in the specified regions of interest

    ! In/output variables:
    character(len=3),           intent(in   ) :: region_name
    real(dp), dimension(:,:  ), intent(in   ) :: poly_mult_sheet
    real(dp), dimension(:,:  ), intent(in   ) :: poly_mult_shelf
    real(dp), dimension(:,:  ), intent(in   ) :: p_line_grounding_line
    real(dp), dimension(:,:  ), intent(in   ) :: p_line_calving_front
    real(dp), dimension(:,:  ), intent(in   ) :: p_line_ice_front
    real(dp), dimension(:,:  ), intent(in   ) :: p_line_coastline
    type(type_mesh),            intent(inout) :: mesh

    ! Local variables:
    character(len=256), parameter                 :: routine_name = 'refine_mesh_in_regions_of_interest'
    integer                                       :: i
    character(len=256)                            :: all_names_ROI, name_ROI
    real(dp), dimension(:,:  ), allocatable       :: poly_ROI
    integer                                       :: n1,n2,nn
    real(dp), dimension(:,:  ), allocatable       :: poly

    ! Add routine to path
    call init_routine( routine_name)

    ! if no regions of interest are specified, do nothing
    if (C%choice_regions_of_interest == '') then
      call finalise_routine( routine_name)
      return
    end if

    all_names_ROI = C%choice_regions_of_interest

    do while (.true.)

      ! == Parse list of input ROIs
      ! ===========================

      ! Get the first region of interest from the list
      i = INDEX( all_names_ROI, '||')
      if (i == 0) then
        ! There is only one left in the list
        name_ROI = TRIM( all_names_ROI)
        all_names_ROI = ''
      else
        ! Get the first first one from the list and remove it
        name_ROI = all_names_ROI( 1:i-1)
        all_names_ROI = all_names_ROI( i+2:LEN_TRIM( all_names_ROI))
      end if

      ! == Check validity of requested ROIs
      ! ===================================

      ! Check if current region is indeed defined in the model
      select case (name_ROI)
        case ('')
          ! No region requested: don't need to do anything
          exit
        case ('PineIsland','Thwaites','Amery','RiiserLarsen','SipleCoast', 'LarsenC','TransMounts','DotsonCrosson', & ! Antarctica
              'Narsarsuaq','Nuuk','Jakobshavn','NGIS','Qaanaaq', &                                                    ! Greenland
              'Patagonia', &                                                                                          ! Patagonia
              'Tijn_test_ISMIP_HOM_A','CalvMIP_quarter')                                                              ! Idealised
          ! List of known regions of interest: these pass the test
        case default
          ! Region not found
          call crash('unknown region of interest "' // TRIM( name_ROI) // '"!')
      end select

      ! == Calculate ROIs
      ! =================

      ! Calculate the polygon describing the specified region of interest
      select case (region_name)
        case ('NAM')
          ! North america

          select case (name_ROI)
            case default
              ! Requested area not in this model domain; skip
              cycle
          end select

        case ('EAS')
          ! Eurasia

          select case (name_ROI)
            case default
              ! Requested area not in this model domain; skip
              cycle
          end select

        case ('GRL')
          ! Greenland

          select case (name_ROI)
            case ('Narsarsuaq')
              call calc_polygon_Narsarsuaq( poly_ROI)
            case ('Nuuk')
              call calc_polygon_Nuuk( poly_ROI)
            case ('Jakobshavn')
              call calc_polygon_Jakobshavn( poly_ROI)
            case ('NGIS')
              call calc_polygon_NGIS( poly_ROI)
            case ('Qaanaaq')
              call calc_polygon_Qaanaaq( poly_ROI)
            case default
              ! Requested area not in this model domain; skip
              cycle
          end select

        case ('ANT')

          select case (name_ROI)
            case ('PineIsland')
              call calc_polygon_Pine_Island_Glacier( poly_ROI)
            case ('Thwaites')
              call calc_polygon_Thwaites_Glacier( poly_ROI)
            case ('Amery')
              call calc_polygon_Amery_ice_shelf( poly_ROI)
            case ('RiiserLarsen')
              call calc_polygon_Riiser_Larsen_ice_shelf( poly_ROI)
            case ('SipleCoast')
              call calc_polygon_Siple_Coast( poly_ROI)
            case ('LarsenC')
              call calc_polygon_Larsen_ice_shelf( poly_ROI)
            case ('TransMounts')
              call calc_polygon_Transantarctic_Mountains( poly_ROI)
            case ('DotsonCrosson')
              call calc_polygon_DotsonCrosson_ice_shelf( poly_ROI)
            case ('Patagonia')
              call calc_polygon_Patagonia( poly_ROI)
            case ('Tijn_test_ISMIP_HOM_A')
              call calc_polygon_Tijn_test_ISMIP_HOM_A( poly_ROI)
            case ('CalvMIP_quarter')
              call calc_polygon_CalvMIP_quarter( poly_ROI)
            case default
              ! Requested area not in this model domain; skip
              cycle
          end select

        case default
          call crash('unknown region name "' // region_name // '"!')
      end select

      ! Refine the mesh in the specified region of interest
      ! ===================================================

      ! Uniform
      call refine_mesh_polygon( mesh, poly_ROI, C%ROI_maximum_resolution_uniform, C%alpha_min)

      ! Polygons: ice sheet, ice shelf

      ! Ice sheet
      ! =========

      n1 = 1
      n2 = 0

      do while (n2 < size( poly_mult_sheet,1))

        ! Copy a single polygon from poly_mult
        nn = nint( poly_mult_sheet( n1,1))
        n2 = n1 + nn
        allocate( poly( nn,2))
        poly = poly_mult_sheet( n1+1:n2,:)
        n1 = n2+1

        ! Refine mesh over this single polygon
        call refine_mesh_polygon_ROI( mesh, poly, C%ROI_maximum_resolution_grounded_ice, C%alpha_min, poly_ROI)

        ! Clean up after yourself
        deallocate( poly)

      end do ! do while (n2 < size( poly_mult_sheet,1))

      ! Ice shelf
      ! =========

      n1 = 1
      n2 = 0

      do while (n2 < size( poly_mult_shelf,1))

        ! Copy a single polygon from poly_mult
        nn = nint( poly_mult_shelf( n1,1))
        n2 = n1 + nn
        allocate( poly( nn,2))
        poly = poly_mult_shelf( n1+1:n2,:)
        n1 = n2+1

        ! Refine mesh over this single polygon. Use the ice sheet polygon set as a
        ! "no-refinement" zone to avoid extreme cases where the ice shelf polygon
        ! encompases the ice sheet one (e.g. in the domains of CalvMIP)
        if (C%choice_refgeo_PD_ANT == 'idealised') then
          call refine_mesh_polygon_ROI( mesh, poly, C%ROI_maximum_resolution_floating_ice, C%alpha_min, poly_ROI, poly_mult_sheet)
        else
          call refine_mesh_polygon_ROI( mesh, poly, C%ROI_maximum_resolution_floating_ice, C%alpha_min, poly_ROI)
        end if

        ! Clean up after yourself
        deallocate( poly)

      end do ! do while (n2 < size( poly_mult_sheet,1))

      ! Lines: grounding line, calving front, ice front, coastline
      call refine_mesh_line_ROI( mesh, p_line_grounding_line, C%ROI_maximum_resolution_grounding_line, C%ROI_grounding_line_width, C%alpha_min, poly_ROI)
      call refine_mesh_line_ROI( mesh, p_line_calving_front , C%ROI_maximum_resolution_calving_front , C%ROI_calving_front_width , C%alpha_min, poly_ROI)
      call refine_mesh_line_ROI( mesh, p_line_ice_front     , C%ROI_maximum_resolution_ice_front     , C%ROI_ice_front_width     , C%alpha_min, poly_ROI)
      call refine_mesh_line_ROI( mesh, p_line_coastline     , C%ROI_maximum_resolution_coastline     , C%ROI_coastline_width     , C%alpha_min, poly_ROI)

      ! Clean up after yourself
      deallocate( poly_ROI)

      ! if no names are left, we are finished
      if (all_names_ROI == '') exit

    end do ! do while (.true.)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine refine_mesh_in_regions_of_interest

! == Some useful tools
! ====================

  ! Mesh creation success message
  subroutine write_mesh_success( mesh)
    ! Write the mesh creation success message to the terminal

    use control_resources_and_error_messaging, only: insert_val_into_string_int, insert_val_into_string_dp

    ! In/output variables:
    type(type_mesh),            intent(in)        :: mesh

    ! Local variables:
    character(len=256), parameter                 :: routine_name = 'write_mesh_success'
    character(len=256)                            :: str

    ! Add routine to path
    call init_routine( routine_name)

    str = '     Set up ' // colour_string( TRIM( mesh%name),'light blue') // ' with {int_01} vertices and {int_02} triangles' // &
      ', with a resolution of {dp_01} to {dp_02} m'
    call insert_val_into_string_int( str, '{int_01}', mesh%nV)
    call insert_val_into_string_int( str, '{int_02}', mesh%nTri)
    call insert_val_into_string_dp(  str, '{dp_01}', MINVAL( mesh%R))
    call insert_val_into_string_dp(  str, '{dp_02}', MAXVAL( mesh%R))

    if (par%master) WRITE(0,'(A)') trim( str)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_mesh_success

end module mesh_creation_main
