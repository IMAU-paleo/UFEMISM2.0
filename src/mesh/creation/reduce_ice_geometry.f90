module reduce_ice_geometry

  ! "Reduce" an ice geometry: calculate the grounding line, calving front,
  ! ice front, ice-sheet polygon and ice-shelf polygon from an ice geometry
  ! defined on either a grid or a mesh; these can then be used to create a new mesh.

  use mpi_f08, only: MPI_COMM_WORLD, MPI_BCAST, MPI_INTEGER, MPI_DOUBLE_PRECISION
  use mpi_basic, only: par, sync
  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use grid_types, only: type_grid
  use mesh_types, only: type_mesh
  use mpi_distributed_memory, only: gather_to_primary
  use grid_basic, only: calc_grid_mask_as_polygons, calc_grid_contour_as_line
  use mpi_distributed_memory_grid, only: gather_gridded_data_to_primary
  use ice_geometry_basics, only: thickness_above_floatation
  use mesh_utilities, only: calc_mesh_mask_as_polygons, calc_mesh_contour_as_line

  implicit none

  private

  public :: reduce_gridded_ice_geometry, reduce_meshed_ice_geometry

contains

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
    integer                                 :: ierr
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
    if (par%primary) then
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

    ! Gather ice geometry data in grid form to the primary
    call gather_gridded_data_to_primary( grid, Hi, Hi_grid)
    call gather_gridded_data_to_primary( grid, Hb, Hb_grid)
    call gather_gridded_data_to_primary( grid, Hs, Hs_grid)
    call gather_gridded_data_to_primary( grid, SL, SL_grid)

    ! Let the primary calculate the polygons and lines

    if (par%primary) then

      ! Allocate memory for thickness above floatation and Hb-SL
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
      allocate( mask_sheet( grid%nx, grid%ny), source = .false.)
      allocate( mask_shelf( grid%nx, grid%ny), source = .false.)

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
      allocate( mask_calc_grounding_line( grid%nx, grid%ny), source = .false.)
      allocate( mask_calc_calving_front(  grid%nx, grid%ny), source = .false.)
      allocate( mask_calc_ice_front(      grid%nx, grid%ny), source = .false.)
      allocate( mask_calc_coastline(      grid%nx, grid%ny), source = .false.)

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

    end if

    ! Broadcast polygon sizes to all processes
    call MPI_BCAST( n_poly_mult_sheet, 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST( n_poly_mult_shelf, 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)

    ! allocate memory for the polygons on the other processes
    if (.not. par%primary) then
      allocate( poly_mult_sheet( n_poly_mult_sheet,2))
      allocate( poly_mult_shelf( n_poly_mult_shelf,2))
    end if

    ! Broadcast polygons to all processes
    call MPI_BCAST( poly_mult_sheet(:,:), n_poly_mult_sheet*2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST( poly_mult_shelf(:,:), n_poly_mult_shelf*2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! Broadcast line sizes to all processes
    call MPI_BCAST( n_line_grounding_line, 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST( n_line_calving_front , 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST( n_line_ice_front     , 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST( n_line_coastline     , 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)

    ! Allocate memory for the lines on the other processes
    if (.not. par%primary) then
      allocate( p_line_grounding_line( n_line_grounding_line,4))
      allocate( p_line_calving_front(  n_line_calving_front ,4))
      allocate( p_line_ice_front(      n_line_ice_front     ,4))
      allocate( p_line_coastline(      n_line_coastline     ,4))
    end if

    ! Broadcast lines to all processes
    call MPI_BCAST( p_line_grounding_line(:,:), n_line_grounding_line*4, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST( p_line_calving_front(:,:) , n_line_calving_front *4, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST( p_line_ice_front(:,:)     , n_line_ice_front     *4, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST( p_line_coastline(:,:)     , n_line_coastline     *4, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

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
    character(len=256), parameter           :: routine_name = 'reduce_meshed_ice_geometry'
    integer                                 :: ierr
    real(dp), dimension(:    ), allocatable :: Hi_tot
    real(dp), dimension(:    ), allocatable :: Hb_tot
    real(dp), dimension(:    ), allocatable :: Hs_tot
    real(dp), dimension(:    ), allocatable :: SL_tot
    real(dp), dimension(:    ), allocatable :: TAF_tot
    real(dp), dimension(:    ), allocatable :: Hb_minus_SL_tot
    integer                                 :: vi
    logical,  dimension(:    ), allocatable :: mask_sheet
    logical,  dimension(:    ), allocatable :: mask_shelf
    logical,  dimension(:    ), allocatable :: mask_calc_grounding_line
    logical,  dimension(:    ), allocatable :: mask_calc_calving_front
    logical,  dimension(:    ), allocatable :: mask_calc_coastline
    integer                                 :: n_poly_mult_sheet
    integer                                 :: n_poly_mult_shelf
    integer                                 :: n_line_grounding_line
    integer                                 :: n_line_calving_front
    integer                                 :: n_line_ice_front
    integer                                 :: n_line_coastline

    ! Add routine to path
    call init_routine( routine_name)

    ! allocate memory for gathered ice geometry data
    if (par%primary) then
      allocate( Hi_tot( mesh%nV))
      allocate( Hb_tot( mesh%nV))
      allocate( Hs_tot( mesh%nV))
      allocate( SL_tot( mesh%nV))
    end if

    ! Gather ice geometry data in grid form to the primary
    call gather_to_primary( Hi, Hi_tot)
    call gather_to_primary( Hb, Hb_tot)
    call gather_to_primary( Hs, Hs_tot)
    call gather_to_primary( SL, SL_tot)

    ! Let the primary calculate the polygons and lines

    if (par%primary) then

      ! allocate memory for thickness above floatation and Hb-SL
      allocate( TAF_tot(         mesh%nV))
      allocate( Hb_minus_SL_tot( mesh%nV))

      ! Calculate thickness above floatation and Hb-SL
      do vi = 1, mesh%nV
        TAF_tot(         vi) = thickness_above_floatation( Hi_tot( vi), Hb_tot( vi), SL_tot( vi))
        Hb_minus_SL_tot( vi) = Hb_tot( vi) - SL_tot( vi)
      end do

      ! Fill in masks for floating/grounded ice
      allocate( mask_sheet( mesh%nV), source = .false.)
      allocate( mask_shelf( mesh%nV), source = .false.)

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
      allocate( mask_calc_grounding_line( mesh%nV), source = .false.)
      allocate( mask_calc_calving_front ( mesh%nV), source = .false.)
      allocate( mask_calc_coastline     ( mesh%nV), source = .false.)

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

    end if

    ! Broadcast polygon sizes to all processes
    call MPI_BCAST( n_poly_mult_sheet, 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST( n_poly_mult_shelf, 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)

    ! allocate memory for the polygons on the other processes
    if (.not. par%primary) then
      allocate( poly_mult_sheet( n_poly_mult_sheet,2))
      allocate( poly_mult_shelf( n_poly_mult_shelf,2))
    end if

    ! Broadcast polygons to all processes
    call MPI_BCAST( poly_mult_sheet(:,:), n_poly_mult_sheet*2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST( poly_mult_shelf(:,:), n_poly_mult_shelf*2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! Broadcast line sizes to all processes
    call MPI_BCAST( n_line_grounding_line, 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST( n_line_calving_front , 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST( n_line_ice_front     , 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST( n_line_coastline     , 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)

    ! allocate memory for the lines on the other processes
    if (.not. par%primary) then
      allocate( p_line_grounding_line( n_line_grounding_line,4))
      allocate( p_line_calving_front(  n_line_calving_front ,4))
      allocate( p_line_ice_front(      n_line_ice_front     ,4))
      allocate( p_line_coastline(      n_line_coastline     ,4))
    end if

    ! Broadcast lines to all processes
    call MPI_BCAST( p_line_grounding_line(:,:), n_line_grounding_line*4, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST( p_line_calving_front(:,:) , n_line_calving_front *4, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST( p_line_ice_front(:,:)     , n_line_ice_front     *4, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST( p_line_coastline(:,:)     , n_line_coastline     *4, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine reduce_meshed_ice_geometry

end module reduce_ice_geometry
