module bedrock_cumulative_density_functions
  !< Routines for calculating sub-grid bedrock cumulative density functions

  ! The cumulative density function (CDF) tells you for a (configurable)
  ! number of bins, which fraction of the high-resolution bedrock data
  ! inside a Voronoi cell/triangle lies below the value of that bin.
  !
  ! E.g., suppose we use 5 bins. Suppose that, for a certain vertex, the CDF is:
  !
  !   bedrock_CDF( vi,:) = [-1252.0, -1231.6, -1211.5, -1188.5, -1183.4]
  !
  ! The value in the first bin indicates the lowest bedrock elevation encountered
  ! within this Voronoi cell, i.e. 0% of the bedrock is below -1252.0 m. The value
  ! in the next bin then indicates the next interval, i.e. 25% of the bedrock is
  ! below -1231.6. and so on until the last, which indicates the highest bedrock
  ! elevation encountered within this Voronoi cell, i.e. 100% of the bedrock is
  ! below -1183.4.
  !
  ! This can then be used to calculate grounded fractions from ice thicknesses. Suppose
  ! we have an ice thickness in this Voronoi cell of 1378.2 m. The bedrock elevation
  ! where this amount of ice would start to float is equal to:
  !
  !   Hb_float = -Hi * ice_density / seawater_density = ... = -1220.0 m
  !
  ! By looking at the CDF, we can determine that 39.4% of the bedrock in this Voronoi
  ! cell is below this elevation, yielding a grounded fraction of 60.6%.
  !
  ! The bedrock CDF is found by scanning all the high-resolution grid cells (of the
  ! original ice geometry dataset, e.g. BedMachine) that overlap with the Voronoi cell.
  ! The bedrock elevations of all these grid cells are listed, and then sorted
  ! ascendingly. For the example of 5 bins (so intervals of 25%), we'd walk through the
  ! list of elevations until we've passed 25% of the numbers; the elevation we're at
  ! by then goes into the second bin. We move forward again until we've passed 50% of
  ! the numbers; the elevation we're at by that goes into the third bin. Et cetera until
  ! all the bins are filled. Clever, eh?

  use mpi_f08, only: MPI_COMM_WORLD, MPI_BCAST, MPI_DOUBLE_PRECISION
  use mpi_basic, only: par
  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash, colour_string
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use reference_geometry_types, only: type_reference_geometry
  use ice_model_types, only: type_ice_model
  use CSR_sparse_matrix_type, only: type_sparse_matrix_CSR_dp
  use remapping_main, only: Atlas
  use remapping_grid_to_mesh_vertices, only: create_map_from_xy_grid_to_mesh_vertices
  use remapping_grid_to_mesh_triangles, only: create_map_from_xy_grid_to_mesh_triangles
  use petsc_basic, only: mat_petsc2CSR
  use mpi_distributed_memory_grid, only: gather_gridded_data_to_primary
  use netcdf_io_main
  use mpi_distributed_memory, only: distribute_from_primary

  implicit none

  private

  public :: initialise_bedrock_CDFs, calc_bedrock_CDFs

contains

  subroutine initialise_bedrock_CDFs( mesh, refgeo, ice, region_name)
    !< Initialise the sub-grid bedrock cumulative density functions

    ! In/output variables:
    type(type_mesh),               intent(in   ) :: mesh
    type(type_reference_geometry), intent(in   ) :: refgeo
    type(type_ice_model),          intent(inout) :: ice
    character(len=3),              intent(in   ) :: region_name

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_bedrock_CDFs'

    ! Add routine to path
    call init_routine( routine_name)

    if (.not. C%choice_subgrid_grounded_fraction == 'bedrock_CDF' .and. &
        .not. C%choice_subgrid_grounded_fraction == 'bilin_interp_TAF+bedrock_CDF') then
      ! Finalise routine path
      call finalise_routine( routine_name)
      return
    end if

    if (par%primary) then
      write(*,"(A)") '     Initialising the sub-grid bedrock cumulative density functions...'
    end if

    if (C%do_read_bedrock_cdf_from_file) then
      ! Read them from the corresponding mesh file
      call initialise_bedrock_CDFs_from_file( mesh, ice, region_name)
    else
      ! Compute them from scratch
      call calc_bedrock_CDFs( mesh, refgeo, ice)
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_bedrock_CDFs

  subroutine initialise_bedrock_CDFs_from_file( mesh, ice, region_name)
    !< Initialise the velocities for the DIVA solver from an external NetCDF file

    ! In/output variables:
    type(type_mesh),      intent(in   ) :: mesh
    type(type_ice_model), intent(inout) :: ice
    character(len=3),     intent(in   ) :: region_name

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_bedrock_CDFs_from_file'
    character(len=256)             :: filename, check

    ! Add routine to path
    call init_routine( routine_name)

    ! Determine the filename to read for this model region
    select case (region_name)
    case default
      call crash('unknown model region "' // region_name // '"!')
    case ('NAM')
      filename  = C%filename_initial_mesh_NAM
      check = C%choice_initial_mesh_NAM
    case ('EAS')
      filename  = C%filename_initial_mesh_EAS
      check = C%choice_initial_mesh_EAS
    case ('GRL')
      filename  = C%filename_initial_mesh_GRL
      check = C%choice_initial_mesh_GRL
    case ('ANT')
      filename  = C%filename_initial_mesh_ANT
      check = C%choice_initial_mesh_ANT
    end select

    ! Write to terminal
    if (par%primary) write(0,*) '      Reading CDF functions from file "' // colour_string( trim( filename),'light blue') // '"...'

    if (.not. check == 'read_from_file') then
      call crash('The initial mesh was not read from a file. Reading a bedrock CDF this way makes no sense!')
    end if

    ! Read meshed data
    call read_field_from_mesh_file_CDF(   filename, 'bedrock_cdf',   ice%bedrock_cdf   )
    call read_field_from_mesh_file_CDF_b( filename, 'bedrock_cdf_b', ice%bedrock_cdf_b )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_bedrock_CDFs_from_file

  subroutine calc_bedrock_CDFs( mesh, refgeo, ice)
    !< Calculate the sub-grid bedrock cumulative density functions

    ! In/output variables:
    type(type_mesh),               intent(in   ) :: mesh
    type(type_reference_geometry), intent(in   ) :: refgeo
    type(type_ice_model),          intent(inout) :: ice

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_bedrock_CDFs'

    ! Add routine to path
    call init_routine( routine_name)

    if (.not. C%choice_subgrid_grounded_fraction == 'bedrock_CDF' .and. &
        .not. C%choice_subgrid_grounded_fraction == 'bilin_interp_TAF+bedrock_CDF') then
      ! Finalise routine path
      call finalise_routine( routine_name)
      return
    end if

    if (par%primary) write(*,"(A)") '       Calculating bedrock CDFs from initial geometry...'

    ! Calculate CDFs separately on the a-grid (vertices) and the b-grid (triangles)
    call calc_bedrock_CDFs_a( mesh, refgeo, ice)
    call calc_bedrock_CDFs_b( mesh, refgeo, ice)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_bedrock_CDFs

  subroutine calc_bedrock_CDFs_a( mesh, refgeo, ice)
    !< Calculate the sub-grid bedrock cumulative density functions on the a-grid (vertices)

    ! In/output variables:
    type(type_mesh),               intent(in   ) :: mesh
    type(type_reference_geometry), intent(in   ) :: refgeo
    type(type_ice_model),          intent(inout) :: ice

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'calc_bedrock_CDFs_a'
    type(type_sparse_matrix_CSR_dp)       :: M_map
    logical                               :: found_map, found_empty_page
    integer                               :: mi, mi_valid, ierr
    real(dp), dimension(:,:), allocatable :: Hb_grid_tot
    real(dp), dimension(:  ), allocatable :: Hb_list
    integer                               :: vi, k, n, i, j
    integer                               :: n_grid_cells, ii0, ii1
    real(dp)                              :: isc, wii0, wii1

    ! Add routine to path
    call init_routine( routine_name)

    ! Browse the Atlas to see if an appropriate mapping object already exists.
    found_map = .false.
    do mi = 1, size( Atlas, 1)
      if (Atlas( mi)%name_src == refgeo%grid_raw%name .and. Atlas( mi)%name_dst == mesh%name) then
        found_map = .true.
        mi_valid  = mi
        exit
      end if
    end do

    ! if no appropriate mapping object could be found, create one.
    if (.not. found_map) then
      found_empty_page = .false.
      do mi = 1, size( Atlas,1)
        if (.not. Atlas( mi)%is_in_use) then
          found_empty_page = .true.
          call create_map_from_xy_grid_to_mesh_vertices( refgeo%grid_raw, mesh, C%output_dir, Atlas( mi))
          mi_valid = mi
          exit
        end if
      end do
      ! Safety
      if (.not. found_empty_page) call crash('No more room in Atlas - assign more memory!')
    end if

    ! Turn map into CSR, so we can use it here
    call mat_petsc2CSR( Atlas( mi_valid)%M, M_map)

    ! Get complete gridded bedrock elevation on all processes
    allocate( Hb_grid_tot( refgeo%grid_raw%nx, refgeo%grid_raw%ny))
    call gather_gridded_data_to_primary( refgeo%grid_raw, refgeo%Hb_grid_raw, Hb_grid_tot)
    call MPI_BCAST( Hb_grid_tot(:,:), refgeo%grid_raw%nx * refgeo%grid_raw%ny, MPI_doUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! allocate memory for list of bedrock elevations
    allocate( Hb_list( refgeo%grid_raw%nx * refgeo%grid_raw%ny ))
    Hb_list = 0._dp

    ! Initialise cumulative density function (CDF)
    ice%bedrock_cdf = 0._dp

    do vi = mesh%vi1, mesh%vi2

      ! Clear the list
      Hb_list = 0._dp

      ! List bedrock elevations from all grid cells overlapping with this vertex's
      ! Voronoi cell (as already determined by the remapping operator)
      ! ==============================================================

      n_grid_cells = 0
      do k = M_map%ptr( vi), M_map%ptr( vi+1)-1
        n_grid_cells = n_grid_cells + 1
        n = M_map%ind( k)
        i = refgeo%grid_raw%n2ij( n,1)
        j = refgeo%grid_raw%n2ij( n,2)
        Hb_list( n_grid_cells) = Hb_grid_tot( i,j)
      end do

      ! Safety
      if (n_grid_cells == 0) call crash('found no overlapping grid cells!')

      ! === Cumulative density function ===
      ! ===================================

      ! Inefficient but easy sorting of Hb_list
      do i = 1, n_grid_cells-1
      do j = i+1, n_grid_cells
        if (Hb_list( i) > Hb_list( j)) then
          Hb_list( i) = Hb_list( i) + Hb_list( j)
          Hb_list( j) = Hb_list( i) - Hb_list( j)
          Hb_list( i) = Hb_list( i) - Hb_list( j)
        end if
      end do
      end do

      ! Set first (0%) and last bins (100%) of the CDF to the minimum
      ! and maximum bedrock elevations scanned, respectively
      ice%bedrock_cdf( vi, 1                          ) = Hb_list( 1)
      ice%bedrock_cdf( vi, C%subgrid_bedrock_cdf_nbins) = Hb_list( n_grid_cells)

      ! Compute the bedrock elevation for each of the other CDF bins
      do i = 2, C%subgrid_bedrock_cdf_nbins - 1
        isc  = 1._dp + (real( n_grid_cells,dp) - 1._dp) * (real( i,dp) - 1._dp) / real( C%subgrid_bedrock_cdf_nbins - 1,dp)
        ii0  = floor(   isc)
        ii1  = ceiling( isc)
        wii0 = real( ii1,dp) - isc
        wii1 = 1.0 - wii0
        ice%bedrock_cdf( vi,i) = wii0 * Hb_list( ii0) + wii1 * Hb_list( ii1)
      end do

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_bedrock_CDFs_a

  subroutine calc_bedrock_CDFs_b( mesh, refgeo, ice)
    !< Calculate the sub-grid bedrock cumulative density functions on the b-grid (triangles)

    ! In/output variables:
    type(type_mesh),               intent(in   ) :: mesh
    type(type_reference_geometry), intent(in   ) :: refgeo
    type(type_ice_model),          intent(inout) :: ice

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'calc_bedrock_CDFs_b'
    type(type_sparse_matrix_CSR_dp)       :: M_map
    logical                               :: found_map, found_empty_page
    integer                               :: mi, mi_valid, ierr
    real(dp), dimension(:,:), allocatable :: Hb_grid_tot
    real(dp), dimension(:  ), allocatable :: Hb_list
    integer                               :: ti, k, n, i, j
    integer                               :: n_grid_cells, ii0, ii1
    real(dp)                              :: isc, wii0, wii1

    ! Add routine to path
    call init_routine( routine_name)

    ! Browse the Atlas to see if an appropriate mapping object already exists.
    found_map = .false.
    do mi = 1, size( Atlas, 1)
      if (Atlas( mi)%name_src == refgeo%grid_raw%name .and. Atlas( mi)%name_dst == (trim( mesh%name) // '_triangles')) then
        found_map = .true.
        mi_valid  = mi
        exit
      end if
    end do

    ! if no appropriate mapping object could be found, create one.
    if (.not. found_map) then
      found_empty_page = .false.
      do mi = 1, size( Atlas,1)
        if (.not. Atlas( mi)%is_in_use) then
          found_empty_page = .true.
          call create_map_from_xy_grid_to_mesh_triangles( refgeo%grid_raw, mesh, C%output_dir, Atlas( mi))
          mi_valid = mi
          exit
        end if
      end do
      ! Safety
      if (.not. found_empty_page) call crash('No more room in Atlas - assign more memory!')
    end if

    ! Turn map into CSR, so we can use it here
    call mat_petsc2CSR( Atlas( mi_valid)%M, M_map)

    ! Get complete gridded bedrock elevation on all processes
    allocate( Hb_grid_tot( refgeo%grid_raw%nx, refgeo%grid_raw%ny))
    call gather_gridded_data_to_primary( refgeo%grid_raw, refgeo%Hb_grid_raw, Hb_grid_tot)
    call MPI_BCAST( Hb_grid_tot(:,:), refgeo%grid_raw%nx * refgeo%grid_raw%ny, MPI_doUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! allocate memory for list of bedrock elevations
    allocate( Hb_list( refgeo%grid_raw%nx * refgeo%grid_raw%ny ))
    Hb_list = 0._dp

    ! Initialise cumulative density function (CDF)
    ice%bedrock_cdf_b = 0._dp

    do ti = mesh%ti1, mesh%ti2

      ! Clear the list
      Hb_list = 0._dp

      ! List bedrock elevations from all grid cells overlapping with this vertex's
      ! Voronoi cell (as already determined by the remapping operator)
      ! ==============================================================

      n_grid_cells = 0
      do k = M_map%ptr( ti), M_map%ptr( ti+1)-1
        n_grid_cells = n_grid_cells + 1
        n = M_map%ind( k)
        i = refgeo%grid_raw%n2ij( n,1)
        j = refgeo%grid_raw%n2ij( n,2)
        Hb_list( n_grid_cells) = Hb_grid_tot( i,j)
      end do

      ! Safety
      if (n_grid_cells == 0) call crash('found no overlapping grid cells!')

      ! === Cumulative density function ===
      ! ===================================

      ! Inefficient but easy sorting of Hb_list
      do i = 1, n_grid_cells-1
      do j = i+1, n_grid_cells
        if (Hb_list( i) > Hb_list( j)) then
          Hb_list( i) = Hb_list( i) + Hb_list( j)
          Hb_list( j) = Hb_list( i) - Hb_list( j)
          Hb_list( i) = Hb_list( i) - Hb_list( j)
        end if
      end do
      end do

      ! Set first (0%) and last bins (100%) of the CDF to the minimum
      ! and maximum bedrock elevations scanned, respectively
      ice%bedrock_cdf_b( ti, 1                          ) = Hb_list( 1)
      ice%bedrock_cdf_b( ti, C%subgrid_bedrock_cdf_nbins) = Hb_list( n_grid_cells)

      ! Compute the bedrock elevation for each of the other CDF bins
      do i = 2, C%subgrid_bedrock_cdf_nbins - 1
        isc  = 1._dp + (real( n_grid_cells,dp) - 1._dp) * (real( i,dp) - 1._dp) / real( C%subgrid_bedrock_cdf_nbins - 1,dp)
        ii0  = floor(   isc)
        ii1  = ceiling( isc)
        wii0 = real( ii1,dp) - isc
        wii1 = 1.0 - wii0
        ice%bedrock_cdf_b( ti,i) = wii0 * Hb_list( ii0) + wii1 * Hb_list( ii1)
      end do

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_bedrock_CDFs_b

  subroutine read_field_from_mesh_file_CDF( filename, field_name_options, &
    d_mesh_partial)
    !< Read a cumulative density function field from a NetCDF file on a mesh

    ! NOTE: the mesh should be read before, and memory allocated for d_mesh_partial!

    ! In/output variables:
    character(len=*),                   intent(in   ) :: filename
    character(len=*),                   intent(in   ) :: field_name_options
    real(dp), dimension(:,:),           intent(  out) :: d_mesh_partial

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'read_field_from_mesh_file_CDF'
    integer                               :: ncid
    integer                               :: id_dim_vi, id_var, id_dim_bins, nbins_loc, nV_loc
    character(len=1024)                   :: var_name
    real(dp), dimension(:,:), allocatable :: d_mesh

    ! Add routine to path
    call init_routine( routine_name)

    ! == Read grid and data from file
    ! ===============================

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_reading( filename, ncid)

    ! Get number of mesh vertices
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nV, id_dim_vi, dim_length = nV_loc)

    ! Check that number of bins in file match the ones in the config file
    call inquire_dim_multopt( filename, ncid, 'bin', id_dim_bins, dim_length = nbins_loc)
    if (nbins_loc /= C%subgrid_bedrock_cdf_nbins) then
      call crash('number of CDF bins in external file ({int_01}) does not match subgrid_bedrock_cdf_nbins in config file!', int_01 = nbins_loc)
    end if

    ! Look for the specified variable in the file
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('couldnt find any of the options "' // trim( field_name_options) // '" in file "' // trim( filename)  // '"!')

    ! allocate memory
    if (par%primary) allocate( d_mesh( nV_loc, C%subgrid_bedrock_cdf_nbins))

    ! Read data from file
    call read_var_primary( filename, ncid, id_var, d_mesh)

    ! Close the NetCDF file
    call close_netcdf_file( ncid)

    ! == Distribute gridded data from the primary to all processes in partial vector form
    ! ==================================================================================

    ! Distribute data
    call distribute_from_primary( d_mesh, d_mesh_partial)

    ! Clean up after yourself
    if (par%primary) deallocate( d_mesh)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_field_from_mesh_file_CDF

  subroutine read_field_from_mesh_file_CDF_b( filename, field_name_options, &
    d_mesh_partial)
    !< Read a cumulative density function field from a NetCDF file on a mesh b-grid

    ! NOTE: the mesh should be read before, and memory allocated for d_mesh_partial!

    ! In/output variables:
    character(len=*),                   intent(in   ) :: filename
    character(len=*),                   intent(in   ) :: field_name_options
    real(dp), dimension(:,:),           intent(  out) :: d_mesh_partial

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'read_field_from_mesh_file_CDF_b'
    integer                                 :: ncid
    integer                                 :: id_dim_ti, id_var, id_dim_bins, nbins_loc, nTri_loc
    character(len=1024)                     :: var_name
    real(dp), dimension(:,:  ), allocatable :: d_mesh

    ! Add routine to path
    call init_routine( routine_name)

    ! == Read grid and data from file
    ! ===============================

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_reading( filename, ncid)

    ! Get number of mesh triangles
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nTri, id_dim_ti, dim_length = nTri_loc)

    ! Check that number of bins in file match the ones in the config file
    call inquire_dim_multopt( filename, ncid, 'bin', id_dim_bins, dim_length = nbins_loc)
    if (nbins_loc /= C%subgrid_bedrock_cdf_nbins) then
      call crash('number of CDF bins in external file ({int_01}) does not match subgrid_bedrock_cdf_nbins in config file!', int_01 = nbins_loc)
    end if

    ! Look for the specified variable in the file
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('couldnt find any of the options "' // trim( field_name_options) // '" in file "' // trim( filename)  // '"!')

    ! allocate memory
    if (par%primary) allocate( d_mesh( nTri_loc, C%subgrid_bedrock_cdf_nbins))

    ! Read data from file
    call read_var_primary( filename, ncid, id_var, d_mesh)

    ! Close the NetCDF file
    call close_netcdf_file( ncid)

    ! == Distribute gridded data from the primary to all processes in partial vector form
    ! ==================================================================================

    ! Distribute data
    call distribute_from_primary( d_mesh, d_mesh_partial)

    ! Clean up after yourself
    if (par%primary) deallocate( d_mesh)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_field_from_mesh_file_CDF_b

end module bedrock_cumulative_density_functions
