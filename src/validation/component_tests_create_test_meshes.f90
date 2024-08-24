module component_tests_create_test_meshes

  ! Create the suite of test meshes for the component tests.

  use mpi
  use precisions                                             , only: dp
  use mpi_basic                                              , only: par, cerr, ierr, recv_status, sync
  use control_resources_and_error_messaging                  , only: warning, crash, happy, init_routine, finalise_routine, colour_string
  use model_configuration                                    , only: C
  use assertions_unit_tests, only: ASSERTION, test_mesh_is_self_consistent
  use mesh_types                                             , only: type_mesh
  use mesh_memory                                            , only: allocate_mesh_primary
  use mesh_creation                                          , only: initialise_dummy_mesh
  use mesh_refinement                                        , only: refine_mesh_uniform, Lloyds_algorithm_single_iteration, refine_mesh_polygon
  use mesh_secondary                                         , only: calc_all_secondary_mesh_data
  use mesh_operators                                         , only: calc_all_matrix_operators_mesh, calc_3D_matrix_operators_mesh, calc_3D_gradient_bk_bk
  use netcdf_basic                                           , only: create_new_netcdf_file_for_writing, close_netcdf_file
  use netcdf_output                                          , only: setup_mesh_in_netcdf_file, add_field_mesh_dp_2D_notime, add_field_mesh_dp_2D_b_notime, &
                                                                     add_field_mesh_dp_2D_c_notime, write_to_field_multopt_mesh_dp_2D_notime, &
                                                                     write_to_field_multopt_mesh_dp_2D_b_notime, write_to_field_multopt_mesh_dp_2D_c_notime, &
                                                                     setup_xy_grid_in_netcdf_file, add_field_grid_dp_2D_notime, &
                                                                     write_to_field_multopt_grid_dp_2D_notime, add_zeta_dimension_to_file, &
                                                                     add_field_mesh_dp_3D_b_notime, write_to_field_multopt_mesh_dp_3D_notime, &
                                                                     write_to_field_multopt_mesh_dp_3D_b_notime

  implicit none

  private

  public :: create_all_test_meshes, read_list_of_test_mesh_filenames

contains

  !> Create the suite of test meshes for the component tests.
  subroutine create_all_test_meshes( foldername_test_meshes, filename_test_meshes_list)

    ! In/output variables:
    character(len=*), intent(out) :: foldername_test_meshes
    character(len=*), intent(out) :: filename_test_meshes_list

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_all_test_meshes'

    ! Add routine to call stack
    call init_routine( routine_name)

    if (par%master) write(0,*) '  Creating suite of test meshes for component tests...'
    if (par%master) write(0,*) ''

    ! Create output folder within the component tests folder to store the meshes
    call create_component_tests_mesh_output_folder( foldername_test_meshes, filename_test_meshes_list)

    ! Create all the test meshes for the Antarctica domain
    call create_all_test_meshes_Antarctica( foldername_test_meshes, filename_test_meshes_list)

    if (par%master) write(0,*) ''

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_all_test_meshes

  !> Create the output folder for the mesh component tests
  subroutine create_component_tests_mesh_output_folder( foldername_test_meshes, filename_test_meshes_list)

    ! In/output variables:
    character(len=*), intent(out) :: foldername_test_meshes
    character(len=*), intent(out) :: filename_test_meshes_list

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_component_tests_mesh_output_folder'
    logical                        :: ex
    integer                        :: io_filename_test_meshes_list, stat
    character(len=1024)            :: msg

    ! Add routine to path
    call init_routine( routine_name)

    foldername_test_meshes = trim(C%output_dir) // '/mesh'

    if (par%master) then

      ! Remove existing folder if necessary
      inquire( file = trim( foldername_test_meshes) // '/.', exist = ex)
      if (ex) then
        call system('rm -rf ' // trim( foldername_test_meshes))
      end if

      ! Create the directory
      call system('mkdir ' // trim( foldername_test_meshes))

      ! Create the text file listing the filenames of the test meshes
      filename_test_meshes_list = trim(foldername_test_meshes) // '/list_of_test_meshes.txt'
      ! Create file
      open(newunit = io_filename_test_meshes_list, file = filename_test_meshes_list, &
        status = "new", action = "write", iostat = stat, iomsg = msg)
      if (stat /= 0) then
        call crash('Could not create text file listing the filenames of the test meshes, error message "' // trim(msg) // '"')
      end if
      close(io_filename_test_meshes_list)

    end if
    call MPI_BCAST( foldername_test_meshes,    len(foldername_test_meshes   ), MPI_CHAR, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST( filename_test_meshes_list, len(filename_test_meshes_list), MPI_CHAR, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_component_tests_mesh_output_folder

  ! Read the list of test mesh filenames from the text file
  subroutine read_list_of_test_mesh_filenames( filename_test_meshes_list, test_mesh_filenames)

    ! In/output variables:
    character(len=*),                               intent(in ) :: filename_test_meshes_list
    character(len=1024), dimension(:), allocatable, intent(out) :: test_mesh_filenames
    integer             :: i, io_filename_test_meshes_list, stat
    character(len=1024) :: msg
    integer             :: number_of_meshes
    character(len=1024) :: mesh_name

    ! Local variables:

    do i = 0, par%n-1
      if (i == par%i) then

        ! Determine number of test meshes
        ! ===============================

        ! Open text file listing the filenames of the test meshes
        open(newunit = io_filename_test_meshes_list, file = filename_test_meshes_list, &
          status = "old", action = "read", iostat = stat, iomsg = msg)
        if (stat /= 0) then
          call crash('Could not create text file listing the filenames of the test meshes, error message "' // trim(msg) // '"')
        end if

        number_of_meshes = 0
        do while (.true.)
          read(unit = io_filename_test_meshes_list, fmt = '(a)', &
            iostat = stat) mesh_name
          if (stat == 0) then
            number_of_meshes = number_of_meshes + 1
          else
            exit
          end if
        end do

        ! Close text file listing the filenames of the test meshes
        close(io_filename_test_meshes_list)

        ! Allocate memory for list of filenames
        ! =====================================

        allocate( test_mesh_filenames( number_of_meshes))

        ! Read the file again, this time actually saving the mesh names
        ! =============================================================

        ! Open text file listing the filenames of the test meshes
        open(newunit = io_filename_test_meshes_list, file = filename_test_meshes_list, &
        status = "old", action = "read", iostat = stat, iomsg = msg)
      if (stat /= 0) then
        call crash('Could not create text file listing the filenames of the test meshes, error message "' // trim(msg) // '"')
      end if

      number_of_meshes = 0
      do while (.true.)
        read(unit = io_filename_test_meshes_list, fmt = '(a)', &
          iostat = stat) mesh_name
        if (stat == 0) then
          number_of_meshes = number_of_meshes + 1
          test_mesh_filenames( number_of_meshes) = trim(adjustl(mesh_name))//'.nc'
        else
          exit
        end if
      end do

      ! Close text file listing the filenames of the test meshes
      close(io_filename_test_meshes_list)

      end if
      call sync
    end do

  end subroutine read_list_of_test_mesh_filenames

  !> Create all the test meshes for the Antarctica domain
  subroutine create_all_test_meshes_Antarctica( foldername_test_meshes, filename_test_meshes_list)

    ! In/output variables:
    character(len=*), intent(in) :: foldername_test_meshes
    character(len=*), intent(in) :: filename_test_meshes_list

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'create_all_test_meshes_Antarctica'
    character(len=1024), parameter      :: domain_name = 'Ant'
    real(dp), parameter                 :: xmin                  = -3040E3_dp
    real(dp), parameter                 :: xmax                  =  3040E3_dp
    real(dp), parameter                 :: ymin                  = -3040E3_dp
    real(dp), parameter                 :: ymax                  =  3040E3_dp
    real(dp), parameter                 :: lambda_M              = 0._dp
    real(dp), parameter                 :: phi_M                 = -90._dp
    real(dp), parameter                 :: beta_stereo           = 71._dp
    real(dp), parameter                 :: alpha_min             = 0.4363_dp
    integer, parameter                  :: nit_Lloyds_algorithm  = 2

    real(dp)                            :: res_min               = 400e3_dp
    real(dp)                            :: res_max               = 75e3_dp
    integer, dimension(4), parameter    :: nits_Lloyds_algorithm = [4,6,8,10]
    real(dp), dimension(6), parameter   :: uniform_resolutions   = [400e3_dp, 300e3_dp, 200e3_dp, 150e3_dp, 100e3_dp, 75e3_dp]
    real(dp), parameter                 :: uniform_resolution    = 150e3_dp

    ! real(dp)                            :: res_min               = 800e3_dp
    ! real(dp)                            :: res_max               = 400e3_dp
    ! integer, dimension(2), parameter    :: nits_Lloyds_algorithm = [4,6]
    ! real(dp), dimension(3), parameter   :: uniform_resolutions   = [800e3_dp, 600e3_dp, 400e3_dp]
    ! real(dp), parameter                 :: uniform_resolution    = 400e3_dp

    integer                             :: i
    character(len=1)                    :: orientation

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Create a set of test meshes covering the unit square with a uniform resolution
    do i = 1, size( uniform_resolutions)
      call create_test_mesh_uniform( xmin, xmax, ymin, ymax, lambda_M, phi_M, beta_stereo, &
        uniform_resolutions( i), alpha_min, nit_Lloyds_algorithm, domain_name, &
        foldername_test_meshes, filename_test_meshes_list)
    end do

    ! Create a set of test meshes covering the unit square with a resolution gradient
    orientation          = 'x'
    call create_test_mesh_gradient( xmin, xmax, ymin, ymax, lambda_M, phi_M, beta_stereo, &
      res_min, res_max, alpha_min, nit_Lloyds_algorithm, orientation, domain_name, &
      foldername_test_meshes, filename_test_meshes_list)
    orientation          = 'y'
    call create_test_mesh_gradient( xmin, xmax, ymin, ymax, lambda_M, phi_M, beta_stereo, &
      res_min, res_max, alpha_min, nit_Lloyds_algorithm, orientation, domain_name, &
      foldername_test_meshes, filename_test_meshes_list)

    ! Create a set of otherwise identical test meshes with increasing smoothness
    do i = 1, size( nits_Lloyds_algorithm)
      call create_test_mesh_uniform( xmin, xmax, ymin, ymax, lambda_M, phi_M, beta_stereo, &
        uniform_resolution, alpha_min, nits_Lloyds_algorithm( i), domain_name, &
        foldername_test_meshes, filename_test_meshes_list)
    end do

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_all_test_meshes_Antarctica

  !> Create a test mesh covering the unit square with a uniform resolution
  subroutine create_test_mesh_uniform( xmin, xmax, ymin, ymax, lambda_M, phi_M, beta_stereo, &
    res_max, alpha_min, nit_Lloyds_algorithm, domain_name, foldername_test_meshes, filename_test_meshes_list)

    ! In/output variables:
    real(dp),         intent(in) :: xmin, xmax, ymin, ymax
    real(dp),         intent(in) :: lambda_M, phi_M, beta_stereo
    real(dp),         intent(in) :: res_max
    real(dp),         intent(in) :: alpha_min
    integer,          intent(in) :: nit_Lloyds_algorithm
    character(len=*), intent(in) :: domain_name
    character(len=*), intent(in) :: foldername_test_meshes
    character(len=*), intent(in) :: filename_test_meshes_list

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_test_mesh_uniform'
    type(type_mesh)                :: mesh
    character(len=14)              :: res_max_str
    character(len=3)               :: nit_Lloyds_algorithm_str
    character(len=1024)            :: mesh_name
    integer                        :: it_Lloyds_algorithm
    integer                        :: ncid
    character(len=1024)            :: filename
    integer                        :: io_filename_test_meshes_list, stat
    character(len=1024)            :: msg

    ! Add routine to call stack
    call init_routine( routine_name)

    write( res_max_str             ,'(es14.4)') res_max
    write( nit_Lloyds_algorithm_str,'(i3)')     nit_Lloyds_algorithm
    mesh_name = 'comp_test_mesh_'//trim(domain_name)//'_uniform_'//trim(adjustl(res_max_str))//'_m'//&
      '_nit_Lloyd_'//trim(adjustl(nit_Lloyds_algorithm_str))

    ! Allocate memory
    call allocate_mesh_primary( mesh, trim(mesh_name), 1000, 2000, 32)

    ! Initialise the dummy mesh
    call initialise_dummy_mesh( mesh, xmin, xmax, ymin, ymax)

    ! Refine the mesh
    call refine_mesh_uniform( mesh, res_max, alpha_min)
    call test_mesh_is_self_consistent( mesh, ASSERTION, &
      trim(mesh_name)//': self_consistency (after refine_mesh_uniform)')

    ! Smooth the mesh
    do it_Lloyds_algorithm = 1, nit_Lloyds_algorithm
      call Lloyds_algorithm_single_iteration( mesh, alpha_min)
      call test_mesh_is_self_consistent( mesh, ASSERTION, &
        trim(mesh_name)//': self_consistency (after Lloyds_algorithm_single_iteration)')
    end do

    ! Calculate secondary geometry data (needed in order to be able to write to NetCDF)
    call calc_all_secondary_mesh_data( mesh, lambda_M, phi_M, beta_stereo)

    ! Write to NetCDF file
    filename = trim( foldername_test_meshes)//'/'//trim( mesh_name)//'.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)
    call close_netcdf_file( ncid)

    ! Add name of mesh to list file
    if (par%master) then
      open(newunit = io_filename_test_meshes_list, file = filename_test_meshes_list, &
        status = "old", action = "write", position = "append", iostat = stat, iomsg = msg)
      if (stat /= 0) then
        call crash('Could not open text file listing the filenames of the test meshes, error message "' // trim(msg) // '"')
      end if
      write(io_filename_test_meshes_list,*) trim(mesh_name)
      close(io_filename_test_meshes_list)
    end if

    if (par%master) write(0,*) '    Created test mesh ', colour_string( trim( mesh_name), 'light blue')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_test_mesh_uniform

  !> Create a test mesh covering the unit square with a resolution gradient
  subroutine create_test_mesh_gradient( xmin, xmax, ymin, ymax, lambda_M, phi_M, beta_stereo, &
    res_min, res_max, alpha_min, nit_Lloyds_algorithm, orientation, domain_name, &
    foldername_test_meshes, filename_test_meshes_list)

    ! In/output variables:
    real(dp),         intent(in) :: xmin, xmax, ymin, ymax
    real(dp),         intent(in) :: lambda_M, phi_M, beta_stereo
    real(dp),         intent(in) :: res_min
    real(dp),         intent(in) :: res_max
    real(dp),         intent(in) :: alpha_min
    integer,          intent(in) :: nit_Lloyds_algorithm
    character(len=*), intent(in) :: orientation
    character(len=*), intent(in) :: domain_name
    character(len=*), intent(in) :: foldername_test_meshes
    character(len=*), intent(in) :: filename_test_meshes_list

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_test_mesh_unit_square_gradient'
    type(type_mesh)                :: mesh
    character(len=14)              :: res_min_str, res_max_str
    character(len=1024)            :: mesh_name
    integer                        :: n,i
    real(dp)                       :: log_res_min, log_res_max, log_res, resolution
    real(dp)                       :: x_lo, y_lo, x_hi, y_hi
    real(dp), dimension(4,2)       :: poly
    integer                        :: it_Lloyds_algorithm
    integer                        :: ncid
    character(len=1024)            :: filename
    integer                        :: io_filename_test_meshes_list, stat
    character(len=1024)            :: msg

    ! Add routine to call stack
    call init_routine( routine_name)

    write( res_min_str,'(es14.4)') res_min
    write( res_max_str,'(es14.4)') res_max
    mesh_name = 'comp_test_mesh_'//trim(domain_name)//'_gradient_'//trim(adjustl(res_min_str))//&
      '-'//trim(adjustl(res_max_str))//'_m_'//orientation

    ! Allocate memory
    call allocate_mesh_primary( mesh, trim(mesh_name), 1000, 2000, 32)

    ! Initialise the dummy mesh
    call initialise_dummy_mesh( mesh, xmin, xmax, ymin, ymax)

    ! Refine the mesh
    n = 15
    do i = 1, n

      ! Smoothly increase resolution
      log_res_min = log( res_min)
      log_res_max = log( res_max)
      log_res = log_res_min + (log_res_max - log_res_min) * real(i-1,dp) / real(n-1,dp)
      resolution = exp( log_res)

      ! Define polygon
      select case (orientation)
      case default
        call crash('orientation should be x or y')
      case ('x')
        x_lo = mesh%xmin + (mesh%xmax - mesh%xmin) * real(i-1,dp) / real(n,dp)
        x_hi = mesh%xmax
        y_lo = mesh%ymin
        y_hi = mesh%ymax
      case ('y')
        x_lo = mesh%xmin
        x_hi = mesh%xmax
        y_lo = mesh%ymin + (mesh%ymax - mesh%ymin) * real(i-1,dp) / real(n,dp)
        y_hi = mesh%ymax
      end select
      poly( 1,:) = [x_lo, y_lo]
      poly( 2,:) = [x_hi, y_lo]
      poly( 3,:) = [x_hi, y_hi]
      poly( 4,:) = [x_lo, y_hi]

      ! Refine mesh over polygon
      call refine_mesh_polygon( mesh, poly, resolution, alpha_min)
      call test_mesh_is_self_consistent( mesh, ASSERTION, &
        trim(mesh_name)//': self_consistency (after refine_mesh_polygon)')

    end do

    ! Smooth the mesh
    do it_Lloyds_algorithm = 1, nit_Lloyds_algorithm
      call Lloyds_algorithm_single_iteration( mesh, alpha_min)
      call test_mesh_is_self_consistent( mesh, ASSERTION, &
        trim(mesh_name)//': self_consistency (after Lloyds_algorithm_single_iteration)')
    end do

    ! Calculate secondary geometry data (needed in order to be able to write to NetCDF)
    call calc_all_secondary_mesh_data( mesh, lambda_M, phi_M, beta_stereo)

    ! Write to NetCDF file
    filename = trim( foldername_test_meshes)//'/'//trim( mesh_name)//'.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)
    call close_netcdf_file( ncid)

    ! Add name of mesh to list file
    if (par%master) then
      open(newunit = io_filename_test_meshes_list, file = filename_test_meshes_list, &
        status = "old", action = "write", position = "append", iostat = stat, iomsg = msg)
      if (stat /= 0) then
        call crash('Could not open text file listing the filenames of the test meshes, error message "' // trim(msg) // '"')
      end if
      write(io_filename_test_meshes_list,*) trim(mesh_name)
      close(io_filename_test_meshes_list)
    end if

    if (par%master) write(0,*) '    Created test mesh ', colour_string( trim( mesh_name), 'light blue')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_test_mesh_gradient

end module component_tests_create_test_meshes
