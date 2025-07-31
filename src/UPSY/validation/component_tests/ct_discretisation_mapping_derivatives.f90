module ct_discretisation_mapping_derivatives

  ! Test the mesh matrix operators for mapping and derivatives

  use mpi_f08, only: MPI_COMM_WORLD, MPI_BCAST, MPI_CHAR
  use precisions, only: dp
  use parameters
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash, colour_string, warning
  use mpi_basic, only: par, sync
  use mesh_types, only: type_mesh
  use netcdf_io_main
  use CSR_matrix_vector_multiplication, only: multiply_CSR_matrix_with_vector_1D_wrapper
  use assertions_basic
  use tests_main

  implicit none

  private

  public :: run_all_map_deriv_tests

  abstract interface
    function test_function(x,y,xmin,xmax,ymin,ymax) result(d_all)
      use precisions, only: dp
      real(dp), intent(in) :: x,y,xmin,xmax,ymin,ymax
      real(dp), dimension(6) :: d_all
    end function test_function
  end interface

contains

  !> Run all mapping/derivative tests.
  subroutine run_all_map_deriv_tests( foldername_discretisation, test_mesh_filenames)

    ! In/output variables:
    character(len=*),               intent(in) :: foldername_discretisation
    character(len=*), dimension(:), intent(in) :: test_mesh_filenames

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_all_map_deriv_tests'
    character(len=1024)            :: foldername_map_deriv
    character(len=1024)            :: test_mesh_filename
    integer                        :: i

    ! Add routine to call stack
    call init_routine( routine_name)

    if (par%primary) write(0,*) '    Running mapping/derivative tests...'
    if (par%primary) write(0,*) ''

    call create_map_deriv_tests_output_folder( foldername_discretisation, foldername_map_deriv)

    do i = 1, size(test_mesh_filenames)
      test_mesh_filename = trim(test_mesh_filenames( i))
      call run_all_map_deriv_tests_on_mesh( foldername_map_deriv, test_mesh_filename)
    end do

    if (par%primary) write(0,*) ''

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_all_map_deriv_tests

  !> Create the output folder for the mappoing/derivative tests
  subroutine create_map_deriv_tests_output_folder( foldername_discretisation, foldername_map_deriv)

    ! In/output variables:
    character(len=*), intent(in)  :: foldername_discretisation
    character(len=*), intent(out) :: foldername_map_deriv

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_map_deriv_tests_output_folder'
    logical                        :: ex
    integer                        :: ierr

    ! Add routine to path
    call init_routine( routine_name)

    foldername_map_deriv = trim(foldername_discretisation) // '/mapping_derivatives'

    if (par%primary) then

      ! Remove existing folder if necessary
      inquire( file = trim( foldername_map_deriv) // '/.', exist = ex)
      if (ex) then
        call system('rm -rf ' // trim( foldername_map_deriv))
      end if

      ! Create the directory
      call system('mkdir ' // trim( foldername_map_deriv))

    end if
    call MPI_BCAST( foldername_map_deriv, len(foldername_map_deriv), MPI_CHAR, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_map_deriv_tests_output_folder

  !> Run all the mapping/derivative tests on a particular mesh
  subroutine run_all_map_deriv_tests_on_mesh( foldername_discretisation, test_mesh_filename)

    ! In/output variables:
    character(len=*), intent(in) :: foldername_discretisation
    character(len=*), intent(in) :: test_mesh_filename

    ! Local variables:
    character(len=1024), parameter    :: routine_name = 'run_all_map_deriv_tests_on_mesh'
    type(type_mesh)                   :: mesh
    integer                           :: ncid
    procedure(test_function), pointer :: test_function_ptr
    character(len=1024)               :: function_name

    ! Add routine to call stack
    call init_routine( routine_name)

    if (par%primary) write(0,*) '      Running mapping/derivative tests on mesh ', &
      colour_string(trim(test_mesh_filename( index( test_mesh_filename,'/',back=.true.)+1:&
      len_trim( test_mesh_filename))),'light blue'), '...'

    ! Set up the mesh from the file (includes calculating secondary geometry data and matrix operators)
    call open_existing_netcdf_file_for_reading( trim(test_mesh_filename), ncid)
    call setup_mesh_from_file( test_mesh_filename, ncid, mesh)
    call close_netcdf_file( ncid)

    ! Check mesh self-consistency
    call assert( test_mesh_is_self_consistent( mesh), 'mesh is not self-consistent')

    ! Run all the mapping tests with different test functions

    test_function_ptr => test_function_linear
    function_name = 'linear'
    call run_all_map_deriv_tests_on_mesh_with_function( foldername_discretisation, &
      test_mesh_filename, mesh, test_function_ptr, function_name)

    test_function_ptr => test_function_quadratic
    function_name = 'quadratic'
    call run_all_map_deriv_tests_on_mesh_with_function( foldername_discretisation, &
      test_mesh_filename, mesh, test_function_ptr, function_name)

    test_function_ptr => test_function_periodic
    function_name = 'periodic'
    call run_all_map_deriv_tests_on_mesh_with_function( foldername_discretisation, &
      test_mesh_filename, mesh, test_function_ptr, function_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_all_map_deriv_tests_on_mesh

  !> Run all the mapping/derivative tests on a particular mesh with a particular test function
  subroutine run_all_map_deriv_tests_on_mesh_with_function( foldername_discretisation, &
    test_mesh_filename, mesh, test_function_ptr, function_name)

    ! In/output variables:
    character(len=*),    intent(in)   :: foldername_discretisation
    character(len=*),    intent(in)   :: test_mesh_filename
    type(type_mesh),     intent(in)   :: mesh
    procedure(test_function), pointer :: test_function_ptr
    character(len=1024), intent(in)   :: function_name

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'run_all_map_deriv_tests_on_mesh_with_function'
    real(dp), dimension(:), allocatable :: d_a_ex, ddx_a_ex, ddy_a_ex
    real(dp), dimension(:), allocatable :: d_b_ex, ddx_b_ex, ddy_b_ex, d2dx2_b_ex, d2dxdy_b_ex, d2dy2_b_ex
    real(dp), dimension(:), allocatable :: d_b_a, d_a_b
    real(dp), dimension(:), allocatable :: ddx_a_a, ddx_a_b
    real(dp), dimension(:), allocatable :: ddx_b_a, ddx_b_b
    real(dp), dimension(:), allocatable :: ddy_a_a, ddy_a_b
    real(dp), dimension(:), allocatable :: ddy_b_a, ddy_b_b
    real(dp), dimension(:), allocatable :: ddx_b_b_2nd, ddy_b_b_2nd, d2dx2_b_b_2nd, d2dxdy_b_b_2nd, d2dy2_b_b_2nd

    ! Add routine to call stack
    call init_routine( routine_name)

    if (par%primary) write(0,*) '        Running all mapping/derivative tests on function ', &
      colour_string(trim(function_name),'light blue'), '...'

    ! Calculate exact solutions
    call run_all_map_deriv_tests_on_mesh_with_function_calc_ex( mesh, test_function_ptr, &
      d_a_ex, ddx_a_ex, ddy_a_ex, &
      d_b_ex, ddx_b_ex, ddy_b_ex, d2dx2_b_ex, d2dxdy_b_ex, d2dy2_b_ex)

    ! Calculate discretised approximations
    call run_all_map_deriv_tests_on_mesh_with_function_calc_disc( mesh, &
      d_a_ex, d_b_ex, &
      d_a_b, d_b_a, &
      ddx_a_a, ddx_a_b, &
      ddx_b_a, ddx_b_b, &
      ddy_a_a, ddy_a_b, &
      ddy_b_a, ddy_b_b, &
      ddx_b_b_2nd, ddy_b_b_2nd, d2dx2_b_b_2nd, d2dxdy_b_b_2nd, d2dy2_b_b_2nd)

    ! Create an output file for the mapping/derivative tests, for a particular mesh with a particular test function
    call  write_map_deriv_test_results_to_file( &
      foldername_discretisation, test_mesh_filename, mesh, function_name, &
      d_a_ex, ddx_a_ex, ddy_a_ex, &
      d_b_ex, ddx_b_ex, ddy_b_ex, d2dx2_b_ex, d2dxdy_b_ex, d2dy2_b_ex, &
      d_a_b, d_b_a, &
      ddx_a_a, ddx_a_b, &
      ddx_b_a, ddx_b_b, &
      ddy_a_a, ddy_a_b, &
      ddy_b_a, ddy_b_b, &
      ddx_b_b_2nd, ddy_b_b_2nd, d2dx2_b_b_2nd, d2dxdy_b_b_2nd, d2dy2_b_b_2nd)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_all_map_deriv_tests_on_mesh_with_function

  !> Calculate the exact solutions for the mapping/derivative tests
  subroutine run_all_map_deriv_tests_on_mesh_with_function_calc_ex( mesh, test_function_ptr, &
    d_a_ex, ddx_a_ex, ddy_a_ex, &
    d_b_ex, ddx_b_ex, ddy_b_ex, d2dx2_b_ex, d2dxdy_b_ex, d2dy2_b_ex)

    ! In/output variables:
    type(type_mesh),                     intent(in)  :: mesh
    procedure(test_function), pointer                :: test_function_ptr
    real(dp), dimension(:), allocatable, intent(out) :: d_a_ex, ddx_a_ex, ddy_a_ex
    real(dp), dimension(:), allocatable, intent(out) :: d_b_ex, ddx_b_ex, ddy_b_ex, d2dx2_b_ex, d2dxdy_b_ex, d2dy2_b_ex

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'run_all_map_deriv_tests_on_mesh_with_function_calc_ex'
    integer                             :: vi,ti,row
    real(dp)                            :: x,y
    real(dp), dimension(6)              :: d_all
    real(dp)                            :: d,ddx,ddy,d2dx2,d2dxdy,d2dy2

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Allocate memory
    allocate( d_a_ex     ( mesh%vi1:mesh%vi2))
    allocate( ddx_a_ex   ( mesh%vi1:mesh%vi2))
    allocate( ddy_a_ex   ( mesh%vi1:mesh%vi2))

    allocate( d_b_ex     ( mesh%ti1:mesh%ti2))
    allocate( ddx_b_ex   ( mesh%ti1:mesh%ti2))
    allocate( ddy_b_ex   ( mesh%ti1:mesh%ti2))
    allocate( d2dx2_b_ex ( mesh%ti1:mesh%ti2))
    allocate( d2dxdy_b_ex( mesh%ti1:mesh%ti2))
    allocate( d2dy2_b_ex ( mesh%ti1:mesh%ti2))

    ! Calculate exact solutions

    ! a-grid (vertices)
    do vi = mesh%vi1, mesh%vi2
      x = mesh%V( vi,1)
      y = mesh%V( vi,2)
      d_all = test_function_ptr( x, y, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax)
      d      = d_all( 1)
      ddx    = d_all( 2)
      ddy    = d_all( 3)
      row = mesh%vi2n( vi)
      d_a_ex(   row) = d
      ddx_a_ex( row) = ddx
      ddy_a_ex( row) = ddy
    end do

    ! b-grid (triangles)
    do ti = mesh%ti1, mesh%ti2
      x = mesh%Trigc( ti,1)
      y = mesh%Trigc( ti,2)
      d_all = test_function_ptr( x, y, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax)
      d      = d_all( 1)
      ddx    = d_all( 2)
      ddy    = d_all( 3)
      d2dx2  = d_all( 4)
      d2dxdy = d_all( 5)
      d2dy2  = d_all( 6)
      row = mesh%ti2n( ti)
      d_b_ex     ( row) = d
      ddx_b_ex   ( row) = ddx
      ddy_b_ex   ( row) = ddy
      d2dx2_b_ex ( row) = d2dx2
      d2dxdy_b_ex( row) = d2dxdy
      d2dy2_b_ex ( row) = d2dy2
    end do

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_all_map_deriv_tests_on_mesh_with_function_calc_ex

  !> Calculate the discretised approximations for the mapping/derivative tests
  subroutine run_all_map_deriv_tests_on_mesh_with_function_calc_disc( mesh, &
    d_a_ex, d_b_ex, d_a_b, d_b_a, &
    ddx_a_a, ddx_a_b, &
    ddx_b_a, ddx_b_b, &
    ddy_a_a, ddy_a_b, &
    ddy_b_a, ddy_b_b, &
    ddx_b_b_2nd, ddy_b_b_2nd, d2dx2_b_b_2nd, d2dxdy_b_b_2nd, d2dy2_b_b_2nd)

    ! In/output variables:
    type(type_mesh),                     intent(in)  :: mesh
    real(dp), dimension(:),              intent(in)  :: d_a_ex, d_b_ex
    real(dp), dimension(:), allocatable, intent(out) :: d_a_b, d_b_a
    real(dp), dimension(:), allocatable, intent(out) :: ddx_a_a, ddx_a_b
    real(dp), dimension(:), allocatable, intent(out) :: ddx_b_a, ddx_b_b
    real(dp), dimension(:), allocatable, intent(out) :: ddy_a_a, ddy_a_b
    real(dp), dimension(:), allocatable, intent(out) :: ddy_b_a, ddy_b_b
    real(dp), dimension(:), allocatable, intent(out) :: ddx_b_b_2nd, ddy_b_b_2nd, d2dx2_b_b_2nd, d2dxdy_b_b_2nd, d2dy2_b_b_2nd

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'run_all_map_deriv_tests_on_mesh_with_function_calc_disc'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Allocate memory

    allocate( d_a_b(       mesh%ti1:mesh%ti2))
    allocate( d_b_a(       mesh%vi1:mesh%vi2))

    allocate( ddx_a_a(     mesh%vi1:mesh%vi2))
    allocate( ddx_a_b(     mesh%ti1:mesh%ti2))

    allocate( ddx_b_a(     mesh%vi1:mesh%vi2))
    allocate( ddx_b_b(     mesh%ti1:mesh%ti2))

    allocate( ddy_a_a(     mesh%vi1:mesh%vi2))
    allocate( ddy_a_b(     mesh%ti1:mesh%ti2))

    allocate( ddy_b_a(     mesh%vi1:mesh%vi2))
    allocate( ddy_b_b(     mesh%ti1:mesh%ti2))

    allocate( ddx_b_b_2nd(     mesh%ti1:mesh%ti2))
    allocate( ddy_b_b_2nd(     mesh%ti1:mesh%ti2))
    allocate( d2dx2_b_b_2nd(   mesh%ti1:mesh%ti2))
    allocate( d2dxdy_b_b_2nd(  mesh%ti1:mesh%ti2))
    allocate( d2dy2_b_b_2nd(   mesh%ti1:mesh%ti2))

    ! Calculate discretised approximations
    call multiply_CSR_matrix_with_vector_1D_wrapper( mesh%M_map_a_b, &
      mesh%pai_V, d_a_ex, mesh%pai_Tri, d_a_b, &
      xx_is_hybrid = .false., yy_is_hybrid = .false., &
      buffer_xx_nih = mesh%buffer1_d_a_nih, buffer_yy_nih = mesh%buffer2_d_b_nih)

    call multiply_CSR_matrix_with_vector_1D_wrapper( mesh%M_map_b_a, &
      mesh%pai_Tri, d_b_ex, mesh%pai_V, d_b_a, &
      xx_is_hybrid = .false., yy_is_hybrid = .false., &
      buffer_xx_nih = mesh%buffer1_d_b_nih, buffer_yy_nih = mesh%buffer2_d_a_nih)



    call multiply_CSR_matrix_with_vector_1D_wrapper( mesh%M_ddx_a_a, &
      mesh%pai_V, d_a_ex, mesh%pai_V, ddx_a_a, &
      xx_is_hybrid = .false., yy_is_hybrid = .false., &
      buffer_xx_nih = mesh%buffer1_d_a_nih, buffer_yy_nih = mesh%buffer2_d_a_nih)
    call multiply_CSR_matrix_with_vector_1D_wrapper( mesh%M_ddx_a_b, &
      mesh%pai_V, d_a_ex, mesh%pai_Tri, ddx_a_b, &
      xx_is_hybrid = .false., yy_is_hybrid = .false., &
      buffer_xx_nih = mesh%buffer1_d_a_nih, buffer_yy_nih = mesh%buffer2_d_b_nih)

    call multiply_CSR_matrix_with_vector_1D_wrapper( mesh%M_ddx_b_a, &
      mesh%pai_Tri, d_b_ex, mesh%pai_V, ddx_b_a, &
      xx_is_hybrid = .false., yy_is_hybrid = .false., &
      buffer_xx_nih = mesh%buffer1_d_b_nih, buffer_yy_nih = mesh%buffer2_d_a_nih)
    call multiply_CSR_matrix_with_vector_1D_wrapper( mesh%M_ddx_b_b, &
      mesh%pai_Tri, d_b_ex, mesh%pai_Tri, ddx_b_b, &
      xx_is_hybrid = .false., yy_is_hybrid = .false., &
      buffer_xx_nih = mesh%buffer1_d_b_nih, buffer_yy_nih = mesh%buffer2_d_b_nih)



    call multiply_CSR_matrix_with_vector_1D_wrapper( mesh%M_ddy_a_a, &
      mesh%pai_V, d_a_ex, mesh%pai_V, ddy_a_a, &
      xx_is_hybrid = .false., yy_is_hybrid = .false., &
      buffer_xx_nih = mesh%buffer1_d_a_nih, buffer_yy_nih = mesh%buffer2_d_a_nih)
    call multiply_CSR_matrix_with_vector_1D_wrapper( mesh%M_ddy_a_b, &
      mesh%pai_V, d_a_ex, mesh%pai_Tri, ddy_a_b, &
      xx_is_hybrid = .false., yy_is_hybrid = .false., &
      buffer_xx_nih = mesh%buffer1_d_a_nih, buffer_yy_nih = mesh%buffer2_d_b_nih)

    call multiply_CSR_matrix_with_vector_1D_wrapper( mesh%M_ddy_b_a, &
      mesh%pai_Tri, d_b_ex, mesh%pai_V, ddy_b_a, &
      xx_is_hybrid = .false., yy_is_hybrid = .false., &
      buffer_xx_nih = mesh%buffer1_d_b_nih, buffer_yy_nih = mesh%buffer2_d_a_nih)
    call multiply_CSR_matrix_with_vector_1D_wrapper( mesh%M_ddy_b_b, &
      mesh%pai_Tri, d_b_ex, mesh%pai_Tri, ddy_b_b, &
      xx_is_hybrid = .false., yy_is_hybrid = .false., &
      buffer_xx_nih = mesh%buffer1_d_b_nih, buffer_yy_nih = mesh%buffer2_d_b_nih)



    call multiply_CSR_matrix_with_vector_1D_wrapper( mesh%M2_ddx_b_b, &
      mesh%pai_Tri, d_b_ex, mesh%pai_Tri, ddx_b_b_2nd, &
      xx_is_hybrid = .false., yy_is_hybrid = .false., &
      buffer_xx_nih = mesh%buffer1_d_b_nih, buffer_yy_nih = mesh%buffer2_d_b_nih)
    call multiply_CSR_matrix_with_vector_1D_wrapper( mesh%M2_ddy_b_b, &
      mesh%pai_Tri, d_b_ex, mesh%pai_Tri, ddy_b_b_2nd, &
      xx_is_hybrid = .false., yy_is_hybrid = .false., &
      buffer_xx_nih = mesh%buffer1_d_b_nih, buffer_yy_nih = mesh%buffer2_d_b_nih)
    call multiply_CSR_matrix_with_vector_1D_wrapper( mesh%M2_d2dx2_b_b, &
      mesh%pai_Tri, d_b_ex, mesh%pai_Tri, d2dx2_b_b_2nd, &
      xx_is_hybrid = .false., yy_is_hybrid = .false., &
      buffer_xx_nih = mesh%buffer1_d_b_nih, buffer_yy_nih = mesh%buffer2_d_b_nih)
    call multiply_CSR_matrix_with_vector_1D_wrapper( mesh%M2_d2dxdy_b_b, &
      mesh%pai_Tri, d_b_ex, mesh%pai_Tri, d2dxdy_b_b_2nd, &
      xx_is_hybrid = .false., yy_is_hybrid = .false., &
      buffer_xx_nih = mesh%buffer1_d_b_nih, buffer_yy_nih = mesh%buffer2_d_b_nih)
    call multiply_CSR_matrix_with_vector_1D_wrapper( mesh%M2_d2dy2_b_b, &
      mesh%pai_Tri, d_b_ex, mesh%pai_Tri, d2dy2_b_b_2nd, &
      xx_is_hybrid = .false., yy_is_hybrid = .false., &
      buffer_xx_nih = mesh%buffer1_d_b_nih, buffer_yy_nih = mesh%buffer2_d_b_nih)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_all_map_deriv_tests_on_mesh_with_function_calc_disc

  !> Write the results of the mapping/derivative tests, for a particular mesh with a particular test function, to a file.
  subroutine write_map_deriv_test_results_to_file( &
    foldername_discretisation, test_mesh_filename, mesh, function_name, &
    d_a_ex, ddx_a_ex, ddy_a_ex, &
    d_b_ex, ddx_b_ex, ddy_b_ex, d2dx2_b_ex, d2dxdy_b_ex, d2dy2_b_ex, &
    d_a_b, d_b_a, &
    ddx_a_a, ddx_a_b, &
    ddx_b_a, ddx_b_b, &
    ddy_a_a, ddy_a_b, &
    ddy_b_a, ddy_b_b, &
    ddx_b_b_2nd, ddy_b_b_2nd, d2dx2_b_b_2nd, d2dxdy_b_b_2nd, d2dy2_b_b_2nd)

    ! In/output variables:
    character(len=*),       intent(in) :: foldername_discretisation
    character(len=*),       intent(in) :: test_mesh_filename
    type(type_mesh),        intent(in) :: mesh
    character(len=1024),    intent(in) :: function_name
    real(dp), dimension(:), intent(in) :: d_a_ex, ddx_a_ex, ddy_a_ex
    real(dp), dimension(:), intent(in) :: d_b_ex, ddx_b_ex, ddy_b_ex, d2dx2_b_ex, d2dxdy_b_ex, d2dy2_b_ex
    real(dp), dimension(:), intent(in) :: d_a_b, d_b_a
    real(dp), dimension(:), intent(in) :: ddx_a_a, ddx_a_b
    real(dp), dimension(:), intent(in) :: ddx_b_a, ddx_b_b
    real(dp), dimension(:), intent(in) :: ddy_a_a, ddy_a_b
    real(dp), dimension(:), intent(in) :: ddy_b_a, ddy_b_b
    real(dp), dimension(:), intent(in) :: ddx_b_b_2nd, ddy_b_b_2nd, d2dx2_b_b_2nd, d2dxdy_b_b_2nd, d2dy2_b_b_2nd

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'write_map_deriv_test_results_to_file'
    integer                        :: ncid
    character(len=1024)            :: mesh_name
    character(len=1024)            :: filename
    integer                        :: i

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Create a file and write the mesh to it
    mesh_name = test_mesh_filename( 1:len_trim( test_mesh_filename)-3)
    i = index( mesh_name, '/', back = .true.)
    mesh_name = mesh_name( i+1:len_trim( mesh_name))
    filename = trim(foldername_discretisation) // '/res_' // &
      trim( mesh_name) // '_' // trim( function_name) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)

    ! Add all the variables
    call add_field_mesh_dp_2D_notime(   filename, ncid, 'd_a_ex'     )
    call add_field_mesh_dp_2D_notime(   filename, ncid, 'ddx_a_ex'   )
    call add_field_mesh_dp_2D_notime(   filename, ncid, 'ddy_a_ex'   )

    call add_field_mesh_dp_2D_b_notime( filename, ncid, 'd_b_ex'     )
    call add_field_mesh_dp_2D_b_notime( filename, ncid, 'ddx_b_ex'   )
    call add_field_mesh_dp_2D_b_notime( filename, ncid, 'ddy_b_ex'   )
    call add_field_mesh_dp_2D_b_notime( filename, ncid, 'd2dx2_b_ex' )
    call add_field_mesh_dp_2D_b_notime( filename, ncid, 'd2dxdy_b_ex')
    call add_field_mesh_dp_2D_b_notime( filename, ncid, 'd2dy2_b_ex' )

    call add_field_mesh_dp_2D_b_notime( filename, ncid, 'd_a_b'      )
    call add_field_mesh_dp_2D_notime(   filename, ncid, 'd_b_a'      )

    call add_field_mesh_dp_2D_notime(   filename, ncid, 'ddx_a_a'    )
    call add_field_mesh_dp_2D_b_notime( filename, ncid, 'ddx_a_b'    )

    call add_field_mesh_dp_2D_notime(   filename, ncid, 'ddx_b_a'    )
    call add_field_mesh_dp_2D_b_notime( filename, ncid, 'ddx_b_b'    )

    call add_field_mesh_dp_2D_notime(   filename, ncid, 'ddy_a_a'    )
    call add_field_mesh_dp_2D_b_notime( filename, ncid, 'ddy_a_b'    )

    call add_field_mesh_dp_2D_notime(   filename, ncid, 'ddy_b_a'    )
    call add_field_mesh_dp_2D_b_notime( filename, ncid, 'ddy_b_b'    )

    call add_field_mesh_dp_2D_b_notime( filename, ncid, 'ddx_b_b_2nd'   )
    call add_field_mesh_dp_2D_b_notime( filename, ncid, 'ddy_b_b_2nd'   )
    call add_field_mesh_dp_2D_b_notime( filename, ncid, 'd2dx2_b_b_2nd' )
    call add_field_mesh_dp_2D_b_notime( filename, ncid, 'd2dxdy_b_b_2nd')
    call add_field_mesh_dp_2D_b_notime( filename, ncid, 'd2dy2_b_b_2nd' )

    ! Write all the variables
    call write_to_field_multopt_mesh_dp_2D_notime(   mesh, filename, ncid, 'd_a_ex'     , d_a_ex     )
    call write_to_field_multopt_mesh_dp_2D_notime(   mesh, filename, ncid, 'ddx_a_ex'   , ddx_a_ex   )
    call write_to_field_multopt_mesh_dp_2D_notime(   mesh, filename, ncid, 'ddy_a_ex'   , ddy_a_ex   )

    call write_to_field_multopt_mesh_dp_2D_b_notime( mesh, filename, ncid, 'd_b_ex'     , d_b_ex     )
    call write_to_field_multopt_mesh_dp_2D_b_notime( mesh, filename, ncid, 'ddx_b_ex'   , ddx_b_ex   )
    call write_to_field_multopt_mesh_dp_2D_b_notime( mesh, filename, ncid, 'ddy_b_ex'   , ddy_b_ex   )
    call write_to_field_multopt_mesh_dp_2D_b_notime( mesh, filename, ncid, 'd2dx2_b_ex' , d2dx2_b_ex )
    call write_to_field_multopt_mesh_dp_2D_b_notime( mesh, filename, ncid, 'd2dxdy_b_ex', d2dxdy_b_ex)
    call write_to_field_multopt_mesh_dp_2D_b_notime( mesh, filename, ncid, 'd2dy2_b_ex' , d2dy2_b_ex )

    call write_to_field_multopt_mesh_dp_2D_b_notime( mesh, filename, ncid, 'd_a_b'      , d_a_b      )
    call write_to_field_multopt_mesh_dp_2D_notime(   mesh, filename, ncid, 'd_b_a'      , d_b_a      )

    call write_to_field_multopt_mesh_dp_2D_notime(   mesh, filename, ncid, 'ddx_a_a'    , ddx_a_a    )
    call write_to_field_multopt_mesh_dp_2D_b_notime( mesh, filename, ncid, 'ddx_a_b'    , ddx_a_b    )

    call write_to_field_multopt_mesh_dp_2D_notime(   mesh, filename, ncid, 'ddx_b_a'    , ddx_b_a    )
    call write_to_field_multopt_mesh_dp_2D_b_notime( mesh, filename, ncid, 'ddx_b_b'    , ddx_b_b    )

    call write_to_field_multopt_mesh_dp_2D_notime(   mesh, filename, ncid, 'ddy_a_a'    , ddy_a_a    )
    call write_to_field_multopt_mesh_dp_2D_b_notime( mesh, filename, ncid, 'ddy_a_b'    , ddy_a_b    )

    call write_to_field_multopt_mesh_dp_2D_notime(   mesh, filename, ncid, 'ddy_b_a'    , ddy_b_a    )
    call write_to_field_multopt_mesh_dp_2D_b_notime( mesh, filename, ncid, 'ddy_b_b'    , ddy_b_b    )

    call write_to_field_multopt_mesh_dp_2D_b_notime( mesh, filename, ncid, 'ddx_b_b_2nd'    , ddx_b_b_2nd    )
    call write_to_field_multopt_mesh_dp_2D_b_notime( mesh, filename, ncid, 'ddy_b_b_2nd'    , ddy_b_b_2nd    )
    call write_to_field_multopt_mesh_dp_2D_b_notime( mesh, filename, ncid, 'd2dx2_b_b_2nd'  , d2dx2_b_b_2nd  )
    call write_to_field_multopt_mesh_dp_2D_b_notime( mesh, filename, ncid, 'd2dxdy_b_b_2nd' , d2dxdy_b_b_2nd )
    call write_to_field_multopt_mesh_dp_2D_b_notime( mesh, filename, ncid, 'd2dy2_b_b_2nd'  , d2dy2_b_b_2nd  )

    ! Close the file
    call close_netcdf_file( ncid)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine write_map_deriv_test_results_to_file

  !> A simple linear test function for the mapping/derivative tests
  function test_function_linear(x,y,xmin,xmax,ymin,ymax) result(d_all)

    ! In/output variables:
    real(dp), intent(in)   :: x,y,xmin,xmax,ymin,ymax
    real(dp), dimension(6) :: d_all

    ! Local variables:
    real(dp) :: c1,c2

    c1 = 2._dp / (xmax - xmin)
    c2 = 3._dp / (ymax - ymin)

    d_all( 1) = c1*x + c2*y ! f
    d_all( 2) = c1          ! df/dx
    d_all( 3) = c2          ! df/dy
    d_all( 4) = 0._dp       ! d2f/dx2
    d_all( 5) = 0._dp       ! d2f/dxdy
    d_all( 6) = 0._dp       ! d2f/dy2

  end function test_function_linear

  !> A simple quadratic test function for the mapping/derivative tests
  function test_function_quadratic(x,y,xmin,xmax,ymin,ymax) result(d_all)

    ! In/output variables:
    real(dp), intent(in)   :: x,y,xmin,xmax,ymin,ymax
    real(dp), dimension(6) :: d_all

    ! Local variables:
    real(dp) :: c1,c2,c3

    c1 = 2._dp / (xmax - xmin)
    c2 = 3._dp / (ymax - ymin)
    c3 = 5._dp / (ymax - ymin)

    d_all( 1) = (c1*x)**2 + (c2*y)**2 + c3*x*y ! f
    d_all( 2) = 2._dp * c1**2 * x  + c3*y      ! df/dx
    d_all( 3) = 2._dp * c2**2 * y  + c3*x      ! df/dy
    d_all( 4) = 2._dp * c1**2                  ! d2f/dx2
    d_all( 5) = c3                             ! d2f/dxdy
    d_all( 6) = 2._dp * c2**2                  ! d2f/dy2

  end function test_function_quadratic

  !> A simple periodic test function for the mapping/derivative tests
  function test_function_periodic(x,y,xmin,xmax,ymin,ymax) result(d_all)

    ! In/output variables:
    real(dp), intent(in)   :: x,y,xmin,xmax,ymin,ymax
    real(dp), dimension(6) :: d_all

    ! Local variables:
    real(dp) :: c1,c2

    c1 = 2._dp * pi / (xmax - xmin)
    c2 = 3._dp * pi / (ymax - ymin)

    d_all( 1) =            sin( c1 * (x - xmin)) *            sin( c2 * (y - ymin)) ! f
    d_all( 2) =   c1     * cos( c1 * (x - xmin)) *            sin( c2 * (y - ymin)) ! df/dx
    d_all( 3) =            sin( c1 * (x - xmin)) *   c2     * cos( c2 * (y - ymin)) ! df/dy
    d_all( 4) = -(c1**2) * sin( c1 * (x - xmin)) *            sin( c2 * (y - ymin)) ! d2f/dx2
    d_all( 5) =   c1     * cos( c1 * (x - xmin)) *   c2     * cos( c2 * (y - ymin)) ! d2f/dxdy
    d_all( 6) =            sin( c1 * (x - xmin)) * -(c2**2) * sin( c2 * (y - ymin)) ! d2f/dy2

  end function test_function_periodic

end module ct_discretisation_mapping_derivatives
