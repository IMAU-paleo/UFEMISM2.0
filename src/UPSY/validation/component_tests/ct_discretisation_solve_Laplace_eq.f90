module ct_discretisation_solve_Laplace_eq

  ! Test the mesh matrix operators by solving a Laplace equation on the mesh

  use assertions_basic
  use tests_main
  use precisions, only: dp
  use mpi_basic, only: par, sync
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash, colour_string, warning
  use mesh_types, only: type_mesh
  use mpi_f08, only: MPI_COMM_WORLD, MPI_BCAST, MPI_CHAR, MPI_WIN
  use mpi_distributed_shared_memory, only: allocate_dist_shared, deallocate_dist_shared
  use CSR_sparse_matrix_type, only: type_sparse_matrix_CSR_dp
  use CSR_matrix_basics, only: allocate_matrix_CSR_dist, add_entry_CSR_dist, read_single_row_CSR_dist, &
  finalise_matrix_CSR_dist, deallocate_matrix_CSR_dist
  use CSR_matrix_solving, only: solve_matrix_equation_CSR_Jacobi
  use netcdf_io_main

  implicit none

  private

  public :: run_all_Laplace_eq_solving_tests

contains

  !> Run all Laplace eq solving tests.
  subroutine run_all_Laplace_eq_solving_tests( foldername_discretisation, test_mesh_filenames)

    ! In/output variables:
    character(len=*),               intent(in) :: foldername_discretisation
    character(len=*), dimension(:), intent(in) :: test_mesh_filenames

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_all_Laplace_eq_solving_tests'
    character(len=1024)            :: foldername_map_deriv
    character(len=1024)            :: test_mesh_filename
    integer                        :: i

    ! Add routine to call stack
    call init_routine( routine_name)

    if (par%primary) write(0,*) '    Running Laplace equation solving tests...'
    if (par%primary) write(0,*) ''

    call create_Laplace_eq_solving_tests_output_folder( foldername_discretisation, foldername_map_deriv)

    do i = 1, size(test_mesh_filenames)
      test_mesh_filename = trim(test_mesh_filenames( i))
      call run_Laplace_eq_solving_test_on_mesh( foldername_map_deriv, test_mesh_filename)
    end do

    if (par%primary) write(0,*) ''

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_all_Laplace_eq_solving_tests

  !> Create the output folder for the Laplace eq tests
  subroutine create_Laplace_eq_solving_tests_output_folder( foldername_discretisation, foldername_map_deriv)

    ! In/output variables:
    character(len=*), intent(in)  :: foldername_discretisation
    character(len=*), intent(out) :: foldername_map_deriv

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_Laplace_eq_tests_output_folder'
    logical                        :: ex
    integer                        :: ierr

    ! Add routine to path
    call init_routine( routine_name)

    foldername_map_deriv = trim(foldername_discretisation) // '/solve_Laplace_eq'

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

  end subroutine create_Laplace_eq_solving_tests_output_folder

  !> Run the Laplace eq solving test on a particular mesh
  subroutine run_Laplace_eq_solving_test_on_mesh( foldername_map_deriv, test_mesh_filename)

    ! In/output variables:
    character(len=*), intent(in) :: foldername_map_deriv
    character(len=*), intent(in) :: test_mesh_filename

    ! Local variables:
    character(len=1024), parameter              :: routine_name = 'run_Laplace_eq_solving_test_on_mesh'
    type(type_mesh)                             :: mesh
    integer                                     :: ncid
    real(dp), dimension(:), contiguous, pointer :: f_nih    => null()
    real(dp), dimension(:), contiguous, pointer :: f_ex_nih => null()
    type(MPI_WIN)                               :: wf_nih, wf_ex_nih
    real(dp)                                    :: c, r0, x, y
    integer                                     :: ti
    integer                                     :: ncols, ncols_loc, nrows, nrows_loc, nnz_est_proc
    type(type_sparse_matrix_CSR_dp)             :: AA
    real(dp), dimension(:), contiguous, pointer :: bb_nih => null()
    type(MPI_WIN)                               :: wbb_nih
    integer,  dimension(:), allocatable         :: single_row_ind
    real(dp), dimension(:), allocatable         :: single_row_d2dx2_val
    real(dp), dimension(:), allocatable         :: single_row_d2dy2_val
    integer                                     :: single_row_nnz
    integer                                     :: k, tj
    real(dp)                                    :: A
    integer                                     :: nit
    real(dp)                                    :: tol

    ! Add routine to call stack
    call init_routine( routine_name)

    if (par%primary) write(0,*) '      Running Laplace equation solving test on mesh ', &
      colour_string(trim(test_mesh_filename( index( test_mesh_filename,'/',back=.true.)+1:&
      len_trim( test_mesh_filename))),'light blue'), '...'

    ! Set up the mesh from the file (includes calculating secondary geometry data and matrix operators)
    call open_existing_netcdf_file_for_reading( trim(test_mesh_filename), ncid)
    call setup_mesh_from_file( test_mesh_filename, ncid, mesh)
    call close_netcdf_file( ncid)

    ! Check mesh self-consistency
    call assert( test_mesh_is_self_consistent( mesh), 'mesh is not self-consistent')

    ! ==   Solve d2f/dx2 + d2f/dy2 = c, f( r >=r0 ) = f_ex   ==
    !       Solution: f = -c/4 r0^2 + c/4 (x^2 + y^2)
    ! =========================================================

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( f_ex_nih, wf_ex_nih, mesh%pai_Tri%n_nih)
    call allocate_dist_shared( f_nih   , wf_nih   , mesh%pai_Tri%n_nih)
    call allocate_dist_shared( bb_nih  , wbb_nih  , mesh%pai_Tri%n_nih)
    f_ex_nih( mesh%pai_Tri%i1_nih:mesh%pai_Tri%i2_nih) => f_ex_nih
    f_nih   ( mesh%pai_Tri%i1_nih:mesh%pai_Tri%i2_nih) => f_nih
    bb_nih  ( mesh%pai_Tri%i1_nih:mesh%pai_Tri%i2_nih) => bb_nih

    ! Calculate exact solution
    c  = -1e-9_dp
    r0 = mesh%xmax * 0.8_dp
    do ti = mesh%ti1, mesh%ti2
      x = mesh%Trigc( ti,1)
      y = mesh%Trigc( ti,2)
      f_ex_nih( ti) = -c/4._dp * r0**2 + c/4._dp * (x**2 + y**2)
    end do

    ! Calculate stiffness matrix A

    ncols           = mesh%nTri          ! from
    ncols_loc       = mesh%nTri_loc
    nrows           = mesh%nTri          ! to
    nrows_loc       = mesh%nTri_loc
    nnz_est_proc    = mesh%M2_ddx_b_b%nnz

    call allocate_matrix_CSR_dist( AA, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc, &
      pai_x = mesh%pai_Tri, pai_y = mesh%pai_Tri)

    allocate( single_row_ind      ( mesh%nC_mem*2))
    allocate( single_row_d2dx2_val( mesh%nC_mem*2))
    allocate( single_row_d2dy2_val( mesh%nC_mem*2))

    do ti = mesh%ti1, mesh%ti2

      x = mesh%Trigc( ti,1)
      y = mesh%Trigc( ti,2)

      if (norm2([x,y]) >= r0) then
        ! f = f_ex

        call add_entry_CSR_dist( AA, ti, ti, 1._dp)

        bb_nih( ti) = f_ex_nih( ti)

      else
        ! d2f/dx2 + d2f/dy2 = c

        ! Read coefficients of the operator matrices
        call read_single_row_CSR_dist( mesh%M2_d2dx2_b_b, ti, single_row_ind, single_row_d2dx2_val, single_row_nnz)
        call read_single_row_CSR_dist( mesh%M2_d2dy2_b_b, ti, single_row_ind, single_row_d2dy2_val, single_row_nnz)

        do k = 1, single_row_nnz
          tj = single_row_ind( k)
          A = single_row_d2dx2_val( k) + single_row_d2dy2_val( k)
          call add_entry_CSR_dist( AA, ti, tj, A)
        end do

        bb_nih( ti) = c

      end if

    end do

    call finalise_matrix_CSR_dist( AA)

    ! Solve matrix equation
    nit = 5000
    tol = 1e-9_dp
    call solve_matrix_equation_CSR_Jacobi( AA, mesh%pai_Tri, f_nih, mesh%pai_Tri, bb_nih, nit, tol)

    ! Write results to NetCDF output file
    call write_Laplace_eq_solving_test_results_to_file( &
      foldername_map_deriv, test_mesh_filename, mesh, f_ex_nih, f_nih)

    ! Clean up after yourself
    call deallocate_dist_shared( f_ex_nih, wf_ex_nih)
    call deallocate_dist_shared( f_nih   , wf_nih   )
    call deallocate_dist_shared( bb_nih  , wbb_nih  )
    call deallocate_matrix_CSR_dist( AA)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_Laplace_eq_solving_test_on_mesh

  !> Write the results of Laplace eq solving test, for a particular mesh with a particular test function, to a file.
  subroutine write_Laplace_eq_solving_test_results_to_file( &
    foldername_map_deriv, test_mesh_filename, mesh, f_ex_nih, f_nih)

    ! In/output variables:
    character(len=*),       intent(in) :: foldername_map_deriv
    character(len=*),       intent(in) :: test_mesh_filename
    type(type_mesh),        intent(in) :: mesh
    real(dp), dimension(mesh%pai_Tri%i1_nih:mesh%pai_Tri%i2_nih), intent(in) :: f_ex_nih, f_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'write_Laplace_eq_solving_test_results_to_file'
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
    filename = trim( foldername_map_deriv) // '/res_' // &
      trim( mesh_name) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)

    ! ! Add all the variables
    call add_field_mesh_dp_2D_b_notime( filename, ncid, 'f_ex')
    call add_field_mesh_dp_2D_b_notime( filename, ncid, 'f'   )

    ! Write all the variables
    call write_to_field_multopt_mesh_dp_2D_b_notime( mesh, filename, ncid, 'f_ex', f_ex_nih, d_is_hybrid = .true.)
    call write_to_field_multopt_mesh_dp_2D_b_notime( mesh, filename, ncid, 'f'   , f_nih   , d_is_hybrid = .true.)

    ! Close the file
    call close_netcdf_file( ncid)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine write_Laplace_eq_solving_test_results_to_file

end module ct_discretisation_solve_Laplace_eq
