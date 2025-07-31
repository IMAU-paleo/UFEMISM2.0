module netcdf_save_single_variables
  !< Contains routines for quickly writing variables to NetCDF to check their content.

#include <petsc/finclude/petscksp.h>
  use petscksp
  use precisions, only: dp
  use mpi_basic, only: par
  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use CSR_sparse_matrix_type, only: type_sparse_matrix_CSR_dp
  use CSR_matrix_basics, only: deallocate_matrix_CSR_dist, gather_CSR_dist_to_primary
  use petsc_basic, only: mat_petsc2CSR
  use mpi_distributed_memory, only: gather_to_primary
  use netcdf, only: NF90_INT, NF90_DOUBLE
  use netcdf_inquire_grid_mesh
  use netcdf_write_var_primary
  use netcdf_basic_wrappers

  implicit none

  private

  public :: write_CSR_matrix_to_NetCDF, write_PETSc_matrix_to_NetCDF, &
    save_variable_as_netcdf_logical_1D, save_variable_as_netcdf_int_1D, save_variable_as_netcdf_int_2D, &
    save_variable_as_netcdf_dp_1D, save_variable_as_netcdf_dp_2D

contains

  subroutine write_PETSc_matrix_to_NetCDF( output_dir, A, filename)
    !< Write a PETSc matrix to a NetCDF file

    ! In- and output variables:
    character(len=*), intent(in) :: output_dir
    type(tMat),       intent(in) :: A
    character(len=*), intent(in) :: filename

    ! Local variables:
    character(len=256), parameter   :: routine_name = 'write_PETSc_matrix_to_NetCDF'
    type(type_sparse_matrix_CSR_dp) :: AA

    ! Add routine to path
    call init_routine( routine_name)

    ! Get matrix in CSR format using native Fortran arrays
    call mat_petsc2CSR( A, AA)

    ! Write the CSR matrix to a file
    call write_CSR_matrix_to_NetCDF( output_dir, AA, filename)

    ! Clean up after yourself
    call deallocate_matrix_CSR_dist( AA)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_PETSc_matrix_to_NetCDF

  subroutine write_CSR_matrix_to_NetCDF( output_dir, AA, filename)
    !< Write a CSR matrix to a NetCDF file

    ! In- and output variables:
    character(len=*),                intent(in) :: output_dir
    type(type_sparse_matrix_CSR_dp), intent(in) :: AA
    character(len=*),                intent(in) :: filename

    ! Local variables:
    character(len=256), parameter   :: routine_name = 'write_CSR_matrix_to_NetCDF'
    character(len=256)              :: filename_applied
    integer                         :: ncid
    integer                         :: id_dim_m, id_dim_mp1, id_dim_n, id_dim_nnz
    integer                         :: id_var_ptr, id_var_ind, id_var_val
    type(type_sparse_matrix_CSR_dp) :: AA_tot

    ! Add routine to path
    call init_routine( routine_name)

    ! Gather distributed matrix to the primary
    call gather_CSR_dist_to_primary( AA, AA_tot)

    ! Append output directory to filename
    filename_applied = trim( output_dir) // '/' // trim( filename) // '.nc'
    call delete_existing_file( filename_applied)

    ! Create a new NetCDF file
    call create_new_netcdf_file_for_writing( filename_applied, ncid)

    ! Create dimensions
    call create_dimension( filename_applied, ncid, 'm'     , AA_tot%m  , id_dim_m  )
    call create_dimension( filename_applied, ncid, 'mplus1', AA_tot%m+1, id_dim_mp1)
    call create_dimension( filename_applied, ncid, 'n'     , AA_tot%n  , id_dim_n  )
    call create_dimension( filename_applied, ncid, 'nnz'   , AA_tot%nnz, id_dim_nnz)

    ! Create variables
    call create_variable( filename_applied, ncid, 'ptr', NF90_INT   , [id_dim_mp1], id_var_ptr)
    call create_variable( filename_applied, ncid, 'ind', NF90_INT   , [id_dim_nnz], id_var_ind)
    call create_variable( filename_applied, ncid, 'val', NF90_DOUBLE, [id_dim_nnz], id_var_val)

    ! Write to NetCDF
    call write_var_primary( filename_applied, ncid, id_var_ptr, AA_tot%ptr               )
    call write_var_primary( filename_applied, ncid, id_var_ind, AA_tot%ind( 1:AA_tot%nnz))
    call write_var_primary(  filename_applied, ncid, id_var_val, AA_tot%val( 1:AA_tot%nnz))

    ! Close the NetCDF file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_CSR_matrix_to_NetCDF

  subroutine save_variable_as_netcdf_logical_1D( output_dir, d_partial, field_name)
    !< Save a single variable to a NetCDF file

    ! Input variables:
    character(len=*),      intent(in) :: output_dir
    logical, dimension(:), intent(in) :: d_partial
    character(len=*),      intent(in) :: field_name

    ! Local variables:
    character(len=256), parameter       :: routine_name = 'save_variable_as_netcdf_logical_1D'
    integer,  dimension(:), allocatable :: d_partial_int

    ! Add routine to path
    call init_routine( routine_name)

    ! allocate memory for integer version
    allocate( d_partial_int( SIZE( d_partial,1)))

    ! Fill integer version
    where (d_partial)
      d_partial_int = 1
    elsewhere
      d_partial_int = 0
    end where

    ! Write integer version to NetCDF
    call save_variable_as_netcdf_int_1D( output_dir, d_partial_int, field_name)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine save_variable_as_netcdf_logical_1D

  subroutine save_variable_as_netcdf_int_1D( output_dir, d_partial, field_name)
    !< Save a single variable to a NetCDF file

    ! Input variables:
    character(len=*),      intent(in) :: output_dir
    integer, dimension(:), intent(in) :: d_partial
    character(len=*),      intent(in) :: field_name

    ! Local variables:
    character(len=256), parameter      :: routine_name = 'save_variable_as_netcdf_int_1D'
    character(len=256)                 :: filename
    integer                            :: ncid
    integer                            :: n_partial, n_tot
    integer                            :: ierr
    integer, dimension(:), allocatable :: d_tot
    integer                            :: id_dim_n1
    integer                            :: id_var

    ! Add routine to path
    call init_routine( routine_name)

    ! Determine file name
    filename = trim( output_dir) // '/' // trim( field_name) // '.nc'
    call delete_existing_file( filename)

    ! Create a new NetCDF file
    call create_new_netcdf_file_for_writing( filename, ncid)

    ! Gather data to the primary
    n_partial = size( d_partial,1)
    call MPI_ALLREDUCE( n_partial, n_tot, 1, MPI_integer, MPI_SUM, MPI_COMM_WORLD, ierr)
    if (par%primary) then
      allocate( d_tot( n_tot))
      call gather_to_primary( d_partial, d_tot)
    else
      call gather_to_primary( d_partial)
    end if

    ! Create dimensions
    call create_dimension( filename, ncid, 'n1', n_tot, id_dim_n1)

    ! Create variable
    call create_variable( filename, ncid, field_name, NF90_INT, (/ id_dim_n1 /), id_var)

    ! Write data
    call write_var_primary( filename, ncid, id_var, d_tot)

    ! Close the NetCDF file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine save_variable_as_netcdf_int_1D

  subroutine save_variable_as_netcdf_int_2D( output_dir, d_partial, field_name)
    !< Save a single variable to a NetCDF file

    ! Input variables:
    character(len=*),        intent(in) :: output_dir
    integer, dimension(:,:), intent(in) :: d_partial
    character(len=*),        intent(in) :: field_name

    ! Local variables:
    character(len=256), parameter        :: routine_name = 'save_variable_as_netcdf_int_2D'
    character(len=256)                   :: filename
    integer                              :: ncid
    integer                              :: n_partial, n_tot, n2
    integer                              :: ierr
    integer, dimension(:,:), allocatable :: d_tot
    integer                              :: id_dim_n1, id_dim_n2
    integer                              :: id_var

    ! Add routine to path
    call init_routine( routine_name)

    ! Determine file name
    filename = trim( output_dir) // '/' // trim( field_name) // '.nc'
    call delete_existing_file( filename)

    ! Create a new NetCDF file
    call create_new_netcdf_file_for_writing( filename, ncid)

    ! Gather data to the primary
    n_partial = size( d_partial,1)
    n2        = size( d_partial,2)
    call MPI_ALLREDUCE( n_partial, n_tot, 1, MPI_integer, MPI_SUM, MPI_COMM_WORLD, ierr)
    if (par%primary) then
      allocate( d_tot( n_tot, n2))
      call gather_to_primary( d_partial, d_tot)
    else
      call gather_to_primary( d_partial)
    end if

    ! Create dimensions
    call create_dimension( filename, ncid, 'n1', n_tot, id_dim_n1)
    call create_dimension( filename, ncid, 'n2', n_tot, id_dim_n2)

    ! Create variable
    call create_variable( filename, ncid, field_name, NF90_INT, (/ id_dim_n1, id_dim_n2 /), id_var)

    ! Write data
    call write_var_primary( filename, ncid, id_var, d_tot)

    ! Close the NetCDF file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine save_variable_as_netcdf_int_2D

  subroutine save_variable_as_netcdf_dp_1D( output_dir, d_partial, field_name)
    !< Save a single variable to a NetCDF file

    ! Input variables:
    character(len=*),       intent(in) :: output_dir
    real(dp), dimension(:), intent(in) :: d_partial
    character(len=*),       intent(in) :: field_name

    ! Local variables:
    character(len=256), parameter       :: routine_name = 'save_variable_as_netcdf_dp_1D'
    character(len=256)                  :: filename
    integer                             :: ncid
    integer                             :: n_partial, n_tot
    integer                             :: ierr
    real(dp), dimension(:), allocatable :: d_tot
    integer                             :: id_dim_n1
    integer                             :: id_var

    ! Add routine to path
    call init_routine( routine_name)

    ! Determine file name
    filename = trim( output_dir) // '/' // trim( field_name) // '.nc'
    call delete_existing_file( filename)

    ! Create a new NetCDF file
    call create_new_netcdf_file_for_writing( filename, ncid)

    ! Gather data to the primary
    n_partial = size( d_partial,1)
    call MPI_ALLREDUCE( n_partial, n_tot, 1, MPI_integer, MPI_SUM, MPI_COMM_WORLD, ierr)
    if (par%primary) then
      allocate( d_tot( n_tot))
      call gather_to_primary( d_partial, d_tot)
    else
      call gather_to_primary( d_partial)
    end if

    ! Create dimensions
    call create_dimension( filename, ncid, 'n1', n_tot, id_dim_n1)

    ! Create variable
    call create_variable( filename, ncid, field_name, NF90_DOUBLE, (/ id_dim_n1 /), id_var)

    ! Write data
    call write_var_primary( filename, ncid, id_var, d_tot)

    ! Close the NetCDF file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine save_variable_as_netcdf_dp_1D

  subroutine save_variable_as_netcdf_dp_2D( output_dir, d_partial, field_name)
    !< Save a single variable to a NetCDF file

    ! Input variables:
    character(len=*),         intent(in) :: output_dir
    real(dp), dimension(:,:), intent(in) :: d_partial
    character(len=*),         intent(in) :: field_name

    ! Local variables:
    character(len=256), parameter         :: routine_name = 'save_variable_as_netcdf_dp_2D'
    character(len=256)                    :: filename
    integer                               :: ncid
    integer                               :: n_partial, n_tot, n2
    integer                               :: ierr
    real(dp), dimension(:,:), allocatable :: d_tot
    integer                               :: id_dim_n1, id_dim_n2
    integer                               :: id_var

    ! Add routine to path
    call init_routine( routine_name)

    ! Determine file name
    filename = trim( output_dir) // '/' // trim( field_name) // '.nc'
    call delete_existing_file( filename)

    ! Create a new NetCDF file
    call create_new_netcdf_file_for_writing( filename, ncid)

    ! Gather data to the primary
    n_partial = size( d_partial,1)
    n2        = size( d_partial,2)
    call MPI_ALLREDUCE( n_partial, n_tot, 1, MPI_integer, MPI_SUM, MPI_COMM_WORLD, ierr)
    if (par%primary) then
      allocate( d_tot( n_tot, n2))
      call gather_to_primary( d_partial, d_tot)
    else
      call gather_to_primary( d_partial)
    end if

    ! Create dimensions
    call create_dimension( filename, ncid, 'n1', n_tot, id_dim_n1)
    call create_dimension( filename, ncid, 'n2', n2   , id_dim_n2)

    ! Create variable
    call create_variable( filename, ncid, field_name, NF90_DOUBLE, (/ id_dim_n1, id_dim_n2 /), id_var)

    ! Write data
    call write_var_primary( filename, ncid, id_var, d_tot)

    ! Close the NetCDF file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine save_variable_as_netcdf_dp_2D

  subroutine delete_existing_file( filename)

    ! In/output variables:
    character(len=*), intent(in) :: filename

    ! Local variables:
    logical :: file_exists

    if (par%primary) then
      inquire( exist = file_exists, file = trim( filename))
      if (file_exists) then
        call system('rm -f ' // filename)
      end if
    end if

  end subroutine delete_existing_file

end module netcdf_save_single_variables