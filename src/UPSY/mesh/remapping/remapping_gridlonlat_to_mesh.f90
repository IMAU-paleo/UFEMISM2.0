module remapping_gridlonlat_to_mesh

  ! Create remapping objects between a lon/lat-grid and a mesh.

  use precisions, only: dp
  use mpi_basic, only: par
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use grid_types, only: type_grid_lonlat
  use mesh_types, only: type_mesh
  use remapping_types, only: type_map
  use CSR_sparse_matrix_type, only: type_sparse_matrix_CSR_dp
  use CSR_matrix_basics, only: allocate_matrix_CSR_dist, finalise_matrix_CSR_dist, &
    deallocate_matrix_CSR_dist, add_entry_CSR_dist
  use petsc_basic, only: mat_CSR2petsc

  implicit none

  private

  public :: create_map_from_lonlat_grid_to_mesh

contains

subroutine create_map_from_lonlat_grid_to_mesh( grid, mesh, map)
  ! Create a new mapping object from a lon/lat-grid to a mesh.
  !
  ! By default uses bilinear interpolation.

  ! In/output variables:
  type(type_grid_lonlat),              intent(in)    :: grid
  type(type_mesh),                     intent(in)    :: mesh
  type(type_map),                      intent(inout) :: map

  ! Local variables:
  character(len=1024), parameter                     :: routine_name = 'create_map_from_lonlat_grid_to_mesh'
  integer                                            :: nrows, ncols, nrows_loc, ncols_loc, nnz_est, nnz_est_proc, nnz_per_row_max
  type(type_sparse_matrix_CSR_dp)                    :: M_CSR
  integer                                            :: vi
  integer                                            :: il,iu,jl,ju
  real(dp)                                           :: wil,wiu,wjl,wju

  ! Add routine to path
  call init_routine( routine_name)

  ! Safety
  if (map%is_in_use) call crash('this map is already in use!')

  ! == Initialise map metadata
  ! ==========================

  map%is_in_use = .true.
  map%name_src  = grid%name
  map%name_dst  = mesh%name
  map%method    = 'bilin'

  ! == Initialise the mapping matrix using the native UFEMISM CSR-matrix format
  ! ===========================================================================

  ! Matrix size
  nrows           = mesh%nV  ! to
  nrows_loc       = mesh%nV_loc
  ncols           = grid%n   ! from
  ncols_loc       = grid%n_loc
  nnz_per_row_max = 4
  nnz_est         = nnz_per_row_max * nrows
  nnz_est_proc    = ceiling( real( nnz_est, dp) / real( par%n, dp))

  call allocate_matrix_CSR_dist( M_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

  ! Fill in the CSR matrix
  do vi = mesh%vi1, mesh%vi2

    ! Find enveloping lat-lon indices
    il  = max(1, min( grid%nlon-1, 1 + floor((mesh%lon( vi) - minval(grid%lon)) / (grid%lon(2) - grid%lon(1)))))
    iu  = il + 1
    wil = (grid%lon(iu) - mesh%lon( vi)) / (grid%lon(2) - grid%lon(1))
    wiu = 1._dp - wil

    ! Exception for pixels near the zero meridian
    if (mesh%lon( vi) < minval(grid%lon)) then
      il  = grid%nlon
      iu  = 1
      wil = (grid%lon( iu) - mesh%lon( vi)) / (grid%lon(2) - grid%lon(1))
      wiu = 1._dp - wil
    elseif (mesh%lon( vi) > maxval(grid%lon)) then
      il  = grid%nlon
      iu  = 1
      wiu = (mesh%lon( vi) - grid%lon( il)) / (grid%lon(2) - grid%lon(1))
      wil = 1._dp - wiu
    end if

    jl  = max(1, min( grid%nlat-1, 1 + floor((mesh%lat( vi) - minval(grid%lat)) / (grid%lat(2) - grid%lat(1)))))
    ju  = jl + 1
    wjl = (grid%lat( ju) - mesh%lat( vi)) / (grid%lat(2) - grid%lat(1))
    wju = 1 - wjl

    ! Add values to the CSR matrix
    call add_entry_CSR_dist( M_CSR, vi, grid%ij2n( il,jl), wil * wjl)
    call add_entry_CSR_dist( M_CSR, vi, grid%ij2n( il,ju), wil * wju)
    call add_entry_CSR_dist( M_CSR, vi, grid%ij2n( iu,jl), wiu * wjl)
    call add_entry_CSR_dist( M_CSR, vi, grid%ij2n( iu,ju), wiu * wju)

  end do

  call finalise_matrix_CSR_dist( M_CSR)

  ! Convert matrices from Fortran to PETSc types
  call mat_CSR2petsc( M_CSR, map%M)

  ! Clean up the Fortran versions
  call deallocate_matrix_CSR_dist( M_CSR)

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine create_map_from_lonlat_grid_to_mesh

end module remapping_gridlonlat_to_mesh
