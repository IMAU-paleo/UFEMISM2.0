module laddie_operators

  ! The main laddie model module.

! ===== Preamble =====
! ====================

  use precisions                                             , only: dp
  use mpi_basic                                              , only: par, sync
  use control_resources_and_error_messaging                  , only: crash, init_routine, finalise_routine, colour_string
  use parameters
  use mesh_types                                             , only: type_mesh
  use laddie_model_types                                     , only: type_laddie_model, type_laddie_timestep
  use mpi_distributed_memory                                 , only: gather_to_all
  use CSR_matrix_basics, only: allocate_matrix_CSR_dist, deallocate_matrix_CSR_dist, &
    add_entry_CSR_dist, add_empty_row_CSR_dist, finalise_matrix_CSR_dist
  use mesh_halo_exchange, only: exchange_halos

  implicit none

contains

! ===== Main routines =====
! =========================

  subroutine update_laddie_operators( mesh, laddie)
    ! Update matrix operators at the start of a new run

    ! In- and output variables

    type(type_mesh),                        intent(in)    :: mesh
    type(type_laddie_model),                intent(inout) :: laddie

    ! Local variables:
    character(len=256), parameter                         :: routine_name = 'update_laddie_operators'
    integer                                               :: ncols, ncols_loc, nrows, nrows_loc, nnz_per_row_est, nnz_est_proc
    integer                                               :: row, ti, n, i, vi, vj, ei, til, tir
    real(dp), dimension(3)                                :: cM_map_H_a_b
    real(dp), dimension(2)                                :: cM_map_H_a_c

    ! Add routine to path
    call init_routine( routine_name)

    call exchange_halos( mesh, laddie%mask_a)
    call exchange_halos( mesh, laddie%mask_b)

    ! Make sure to deallocate before allocating
    call deallocate_matrix_CSR_dist( laddie%M_map_H_a_b)
    call deallocate_matrix_CSR_dist( laddie%M_map_H_a_c)
    call deallocate_matrix_CSR_dist( laddie%M_map_UV_b_c)

    ! == Initialise the matrix using the native UFEMISM CSR-matrix format
    ! ===================================================================

    ! Matrix size
    ncols           = mesh%nV        ! from
    ncols_loc       = mesh%nV_loc
    nrows           = mesh%nTri      ! to
    nrows_loc       = mesh%nTri_loc
    nnz_per_row_est = 3
    nnz_est_proc    = nrows_loc * nnz_per_row_est

    call allocate_matrix_CSR_dist( laddie%M_map_H_a_b, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc, &
      pai_x = mesh%pai_V, pai_y = mesh%pai_Tri)

    ! == Calculate coefficients
    ! =========================

    do row = laddie%M_map_H_a_b%i1, laddie%M_map_H_a_b%i2

      ! The vertex represented by this matrix row
      ti = mesh%n2ti( row)

      if (laddie%mask_b( ti)) then

        ! Initialise
        cM_map_H_a_b = 0._dp

        ! Set counter of vertices to average over
        n = 0

        ! Loop over vertices
        do i = 1, 3
          vi = mesh%Tri( ti, i)
          ! Only add vertex if in mask_a
          if (laddie%mask_a( vi)) then
            ! Set weight factor
            cM_map_H_a_b( i) = 1._dp
            n = n + 1
          end if
        end do

        ! Scale weight factor by number of available vertices
        cM_map_H_a_b = cM_map_H_a_b / real( n, dp)

        ! Add weight to matrix
        do i = 1, 3
          vi = mesh%Tri( ti, i)
          call add_entry_CSR_dist( laddie%M_map_H_a_b, ti, vi, cM_map_H_a_b( i))
        end do

      else
        ! Outside laddie domain, so skip
        call add_empty_row_CSR_dist( laddie%M_map_H_a_b, ti)

      end if

    end do

    ! == Initialise the matrix using the native UFEMISM CSR-matrix format
    ! ===================================================================

    ! Matrix size
    ncols           = mesh%nV        ! from
    ncols_loc       = mesh%nV_loc
    nrows           = mesh%nE        ! to
    nrows_loc       = mesh%nE_loc
    nnz_per_row_est = 2
    nnz_est_proc    = nrows_loc * nnz_per_row_est

    call allocate_matrix_CSR_dist( laddie%M_map_H_a_c, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc, &
      pai_x = mesh%pai_V, pai_y = mesh%pai_E)

    ! == Calculate coefficients
    ! =========================

    do row = laddie%M_map_H_a_c%i1, laddie%M_map_H_a_c%i2

      ! The vertex represented by this matrix row
      ei = mesh%n2ei( row)

      ! Get neighbouring vertices
      vi = mesh%EV( ei, 1)
      vj = mesh%EV( ei, 2)

      ! Get masked average between the two vertices
      if (laddie%mask_a( vi) .and. laddie%mask_a( vj)) then
        cM_map_H_a_c = [0.5_dp, 0.5_dp]
      elseif (laddie%mask_a( vi)) then
        cM_map_H_a_c = [1._dp, 0._dp]
      elseif (laddie%mask_a( vj)) then
        cM_map_H_a_c = [0._dp, 1._dp]
      else
        cM_map_H_a_c = 0._dp
      end if

      ! Add weight to matrix
      call add_entry_CSR_dist( laddie%M_map_H_a_c, ei, vi, cM_map_H_a_c( 1))
      call add_entry_CSR_dist( laddie%M_map_H_a_c, ei, vj, cM_map_H_a_c( 2))

    end do

    ! == Initialise the matrix using the native UFEMISM CSR-matrix format
    ! ===================================================================

    ! Matrix size
    ncols           = mesh%nTri        ! from
    ncols_loc       = mesh%nTri_loc
    nrows           = mesh%nE        ! to
    nrows_loc       = mesh%nE_loc
    nnz_per_row_est = 2
    nnz_est_proc    = nrows_loc * nnz_per_row_est

    call allocate_matrix_CSR_dist( laddie%M_map_UV_b_c, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc, &
      pai_x = mesh%pai_Tri, pai_y = mesh%pai_E)

    ! == Calculate coefficients
    ! =========================

    do row = laddie%M_map_UV_b_c%i1, laddie%M_map_UV_b_c%i2

      ! The vertex represented by this matrix row
      ei = mesh%n2ei( row)

      ! Get neighbouring triangles
      til = mesh%ETri( ei, 1)
      tir = mesh%ETri( ei, 2)

      if (til == 0 .and. tir > 0) then
        ! Only triangle on right side exists
        if (laddie%mask_b( tir)) then
          ! Within laddie domain, so add
          call add_entry_CSR_dist( laddie%M_map_UV_b_c, ei, tir, 1._dp)
        else
          ! Outside laddie domain, so omit
          call add_empty_row_CSR_dist( laddie%M_map_UV_b_c, ei)
        end if
      elseif (tir == 0 .and. til > 0) then
        ! Only triangle on left side exists
        if (laddie%mask_b( til)) then
          ! Within laddie domain, so add
          call add_entry_CSR_dist( laddie%M_map_UV_b_c, ei, til, 1._dp)
        else
          ! Outside laddie domain, so omit
          call add_empty_row_CSR_dist( laddie%M_map_UV_b_c, ei)
        end if
      elseif (til > 0 .and. tir > 0) then
        ! Both triangles exist
        if (laddie%mask_b( til) .or. laddie%mask_b( tir)) then
          ! At least one traingle in laddie domain, so add average
          call add_entry_CSR_dist( laddie%M_map_UV_b_c, ei, til, 0.5_dp)
          call add_entry_CSR_dist( laddie%M_map_UV_b_c, ei, tir, 0.5_dp)
        else
          ! Both outside laddie domain, so omit
          call add_empty_row_CSR_dist( laddie%M_map_UV_b_c, ei)
        end if
      else
          call crash('something is seriously wrong with the ETri array of this mesh!')
      end if

    end do

    ! Crop matrix memory
    call finalise_matrix_CSR_dist( laddie%M_map_H_a_b)
    call finalise_matrix_CSR_dist( laddie%M_map_H_a_c)
    call finalise_matrix_CSR_dist( laddie%M_map_UV_b_c)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine update_laddie_operators

end module laddie_operators
