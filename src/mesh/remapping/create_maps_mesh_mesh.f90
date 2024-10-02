module create_maps_mesh_mesh

  ! Create remapping objects between two meshes.

#include <petsc/finclude/petscksp.h>
  use petscksp
  use precisions, only: dp
  use mpi_basic, only: par, sync
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use mesh_types, only: type_mesh
  use remapping_types, only: type_map, type_single_row_mapping_matrices
  use CSR_sparse_matrix_utilities, only: type_sparse_matrix_CSR_dp, allocate_matrix_CSR_dist, &
    add_empty_row_CSR_dist, add_entry_CSR_dist, deallocate_matrix_CSR_dist
  use math_utilities, only: triangle_area
  use mesh_utilities, only: calc_Voronoi_cell, find_containing_triangle, find_containing_vertex
  use petsc_basic, only: mat_CSR2petsc, mat_petsc2CSR, MatConvert
  use line_tracing_triangles, only: trace_line_tri
  use line_tracing_Voronoi, only: trace_line_Vor

  implicit none

  private

  public :: create_map_from_mesh_to_mesh_nearest_neighbour
  public :: create_map_from_mesh_to_mesh_trilin
  public :: create_map_from_mesh_to_mesh_2nd_order_conservative

contains

  !> Create a new mapping object from a mesh to a mesh using nearest-neighbour interpolation.
  subroutine create_map_from_mesh_to_mesh_nearest_neighbour( mesh_src, mesh_dst, map)

    ! In/output variables
    type(type_mesh), intent(in)    :: mesh_src
    type(type_mesh), intent(in)    :: mesh_dst
    type(type_map),  intent(inout) :: map

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'create_map_from_mesh_to_mesh_nearest_neighbour'
    integer                         :: ncols, nrows, nrows_loc, ncols_loc, nnz_per_row_max, nnz_est_proc
    type(type_sparse_matrix_CSR_dp) :: M_CSR
    integer                         :: row, vi_dst
    real(dp), dimension(2)          :: p
    integer                         :: vi_src, col

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (map%is_in_use) call crash('this map is already in use!')

    ! == Initialise map metadata
    ! ==========================

    map%is_in_use = .true.
    map%name_src  = mesh_src%name
    map%name_dst  = mesh_dst%name
    map%method    = 'nearest_neighbour'

    ! == Initialise the matrix using the native UFEMISM CSR-matrix format
    ! ===================================================================

    ! Matrix size
    nrows           = mesh_dst%nV   ! to
    nrows_loc       = mesh_dst%nV_loc
    ncols           = mesh_src%nV   ! from
    ncols_loc       = mesh_src%nV_loc
    nnz_per_row_max = 1
    nnz_est_proc    = nnz_per_row_max * nrows_loc

    call allocate_matrix_CSR_dist( M_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! For all mesh_dst vertices, find the mesh_src triangle containing them
    vi_src = 1
    do row = mesh_dst%vi1, mesh_dst%vi2

      vi_dst = mesh_dst%n2vi( row)

      p = mesh_dst%V( vi_dst,:)
      call find_containing_vertex( mesh_src, p, vi_src)

      col = mesh_src%vi2n( vi_src)

      ! Add to the matrix
      call add_entry_CSR_dist( M_CSR, row, col, 1._dp)

    end do ! do row = mesh_dst%vi1, mesh_dst%vi2
    call sync

    ! Convert matrices from Fortran to PETSc types
    call mat_CSR2petsc( M_CSR, map%M)

    ! Clean up after yourself
    call deallocate_matrix_CSR_dist( M_CSR)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_map_from_mesh_to_mesh_nearest_neighbour

  !> Create a new mapping object from a mesh to a mesh using trilinear interpolation.
  subroutine create_map_from_mesh_to_mesh_trilin( mesh_src, mesh_dst, map)

    ! In/output variables
    type(type_mesh), intent(in)    :: mesh_src
    type(type_mesh), intent(in)    :: mesh_dst
    type(type_map),  intent(inout) :: map

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'create_map_from_mesh_to_mesh_trilin'
    integer                         :: ncols, nrows, nrows_loc, ncols_loc, nnz_per_row_max, nnz_est_proc
    type(type_sparse_matrix_CSR_dp) :: M_CSR
    integer                         :: row, vi_dst
    real(dp), dimension(2)          :: p
    integer                         :: ti_src, via, vib, vic
    real(dp), dimension(2)          :: pa, pb, pc
    real(dp)                        :: Atri_abp, Atri_bcp, Atri_cap, Atri_abc, wa, wb, wc
    integer                         :: cola, colb, colc

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (map%is_in_use) call crash('this map is already in use!')

    ! == Initialise map metadata
    ! ==========================

    map%is_in_use = .true.
    map%name_src  = mesh_src%name
    map%name_dst  = mesh_dst%name
    map%method    = 'trilin'

    ! == Initialise the matrix using the native UFEMISM CSR-matrix format
    ! ===================================================================

    ! Matrix size
    nrows           = mesh_dst%nV   ! to
    nrows_loc       = mesh_dst%nV_loc
    ncols           = mesh_src%nV   ! from
    ncols_loc       = mesh_src%nV_loc
    nnz_per_row_max = 3
    nnz_est_proc    = nnz_per_row_max * nrows_loc

    call allocate_matrix_CSR_dist( M_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! For all mesh_dst vertices, find the mesh_src triangle containing them
    ti_src = 1
    do row = mesh_dst%vi1, mesh_dst%vi2

      vi_dst = mesh_dst%n2vi( row)

      p = mesh_dst%V( vi_dst,:)
      call find_containing_triangle( mesh_src, p, ti_src)

      ! Calculate the trilinear interpolation weights
      via = mesh_src%Tri( ti_src,1)
      vib = mesh_src%Tri( ti_src,2)
      vic = mesh_src%Tri( ti_src,3)

      pa  = mesh_src%V( via,:)
      pb  = mesh_src%V( vib,:)
      pc  = mesh_src%V( vic,:)

      Atri_abp = triangle_area( pa, pb, p)
      Atri_bcp = triangle_area( pb, pc, p)
      Atri_cap = triangle_area( pc, pa, p)
      Atri_abc = Atri_abp + Atri_bcp + Atri_cap

      wa = Atri_bcp / Atri_abc
      wb = Atri_cap / Atri_abc
      wc = Atri_abp / Atri_abc

      ! Matrix columns corresponding to these three vertices
      cola = mesh_src%vi2n( via)
      colb = mesh_src%vi2n( vib)
      colc = mesh_src%vi2n( vic)

      ! Add to the matrix
      call add_entry_CSR_dist( M_CSR, row, cola, wa)
      call add_entry_CSR_dist( M_CSR, row, colb, wb)
      call add_entry_CSR_dist( M_CSR, row, colc, wc)

    end do ! do row = mesh_dst%vi1, mesh_dst%vi2
    call sync

    ! Convert matrices from Fortran to PETSc types
    call mat_CSR2petsc( M_CSR, map%M)

    ! Clean up after yourself
    call deallocate_matrix_CSR_dist( M_CSR)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_map_from_mesh_to_mesh_trilin

  !> Create a new mapping object from a mesh to a mesh using 2nd-order conservative interpolation.
  subroutine create_map_from_mesh_to_mesh_2nd_order_conservative( mesh_src, mesh_dst, map)

    ! In/output variables
    type(type_mesh), intent(in)    :: mesh_src
    type(type_mesh), intent(in)    :: mesh_dst
    type(type_map),  intent(inout) :: map

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'create_map_from_mesh_to_mesh_2nd_order_conservative'
    type(PetscErrorCode)                    :: perr
    logical                                 :: count_coincidences
    integer                                 :: nnz_per_row_max
    type(tMat)                              :: B_xdy_b_a  , B_mxydx_b_a  , B_xydy_b_a
    type(tMat)                              :: B_xdy_a_b  , B_mxydx_a_b  , B_xydy_a_b
    type(tMat)                              :: B_xdy_b_a_T, B_mxydx_b_a_T, B_xydy_b_a_T
    type(tMat)                              :: w0, w1x, w1y
    integer                                 :: istart, iend, n, k, ti
    integer                                 :: ncols
    integer,  dimension(:    ), allocatable :: cols
    real(dp), dimension(:    ), allocatable :: vals, w0_row, w1x_row, w1y_row
    real(dp)                                :: A_overlap_tot
    type(tMat)                              :: M_map_a_b, M_ddx_a_b, M_ddy_a_b
    type(tMat)                              :: M1, M2, M_cons_1st_order

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (map%is_in_use) call crash('this map is already in use!')

    ! == Initialise map metadata
    ! ==========================

    map%is_in_use = .true.
    map%name_src  = mesh_src%name
    map%name_dst  = mesh_dst%name
    map%method    = '2nd_order_conservative'

    ! Integrate around the Voronoi cells of the destination mesh through the triangles of the source mesh
    count_coincidences = .true.
    call integrate_Voronoi_cells_through_triangles( mesh_dst, mesh_src, B_xdy_a_b, B_mxydx_a_b, B_xydy_a_b, count_coincidences)

    ! Integrate around the triangles of the source mesh through the Voronoi cells of the destination mesh
    count_coincidences = .false.
    call integrate_triangles_through_Voronoi_cells( mesh_src, mesh_dst, B_xdy_b_a, B_mxydx_b_a, B_xydy_b_a, count_coincidences)

    ! Transpose line integral matrices
    !if (par%master) WRITE(0,*) 'calc_remapping_operators_mesh_mesh_conservative - transposing line integral matrices...'
    call MatCreateTranspose( B_xdy_b_a  , B_xdy_b_a_T  , perr)
    call MatCreateTranspose( B_mxydx_b_a, B_mxydx_b_a_T, perr)
    call MatCreateTranspose( B_xydy_b_a , B_xydy_b_a_T , perr)

    ! Combine line integrals around areas of overlap to get surface integrals over areas of overlap
    call MatAXPY( B_xdy_a_b  , 1._dp, B_xdy_b_a_T  , UNKNOWN_NONZERO_PATTERN, perr)
    call MatAXPY( B_mxydx_a_b, 1._dp, B_mxydx_b_a_T, UNKNOWN_NONZERO_PATTERN, perr)
    call MatAXPY( B_xydy_a_b , 1._dp, B_xydy_b_a_T , UNKNOWN_NONZERO_PATTERN, perr)

    call MatDestroy( B_xdy_b_a_T  , perr)
    call MatDestroy( B_mxydx_b_a_T, perr)
    call MatDestroy( B_xydy_b_a_T , perr)

    ! Calculate w0, w1x, w1y for the mesh-to-grid remapping operator
    call MatConvert( B_xdy_a_b, MATAIJ, MAT_INITIAL_MATRIX, w0, perr)
    call MatConvert( B_xdy_a_b, MATAIJ, MAT_INITIAL_MATRIX, w1x, perr)
    call MatConvert( B_xdy_a_b, MATAIJ, MAT_INITIAL_MATRIX, w1y, perr)

    ! Estimate maximum number of non-zeros per row (i.e. maximum number of grid cells overlapping with a mesh triangle)
    nnz_per_row_max = max( 32, max( ceiling( 2._dp * maxval( mesh_src%TriA) / minval( mesh_dst%A   )), &
                                    ceiling( 2._dp * maxval( mesh_dst%A   ) / minval( mesh_src%TriA))) )

    ! allocate memory for a single matrix row
    allocate( cols(    nnz_per_row_max))
    allocate( vals(    nnz_per_row_max))
    allocate( w0_row(  nnz_per_row_max))
    allocate( w1x_row( nnz_per_row_max))
    allocate( w1y_row( nnz_per_row_max))

    call MatGetOwnershipRange( B_xdy_a_b  , istart, iend, perr)

    do n = istart+1, iend ! +1 because PETSc indexes from 0

      ! w0
      call MatGetRow( B_xdy_a_b, n-1, ncols, cols, vals, perr)
      A_overlap_tot = sum( vals( 1:ncols))
      do k = 1, ncols
        w0_row( k) = vals( k) / A_overlap_tot
        call MatSetValues( w0, 1, n-1, 1, cols( k), w0_row( k), INSERT_VALUES, perr)
      end do
      call MatRestoreRow( B_xdy_a_b, n-1, ncols, cols, vals, perr)

      ! w1x
      call MatGetRow( B_mxydx_a_b, n-1, ncols, cols, vals, perr)
      do k = 1, ncols
        ti = cols( k)+1
        w1x_row( k) = (vals( k) / A_overlap_tot) - (mesh_src%TriGC( ti,1) * w0_row( k))
        call MatSetValues( w1x, 1, n-1, 1, cols( k), w1x_row( k), INSERT_VALUES, perr)
      end do
      call MatRestoreRow( B_mxydx_a_b, n-1, ncols, cols, vals, perr)

      ! w1y
      call MatGetRow( B_xydy_a_b, n-1, ncols, cols, vals, perr)
      do k = 1, ncols
        ti = cols( k)+1
        w1y_row( k) = (vals( k) / A_overlap_tot) - (mesh_src%TriGC( ti,2) * w0_row( k))
        call MatSetValues( w1y, 1, n-1, 1, cols( k), w1y_row( k), INSERT_VALUES, perr)
      end do
      call MatRestoreRow( B_xydy_a_b, n-1, ncols, cols, vals, perr)

    end do
    call MatAssemblyBegin( w0, MAT_FINAL_ASSEMBLY, perr)
    call MatAssemblyEnd(   w0, MAT_FINAL_ASSEMBLY, perr)
    call MatAssemblyBegin( w1x, MAT_FINAL_ASSEMBLY, perr)
    call MatAssemblyEnd(   w1x, MAT_FINAL_ASSEMBLY, perr)
    call MatAssemblyBegin( w1y, MAT_FINAL_ASSEMBLY, perr)
    call MatAssemblyEnd(   w1y, MAT_FINAL_ASSEMBLY, perr)
    call sync

    call MatDestroy( B_xdy_a_b  , perr)
    call MatDestroy( B_mxydx_a_b, perr)
    call MatDestroy( B_xydy_a_b , perr)

    ! == Calculate the remapping matrices
    ! ===================================

    ! Safety
    if (.not. allocateD( mesh_src%vi2n)) then
      call crash('matrix operators for mesh "' // trim( mesh_src%name) // '" have not been calculated!')
    end if

    ! Convert matrices to PETSc format
    call mat_CSR2petsc( mesh_src%M_map_a_b, M_map_a_b)
    call mat_CSR2petsc( mesh_src%M_ddx_a_b, M_ddx_a_b)
    call mat_CSR2petsc( mesh_src%M_ddy_a_b, M_ddy_a_b)

    ! 1st-order = w0 * map_a_b
    call MatMatMult( w0 , M_map_a_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_real, M_cons_1st_order, perr)

    ! 2nd-order = 1st-order + w1x * ddx_a_b + w1y * ddy_a_b
    call MatMatMult( w1x, M_ddx_a_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_real, M1, perr)  ! This can be done more efficiently now that the non-zero structure is known...
    call MatMatMult( w1y, M_ddy_a_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_real, M2, perr)

    call MatConvert( M_cons_1st_order, MATAIJ, MAT_INITIAL_MATRIX, map%M, perr)
    call MatAXPY( map%M, 1._dp, M1, DifFERENT_NONZERO_PATTERN, perr)
    call MatAXPY( map%M, 1._dp, M2, DifFERENT_NONZERO_PATTERN, perr)

    call MatDestroy( w0       , perr)
    call MatDestroy( w1x      , perr)
    call MatDestroy( w1y      , perr)
    call MatDestroy( M_map_a_b, perr)
    call MatDestroy( M_ddx_a_b, perr)
    call MatDestroy( M_ddy_a_b, perr)

    call MatDestroy( M1              , perr)
    call MatDestroy( M2              , perr)

    ! == Apply some final corrections
    ! ===============================

    call correct_mesh_to_mesh_map( mesh_src, mesh_dst, M_cons_1st_order, map%M)

    call MatDestroy( M_cons_1st_order, perr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_map_from_mesh_to_mesh_2nd_order_conservative

  !> Apply some final corrections to the 2nd-order conservative mesh-to-mesh remapping operator:
  !> - set remapped data to zero on the domain border
  !> - use direct copying for identical vertices
  subroutine correct_mesh_to_mesh_map( mesh_src, mesh_dst, M_cons_1st_order, M_cons_2nd_order)

    ! In/output variables
    type(type_mesh), intent(in)    :: mesh_src
    type(type_mesh), intent(in)    :: mesh_dst
    type(tMat),      intent(in)    :: M_cons_1st_order
    type(tMat),      intent(inout) :: M_cons_2nd_order

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'correct_mesh_to_mesh_map'
    integer                                 :: perr
    type(type_sparse_matrix_CSR_dp)         :: M_cons_1st_order_CSR
    type(type_sparse_matrix_CSR_dp)         :: M_cons_2nd_order_CSR
    integer                                 :: i, vi_dst, k1, k2, k
    integer                                 :: j, vi_src
    logical                                 :: do_direct_copy
    integer                                 :: vi_src_copy
    logical                                 :: Voronoi_cells_are_identical
    real(dp), dimension( mesh_src%nC_mem,2) :: Vor_src   , Vor_dst
    integer,  dimension( mesh_src%nC_mem  ) :: Vor_src_vi, Vor_dst_vi
    integer,  dimension( mesh_src%nC_mem  ) :: Vor_src_ti, Vor_dst_ti
    integer                                 :: nVor_src  , nVor_dst
    integer                                 :: vori
    type(type_map)                          :: map_trilin
    type(type_sparse_matrix_CSR_dp)         :: M_trilin_CSR
    logical,  dimension( mesh_dst%nV)       :: isgood_1st_order
    logical,  dimension( mesh_dst%nV)       :: isgood_2nd_order
    integer                                 :: kk1,kk2,kk

    ! Add routine to path
    call init_routine( routine_name)

    ! Convert matrices to CSR format for easier handling
    call mat_petsc2CSR( M_cons_1st_order, M_cons_1st_order_CSR)
    call mat_petsc2CSR( M_cons_2nd_order, M_cons_2nd_order_CSR)

    ! == Set to zero
    !
    ! 2nd-order conservative doesn't work all that well on the
    ! domain border; just set the result to zero there.

    do i = M_cons_2nd_order_CSR%i1, M_cons_2nd_order_CSR%i2

      k1 = M_cons_2nd_order_CSR%ptr( i)
      k2 = M_cons_2nd_order_CSR%ptr( i+1) - 1

      vi_dst = mesh_dst%n2vi( i)

      if (mesh_dst%VBI( vi_dst) > 0) then
        do k = k1, k2
          M_cons_2nd_order_CSR%val( k) = 0._dp
        end do
      end if

    end do ! do i = M_cons_2nd_order_CSR%i1, M_cons_2nd_order_CSR%i2

    ! == Direct copying
    !
    ! With the new mesh generation code (of UFE2.0), many vertices away from the grounding line
    ! remain unchanged after a mesh update. The vertex-to-triangle-to-vertex remapping is slightly
    ! diffusive, so instead we can just copy data directly for those vertices.

    do i = M_cons_2nd_order_CSR%i1, M_cons_2nd_order_CSR%i2

      k1 = M_cons_2nd_order_CSR%ptr( i)
      k2 = M_cons_2nd_order_CSR%ptr( i+1) - 1

      vi_dst = mesh_dst%n2vi( i)

      do_direct_copy = .false.
      vi_src_copy    = 0

      ! Loop over all src vertices contributing to this dst vertex
      do k = k1, k2

        j = M_cons_2nd_order_CSR%ind( k)
        vi_src = mesh_src%n2vi( j)

        if (norm2( mesh_dst%V( vi_dst,:) - mesh_src%V( vi_src,:)) < mesh_dst%R( vi_dst) / 1E2_dp) then
          ! These vertices coincide; check if their Voronoi cells are identical

          Voronoi_cells_are_identical = .true.

          call calc_Voronoi_cell( mesh_src, vi_src, 0._dp, Vor_src, Vor_src_vi, Vor_src_ti, nVor_src)
          call calc_Voronoi_cell( mesh_dst, vi_dst, 0._dp, Vor_dst, Vor_dst_vi, Vor_dst_ti, nVor_dst)

          if (nVor_src /= nVor_dst) then
            Voronoi_cells_are_identical = .false.
          else
            do vori = 1, nVor_src
              if (norm2( Vor_src( vori,:) - Vor_dst( vori,:)) > mesh_dst%R( vi_dst) / 1E2_dp) then
                Voronoi_cells_are_identical = .false.
              end if
            end do
          end if

          if (Voronoi_cells_are_identical) then
            ! These two vertices have identical Voronoi cells; use direct copying
            do_direct_copy = .true.
            vi_src_copy    = vi_src
            exit
          end if ! if (Voronoi_cells_are_identical) then

        end if ! if (norm2( mesh_dst%V( vi_dst,:) - mesh_src%V( vi_src,:)) < mesh_dst%tol_dist) then

      end do ! do k = k1, k2

      ! if a source vertex with an identical Voronoi cell to this dst vertex was
      ! found, copy data from that vertex directly
      if (do_direct_copy) then
        ! Loop over all src vertices contributing to this dst vertex; set all
        ! contributions to zero except the one we're copying (which is set to 1)

        do k = k1, k2

          j = M_cons_2nd_order_CSR%ind( k)
          vi_src = mesh_src%n2vi( j)

          if (vi_src == vi_src_copy) then
            M_cons_2nd_order_CSR%val( k) = 1._dp
          else
            M_cons_2nd_order_CSR%val( k) = 0._dp
          end if

        end do ! do k = k1, k2

      end if ! if (do_direct_copy) then

    end do ! do i = M_cons_2nd_order_CSR%i1, M_cons_2nd_order_CSR%i2

    ! == Remapping errors
    !
    ! On very rare occasions, the remapping operator is just wrong, likely due to round-off
    ! errors in determining if vertices coincide or not. Usually, the 1st-order operator
    ! is fine. if that one fails too, just replace the answer with trilinear interpolation.
    !
    ! Faulty operators can be detected by negative coefficients in the remapping matrix,
    ! which automatically violate conservation of extreme values.

    ! Calculate the trilinear interpolation operator to serve as a back-up
    call create_map_from_mesh_to_mesh_trilin( mesh_src, mesh_dst, map_trilin)
    call mat_petsc2CSR( map_trilin%M, M_trilin_CSR)

    ! Find faulty operators in the 1st-order conservative remapping operator
    isgood_1st_order = .true.
    do i = M_cons_1st_order_CSR%i1, M_cons_1st_order_CSR%i2

      k1 = M_cons_1st_order_CSR%ptr( i)
      k2 = M_cons_1st_order_CSR%ptr( i+1) - 1

      vi_dst = mesh_dst%n2vi( i)

      do k = k1, k2
        if (M_cons_1st_order_CSR%val( k) < 0._dp) then
          isgood_1st_order( vi_dst) = .false.
        end if
      end do

    end do ! do i = M_cons_1st_order_CSR%i1, M_cons_1st_order_CSR%i2

    ! Find faulty operators in the 2nd-order conservative remapping operator
    isgood_2nd_order = .true.
    do i = M_cons_2nd_order_CSR%i1, M_cons_2nd_order_CSR%i2

      k1 = M_cons_2nd_order_CSR%ptr( i)
      k2 = M_cons_2nd_order_CSR%ptr( i+1) - 1

      vi_dst = mesh_dst%n2vi( i)

      do k = k1, k2
        if (M_cons_2nd_order_CSR%val( k) < 0._dp) then
          isgood_2nd_order( vi_dst) = .false.
        end if
      end do

    end do ! do i = M_cons_2nd_order_CSR%i1, M_cons_2nd_order_CSR%i2

    ! Replace faulty operators in the 2nd-order conservative remapping operator

    do i = M_cons_2nd_order_CSR%i1, M_cons_2nd_order_CSR%i2

      vi_dst = mesh_dst%n2vi( i)

      if (.not. isgood_2nd_order( vi_dst)) then
        ! Replace this faulty operator

        if (isgood_1st_order( vi_dst)) then
          ! Replace with the 1st-order conservative remapping operator

          ! First set all coefficients of the 2nd-order operator to zero
          k1 = M_cons_2nd_order_CSR%ptr( i)
          k2 = M_cons_2nd_order_CSR%ptr( i+1) - 1
          do k = k1, k2
            M_cons_2nd_order_CSR%val( k) = 0._dp
          end do

          ! then copy values from the 1st-order operator
          kk1 = M_cons_1st_order_CSR%ptr( i)
          kk2 = M_cons_1st_order_CSR%ptr( i+1) - 1
          do kk = kk1, kk2
            k = kk - kk1 + k1
            M_cons_2nd_order_CSR%ind( k) = M_cons_1st_order_CSR%ind( kk)
            M_cons_2nd_order_CSR%val( k) = M_cons_1st_order_CSR%val( kk)
          end do

        else ! if (isgood_1st_order( vi_dst)) then
          ! Replace with the trilinear interpolation operator

          ! First set all coefficients of the 2nd-order operator to zero
          k1 = M_cons_2nd_order_CSR%ptr( i)
          k2 = M_cons_2nd_order_CSR%ptr( i+1) - 1
          do k = k1, k2
            M_cons_2nd_order_CSR%val( k) = 0._dp
          end do

          ! then copy values from the trilinear interpolation operator
          kk1 = M_trilin_CSR%ptr( i)
          kk2 = M_trilin_CSR%ptr( i+1) - 1
          do kk = kk1, kk2
            k = kk - kk1 + k1
            M_cons_2nd_order_CSR%ind( k) = M_trilin_CSR%ind( kk)
            M_cons_2nd_order_CSR%val( k) = M_trilin_CSR%val( kk)
          end do

        end if ! if (isgood_1st_order( vi_dst)) then

      end if ! if (.not. isgood_2nd_order( vi_dst)) then

    end do ! do i = M_cons_2nd_order_CSR%i1, M_cons_2nd_order_CSR%i2

    ! Convert back to PETSc format
    call MatDestroy( M_cons_2nd_order, perr)
    call mat_CSR2petsc( M_cons_2nd_order_CSR, M_cons_2nd_order)

    ! Clean up after yourself
    call deallocate_matrix_CSR_dist( M_cons_1st_order_CSR)
    call deallocate_matrix_CSR_dist( M_cons_2nd_order_CSR)
    call deallocate_matrix_CSR_dist( M_trilin_CSR)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine correct_mesh_to_mesh_map

  !> Integrate around the triangles of mesh_tri through the Voronoi cells of mesh_Vor
  subroutine integrate_triangles_through_Voronoi_cells( mesh_tri, mesh_Vor, &
    B_xdy_b_a, B_mxydx_b_a, B_xydy_b_a, count_coincidences)

    ! In/output variables
    type(type_mesh), intent(in)    :: mesh_tri
    type(type_mesh), intent(in)    :: mesh_Vor
    type(tMat),      intent(out)   :: B_xdy_b_a
    type(tMat),      intent(out)   :: B_mxydx_b_a
    type(tMat),      intent(out)   :: B_xydy_b_a
    logical,         intent(in)    :: count_coincidences

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'integrate_triangles_through_Voronoi_cells'
    integer                                :: nrows, ncols, nrows_loc, ncols_loc, nnz_est, nnz_est_proc, nnz_per_row_max
    type(type_sparse_matrix_CSR_dp)        :: B_xdy_b_a_CSR, B_mxydx_b_a_CSR, B_xydy_b_a_CSR
    type(type_single_row_mapping_matrices) :: single_row
    integer                                :: via, vib, vic, ti, vi_hint, k
    real(dp), dimension(2)                 :: p, q

    ! Add routine to path
    call init_routine( routine_name)

    ! == Initialise the three matrices using the native UFEMISM CSR-matrix format
    ! ===========================================================================

    ! Matrix sise
    nrows           = mesh_tri%nTri  ! to
    nrows_loc       = mesh_tri%nTri_loc
    ncols           = mesh_Vor%nV    ! from
    ncols_loc       = mesh_Vor%nV_loc
    nnz_est         = 4 * max( nrows, ncols)
    nnz_est_proc    = ceiling( real( nnz_est, dp) / real( par%n, dp))
    nnz_per_row_max = max( 32, max( ceiling( 2._dp * maxval( mesh_tri%TriA) / minval( mesh_Vor%A   )), &
                                    ceiling( 2._dp * maxval( mesh_Vor%A   ) / minval( mesh_tri%TriA)) ))

    call allocate_matrix_CSR_dist( B_xdy_b_a_CSR  , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( B_mxydx_b_a_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( B_xydy_b_a_CSR , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! Initialise results from integrating a single triangle through the Voronoi cells
    single_row%n_max = 100
    single_row%n     = 0
    allocate( single_row%index_left( single_row%n_max))
    allocate( single_row%LI_xdy(     single_row%n_max))
    allocate( single_row%LI_mxydx(   single_row%n_max))
    allocate( single_row%LI_xydy(    single_row%n_max))

    ! == Trace all the line segments to fill the matrices
    ! ===================================================

    vi_hint = 1

    do ti = mesh_tri%ti1, mesh_tri%ti2

      ! Clean up single row results
      single_row%n            = 0
      single_row%index_left   = 0
      single_row%LI_xdy       = 0
      single_row%LI_mxydx     = 0
      single_row%LI_xydy      = 0

      ! The three vertices spanning this triangle
      via = mesh_tri%Tri( ti,1)
      vib = mesh_tri%Tri( ti,2)
      vic = mesh_tri%Tri( ti,3)

      ! Integrate over the three triangle sides
      p = mesh_tri%V( via,:)
      q = mesh_tri%V( vib,:)
      call trace_line_Vor( mesh_Vor, p, q, single_row, count_coincidences, vi_hint)

      p = mesh_tri%V( vib,:)
      q = mesh_tri%V( vic,:)
      call trace_line_Vor( mesh_Vor, p, q, single_row, count_coincidences, vi_hint)

      p = mesh_tri%V( vic,:)
      q = mesh_tri%V( via,:)
      call trace_line_Vor( mesh_Vor, p, q, single_row, count_coincidences, vi_hint)

      ! Add the results for this triangle to the sparse matrix
      if (single_row%n == 0) then
        call add_empty_row_CSR_dist( B_xdy_b_a_CSR  , ti)
        call add_empty_row_CSR_dist( B_mxydx_b_a_CSR, ti)
        call add_empty_row_CSR_dist( B_xydy_b_a_CSR , ti)
      else
        do k = 1, single_row%n
          call add_entry_CSR_dist( B_xdy_b_a_CSR  , ti, single_row%index_left( k), single_row%LI_xdy(   k))
          call add_entry_CSR_dist( B_mxydx_b_a_CSR, ti, single_row%index_left( k), single_row%LI_mxydx( k))
          call add_entry_CSR_dist( B_xydy_b_a_CSR , ti, single_row%index_left( k), single_row%LI_xydy(  k))
        end do
      end if

    end do ! do ti = mesh_tri%ti1, mesh_tri%ti2
    call sync

    ! Convert matrices from Fortran to PETSc types
    call mat_CSR2petsc( B_xdy_b_a_CSR  , B_xdy_b_a  )
    call mat_CSR2petsc( B_mxydx_b_a_CSR, B_mxydx_b_a)
    call mat_CSR2petsc( B_xydy_b_a_CSR , B_xydy_b_a )

    ! Clean up the Fortran versions
    call deallocate_matrix_CSR_dist( B_xdy_b_a_CSR  )
    call deallocate_matrix_CSR_dist( B_mxydx_b_a_CSR)
    call deallocate_matrix_CSR_dist( B_xydy_b_a_CSR )

    ! Clean up after yourself
    deallocate( single_row%index_left )
    deallocate( single_row%LI_xdy     )
    deallocate( single_row%LI_mxydx   )
    deallocate( single_row%LI_xydy    )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine integrate_triangles_through_Voronoi_cells

  !> Integrate around the grid cells of the grid through the triangles of the mesh
  subroutine integrate_Voronoi_cells_through_triangles( mesh_Vor, mesh_tri, &
    B_xdy_a_b, B_mxydx_a_b, B_xydy_a_b, count_coincidences)

    ! In/output variables
    type(type_mesh), intent(in)    :: mesh_Vor
    type(type_mesh), intent(in)    :: mesh_tri
    type(tMat),      intent(out)   :: B_xdy_a_b
    type(tMat),      intent(out)   :: B_mxydx_a_b
    type(tMat),      intent(out)   :: B_xydy_a_b
    logical,         intent(in)    :: count_coincidences

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'integrate_Voronoi_cells_through_triangles'
    integer                                 :: nrows, ncols, nrows_loc, ncols_loc, nnz_est, nnz_est_proc, nnz_per_row_max
    type(type_sparse_matrix_CSR_dp)         :: B_xdy_a_b_CSR, B_mxydx_a_b_CSR, B_xydy_a_b_CSR
    type(type_single_row_mapping_matrices)  :: single_row
    integer                                 :: vi, vori1, vori2, k, ti_hint
    real(dp), dimension( mesh_Vor%nC_mem,2) :: Vor
    integer,  dimension( mesh_Vor%nC_mem  ) :: Vor_vi
    integer,  dimension( mesh_Vor%nC_mem  ) :: Vor_ti
    integer                                 :: nVor
    real(dp), dimension(2)                  :: p, q

    ! Add routine to path
    call init_routine( routine_name)

    ! == Initialise the three matrices using the native UFEMISM CSR-matrix format
    ! ===========================================================================

    ! Matrix sise
    nrows           = mesh_Vor%nV    ! to
    nrows_loc       = mesh_Vor%nV_loc
    ncols           = mesh_tri%nTri  ! from
    ncols_loc       = mesh_tri%nTri_loc
    nnz_est         = 4 * max( nrows, ncols)
    nnz_est_proc    = ceiling( real( nnz_est, dp) / real( par%n, dp))
    nnz_per_row_max = max( 32, max( ceiling( 2._dp * maxval( mesh_tri%TriA) / minval( mesh_vor%A   )), &
                                    ceiling( 2._dp * maxval( mesh_vor%A   ) / minval( mesh_tri%TriA)) ))

    call allocate_matrix_CSR_dist( B_xdy_a_b_CSR  , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( B_mxydx_a_b_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( B_xydy_a_b_CSR , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! Initialise results from integrating a single triangle through the Voronoi cells
    single_row%n_max = 100
    single_row%n     = 0
    allocate( single_row%index_left( single_row%n_max))
    allocate( single_row%LI_xdy(     single_row%n_max))
    allocate( single_row%LI_mxydx(   single_row%n_max))
    allocate( single_row%LI_xydy(    single_row%n_max))

    ! == Trace all the line segments to fill the matrices
    ! ===================================================

    ti_hint = 1

    do vi = mesh_Vor%vi1, mesh_Vor%vi2 ! +1 because PETSc indexes from 0

      ! Clean up single row results
      single_row%n            = 0
      single_row%index_left   = 0
      single_row%LI_xdy       = 0
      single_row%LI_mxydx     = 0
      single_row%LI_xydy      = 0

      ! Integrate over the complete Voronoi cell boundary
      call calc_Voronoi_cell( mesh_Vor, vi, 0._dp, Vor, Vor_vi, Vor_ti, nVor)
      do vori1 = 1, nVor
        vori2 = vori1 + 1
        if (vori2 > nVor) vori2 = 1
        p = Vor( vori1,:)
        q = Vor( vori2,:)
        call trace_line_tri( mesh_tri, p, q, single_row, count_coincidences, ti_hint)
      end do

      ! Add the results for this triangle to the sparse matrix
      if (single_row%n == 0) then
        call add_empty_row_CSR_dist( B_xdy_a_b_CSR  , vi)
        call add_empty_row_CSR_dist( B_mxydx_a_b_CSR, vi)
        call add_empty_row_CSR_dist( B_xydy_a_b_CSR , vi)
      else
        do k = 1, single_row%n
          call add_entry_CSR_dist( B_xdy_a_b_CSR  , vi, single_row%index_left( k), single_row%LI_xdy(   k))
          call add_entry_CSR_dist( B_mxydx_a_b_CSR, vi, single_row%index_left( k), single_row%LI_mxydx( k))
          call add_entry_CSR_dist( B_xydy_a_b_CSR , vi, single_row%index_left( k), single_row%LI_xydy(  k))
        end do
      end if

    end do ! do vi = mesh_Vor%vi1, mesh_Vor%vi2
    call sync

    ! Convert matrices from Fortran to PETSc types
    call mat_CSR2petsc( B_xdy_a_b_CSR  , B_xdy_a_b  )
    call mat_CSR2petsc( B_mxydx_a_b_CSR, B_mxydx_a_b)
    call mat_CSR2petsc( B_xydy_a_b_CSR , B_xydy_a_b )

    ! Clean up the Fortran versions
    call deallocate_matrix_CSR_dist( B_xdy_a_b_CSR  )
    call deallocate_matrix_CSR_dist( B_mxydx_a_b_CSR)
    call deallocate_matrix_CSR_dist( B_xydy_a_b_CSR )

    ! Clean up after yourself
    deallocate( single_row%index_left )
    deallocate( single_row%LI_xdy     )
    deallocate( single_row%LI_mxydx   )
    deallocate( single_row%LI_xydy    )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine integrate_Voronoi_cells_through_triangles

end module create_maps_mesh_mesh
