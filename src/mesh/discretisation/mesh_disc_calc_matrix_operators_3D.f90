module mesh_disc_calc_matrix_operators_3D

  ! Routines for calculating 3-D matrix operators on the mesh.

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use mesh_types, only: type_mesh
  use CSR_matrix_basics, only: deallocate_matrix_CSR_dist, allocate_matrix_CSR_dist, &
    read_single_row_CSR_dist, add_entry_CSR_dist, finalise_matrix_CSR_dist

  implicit none

contains

subroutine calc_3D_matrix_operators_mesh( mesh, &
  dzeta_dx_ak, dzeta_dy_ak, dzeta_dx_bk, dzeta_dy_bk, &
  dzeta_dz_bk, dzeta_dz_bks, &
  d2zeta_dx2_bk, d2zeta_dxdy_bk, d2zeta_dy2_bk)
  ! Calculate all 3-D gradient operators in Cartesian coordinates

  ! In/output variables:
  type(type_mesh),                                    intent(inout) :: mesh
  real(dp), dimension(mesh%vi1:mesh%vi2,1:mesh%nz  ), intent(in   ) :: dzeta_dx_ak, dzeta_dy_ak
  real(dp), dimension(mesh%ti1:mesh%ti2,1:mesh%nz  ), intent(in   ) :: dzeta_dx_bk, dzeta_dy_bk
  real(dp), dimension(mesh%ti1:mesh%ti2,1:mesh%nz  ), intent(in   ) :: dzeta_dz_bk
  real(dp), dimension(mesh%ti1:mesh%ti2,1:mesh%nz-1), intent(in   ) :: dzeta_dz_bks
  real(dp), dimension(mesh%ti1:mesh%ti2,1:mesh%nz  ), intent(in   ) :: d2zeta_dx2_bk
  real(dp), dimension(mesh%ti1:mesh%ti2,1:mesh%nz  ), intent(in   ) :: d2zeta_dxdy_bk
  real(dp), dimension(mesh%ti1:mesh%ti2,1:mesh%nz  ), intent(in   ) :: d2zeta_dy2_bk

  ! Local variables:
  character(len=256), parameter :: routine_name = 'calc_3D_matrix_operators_mesh'

  ! Add routine to path
  call init_routine( routine_name)

  ! bk to ak (for calculating the horizontal stretch/shear strain rates in the BPA)
  call calc_3D_matrix_operators_mesh_bk_ak( mesh, dzeta_dx_ak, dzeta_dy_ak)

  ! ak to bk (for calculating the horizontal gradients of the effective viscosity in the BPA)
  call calc_3D_matrix_operators_mesh_ak_bk( mesh, dzeta_dx_bk, dzeta_dy_bk)

  ! bk to bks (for calculating the vertical shear strain rates in the BPA)
  call calc_3D_matrix_operators_mesh_bk_bks( mesh, dzeta_dz_bks)

  ! bks to bk (for calculating the vertical gradient of the effective viscosity in the BPA)
  call calc_3D_matrix_operators_mesh_bks_bk( mesh, dzeta_dz_bk)

  ! Map between the bks-grid and the ak-grid (for calculating strain rates in the BPA)
  call calc_3D_mapping_operator_mesh_bks_ak( mesh)
  call calc_3D_mapping_operator_mesh_ak_bks( mesh)

  ! bk to bk (for constructing the BPA stiffness matrix)
  call calc_3D_matrix_operators_mesh_bk_bk( mesh, dzeta_dx_bk, dzeta_dy_bk, dzeta_dz_bk, d2zeta_dx2_bk, d2zeta_dxdy_bk, d2zeta_dy2_bk)

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_3D_matrix_operators_mesh

subroutine calc_3D_matrix_operators_mesh_bk_ak( mesh, dzeta_dx_ak, dzeta_dy_ak)
  ! Calculate all 3-D gradient operators in Cartesian coordinates
  !
  ! The basic operators are defined in transformed coordinates [xh, yh, zeta], which
  ! are defined as:
  !
  !   xh   = x
  !   yh   = y
  !   zeta = (Hs - z) / Hi
  !
  ! Applying the chain rule to the gradient operators d/dx, d/dy yields:
  !
  !    d/dx = d/dxh + dzeta/dx d/dzeta
  !    d/dy = d/dyh + dzeta/dy d/dzeta
  !
  ! Theoretically, we could convert the basic d/dx, d/dy, d/dzeta matrix operators to PETSc format, perform
  ! matrix multiplications on them (so e.g. M_ddx = M_ddxh + D( dzeta/dx) M_ddzeta), and then convert
  ! the result back to CSR format. However, this is rather cumbersome to do (especially because
  ! the basic operators are defined in 2-D for the horizontal ones and in 1-D for the vertical,
  ! so we'd need to convert them to act on the 3-D mesh first).
  !
  ! Since the resulting operators  are relatively easy to interpret, we can just calculate
  ! their coefficients directly, which is done here.

  ! In/output variables:
  type(type_mesh),                                  intent(inout) :: mesh
  real(dp), dimension(mesh%vi1:mesh%vi2,1:mesh%nz), intent(in   ) :: dzeta_dx_ak, dzeta_dy_ak

  ! Local variables:
  character(len=256), parameter           :: routine_name = 'calc_3D_matrix_operators_mesh_bk_ak'
  integer                                 :: ncols, nrows, ncols_loc, nrows_loc, nnz_est_proc
  integer                                 :: vi
  integer,  dimension(:    ), allocatable :: single_row_vi_ind
  integer                                 :: single_row_vi_nnz
  real(dp), dimension(:    ), allocatable :: single_row_map_val
  real(dp), dimension(:    ), allocatable :: single_row_ddxh_val
  real(dp), dimension(:    ), allocatable :: single_row_ddyh_val
  integer                                 :: k
  integer,  dimension(:    ), allocatable :: single_row_k_ind
  integer                                 :: single_row_k_nnz
  real(dp), dimension(:    ), allocatable :: single_row_ddzeta_val
  integer                                 :: row_vik
  integer                                 :: ii,col_ti,ti,jj,kk,col_tikk
  real(dp)                                :: dzeta_dx, dzeta_dy
  real(dp)                                :: c_map, c_ddxh, c_ddyh, c_ddzeta
  real(dp)                                :: c_ddx, c_ddy

  ! Add routine to path
  call init_routine( routine_name)

  ! == Initialise the matrices using the native UFEMISM CSR-matrix format
  ! =====================================================================

  ! Deallocate existing matrices if necessary
  if (allocateD( mesh%M_ddx_bk_ak%ptr)) call deallocate_matrix_CSR_dist( mesh%M_ddx_bk_ak)
  if (allocateD( mesh%M_ddy_bk_ak%ptr)) call deallocate_matrix_CSR_dist( mesh%M_ddy_bk_ak)

  ! Matrix size
  ncols           = mesh%nTri     * mesh%nz ! from
  ncols_loc       = mesh%nTri_loc * mesh%nz
  nrows           = mesh%nV       * mesh%nz ! to
  nrows_loc       = mesh%nV_loc   * mesh%nz
  nnz_est_proc    = mesh%M_map_b_a%nnz * mesh%nz * 3

  call allocate_matrix_CSR_dist( mesh%M_ddx_bk_ak, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
  call allocate_matrix_CSR_dist( mesh%M_ddy_bk_ak, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

  ! Calculate shape functions and fill them into the matrices
  ! =========================================================

  ! allocate memory for single matrix rows
  allocate( single_row_vi_ind(   mesh%nC_mem*2))
  allocate( single_row_map_val(  mesh%nC_mem*2))
  allocate( single_row_ddxh_val( mesh%nC_mem*2))
  allocate( single_row_ddyh_val( mesh%nC_mem*2))

  allocate( single_row_k_ind(      mesh%nz))
  allocate( single_row_ddzeta_val( mesh%nz))

  ! Loop over all vertices
  do vi = mesh%vi1, mesh%vi2

    ! Read coefficients from the 2-D gradient operators for this triangle
    call read_single_row_CSR_dist( mesh%M_map_b_a, vi, single_row_vi_ind, single_row_map_val , single_row_vi_nnz)
    call read_single_row_CSR_dist( mesh%M_ddx_b_a, vi, single_row_vi_ind, single_row_ddxh_val, single_row_vi_nnz)
    call read_single_row_CSR_dist( mesh%M_ddy_b_a, vi, single_row_vi_ind, single_row_ddyh_val, single_row_vi_nnz)

    ! Loop over all layers
    do k = 1, mesh%nz

      ! Vertex vi, layer k corresponds to this matrix row
      row_vik = mesh%vik2n( vi,k)

      ! Read coefficients from the zeta gradient operators for this layer
      call read_single_row_CSR_dist( mesh%M_ddzeta_k_k_1D, k, single_row_k_ind, single_row_ddzeta_val, single_row_k_nnz)

      ! Gradients of zeta at vertex vi, layer k
      dzeta_dx = dzeta_dx_ak( vi,k)
      dzeta_dy = dzeta_dy_ak( vi,k)

      ! Loop over the entire 3-D local neighbourhood, calculate
      ! coefficients for all 3-D matrix operators

      do ii = 1, single_row_vi_nnz

        col_ti = single_row_vi_ind( ii)
        ti = mesh%n2ti( col_ti)

        ! Coefficients for horizontal gradient matrix operators
        c_map  = single_row_map_val(  ii)
        c_ddxh = single_row_ddxh_val( ii)
        c_ddyh = single_row_ddyh_val( ii)

        do jj = 1, single_row_k_nnz

          kk = single_row_k_ind( jj)

          ! Triangle ti, layer kk corresponds to this matrix row
          col_tikk = mesh%tik2n( ti,kk)

          ! Coefficients for vertical gradient matrix operators
          c_ddzeta = single_row_ddzeta_val( jj)

          ! Calculate coefficients
          c_ddx = 0._dp
          c_ddy = 0._dp

          ! Horizontal-only part
          if (kk == k) then
            c_ddx = c_ddx + c_ddxh                      ! Now:  d/dx   =  d/dxh ...
            c_ddy = c_ddy + c_ddyh                      ! Now:  d/dy   =  d/dyh ...
          end if ! if (kk == k) then

          ! Mixed part
          c_ddx = c_ddx + dzeta_dx * c_map * c_ddzeta   ! Now:  d/dx   =  d/dxh    +  dzeta/dx   d/dzeta
          c_ddy = c_ddy + dzeta_dy * c_map * c_ddzeta   ! Now:  d/dy   =  d/dyh    +  dzeta/dy   d/dzeta

          ! Add to CSR matrices
          call add_entry_CSR_dist( mesh%M_ddx_bk_ak, row_vik, col_tikk, c_ddx)
          call add_entry_CSR_dist( mesh%M_ddy_bk_ak, row_vik, col_tikk, c_ddy)

        end do ! do jj = 1, single_row_k_nnz
      end do ! do ii = 1, single_row_vi_nnz

    end do ! do k = 1, mesh%nz
  end do ! do vi = mesh%vi1, mesh%vi2

  ! Crop matrix memory
  call finalise_matrix_CSR_dist( mesh%M_ddx_bk_ak)
  call finalise_matrix_CSR_dist( mesh%M_ddy_bk_ak)

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_3D_matrix_operators_mesh_bk_ak

subroutine calc_3D_matrix_operators_mesh_ak_bk( mesh, dzeta_dx_bk, dzeta_dy_bk)
  ! Calculate all 3-D gradient operators in Cartesian coordinates
  !
  ! The basic operators are defined in transformed coordinates [xh, yh, zeta], which
  ! are defined as:
  !
  !   xh   = x
  !   yh   = y
  !   zeta = (Hs - z) / Hi
  !
  ! Applying the chain rule to the gradient operators d/dx, d/dy yields:
  !
  !    d/dx = d/dxh + dzeta/dx d/dzeta
  !    d/dy = d/dyh + dzeta/dy d/dzeta
  !
  ! Theoretically, we could convert the basic d/dx, d/dy, d/dzeta matrix operators to PETSc format, perform
  ! matrix multiplications on them (so e.g. M_ddx = M_ddxh + D( dzeta/dx) M_ddzeta), and then convert
  ! the result back to CSR format. However, this is rather cumbersome to do (especially because
  ! the basic operators are defined in 2-D for the horizontal ones and in 1-D for the vertical,
  ! so we'd need to convert them to act on the 3-D mesh first).
  !
  ! Since the resulting operators  are relatively easy to interpret, we can just calculate
  ! their coefficients directly, which is done here.

  ! In/output variables:
  type(type_mesh),                                  intent(inout) :: mesh
  real(dp), dimension(mesh%ti1:mesh%ti2,1:mesh%nz), intent(in   ) :: dzeta_dx_bk, dzeta_dy_bk

  ! Local variables:
  character(len=256), parameter           :: routine_name = 'calc_3D_matrix_operators_mesh_ak_bk'
  integer                                 :: ncols, nrows, ncols_loc, nrows_loc, nnz_est_proc
  integer                                 :: ti
  integer,  dimension(:    ), allocatable :: single_row_ti_ind
  integer                                 :: single_row_ti_nnz
  real(dp), dimension(:    ), allocatable :: single_row_map_val
  real(dp), dimension(:    ), allocatable :: single_row_ddxh_val
  real(dp), dimension(:    ), allocatable :: single_row_ddyh_val
  integer                                 :: k
  integer,  dimension(:    ), allocatable :: single_row_k_ind
  integer                                 :: single_row_k_nnz
  real(dp), dimension(:    ), allocatable :: single_row_ddzeta_val
  integer                                 :: row_tik
  integer                                 :: ii,col_vi,vi,jj,kk,col_vikk
  real(dp)                                :: dzeta_dx, dzeta_dy
  real(dp)                                :: c_map, c_ddxh, c_ddyh, c_ddzeta
  real(dp)                                :: c_ddx, c_ddy

  ! Add routine to path
  call init_routine( routine_name)

  ! == Initialise the matrices using the native UFEMISM CSR-matrix format
  ! =====================================================================

  ! Deallocate existing matrices if necessary
  if (allocateD( mesh%M_ddx_ak_bk%ptr)) call deallocate_matrix_CSR_dist( mesh%M_ddx_ak_bk)
  if (allocateD( mesh%M_ddy_ak_bk%ptr)) call deallocate_matrix_CSR_dist( mesh%M_ddy_ak_bk)

  ! Matrix size
  ncols           = mesh%nV       * mesh%nz ! from
  ncols_loc       = mesh%nV_loc   * mesh%nz
  nrows           = mesh%nTri     * mesh%nz ! to
  nrows_loc       = mesh%nTri_loc * mesh%nz
  nnz_est_proc    = mesh%M_map_a_b%nnz * mesh%nz * 3

  call allocate_matrix_CSR_dist( mesh%M_ddx_ak_bk, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
  call allocate_matrix_CSR_dist( mesh%M_ddy_ak_bk, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

  ! Calculate shape functions and fill them into the matrices
  ! =========================================================

  ! allocate memory for single matrix rows
  allocate( single_row_ti_ind(   mesh%nC_mem*2))
  allocate( single_row_map_val(  mesh%nC_mem*2))
  allocate( single_row_ddxh_val( mesh%nC_mem*2))
  allocate( single_row_ddyh_val( mesh%nC_mem*2))

  allocate( single_row_k_ind(      mesh%nz))
  allocate( single_row_ddzeta_val( mesh%nz))

  ! Loop over all triangles
  do ti = mesh%ti1, mesh%ti2

    ! Read coefficients from the 2-D gradient operators for this triangle
    call read_single_row_CSR_dist( mesh%M_map_a_b, ti, single_row_ti_ind, single_row_map_val , single_row_ti_nnz)
    call read_single_row_CSR_dist( mesh%M_ddx_a_b, ti, single_row_ti_ind, single_row_ddxh_val, single_row_ti_nnz)
    call read_single_row_CSR_dist( mesh%M_ddy_a_b, ti, single_row_ti_ind, single_row_ddyh_val, single_row_ti_nnz)

    ! Loop over all layers
    do k = 1, mesh%nz

      ! Triangle ti, layer k corresponds to this matrix row
      row_tik = mesh%tik2n( ti,k)

      ! Read coefficients from the zeta gradient operators for this layer
      call read_single_row_CSR_dist( mesh%M_ddzeta_k_k_1D, k, single_row_k_ind, single_row_ddzeta_val, single_row_k_nnz)

      ! Gradients of zeta at triangle ti, layer k
      dzeta_dx = dzeta_dx_bk( ti,k)
      dzeta_dy = dzeta_dy_bk( ti,k)

      ! Loop over the entire 3-D local neighbourhood, calculate
      ! coefficients for all 3-D matrix operators

      do ii = 1, single_row_ti_nnz

        col_vi = single_row_ti_ind( ii)
        vi = mesh%n2vi( col_vi)

        ! Coefficients for horizontal gradient matrix operators
        c_map  = single_row_map_val(  ii)
        c_ddxh = single_row_ddxh_val( ii)
        c_ddyh = single_row_ddyh_val( ii)

        do jj = 1, single_row_k_nnz

          kk = single_row_k_ind( jj)

          ! Vertex vi, layer kk corresponds to this matrix row
          col_vikk = mesh%vik2n( vi,kk)

          ! Coefficients for vertical gradient matrix operators
          c_ddzeta = single_row_ddzeta_val( jj)

          ! Calculate coefficients
          c_ddx = 0._dp
          c_ddy = 0._dp

          ! Horizontal-only part
          if (kk == k) then
            c_ddx = c_ddx + c_ddxh                      ! Now:  d/dx   =  d/dxh ...
            c_ddy = c_ddy + c_ddyh                      ! Now:  d/dy   =  d/dyh ...
          end if ! if (kk == k) then

          ! Mixed part
          c_ddx = c_ddx + dzeta_dx * c_map * c_ddzeta   ! Now:  d/dx   =  d/dxh    +  dzeta/dx   d/dzeta
          c_ddy = c_ddy + dzeta_dy * c_map * c_ddzeta   ! Now:  d/dy   =  d/dyh    +  dzeta/dy   d/dzeta

          ! Add to CSR matrices
          call add_entry_CSR_dist( mesh%M_ddx_ak_bk, row_tik, col_vikk, c_ddx)
          call add_entry_CSR_dist( mesh%M_ddy_ak_bk, row_tik, col_vikk, c_ddy)

        end do ! do jj = 1, single_row_k_nnz
      end do ! do ii = 1, single_row_ti_nnz

    end do ! do k = 1, mesh%nz
  end do ! do ti = mesh%ti1, mesh%ti2

  ! Crop matrix memory
  call finalise_matrix_CSR_dist( mesh%M_ddx_ak_bk)
  call finalise_matrix_CSR_dist( mesh%M_ddy_ak_bk)

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_3D_matrix_operators_mesh_ak_bk

subroutine calc_3D_matrix_operators_mesh_bk_bks( mesh, dzeta_dz_bks)
  ! Calculate all 3-D gradient operators in Cartesian coordinates
  !
  ! The basic operators are defined in transformed coordinates [xh, yh, zeta], which
  ! are defined as:
  !
  !   xh   = x
  !   yh   = y
  !   zeta = (Hs - z) / Hi
  !
  ! Applying the chain rule to the gradient operators d/dz, d2/dz2 yields:
  !
  !    d/dz   =  dzeta/dz     d/dzeta
  !   d2/dz2  = (dzeta/dz)^2 d2/dzeta2
  !
  ! Theoretically, we could convert all these basic matrix operators to PETSc format, perform
  ! matrix multiplications on them (so e.g. M_d2dxhdzeta = M_ddxh * M_ddzeta), and then convert
  ! the result back to CSR format. However, this is rather cumbersome to do (especially because
  ! the basic operators are defined in 2-D for the horizontal ones and in 1-D for the vertical,
  ! so we'd need to convert them to act on the 3-D mesh first).
  !
  ! Since the resulting operators  are relatively easy to interpret, we can just calculate
  ! their coefficients directly, which is done here.

  ! In/output variables:
  type(type_mesh),                                    intent(inout) :: mesh
  real(dp), dimension(mesh%ti1:mesh%ti2,1:mesh%nz-1), intent(in   ) :: dzeta_dz_bks

  ! Local variables:
  character(len=256), parameter           :: routine_name = 'calc_3D_matrix_operators_mesh_bk_bks'
  integer                                 :: ncols, nrows, ncols_loc, nrows_loc, nnz_per_row_est, nnz_est_proc
  integer                                 :: ti
  integer                                 :: ks
  integer,  dimension(:    ), allocatable :: single_row_ks_ind
  integer                                 :: single_row_ks_nnz
  real(dp), dimension(:    ), allocatable :: single_row_ddzeta_val
  integer                                 :: row_tiks
  integer                                 :: jj,k,col_tik
  real(dp)                                :: dzeta_dz
  real(dp)                                :: c_ddzeta
  real(dp)                                :: c_ddz

  ! Add routine to path
  call init_routine( routine_name)

  ! == Initialise the matrices using the native UFEMISM CSR-matrix format
  ! =====================================================================

  ! Deallocate existing matrices if necessary
  if (allocateD( mesh%M_ddz_bk_bks%ptr)) call deallocate_matrix_CSR_dist( mesh%M_ddz_bk_bks)

  ! Matrix size
  ncols           = mesh%nTri     *  mesh%nz    ! from
  ncols_loc       = mesh%nTri_loc *  mesh%nz
  nrows           = mesh%nTri     * (mesh%nz-1) ! to
  nrows_loc       = mesh%nTri_loc * (mesh%nz-1)
  nnz_per_row_est = 2
  nnz_est_proc    = nrows_loc * nnz_per_row_est

  call allocate_matrix_CSR_dist( mesh%M_ddz_bk_bks, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

  ! Calculate shape functions and fill them into the matrices
  ! =========================================================

  ! allocate memory for single matrix rows
  allocate( single_row_ks_ind(     mesh%nz))
  allocate( single_row_ddzeta_val( mesh%nz))

  ! Loop over all triangles
  do ti = mesh%ti1, mesh%ti2

    ! Loop over all staggered layers
    do ks = 1, mesh%nz-1

      ! Triangle ti, staggered layer ks corresponds to this matrix row
      row_tiks = mesh%tiks2n( ti,ks)

      ! Read coefficients from the zeta gradient operators for this staggered layer
      call read_single_row_CSR_dist( mesh%M_ddzeta_k_ks_1D, ks, single_row_ks_ind, single_row_ddzeta_val, single_row_ks_nnz)

      ! Gradients of zeta at triangle ti, staggered layer ks
      dzeta_dz = dzeta_dz_bks( ti,ks)

      ! Loop over the vertical local neighbourhood, calculate
      ! coefficients for all 3-D matrix operators

      do jj = 1, single_row_ks_nnz

        k = single_row_ks_ind( jj)

        ! Triangle tj, layer kk corresponds to this matrix row
        col_tik = mesh%tik2n( ti,k)

        ! Coefficients for vertical gradient matrix operators
        c_ddzeta = single_row_ddzeta_val( jj)

        ! Calculate coefficients
        c_ddz = dzeta_dz * c_ddzeta

        ! Add to CSR matrices
        call add_entry_CSR_dist( mesh%M_ddz_bk_bks, row_tiks, col_tik, c_ddz)

      end do ! do jj = 1, single_row_ks_nnz

    end do ! do ks = 1, mesh%nz-1
  end do ! do ti = mesh%ti1, mesh%ti2

  ! Clean up after yourself
  DEallocate( single_row_ks_ind)
  DEallocate( single_row_ddzeta_val)

  ! Crop matrix memory
  call finalise_matrix_CSR_dist( mesh%M_ddz_bk_bks)

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_3D_matrix_operators_mesh_bk_bks

subroutine calc_3D_matrix_operators_mesh_bks_bk( mesh, dzeta_dz_bk)
  ! Calculate all 3-D gradient operators in Cartesian coordinates
  !
  ! The basic operators are defined in transformed coordinates [xh, yh, zeta], which
  ! are defined as:
  !
  !   xh   = x
  !   yh   = y
  !   zeta = (Hs - z) / Hi
  !
  ! Applying the chain rule to the gradient operators d/dz, d2/dz2 yields:
  !
  !    d/dz   =  dzeta/dz     d/dzeta
  !   d2/dz2  = (dzeta/dz)^2 d2/dzeta2
  !
  ! Theoretically, we could convert all these basic matrix operators to PETSc format, perform
  ! matrix multiplications on them (so e.g. M_d2dxhdzeta = M_ddxh * M_ddzeta), and then convert
  ! the result back to CSR format. However, this is rather cumbersome to do (especially because
  ! the basic operators are defined in 2-D for the horizontal ones and in 1-D for the vertical,
  ! so we'd need to convert them to act on the 3-D mesh first).
  !
  ! Since the resulting operators  are relatively easy to interpret, we can just calculate
  ! their coefficients directly, which is done here.

  ! In/output variables:
  type(type_mesh),                                  intent(inout) :: mesh
  real(dp), dimension(mesh%ti1:mesh%ti2,1:mesh%nz), intent(in   ) :: dzeta_dz_bk

  ! Local variables:
  character(len=256), parameter           :: routine_name = 'calc_3D_matrix_operators_mesh_bks_bk'
  integer                                 :: ncols, nrows, ncols_loc, nrows_loc, nnz_per_row_est, nnz_est_proc
  integer                                 :: ti
  integer                                 :: k
  integer,  dimension(:    ), allocatable :: single_row_k_ind
  integer                                 :: single_row_k_nnz
  real(dp), dimension(:    ), allocatable :: single_row_map_val
  real(dp), dimension(:    ), allocatable :: single_row_ddzeta_val
  integer                                 :: row_tik
  integer                                 :: jj,ks,col_tiks
  real(dp)                                :: dzeta_dz
  real(dp)                                :: c_ddzeta
  real(dp)                                :: c_ddz, c_map

  ! Add routine to path
  call init_routine( routine_name)

  ! == Initialise the matrices using the native UFEMISM CSR-matrix format
  ! =====================================================================

  ! Deallocate existing matrices if necessary
  if (allocateD( mesh%M_map_bks_bk%ptr)) call deallocate_matrix_CSR_dist( mesh%M_map_bks_bk)
  if (allocateD( mesh%M_ddz_bks_bk%ptr)) call deallocate_matrix_CSR_dist( mesh%M_ddz_bks_bk)

  ! Matrix size
  ncols           = mesh%nTri     * (mesh%nz-1) ! from
  ncols_loc       = mesh%nTri_loc * (mesh%nz-1)
  nrows           = mesh%nTri     *  mesh%nz    ! to
  nrows_loc       = mesh%nTri_loc *  mesh%nz
  nnz_per_row_est = 2
  nnz_est_proc    = nrows_loc * nnz_per_row_est

  call allocate_matrix_CSR_dist( mesh%M_map_bks_bk, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
  call allocate_matrix_CSR_dist( mesh%M_ddz_bks_bk, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

  ! Calculate shape functions and fill them into the matrices
  ! =========================================================

  ! allocate memory for single matrix rows
  allocate( single_row_k_ind(      mesh%nz))
  allocate( single_row_map_val(    mesh%nz))
  allocate( single_row_ddzeta_val( mesh%nz))

  ! Loop over all triangles
  do ti = mesh%ti1, mesh%ti2

    ! Loop over all layers
    do k = 1, mesh%nz

      ! Triangle ti, layer k corresponds to this matrix row
      row_tik = mesh%tik2n( ti,k)

      ! Read coefficients from the zeta gradient operators for this staggered layer
      call read_single_row_CSR_dist( mesh%M_map_ks_k_1D   , k, single_row_k_ind, single_row_map_val   , single_row_k_nnz)
      call read_single_row_CSR_dist( mesh%M_ddzeta_ks_k_1D, k, single_row_k_ind, single_row_ddzeta_val, single_row_k_nnz)

      ! Gradients of zeta at triangle ti, layer k
      dzeta_dz = dzeta_dz_bk( ti,k)

      ! Loop over the vertical local neighbourhood, calculate
      ! coefficients for all 3-D matrix operators

      do jj = 1, single_row_k_nnz

        ks = single_row_k_ind( jj)

        ! Triangle tj, layer kk corresponds to this matrix row
        col_tiks = mesh%tiks2n( ti,ks)

        ! Coefficients for vertical gradient matrix operators
        c_map    = single_row_map_val(    jj)
        c_ddzeta = single_row_ddzeta_val( jj)

        ! Calculate coefficients
        c_ddz = dzeta_dz * c_ddzeta

        ! Add to CSR matrices
        call add_entry_CSR_dist( mesh%M_map_bks_bk, row_tik, col_tiks, c_map)
        call add_entry_CSR_dist( mesh%M_ddz_bks_bk, row_tik, col_tiks, c_ddz)

      end do ! do jj = 1, single_row_k_nnz

    end do ! do ks = 1, mesh%nz
  end do ! do ti = mesh%ti1, mesh%ti2

  ! Crop matrix memory
  call finalise_matrix_CSR_dist( mesh%M_map_bks_bk)
  call finalise_matrix_CSR_dist( mesh%M_ddz_bks_bk)

  ! Clean up after yourself
  DEallocate( single_row_k_ind)
  DEallocate( single_row_map_val)
  DEallocate( single_row_ddzeta_val)

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_3D_matrix_operators_mesh_bks_bk

subroutine calc_3D_mapping_operator_mesh_bks_ak( mesh)
  ! Calculate mapping operator from the bks-grid to the ak-grid

  ! In/output variables:
  type(type_mesh), intent(inout) :: mesh

  ! Local variables:
  character(len=256), parameter           :: routine_name = 'calc_3D_mapping_operator_mesh_bks_ak'
  integer                                 :: ncols, nrows, ncols_loc, nrows_loc, nnz_est_proc
  integer                                 :: vi
  integer,  dimension(:    ), allocatable :: single_row_vi_ind
  integer                                 :: single_row_vi_nnz
  real(dp), dimension(:    ), allocatable :: single_row_map_b_a_val
  integer                                 :: k
  integer,  dimension(:    ), allocatable :: single_row_k_ind
  integer                                 :: single_row_k_nnz
  real(dp), dimension(:    ), allocatable :: single_row_map_ks_k_val
  integer                                 :: row_vik
  integer                                 :: ii,col_ti,ti,jj,ks,col_tiks
  real(dp)                                :: c_map_b_a, c_map_ks_k
  real(dp)                                :: c_map

  ! Add routine to path
  call init_routine( routine_name)

  ! == Initialise the matrices using the native UFEMISM CSR-matrix format
  ! =====================================================================

  ! Deallocate existing matrices if necessary
  if (allocateD( mesh%M_map_bks_ak%ptr)) call deallocate_matrix_CSR_dist( mesh%M_map_bks_ak)

  ! Matrix size
  ncols           = mesh%nTri     * (mesh%nz-1) ! from
  ncols_loc       = mesh%nTri_loc * (mesh%nz-1)
  nrows           = mesh%nV       *  mesh%nz    ! to
  nrows_loc       = mesh%nV_loc   *  mesh%nz
  nnz_est_proc    = mesh%M_map_b_a%nnz * mesh%nz * 2

  call allocate_matrix_CSR_dist( mesh%M_map_bks_ak, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

  ! Calculate shape functions and fill them into the matrices
  ! =========================================================

  ! allocate memory for single matrix rows
  allocate( single_row_vi_ind(      mesh%nC_mem*2))
  allocate( single_row_map_b_a_val( mesh%nC_mem*2))

  allocate( single_row_k_ind(        mesh%nz      ))
  allocate( single_row_map_ks_k_val( mesh%nz      ))

  ! Loop over all vertices
  do vi = mesh%vi1, mesh%vi2

    ! Read coefficients from the 2-D gradient operators for this vertex
    call read_single_row_CSR_dist( mesh%M_map_b_a, vi, single_row_vi_ind, single_row_map_b_a_val, single_row_vi_nnz)

    ! Loop over all layers
    do k = 1, mesh%nz

      ! Vertex vi, layer k corresponds to this matrix row
      row_vik = mesh%vik2n( vi,k)

      ! Read coefficients from the zeta gradient operators for this layer
      call read_single_row_CSR_dist( mesh%M_map_ks_k_1D, k, single_row_k_ind, single_row_map_ks_k_val, single_row_k_nnz)

      ! Loop over the entire 3-D local neighbourhood, calculate
      ! coefficients for all 3-D matrix operators

      do ii = 1, single_row_vi_nnz

        col_ti = single_row_vi_ind( ii)
        ti = mesh%n2ti( col_ti)

        ! Coefficients for horizontal gradient matrix operators
        c_map_b_a = single_row_map_b_a_val( ii)

        do jj = 1, single_row_k_nnz

          ks = single_row_k_ind( jj)

          ! Triangle ti, staggered layer ks corresponds to this matrix row
          col_tiks = mesh%tiks2n( ti,ks)

          ! Coefficients for vertical gradient matrix operators
          c_map_ks_k = single_row_map_ks_k_val( jj)

          ! Calculate coefficients
          c_map = c_map_b_a * c_map_ks_k

          ! Add to CSR matrices
          call add_entry_CSR_dist( mesh%M_map_bks_ak, row_vik, col_tiks, c_map)

        end do ! do jj = 1, single_row_k_nnz
      end do ! do ii = 1, single_row_vi_nnz

    end do ! do k = 1, mesh%nz
  end do ! do vi = mesh%vi1, mesh%vi2

  ! Crop matrix memory
  call finalise_matrix_CSR_dist( mesh%M_map_bks_ak)

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_3D_mapping_operator_mesh_bks_ak

subroutine calc_3D_mapping_operator_mesh_ak_bks( mesh)
  ! Calculate mapping operator from the ak-grid to the bks-grid

  ! In/output variables:
  type(type_mesh), intent(inout) :: mesh

  ! Local variables:
  character(len=256), parameter           :: routine_name = 'calc_3D_mapping_operator_mesh_ak_bks'
  integer                                 :: ncols, nrows, ncols_loc, nrows_loc, nnz_est_proc
  integer                                 :: ti
  integer,  dimension(:    ), allocatable :: single_row_ti_ind
  integer                                 :: single_row_ti_nnz
  real(dp), dimension(:    ), allocatable :: single_row_map_a_b_val
  integer                                 :: ks
  integer,  dimension(:    ), allocatable :: single_row_ks_ind
  integer                                 :: single_row_ks_nnz
  real(dp), dimension(:    ), allocatable :: single_row_map_k_ks_val
  integer                                 :: row_tiks
  integer                                 :: ii,col_vi,vi,jj,k,col_vik
  real(dp)                                :: c_map_a_b, c_map_k_ks
  real(dp)                                :: c_map

  ! Add routine to path
  call init_routine( routine_name)

  ! == Initialise the matrices using the native UFEMISM CSR-matrix format
  ! =====================================================================

  ! Deallocate existing matrices if necessary
  if (allocateD( mesh%M_map_ak_bks%ptr)) call deallocate_matrix_CSR_dist( mesh%M_map_ak_bks)

  ! Matrix size
  ncols           = mesh%nV       *  mesh%nz    ! from
  ncols_loc       = mesh%nV_loc   *  mesh%nz
  nrows           = mesh%nTri     * (mesh%nz-1) ! to
  nrows_loc       = mesh%nTri_loc * (mesh%nz-1)
  nnz_est_proc    = mesh%M_map_a_b%nnz * (mesh%nz-1) * 2

  call allocate_matrix_CSR_dist( mesh%M_map_ak_bks, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

  ! Calculate shape functions and fill them into the matrices
  ! =========================================================

  ! allocate memory for single matrix rows
  allocate( single_row_ti_ind(      mesh%nC_mem*2))
  allocate( single_row_map_a_b_val( mesh%nC_mem*2))

  allocate( single_row_ks_ind(       mesh%nz      ))
  allocate( single_row_map_k_ks_val( mesh%nz      ))

  ! Loop over all triangles
  do ti = mesh%ti1, mesh%ti2

    ! Read coefficients from the 2-D gradient operators for this triangle
    call read_single_row_CSR_dist( mesh%M_map_a_b, ti, single_row_ti_ind, single_row_map_a_b_val, single_row_ti_nnz)

    ! Loop over all staggered layers
    do ks = 1, mesh%nz-1

      ! Triangle ti, staggered layer ks corresponds to this matrix row
      row_tiks = mesh%tiks2n( ti,ks)

      ! Read coefficients from the zeta gradient operators for this staggered layer
      call read_single_row_CSR_dist( mesh%M_map_k_ks_1D, ks, single_row_ks_ind, single_row_map_k_ks_val, single_row_ks_nnz)

      ! Loop over the entire 3-D local neighbourhood, calculate
      ! coefficients for all 3-D matrix operators

      do ii = 1, single_row_ti_nnz

        col_vi = single_row_ti_ind( ii)
        vi = mesh%n2vi( col_vi)

        ! Coefficients for horizontal gradient matrix operators
        c_map_a_b = single_row_map_a_b_val( ii)

        do jj = 1, single_row_ks_nnz

          k = single_row_ks_ind( jj)

          ! Vertex vi, layer k corresponds to this matrix row
          col_vik = mesh%vik2n( vi,k)

          ! Coefficients for vertical gradient matrix operators
          c_map_k_ks = single_row_map_k_ks_val( jj)

          ! Calculate coefficients
          c_map = c_map_a_b * c_map_k_ks

          ! Add to CSR matrices
          call add_entry_CSR_dist( mesh%M_map_ak_bks, row_tiks, col_vik, c_map)

        end do ! do jj = 1, single_row_ks_nnz
      end do ! do ii = 1, single_row_ti_nnz

    end do ! do ks = 1, mesh%nz-1
  end do !     do ti = mesh%ti1, mesh%ti2

  ! Crop matrix memory
  call finalise_matrix_CSR_dist( mesh%M_map_ak_bks)

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_3D_mapping_operator_mesh_ak_bks

subroutine calc_3D_matrix_operators_mesh_bk_bk( mesh, &
  dzeta_dx_bk, dzeta_dy_bk, dzeta_dz_bk, d2zeta_dx2_bk, d2zeta_dxdy_bk, d2zeta_dy2_bk)
  ! Calculate all 3-D gradient operators in Cartesian coordinates
  !
  ! The basic operators are defined in transformed coordinates [xh, yh, zeta], which
  ! are defined as:
  !
  !   xh   = x
  !   yh   = y
  !   zeta = (Hs - z) / Hi
  !
  ! Applying the chain rule to the gradient operators d/dx, d/dy, d/dz, d2/dx2, d2/dxdy,
  ! d2/dy2, d2/dz2 yields:
  !
  !    d/dx   =  d/dxh    +  dzeta/dx   d/dzeta
  !    d/dy   =  d/dyh    +  dzeta/dy   d/dzeta
  !    d/dz   =              dzeta/dz   d/dzeta
  !   d2/dx2  = d2/dxh2   + d2zeta/dx2  d/dzeta + (dzeta/dx)^2       d2/dzeta2 + 2 dzeta/dx d2/dxhdzeta
  !   d2/dxdy = d2/dxhdyh + d2zeta/dxdy d/dzeta +  dzeta/dx dzeta/dy d2/dzeta2 +   dzeta/dx d2/dyhdzeta + dzeta/dy dxhdzeta
  !   d2/dy2  = d2/dyh2   + d2zeta/dy2  d/dzeta + (dzeta/dy)^2       d2/dzeta2 + 2 dzeta/dy d2/dyhdzeta
  !   d2/dz2  =                                   (dzeta/dz)^2       d2/dzeta2
  !
  ! The d/dxh, d/dyh, d2/dxh2, d2/dxhdyh, d2/dyh2 operators all act in the horizontal plane, so
  ! calculating them requires only information from neighbouring triangles in the same horizontal
  ! layer k. The d/dzeta, d2/dzeta2 operators only act in the vertical column, so they only
  ! require information from the two adjacent horizontal layers k-1, k+1 at the same triangle ti.
  ! Only the two "mixed" gradient operators, d2/dxhdzeta, d2/dyhdzeta, require information from
  ! all neighbouring triangles in both adjacent layers.
  !
  ! Theoretically, we could convert all these basic matrix operators to PETSc format, perform
  ! matrix multiplications on them (so e.g. M_d2dxhdzeta = M_ddxh * M_ddzeta), and then convert
  ! the result back to CSR format. However, this is rather cumbersome to do (especially because
  ! the basic operators are defined in 2-D for the horizontal ones and in 1-D for the vertical,
  ! so we'd need to convert them to act on the 3-D mesh first).
  !
  ! Since the resulting operators  are relatively easy to interpret, we can just calculate
  ! their coefficients directly, which is done here.

  ! In/output variables:
  type(type_mesh),                                  intent(inout) :: mesh
  real(dp), dimension(mesh%ti1:mesh%ti2,1:mesh%nz), intent(in   ) :: dzeta_dx_bk
  real(dp), dimension(mesh%ti1:mesh%ti2,1:mesh%nz), intent(in   ) :: dzeta_dy_bk
  real(dp), dimension(mesh%ti1:mesh%ti2,1:mesh%nz), intent(in   ) :: dzeta_dz_bk
  real(dp), dimension(mesh%ti1:mesh%ti2,1:mesh%nz), intent(in   ) :: d2zeta_dx2_bk
  real(dp), dimension(mesh%ti1:mesh%ti2,1:mesh%nz), intent(in   ) :: d2zeta_dxdy_bk
  real(dp), dimension(mesh%ti1:mesh%ti2,1:mesh%nz), intent(in   ) :: d2zeta_dy2_bk

  ! Local variables:
  character(len=256), parameter           :: routine_name = 'calc_3D_matrix_operators_mesh_bk_bk'
  integer                                 :: ncols, nrows, ncols_loc, nrows_loc, nnz_est_proc
  integer                                 :: ti
  integer,  dimension(:    ), allocatable :: single_row_ti_ind
  integer                                 :: single_row_ti_nnz
  real(dp), dimension(:    ), allocatable :: single_row_ddxh_val
  real(dp), dimension(:    ), allocatable :: single_row_ddyh_val
  real(dp), dimension(:    ), allocatable :: single_row_d2dxh2_val
  real(dp), dimension(:    ), allocatable :: single_row_d2dxhdyh_val
  real(dp), dimension(:    ), allocatable :: single_row_d2dyh2_val
  integer                                 :: k
  integer,  dimension(:    ), allocatable :: single_row_k_ind
  integer                                 :: single_row_k_nnz
  real(dp), dimension(:    ), allocatable :: single_row_ddzeta_val
  real(dp), dimension(:    ), allocatable :: single_row_d2dzeta2_val
  integer                                 :: row_tik
  integer                                 :: ii,col_tj,tj,jj,kk,col_tjkk
  real(dp)                                :: dzeta_dx, dzeta_dy, dzeta_dz, d2zeta_dx2, d2zeta_dxdy, d2zeta_dy2
  real(dp)                                :: c_ddxh, c_ddyh, c_d2dxh2, c_d2dxhdyh, c_d2dyh2, c_ddzeta, c_d2dzeta2
  real(dp)                                :: c_ddx, c_ddy, c_ddz, c_d2dx2, c_d2dxdy, c_d2dy2, c_d2dz2

  ! Add routine to path
  call init_routine( routine_name)

  ! == Initialise the matrices using the native UFEMISM CSR-matrix format
  ! =====================================================================

  ! Deallocate existing matrices if necessary
  if (allocateD( mesh%M2_ddx_bk_bk%ptr   )) call deallocate_matrix_CSR_dist( mesh%M2_ddx_bk_bk)
  if (allocateD( mesh%M2_ddy_bk_bk%ptr   )) call deallocate_matrix_CSR_dist( mesh%M2_ddy_bk_bk)
  if (allocateD( mesh%M2_ddz_bk_bk%ptr   )) call deallocate_matrix_CSR_dist( mesh%M2_ddz_bk_bk)
  if (allocateD( mesh%M2_d2dx2_bk_bk%ptr )) call deallocate_matrix_CSR_dist( mesh%M2_d2dx2_bk_bk)
  if (allocateD( mesh%M2_d2dxdy_bk_bk%ptr)) call deallocate_matrix_CSR_dist( mesh%M2_d2dxdy_bk_bk)
  if (allocateD( mesh%M2_d2dy2_bk_bk%ptr )) call deallocate_matrix_CSR_dist( mesh%M2_d2dy2_bk_bk)
  if (allocateD( mesh%M2_d2dz2_bk_bk%ptr )) call deallocate_matrix_CSR_dist( mesh%M2_d2dz2_bk_bk)

  ! Matrix size
  ncols           = mesh%nTri     * mesh%nz ! from
  ncols_loc       = mesh%nTri_loc * mesh%nz
  nrows           = mesh%nTri     * mesh%nz ! to
  nrows_loc       = mesh%nTri_loc * mesh%nz
  nnz_est_proc    = mesh%M2_ddx_b_b%nnz * mesh%nz * 3

  call allocate_matrix_CSR_dist( mesh%M2_ddx_bk_bk   , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
  call allocate_matrix_CSR_dist( mesh%M2_ddy_bk_bk   , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
  call allocate_matrix_CSR_dist( mesh%M2_ddz_bk_bk   , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
  call allocate_matrix_CSR_dist( mesh%M2_d2dx2_bk_bk , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
  call allocate_matrix_CSR_dist( mesh%M2_d2dxdy_bk_bk, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
  call allocate_matrix_CSR_dist( mesh%M2_d2dy2_bk_bk , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
  call allocate_matrix_CSR_dist( mesh%M2_d2dz2_bk_bk , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

  ! Calculate shape functions and fill them into the matrices
  ! =========================================================

  ! allocate memory for single matrix rows
  allocate( single_row_ti_ind(       mesh%nC_mem*2))
  allocate( single_row_ddxh_val(     mesh%nC_mem*2))
  allocate( single_row_ddyh_val(     mesh%nC_mem*2))
  allocate( single_row_d2dxh2_val(   mesh%nC_mem*2))
  allocate( single_row_d2dxhdyh_val( mesh%nC_mem*2))
  allocate( single_row_d2dyh2_val(   mesh%nC_mem*2))

  allocate( single_row_k_ind(        mesh%nz      ))
  allocate( single_row_ddzeta_val(   mesh%nz      ))
  allocate( single_row_d2dzeta2_val( mesh%nz      ))

  ! Loop over all triangles
  do ti = mesh%ti1, mesh%ti2

    ! Read coefficients from the 2-D gradient operators for this triangle
    call read_single_row_CSR_dist( mesh%M2_ddx_b_b   , ti, single_row_ti_ind, single_row_ddxh_val    , single_row_ti_nnz)
    call read_single_row_CSR_dist( mesh%M2_ddy_b_b   , ti, single_row_ti_ind, single_row_ddyh_val    , single_row_ti_nnz)
    call read_single_row_CSR_dist( mesh%M2_d2dx2_b_b , ti, single_row_ti_ind, single_row_d2dxh2_val  , single_row_ti_nnz)
    call read_single_row_CSR_dist( mesh%M2_d2dxdy_b_b, ti, single_row_ti_ind, single_row_d2dxhdyh_val, single_row_ti_nnz)
    call read_single_row_CSR_dist( mesh%M2_d2dy2_b_b , ti, single_row_ti_ind, single_row_d2dyh2_val  , single_row_ti_nnz)

    ! Loop over all layers
    do k = 1, mesh%nz

      ! Triangle ti, layer k corresponds to this matrix row
      row_tik = mesh%tik2n( ti,k)

      ! Read coefficients from the zeta gradient operators for this layer
      call read_single_row_CSR_dist( mesh%M_ddzeta_k_k_1D  , k, single_row_k_ind, single_row_ddzeta_val  , single_row_k_nnz)
      call read_single_row_CSR_dist( mesh%M_d2dzeta2_k_k_1D, k, single_row_k_ind, single_row_d2dzeta2_val, single_row_k_nnz)

      ! Gradients of zeta at triangle ti, layer k
      dzeta_dx    = dzeta_dx_bk(    ti,k)
      dzeta_dy    = dzeta_dy_bk(    ti,k)
      dzeta_dz    = dzeta_dz_bk(    ti,k)
      d2zeta_dx2  = d2zeta_dx2_bk(  ti,k)
      d2zeta_dxdy = d2zeta_dxdy_bk( ti,k)
      d2zeta_dy2  = d2zeta_dy2_bk(  ti,k)

      ! Loop over the entire 3-D local neighbourhood, calculate
      ! coefficients for all 3-D matrix operators

      do ii = 1, single_row_ti_nnz

        col_tj = single_row_ti_ind( ii)
        tj = mesh%n2ti( col_tj)

        ! Coefficients for horizontal gradient matrix operators
        c_ddxh     = single_row_ddxh_val(     ii)
        c_ddyh     = single_row_ddyh_val(     ii)
        c_d2dxh2   = single_row_d2dxh2_val(   ii)
        c_d2dxhdyh = single_row_d2dxhdyh_val( ii)
        c_d2dyh2   = single_row_d2dyh2_val(   ii)

        do jj = 1, single_row_k_nnz

          kk = single_row_k_ind( jj)

          ! Triangle tj, layer kk corresponds to this matrix row
          col_tjkk = mesh%tik2n( tj,kk)

          ! Coefficients for vertical gradient matrix operators
          c_ddzeta   = single_row_ddzeta_val(   jj)
          c_d2dzeta2 = single_row_d2dzeta2_val( jj)

          ! Calculate coefficients
          c_ddx    = 0._dp
          c_ddy    = 0._dp
          c_ddz    = 0._dp
          c_d2dx2  = 0._dp
          c_d2dxdy = 0._dp
          c_d2dy2  = 0._dp
          c_d2dz2  = 0._dp

          ! Horizontal-only part
          if (kk == k) then
            c_ddx    = c_ddx    + c_ddxh                                                            ! Now:  d/dx   =  d/dxh ...
            c_ddy    = c_ddy    + c_ddyh                                                            ! Now:  d/dy   =  d/dyh ...
            c_d2dx2  = c_d2dx2  + c_d2dxh2                                                          ! Now: d2/dx2  = d2/dxh2 ...
            c_d2dxdy = c_d2dxdy + c_d2dxhdyh                                                        ! Now: d2/dxdy = d2/dxhdyh ...
            c_d2dy2  = c_d2dy2  + c_d2dyh2                                                          ! Now: d2/dy2  = d2/dyh2 ...
          end if ! if (kk == k) then

          ! Vertical-only part
          if (tj == ti) then
            c_ddx    = c_ddx    + dzeta_dx    * c_ddzeta                                            ! Now:  d/dx   =  d/dxh    +  dzeta/dx   d/dzeta
            c_ddy    = c_ddy    + dzeta_dy    * c_ddzeta                                            ! Now:  d/dy   =  d/dyh    +  dzeta/dy   d/dzeta
            c_ddz    = c_ddz    + dzeta_dz    * c_ddzeta                                            ! Now:  d/dz   =              dzeta/dz   d/dzeta
            c_d2dx2  = c_d2dx2  + d2zeta_dx2  * c_ddzeta + dzeta_dx * dzeta_dx * c_d2dzeta2         ! Now: d2/dx2  = d2/dxh2   + d2zeta/dx2  d/dzeta + (dzeta/dx)^2       d2/dzeta2 + ...
            c_d2dxdy = c_d2dxdy + d2zeta_dxdy * c_ddzeta + dzeta_dx * dzeta_dy * c_d2dzeta2         ! Now: d2/dxdy = d2/dxhdyh + d2zeta/dxdy d/dzeta +  dzeta/dx dzeta/dy d2/dzeta2 + ...
            c_d2dy2  = c_d2dy2  + d2zeta_dy2  * c_ddzeta + dzeta_dy * dzeta_dy * c_d2dzeta2         ! Now: d2/dy2  = d2/dyh2   + d2zeta/dy2  d/dzeta + (dzeta/dy)^2       d2/dzeta2 + ...
            c_d2dz2  = c_d2dz2  + dzeta_dz**2 * c_d2dzeta2                                          ! Now: d2/dz2  =                                   (dzeta/dz)^2       d2/dzeta2
          end if ! if (tj == ti) then

          ! Mixed part
          c_d2dx2  = c_d2dx2  + 2._dp * dzeta_dx * c_ddxh * c_ddzeta                                ! Now: d2/dx2  = d2/dxh2   + d2zeta/dx2  d/dzeta + (dzeta/dx)^2       d2/dzeta2 + 2 dzeta/dx d2/dxhdzeta
          c_d2dxdy = c_d2dxdy +         dzeta_dx * c_ddyh * c_ddzeta + dzeta_dy * c_ddxh * c_ddzeta ! Now: d2/dxdy = d2/dxhdyh + d2zeta/dxdy d/dzeta +  dzeta/dx dzeta/dy d2/dzeta2 +   dzeta/dx d2/dyhdzeta + dzeta/dy dxhdzeta
          c_d2dy2  = c_d2dy2  + 2._dp * dzeta_dy * c_ddyh * c_ddzeta                                ! Now: d2/dy2  = d2/dyh2   + d2zeta/dy2  d/dzeta + (dzeta/dy)^2       d2/dzeta2 + 2 dzeta/dy d2/dyhdzeta

          ! Add to CSR matrices
          call add_entry_CSR_dist( mesh%M2_ddx_bk_bk   , row_tik, col_tjkk, c_ddx   )
          call add_entry_CSR_dist( mesh%M2_ddy_bk_bk   , row_tik, col_tjkk, c_ddy   )
          call add_entry_CSR_dist( mesh%M2_ddz_bk_bk   , row_tik, col_tjkk, c_ddz   )
          call add_entry_CSR_dist( mesh%M2_d2dx2_bk_bk , row_tik, col_tjkk, c_d2dx2 )
          call add_entry_CSR_dist( mesh%M2_d2dxdy_bk_bk, row_tik, col_tjkk, c_d2dxdy)
          call add_entry_CSR_dist( mesh%M2_d2dy2_bk_bk , row_tik, col_tjkk, c_d2dy2 )
          call add_entry_CSR_dist( mesh%M2_d2dz2_bk_bk , row_tik, col_tjkk, c_d2dz2 )

        end do ! do jj = 1, single_row_k_nnz
      end do ! do ii = 1, single_row_ti_nnz

    end do ! do k = 1, mesh%nz
  end do ! do ti = mesh%ti1, mesh%ti2

  ! Crop matrix memory
  call finalise_matrix_CSR_dist( mesh%M2_ddx_bk_bk)
  call finalise_matrix_CSR_dist( mesh%M2_ddy_bk_bk)
  call finalise_matrix_CSR_dist( mesh%M2_ddz_bk_bk)
  call finalise_matrix_CSR_dist( mesh%M2_d2dx2_bk_bk)
  call finalise_matrix_CSR_dist( mesh%M2_d2dxdy_bk_bk)
  call finalise_matrix_CSR_dist( mesh%M2_d2dy2_bk_bk)
  call finalise_matrix_CSR_dist( mesh%M2_d2dz2_bk_bk)

  ! Clean up after yourself
  DEallocate( single_row_ti_ind)
  DEallocate( single_row_ddxh_val)
  DEallocate( single_row_ddyh_val)
  DEallocate( single_row_d2dxh2_val)
  DEallocate( single_row_d2dxhdyh_val)
  DEallocate( single_row_d2dyh2_val)
  DEallocate( single_row_k_ind)
  DEallocate( single_row_ddzeta_val)
  DEallocate( single_row_d2dzeta2_val)

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_3D_matrix_operators_mesh_bk_bk

end module mesh_disc_calc_matrix_operators_3D
