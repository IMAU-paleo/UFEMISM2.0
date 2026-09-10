submodule(momentum_balance_solver_BPA) solve_linearised_BPA

contains

  subroutine solve_BPA_linearised( self, ice, n_Axb_its, &
    BC_prescr_mask_bk, BC_prescr_u_bk, BC_prescr_v_bk)

    ! In/output variables:
    class(type_momentum_balance_solver_BPA),                         intent(inout) :: self
    class(atype_ice_model_data),                                     intent(in   ) :: ice
    integer,                                                         intent(  out) :: n_Axb_its              ! Number of iterations used in the iterative solver
    integer,  dimension(self%mesh%ti1:self%mesh%ti2,1:self%mesh%nz), intent(in   ) :: BC_prescr_mask_bk      ! Mask of triangles where velocity is prescribed
    real(dp), dimension(self%mesh%ti1:self%mesh%ti2,1:self%mesh%nz), intent(in   ) :: BC_prescr_u_bk         ! Prescribed velocities in the x-direction
    real(dp), dimension(self%mesh%ti1:self%mesh%ti2,1:self%mesh%nz), intent(in   ) :: BC_prescr_v_bk         ! Prescribed velocities in the y-direction

    ! Local variables:
    character(len=*), parameter         :: routine_name = 'solve_BPA_linearised'
    type(type_CSR_matrix_dp)            :: A_CSR
    real(dp), dimension(:), allocatable :: bb
    real(dp), dimension(:), allocatable :: uv_bkuv
    integer                             :: row_tikuv,ti,k

    ! Add routine to path
    call init_routine( routine_name)

    ! Store the previous solution
    call gather_to_all( self%u_bk, self%u_bk_prev)
    call gather_to_all( self%v_bk, self%v_bk_prev)

    ! Assemble the matrix equation
    call self%assemble_BPA_linearised_matrix_eq( ice, &
      BC_prescr_mask_bk, BC_prescr_u_bk, BC_prescr_v_bk, &
      A_CSR, bb, uv_bkuv)

    ! Use PETSc to solve the matrix equation
    call solve_matrix_equation_CSR_PETSc( A_CSR, bb, uv_bkuv, &
      self%PETSc_rtol, self%PETSc_abstol, n_Axb_its, &
      PETSc_KSPtype = C%stress_balance_PETSc_KSPtype, PETSc_PCtype = C%stress_balance_PETSc_PCtype)

    ! Disentangle the u and v components of the velocity solution
    do ti = self%mesh%ti1, self%mesh%ti2
    do k  = 1, self%mesh%nz

      ! u
      row_tikuv = self%mesh%tikuv2n( ti,k,1)
      self%u_bk( ti,k) = uv_bkuv( row_tikuv)

      ! v
      row_tikuv = self%mesh%tikuv2n( ti,k,2)
      self%v_bk( ti,k) = uv_bkuv( row_tikuv)

    end do
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine solve_BPA_linearised

  subroutine assemble_BPA_linearised_matrix_eq( self, ice, &
    BC_prescr_mask_bk, BC_prescr_u_bk, BC_prescr_v_bk, &
    A_CSR, bb, uv_bkuv)

    ! In/output variables:
    class(type_momentum_balance_solver_BPA),                         intent(in   ) :: self
    class(atype_ice_model_data),                                     intent(in   ) :: ice
    integer,  dimension(self%mesh%ti1:self%mesh%ti2,1:self%mesh%nz), intent(in   ) :: BC_prescr_mask_bk      ! Mask of triangles where velocity is prescribed
    real(dp), dimension(self%mesh%ti1:self%mesh%ti2,1:self%mesh%nz), intent(in   ) :: BC_prescr_u_bk         ! Prescribed velocities in the x-direction
    real(dp), dimension(self%mesh%ti1:self%mesh%ti2,1:self%mesh%nz), intent(in   ) :: BC_prescr_v_bk         ! Prescribed velocities in the y-direction
    type(type_CSR_matrix_dp),                                        intent(inout) :: A_CSR
    real(dp), dimension(:), allocatable,                             intent(inout) :: bb
    real(dp), dimension(:), allocatable,                             intent(inout) :: uv_bkuv

    ! Local variables:
    character(len=*), parameter :: routine_name = 'assemble_BPA_linearised_matrix_eq'
    integer                     :: ncols, ncols_loc, nrows, nrows_loc, nnz_est_proc
    integer                     :: row_tikuv,ti,k,uv

    ! Add routine to path
    call init_routine( routine_name)

    ! == Initialise the stiffness matrix using the native UFEMISM CSR-matrix format
    ! =============================================================================

    ! Matrix size
    ncols           = self%mesh%nTri     * self%mesh%nz * 2      ! from
    ncols_loc       = self%mesh%nTri_loc * self%mesh%nz * 2
    nrows           = self%mesh%nTri     * self%mesh%nz * 2      ! to
    nrows_loc       = self%mesh%nTri_loc * self%mesh%nz * 2
    nnz_est_proc    = self%mesh%M2_ddx_bk_bk%nnz   * 4

    call A_CSR%allocate( nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! allocate memory for the load vector and the solution
    allocate( bb(      A_CSR%i1:A_CSR%i2))
    allocate( uv_bkuv( A_CSR%i1:A_CSR%i2))

    ! Fill in the current velocity solution
    do ti = self%mesh%ti1, self%mesh%ti2
    do k  = 1, self%mesh%nz

      ! u
      row_tikuv = self%mesh%tikuv2n( ti,k,1)
      uv_bkuv( row_tikuv) = self%u_bk( ti,k)

      ! v
      row_tikuv = self%mesh%tikuv2n( ti,k,2)
      uv_bkuv( row_tikuv) = self%v_bk( ti,k)

    end do
    end do

    ! == Construct the stiffness matrix for the linearised BPA
    ! ========================================================

    do row_tikuv = A_CSR%i1, A_CSR%i2

      ti = self%mesh%n2tikuv( row_tikuv,1)
      k  = self%mesh%n2tikuv( row_tikuv,2)
      uv = self%mesh%n2tikuv( row_tikuv,3)

      if (BC_prescr_mask_bk( ti,k) == 1) then
        ! Dirichlet boundary condition; velocities are prescribed for this triangle

        ! Stiffness matrix: diagonal element set to 1
        call A_CSR%add_entry( row_tikuv, row_tikuv, 1._dp)

        ! Load vector: prescribed velocity
        if     (uv == 1) then
          bb( row_tikuv) = BC_prescr_u_bk( ti,k)
        elseif (uv == 2) then
          bb( row_tikuv) = BC_prescr_v_bk( ti,k)
        else
          call crash('uv can only be 1 or 2!')
        end if

      elseif (self%mesh%TriBI( ti) == 1 .or. self%mesh%TriBI( ti) == 2) then
        ! Northern domain border

        call self%calc_BPA_stiffness_matrix_row_BC_north( A_CSR, bb, row_tikuv)

      elseif ( self%mesh%TriBI( ti) == 3 .or.  self%mesh%TriBI( ti) == 4) then
        ! Eastern domain border

        call self%calc_BPA_stiffness_matrix_row_BC_east( A_CSR, bb, row_tikuv)

      elseif (self%mesh%TriBI( ti) == 5 .or. self%mesh%TriBI( ti) == 6) then
        ! Southern domain border

        call self%calc_BPA_stiffness_matrix_row_BC_south( A_CSR, bb, row_tikuv)

      elseif (self%mesh%TriBI( ti) == 7 .or. self%mesh%TriBI( ti) == 8) then
        ! Western domain border

        call self%calc_BPA_stiffness_matrix_row_BC_west( A_CSR, bb, row_tikuv)

      elseif (k == 1) then
        ! Ice surface

        call self%calc_BPA_stiffness_matrix_row_BC_surf( ice, A_CSR, bb, row_tikuv)

      elseif (k == self%mesh%nz) then
        ! Ice base

        call self%calc_BPA_stiffness_matrix_row_BC_base( ice, A_CSR, bb, row_tikuv)

      else
        ! No boundary conditions apply; solve the BPA

        call self%calc_BPA_stiffness_matrix_row_free( A_CSR, bb, row_tikuv)

      end if

    end do

    call A_CSR%finalise()

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine assemble_BPA_linearised_matrix_eq

  subroutine calc_BPA_stiffness_matrix_row_free( self, A_CSR, bb, row_tikuv)
    !< Add coefficients to this matrix row to represent the linearised BPA

    ! The BPA reads;
    !
    !   d/dx [ 2 eta ( 2 du/dx + dv/dy )] + d/dy [ eta ( du/dy + dv/dx)] + d/dz [ eta du/dz] = -tau_dx
    !
    !   d/dy [ 2 eta ( 2 dv/dy + du/dx )] + d/dx [ eta ( dv/dx + du/dy)] + d/dz [ eta dv/dz] = -tau_dy
    !
    ! Using the chain rule, this expands to read:
    !
    !   4 eta d2u/dx2 + 4 deta/dx du/dx + 2 eta d2v/dxdy + 2 deta/dx dv/dy + ...
    !     eta d2u/dy2 +   deta/dy du/dy +   eta d2v/dxdy +   deta/dy dv/dx + ...
    !     eta d2u/dz2 +   deta/dz du/dz = -tau_dx
    !
    !   4 eta d2v/dy2 + 4 deta/dy dv/dy + 2 eta d2u/dxdy + 2 deta/dy du/dx + ...
    !     eta d2v/dx2 +   deta/dx dv/dx +   eta d2u/dxdy +   deta/dx du/dy + ...
    !     eta d2v/dz2 +   deta/dz dv/dz = -tau_dy
    !
    ! Rearranging to gather the terms involving u and v gives:
    !
    !   4 eta d2u/dx2  + 4 deta/dx du/dx + eta d2u/dy2 + deta/dy du/dy + eta d2u/dz2 + deta/dz du/dz + ...
    !   3 eta d2v/dxdy + 2 deta/dx dv/dy +               deta/dy dv/dx = -tau_dx
    !
    !   4 eta d2v/dy2  + 4 deta/dy dv/dy + eta d2v/dx2 + deta/dx dv/dx + eta d2v/dz2 + deta/dz dv/dz + ...
    !   3 eta d2u/dxdy + 2 deta/dy du/dx +               deta/dx du/dy = -tau_dy
    !
    ! We define the velocities u,v, and the driving stress tau_d on the bk-grid
    ! (triangles, regular vertical), and the effective viscosity eta on the ak-grid
    ! (vertices, regular vertical) and the bks-grid (triangles, staggered vertical).
    ! From this, we then calculate the gradients of eta on the bk-grid.

    ! In/output variables:
    class(type_momentum_balance_solver_BPA), intent(in   ) :: self
    type(type_CSR_matrix_dp),                intent(inout) :: A_CSR
    real(dp), dimension(A_CSR%i1:A_CSR%i2),  intent(inout) :: bb
    integer,                                 intent(in   ) :: row_tikuv

    ! Local variables:
    integer                             :: ti, k, uv
    real(dp)                            :: eta, deta_dx, deta_dy, deta_dz, tau_dx, tau_dy
    integer,  dimension(:), allocatable :: single_row_ind
    real(dp), dimension(:), allocatable :: single_row_ddx_val
    real(dp), dimension(:), allocatable :: single_row_ddy_val
    real(dp), dimension(:), allocatable :: single_row_ddz_val
    real(dp), dimension(:), allocatable :: single_row_d2dx2_val
    real(dp), dimension(:), allocatable :: single_row_d2dxdy_val
    real(dp), dimension(:), allocatable :: single_row_d2dy2_val
    real(dp), dimension(:), allocatable :: single_row_d2dz2_val
    integer                             :: single_row_nnz
    integer                             :: row_tik
    real(dp)                            :: Au, Av
    integer                             :: n, row_tjkk, tj, kk, col_tjkku, col_tjkkv

    ! Relevant indices for this triangle and layer
    ti = self%mesh%n2tikuv( row_tikuv,1)
    k  = self%mesh%n2tikuv( row_tikuv,2)
    uv = self%mesh%n2tikuv( row_tikuv,3)

    ! eta, deta/dx, deta/dy, deta/dz, tau_dx, and tau_dy on this triangle and layer
    eta     = self%eta_bk(     ti,k)
    deta_dx = self%deta_dx_bk( ti,k)
    deta_dy = self%deta_dy_bk( ti,k)
    deta_dz = self%deta_dz_bk( ti,k)
    tau_dx  = self%tau_dx_b(   ti  )
    tau_dy  = self%tau_dy_b(   ti  )

    ! allocate memory for single matrix rows
    allocate( single_row_ind(        self%mesh%nC_mem*3*2))
    allocate( single_row_ddx_val(    self%mesh%nC_mem*3*2))
    allocate( single_row_ddy_val(    self%mesh%nC_mem*3*2))
    allocate( single_row_ddz_val(    self%mesh%nC_mem*3*2))
    allocate( single_row_d2dx2_val(  self%mesh%nC_mem*3*2))
    allocate( single_row_d2dxdy_val( self%mesh%nC_mem*3*2))
    allocate( single_row_d2dy2_val(  self%mesh%nC_mem*3*2))
    allocate( single_row_d2dz2_val(  self%mesh%nC_mem*3*2))

    ! Read coefficients of the operator matrices
    row_tik = self%mesh%tik2n( ti,k)
    call self%mesh%M2_ddx_bk_bk%read_single_row(    row_tik, single_row_ind, single_row_ddx_val   , single_row_nnz)
    call self%mesh%M2_ddy_bk_bk%read_single_row(    row_tik, single_row_ind, single_row_ddy_val   , single_row_nnz)
    call self%mesh%M2_ddz_bk_bk%read_single_row(    row_tik, single_row_ind, single_row_ddz_val   , single_row_nnz)
    call self%mesh%M2_d2dx2_bk_bk%read_single_row(  row_tik, single_row_ind, single_row_d2dx2_val , single_row_nnz)
    call self%mesh%M2_d2dxdy_bk_bk%read_single_row( row_tik, single_row_ind, single_row_d2dxdy_val, single_row_nnz)
    call self%mesh%M2_d2dy2_bk_bk%read_single_row(  row_tik, single_row_ind, single_row_d2dy2_val , single_row_nnz)
    call self%mesh%M2_d2dz2_bk_bk%read_single_row(  row_tik, single_row_ind, single_row_d2dz2_val , single_row_nnz)

    if (uv == 1) then
      ! x-component

      do n = 1, single_row_nnz

        ! Relevant indices for this neighbouring triangle
        row_tjkk  = single_row_ind( n)
        tj        = self%mesh%n2tik( row_tjkk,1)
        kk        = self%mesh%n2tik( row_tjkk,2)
        col_tjkku = self%mesh%tikuv2n( tj,kk,1)
        col_tjkkv = self%mesh%tikuv2n( tj,kk,2)

        !   4 eta d2u/dx2  + 4 deta/dx du/dx + eta d2u/dy2 + deta/dy du/dy + eta d2u/dz2 + deta/dz du/dz + ...
        !   3 eta d2v/dxdy + 2 deta/dx dv/dy +               deta/dy dv/dx = -tau_dx

        ! Combine the mesh operators
        Au = 4._dp * eta     * single_row_d2dx2_val(  n) + &  ! 4  eta    d2u/dx2
             4._dp * deta_dx * single_row_ddx_val(    n) + &  ! 4 deta/dx  du/dx
                     eta     * single_row_d2dy2_val(  n) + &  !    eta    d2u/dy2
                     deta_dy * single_row_ddy_val(    n) + &  !   deta/dy  du/dy
                     eta     * single_row_d2dz2_val(  n) + &  !    eta    d2u/dz2
                     deta_dz * single_row_ddz_val(    n)      !   deta/dz  du/dz

        Av = 3._dp * eta     * single_row_d2dxdy_val( n) + &  ! 3  eta    d2v/dxdy
             2._dp * deta_dx * single_row_ddy_val(    n) + &  ! 2 deta/dx  dv/dy
                     deta_dy * single_row_ddx_val(    n)      !   deta/dy  dv/dx

        ! Add coefficients to the stiffness matrix
        call A_CSR%add_entry( row_tikuv, col_tjkku, Au)
        call A_CSR%add_entry( row_tikuv, col_tjkkv, Av)

      end do

      ! Load vector
      bb( row_tikuv) = -tau_dx

    elseif (uv == 2) then
      ! y-component

      do n = 1, single_row_nnz

        ! Releuant indices for this neighbouring triangle
        row_tjkk  = single_row_ind( n)
        tj        = self%mesh%n2tik( row_tjkk,1)
        kk        = self%mesh%n2tik( row_tjkk,2)
        col_tjkku = self%mesh%tikuv2n( tj,kk,1)
        col_tjkkv = self%mesh%tikuv2n( tj,kk,2)

        !   4 eta d2v/dy2  + 4 deta/dy dv/dy + eta d2v/dx2 + deta/dx dv/dx + eta d2v/dz2 + deta/dz dv/dz + ...
        !   3 eta d2u/dxdy + 2 deta/dy du/dx +               deta/dx du/dy = -tau_dy

        ! Combine the mesh operators
        Av = 4._dp * eta     * single_row_d2dy2_val(  n) + &  ! 4  eta    d2v/dy2
             4._dp * deta_dy * single_row_ddy_val(    n) + &  ! 4 deta/dy  dv/dy
                     eta     * single_row_d2dx2_val(  n) + &  !    eta    d2v/dx2
                     deta_dx * single_row_ddx_val(    n) + &  !   deta/dx  dv/dx
                     eta     * single_row_d2dz2_val(  n) + &  !    eta    d2v/dz2
                     deta_dz * single_row_ddz_val(    n)      !   deta/dz  dv/dz

        Au = 3._dp * eta     * single_row_d2dxdy_val( n) + &  ! 3  eta    d2u/dxdy
             2._dp * deta_dy * single_row_ddx_val(    n) + &  ! 2 deta/dy  du/dx
                     deta_dx * single_row_ddy_val(    n)      !   deta/dx  du/dy

        ! Add coefficients to the stiffness matrix
        call A_CSR%add_entry( row_tikuv, col_tjkkv, Av)
        call A_CSR%add_entry( row_tikuv, col_tjkku, Au)

      end do

      ! Load uector
      bb( row_tikuv) = -tau_dy

    else
      call crash('uv can only be 1 or 2!')
    end if

  end subroutine calc_BPA_stiffness_matrix_row_free

  subroutine calc_BPA_stiffness_matrix_row_BC_surf( self, ice, A_CSR, bb, row_tikuv)
    !< Add coefficients to this matrix row to represent the boundary conditions to the BPA at the ice surface

    ! At the ice surface (k=1), the zero-stress boundary condition implies that:
    !
    ! [1]     2 dh/dx (2 du/dx + dv/dy) + dh/dy (du/dy + dv/dx) - du/dz = 0
    !
    ! The two-sided differencing schemes for the first and second derivatives du/dz, d2u/dz2 read:
    !
    ! [2]     du/dz  =  dzeta/dz    (u( k+1) - u( k-1)) / (2 dzeta)
    ! [3]    d2u/dz2 = (dzeta/dz)^2 (u( k+1) + u( k-1) - 2 u( k)) / dzeta^2
    !
    ! A the ice surface, u( k-1) doesn't actually exist, but we can treat it as a "ghost point".
    ! Since, at the ice surface, we know du/dz from the zero-stress boundary condition [1]:
    !
    ! [4]     du/dz = 2 dh/dx (2 du/dx + dv/dy) + dh/dy (du/dy + dv/dx)
    !
    ! Substituting [4] into [2] yields:
    !
    !         dzeta/dz (u( k+1) - u( k-1)) / (2 dzeta) = 2 dh/dx (2 du/dx + dv/dy) + dh/dy (du/dy + dv/dx)
    !         u( k+1) - u( k-1) = 2 dzeta / (dzeta/dz)  (2 dh/dx (2 du/dx + dv/dy) + dh/dy (du/dy + dv/dx))
    ! [5]     u( k-1) = u( k+1) - 2 dzeta / (dzeta/dz)  (2 dh/dx (2 du/dx + dv/dy) + dh/dy (du/dy + dv/dx))
    !
    ! Substituting [5] into [3] yields:
    !
    !         d2u/dz2 = (dzeta/dz)^2 1/dzeta^2 ( -2 u( k) + u( k+1) + u( k-1))
    !                 = (dzeta/dz)^2 1/dzeta^2 ( -2 u( k) + 2 u( k+1) - 2 dzeta / (dzeta/dz) (2 dh/dx (2 du/dx + dv/dy) + dh/dy (du/dy + dv/dx)))
    ! [6]             = (dzeta/dz)^2 2/dzeta^2 ( u( k+1) - u( k) - dzeta / (dzeta/dz) (2 dh/dx (2 du/dx + dv/dy) + dh/dy (du/dy + dv/dx)))
    !
    ! The product-rule-expanded form of the BPA reads:
    !
    ! [7]     4 eta d2u/dx2  + 4 deta/dx du/dx + ...
    !           eta d2u/dy2  +   deta/dy du/dy + ...
    !         3 eta d2v/dxdy + 2 deta/dx dv/dy + ...
    !                            deta/dy dv/dx + ...
    !           eta d2u/dz2  +   deta/dz du/dz = - tau_d,x
    !
    ! Substituting [4] and [6] into [7] yields:
    !
    !         4 eta d2u/dx2  + 4 deta/dx du/dx + ...
    !           eta d2u/dy2  +   deta/dy du/dy + ...
    !         3 eta d2v/dxdy + 2 deta/dx dv/dy + ...
    !                            deta/dy dv/dx + ...
    !           eta (dzeta/dz)^2 2/dzeta^2 ( u( k+1) - u( k) - dzeta / (dzeta/dz) (2 dh/dx (2 du/dx + dv/dy) + dh/dy (du/dy + dv/dx))) + ...
    !          deta/dx 2 dh/dx (2 du/dx + dv/dy) + dh/dy (du/dy + dv/dx) = -tau_d,x
    !
    ! Rearranging to group the terms involving du/dx, du/dy, dv/dx, du/dy yields:
    !
    !         4 eta d2u/dx2 + eta d2u/dy2 + 3 eta d2v/dxdy + ...
    !         du/dx   [ 4 deta/dx + eta (dzeta/dz)^2 ( 2 / dzeta^2) (-dzeta / (dzeta/dz)) 4 dh/dx + 4 deta/dz dh/dx] + ...
    !         du/dy   [   deta/dy + eta (dzeta/dz)^2 ( 2 / dzeta^2) (-dzeta / (dzeta/dz))   dh/dy +   deta/dz dh/dy] + ...
    !         dv/dy   [ 2 deta/dx + eta (dzeta/dz)^2 ( 2 / dzeta^2) (-dzeta / (dzeta/dz)) 2 dh/dx + 2 deta/dz dh/dx] + ...
    !         dv/dx   [   deta/dy + eta (dzeta/dz)^2 ( 2 / dzeta^2) (-dzeta / (dzeta/dz)    dh/dy +   deta/dz dh/dy] + ...
    !         u( k+1) [             eta (dzeta/dz)^2 ( 2 / dzeta^2)] + ...
    !         u( k  ) [            -eta (dzeta/dz)^2 ( 2 / dzeta^2)] = -tau_d,x
    !
    ! Rearranging some of the cancelling dzeta/dz and dzeta terms yields:
    !
    ! [8]     4 eta d2u/dx2 + eta d2u/dy2 + 3 eta d2v/dxdy + ...
    !         du/dx   [ 4 deta/dx - 2 eta / dzeta    dzeta/dz   4 dh/dx + 4 deta/dz dh/dx] + ...
    !         du/dy   [   deta/dy - 2 eta / dzeta    dzeta/dz     dh/dy +   deta/dz dh/dy] + ...
    !         dv/dy   [ 2 deta/dx - 2 eta / dzeta    dzeta/dz   2 dh/dx + 2 deta/dz dh/dx] + ...
    !         dv/dx   [   deta/dy - 2 eta / dzeta    dzeta/dz     dh/dy +   deta/dz dh/dy] + ...
    !         u( k+1) [             2 eta / dzeta^2 (dzeta/dz)^2] + ...
    !         u( k  ) [            -2 eta / dzeta^2 (dzeta/dz)^2] = -tau_d,x

    ! In/output variables:
    class(type_momentum_balance_solver_BPA), intent(in   ) :: self
    class(atype_ice_model_data),             intent(in   ) :: ice
    type(type_CSR_matrix_dp),                intent(inout) :: A_CSR
    real(dp), dimension(A_CSR%i1:A_CSR%i2),  intent(inout) :: bb
    integer,                                 intent(in   ) :: row_tikuv

    ! Local variables:
    integer                             :: ti, k, uv
    real(dp)                            :: eta, deta_dx, deta_dy, deta_dz, tau_dx, tau_dy, dh_dx, dh_dy, dzeta_dz, dzeta
    integer,  dimension(:), allocatable :: single_row_ind
    real(dp), dimension(:), allocatable :: single_row_ddx_val
    real(dp), dimension(:), allocatable :: single_row_ddy_val
    real(dp), dimension(:), allocatable :: single_row_d2dx2_val
    real(dp), dimension(:), allocatable :: single_row_d2dxdy_val
    real(dp), dimension(:), allocatable :: single_row_d2dy2_val
    integer                             :: single_row_nnz
    integer                             :: row_tik
    real(dp)                            :: cu_dudx, cu_dudy, cu_d2udx2, cu_d2udy2, cu_dvdx, cu_dvdy, cu_d2vdxdy, cu_uk, cu_ukp1
    real(dp)                            :: cv_dvdy, cv_dvdx, cv_d2vdy2, cv_d2vdx2, cv_dudy, cv_dudx, cv_d2udxdy, cv_vk, cv_vkp1
    real(dp)                            :: Au, Av
    integer                             :: n, row_tjkk, tj, kk, col_tjkku, col_tjkkv

    ! Relevant indices for this triangle and layer
    ti = self%mesh%n2tikuv( row_tikuv,1)
    k  = self%mesh%n2tikuv( row_tikuv,2)
    uv = self%mesh%n2tikuv( row_tikuv,3)

    ! Safety
    if (k /= 1) call crash('Received k = {int_01}; only applicable at ice surface!', int_01 = k)

    ! eta, deta/dx, deta/dy, deta/dz, tau_dx, and tau_dy on this triangle and layer
    eta      = self%eta_bks(     ti,1)
    deta_dx  = self%deta_dx_bk(  ti,1)
    deta_dy  = self%deta_dy_bk(  ti,1)
    deta_dz  = self%deta_dz_bk(  ti,1)
    tau_dx   = self%tau_dx_b(    ti  )
    tau_dy   = self%tau_dy_b(    ti  )
    dh_dx    = self%dh_dx_b(     ti  )
    dh_dy    = self%dh_dy_b(     ti  )
    dzeta_dz = ice%dzeta_dz_bk( ti,1)
    dzeta    = self%mesh%zeta( 2) - self%mesh%zeta( 1)

    ! Calculate coefficients for the different gradients of u and v
    if (uv == 1) then
      ! x-component

      ! 4 eta d2u/dx2 + eta d2u/dy2 + 3 eta d2v/dxdy + ...
      ! du/dx   [ 4 deta/dx - 2 eta / dzeta    dzeta/dz   4 dh/dx + 4 deta/dz dh/dx] + ...
      ! du/dy   [   deta/dy - 2 eta / dzeta    dzeta/dz     dh/dy +   deta/dz dh/dy] + ...
      ! dv/dy   [ 2 deta/dx - 2 eta / dzeta    dzeta/dz   2 dh/dx + 2 deta/dz dh/dx] + ...
      ! dv/dx   [   deta/dy - 2 eta / dzeta    dzeta/dz     dh/dy +   deta/dz dh/dy] + ...
      ! u( k+1) [             2 eta / dzeta^2 (dzeta/dz)^2] + ...
      ! u( k  ) [            -2 eta / dzeta^2 (dzeta/dz)^2] = -tau_d,x

      cu_dudx    =  4._dp * deta_dx - 2._dp * eta / dzeta * dzeta_dz * 4._dp * dh_dx + 4._dp * deta_dz * dh_dx
      cu_dudy    =          deta_dy - 2._dp * eta / dzeta * dzeta_dz *         dh_dy +         deta_dz * dh_dy
      cu_dvdy    =  2._dp * deta_dx - 2._dp * eta / dzeta * dzeta_dz * 2._dp * dh_dx + 2._dp * deta_dz * dh_dx
      cu_dvdx    =          deta_dy - 2._dp * eta / dzeta * dzeta_dz *         dh_dy +         deta_dz * dh_dy
      cu_d2udx2  =  4._dp * eta
      cu_d2udy2  =          eta
      cu_d2vdxdy =  3._dp * eta
      cu_uk      = -2._dp * eta / dzeta**2 * dzeta_dz**2
      cu_ukp1    =  2._dp * eta / dzeta**2 * dzeta_dz**2

    elseif (uv == 2) then
      ! y-component

      ! 4 eta d2v/dy2 + eta d2v/dx2 + 3 eta d2u/dxdy + ...
      ! dv/dy   [ 4 deta/dy - 2 eta / dzeta    dzeta/dz   4 dh/dy + 4 deta/dz dh/dy] + ...
      ! dv/dx   [   deta/dx - 2 eta / dzeta    dzeta/dz     dh/dx +   deta/dz dh/dx] + ...
      ! du/dx   [ 2 deta/dy - 2 eta / dzeta    dzeta/dz   2 dh/dy + 2 deta/dz dh/dy] + ...
      ! du/dy   [   deta/dx - 2 eta / dzeta    dzeta/dz     dh/dx +   deta/dz dh/dx] + ...
      ! v( k+1) [             2 eta / dzeta^2 (dzeta/dz)^2] + ...
      ! v( k  ) [            -2 eta / dzeta^2 (dzeta/dz)^2] = -tau_d,y

      cv_dvdy    =  4._dp * deta_dy - 2._dp * eta / dzeta * dzeta_dz * 4._dp * dh_dy + 4._dp * deta_dz * dh_dy
      cv_dvdx    =          deta_dx - 2._dp * eta / dzeta * dzeta_dz *         dh_dx +         deta_dz * dh_dx
      cv_dudx    =  2._dp * deta_dy - 2._dp * eta / dzeta * dzeta_dz * 2._dp * dh_dy + 2._dp * deta_dz * dh_dy
      cv_dudy    =          deta_dx - 2._dp * eta / dzeta * dzeta_dz *         dh_dx +         deta_dz * dh_dx
      cv_d2vdy2  =  4._dp * eta
      cv_d2vdx2  =          eta
      cv_d2udxdy =  3._dp * eta
      cv_vk      = -2._dp * eta / dzeta**2 * dzeta_dz**2
      cv_vkp1    =  2._dp * eta / dzeta**2 * dzeta_dz**2

    else
      call crash('uv can only be 1 or 2!')
    end if

    ! allocate memory for single matrix rows
    allocate( single_row_ind(        self%mesh%nC_mem*3*2))
    allocate( single_row_ddx_val(    self%mesh%nC_mem*3*2))
    allocate( single_row_ddy_val(    self%mesh%nC_mem*3*2))
    allocate( single_row_d2dx2_val(  self%mesh%nC_mem*3*2))
    allocate( single_row_d2dxdy_val( self%mesh%nC_mem*3*2))
    allocate( single_row_d2dy2_val(  self%mesh%nC_mem*3*2))

    ! Read coefficients of the operator matrices
    row_tik = self%mesh%tik2n( ti,k)
    call self%mesh%M2_ddx_bk_bk%read_single_row(    row_tik, single_row_ind, single_row_ddx_val   , single_row_nnz)
    call self%mesh%M2_ddy_bk_bk%read_single_row(    row_tik, single_row_ind, single_row_ddy_val   , single_row_nnz)
    call self%mesh%M2_d2dx2_bk_bk%read_single_row(  row_tik, single_row_ind, single_row_d2dx2_val , single_row_nnz)
    call self%mesh%M2_d2dxdy_bk_bk%read_single_row( row_tik, single_row_ind, single_row_d2dxdy_val, single_row_nnz)
    call self%mesh%M2_d2dy2_bk_bk%read_single_row(  row_tik, single_row_ind, single_row_d2dy2_val , single_row_nnz)

    if (uv == 1) then
      ! x-component

      do n = 1, single_row_nnz

        ! Relevant indices for this neighbouring triangle
        row_tjkk  = single_row_ind( n)
        tj        = self%mesh%n2tik( row_tjkk,1)
        kk        = self%mesh%n2tik( row_tjkk,2)
        col_tjkku = self%mesh%tikuv2n( tj,kk,1)
        col_tjkkv = self%mesh%tikuv2n( tj,kk,2)

        ! Combine coefficients
        Au = cu_dudx    * single_row_ddx_val(    n) + &
             cu_dudy    * single_row_ddy_val(    n) + &
             cu_d2udx2  * single_row_d2dx2_val(  n) + &
             cu_d2udy2  * single_row_d2dy2_val(  n)
        if (tj == ti .and. kk == k  ) Au = Au + cu_uk
        if (tj == ti .and. kk == k+1) Au = Au + cu_ukp1

        Av = cu_dvdy    * single_row_ddy_val(    n) + &
             cu_dvdx    * single_row_ddx_val(    n) + &
             cu_d2vdxdy * single_row_d2dxdy_val( n)

        ! Add coefficients to the stiffness matrix
        call A_CSR%add_entry( row_tikuv, col_tjkku, Au)
        call A_CSR%add_entry( row_tikuv, col_tjkkv, Av)

      end do

      ! Load vector
      bb( row_tikuv) = -tau_dx

    elseif (uv == 2) then
      ! y-component

      do n = 1, single_row_nnz

        ! Relevant indices for this neighbouring triangle
        row_tjkk  = single_row_ind( n)
        tj        = self%mesh%n2tik( row_tjkk,1)
        kk        = self%mesh%n2tik( row_tjkk,2)
        col_tjkkv = self%mesh%tikuv2n( tj,kk,2)
        col_tjkku = self%mesh%tikuv2n( tj,kk,1)

        ! Combine coefficients
        Av = cv_dvdy    * single_row_ddy_val(    n) + &
             cv_dvdx    * single_row_ddx_val(    n) + &
             cv_d2vdy2  * single_row_d2dy2_val(  n) + &
             cv_d2vdx2  * single_row_d2dx2_val(  n)
        if (tj == ti .and. kk == k  ) Av = Av + cv_vk
        if (tj == ti .and. kk == k+1) Av = Av + cv_vkp1

        Au = cv_dudx    * single_row_ddx_val(    n) + &
             cv_dudy    * single_row_ddy_val(    n) + &
             cv_d2udxdy * single_row_d2dxdy_val( n)

        ! Add coefficients to the stiffness matrix
        call A_CSR%add_entry( row_tikuv, col_tjkkv, Av)
        call A_CSR%add_entry( row_tikuv, col_tjkku, Au)

      end do

      ! Load uector
      bb( row_tikuv) = -tau_dy

    else
      call crash('uv can only be 1 or 2!')
    end if

  end subroutine calc_BPA_stiffness_matrix_row_BC_surf

  subroutine calc_BPA_stiffness_matrix_row_BC_base( self, ice, A_CSR, bb, row_tikuv)
    !< Add coefficients to this matrix row to represent the boundary conditions to the BPA at the ice base

    ! At the ice surface (k=1), the zero-stress boundary condition implies that:
    !
    ! [1]     2 db/dx (2 du/dx + dv/dy) + db/dy (du/dy + dv/dx) - du/dz + beta_b/eta u = 0
    !
    ! As a shorthand, define:
    !
    ! [2]     P = du/dz = 2 db/dx (2 du/dx + dv/dy) + db/dy (du/dy + dv/dx) + beta_b/eta u
    !
    ! The two-sided differencing schemes for the first and second derivatives du/dz, d2u/dz2 read:
    !
    ! [3]     du/dz  =  dzeta/dz    (u( k+1) - u( k-1)          ) / (2 dzeta)
    ! [4]    d2u/dz2 = (dzeta/dz)^2 (u( k+1) + u( k-1) - 2 u( k)) /    dzeta^2
    !
    ! Substituting [3] into [2] and rearranging yields an expression for the ghost node u( k+1):
    !
    !         dzeta/dz (u( k+1) - u( k-1)) / (2 dzeta) = P
    !         u( k+1) - u( k-1) = 2 P dzeta / dzeta_dz
    ! [5]     u( k+1) = u( k-1) + 2 P dzeta / dzeta_dz
    !
    ! We then substitute [5] into [4] to find an expression for d2u/dz2 that no longer depends
    ! on the ghost node u( k+1):
    !
    !         d2u/dz2 = (dzeta/dz)^2 1 / dzeta^2 (u( k+1) + u( k-1) - 2 u( k))
    !                 = (dzeta/dz)^2 1 / dzeta^2 (u( k-1) + 2 P dzeta / dzeta_dz + u( k-1) - 2 u( k))
    !                 = (dzeta/dz)^2 1 / dzeta^2 (2 u( k-1) - 2 u( k) + 2 P dzeta / dzeta_dz)
    ! [6]             = (dzeta/dz)^2 2 / dzeta^2 (u( k-1) - u( k) + P dzeta / dzeta_dz)
    !
    ! The product-rule-expanded form of the BPA reads:
    !
    ! [7]     4 eta d2u/dx2  + 4 deta/dx du/dx + eta d2u/dy2  +   deta/dy du/dy + ...
    !         3 eta d2v/dxdy + 2 deta/dx dv/dy + deta/dy dv/dx + ...
    !           eta d2u/dz2  +   deta/dz du/dz = -tau_d,x
    !
    ! Substituting [6] and [2] into [7] yields:
    !
    ! [8]     4 eta d2u/dx2  + 4 deta/dx du/dx + eta d2u/dy2  +   deta/dy du/dy + ...
    !         3 eta d2v/dxdy + 2 deta/dx dv/dy + deta/dy dv/dx + ...
    !           eta [(dzeta/dz)^2 2 / dzeta^2 (u( k-1) - u( k) + P dzeta / dzeta_dz)] + ...
    !          deta/dz P = -tau_d,x
    !
    ! As more shorthand, define:
    !
    ! [9]     Q = eta (dzeta/dz)^2 2 / dzeta^2 = 2 eta / dzeta^2 (dzeta/dz)^2
    !
    ! Substituting [9] into [8] yields:
    !
    ! [10]    4 eta d2u/dx2  + 4 deta/dx du/dx + eta d2u/dy2  +   deta/dy du/dy + ...
    !         3 eta d2v/dxdy + 2 deta/dx dv/dy + deta/dy dv/dx + ...
    !         Q u( k-1) - Q u( k) + P Q dzeta / dzeta_dz)] + deta/dz P = -tau_d,x
    !
    ! As even more shorthand, define:
    !
    ! [11]    R = Q dzeta / (dzeta/dz) + deta/dz = 2 eta / dzeta dzeta/dz + deta/dz
    !
    ! Substituting [11] into [10] yields:
    !
    ! [12]    4 eta d2u/dx2  + 4 deta/dx du/dx + eta d2u/dy2  +   deta/dy du/dy + ...
    !         3 eta d2v/dxdy + 2 deta/dx dv/dy + deta/dy dv/dx + ...
    !         Q u( k-1) - Q u( k) + R P = -tau_d,x
    !
    ! Substituting [2] back into [12] yields:
    !
    ! [12]    4 eta d2u/dx2  + 4 deta/dx du/dx + eta d2u/dy2  +   deta/dy du/dy + ...
    !         3 eta d2v/dxdy + 2 deta/dx dv/dy + deta/dy dv/dx + ...
    !         Q u( k-1) - Q u( k) + ...
    !         R [2 db/dx (2 du/dx + dv/dy) + db/dy (du/dy + dv/dx) + beta_b/eta u] = -tau_d,x
    !
    ! Rearranging to collect the terms involving du/dx, du/dy, dv/dx, dv/dy yields:
    !
    ! [13]    4 eta d2u/dx2  + eta d2u/dy2 + 3 eta d2v/dxdy + ...
    !         du/dx   [ 4 deta/dx + 4 R db/dx] + ...
    !         du/dy   [   deta/dy +   R db/dy] + ...
    !         dv/dy   [ 2 deta/dx + 2 R db/dx] + ...
    !         dv/dx   [   deta/dy +   R db/dy] + ...
    !         u( k  ) [R beta_b/eta - Q] + ...
    !         u( k-1) [Q] = -tau_d,x

    ! In/output variables:
    class(type_momentum_balance_solver_BPA), intent(in   ) :: self
    class(atype_ice_model_data),             intent(in   ) :: ice
    type(type_CSR_matrix_dp),                intent(inout) :: A_CSR
    real(dp), dimension(A_CSR%i1:A_CSR%i2),  intent(inout) :: bb
    integer,                                 intent(in   ) :: row_tikuv

    ! Local variables:
    integer                             :: ti, k, uv
    real(dp)                            :: eta, deta_dx, deta_dy, deta_dz, tau_dx, tau_dy, db_dx, db_dy, dzeta_dz, dzeta, basal_friction_coefficient, Q, R
    integer,  dimension(:), allocatable :: single_row_ind
    real(dp), dimension(:), allocatable :: single_row_ddx_val
    real(dp), dimension(:), allocatable :: single_row_ddy_val
    real(dp), dimension(:), allocatable :: single_row_d2dx2_val
    real(dp), dimension(:), allocatable :: single_row_d2dxdy_val
    real(dp), dimension(:), allocatable :: single_row_d2dy2_val
    integer                             :: single_row_nnz
    integer                             :: row_tik
    real(dp)                            :: cu_dudx, cu_dudy, cu_d2udx2, cu_d2udy2, cu_dvdx, cu_dvdy, cu_d2vdxdy, cu_uk, cu_ukm1
    real(dp)                            :: cv_dvdy, cv_dvdx, cv_d2vdy2, cv_d2vdx2, cv_dudy, cv_dudx, cv_d2udxdy, cv_vk, cv_vkm1
    real(dp)                            :: Au, Av
    integer                             :: n, row_tjkk, tj, kk, col_tjkku, col_tjkkv

    ! Relevant indices for this triangle and layer
    ti = self%mesh%n2tikuv( row_tikuv,1)
    k  = self%mesh%n2tikuv( row_tikuv,2)
    uv = self%mesh%n2tikuv( row_tikuv,3)

    ! Safety
    if (k /= self%mesh%nz) call crash('Received k = {int_01}; only applicable at ice base!', int_01 = k)

    ! Exception for the case of no sliding
    if (C%choice_sliding_law == 'no_sliding') then
      ! u = v = 0
      call A_CSR%add_entry( row_tikuv, row_tikuv, 1._dp)
      bb( row_tikuv) = 0._dp
      return
    end if

    ! eta, deta/dx, deta/dy, deta/dz, tau_dx, and tau_dy on this triangle and layer
    eta                        = self%eta_bks(                      ti,self%mesh%nz-1)
    deta_dx                    = self%deta_dx_bk(                   ti,self%mesh%nz  )
    deta_dy                    = self%deta_dy_bk(                   ti,self%mesh%nz  )
    deta_dz                    = self%deta_dz_bk(                   ti,self%mesh%nz  )
    tau_dx                     = self%tau_dx_b(                     ti               )
    tau_dy                     = self%tau_dy_b(                     ti               )
    db_dx                      = self%db_dx_b(                      ti               )
    db_dy                      = self%db_dy_b(                      ti               )
    basal_friction_coefficient = self%basal_friction_coefficient_b( ti               )
    dzeta_dz                   = ice%dzeta_dz_bk(                   ti,self%mesh%nz  )
    dzeta                      = self%mesh%zeta( self%mesh%nz) - self%mesh%zeta( self%mesh%nz-1)

    ! [9]     Q = eta (dzeta/dz)^2 2 / dzeta^2 = 2 eta / dzeta^2 (dzeta/dz)^2
    Q = 2._dp * eta / dzeta**2 * dzeta_dz**2

    ! [11]    R = Q dzeta / (dzeta/dz) + deta/dz = 2 eta / dzeta dzeta/dz + deta/dz
    R = 2._dp * eta / dzeta    * dzeta_dz    + deta_dz

    ! Calculate coefficients for the different gradients of u and v
    if (uv == 1) then
      ! x-component

      ! [13]    4 eta d2u/dx2  + eta d2u/dy2 + 3 eta d2v/dxdy + ...
      !         du/dx   [ 4 deta/dx + 4 R db/dx] + ...
      !         du/dy   [   deta/dy +   R db/dy] + ...
      !         dv/dy   [ 2 deta/dx + 2 R db/dx] + ...
      !         dv/dx   [   deta/dy +   R db/dy] + ...
      !         u( k  ) [R beta_b/eta - Q] + ...
      !         u( k-1) [Q] = -tau_d,x

      cu_dudx    =  4._dp * deta_dx + 4._dp * R * db_dx
      cu_dudy    =          deta_dy +         R * db_dy
      cu_dvdy    =  2._dp * deta_dx + 2._dp * R * db_dx
      cu_dvdx    =          deta_dy +         R * db_dy
      cu_d2udx2  =  4._dp * eta
      cu_d2udy2  =          eta
      cu_d2vdxdy =  3._dp * eta
      cu_uk      = ( R * basal_friction_coefficient / eta) - Q
      cu_ukm1    = Q

    elseif (uv == 2) then
      ! y-component

      ! [13]    4 eta d2v/dy2  + eta d2v/dx2 + 3 eta d2u/dydx + ...
      !         dv/dy   [ 4 deta/dy + 4 R db/dy] + ...
      !         dv/dx   [   deta/dx +   R db/dx] + ...
      !         du/dx   [ 2 deta/dy + 2 R db/dy] + ...
      !         du/dy   [   deta/dx +   R db/dx] + ...
      !         v( k  ) [R beta_b/eta - Q] + ...
      !         v( k-1) [Q] = -tau_d,y

      cv_dvdy    =  4._dp * deta_dy + 4._dp * R * db_dy
      cv_dvdx    =          deta_dx +         R * db_dx
      cv_dudx    =  2._dp * deta_dy + 2._dp * R * db_dy
      cv_dudy    =          deta_dx +         R * db_dx
      cv_d2vdy2  =  4._dp * eta
      cv_d2vdx2  =          eta
      cv_d2udxdy =  3._dp * eta
      cv_vk      = ( R * basal_friction_coefficient / eta) - Q
      cv_vkm1    = Q

    else
      call crash('uv can only be 1 or 2!')
    end if

    ! allocate memory for single matrix rows
    allocate( single_row_ind(        self%mesh%nC_mem*3*2))
    allocate( single_row_ddx_val(    self%mesh%nC_mem*3*2))
    allocate( single_row_ddy_val(    self%mesh%nC_mem*3*2))
    allocate( single_row_d2dx2_val(  self%mesh%nC_mem*3*2))
    allocate( single_row_d2dxdy_val( self%mesh%nC_mem*3*2))
    allocate( single_row_d2dy2_val(  self%mesh%nC_mem*3*2))

    ! Read coefficients of the operator matrices
    row_tik = self%mesh%tik2n( ti,k)
    call self%mesh%M2_ddx_bk_bk%read_single_row(    row_tik, single_row_ind, single_row_ddx_val   , single_row_nnz)
    call self%mesh%M2_ddy_bk_bk%read_single_row(    row_tik, single_row_ind, single_row_ddy_val   , single_row_nnz)
    call self%mesh%M2_d2dx2_bk_bk%read_single_row(  row_tik, single_row_ind, single_row_d2dx2_val , single_row_nnz)
    call self%mesh%M2_d2dxdy_bk_bk%read_single_row( row_tik, single_row_ind, single_row_d2dxdy_val, single_row_nnz)
    call self%mesh%M2_d2dy2_bk_bk%read_single_row(  row_tik, single_row_ind, single_row_d2dy2_val , single_row_nnz)

    if (uv == 1) then
      ! x-component

      do n = 1, single_row_nnz

        ! Relevant indices for this neighbouring triangle
        row_tjkk  = single_row_ind( n)
        tj        = self%mesh%n2tik( row_tjkk,1)
        kk        = self%mesh%n2tik( row_tjkk,2)
        col_tjkku = self%mesh%tikuv2n( tj,kk,1)
        col_tjkkv = self%mesh%tikuv2n( tj,kk,2)

        ! Combine coefficients
        Au = cu_dudx    * single_row_ddx_val(    n) + &
             cu_dudy    * single_row_ddy_val(    n) + &
             cu_d2udx2  * single_row_d2dx2_val(  n) + &
             cu_d2udy2  * single_row_d2dy2_val(  n)
        if (tj == ti .and. kk == k  ) Au = Au + cu_uk
        if (tj == ti .and. kk == k-1) Au = Au + cu_ukm1

        Av = cu_dvdy    * single_row_ddy_val(    n) + &
             cu_dvdx    * single_row_ddx_val(    n) + &
             cu_d2vdxdy * single_row_d2dxdy_val( n)

        ! Add coefficients to the stiffness matrix
        call A_CSR%add_entry( row_tikuv, col_tjkku, Au)
        call A_CSR%add_entry( row_tikuv, col_tjkkv, Av)

      end do

      ! Load vector
      bb( row_tikuv) = -tau_dx

    elseif (uv == 2) then
      ! y-component

      do n = 1, single_row_nnz

        ! Relevant indices for this neighbouring triangle
        row_tjkk  = single_row_ind( n)
        tj        = self%mesh%n2tik( row_tjkk,1)
        kk        = self%mesh%n2tik( row_tjkk,2)
        col_tjkkv = self%mesh%tikuv2n( tj,kk,2)
        col_tjkku = self%mesh%tikuv2n( tj,kk,1)

        ! Combine coefficients
        Av = cv_dvdy    * single_row_ddy_val(    n) + &
             cv_dvdx    * single_row_ddx_val(    n) + &
             cv_d2vdy2  * single_row_d2dy2_val(  n) + &
             cv_d2vdx2  * single_row_d2dx2_val(  n)
        if (tj == ti .and. kk == k  ) Av = Av + cv_vk
        if (tj == ti .and. kk == k-1) Av = Av + cv_vkm1

        Au = cv_dudx    * single_row_ddx_val(    n) + &
             cv_dudy    * single_row_ddy_val(    n) + &
             cv_d2udxdy * single_row_d2dxdy_val( n)

        ! Add coefficients to the stiffness matrix
        call A_CSR%add_entry( row_tikuv, col_tjkkv, Av)
        call A_CSR%add_entry( row_tikuv, col_tjkku, Au)

      end do

      ! Load uector
      bb( row_tikuv) = -tau_dy

    else
      call crash('uv can only be 1 or 2!')
    end if

  end subroutine calc_BPA_stiffness_matrix_row_BC_base

  subroutine calc_BPA_stiffness_matrix_row_BC_west( self, A_CSR, bb, row_tikuv)
    !< Add coefficients to this matrix row to represent boundary conditions at the
    !< western domain border.

    ! In/output variables:
    class(type_momentum_balance_solver_BPA), intent(in   ) :: self
    type(type_CSR_matrix_dp),                intent(inout) :: A_CSR
    real(dp), dimension(A_CSR%i1:A_CSR%i2),  intent(inout) :: bb
    integer,                                 intent(in   ) :: row_tikuv

    ! Local variables:
    integer                               :: ti,k,uv,row_ti
    integer                               :: tj, col_tjkuv
    integer,  dimension(self%mesh%nC_mem) :: ti_copy
    real(dp), dimension(self%mesh%nC_mem) :: wti_copy
    real(dp)                              :: u_fixed, v_fixed
    integer                               :: n, n_neighbours

    ti = self%mesh%n2tikuv( row_tikuv,1)
    k  = self%mesh%n2tikuv( row_tikuv,2)
    uv = self%mesh%n2tikuv( row_tikuv,3)
    row_ti = self%mesh%ti2n( ti)

    if (uv == 1) then
      ! x-component

      if     (C%BC_u_west == 'infinite') then
        ! du/dx = 0
        !
        ! notE: using the d/dx operator matrix doesn't always work well, not sure why...

        ! Set u on this triangle equal to the average value on its neighbours
        n_neighbours = 0
        do n = 1, 3
          tj = self%mesh%TriC( ti,n)
          if (tj == 0) CYCLE
          n_neighbours = n_neighbours + 1
          col_tjkuv = self%mesh%tikuv2n( tj,k,uv)
          call A_CSR%add_entry( row_tikuv, col_tjkuv, 1._dp)
        end do
        if (n_neighbours == 0) call crash('whaa!')
        call A_CSR%add_entry( row_tikuv, row_tikuv, -1._dp * real( n_neighbours,dp))

        ! Load vector
        bb( row_tikuv) = 0._dp

      elseif (C%BC_u_west == 'zero') then
        ! u = 0

        ! Stiffness matrix
        call A_CSR%add_entry( row_tikuv, row_tikuv, 1._dp)

        ! Load vector
        bb( row_tikuv) = 0._dp

      elseif (C%BC_u_west == 'periodic_ISMIP-HOM') then
        ! u(x,y) = u(x+-L/2,y+-L/2)

        ! Find the triangle ti_copy that is displaced by [x+-L/2,y+-L/2] relative to ti
        call find_ti_copy_ISMIP_HOM_periodic( self%mesh, C%refgeo_idealised_ISMIP_HOM_L, ti, ti_copy, wti_copy)

        ! Set value at ti equal to value at ti_copy
        call A_CSR%add_entry( row_tikuv, row_tikuv,  1._dp)
        u_fixed = 0._dp
        do n = 1, self%mesh%nC_mem
          tj = ti_copy( n)
          if (tj == 0) CYCLE
          u_fixed = u_fixed + wti_copy( n) * self%u_bk_prev( tj,k)
        end do
        ! Relax solution to improve stability
        u_fixed = (C%visc_it_relax * u_fixed) + ((1._dp - C%visc_it_relax) * self%u_bk_prev( ti,k))
        ! Set load vector
        bb( row_tikuv) = u_fixed

      else
        call crash('unknown BC_u_west "' // trim( C%BC_u_west) // '"!')
      end if

    elseif (uv == 2) then
      ! y-component

      if     (C%BC_v_west == 'infinite') then
        ! dv/dx = 0
        !
        ! notE: using the d/dx operator matrix doesn't always work well, not sure why...

        ! Set v on this triangle equal to the average value on its neighbours
        n_neighbours = 0
        do n = 1, 3
          tj = self%mesh%TriC( ti,n)
          if (tj == 0) CYCLE
          n_neighbours = n_neighbours + 1
          col_tjkuv = self%mesh%tikuv2n( tj,k,uv)
          call A_CSR%add_entry( row_tikuv, col_tjkuv, 1._dp)
        end do
        if (n_neighbours == 0) call crash('whaa!')
        call A_CSR%add_entry( row_tikuv, row_tikuv, -1._dp * real( n_neighbours,dp))

        ! Load vector
        bb( row_tikuv) = 0._dp

      elseif (C%BC_v_west == 'zero') then
        ! v = 0

        ! Stiffness matrix
        call A_CSR%add_entry( row_tikuv, row_tikuv, 1._dp)

        ! Load vector
        bb( row_tikuv) = 0._dp

      elseif (C%BC_v_west == 'periodic_ISMIP-HOM') then
        ! v(x,y) = v(x+-L/2,y+-L/2)

        ! Find the triangle ti_copy that is displaced by [x+-L/2,y+-L/2] relative to ti
        call find_ti_copy_ISMIP_HOM_periodic( self%mesh, C%refgeo_idealised_ISMIP_HOM_L, ti, ti_copy, wti_copy)

        ! Set value at ti equal to value at ti_copy
        call A_CSR%add_entry( row_tikuv, row_tikuv,  1._dp)
        v_fixed = 0._dp
        do n = 1, self%mesh%nC_mem
          tj = ti_copy( n)
          if (tj == 0) CYCLE
          v_fixed = v_fixed + wti_copy( n) * self%v_bk_prev( tj,k)
        end do
        ! Relax solution to improve stability
        v_fixed = (C%visc_it_relax * v_fixed) + ((1._dp - C%visc_it_relax) * self%v_bk_prev( ti,k))
        ! Set load vector
        bb( row_tikuv) = v_fixed

      else
        call crash('unknown BC_u_west "' // trim( C%BC_u_west) // '"!')
      end if

    else
      call crash('uv can only be 1 or 2!')
    end if

  end subroutine calc_BPA_stiffness_matrix_row_BC_west

  subroutine calc_BPA_stiffness_matrix_row_BC_east( self, A_CSR, bb, row_tikuv)
    !< Add coefficients to this matrix row to represent boundary conditions at the
    !< eastern domain border.

    ! In/output variables:
    class(type_momentum_balance_solver_BPA), intent(in   ) :: self
    type(type_CSR_matrix_dp),                intent(inout) :: A_CSR
    real(dp), dimension(A_CSR%i1:A_CSR%i2),  intent(inout) :: bb
    integer,                                 intent(in   ) :: row_tikuv

    ! Local variables:
    integer                               :: ti,k,uv,row_ti
    integer                               :: tj, col_tjkuv
    integer,  dimension(self%mesh%nC_mem) :: ti_copy
    real(dp), dimension(self%mesh%nC_mem) :: wti_copy
    real(dp)                              :: u_fixed, v_fixed
    integer                               :: n, n_neighbours

    ti = self%mesh%n2tikuv( row_tikuv,1)
    k  = self%mesh%n2tikuv( row_tikuv,2)
    uv = self%mesh%n2tikuv( row_tikuv,3)
    row_ti = self%mesh%ti2n( ti)

    if (uv == 1) then
      ! x-component

      if     (C%BC_u_east == 'infinite') then
        ! du/dx = 0
        !
        ! notE: using the d/dx operator matrix doesn't always work well, not sure why...

        ! Set u on this triangle equal to the average value on its neighbours
        n_neighbours = 0
        do n = 1, 3
          tj = self%mesh%TriC( ti,n)
          if (tj == 0) CYCLE
          n_neighbours = n_neighbours + 1
          col_tjkuv = self%mesh%tikuv2n( tj,k,uv)
          call A_CSR%add_entry( row_tikuv, col_tjkuv, 1._dp)
        end do
        if (n_neighbours == 0) call crash('whaa!')
        call A_CSR%add_entry( row_tikuv, row_tikuv, -1._dp * real( n_neighbours,dp))

        ! Load vector
        bb( row_tikuv) = 0._dp

      elseif (C%BC_u_east == 'zero') then
        ! u = 0

        ! Stiffness matrix
        call A_CSR%add_entry( row_tikuv, row_tikuv, 1._dp)

        ! Load vector
        bb( row_tikuv) = 0._dp

      elseif (C%BC_u_east == 'periodic_ISMIP-HOM') then
        ! u(x,y) = u(x+-L/2,y+-L/2)

        ! Find the triangle ti_copy that is displaced by [x+-L/2,y+-L/2] relative to ti
        call find_ti_copy_ISMIP_HOM_periodic( self%mesh, C%refgeo_idealised_ISMIP_HOM_L, ti, ti_copy, wti_copy)

        ! Set value at ti equal to value at ti_copy
        call A_CSR%add_entry( row_tikuv, row_tikuv,  1._dp)
        u_fixed = 0._dp
        do n = 1, self%mesh%nC_mem
          tj = ti_copy( n)
          if (tj == 0) CYCLE
          u_fixed = u_fixed + wti_copy( n) * self%u_bk_prev( tj,k)
        end do
        ! Relax solution to improve stability
        u_fixed = (C%visc_it_relax * u_fixed) + ((1._dp - C%visc_it_relax) * self%u_bk_prev( ti,k))
        ! Set load vector
        bb( row_tikuv) = u_fixed

      else
        call crash('unknown BC_u_east "' // trim( C%BC_u_east) // '"!')
      end if

    elseif (uv == 2) then
      ! y-component

      if     (C%BC_v_east == 'infinite') then
        ! dv/dx = 0
        !
        ! notE: using the d/dx operator matrix doesn't always work well, not sure why...

        ! Set v on this triangle equal to the average value on its neighbours
        n_neighbours = 0
        do n = 1, 3
          tj = self%mesh%TriC( ti,n)
          if (tj == 0) CYCLE
          n_neighbours = n_neighbours + 1
          col_tjkuv = self%mesh%tikuv2n( tj,k,uv)
          call A_CSR%add_entry( row_tikuv, col_tjkuv, 1._dp)
        end do
        if (n_neighbours == 0) call crash('whaa!')
        call A_CSR%add_entry( row_tikuv, row_tikuv, -1._dp * real( n_neighbours,dp))

        ! Load vector
        bb( row_tikuv) = 0._dp

      elseif (C%BC_v_east == 'zero') then
        ! v = 0

        ! Stiffness matrix
        call A_CSR%add_entry( row_tikuv, row_tikuv, 1._dp)

        ! Load vector
        bb( row_tikuv) = 0._dp

      elseif (C%BC_v_east == 'periodic_ISMIP-HOM') then
        ! v(x,y) = v(x+-L/2,y+-L/2)

        ! Find the triangle ti_copy that is displaced by [x+-L/2,y+-L/2] relative to ti
        call find_ti_copy_ISMIP_HOM_periodic( self%mesh, C%refgeo_idealised_ISMIP_HOM_L, ti, ti_copy, wti_copy)

        ! Set value at ti equal to value at ti_copy
        call A_CSR%add_entry( row_tikuv, row_tikuv,  1._dp)
        v_fixed = 0._dp
        do n = 1, self%mesh%nC_mem
          tj = ti_copy( n)
          if (tj == 0) CYCLE
          v_fixed = v_fixed + wti_copy( n) * self%v_bk_prev( tj,k)
        end do
        ! Relax solution to improve stability
        v_fixed = (C%visc_it_relax * v_fixed) + ((1._dp - C%visc_it_relax) * self%v_bk_prev( ti,k))
        ! Set load vector
        bb( row_tikuv) = v_fixed

      else
        call crash('unknown BC_u_east "' // trim( C%BC_u_east) // '"!')
      end if

    else
      call crash('uv can only be 1 or 2!')
    end if

  end subroutine calc_BPA_stiffness_matrix_row_BC_east

  subroutine calc_BPA_stiffness_matrix_row_BC_south( self, A_CSR, bb, row_tikuv)
    !< Add coefficients to this matrix row to represent boundary conditions at the
    !< southern domain border.

    ! In/output variables:
    class(type_momentum_balance_solver_BPA), intent(in   ) :: self
    type(type_CSR_matrix_dp),                intent(inout) :: A_CSR
    real(dp), dimension(A_CSR%i1:A_CSR%i2),  intent(inout) :: bb
    integer,                                 intent(in   ) :: row_tikuv

    ! Local variables:
    integer                               :: ti,k,uv,row_ti
    integer                               :: tj, col_tjkuv
    integer,  dimension(self%mesh%nC_mem) :: ti_copy
    real(dp), dimension(self%mesh%nC_mem) :: wti_copy
    real(dp)                              :: u_fixed, v_fixed
    integer                               :: n, n_neighbours

    ti = self%mesh%n2tikuv( row_tikuv,1)
    k  = self%mesh%n2tikuv( row_tikuv,2)
    uv = self%mesh%n2tikuv( row_tikuv,3)
    row_ti = self%mesh%ti2n( ti)

    if (uv == 1) then
      ! x-component

      if     (C%BC_u_south == 'infinite') then
        ! du/dy = 0
        !
        ! notE: using the d/dy operator matrix doesn't always work well, not sure why...

        ! Set u on this triangle equal to the average value on its neighbours
        n_neighbours = 0
        do n = 1, 3
          tj = self%mesh%TriC( ti,n)
          if (tj == 0) CYCLE
          n_neighbours = n_neighbours + 1
          col_tjkuv = self%mesh%tikuv2n( tj,k,uv)
          call A_CSR%add_entry( row_tikuv, col_tjkuv, 1._dp)
        end do
        if (n_neighbours == 0) call crash('whaa!')
        call A_CSR%add_entry( row_tikuv, row_tikuv, -1._dp * real( n_neighbours,dp))

        ! Load vector
        bb( row_tikuv) = 0._dp

      elseif (C%BC_u_south == 'zero') then
        ! u = 0

        ! Stiffness matrix
        call A_CSR%add_entry( row_tikuv, row_tikuv, 1._dp)

        ! Load vector
        bb( row_tikuv) = 0._dp

      elseif (C%BC_u_south == 'periodic_ISMIP-HOM') then
        ! u(x,y) = u(x+-L/2,y+-L/2)

        ! Find the triangle ti_copy that is displaced by [x+-L/2,y+-L/2] relative to ti
        call find_ti_copy_ISMIP_HOM_periodic( self%mesh, C%refgeo_idealised_ISMIP_HOM_L, ti, ti_copy, wti_copy)

        ! Set value at ti equal to value at ti_copy
        call A_CSR%add_entry( row_tikuv, row_tikuv,  1._dp)
        u_fixed = 0._dp
        do n = 1, self%mesh%nC_mem
          tj = ti_copy( n)
          if (tj == 0) CYCLE
          u_fixed = u_fixed + wti_copy( n) * self%u_bk_prev( tj,k)
        end do
        ! Relax solution to improve stability
        u_fixed = (C%visc_it_relax * u_fixed) + ((1._dp - C%visc_it_relax) * self%u_bk_prev( ti,k))
        ! Set load vector
        bb( row_tikuv) = u_fixed

      else
        call crash('unknown BC_u_south "' // trim( C%BC_u_south) // '"!')
      end if

    elseif (uv == 2) then
      ! y-component

      if     (C%BC_v_south == 'infinite') then
        ! dv/dy = 0
        !
        ! notE: using the d/dy operator matrix doesn't always work well, not sure why...

        ! Set v on this triangle equal to the average value on its neighbours
        n_neighbours = 0
        do n = 1, 3
          tj = self%mesh%TriC( ti,n)
          if (tj == 0) CYCLE
          n_neighbours = n_neighbours + 1
          col_tjkuv = self%mesh%tikuv2n( tj,k,uv)
          call A_CSR%add_entry( row_tikuv, col_tjkuv, 1._dp)
        end do
        if (n_neighbours == 0) call crash('whaa!')
        call A_CSR%add_entry( row_tikuv, row_tikuv, -1._dp * real( n_neighbours,dp))

        ! Load vector
        bb( row_tikuv) = 0._dp

      elseif (C%BC_v_south == 'zero') then
        ! v = 0

        ! Stiffness matrix
        call A_CSR%add_entry( row_tikuv, row_tikuv, 1._dp)

        ! Load vector
        bb( row_tikuv) = 0._dp

      elseif (C%BC_v_south == 'periodic_ISMIP-HOM') then
        ! v(x,y) = v(x+-L/2,y+-L/2)

        ! Find the triangle ti_copy that is displaced by [x+-L/2,y+-L/2] relative to ti
        call find_ti_copy_ISMIP_HOM_periodic( self%mesh, C%refgeo_idealised_ISMIP_HOM_L, ti, ti_copy, wti_copy)

        ! Set value at ti equal to value at ti_copy
        call A_CSR%add_entry( row_tikuv, row_tikuv,  1._dp)
        v_fixed = 0._dp
        do n = 1, self%mesh%nC_mem
          tj = ti_copy( n)
          if (tj == 0) CYCLE
          v_fixed = v_fixed + wti_copy( n) * self%v_bk_prev( tj,k)
        end do
        ! Relax solution to improve stability
        v_fixed = (C%visc_it_relax * v_fixed) + ((1._dp - C%visc_it_relax) * self%v_bk_prev( ti,k))
        ! Set load vector
        bb( row_tikuv) = v_fixed

      else
        call crash('unknown BC_u_south "' // trim( C%BC_u_south) // '"!')
      end if

    else
      call crash('uv can only be 1 or 2!')
    end if

  end subroutine calc_BPA_stiffness_matrix_row_BC_south

  subroutine calc_BPA_stiffness_matrix_row_BC_north( self, A_CSR, bb, row_tikuv)
    !< Add coefficients to this matrix row to represent boundary conditions at the
    !< northern domain border.

    ! In/output variables:
    class(type_momentum_balance_solver_BPA), intent(in   ) :: self
    type(type_CSR_matrix_dp),                intent(inout) :: A_CSR
    real(dp), dimension(A_CSR%i1:A_CSR%i2),  intent(inout) :: bb
    integer,                                 intent(in   ) :: row_tikuv

    ! Local variables:
    integer                               :: ti,k,uv,row_ti
    integer                               :: tj, col_tjkuv
    integer,  dimension(self%mesh%nC_mem) :: ti_copy
    real(dp), dimension(self%mesh%nC_mem) :: wti_copy
    real(dp)                              :: u_fixed, v_fixed
    integer                               :: n, n_neighbours

    ti = self%mesh%n2tikuv( row_tikuv,1)
    k  = self%mesh%n2tikuv( row_tikuv,2)
    uv = self%mesh%n2tikuv( row_tikuv,3)
    row_ti = self%mesh%ti2n( ti)

    if (uv == 1) then
      ! x-component

      if     (C%BC_u_north == 'infinite') then
        ! du/dy = 0
        !
        ! notE: using the d/dy operator matrix doesn't always work well, not sure why...

        ! Set u on this triangle equal to the average value on its neighbours
        n_neighbours = 0
        do n = 1, 3
          tj = self%mesh%TriC( ti,n)
          if (tj == 0) CYCLE
          n_neighbours = n_neighbours + 1
          col_tjkuv = self%mesh%tikuv2n( tj,k,uv)
          call A_CSR%add_entry( row_tikuv, col_tjkuv, 1._dp)
        end do
        if (n_neighbours == 0) call crash('whaa!')
        call A_CSR%add_entry( row_tikuv, row_tikuv, -1._dp * real( n_neighbours,dp))

        ! Load vector
        bb( row_tikuv) = 0._dp

      elseif (C%BC_u_north == 'zero') then
        ! u = 0

        ! Stiffness matrix
        call A_CSR%add_entry( row_tikuv, row_tikuv, 1._dp)

        ! Load vector
        bb( row_tikuv) = 0._dp

      elseif (C%BC_u_north == 'periodic_ISMIP-HOM') then
        ! u(x,y) = u(x+-L/2,y+-L/2)

        ! Find the triangle ti_copy that is displaced by [x+-L/2,y+-L/2] relative to ti
        call find_ti_copy_ISMIP_HOM_periodic( self%mesh, C%refgeo_idealised_ISMIP_HOM_L, ti, ti_copy, wti_copy)

        ! Set value at ti equal to value at ti_copy
        call A_CSR%add_entry( row_tikuv, row_tikuv,  1._dp)
        u_fixed = 0._dp
        do n = 1, self%mesh%nC_mem
          tj = ti_copy( n)
          if (tj == 0) CYCLE
          u_fixed = u_fixed + wti_copy( n) * self%u_bk_prev( tj,k)
        end do
        ! Relax solution to improve stability
        u_fixed = (C%visc_it_relax * u_fixed) + ((1._dp - C%visc_it_relax) * self%u_bk_prev( ti,k))
        ! Set load vector
        bb( row_tikuv) = u_fixed

      else
        call crash('unknown BC_u_north "' // trim( C%BC_u_north) // '"!')
      end if

    elseif (uv == 2) then
      ! y-component

      if     (C%BC_v_north == 'infinite') then
        ! dv/dy = 0
        !
        ! notE: using the d/dy operator matrix doesn't always work well, not sure why...

        ! Set v on this triangle equal to the average value on its neighbours
        n_neighbours = 0
        do n = 1, 3
          tj = self%mesh%TriC( ti,n)
          if (tj == 0) CYCLE
          n_neighbours = n_neighbours + 1
          col_tjkuv = self%mesh%tikuv2n( tj,k,uv)
          call A_CSR%add_entry( row_tikuv, col_tjkuv, 1._dp)
        end do
        if (n_neighbours == 0) call crash('whaa!')
        call A_CSR%add_entry( row_tikuv, row_tikuv, -1._dp * real( n_neighbours,dp))

        ! Load vector
        bb( row_tikuv) = 0._dp

      elseif (C%BC_v_north == 'zero') then
        ! v = 0

        ! Stiffness matrix
        call A_CSR%add_entry( row_tikuv, row_tikuv, 1._dp)

        ! Load vector
        bb( row_tikuv) = 0._dp

      elseif (C%BC_v_north == 'periodic_ISMIP-HOM') then
        ! v(x,y) = v(x+-L/2,y+-L/2)

        ! Find the triangle ti_copy that is displaced by [x+-L/2,y+-L/2] relative to ti
        call find_ti_copy_ISMIP_HOM_periodic( self%mesh, C%refgeo_idealised_ISMIP_HOM_L, ti, ti_copy, wti_copy)

        ! Set value at ti equal to value at ti_copy
        call A_CSR%add_entry( row_tikuv, row_tikuv,  1._dp)
        v_fixed = 0._dp
        do n = 1, self%mesh%nC_mem
          tj = ti_copy( n)
          if (tj == 0) CYCLE
          v_fixed = v_fixed + wti_copy( n) * self%v_bk_prev( tj,k)
        end do
        ! Relax solution to improve stability
        v_fixed = (C%visc_it_relax * v_fixed) + ((1._dp - C%visc_it_relax) * self%v_bk_prev( ti,k))
        ! Set load vector
        bb( row_tikuv) = v_fixed

      else
        call crash('unknown BC_u_north "' // trim( C%BC_u_north) // '"!')
      end if

    else
      call crash('uv can only be 1 or 2!')
    end if

  end subroutine calc_BPA_stiffness_matrix_row_BC_north

end submodule solve_linearised_BPA
