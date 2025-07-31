module mesh_disc_apply_operators

  ! Apply the matrix operators to calculate gradients of functions on the mesh.

  use assertions_basic, only: assert
  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use mesh_types, only: type_mesh
  use CSR_sparse_matrix_type, only: type_sparse_matrix_CSR_dp
  use mpi_distributed_memory, only: gather_to_all
  use CSR_matrix_vector_multiplication, only: multiply_CSR_matrix_with_vector_1D_wrapper, &
    multiply_CSR_matrix_with_vector_2D_wrapper

  implicit none

contains

  subroutine ddx_a_a_2D( mesh, d_a, ddx_a, d_a_is_hybrid, ddx_a_is_hybrid)
    ! ddx a 2-D data field from the a-grid (vertices) to the a-grid (vertices)

    ! In/output variables:
    type(type_mesh),            intent(in   ) :: mesh
    real(dp), dimension(:    ), intent(in   ) :: d_a
    real(dp), dimension(:    ), intent(  out) :: ddx_a
    logical, optional,          intent(in   ) :: d_a_is_hybrid, ddx_a_is_hybrid

    ! Local variables:
    character(len=256), parameter :: routine_name = 'ddx_a_a_2D'

    ! Add routine to path
    call init_routine( routine_name)

    ! Perform the ddxping operation as a matrix multiplication
    call multiply_CSR_matrix_with_vector_1D_wrapper( mesh%M_ddx_a_a, &
      mesh%pai_V, d_a, mesh%pai_V, ddx_a, &
      xx_is_hybrid = d_a_is_hybrid, yy_is_hybrid = ddx_a_is_hybrid, &
      buffer_xx_nih = mesh%buffer1_d_a_nih, buffer_yy_nih = mesh%buffer2_d_a_nih)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine ddx_a_a_2D

  subroutine ddx_a_a_3D( mesh, d_a, ddx_a, d_a_is_hybrid, ddx_a_is_hybrid)
    ! ddx a 3-D data field from the a-grid (vertices) to the a-grid (vertices)

    ! In/output variables:
    type(type_mesh),            intent(in   ) :: mesh
    real(dp), dimension(:,:  ), intent(in   ) :: d_a
    real(dp), dimension(:,:  ), intent(out  ) :: ddx_a
    logical, optional,          intent(in   ) :: d_a_is_hybrid, ddx_a_is_hybrid

    ! Local variables:
    character(len=256), parameter :: routine_name = 'ddx_a_a_3D'

    ! Add routine to path
    call init_routine( routine_name)

    ! Perform the ddxping operation as a matrix multiplication
    call multiply_CSR_matrix_with_vector_2D_wrapper( mesh%M_ddx_a_a, &
      mesh%pai_V, d_a, mesh%pai_V, ddx_a, &
      xx_is_hybrid = d_a_is_hybrid, yy_is_hybrid = ddx_a_is_hybrid, &
      buffer_xx_nih = mesh%buffer1_d_ak_nih, buffer_yy_nih = mesh%buffer2_d_ak_nih)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine ddx_a_a_3D

  subroutine ddy_a_a_2D( mesh, d_a, ddy_a, d_a_is_hybrid, ddy_a_is_hybrid)
    ! ddy a 2-D data field from the a-grid (vertices) to the a-grid (vertices)

    ! In/output variables:
    type(type_mesh),            intent(in   ) :: mesh
    real(dp), dimension(:    ), intent(in   ) :: d_a
    real(dp), dimension(:    ), intent(out  ) :: ddy_a
    logical, optional,          intent(in   ) :: d_a_is_hybrid, ddy_a_is_hybrid

    ! Local variables:
    character(len=256), parameter :: routine_name = 'ddy_a_a_2D'

    ! Add routine to path
    call init_routine( routine_name)

    ! Perform the ddyping operation as a matrix multiplication
    call multiply_CSR_matrix_with_vector_1D_wrapper( mesh%M_ddy_a_a, &
      mesh%pai_V, d_a, mesh%pai_V, ddy_a, &
      xx_is_hybrid = d_a_is_hybrid, yy_is_hybrid = ddy_a_is_hybrid, &
      buffer_xx_nih = mesh%buffer1_d_a_nih, buffer_yy_nih = mesh%buffer2_d_a_nih)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine ddy_a_a_2D

  subroutine ddy_a_a_3D( mesh, d_a, ddy_a, d_a_is_hybrid, ddy_a_is_hybrid)
    ! ddy a 3-D data field from the a-grid (vertices) to the a-grid (vertices)

    ! In/output variables:
    type(type_mesh),              intent(in   ) :: mesh
    real(dp), dimension(:,:  ),   intent(in   ) :: d_a
    real(dp), dimension(:,:  ),   intent(out  ) :: ddy_a
    logical, optional,            intent(in   ) :: d_a_is_hybrid, ddy_a_is_hybrid

    ! Local variables:
    character(len=256), parameter :: routine_name = 'ddy_a_a_3D'

    ! Add routine to path
    call init_routine( routine_name)

    call multiply_CSR_matrix_with_vector_2D_wrapper( mesh%M_ddy_a_a, &
      mesh%pai_V, d_a, mesh%pai_V, ddy_a, &
      xx_is_hybrid = d_a_is_hybrid, yy_is_hybrid = ddy_a_is_hybrid, &
      buffer_xx_nih = mesh%buffer1_d_ak_nih, buffer_yy_nih = mesh%buffer2_d_ak_nih)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine ddy_a_a_3D

  subroutine map_a_b_2D( mesh, d_a, d_b, d_a_is_hybrid, d_b_is_hybrid)
    ! Map a 2-D data field from the a-grid (vertices) to the b-grid (triangles)

    ! In/output variables:
    type(type_mesh),            intent(in   ) :: mesh
    real(dp), dimension(:    ), intent(in   ) :: d_a
    real(dp), dimension(:    ), intent(out  ) :: d_b
    logical, optional,          intent(in   ) :: d_a_is_hybrid, d_b_is_hybrid

    ! Local variables:
    character(len=256), parameter :: routine_name = 'map_a_b_2D'

    ! Add routine to path
    call init_routine( routine_name)

    ! Perform the mapping operation as a matrix multiplication
    call multiply_CSR_matrix_with_vector_1D_wrapper( mesh%M_map_a_b, &
      mesh%pai_V, d_a, mesh%pai_Tri, d_b, &
      xx_is_hybrid = d_a_is_hybrid, yy_is_hybrid = d_b_is_hybrid, &
      buffer_xx_nih = mesh%buffer1_d_a_nih, buffer_yy_nih = mesh%buffer2_d_b_nih)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_a_b_2D

  subroutine map_a_b_3D( mesh, d_a, d_b, d_a_is_hybrid, d_b_is_hybrid)
    ! Map a 3-D data field from the a-grid (vertices) to the b-grid (triangles)

    ! In/output variables:
    type(type_mesh),            intent(in   ) :: mesh
    real(dp), dimension(:,:  ), intent(in   ) :: d_a
    real(dp), dimension(:,:  ), intent(out  ) :: d_b
    logical, optional,          intent(in   ) :: d_a_is_hybrid, d_b_is_hybrid

    ! Local variables:
    character(len=256), parameter :: routine_name = 'map_a_b_3D'

    ! Add routine to path
    call init_routine( routine_name)

    ! Perform the mapping operation as a matrix multiplication
    call multiply_CSR_matrix_with_vector_2D_wrapper( mesh%M_map_a_b, &
      mesh%pai_V, d_a, mesh%pai_Tri, d_b, &
      xx_is_hybrid = d_a_is_hybrid, yy_is_hybrid = d_b_is_hybrid, &
      buffer_xx_nih = mesh%buffer1_d_ak_nih, buffer_yy_nih = mesh%buffer2_d_bk_nih)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_a_b_3D

  subroutine map_b_a_2D( mesh, d_b, d_a, d_b_is_hybrid, d_a_is_hybrid)
    ! Map a 2-D data field from the b-grid (triangles) to the a-grid (vertices)

    ! In/output variables:
    type(type_mesh),            intent(in   ) :: mesh
    real(dp), dimension(:    ), intent(in   ) :: d_b
    real(dp), dimension(:    ), intent(out  ) :: d_a
    logical, optional,          intent(in   ) :: d_b_is_hybrid, d_a_is_hybrid

    ! Local variables:
    character(len=256), parameter :: routine_name = 'map_b_a_2D'

    ! Add routine to path
    call init_routine( routine_name)

    ! Perform the mapping operation as a matrix multiplication
    call multiply_CSR_matrix_with_vector_1D_wrapper( mesh%M_map_b_a, &
      mesh%pai_Tri, d_b, mesh%pai_V, d_a, &
      xx_is_hybrid = d_b_is_hybrid, yy_is_hybrid = d_a_is_hybrid, &
      buffer_xx_nih = mesh%buffer1_d_b_nih, buffer_yy_nih = mesh%buffer2_d_a_nih)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_b_a_2D

  subroutine map_b_a_3D( mesh, d_b, d_a, d_b_is_hybrid, d_a_is_hybrid)
    ! Map a 3-D data field from the b-grid (triangles) to the a-grid (vertices)

    ! In/output variables:
    type(type_mesh),            intent(in   ) :: mesh
    real(dp), dimension(:,:  ), intent(in   ) :: d_b
    real(dp), dimension(:,:  ), intent(out  ) :: d_a
    logical, optional,          intent(in   ) :: d_b_is_hybrid, d_a_is_hybrid

    ! Local variables:
    character(len=256), parameter :: routine_name = 'map_b_a_3D'

    ! Add routine to path
    call init_routine( routine_name)

    ! Perform the mapping operation as a matrix multiplication
    call multiply_CSR_matrix_with_vector_2D_wrapper( mesh%M_map_b_a, &
      mesh%pai_Tri, d_b, mesh%pai_V, d_a, &
      xx_is_hybrid = d_b_is_hybrid, yy_is_hybrid = d_a_is_hybrid, &
      buffer_xx_nih = mesh%buffer1_d_bk_nih, buffer_yy_nih = mesh%buffer2_d_ak_nih)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_b_a_3D

  subroutine ddx_a_b_2D( mesh, d_a, ddx_b, d_a_is_hybrid, ddx_b_is_hybrid)
    ! ddx a 2-D data field from the a-grid (vertices) to the b-grid (triangles)

    ! In/output variables:
    type(type_mesh),            intent(in   ) :: mesh
    real(dp), dimension(:    ), intent(in   ) :: d_a
    real(dp), dimension(:    ), intent(out  ) :: ddx_b
    logical, optional,          intent(in   ) :: d_a_is_hybrid, ddx_b_is_hybrid

    ! Local variables:
    character(len=256), parameter :: routine_name = 'ddx_a_b_2D'

    ! Add routine to path
    call init_routine( routine_name)

    ! Perform the ddxping operation as a matrix multiplication
    call multiply_CSR_matrix_with_vector_1D_wrapper( mesh%M_ddx_a_b, &
      mesh%pai_V, d_a, mesh%pai_Tri, ddx_b, &
      xx_is_hybrid = d_a_is_hybrid, yy_is_hybrid = ddx_b_is_hybrid, &
      buffer_xx_nih = mesh%buffer1_d_a_nih, buffer_yy_nih = mesh%buffer2_d_b_nih)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine ddx_a_b_2D

  subroutine ddx_a_b_3D( mesh, d_a, ddx_b, d_a_is_hybrid, ddx_b_is_hybrid)
    ! ddx a 3-D data field from the a-grid (vertices) to the b-grid (triangles)

    ! In/output variables:
    type(type_mesh),            intent(in   ) :: mesh
    real(dp), dimension(:,:  ), intent(in   ) :: d_a
    real(dp), dimension(:,:  ), intent(out  ) :: ddx_b
    logical, optional,          intent(in   ) :: d_a_is_hybrid, ddx_b_is_hybrid

    ! Local variables:
    character(len=256), parameter :: routine_name = 'ddx_a_b_3D'

    ! Add routine to path
    call init_routine( routine_name)

    ! Perform the ddxping operation as a matrix multiplication
    call multiply_CSR_matrix_with_vector_2D_wrapper( mesh%M_ddx_a_b, &
      mesh%pai_V, d_a, mesh%pai_Tri, ddx_b, &
      xx_is_hybrid = d_a_is_hybrid, yy_is_hybrid = ddx_b_is_hybrid, &
      buffer_xx_nih = mesh%buffer1_d_ak_nih, buffer_yy_nih = mesh%buffer2_d_bk_nih)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine ddx_a_b_3D

  subroutine ddx_b_a_2D( mesh, d_b, ddx_a, d_b_is_hybrid, ddx_a_is_hybrid)
    ! ddx a 2-D data field from the b-grid (triangles) to the a-grid (vertices)

    ! In/output variables:
    type(type_mesh),            intent(in   ) :: mesh
    real(dp), dimension(:    ), intent(in   ) :: d_b
    real(dp), dimension(:    ), intent(out  ) :: ddx_a
    logical, optional,          intent(in   ) :: d_b_is_hybrid, ddx_a_is_hybrid

    ! Local variables:
    character(len=256), parameter :: routine_name = 'ddx_b_a_2D'

    ! Add routine to path
    call init_routine( routine_name)

    ! Perform the ddxping operation as a matrix multiplication
    call multiply_CSR_matrix_with_vector_1D_wrapper( mesh%M_ddx_b_a, &
      mesh%pai_Tri, d_b, mesh%pai_V, ddx_a, &
      xx_is_hybrid = d_b_is_hybrid, yy_is_hybrid = ddx_a_is_hybrid, &
      buffer_xx_nih = mesh%buffer1_d_b_nih, buffer_yy_nih = mesh%buffer2_d_a_nih)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine ddx_b_a_2D

  subroutine ddx_b_a_3D( mesh, d_b, ddx_a, d_b_is_hybrid, ddx_a_is_hybrid)
    ! ddx a 3-D data field from the b-grid (triangles) to the a-grid (vertices)

    ! In/output variables:
    type(type_mesh),            intent(in   ) :: mesh
    real(dp), dimension(:,:  ), intent(in   ) :: d_b
    real(dp), dimension(:,:  ), intent(out  ) :: ddx_a
    logical, optional,          intent(in   ) :: d_b_is_hybrid, ddx_a_is_hybrid

    ! Local variables:
    character(len=256), parameter :: routine_name = 'ddx_b_a_3D'

    ! Add routine to path
    call init_routine( routine_name)

    ! Perform the ddxping operation as a matrix multiplication
    call multiply_CSR_matrix_with_vector_2D_wrapper( mesh%M_ddx_b_a, &
      mesh%pai_Tri, d_b, mesh%pai_V, ddx_a, &
      xx_is_hybrid = d_b_is_hybrid, yy_is_hybrid = ddx_a_is_hybrid, &
      buffer_xx_nih = mesh%buffer1_d_bk_nih, buffer_yy_nih = mesh%buffer2_d_ak_nih)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine ddx_b_a_3D

  subroutine ddy_a_b_2D( mesh, d_a, ddy_b, d_a_is_hybrid, ddy_b_is_hybrid)
    ! ddy a 2-D data field from the a-grid (vertices) to the b-grid (triangles)

    ! In/output variables:
    type(type_mesh),            intent(in   ) :: mesh
    real(dp), dimension(:    ), intent(in   ) :: d_a
    real(dp), dimension(:    ), intent(out  ) :: ddy_b
    logical, optional,          intent(in   ) :: d_a_is_hybrid, ddy_b_is_hybrid

    ! Local variables:
    character(len=256), parameter :: routine_name = 'ddy_a_b_2D'

    ! Add routine to path
    call init_routine( routine_name)

    ! Perform the ddyping operation as a matrix multiplication
    call multiply_CSR_matrix_with_vector_1D_wrapper( mesh%M_ddy_a_b, &
      mesh%pai_V, d_a, mesh%pai_Tri, ddy_b, &
      xx_is_hybrid = d_a_is_hybrid, yy_is_hybrid = ddy_b_is_hybrid, &
      buffer_xx_nih = mesh%buffer1_d_a_nih, buffer_yy_nih = mesh%buffer2_d_b_nih)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine ddy_a_b_2D

  subroutine ddy_a_b_3D( mesh, d_a, ddy_b, d_a_is_hybrid, ddy_b_is_hybrid)
    ! ddy a 3-D data field from the a-grid (vertices) to the b-grid (triangles)

    ! In/output variables:
    type(type_mesh),            intent(in   ) :: mesh
    real(dp), dimension(:,:  ), intent(in   ) :: d_a
    real(dp), dimension(:,:  ), intent(out  ) :: ddy_b
    logical, optional,          intent(in   ) :: d_a_is_hybrid, ddy_b_is_hybrid

    ! Local variables:
    character(len=256), parameter :: routine_name = 'ddy_a_b_3D'

    ! Add routine to path
    call init_routine( routine_name)

    ! Perform the ddyping operation as a matrix multiplication
    call multiply_CSR_matrix_with_vector_2D_wrapper( mesh%M_ddy_a_b, &
      mesh%pai_V, d_a, mesh%pai_Tri, ddy_b, &
      xx_is_hybrid = d_a_is_hybrid, yy_is_hybrid = ddy_b_is_hybrid, &
      buffer_xx_nih = mesh%buffer1_d_ak_nih, buffer_yy_nih = mesh%buffer2_d_bk_nih)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine ddy_a_b_3D

  subroutine ddy_b_a_2D( mesh, d_b, ddy_a, d_b_is_hybrid, ddy_a_is_hybrid)
    ! ddy a 2-D data field from the b-grid (triangles) to the a-grid (vertices)

    ! In/output variables:
    type(type_mesh),            intent(in   ) :: mesh
    real(dp), dimension(:    ), intent(in   ) :: d_b
    real(dp), dimension(:    ), intent(out  ) :: ddy_a
    logical, optional,          intent(in   ) :: d_b_is_hybrid, ddy_a_is_hybrid

    ! Local variables:
    character(len=256), parameter :: routine_name = 'ddy_b_a_2D'

    ! Add routine to path
    call init_routine( routine_name)

    ! Perform the ddyping operation as a matrix multiplication
    call multiply_CSR_matrix_with_vector_1D_wrapper( mesh%M_ddy_b_a, &
      mesh%pai_Tri, d_b, mesh%pai_V, ddy_a, &
      xx_is_hybrid = d_b_is_hybrid, yy_is_hybrid = ddy_a_is_hybrid, &
      buffer_xx_nih = mesh%buffer1_d_b_nih, buffer_yy_nih = mesh%buffer2_d_a_nih)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine ddy_b_a_2D

  subroutine ddy_b_a_3D( mesh, d_b, ddy_a, d_b_is_hybrid, ddy_a_is_hybrid)
    ! ddy a 3-D data field from the b-grid (triangles) to the a-grid (vertices)

    ! In/output variables:
    type(type_mesh),            intent(in   ) :: mesh
    real(dp), dimension(:,:  ), intent(in   ) :: d_b
    real(dp), dimension(:,:  ), intent(out  ) :: ddy_a
    logical, optional,          intent(in   ) :: d_b_is_hybrid, ddy_a_is_hybrid

    ! Local variables:
    character(len=256), parameter :: routine_name = 'ddy_b_a_3D'

    ! Add routine to path
    call init_routine( routine_name)

    ! Perform the ddyping operation as a matrix multiplication
    call multiply_CSR_matrix_with_vector_2D_wrapper( mesh%M_ddy_b_a, &
      mesh%pai_Tri, d_b, mesh%pai_V, ddy_a, &
      xx_is_hybrid = d_b_is_hybrid, yy_is_hybrid = ddy_a_is_hybrid, &
      buffer_xx_nih = mesh%buffer1_d_bk_nih, buffer_yy_nih = mesh%buffer2_d_ak_nih)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine ddy_b_a_3D

  subroutine calc_3D_gradient_bk_ak( mesh, AA, d_bk, grad_d_ak)
    ! Apply a 3-D gradient operator to a 3-D data field

    ! In- and output variables:
    type(type_mesh),                                intent(in   ) :: mesh
    type(type_sparse_matrix_CSR_dp),                intent(in   ) :: AA
    real(dp), dimension(mesh%ti1:mesh%ti2,mesh%nz), intent(in   ) :: d_bk
    real(dp), dimension(mesh%vi1:mesh%vi2,mesh%nz), intent(out  ) :: grad_d_ak

    ! Local variables:
    character(len=256), parameter           :: routine_name = 'calc_3D_gradient_bk_ak'
    real(dp), dimension(:,:  ), allocatable :: d_bk_tot
    integer                                 :: vi,k,row_vik,ii1,ii2,ii,col_tikk,ti,kk
    real(dp)                                :: cAA

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety: check sizes
    if (AA%m_loc /= mesh%nV_loc   * mesh%nz .or. &
        AA%m     /= mesh%nV       * mesh%nz .or. &
        AA%n_loc /= mesh%nTri_loc * mesh%nz .or. &
        AA%n     /= mesh%nTri     * mesh%nz .or. &
        size(      d_bk,1) /= mesh%nTri_loc .or. size(      d_bk,2) /= mesh%nz .or. &
        size( grad_d_ak,1) /= mesh%nV_loc   .or. size( grad_d_ak,2) /= mesh%nz) then
      call crash('matrix and vector sizes dont match')
    end if

    ! allocate memory for gathered vector x
    allocate( d_bk_tot( mesh%nTri, mesh%nz))

    ! Gather data
    call gather_to_all( d_bk, d_bk_tot)

    ! Calculate gradient
    do vi = mesh%vi1, mesh%vi2
    do k  = 1, mesh%nz

      ! Vertex vi, layer k corresponds to this matrix row
      row_vik = mesh%vik2n( vi,k)

      ! Initialise
      grad_d_ak( vi,k) = 0._dp

      ! Loop over all contributing 3-D neighbours
      ii1 = AA%ptr( row_vik)
      ii2 = AA%ptr( row_vik+1) - 1

      do ii = ii1, ii2

        ! Read matrix coefficient
        col_tikk = AA%ind( ii)
        cAA      = AA%val( ii)

        ! This matrix column corresponds to triangle ti, layer kk
        ti = mesh%n2tik( col_tikk,1)
        kk = mesh%n2tik( col_tikk,2)

        ! Add contribution
        grad_d_ak( vi,k) = grad_d_ak( vi,k) + cAA * d_bk_tot( ti,kk)

      end do ! do ii = ii1, ii2

    end do ! do k  = 1, mesh%nz
    end do ! do vi = mesh%vi1, mesh%vi2

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_3D_gradient_bk_ak

  subroutine calc_3D_gradient_ak_bk( mesh, AA, d_ak, grad_d_bk)
    ! Apply a 3-D gradient operator to a 3-D data field

    ! In- and output variables:
    type(type_mesh),                                intent(in   ) :: mesh
    type(type_sparse_matrix_CSR_dp),                intent(in   ) :: AA
    real(dp), dimension(mesh%vi1:mesh%vi2,mesh%nz), intent(in   ) :: d_ak
    real(dp), dimension(mesh%ti1:mesh%ti2,mesh%nz), intent(out  ) :: grad_d_bk

    ! Local variables:
    character(len=256), parameter           :: routine_name = 'calc_3D_gradient_ak_bk'
    real(dp), dimension(:,:  ), allocatable :: d_ak_tot
    integer                                 :: ti,k,row_tik,ii1,ii2,ii,col_vikk,vi,kk
    real(dp)                                :: cAA

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety: check sizes
    if (AA%m_loc /= mesh%nTri_loc * mesh%nz .or. &
        AA%m     /= mesh%nTri     * mesh%nz .or. &
        AA%n_loc /= mesh%nV_loc   * mesh%nz .or. &
        AA%n     /= mesh%nV       * mesh%nz .or. &
        size(      d_ak,1) /= mesh%nV_loc   .or. size(      d_ak,2) /= mesh%nz .or. &
        size( grad_d_bk,1) /= mesh%nTri_loc .or. size( grad_d_bk,2) /= mesh%nz) then
      call crash('matrix and vector sizes dont match')
    end if

    ! allocate memory for gathered vector x
    allocate( d_ak_tot( mesh%nV, mesh%nz))

    ! Gather data
    call gather_to_all( d_ak, d_ak_tot)

    ! Calculate gradient
    do ti = mesh%ti1, mesh%ti2
    do k  = 1, mesh%nz

      ! Triangle ti, layer k corresponds to this matrix row
      row_tik = mesh%tik2n( ti,k)

      ! Initialise
      grad_d_bk( ti,k) = 0._dp

      ! Loop over all contributing 3-D neighbours
      ii1 = AA%ptr( row_tik)
      ii2 = AA%ptr( row_tik+1) - 1

      do ii = ii1, ii2

        ! Read matrix coefficient
        col_vikk = AA%ind( ii)
        cAA      = AA%val( ii)

        ! This matrix column corresponds to vertex vi, layer kk
        vi = mesh%n2vik( col_vikk,1)
        kk = mesh%n2vik( col_vikk,2)

        ! Add contribution
        grad_d_bk( ti,k) = grad_d_bk( ti,k) + cAA * d_ak_tot( vi,kk)

      end do ! do ii = ii1, ii2

    end do ! do k  = 1, mesh%nz
    end do ! do ti = mesh%ti1, mesh%ti2

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_3D_gradient_ak_bk

  subroutine calc_3D_gradient_bk_bks( mesh, AA, d_bk, grad_d_bks)
    ! Apply a 3-D gradient operator to a 3-D data field

    ! In- and output variables:
    type(type_mesh),                                  intent(in   ) :: mesh
    type(type_sparse_matrix_CSR_dp),                  intent(in   ) :: AA
    real(dp), dimension(mesh%ti1:mesh%ti2,mesh%nz),   intent(in   ) :: d_bk
    real(dp), dimension(mesh%ti1:mesh%ti2,mesh%nz-1), intent(out  ) :: grad_d_bks

    ! Local variables:
    character(len=256), parameter           :: routine_name = 'calc_3D_gradient_bk_bks'
    real(dp), dimension(:,:  ), allocatable :: d_bk_tot
    integer                                 :: ti,ks,row_tiks,ii1,ii2,ii,col_tjk,tj,k
    real(dp)                                :: cAA

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety: check sizes
    if (AA%m_loc /= mesh%nTri_loc * (mesh%nz-1) .or. &
        AA%m     /= mesh%nTri     * (mesh%nz-1) .or. &
        AA%n_loc /= mesh%nTri_loc *  mesh%nz    .or. &
        AA%n     /= mesh%nTri     *  mesh%nz    .or. &
        size(      d_bk ,1) /= mesh%nTri_loc .or. size(      d_bk ,2) /= mesh%nz .or. &
        size( grad_d_bks,1) /= mesh%nTri_loc .or. size( grad_d_bks,2) /= mesh%nz-1) then
      call crash('matrix and vector sizes dont match')
    end if

    ! allocate memory for gathered vector x
    allocate( d_bk_tot( mesh%nTri, mesh%nz))

    ! Gather data
    call gather_to_all( d_bk, d_bk_tot)

    ! Calculate gradient
    do ti = mesh%ti1, mesh%ti2
    do ks = 1, mesh%nz-1

      ! Triangle ti, staggered layer ks corresponds to this matrix row
      row_tiks = mesh%tiks2n( ti,ks)

      ! Initialise
      grad_d_bks( ti,ks) = 0._dp

      ! Loop over all contributing 3-D neighbours
      ii1 = AA%ptr( row_tiks)
      ii2 = AA%ptr( row_tiks+1) - 1

      do ii = ii1, ii2

        ! Read matrix coefficient
        col_tjk = AA%ind( ii)
        cAA     = AA%val( ii)

        ! This matrix column corresponds to triangle tj, layer k
        tj = mesh%n2tik( col_tjk,1)
        k  = mesh%n2tik( col_tjk,2)

        ! Add contribution
        grad_d_bks( ti,ks) = grad_d_bks( ti,ks) + cAA * d_bk_tot( tj,k)

      end do ! do ii = ii1, ii2

    end do ! do ks  = 1, mesh%nz-1
    end do ! do ti = mesh%ti1, mesh%ti2

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_3D_gradient_bk_bks

  subroutine calc_3D_gradient_bks_bk( mesh, AA, d_bks, grad_d_bk)
    ! Apply a 3-D gradient operator to a 3-D data field

    ! In- and output variables:
    type(type_mesh),                                  intent(in   ) :: mesh
    type(type_sparse_matrix_CSR_dp),                  intent(in   ) :: AA
    real(dp), dimension(mesh%ti1:mesh%ti2,mesh%nz-1), intent(in   ) :: d_bks
    real(dp), dimension(mesh%ti1:mesh%ti2,mesh%nz),   intent(out  ) :: grad_d_bk

    ! Local variables:
    character(len=256), parameter           :: routine_name = 'calc_3D_gradient_bks_bk'
    real(dp), dimension(:,:  ), allocatable :: d_bks_tot
    integer                                 :: ti,k,row_tik,ii1,ii2,ii,col_tjks,tj,ks
    real(dp)                                :: cAA

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety: check sizes
    if (AA%m_loc /= mesh%nTri_loc *  mesh%nz    .or. &
        AA%m     /= mesh%nTri     *  mesh%nz    .or. &
        AA%n_loc /= mesh%nTri_loc * (mesh%nz-1) .or. &
        AA%n     /= mesh%nTri     * (mesh%nz-1) .or. &
        size(      d_bks,1) /= mesh%nTri_loc .or. size(      d_bks,2) /= mesh%nz-1 .or. &
        size( grad_d_bk ,1) /= mesh%nTri_loc .or. size( grad_d_bk ,2) /= mesh%nz) then
      call crash('matrix and vector sizes dont match')
    end if

    ! allocate memory for gathered vector x
    allocate( d_bks_tot( mesh%nTri, mesh%nz-1))

    ! Gather data
    call gather_to_all( d_bks, d_bks_tot)

    ! Calculate gradient
    do ti = mesh%ti1, mesh%ti2
    do k = 1, mesh%nz

      ! Triangle ti, layer k corresponds to this matrix row
      row_tik = mesh%tik2n( ti,k)

      ! Initialise
      grad_d_bk( ti,k) = 0._dp

      ! Loop over all contributing 3-D neighbours
      ii1 = AA%ptr( row_tik)
      ii2 = AA%ptr( row_tik+1) - 1

      do ii = ii1, ii2

        ! Read matrix coefficient
        col_tjks = AA%ind( ii)
        cAA      = AA%val( ii)

        ! This matrix column corresponds to triangle tj, staggered layer ks
        tj = mesh%n2tiks( col_tjks,1)
        ks = mesh%n2tiks( col_tjks,2)

        ! Add contribution
        grad_d_bk( ti,k) = grad_d_bk( ti,k) + cAA * d_bks_tot( tj,ks)

      end do ! do ii = ii1, ii2

    end do ! do k  = 1, mesh%nz
    end do ! do ti = mesh%ti1, mesh%ti2

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_3D_gradient_bks_bk

  subroutine map_bks_ak( mesh, AA, d_bks, grad_d_ak)
    ! Apply a 3-D gradient operator to a 3-D data field

    ! In- and output variables:
    type(type_mesh),                     intent(in   )    :: mesh
    type(type_sparse_matrix_CSR_dp),     intent(in   )    :: AA
    real(dp), dimension(mesh%ti1:mesh%ti2,mesh%nz-1),        intent(in   )    :: d_bks
    real(dp), dimension(mesh%vi1:mesh%vi2,mesh%nz),          intent(out  )    :: grad_d_ak

    ! Local variables:
    character(len=256), parameter                      :: routine_name = 'map_bks_ak'
    real(dp), dimension(:,:  ), allocatable            :: d_bks_tot
    integer                                            :: vi,k,row_vik,ii1,ii2,ii,col_tiks,ti,ks
    real(dp)                                           :: cAA

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety: check sizes
    if (AA%m_loc /= mesh%nV_loc   *  mesh%nz    .or. &
        AA%m     /= mesh%nV       *  mesh%nz    .or. &
        AA%n_loc /= mesh%nTri_loc * (mesh%nz-1) .or. &
        AA%n     /= mesh%nTri     * (mesh%nz-1) .or. &
        size(      d_bks,1) /= mesh%nTri_loc .or. size(      d_bks,2) /= mesh%nz-1 .or. &
        size( grad_d_ak ,1) /= mesh%nV_loc   .or. size( grad_d_ak ,2) /= mesh%nz) then
      call crash('matrix and vector sizes dont match')
    end if

    ! allocate memory for gathered vector x
    allocate( d_bks_tot( mesh%nTri, mesh%nz-1))

    ! Gather data
    call gather_to_all( d_bks, d_bks_tot)

    ! Calculate gradient
    do vi = mesh%vi1, mesh%vi2
    do k = 1, mesh%nz

      ! Vertex vi, layer k corresponds to this matrix row
      row_vik = mesh%vik2n( vi,k)

      ! Initialise
      grad_d_ak( vi,k) = 0._dp

      ! Loop over all contributing 3-D neighbours
      ii1 = AA%ptr( row_vik)
      ii2 = AA%ptr( row_vik+1) - 1

      do ii = ii1, ii2

        ! Read matrix coefficient
        col_tiks = AA%ind( ii)
        cAA      = AA%val( ii)

        ! This matrix column corresponds to triangle ti, staggered layer ks
        ti = mesh%n2tiks( col_tiks,1)
        ks = mesh%n2tiks( col_tiks,2)

        ! Add contribution
        grad_d_ak( vi,k) = grad_d_ak( vi,k) + cAA * d_bks_tot( ti,ks)

      end do ! do ii = ii1, ii2

    end do ! do k  = 1, mesh%nz
    end do ! do vi = mesh%vi1, mesh%vi2

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_bks_ak

  subroutine map_ak_bks( mesh, AA, d_ak, grad_d_bks)
    ! Apply a 3-D gradient operator to a 3-D data field

    ! In- and output variables:
    type(type_mesh),                                  intent(in   ) :: mesh
    type(type_sparse_matrix_CSR_dp),                  intent(in   ) :: AA
    real(dp), dimension(mesh%vi1:mesh%vi2,mesh%nz),   intent(in   ) :: d_ak
    real(dp), dimension(mesh%ti1:mesh%ti2,mesh%nz-1), intent(out  ) :: grad_d_bks

    ! Local variables:
    character(len=256), parameter           :: routine_name = 'map_ak_bks'
    real(dp), dimension(:,:  ), allocatable :: d_ak_tot
    integer                                 :: ti,ks,row_tiks,ii1,ii2,ii,col_vik,vi,k
    real(dp)                                :: cAA

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety: check sizes
    if (AA%m_loc /= mesh%nTri_loc * (mesh%nz-1) .or. &
        AA%m     /= mesh%nTri     * (mesh%nz-1) .or. &
        AA%n_loc /= mesh%nV_loc   *  mesh%nz    .or. &
        AA%n     /= mesh%nV       *  mesh%nz    .or. &
        size(      d_ak ,1) /= mesh%nV_loc   .or. size(      d_ak ,2) /= mesh%nz .or. &
        size( grad_d_bks,1) /= mesh%nTri_loc .or. size( grad_d_bks,2) /= mesh%nz-1) then
      call crash('matrix and vector sizes dont match')
    end if

    ! allocate memory for gathered vector x
    allocate( d_ak_tot( mesh%nV, mesh%nz))

    ! Gather data
    call gather_to_all( d_ak, d_ak_tot)

    ! Calculate gradient
    do ti = mesh%ti1, mesh%ti2
    do ks = 1, mesh%nz-1

      ! Triangle ti, staggered layer ks corresponds to this matrix row
      row_tiks = mesh%tiks2n( ti,ks)

      ! Initialise
      grad_d_bks( ti,ks) = 0._dp

      ! Loop over all contributing 3-D neighbours
      ii1 = AA%ptr( row_tiks)
      ii2 = AA%ptr( row_tiks+1) - 1

      do ii = ii1, ii2

        ! Read matrix coefficient
        col_vik = AA%ind( ii)
        cAA     = AA%val( ii)

        ! This matrix column corresponds to vertex vi, layer k
        vi = mesh%n2vik( col_vik,1)
        k  = mesh%n2vik( col_vik,2)

        ! Add contribution
        grad_d_bks( ti,ks) = grad_d_bks( ti,ks) + cAA * d_ak_tot( vi,k)

      end do ! do ii = ii1, ii2

    end do ! do ks = 1, mesh%nz-1
    end do ! do ti = mesh%ti1, mesh%ti2

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_ak_bks

  subroutine calc_3D_gradient_bk_bk( mesh, AA, d_bk, grad_d_bk)
    ! Apply a 3-D gradient operator to a 3-D data field

    ! In- and output variables:
    type(type_mesh),                                intent(in   ) :: mesh
    type(type_sparse_matrix_CSR_dp),                intent(in   ) :: AA
    real(dp), dimension(mesh%ti1:mesh%ti2,mesh%nz), intent(in   ) :: d_bk
    real(dp), dimension(mesh%ti1:mesh%ti2,mesh%nz), intent(out  ) :: grad_d_bk

    ! Local variables:
    character(len=256), parameter           :: routine_name = 'calc_3D_gradient_bk_bk'
    real(dp), dimension(:,:  ), allocatable :: d_bk_tot
    integer                                 :: ti,k,row_tik,ii1,ii2,ii,col_tjkk,tj,kk
    real(dp)                                :: cAA

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety: check sizes
    if (AA%m_loc /= mesh%nTri_loc * mesh%nz .or. &
        AA%m     /= mesh%nTri     * mesh%nz .or. &
        AA%n_loc /= mesh%nTri_loc * mesh%nz .or. &
        AA%n     /= mesh%nTri     * mesh%nz .or. &
        size(      d_bk,1) /= mesh%nTri_loc .or. size(      d_bk,2) /= mesh%nz .or. &
        size( grad_d_bk,1) /= mesh%nTri_loc .or. size( grad_d_bk,2) /= mesh%nz) then
      call crash('matrix and vector sizes dont match')
    end if

    ! allocate memory for gathered vector x
    allocate( d_bk_tot( mesh%nTri, mesh%nz))

    ! Gather data
    call gather_to_all( d_bk, d_bk_tot)

    ! Calculate gradient
    do ti = mesh%ti1, mesh%ti2
    do k  = 1, mesh%nz

      ! Triangle ti, layer k corresponds to this matrix row
      row_tik = mesh%tik2n( ti,k)

      ! Initialise
      grad_d_bk( ti,k) = 0._dp

      ! Loop over all contributing 3-D neighbours
      ii1 = AA%ptr( row_tik)
      ii2 = AA%ptr( row_tik+1) - 1

      do ii = ii1, ii2

        ! Read matrix coefficient
        col_tjkk = AA%ind( ii)
        cAA      = AA%val( ii)

        ! This matrix column corresponds to triangle tj, layer kk
        tj = mesh%n2tik( col_tjkk,1)
        kk = mesh%n2tik( col_tjkk,2)

        ! Add contribution
        grad_d_bk( ti,k) = grad_d_bk( ti,k) + cAA * d_bk_tot( tj,kk)

      end do ! do ii = ii1, ii2

    end do ! do k  = 1, mesh%nz
    end do ! do ti = mesh%ti1, mesh%ti2

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_3D_gradient_bk_bk

end module mesh_disc_apply_operators
