module zeta_gradients
  !< Calculate different terms in the Jacobian of the zeta coordinate transformaton

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model
  use mesh_disc_apply_operators, only: ddx_a_a_2D, ddy_a_a_2D, map_a_b_2D, ddx_a_b_2D, ddy_a_b_2D, &
    ddx_b_a_2D, ddy_b_a_2D

  implicit none

  private

  public :: calc_zeta_gradients

contains

  subroutine calc_zeta_gradients( mesh, ice)
    !< Calculate all the gradients of zeta,
    !< needed to perform the scaled vertical coordinate transformation

    ! In/output variables:
    type(type_mesh),      intent(in   ) :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=1024), parameter    :: routine_name = 'calc_zeta_gradients'
    ! real(dp), dimension(:), pointer :: Hi_a
    real(dp), dimension(:), pointer   :: dHi_dx_a
    real(dp), dimension(:), pointer   :: dHi_dy_a
    real(dp), dimension(:), pointer   :: d2Hi_dx2_a
    real(dp), dimension(:), pointer   :: d2Hi_dxdy_a
    real(dp), dimension(:), pointer   :: d2Hi_dy2_a
    real(dp), dimension(:), pointer   :: Hi_b
    real(dp), dimension(:), pointer   :: dHi_dx_b
    real(dp), dimension(:), pointer   :: dHi_dy_b
    real(dp), dimension(:), pointer   :: d2Hi_dx2_b
    real(dp), dimension(:), pointer   :: d2Hi_dxdy_b
    real(dp), dimension(:), pointer   :: d2Hi_dy2_b
    real(dp), dimension(:), pointer   :: Hs_b
    ! real(dp), dimension(:), pointer :: Hs_a
    real(dp), dimension(:), pointer   :: dHs_dx_a
    real(dp), dimension(:), pointer   :: dHs_dy_a
    real(dp), dimension(:), pointer   :: d2Hs_dx2_a
    real(dp), dimension(:), pointer   :: d2Hs_dxdy_a
    real(dp), dimension(:), pointer   :: d2Hs_dy2_a
    real(dp), dimension(:), pointer   :: dHs_dx_b
    real(dp), dimension(:), pointer   :: dHs_dy_b
    real(dp), dimension(:), pointer   :: d2Hs_dx2_b
    real(dp), dimension(:), pointer   :: d2Hs_dxdy_b
    real(dp), dimension(:), pointer   :: d2Hs_dy2_b
    integer                           :: vi,ti,k,ks
    real(dp)                          :: Hi, zeta

    ! Add routine to path
    call init_routine( routine_name)

    ! allocate shared memory
    ! allocate( Hi_a(        mesh%vi1:mesh%vi2))
    allocate( dHi_dx_a(    mesh%vi1:mesh%vi2))
    allocate( dHi_dy_a(    mesh%vi1:mesh%vi2))
    allocate( d2Hi_dx2_a(  mesh%vi1:mesh%vi2))
    allocate( d2Hi_dxdy_a( mesh%vi1:mesh%vi2))
    allocate( d2Hi_dy2_a(  mesh%vi1:mesh%vi2))

    allocate( Hi_b(        mesh%ti1:mesh%ti2))
    allocate( dHi_dx_b(    mesh%ti1:mesh%ti2))
    allocate( dHi_dy_b(    mesh%ti1:mesh%ti2))
    allocate( d2Hi_dx2_b(  mesh%ti1:mesh%ti2))
    allocate( d2Hi_dxdy_b( mesh%ti1:mesh%ti2))
    allocate( d2Hi_dy2_b(  mesh%ti1:mesh%ti2))

    ! allocate( Hs_a(        mesh%vi1:mesh%vi2))
    allocate( dHs_dx_a(    mesh%vi1:mesh%vi2))
    allocate( dHs_dy_a(    mesh%vi1:mesh%vi2))
    allocate( d2Hs_dx2_a(  mesh%vi1:mesh%vi2))
    allocate( d2Hs_dxdy_a( mesh%vi1:mesh%vi2))
    allocate( d2Hs_dy2_a(  mesh%vi1:mesh%vi2))

    allocate( Hs_b(        mesh%ti1:mesh%ti2))
    allocate( dHs_dx_b(    mesh%ti1:mesh%ti2))
    allocate( dHs_dy_b(    mesh%ti1:mesh%ti2))
    allocate( d2Hs_dx2_b(  mesh%ti1:mesh%ti2))
    allocate( d2Hs_dxdy_b( mesh%ti1:mesh%ti2))
    allocate( d2Hs_dy2_b(  mesh%ti1:mesh%ti2))

    ! Calculate gradients of Hi and Hs on both grids

    ! call map_a_a_2D( mesh, ice%Hi_a, Hi_a       )
    call ddx_a_a_2D( mesh, ice%Hi  , dHi_dx_a   )
    call ddy_a_a_2D( mesh, ice%Hi  , dHi_dy_a   )
    call map_a_b_2D( mesh, ice%Hi  , Hi_b       )
    call ddx_a_b_2D( mesh, ice%Hi  , dHi_dx_b   )
    call ddy_a_b_2D( mesh, ice%Hi  , dHi_dy_b   )
    call ddx_b_a_2D( mesh, dHi_dx_b, d2Hi_dx2_a )
    call ddy_b_a_2D( mesh, dHi_dx_b, d2Hi_dxdy_a)
    call ddy_b_a_2D( mesh, dHi_dy_b, d2Hi_dy2_a )
    call ddx_a_b_2D( mesh, dHi_dx_a, d2Hi_dx2_b )
    call ddy_a_b_2D( mesh, dHi_dx_a, d2Hi_dxdy_b)
    call ddy_a_b_2D( mesh, dHi_dy_a, d2Hi_dy2_b )

    ! call map_a_a_2D( mesh, ice%Hs_a, Hs_a       )
    call ddx_a_a_2D( mesh, ice%Hs  , dHs_dx_a   )
    call ddy_a_a_2D( mesh, ice%Hs  , dHs_dy_a   )
    call map_a_b_2D( mesh, ice%Hs  , Hs_b       )
    call ddx_a_b_2D( mesh, ice%Hs  , dHs_dx_b   )
    call ddy_a_b_2D( mesh, ice%Hs  , dHs_dy_b   )
    call ddx_b_a_2D( mesh, dHs_dx_b, d2Hs_dx2_a )
    call ddy_b_a_2D( mesh, dHs_dx_b, d2Hs_dxdy_a)
    call ddy_b_a_2D( mesh, dHs_dy_b, d2Hs_dy2_a )
    call ddx_a_b_2D( mesh, dHs_dx_a, d2Hs_dx2_b )
    call ddy_a_b_2D( mesh, dHs_dx_a, d2Hs_dxdy_b)
    call ddy_a_b_2D( mesh, dHs_dy_a, d2Hs_dy2_b )

    ! Calculate zeta gradients on all grids

    ! ak
    do vi = mesh%vi1, mesh%vi2

      Hi = max( 10._dp, ice%Hi( vi))

      do k = 1, mesh%nz

        zeta = mesh%zeta( k)

        ice%dzeta_dt_ak(    vi,k) = ( 1._dp / Hi) * (ice%dHs_dt( vi) - zeta * ice%dHi_dt( vi))

        ice%dzeta_dx_ak(    vi,k) = ( 1._dp / Hi) * (dHs_dx_a( vi) - zeta * dHi_dx_a( vi))
        ice%dzeta_dy_ak(    vi,k) = ( 1._dp / Hi) * (dHs_dy_a( vi) - zeta * dHi_dy_a( vi))
        ice%dzeta_dz_ak(    vi,k) = (-1._dp / Hi)

        ice%d2zeta_dx2_ak(  vi,k) = (dHi_dx_a( vi) * -1._dp / Hi) * ice%dzeta_dx_ak( vi,k) + (1._dp / Hi) * (d2Hs_dx2_a(  vi) - zeta * d2Hi_dx2_a(  vi))
        ice%d2zeta_dxdy_ak( vi,k) = (dHi_dy_a( vi) * -1._dp / Hi) * ice%dzeta_dx_ak( vi,k) + (1._dp / Hi) * (d2Hs_dxdy_a( vi) - zeta * d2Hi_dxdy_a( vi))
        ice%d2zeta_dy2_ak(  vi,k) = (dHi_dy_a( vi) * -1._dp / Hi) * ice%dzeta_dy_ak( vi,k) + (1._dp / Hi) * (d2Hs_dy2_a(  vi) - zeta * d2Hi_dy2_a(  vi))

      end do
    end do

    ! bk
    do ti = mesh%ti1, mesh%ti2

      Hi = max( 10._dp, Hi_b( ti))

      do k = 1, mesh%nz

        zeta = mesh%zeta( k)

        ice%dzeta_dx_bk(    ti,k) = ( 1._dp / Hi) * (dHs_dx_b( ti) - zeta * dHi_dx_b( ti))
        ice%dzeta_dy_bk(    ti,k) = ( 1._dp / Hi) * (dHs_dy_b( ti) - zeta * dHi_dy_b( ti))
        ice%dzeta_dz_bk(    ti,k) = (-1._dp / Hi)

        ice%d2zeta_dx2_bk(  ti,k) = (dHi_dx_b( ti) * -1._dp / Hi) * ice%dzeta_dx_bk( ti,k) + (1._dp / Hi) * (d2Hs_dx2_b(  ti) - zeta * d2Hi_dx2_b(  ti))
        ice%d2zeta_dxdy_bk( ti,k) = (dHi_dy_b( ti) * -1._dp / Hi) * ice%dzeta_dx_bk( ti,k) + (1._dp / Hi) * (d2Hs_dxdy_b( ti) - zeta * d2Hi_dxdy_b( ti))
        ice%d2zeta_dy2_bk(  ti,k) = (dHi_dy_b( ti) * -1._dp / Hi) * ice%dzeta_dy_bk( ti,k) + (1._dp / Hi) * (d2Hs_dy2_b(  ti) - zeta * d2Hi_dy2_b(  ti))

      end do
    end do

    ! bks
    do ti = mesh%ti1, mesh%ti2

      Hi = max( 10._dp, Hi_b( ti))

      do ks = 1, mesh%nz-1

        zeta = mesh%zeta_stag( ks)

        ice%dzeta_dx_bks(    ti,ks) = ( 1._dp / Hi) * (dHs_dx_b( ti) - zeta * dHi_dx_b( ti))
        ice%dzeta_dy_bks(    ti,ks) = ( 1._dp / Hi) * (dHs_dy_b( ti) - zeta * dHi_dy_b( ti))
        ice%dzeta_dz_bks(    ti,ks) = (-1._dp / Hi)

        ice%d2zeta_dx2_bks(  ti,ks) = (dHi_dx_b( ti) * -1._dp / Hi) * ice%dzeta_dx_bks( ti,ks) + (1._dp / Hi) * (d2Hs_dx2_b(  ti) - zeta * d2Hi_dx2_b(  ti))
        ice%d2zeta_dxdy_bks( ti,ks) = (dHi_dy_b( ti) * -1._dp / Hi) * ice%dzeta_dx_bks( ti,ks) + (1._dp / Hi) * (d2Hs_dxdy_b( ti) - zeta * d2Hi_dxdy_b( ti))
        ice%d2zeta_dy2_bks(  ti,ks) = (dHi_dy_b( ti) * -1._dp / Hi) * ice%dzeta_dy_bks( ti,ks) + (1._dp / Hi) * (d2Hs_dy2_b(  ti) - zeta * d2Hi_dy2_b(  ti))

      end do
    end do

    ! Clean after yourself
    ! deallocate( Hi_a        )
    deallocate( dHi_dx_a    )
    deallocate( dHi_dy_a    )
    deallocate( d2Hi_dx2_a  )
    deallocate( d2Hi_dxdy_a )
    deallocate( d2Hi_dy2_a  )

    deallocate( Hi_b        )
    deallocate( dHi_dx_b    )
    deallocate( dHi_dy_b    )
    deallocate( d2Hi_dx2_b  )
    deallocate( d2Hi_dxdy_b )
    deallocate( d2Hi_dy2_b  )

    ! deallocate( Hs_a        )
    deallocate( dHs_dx_a    )
    deallocate( dHs_dy_a    )
    deallocate( d2Hs_dx2_a  )
    deallocate( d2Hs_dxdy_a )
    deallocate( d2Hs_dy2_a  )

    deallocate( Hs_b        )
    deallocate( dHs_dx_b    )
    deallocate( dHs_dy_b    )
    deallocate( d2Hs_dx2_b  )
    deallocate( d2Hs_dxdy_b )
    deallocate( d2Hs_dy2_b  )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_zeta_gradients

end module zeta_gradients
