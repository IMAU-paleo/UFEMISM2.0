module SIA_main

  !< Routines for calculating ice velocities using the Shallow Ice Approximation (SIA)

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model, type_ice_velocity_solver_SIA
  use parameters, only: grav, ice_density
  use reallocate_mod, only: reallocate_bounds
  use constitutive_equation, only: calc_ice_rheology_Glen
  use mesh_disc_apply_operators, only: ddx_a_b_2D, ddy_a_b_2D, map_a_b_2D, map_a_b_3D, ddx_a_a_2D, ddy_a_a_2D
  use mesh_zeta, only: integrate_from_zeta_is_one_to_zeta_is_zetap

  implicit none

contains

  subroutine initialise_SIA_solver( mesh, SIA)
    !< Initialise the SIA solver

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_ice_velocity_solver_SIA), intent(  out) :: SIA

    ! Local variables:
    character(len=1024) :: routine_name = 'initialise_SIA_solver'

    ! Add routine to path
    call init_routine( routine_name)

    ! allocate memory
    call allocate_SIA_solver( mesh, SIA)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_SIA_solver

  subroutine allocate_SIA_solver( mesh, SIA)
    !< allocate memory for the SIA solver

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_ice_velocity_solver_SIA), intent(  out) :: SIA

    ! Local variables:
    character(len=1024) :: routine_name = 'allocate_SIA_solver'

    ! Add routine to path
    call init_routine( routine_name)

    allocate( SIA%u_3D_b  ( mesh%ti1:mesh%ti2, mesh%nz))
    allocate( SIA%v_3D_b  ( mesh%ti1:mesh%ti2, mesh%nz))
    allocate( SIA%du_dz_3D( mesh%vi1:mesh%vi2, mesh%nz))
    allocate( SIA%dv_dz_3D( mesh%vi1:mesh%vi2, mesh%nz))
    allocate( SIA%D_3D_b  ( mesh%ti1:mesh%ti2, mesh%nz))

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine allocate_SIA_solver

  subroutine solve_SIA( mesh, ice, SIA)
    !< Calculate ice velocities by solving the Shallow Ice Approximation with Glen's flow law

    ! In/output variables:
    type(type_mesh),                     intent(in   ) :: mesh
    type(type_ice_model),                intent(inout) :: ice
    type(type_ice_velocity_solver_SIA),  intent(inout) :: SIA

    ! Local variables:
    character(len=1024)                     :: routine_name = 'solve_SIA'
    real(dp), dimension(:    ), allocatable :: Hi_b, Hs_b, dHs_dx, dHs_dy, dHs_dx_b, dHs_dy_b
    real(dp), dimension(:,:  ), allocatable :: A_flow_b
    integer                                 :: vi,ti,k
    real(dp)                                :: abs_grad_Hs
    real(dp), dimension(mesh%nz)            :: z, int_A_hminzetan

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (.not. C%choice_flow_law == 'Glen') then
      call crash('the analytical solution to the SIA is only valid when using Glens flow law!')
    end if

    ! Allocate memory
    allocate( Hi_b(     mesh%ti1:mesh%ti2         ))
    allocate( Hs_b(     mesh%ti1:mesh%ti2         ))
    allocate( dHs_dx(   mesh%vi1:mesh%vi2         ))
    allocate( dHs_dy(   mesh%vi1:mesh%vi2         ))
    allocate( dHs_dx_b( mesh%ti1:mesh%ti2         ))
    allocate( dHs_dy_b( mesh%ti1:mesh%ti2         ))
    allocate( A_flow_b( mesh%ti1:mesh%ti2, mesh%nz))

    ! Calculate flow factors
    call calc_ice_rheology_Glen( mesh, ice)

    ! Calculate ice thickness, surface elevation, surface slopes, and ice flow factor on the b-grid
    call map_a_b_2D( mesh, ice%Hi    , Hi_b    )
    call map_a_b_2D( mesh, ice%Hs    , Hs_b    )
    call ddx_a_a_2D( mesh, ice%Hs    , dHs_dx  )
    call ddy_a_a_2D( mesh, ice%Hs    , dHs_dy  )
    call ddx_a_b_2D( mesh, ice%Hs    , dHs_dx_b)
    call ddy_a_b_2D( mesh, ice%Hs    , dHs_dy_b)
    call map_a_b_3D( mesh, ice%A_flow, A_flow_b)

    ! Calculate velocities and strain rates according to the analytical solution of the SIA:
    ! (see also Bueler and Brown, 2009, Eqs. 12-13)
    !
    !   D( z) = -2 (rho g)^n (abs(grad H))^(n-1) int_b_z( A(T*) (h - zeta)^n ) dzeta
    !   u( z) = dh/dx D( z)
    !   v( z) = dh/dy D( z)
    !
    !   du/dz( z) = -2 (rho g)^n (abs(grad H))^(n-1) A(T*) (h - z)^n dh/dx
    !   dv/dz( z) = -2 (rho g)^n (abs(grad H))^(n-1) A(T*) (h - z)^n dh/dy

    ! Calculate velocities
    do ti = mesh%ti1, mesh%ti2

      ! Calculate the integral from b to z of (A_flow * (h - zeta)^n) dzeta
      z = Hs_b( ti) - mesh%zeta * Hi_b( ti)
      int_A_hminzetan = integrate_from_zeta_is_one_to_zeta_is_zetap( z, A_flow_b( ti,:) * (Hs_b( ti) - z)**C%Glens_flow_law_exponent)

      ! Calculate the diffusivity term
      abs_grad_Hs = SQRT( dHs_dx_b( ti)**2 + dHs_dy_b( ti)**2)
      SIA%D_3D_b( ti,:) = -2._dp * (ice_density * grav)**C%Glens_flow_law_exponent * abs_grad_Hs**(C%Glens_flow_law_exponent - 1._dp) * int_A_hminzetan

      ! Safety
      SIA%D_3D_b( ti,:) = MAX( -C%SIA_maximum_diffusivity, SIA%D_3D_b( ti,:))

      ! Calculate the velocities
      SIA%u_3D_b( ti,:) = SIA%D_3D_b( ti,:) * dHs_dx_b( ti)
      SIA%v_3D_b( ti,:) = SIA%D_3D_b( ti,:) * dHs_dy_b( ti)

    end do

    ! Calculate vertical shear strain rates (needed later to calculate strain heating in thermodynamics)
    do vi = mesh%vi1, mesh%vi2

      abs_grad_Hs = SQRT( dHs_dx( vi)**2 + dHs_dy( vi)**2)
      z = ice%Hs( vi) - mesh%zeta * ice%Hi( vi)

      do k = 1, mesh%nz
        SIA%du_dz_3D( vi,k) = -2._dp * (ice_density * grav)**C%Glens_flow_law_exponent * abs_grad_Hs**(C%Glens_flow_law_exponent - 1._dp) * &
          ice%A_flow( vi,k) * (ice%Hs( vi) - z( k))**C%Glens_flow_law_exponent * dHs_dx( vi)
        SIA%dv_dz_3D( vi,k) = -2._dp * (ice_density * grav)**C%Glens_flow_law_exponent * abs_grad_Hs**(C%Glens_flow_law_exponent - 1._dp) * &
          ice%A_flow( vi,k) * (ice%Hs( vi) - z( k))**C%Glens_flow_law_exponent * dHs_dy( vi)
      end do

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine solve_SIA

  subroutine remap_SIA_solver( mesh_old, mesh_new, SIA)
    ! Remap the SIA solver

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh_old
    type(type_mesh),                    intent(in   ) :: mesh_new
    type(type_ice_velocity_solver_SIA), intent(inout) :: SIA

    ! Local variables:
    character(len=1024) :: routine_name = 'remap_SIA_solver'
    real(dp)            :: dp_dummy

    ! Add routine to path
    call init_routine( routine_name)

    ! To prevent compiler warnings
    dp_dummy = mesh_old%V( 1,1)
    dp_dummy = mesh_new%V( 1,1)

    ! Solution: reallocate
    call reallocate_bounds( SIA%u_3D_b  , mesh_new%ti1, mesh_new%ti2, mesh_new%nz)
    call reallocate_bounds( SIA%v_3D_b  , mesh_new%ti1, mesh_new%ti2, mesh_new%nz)
    call reallocate_bounds( SIA%du_dz_3D, mesh_new%vi1, mesh_new%vi2, mesh_new%nz)
    call reallocate_bounds( SIA%dv_dz_3D, mesh_new%vi1, mesh_new%vi2, mesh_new%nz)

    ! Intermediate data fields: reallocate
    call reallocate_bounds( SIA%D_3D_b  , mesh_new%ti1, mesh_new%ti2, mesh_new%nz)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remap_SIA_solver

end module SIA_main
