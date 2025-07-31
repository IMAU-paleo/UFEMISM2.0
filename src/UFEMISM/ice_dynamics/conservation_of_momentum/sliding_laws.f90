module sliding_laws

  ! Contains all the different sliding laws.

  use precisions, only: dp
  use control_resources_and_error_messaging, only: crash, init_routine, finalise_routine
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model
  use bed_roughness_model_types, only: type_bed_roughness_model
  use parameters
  use mesh_disc_apply_operators, only: map_b_a_2D
  use mesh_utilities, only: extrapolate_Gaussian
  use mpi_distributed_memory, only: gather_to_all
  use Schoof_SSA_solution, only: Schoof2006_icestream

  implicit none

  private

  public :: calc_basal_friction_coefficient

contains

  subroutine calc_basal_friction_coefficient( mesh, ice, bed_roughness, u_b, v_b)
    ! Calculate the effective basal friction coefficient using the specified sliding law

    ! In/output variables:
    type(type_mesh),                intent(in   ) :: mesh
    type(type_ice_model),           intent(inout) :: ice
    type(type_bed_roughness_model), intent(in   ) :: bed_roughness
    real(dp), dimension(:),         intent(in   ) :: u_b
    real(dp), dimension(:),         intent(in   ) :: v_b

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'calc_basal_friction_coefficient'
    real(dp), dimension(mesh%vi1:mesh%vi2) :: u_a, v_a

    ! Add routine to path
    call init_routine( routine_name)

    ! Map velocities to the a-grid
    call map_b_a_2D( mesh, u_b, u_a)
    call map_b_a_2D( mesh, v_b, v_a)

    select case (C%choice_sliding_law)
    case default
      call crash('unknown choice_sliding_law "' // trim( C%choice_sliding_law) // '"')
    case ('no_sliding')
      ! No sliding allowed (implemented directly in linear equations of momentum balance; choice of beta is trivial)
      ice%basal_friction_coefficient = 0._dp
    case ('idealised')
      ! Sliding laws for some idealised experiments
      call calc_sliding_law_idealised(   mesh, ice, u_a, v_a)
    case ('Weertman')
      ! Weertman-type ("power law") sliding law
      call calc_sliding_law_Weertman(    mesh, ice, bed_roughness, u_a, v_a)
    case ('Coulomb')
      ! Coulomb-type sliding law
      call calc_sliding_law_Coulomb(     mesh, ice, bed_roughness, u_a, v_a)
    case ('Budd')
      ! Regularised Coulomb-type sliding law
      call calc_sliding_law_Budd(        mesh, ice, bed_roughness, u_a, v_a)
    case ('Tsai2015')
      ! Modified power-law relation according to Tsai et al. (2015)
      call calc_sliding_law_Tsai2015(    mesh, ice, bed_roughness, u_a, v_a)
    case ('Schoof2005')
      ! Modified power-law relation according to Schoof (2005)
      call calc_sliding_law_Schoof2005(  mesh, ice, bed_roughness, u_a, v_a)
    case ('Zoet-Iverson')
      ! Zoet-Iverson sliding law (Zoet & Iverson, 2020)
      call calc_sliding_law_ZoetIverson( mesh, ice, bed_roughness, u_a, v_a)
    end select

    ! Limit basal friction coefficient to improve stability
    ice%basal_friction_coefficient = min( C%slid_beta_max, ice%basal_friction_coefficient)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_basal_friction_coefficient

! == Different sliding laws
! =========================

  subroutine calc_sliding_law_Weertman( mesh, ice, bed_roughness, u_a, v_a)
    ! Weertman-type ("power law") sliding law

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    type(type_ice_model),                   intent(inout) :: ice
    type(type_bed_roughness_model),         intent(in   ) :: bed_roughness
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   ) :: u_a, v_a

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'calc_sliding_law_Weertman'
    real(dp), dimension(mesh%vi1:mesh%vi2) :: bed_roughness_applied
    integer                                :: vi
    real(dp)                               :: uabs

    ! Add routine to path
    call init_routine( routine_name)

    ! == Bed roughness
    ! ================

    ! Assume that bed roughness is represented by the beta_sq field

    ! Scale bed roughness based on grounded area fractions
    call apply_grounded_fractions_to_bed_roughness( mesh, ice, &
      bed_roughness%beta_sq, bed_roughness_applied)

    ! == Basal friction field
    ! =======================

    ! Calculate beta
    do vi = mesh%vi1, mesh%vi2

      ! Include a normalisation term following Bueler & Brown (2009) to prevent divide-by-zero errors.
      uabs = sqrt( C%slid_delta_v**2 + u_a( vi)**2 + v_a( vi)**2)

      ! Asay-Davis et al. (2016), Eq. 6; replacing beta_sq by the bed roughness field from above
      ice%basal_friction_coefficient( vi) = bed_roughness_applied( vi) * uabs ** (1._dp / C%slid_Weertman_m - 1._dp)

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_sliding_law_Weertman

  subroutine calc_sliding_law_Coulomb( mesh, ice, bed_roughness, u_a, v_a)
    ! Coulomb-type sliding law

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    type(type_ice_model),                   intent(inout) :: ice
    type(type_bed_roughness_model),         intent(in   ) :: bed_roughness
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   ) :: u_a, v_a

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'calc_sliding_law_Coulomb'
    real(dp), dimension(mesh%vi1:mesh%vi2) :: bed_roughness_applied
    integer                                :: vi
    real(dp)                               :: uabs

    ! Add routine to path
    call init_routine( routine_name)

    ! == Bed roughness
    ! ================

    ! Compute bed roughness based on till friction angle

    ! Scale bed roughness based on grounded area fractions
    call apply_grounded_fractions_to_bed_roughness( mesh, ice, &
      bed_roughness%till_friction_angle, bed_roughness_applied)

    ! == Till yield stress
    ! ====================

    ! Calculate the till yield stress from the effective pressure and bed roughness
    ice%till_yield_stress = ice%effective_pressure * tan(pi / 180._dp) * bed_roughness_applied

    ! == Extend till yield stress over ice-free land neighbours
    ! =========================================================

    call extend_till_yield_stress_to_neighbours( mesh, ice)

    ! == Basal friction field
    ! =======================

    ! Calculate beta
    do vi = mesh%vi1, mesh%vi2

      ! Include a normalisation term following Bueler & Brown (2009) to prevent divide-by-zero errors.
      uabs = sqrt( C%slid_delta_v**2 + u_a( vi)**2 + v_a( vi)**2)

      ice%basal_friction_coefficient( vi) = ice%till_yield_stress( vi) / uabs

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_sliding_law_Coulomb

  subroutine calc_sliding_law_Budd( mesh, ice, bed_roughness, u_a, v_a)
    ! Budd-type sliding law (formerly known as regularised Coulomb-type sliding law)

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    type(type_ice_model),                   intent(inout) :: ice
    type(type_bed_roughness_model),         intent(in   ) :: bed_roughness
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   ) :: u_a, v_a

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'calc_sliding_law_Budd'
    real(dp), dimension(mesh%vi1:mesh%vi2) :: bed_roughness_applied
    integer                                :: vi
    real(dp)                               :: uabs

    ! Add routine to path
    call init_routine( routine_name)

    ! == Bed roughness
    ! ================

    ! Compute bed roughness based on till friction angle

    ! Scale bed roughness based on grounded area fractions
    call apply_grounded_fractions_to_bed_roughness( mesh, ice, &
      bed_roughness%till_friction_angle, bed_roughness_applied)

    ! == Till yield stress
    ! ====================

    ! Calculate the till yield stress from the effective pressure and bed roughness
    ice%till_yield_stress = ice%effective_pressure * tan(pi / 180._dp) * bed_roughness_applied

    ! == Extend till yield stress over ice-free land neighbours
    ! =========================================================

    call extend_till_yield_stress_to_neighbours( mesh, ice)

    ! == Basal friction field
    ! =======================

    ! Calculate beta
    do vi = mesh%vi1, mesh%vi2

      ! Include a normalisation term following Bueler & Brown (2009) to prevent divide-by-zero errors.
      uabs = sqrt( C%slid_delta_v**2 + u_a( vi)**2 + v_a( vi)**2)

      ice%basal_friction_coefficient( vi) = ice%till_yield_stress( vi) * uabs ** (C%slid_Budd_q_plastic - 1._dp) / (C%slid_Budd_u_threshold ** C%slid_Budd_q_plastic)

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_sliding_law_Budd

  subroutine calc_sliding_law_Tsai2015(  mesh, ice, bed_roughness, u_a, v_a)
    ! Modified power-law relation according to Tsai et al. (2015)
    ! (implementation based on equations provided by Asay-Davis et al., 2016)
    !
    ! Asay-Davis et al.: Experimental design for three interrelated marine ice sheet and ocean model
    ! intercomparison projects: MISMIP v. 3 (MISMIP+), ISOMIP v. 2 (ISOMIP+) and MISOMIP v. 1 (MISOMIP1),
    ! Geoscientific Model Development 9, 2471-2497, 2016
    !
    ! Tsai et al.: Marine ice-sheet profiles and stability under Coulomb basal conditions,
    ! Journal of Glaciology 61, 205–215, doi:10.3189/2015JoG14J221, 2015.

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    type(type_ice_model),                   intent(inout) :: ice
    type(type_bed_roughness_model),         intent(in   ) :: bed_roughness
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   ) :: u_a, v_a

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'calc_sliding_law_Tsai2015'
    real(dp), dimension(mesh%vi1:mesh%vi2) :: bed_roughness_applied
    integer                                :: vi
    real(dp)                               :: uabs

    ! Add routine to path
    call init_routine( routine_name)

    ! == Bed roughness
    ! ================

    ! Assume that bed roughness is represented by the beta_sq field

    ! Scale bed roughness based on grounded area fractions
    call apply_grounded_fractions_to_bed_roughness( mesh, ice, &
      bed_roughness%beta_sq, bed_roughness_applied)

    ! == Basal friction field
    ! =======================

    ! Calculate beta
    do vi = mesh%vi1, mesh%vi2

      ! Include a normalisation term following Bueler & Brown (2009) to prevent divide-by-zero errors.
      uabs = sqrt( C%slid_delta_v**2 + u_a( vi)**2 + v_a( vi)**2)

      ! Asay-Davis et al. (2016), Eq. 7; replacing beta_sq by the bed roughness field from above
      ice%basal_friction_coefficient( vi) = min( bed_roughness%alpha_sq( vi) * ice%effective_pressure( vi), &
        bed_roughness_applied( vi) * uabs ** (1._dp / C%slid_Weertman_m)) * uabs**(-1._dp)

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_sliding_law_Tsai2015

  subroutine calc_sliding_law_Schoof2005(  mesh, ice, bed_roughness, u_a, v_a)
    ! Modified power-law relation according to Tsai et al. (2015)
    ! (implementation based on equations provided by Asay-Davis et al., 2016)
    !
    ! Asay-Davis et al.: Experimental design for three interrelated marine ice sheet and ocean model
    ! intercomparison projects: MISMIP v. 3 (MISMIP+), ISOMIP v. 2 (ISOMIP+) and MISOMIP v. 1 (MISOMIP1),
    ! Geoscientific Model Development 9, 2471-2497, 2016
    !
    ! Schoof: The effect of cavitation on glacier sliding, P. Roy. Soc. A-Math. Phy., 461, 609–627, doi:10.1098/rspa.2004.1350, 2005

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    type(type_ice_model),                   intent(inout) :: ice
    type(type_bed_roughness_model),         intent(in   ) :: bed_roughness
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   ) :: u_a, v_a

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'calc_sliding_law_Schoof2005'
    real(dp), dimension(mesh%vi1:mesh%vi2) :: bed_roughness_applied
    integer                                :: vi
    real(dp)                               :: uabs

    ! Add routine to path
    call init_routine( routine_name)

    ! == Bed roughness
    ! ================

    ! Assume that bed roughness is represented by the beta_sq field

    ! Scale bed roughness based on grounded area fractions
    call apply_grounded_fractions_to_bed_roughness( mesh, ice, &
      bed_roughness%beta_sq, bed_roughness_applied)

    ! == Basal friction field
    ! =======================

    ! Calculate beta
    do vi = mesh%vi1, mesh%vi2

      ! Include a normalisation term following Bueler & Brown (2009) to prevent divide-by-zero errors.
      uabs = sqrt( C%slid_delta_v**2 + u_a( vi)**2 + v_a( vi)**2)

      ! Asay-Davis et al. (2016), Eq. 11; replacing beta_sq by the bed roughness field from above
      ice%basal_friction_coefficient( vi) = ((bed_roughness_applied( vi) * uabs**(1._dp / C%slid_Weertman_m) * bed_roughness%alpha_sq( vi) * ice%effective_pressure( vi)) / &
        ((bed_roughness_applied( vi)**C%slid_Weertman_m * uabs + (bed_roughness%alpha_sq( vi) * ice%effective_pressure( vi))**C%slid_Weertman_m)**(1._dp / C%slid_Weertman_m))) * uabs**(-1._dp)

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_sliding_law_Schoof2005

  subroutine calc_sliding_law_ZoetIverson( mesh, ice, bed_roughness, u_a, v_a)
    ! Zoet-Iverson sliding law (Zoet & Iverson, 2020)

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    type(type_ice_model),                   intent(inout) :: ice
    type(type_bed_roughness_model),         intent(in   ) :: bed_roughness
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   ) :: u_a, v_a

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'calc_sliding_law_ZoetIverson'
    real(dp), dimension(mesh%vi1:mesh%vi2) :: bed_roughness_applied
    integer                                :: vi
    real(dp)                               :: uabs

    ! Add routine to path
    call init_routine( routine_name)

    ! == Bed roughness
    ! ================

    ! Compute bed roughness based on till friction angle

    ! Scale bed roughness based on grounded area fractions
    call apply_grounded_fractions_to_bed_roughness( mesh, ice, &
      bed_roughness%till_friction_angle, bed_roughness_applied)

    ! == Till yield stress
    ! ====================

    ! Calculate the till yield stress from the effective pressure and bed roughness
    ice%till_yield_stress = ice%effective_pressure * tan(pi / 180._dp) * bed_roughness_applied

    ! == Extend till yield stress over ice-free land neighbours
    ! =========================================================

    call extend_till_yield_stress_to_neighbours( mesh, ice)

    ! == Basal friction field
    ! =======================

    ! Calculate beta
    do vi = mesh%vi1, mesh%vi2

      ! Include a normalisation term following Bueler & Brown (2009) to prevent divide-by-zero errors.
      uabs = sqrt( C%slid_delta_v**2 + u_a( vi)**2 + v_a( vi)**2)

      ! Zoet & Iverson (2020), Eq. (3) (divided by u to give beta = tau_b / u)
      ice%basal_friction_coefficient( vi) = ice%till_yield_stress( vi) * &
        (uabs**(1._dp / C%slid_ZI_p - 1._dp)) * ((uabs + C%slid_ZI_ut)**(-1._dp / C%slid_ZI_p))

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_sliding_law_ZoetIverson

  subroutine calc_sliding_law_idealised(  mesh, ice, u_a, v_a)
    ! Sliding laws for some idealised experiments

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    type(type_ice_model),                   intent(inout) :: ice
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   ) :: u_a, v_a

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_sliding_law_idealised'

    ! Add routine to path
    call init_routine( routine_name)

    ! Use the specified idealised sliding law
    select case (C%choice_idealised_sliding_law)
    case default
      call crash('unknown choice_idealised_sliding_law "' // trim( C%choice_idealised_sliding_law) // '"')
    case ('SSA_icestream')
      call calc_sliding_law_idealised_SSA_icestream( mesh, ice, u_a, v_a)
    case ('ISMIP-HOM_C')
      ! ISMIP-HOM experiment C
      call calc_sliding_law_idealised_ISMIP_HOM_C( mesh, ice)
    case ('ISMIP-HOM_D')
      ! ISMIP-HOM experiment D
      call calc_sliding_law_idealised_ISMIP_HOM_D( mesh, ice)
    case ('ISMIP-HOM_E')
      ! ISMIP-HOM experiment E
      call crash('the Glacier Arolla experiment is not implemented in UFEMISM')
    case ('ISMIP-HOM_F')
      ! ISMIP-HOM experiment F
      call calc_sliding_law_idealised_ISMIP_HOM_F( mesh, ice)
    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_sliding_law_idealised

  subroutine calc_sliding_law_idealised_SSA_icestream( mesh, ice, u_a, v_a)
    ! Sliding laws for some idealised experiments
    !
    ! SSA_icestream

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    type(type_ice_model),                   intent(inout) :: ice
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   ) :: u_a, v_a

    ! Local variables:
    character(len=1024), parameter:: routine_name = 'calc_sliding_law_idealised_SSA_icestream'
    integer                       :: vi
    real(dp)                      :: y, u, uabs

    ! Add routine to path
    call init_routine( routine_name)

    do vi = mesh%vi1, mesh%vi2

      y = mesh%V( vi,2)
      call Schoof2006_icestream( C%uniform_Glens_flow_factor, C%Glens_flow_law_exponent, C%refgeo_idealised_SSA_icestream_Hi, &
        C%refgeo_idealised_SSA_icestream_dhdx, C%refgeo_idealised_SSA_icestream_L, C%refgeo_idealised_SSA_icestream_m, &
        y, u, ice%till_yield_stress( vi))

      ! Include a normalisation term following Bueler & Brown (2009) to prevent divide-by-zero errors.
      uabs = sqrt( C%slid_delta_v**2 + u_a( vi)**2 + v_a( vi)**2)

      ice%basal_friction_coefficient( vi) = ice%till_yield_stress( vi) / uabs

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_sliding_law_idealised_SSA_icestream

  subroutine calc_sliding_law_idealised_ISMIP_HOM_C( mesh, ice)
    ! Sliding laws for some idealised experiments
    !
    ! ISMIP-HOM experiment C

    ! In/output variables:
    type(type_mesh),      intent(in   ) :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=1024), parameter                     :: routine_name = 'calc_sliding_law_idealised_ISMIP_HOM_C'
    integer                                            :: vi
    real(dp)                                           :: x,y

    ! Add routine to path
    call init_routine( routine_name)

    do vi = mesh%vi1, mesh%vi2
      x = mesh%V( vi,1)
      y = mesh%V( vi,2)
      ice%basal_friction_coefficient( vi) = 1000._dp + &
        1000._dp * sin( 2._dp * pi * x / C%refgeo_idealised_ISMIP_HOM_L) * &
                   sin( 2._dp * pi * y / C%refgeo_idealised_ISMIP_HOM_L)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_sliding_law_idealised_ISMIP_HOM_C

  subroutine calc_sliding_law_idealised_ISMIP_HOM_D( mesh, ice)
    ! Sliding laws for some idealised experiments
    !
    ! ISMIP-HOM experiment D

    ! In/output variables:
    type(type_mesh),      intent(in   ) :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_sliding_law_idealised_ISMIP_HOM_D'
    integer                        :: vi
    real(dp)                       :: x

    ! Add routine to path
    call init_routine( routine_name)

    do vi = mesh%vi1, mesh%vi2
      x = mesh%V( vi,1)
      ice%basal_friction_coefficient( vi) = 1000._dp + &
        1000._dp * sin( 2._dp * pi * x / C%refgeo_idealised_ISMIP_HOM_L)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_sliding_law_idealised_ISMIP_HOM_D

  subroutine calc_sliding_law_idealised_ISMIP_HOM_F( mesh, ice)
    ! Sliding laws for some idealised experiments
    !
    ! ISMIP-HOM experiment F

    ! In/output variables:
    type(type_mesh),      intent(in   ) :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_sliding_law_idealised_ISMIP_HOM_F'
    integer                        :: vi

    ! Add routine to path
    call init_routine( routine_name)

    do vi = mesh%vi1, mesh%vi2
      ice%basal_friction_coefficient( vi) = (C%uniform_Glens_flow_factor * 1000._dp)**(-1._dp)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_sliding_law_idealised_ISMIP_HOM_F

! == Utilities
! ============

  subroutine apply_grounded_fractions_to_bed_roughness( mesh, ice, bed_roughness, bed_roughness_applied)
    ! Scale bed roughness based on grounded area fractions

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    type(type_ice_model),                   intent(in   ) :: ice
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   ) :: bed_roughness
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(  out) :: bed_roughness_applied

    ! Local variables:
    character(len=1024), parameter:: routine_name = 'apply_grounded_fractions_to_bed_roughness'
    integer                       :: vi
    real(dp)                      :: weight_gr, exponent_hi, exponent_hs, exponent_gr

    ! Add routine to path
    call init_routine( routine_name)

    ! Check if this is actually wanted
    if (.not. C%do_subgrid_friction_on_A_grid) then
      bed_roughness_applied = bed_roughness
      call finalise_routine( routine_name)
      return
    end if

    do vi = mesh%vi1, mesh%vi2

      ! Initialise grounded area fraction weight
      weight_gr = 1._dp

      ! Compute exponent for this vertex's weight based on ice thickness
      exponent_hi = log10( max( 1._dp, ice%Hi( vi)))
      ! Compute exponent for this vertex's weight based on surface gradients
      exponent_hs = ice%Hs_slope( vi) / 0.005_dp
      ! Compute final exponent for this vertex's weight
      exponent_gr = max( 0._dp, exponent_hi - exponent_hs)

      ! Compute a weight based on the grounded area fractions
      if (ice%mask_gl_gr( vi)) then
        weight_gr = ice%fraction_gr( vi)**exponent_gr

      elseif (ice%mask_cf_gr( vi)) then
        weight_gr = ice%fraction_gr( vi)**exponent_gr

      elseif (ice%mask_gl_fl( vi)) then
        weight_gr = ice%fraction_gr( vi)**exponent_gr

      elseif (ice%mask_grounded_ice( vi)) then
        weight_gr = 1._dp

      elseif (ice%mask_floating_ice( vi)) then
        weight_gr = 0._dp

      elseif (ice%mask_icefree_ocean( vi)) then
        weight_gr = 0._dp

      end if

      ! Just in case
      weight_gr = min( 1._dp, max( 0._dp, weight_gr))

      ! Compute till friction angle accounting for grounded area fractions
      bed_roughness_applied( vi) = bed_roughness(vi) * weight_gr

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine apply_grounded_fractions_to_bed_roughness

  subroutine extend_till_yield_stress_to_neighbours( mesh, ice)
    ! Extend till yield stress over ice-free land neighbours

    ! In/output variables:
    type(type_mesh),      intent(in   ) :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'extend_till_yield_stress_to_neighbours'
    integer                        :: vi, ci, vc
    real(dp)                       :: min_neighbour
    logical,  dimension(mesh%nV)   :: mask_grounded_ice_tot
    real(dp), dimension(mesh%nV)   :: till_yield_stress_tot
    logical                        :: found_grounded_neighbour

    ! Add routine to path
    call init_routine( routine_name)

    ! Gather data from all processes
    call gather_to_all( ice%mask_grounded_ice, mask_grounded_ice_tot)
    call gather_to_all( ice%till_yield_stress, till_yield_stress_tot)

    do vi = mesh%vi1, mesh%vi2

      ! Skip if not ice-free land
      if (.not. ice%mask_icefree_land( vi)) cycle

      ! Initialise
      found_grounded_neighbour = .false.

      ! Initialise with default values and assume a till yield
      ! stress equivalent to a column of ice of 1000 metres
      ! with no hydrology and maximum bed roughness.
      min_neighbour = 1000._dp * ice_density * grav

      ! Check grounded neighbours
      do ci = 1, mesh%nC( vi)
        vc = mesh%C( vi, ci)
        if (mask_grounded_ice_tot( vc)) then
          min_neighbour = min( min_neighbour, till_yield_stress_tot( vc))
          found_grounded_neighbour = .true.
        end if
      end do

      if (found_grounded_neighbour) then
        ! Use the minimum value among neighbours
        ice%till_yield_stress( vi) = min_neighbour
      else
        ! Use a default minimum value to avoid 0 friction
        ice%till_yield_stress( vi) = C%Hi_min * ice_density * grav
      end if

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine extend_till_yield_stress_to_neighbours

end module sliding_laws
