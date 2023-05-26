MODULE basal_conditions_main

  ! The main ice-dynamical model module.

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE control_resources_and_error_messaging                  , ONLY: crash, init_routine, finalise_routine
  USE model_configuration                                    , ONLY: C
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model

  IMPLICIT NONE

CONTAINS

! ===== Main routines =====
! =========================

  subroutine initialise_basal_conditions( mesh, ice)
    ! Allocation and initialisation

    implicit none

    ! Input variables:
    type(type_mesh),      intent(in)    :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=256), parameter       :: routine_name = 'initialise_basal_conditions'

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    call init_routine( routine_name)

    if (par%master) then
      write(*,"(A)") '   Initialising basal conditions...'
    end if
    call sync

    ! === Hydrology ===
    ! =================

    call initialise_basal_hydrology( mesh, ice)

    ! === Bed roughness ===
    ! =====================

    call initialise_bed_roughness( mesh, ice)

    ! === Geothermal heat ===
    ! =======================

    call initialise_geothermal_heat_flux( mesh, ice)

    ! === Finalisation ===
    ! ====================

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_basal_conditions

! ===== Hydrology =====
! =====================

  subroutine initialise_basal_hydrology( mesh, ice)
    ! Initialise basal hydrology

    implicit none

    ! Input variables:
    type(type_mesh),      intent(in)    :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=256), parameter       :: routine_name = 'initialise_basal_hydrology'

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    call init_routine( routine_name)

    ! === Stuff ===
    ! =============

    ! Nothing needed for now

    ! === Finalisation ===
    ! ====================

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_basal_hydrology

! ===== Bed roughness =====
! =========================

  subroutine initialise_bed_roughness( mesh, ice)
    ! Initialise basal friction coefficients

    implicit none

    ! Input variables:
    type(type_mesh),      intent(in)    :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=256), parameter       :: routine_name = 'initialise_bed_roughness'

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    call init_routine( routine_name)

    ! === Sliding law choice ===
    ! ==========================

    select case (C%choice_sliding_law)

      case ('no_sliding','idealised')
        ! Nothing to do; return

        call finalise_routine( routine_name)
        return

      case ('Weertman')
        ! Weertman-type sliding law

        ! Initialise with uniform friction coefficient
        ice%friction_coef_1 = C%slid_Weertman_beta_sq_uniform

      case ('Coulomb')
        ! Coulomb-type sliding law

        ! Initialise with uniform friction coefficient
        ice%friction_coef_1 = C%slid_Coulomb_phi_fric_uniform

      case ('Budd')
        ! Regularised Coulomb-type sliding law

        ! Initialise with uniform friction coefficient
        ice%friction_coef_1 = C%slid_Budd_phi_fric_uniform

      case ('Schoof2005')
        ! Schoof (2005) sliding law

        ! Initialise with uniform friction coefficients
        ice%friction_coef_1 = C%slid_Schoof2005_alpha_sq_uniform
        ice%friction_coef_2 = C%slid_Schoof2005_beta_sq_uniform

      case ('Tsai2015')
        ! Tsai (2015) sliding law

        ! Initialise with uniform friction coefficients
        ice%friction_coef_1 = C%slid_Tsai2015_alpha_sq_uniform
        ice%friction_coef_2 = C%slid_Tsai2015_beta_sq_uniform

      case ('Zoet-Iverson')
        ! Zoet-Iverson sliding law (Zoet & Iverson, 2020)

        ! Initialise with uniform friction coefficient
        ice%friction_coef_1 = C%slid_ZI_phi_fric_uniform

      case default
        ! Unknown case
        call crash('unknown choice_sliding_law "' // &
                    trim( C%choice_sliding_law) // '"!')

    end select

    ! === Friction coefficients choice ===
    ! ====================================

    select case (C%choice_bed_roughness)

      case ('uniform')
        ! Uniform values already assigned above

      case ('parameterised')
        ! Apply the chosen parameterisation for friction coefficients
        ! call calc_friction_coefficients( mesh, ice)
        call crash('parameterised friction coefficients not yet implemented!')

      case ('read_from_file')
        ! If friction coefficients are prescribed, read them from the provided file
        ! call read_friction_coefficients_from_file( mesh, ice)
        call crash('friction coefficients from file not yet implemented!')

      case default
        ! Unknown case
        call crash('unknown choice_bed_roughness "' // &
                    trim( C%choice_bed_roughness) // '"!')

    end select

    ! === Finalisation ===
    ! ====================

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_bed_roughness

! ===== Geothermal heat =====
! ===========================

  subroutine initialise_geothermal_heat_flux( mesh, ice)
    ! Initialise geothermal heat flux

    implicit none

    ! Input variables:
    type(type_mesh),      intent(in)    :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=256), parameter       :: routine_name = 'initialise_geothermal_heat_flux'

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    call init_routine( routine_name)

    ! === GHF source ===
    ! ==================

    select case (C%choice_geothermal_heat_flux)

      case ('constant')
        ! Uniform value over whole domain
        ice%geothermal_heat_flux = C%constant_geothermal_heat_flux

      case ('read_from_file')
        ! Spatially variable field
        call crash ('GHF read from file not yet implemented!')

      case default
        ! Unknown case
        call crash('unknown choice_geothermal_heat_flux "' // &
                    trim( C%choice_geothermal_heat_flux) // '"!')

    end select

    ! === Finalisation ===
    ! ====================

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_geothermal_heat_flux

END MODULE basal_conditions_main
