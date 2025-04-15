MODULE laddie_main

  ! The main laddie model module.

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE control_resources_and_error_messaging                  , ONLY: crash, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  USE laddie_model_types                                     , ONLY: type_laddie_model, type_laddie_timestep
  USE ocean_model_types                                      , ONLY: type_ocean_model
  USE BMB_model_types                                        , ONLY: type_BMB_model
  USE reallocate_mod                                         , ONLY: reallocate_bounds
  USE remapping_main                                         , ONLY: map_from_mesh_to_mesh_with_reallocation_2D
  USE laddie_utilities                                       , ONLY: compute_ambient_TS, allocate_laddie_model, &
                                                                     allocate_laddie_timestep, map_H_a_b, map_H_a_c, &
                                                                     print_diagnostics
  USE laddie_thickness                                       , ONLY: compute_H_npx
  USE laddie_velocity                                        , ONLY: compute_UV_npx, compute_viscUV
  USE laddie_tracers                                         , ONLY: compute_TS_npx, compute_diffTS
  use laddie_operators                                       , only: update_laddie_operators
  USE mesh_utilities                                         , ONLY: extrapolate_Gaussian
  USE mpi_distributed_memory                                 , ONLY: gather_to_all
  use mesh_integrate_over_domain, only: calc_and_print_min_mean_max

  IMPLICIT NONE

CONTAINS

! ===== Main routines =====
! =========================

  SUBROUTINE run_laddie_model( mesh, ice, ocean, laddie, time, duration)
    ! Run the laddie model

    ! In- and output variables

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_ocean_model),                 INTENT(IN)    :: ocean
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie
    REAL(dp),                               INTENT(IN)    :: time
    REAL(dp),                               INTENT(IN)    :: duration

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'run_laddie_model'
    INTEGER                                               :: vi, ti, ei
    REAL(dp)                                              :: tl               ! [s] Laddie time
    REAL(dp)                                              :: dt               ! [s] Laddie time step
    REAL(dp), PARAMETER                                   :: time_relax_laddie = 0.02_dp ! [days]
    REAL(dp), PARAMETER                                   :: fac_dt_relax = 3.0_dp ! Reduction factor of time step


    ! Add routine to path
    CALL init_routine( routine_name)

    ! == Preparation ==
    ! =================

    ! Extrapolate data into new cells
    CALL extrapolate_laddie_variables( mesh, ice, laddie)

    ! == Update masks ==
    CALL update_laddie_masks( mesh, ice, laddie)

    ! Set values to zero if outside laddie mask
    DO vi = mesh%vi1, mesh%vi2
      IF (.NOT. laddie%mask_a( vi)) THEN
        laddie%now%H( vi)     = 0.0_dp
        laddie%now%T( vi)     = 0.0_dp
        laddie%now%S( vi)     = 0.0_dp
        laddie%melt( vi)  = 0.0_dp
        laddie%entr( vi)  = 0.0_dp
      END IF
    END DO

    DO ti = mesh%ti1, mesh%ti2
      IF (.NOT. laddie%mask_b( ti)) THEN
        laddie%now%U( ti)     = 0.0_dp
        laddie%now%V( ti)     = 0.0_dp
        laddie%now%H_b( ti)   = 0.0_dp
      END IF
    END DO

    ! Simply set H_c zero everywhere, will be recomputed through mapping later
    laddie%now%H_c = 0.0_dp

    ! == Update operators ==
    CALL update_laddie_operators( mesh, ice, laddie)

    ! == Main time loop ==
    ! ====================

    tl = 0.0_dp

    DO WHILE (tl < duration * sec_per_day)

      ! Set time step
      IF (tl < time_relax_laddie * sec_per_day) THEN
        ! Relaxation, take short time step
        dt = C%dt_laddie / fac_dt_relax
      ELSE
        ! Regular timestep
        dt = C%dt_laddie
      END IF

      SELECT CASE(C%choice_laddie_integration_scheme)
        CASE DEFAULT
          CALL crash('unknown choice_laddie_integration_scheme "' // TRIM( C%choice_laddie_integration_scheme) // '"')
        CASE ('euler')
          CALL integrate_euler( mesh, ice, ocean, laddie, tl, time, dt)
        CASE ('fbrk3')
          CALL integrate_fbrk3( mesh, ice, ocean, laddie, tl, time, dt)
      END SELECT

      ! Display or save fields
      CALL print_diagnostics( mesh, laddie, tl)

      call crash('whoopsiedaisy')

    END DO !DO WHILE (tl < C%time_duration_laddie)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_laddie_model

  SUBROUTINE initialise_laddie_model( mesh, laddie, ocean, ice)
    ! Initialise the laddie model

    ! In- and output variables

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie
    TYPE(type_ocean_model),                 INTENT(IN)    :: ocean
    TYPE(type_ice_model),                   INTENT(IN)    :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'initialise_laddie_model'
    INTEGER                                               :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to terminal
    IF (par%primary)  WRITE(*,"(A)") '   Initialising LADDIE model...'

    ! Allocate variables
    CALL allocate_laddie_model( mesh, laddie)

    ! Mask on a grid
    DO vi = mesh%vi1, mesh%vi2
      laddie%mask_a( vi)  = ice%mask_floating_ice( vi)
    END DO

    ! == Update operators ==
    CALL update_laddie_operators( mesh, ice, laddie)

    ! Initialise requested timesteps
    CALL initialise_laddie_model_timestep( mesh, laddie, ocean, ice, laddie%now)

    SELECT CASE(C%choice_laddie_integration_scheme)
      CASE DEFAULT
        CALL crash('unknown choice_laddie_integration_scheme "' // TRIM( C%choice_laddie_integration_scheme) // '"')
      CASE ('euler')
        CALL initialise_laddie_model_timestep( mesh, laddie, ocean, ice, laddie%np1)
      CASE ('fbrk3')
        CALL initialise_laddie_model_timestep( mesh, laddie, ocean, ice, laddie%np13)
        CALL initialise_laddie_model_timestep( mesh, laddie, ocean, ice, laddie%np12)
        CALL initialise_laddie_model_timestep( mesh, laddie, ocean, ice, laddie%np1)
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_laddie_model

  SUBROUTINE initialise_laddie_model_timestep( mesh, laddie, ocean, ice, npx)
    ! Initialise the laddie model for given timestep

    ! In- and output variables

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie
    TYPE(type_ocean_model),                 INTENT(IN)    :: ocean
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_laddie_timestep),             INTENT(INOUT) :: npx

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'initialise_laddie_model_timestep'
    INTEGER                                               :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate timestep
    CALL allocate_laddie_timestep( mesh, npx)

    ! Layer thickness
    DO vi = mesh%vi1, mesh%vi2
       IF (laddie%mask_a( vi)) THEN
         npx%H( vi)      = C%laddie_initial_thickness
       END IF
    END DO

    ! Layer thickness on b and c grid
    CALL map_H_a_b( mesh, laddie, npx%H, npx%H_b)
    CALL map_H_a_c( mesh, laddie, npx%H, npx%H_c)

    ! Initialise ambient T and S
    CALL compute_ambient_TS( mesh, ice, ocean, laddie, npx%H)

    ! Initialise main T and S
    DO vi = mesh%vi1, mesh%vi2
       IF (laddie%mask_a( vi)) THEN
         npx%T( vi)      = laddie%T_amb( vi) + C%laddie_initial_T_offset
         npx%S( vi)      = laddie%S_amb( vi) + C%laddie_initial_S_offset
       END IF
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_laddie_model_timestep

  SUBROUTINE integrate_euler( mesh, ice, ocean, laddie, tl, time, dt)
    ! Integrate 1 timestep Euler scheme

    ! In- and output variables

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_ocean_model),                 INTENT(IN)    :: ocean
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie
    REAL(dp),                               INTENT(INOUT) :: tl
    REAL(dp),                               INTENT(IN)    :: time
    REAL(dp),                               INTENT(IN)    :: dt

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'integrate_euler'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Integrate H 1 time step
    CALL compute_H_npx( mesh, ice, ocean, laddie, laddie%now, laddie%np1, time, dt)

    ! Update diffusive terms based on now time step
    CALL update_diffusive_terms( mesh, ice, laddie, laddie%now)

    ! Integrate U and V 1 time step
    CALL compute_UV_npx( mesh, ice, ocean, laddie, laddie%now, laddie%np1, laddie%now%H, dt, .true.)

    ! Integrate T and S 1 time step
    CALL compute_TS_npx( mesh, ice, laddie, laddie%now, laddie%np1, laddie%now%H, dt, .true.)

    ! == Move time ==
    CALL move_laddie_timestep( laddie, tl, dt)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE integrate_euler

  SUBROUTINE integrate_fbrk3( mesh, ice, ocean, laddie, tl, time, dt)
    ! Integrate 1 timestep Forward-Backward Runge Kutta 3 scheme

    ! Based on Lilly et al (2023, MWR) doi:10.1175/MWR-D-23-0113.1

    ! In- and output variables

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_ocean_model),                 INTENT(IN)    :: ocean
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie
    REAL(dp),                               INTENT(INOUT) :: tl
    REAL(dp),                               INTENT(IN)    :: time
    REAL(dp),                               INTENT(IN)    :: dt

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'integrate_fbrk3'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! == Stage 1: explicit 1/3 timestep ==
    ! == RHS terms defined at n ==========
    ! ====================================

    ! Integrate H 1/3 time step
    CALL compute_H_npx( mesh, ice, ocean, laddie, laddie%now, laddie%np13, time, dt/3)

    ! Compute Hstar
    laddie%Hstar = C%laddie_fbrk3_beta1 * laddie%np13%H + (1-C%laddie_fbrk3_beta1) * laddie%now%H
    call calc_and_print_min_mean_max( mesh, laddie%Hstar, 'laddie%Hstar')

    ! Update diffusive terms
    CALL update_diffusive_terms( mesh, ice, laddie, laddie%now)

    ! Integrate U and V 1/3 time step
    CALL compute_UV_npx( mesh, ice, ocean, laddie, laddie%now, laddie%np13, laddie%Hstar, dt/3, .false.)

    ! Integrate T and S 1/3 time step
    CALL compute_TS_npx( mesh, ice, laddie, laddie%now, laddie%np13, laddie%now%H, dt/3, .false.)

    ! == Stage 2: explicit 1/2 timestep ==
    ! == RHS terms defined at n + 1/3 ====
    ! ====================================

    ! Integrate H 1/2 time step
    CALL compute_H_npx( mesh, ice, ocean, laddie, laddie%np13, laddie%np12, time, dt/2)

    ! Compute new Hstar
    laddie%Hstar = C%laddie_fbrk3_beta2 * laddie%np12%H + (1-C%laddie_fbrk3_beta2) * laddie%now%H
    call calc_and_print_min_mean_max( mesh, laddie%Hstar, 'laddie%Hstar')

    ! Update diffusive terms
    !CALL update_diffusive_terms( mesh, ice, laddie, laddie%np13)

    ! Integrate U and V 1/2 time step
    CALL compute_UV_npx( mesh, ice, ocean, laddie, laddie%np13, laddie%np12, laddie%Hstar, dt/2, .false.)

    ! Integrate T and S 1/2 time step
    CALL compute_TS_npx( mesh, ice, laddie, laddie%np13, laddie%np12, laddie%np13%H, dt/2, .false.)

    ! == Stage 3: explicit 1 timestep ====
    ! == RHS terms defined at n + 1/2 ====
    ! ====================================

    ! Integrate H 1 time step
    CALL compute_H_npx( mesh, ice, ocean, laddie, laddie%np12, laddie%np1, time, dt)

    ! Compute new Hstar
    laddie%Hstar = C%laddie_fbrk3_beta3 * laddie%np1%H + (1-2*C%laddie_fbrk3_beta3) * laddie%np12%H + C%laddie_fbrk3_beta3 * laddie%now%H
    call calc_and_print_min_mean_max( mesh, laddie%Hstar, 'laddie%Hstar')

    ! Update diffusive terms
    !CALL update_diffusive_terms( mesh, ice, laddie, laddie%np12)

    ! Integrate U and V 1 time step
    CALL compute_UV_npx( mesh, ice, ocean, laddie, laddie%np12, laddie%np1, laddie%Hstar, dt, .true.)

    ! Integrate T and S 1 time step
    CALL compute_TS_npx( mesh, ice, laddie, laddie%np12, laddie%np1, laddie%np12%H, dt, .true.)

    ! ===============
    ! == Move time ==
    CALL move_laddie_timestep( laddie, tl, dt)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE integrate_fbrk3

  SUBROUTINE move_laddie_timestep( laddie, tl, dt)
    ! Increase laddie time tl by timestep dt and overwrite now timestep

    ! In- and output variables

    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie
    REAL(dp),                               INTENT(INOUT) :: tl
    REAL(dp),                               INTENT(IN)    :: dt

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'move_laddie_timestep'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Increase laddie time
    tl = tl + dt

    ! Move main variables by 1 time step
    laddie%now%H = laddie%np1%H
    laddie%now%T = laddie%np1%T
    laddie%now%S = laddie%np1%S
    laddie%now%U = laddie%np1%U
    laddie%now%V = laddie%np1%V
    laddie%now%H_b = laddie%np1%H_b
    laddie%now%H_c = laddie%np1%H_c
    laddie%now%U_a = laddie%np1%U_a
    laddie%now%U_c = laddie%np1%U_c
    laddie%now%V_a = laddie%np1%V_a
    laddie%now%V_c = laddie%np1%V_c

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE move_laddie_timestep

  SUBROUTINE update_diffusive_terms( mesh, ice, laddie, npxref)
    ! Update diffusivity and viscosity. Based on reference timestep npxref

    ! For stability, most studies base diffusive terms on the now timestep

    ! In- and output variables

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie
    TYPE(type_laddie_timestep),             INTENT(IN)    :: npxref

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'update_diffusive_terms'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Compute diffusivities
    CALL compute_diffTS( mesh, ice, laddie, npxref)

    ! Compute viscosities
    CALL compute_viscUV( mesh, ice, laddie, npxref)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE update_diffusive_terms

  SUBROUTINE update_laddie_masks( mesh, ice, laddie)
    ! Update bunch of masks for laddie at the start of a new run

    ! In- and output variables

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'update_laddie_masks'
    INTEGER                                               :: vi, ti, i, no
    LOGICAL, DIMENSION(mesh%nV)                           :: mask_a_tot, mask_gr_a_tot, mask_oc_a_tot

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Mask on a grid
    DO vi = mesh%vi1, mesh%vi2
      ! Check whether vertex on border
      IF (mesh%VBI( vi) > 0) THEN
        laddie%mask_a( vi)    = .false.
        laddie%mask_gr_a( vi) = .true.
      ELSE IF (ice%Hib( vi) - ice%Hb( vi) < 2*C%laddie_thickness_minimum) THEN
        laddie%mask_a( vi)    = .false.
        laddie%mask_gr_a( vi) = .true.
      ELSE IF (ice%Hi( vi) < 1.0 .and. ice%mask_floating_ice( vi)) THEN
        laddie%mask_a( vi)    = .false.
        laddie%mask_oc_a( vi) = .true.
      ELSE
        ! Inherit regular masks
        laddie%mask_a( vi)    = ice%mask_floating_ice( vi)
        laddie%mask_gr_a( vi) = ice%mask_grounded_ice( vi) .OR. ice%mask_icefree_land( vi)
        laddie%mask_oc_a( vi) = ice%mask_icefree_ocean( vi)
      END IF
    END DO

    ! Mask on b grid
    CALL gather_to_all( laddie%mask_a, mask_a_tot)
    CALL gather_to_all( laddie%mask_gr_a, mask_gr_a_tot)
    CALL gather_to_all( laddie%mask_oc_a, mask_oc_a_tot)

    DO ti = mesh%ti1, mesh%ti2
      ! Initialise as false to overwrite previous mask
      laddie%mask_b( ti)    = .false.
      laddie%mask_gl_b( ti) = .false.
      laddie%mask_cf_b( ti) = .false.
      laddie%mask_oc_b( ti) = .false.

      ! Define floating mask if any of the three vertices is floating
      DO i = 1, 3
        vi = mesh%Tri( ti, i)
        IF (mask_a_tot( vi)) THEN
          ! Set true if any of the three vertices is floating
          laddie%mask_b( ti) = .true.
        END IF
      END DO

      ! Define grounding line triangles
      DO i = 1, 3
        vi = mesh%Tri( ti, i)
        ! Check if any connected vertex is grounded
        IF (mask_gr_a_tot( vi)) THEN
          ! Omit triangle from floating mask. Adjust for no slip conditions
          laddie%mask_b( ti) = .false.
          ! Define as grounding line triangle
          laddie%mask_gl_b( ti) = .true.
        END IF
      END DO

      ! Also define border triangles as grounding line
      IF (mesh%TriBI( ti) > 0) THEN
        ! Omit triangle from floating mask. Adjust for no slip conditions
        laddie%mask_b( ti) = .false.
        ! Define as grounding line triangle
        laddie%mask_gl_b( ti) = .true.
      END IF

      ! For non-grounding line triangles:
      IF (.NOT. laddie%mask_gl_b( ti)) THEN
        ! Define calving front triangles
        DO i = 1, 3
          vi = mesh%Tri( ti, i)
          ! Check if any vertex is icefree ocean
          IF (mask_oc_a_tot( vi)) THEN
            ! Define as calving front triangle
            laddie%mask_cf_b( ti) = .true.
          END IF
        END DO
      END IF

      ! Define ocean triangles
      no = 0 ! Number of ice free ocean vertices
      DO i = 1, 3
        vi = mesh%Tri( ti, i)
        ! Check if vertex is icefree ocean
        IF (mask_oc_a_tot( vi)) THEN
          no = no + 1
        END IF
      END DO
      ! Check whether all vertices are icefree ocean
      IF (no == 3) THEN
        ! Define as ocean triangle
        laddie%mask_oc_b( ti) = .true.
      END IF

    END DO !ti = mesh%ti1, mesh%ti2

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE update_laddie_masks

  SUBROUTINE extrapolate_laddie_variables( mesh, ice, laddie)
    ! Update bunch of masks for laddie at the start of a new run

    ! In- and output variables

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'extrapolate_laddie_variables'
    INTEGER                                               :: vi
    INTEGER, DIMENSION(mesh%vi1: mesh%vi2)                :: mask
    REAL(dp), PARAMETER                                   :: sigma = 16000.0_dp

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Initialise mask
    mask = 0

    ! Determine mask for seed (2: previously floating cells), fill (1: new floating cells), or ignore (0: grounded/ocean)
    DO vi = mesh%vi1, mesh%vi2
      ! Skip if vertex is at border
      IF (mesh%VBI( vi) > 0) CYCLE

      ! Skip if water column thickness is insufficient, treated as grounded for now
      IF (ice%Hib( vi) - ice%Hb( vi) < 2*C%laddie_thickness_minimum) CYCLE

      IF (ice%Hi( vi) < 1.0 .and. ice%mask_floating_ice( vi)) CYCLE

      ! Currently floating ice, so either seed or fill here
      IF (ice%mask_floating_ice( vi)) THEN
        IF (laddie%mask_a( vi)) THEN
          ! Data already available here, so use as seed
          mask( vi) = 2
        ELSE
          ! New floating cells, so fill here
          mask (vi) = 1
        END IF
      END IF
    END DO

    ! Apply extrapolation to H, T and S
    CALL extrapolate_Gaussian( mesh, mask, laddie%now%H, sigma)
    CALL extrapolate_Gaussian( mesh, mask, laddie%now%T, sigma)
    CALL extrapolate_Gaussian( mesh, mask, laddie%now%S, sigma)

    ! The above should ensure that all (newly) floating vertices have a non-zero thickness
    ! In case the extrapolation did not cover this, apply a backup check to set values
    ! at non-zero initialisation
    DO vi = mesh%vi1, mesh%vi2
      ! Skip if vertex is at border
      IF (mesh%VBI( vi) > 0) CYCLE

      ! Skip if water column thickness is insufficient, treated as grounded for now
      IF (ice%Hib( vi) - ice%Hb( vi) < 2*C%laddie_thickness_minimum) CYCLE

      IF (ice%Hi( vi) < 1.0 .and. ice%mask_floating_ice( vi)) CYCLE

      ! Currently floating ice, so either seed or fill here
      IF (ice%mask_floating_ice( vi)) THEN
        IF (laddie%now%H( vi) == 0.0_dp) THEN
          laddie%now%H( vi) = C%laddie_thickness_minimum
          laddie%now%T( vi) = laddie%T_amb( vi) + C%laddie_initial_T_offset
          laddie%now%S( vi) = laddie%S_amb( vi) + C%laddie_initial_S_offset
        END IF
      END IF
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE extrapolate_laddie_variables

  SUBROUTINE remap_laddie_model( mesh_old, mesh_new, ice, ocean, laddie, time)
    ! Reallocate and remap laddie variables

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                        INTENT(IN)    :: mesh_new
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_ocean_model),                 INTENT(IN)    :: ocean
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'remap_laddie_model'
    INTEGER                                               :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    ! ! == Regular variables ==

    ! ! Thickness
    ! CALL reallocate_bounds( laddie%dH_dt,                mesh_new%vi1, mesh_new%vi2)

    ! ! Temperatures
    ! CALL reallocate_bounds( laddie%T_amb,                mesh_new%vi1, mesh_new%vi2)
    ! CALL reallocate_bounds( laddie%T_base,               mesh_new%vi1, mesh_new%vi2)
    ! CALL reallocate_bounds( laddie%T_freeze,             mesh_new%vi1, mesh_new%vi2)

    ! ! Salinities
    ! CALL reallocate_bounds( laddie%S_amb,                mesh_new%vi1, mesh_new%vi2)
    ! CALL reallocate_bounds( laddie%S_base,               mesh_new%vi1, mesh_new%vi2)

    ! ! Densities and buoyancies
    ! CALL reallocate_bounds( laddie%rho,                  mesh_new%vi1, mesh_new%vi2)
    ! CALL reallocate_bounds( laddie%rho_amb,              mesh_new%vi1, mesh_new%vi2)
    ! CALL reallocate_bounds( laddie%drho_amb,             mesh_new%vi1, mesh_new%vi2)
    ! CALL reallocate_bounds( laddie%Hdrho_amb,            mesh_new%vi1, mesh_new%vi2)
    ! CALL reallocate_bounds( laddie%Hdrho_amb_b,          mesh_new%ti1, mesh_new%ti2)
    ! CALL reallocate_bounds( laddie%drho_base,            mesh_new%vi1, mesh_new%vi2)

    ! ! Friction velocity
    ! CALL reallocate_bounds( laddie%u_star,               mesh_new%vi1, mesh_new%vi2)

    ! ! Physical parameter fields
    ! CALL reallocate_bounds( laddie%gamma_T,              mesh_new%vi1, mesh_new%vi2)
    ! CALL reallocate_bounds( laddie%gamma_S,              mesh_new%vi1, mesh_new%vi2)
    ! CALL reallocate_bounds( laddie%A_h,                  mesh_new%ti1, mesh_new%ti2)
    ! CALL reallocate_bounds( laddie%K_h,                  mesh_new%vi1, mesh_new%vi2)

    ! ! Vertical rates
    ! CALL reallocate_bounds( laddie%melt,                 mesh_new%vi1, mesh_new%vi2)
    ! CALL reallocate_bounds( laddie%entr,                 mesh_new%vi1, mesh_new%vi2)
    ! CALL reallocate_bounds( laddie%entr_dmin,            mesh_new%vi1, mesh_new%vi2)
    ! CALL reallocate_bounds( laddie%detr,                 mesh_new%vi1, mesh_new%vi2)
    ! CALL reallocate_bounds( laddie%entr_tot,             mesh_new%vi1, mesh_new%vi2)

    ! ! Horizontal fluxes
    ! CALL reallocate_bounds( laddie%divQH,                mesh_new%vi1, mesh_new%vi2)
    ! CALL reallocate_bounds( laddie%divQU,                mesh_new%ti1, mesh_new%ti2)
    ! CALL reallocate_bounds( laddie%divQV,                mesh_new%ti1, mesh_new%ti2)
    ! CALL reallocate_bounds( laddie%divQT,                mesh_new%vi1, mesh_new%vi2)
    ! CALL reallocate_bounds( laddie%divQS,                mesh_new%vi1, mesh_new%vi2)

    ! ! Viscosities
    ! CALL reallocate_bounds( laddie%viscU,                mesh_new%ti1, mesh_new%ti2)
    ! CALL reallocate_bounds( laddie%viscV,                mesh_new%ti1, mesh_new%ti2)

    ! ! Diffusivities
    ! CALL reallocate_bounds( laddie%diffT,                mesh_new%vi1, mesh_new%vi2)
    ! CALL reallocate_bounds( laddie%diffS,                mesh_new%vi1, mesh_new%vi2)

    ! ! RHS terms
    ! CALL reallocate_bounds( laddie%ddrho_amb_dx_b,       mesh_new%ti1, mesh_new%ti2)
    ! CALL reallocate_bounds( laddie%ddrho_amb_dy_b,       mesh_new%ti1, mesh_new%ti2)
    ! CALL reallocate_bounds( laddie%dH_dx_b,              mesh_new%ti1, mesh_new%ti2)
    ! CALL reallocate_bounds( laddie%dH_dy_b,              mesh_new%ti1, mesh_new%ti2)
    ! CALL reallocate_bounds( laddie%detr_b,               mesh_new%ti1, mesh_new%ti2)

    ! ! Masks
    ! CALL reallocate_bounds( laddie%mask_a,               mesh_new%vi1, mesh_new%vi2)
    ! CALL reallocate_bounds( laddie%mask_gr_a,            mesh_new%vi1, mesh_new%vi2)
    ! CALL reallocate_bounds( laddie%mask_oc_a,            mesh_new%vi1, mesh_new%vi2)
    ! CALL reallocate_bounds( laddie%mask_b,               mesh_new%ti1, mesh_new%ti2)
    ! CALL reallocate_bounds( laddie%mask_gl_b,            mesh_new%ti1, mesh_new%ti2)
    ! CALL reallocate_bounds( laddie%mask_cf_b,            mesh_new%ti1, mesh_new%ti2)
    ! CALL reallocate_bounds( laddie%mask_oc_b,            mesh_new%ti1, mesh_new%ti2)

    ! ! == Re-initialise masks ==
    ! CALL update_laddie_masks( mesh_new, ice, laddie)

    ! ! == Update operators ==
    ! CALL update_laddie_operators( mesh_new, ice, laddie)

    ! ! == Timestep variables ==
    ! CALL remap_laddie_timestep( mesh_old, mesh_new, ice, ocean, laddie, laddie%now)

    ! SELECT CASE(C%choice_laddie_integration_scheme)
    !   CASE DEFAULT
    !     CALL crash('unknown choice_laddie_integration_scheme "' // TRIM( C%choice_laddie_integration_scheme) // '"')
    !   CASE ('euler')
    !     CALL remap_laddie_timestep( mesh_old, mesh_new, ice, ocean, laddie, laddie%np1)
    !   CASE ('fbrk3')
    !     CALL remap_laddie_timestep( mesh_old, mesh_new, ice, ocean, laddie, laddie%np13)
    !     CALL remap_laddie_timestep( mesh_old, mesh_new, ice, ocean, laddie, laddie%np12)
    !     CALL remap_laddie_timestep( mesh_old, mesh_new, ice, ocean, laddie, laddie%np1)
    ! END SELECT

    ! ! == Re-initialise ==
    ! CALL run_laddie_model( mesh_new, ice, ocean, laddie, time, C%time_duration_laddie_init)

    call crash('fixme!')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_laddie_model

  SUBROUTINE remap_laddie_timestep( mesh_old, mesh_new, ice, ocean, laddie, npx)
    ! Remap laddie timestep

    ! In- and output variables

    TYPE(type_mesh),                        INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                        INTENT(IN)    :: mesh_new
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_ocean_model),                 INTENT(IN)    :: ocean
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie
    TYPE(type_laddie_timestep),             INTENT(INOUT) :: npx

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'remap_laddie_timestep'
    INTEGER                                               :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    ! ! Reallocate
    ! CALL map_from_mesh_to_mesh_with_reallocation_2D( mesh_old, mesh_new, npx%H, '2nd_order_conservative')
    ! CALL reallocate_bounds( npx%H_b,                     mesh_new%ti1, mesh_new%ti2)
    ! CALL reallocate_bounds( npx%H_c,                     mesh_new%ei1, mesh_new%ei2)
    ! CALL reallocate_bounds( npx%U,                       mesh_new%ti1, mesh_new%ti2)
    ! CALL reallocate_bounds( npx%U_a,                     mesh_new%vi1, mesh_new%vi2)
    ! CALL reallocate_bounds( npx%U_c,                     mesh_new%ei1, mesh_new%ei2)
    ! CALL reallocate_bounds( npx%V,                       mesh_new%ti1, mesh_new%ti2)
    ! CALL reallocate_bounds( npx%V_a,                     mesh_new%vi1, mesh_new%vi2)
    ! CALL reallocate_bounds( npx%V_c,                     mesh_new%ei1, mesh_new%ei2)
    ! CALL map_from_mesh_to_mesh_with_reallocation_2D( mesh_old, mesh_new, npx%T, '2nd_order_conservative')
    ! CALL map_from_mesh_to_mesh_with_reallocation_2D( mesh_old, mesh_new, npx%S, '2nd_order_conservative')

    ! ! == Re-initialise ==

    ! ! Layer thickness on b and c grid
    ! CALL map_H_a_b( mesh_new, laddie, npx%H, npx%H_b)
    ! CALL map_H_a_c( mesh_new, laddie, npx%H, npx%H_c)

    call crash('fixme!')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_laddie_timestep

END MODULE laddie_main

