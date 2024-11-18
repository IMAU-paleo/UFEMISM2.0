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
  USE laddie_utilities                                       , ONLY: compute_ambient_TS, allocate_laddie_model, &
                                                                     allocate_laddie_timestep, map_H_a_b, map_H_a_c, &
                                                                     print_diagnostics
  USE laddie_thickness                                       , ONLY: compute_H_npx
  USE laddie_velocity                                        , ONLY: compute_UV_npx, compute_viscUV
  USE laddie_tracers                                         , ONLY: compute_TS_npx, compute_diffTS
  USE mpi_distributed_memory                                 , ONLY: gather_to_all_logical_1D

  IMPLICIT NONE
    
CONTAINS
    
! ===== Main routines =====
! =========================

  SUBROUTINE run_laddie_model( mesh, ice, ocean, laddie)
    ! Run the laddie model

    ! In- and output variables

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_ocean_model),                 INTENT(IN)    :: ocean
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'run_laddie_model'
    INTEGER                                               :: vi, ti
    REAL(dp)                                              :: tl               ! [s] Laddie time
    REAL(dp)                                              :: dt               ! [s] Laddie time step
 
    ! Add routine to path
    CALL init_routine( routine_name)

    ! == Preparation ==
    ! =================

    ! Get time step
    tl = 0.0_dp
    dt = C%dt_laddie

    ! == Update masks ==
    CALL update_laddie_masks( mesh, ice, laddie)

    ! Extrapolate new cells
    ! TODO, use Gaussian extrap routine

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

    laddie%H_c = 0.0_dp

    ! == Main time loop ==
    ! ====================

    DO WHILE (tl < C%time_duration_laddie * sec_per_day)

      SELECT CASE(C%choice_laddie_integration_scheme)
        CASE DEFAULT
          CALL crash('unknown choice_laddie_integration_scheme "' // TRIM( C%choice_laddie_integration_scheme) // '"')
        CASE ('euler')
          CALL integrate_euler( mesh, ice, ocean, laddie, tl, dt)  
        CASE ('fbrk3')
          CALL integrate_fbrk3( mesh, ice, ocean, laddie, tl, dt)  
      END SELECT

      ! Display or save fields
      CALL print_diagnostics( laddie, tl)

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
    IF (par%master)  WRITE(*,"(A)") '   Initialising LADDIE model...'

    ! Allocate variables
    CALL allocate_laddie_model( mesh, laddie)
    CALL allocate_laddie_timestep( mesh, laddie%now)

    ! Allocate timestep
    SELECT CASE(C%choice_laddie_integration_scheme)
      CASE DEFAULT
        CALL crash('unknown choice_laddie_integration_scheme "' // TRIM( C%choice_laddie_integration_scheme) // '"')
      CASE ('euler')
        CALL allocate_laddie_timestep( mesh, laddie%np1)
      CASE ('fbrk3')
        CALL allocate_laddie_timestep( mesh, laddie%np13)
        CALL allocate_laddie_timestep( mesh, laddie%np12)
        CALL allocate_laddie_timestep( mesh, laddie%np1)
    END SELECT

    ! Mask on a grid
    DO vi = mesh%vi1, mesh%vi2
      laddie%mask_a( vi)  = ice%mask_floating_ice( vi)
    END DO

    ! Layer thickness 
    DO vi = mesh%vi1, mesh%vi2
       IF (laddie%mask_a( vi)) THEN
         laddie%now%H( vi)      = C%laddie_initial_thickness
         SELECT CASE(C%choice_laddie_integration_scheme)
           CASE DEFAULT
             CALL crash('unknown choice_laddie_integration_scheme "' // TRIM( C%choice_laddie_integration_scheme) // '"')
           CASE ('euler')
             laddie%np1%H( vi)   = C%laddie_initial_thickness
           CASE ('fbrk3')
             laddie%np13%H( vi)  = C%laddie_initial_thickness
             laddie%np12%H( vi)  = C%laddie_initial_thickness
             laddie%np1%H( vi)   = C%laddie_initial_thickness
         END SELECT
       END IF
    END DO

    ! Layer thickness on b grid
    CALL map_H_a_b( mesh, laddie, laddie%now%H, laddie%now%H_b)
    SELECT CASE(C%choice_laddie_integration_scheme)
      CASE DEFAULT
        CALL crash('unknown choice_laddie_integration_scheme "' // TRIM( C%choice_laddie_integration_scheme) // '"')
      CASE ('euler')
        CALL map_H_a_b( mesh, laddie, laddie%np1%H, laddie%np1%H_b)
      CASE ('fbrk3')
        CALL map_H_a_b( mesh, laddie, laddie%np13%H, laddie%np13%H_b)
        CALL map_H_a_b( mesh, laddie, laddie%np12%H, laddie%np12%H_b)
        CALL map_H_a_b( mesh, laddie, laddie%np1%H, laddie%np1%H_b)
    END SELECT

    ! Layer thickness on c grid
    CALL map_H_a_c( mesh, laddie, laddie%now%H, laddie%now%H_c)
    SELECT CASE(C%choice_laddie_integration_scheme)
      CASE DEFAULT
        CALL crash('unknown choice_laddie_integration_scheme "' // TRIM( C%choice_laddie_integration_scheme) // '"')
      CASE ('euler')
        CALL map_H_a_c( mesh, laddie, laddie%np1%H, laddie%np1%H_c)
      CASE ('fbrk3')
        CALL map_H_a_c( mesh, laddie, laddie%np13%H, laddie%np13%H_c)
        CALL map_H_a_c( mesh, laddie, laddie%np12%H, laddie%np12%H_c)
        CALL map_H_a_c( mesh, laddie, laddie%np1%H, laddie%np1%H_c)
    END SELECT

    ! Initialise ambient T and S
    CALL compute_ambient_TS( mesh, ice, ocean, laddie, laddie%now%H)

    ! Initialise main T and S
    DO vi = mesh%vi1, mesh%vi2
       IF (laddie%mask_a( vi)) THEN
         laddie%now%T( vi)      = laddie%T_amb( vi) + C%laddie_initial_T_offset 
         laddie%now%S( vi)      = laddie%S_amb( vi) + C%laddie_initial_S_offset
         SELECT CASE(C%choice_laddie_integration_scheme)
           CASE DEFAULT
             CALL crash('unknown choice_laddie_integration_scheme "' // TRIM( C%choice_laddie_integration_scheme) // '"')
           CASE ('euler')
             laddie%np1%T( vi)   = laddie%T_amb( vi) + C%laddie_initial_T_offset
             laddie%np1%S( vi)   = laddie%S_amb( vi) + C%laddie_initial_S_offset
           CASE ('fbrk3')
             laddie%np13%T( vi)  = laddie%T_amb( vi) + C%laddie_initial_T_offset
             laddie%np13%S( vi)  = laddie%S_amb( vi) + C%laddie_initial_S_offset
             laddie%np12%T( vi)  = laddie%T_amb( vi) + C%laddie_initial_T_offset
             laddie%np12%S( vi)  = laddie%S_amb( vi) + C%laddie_initial_S_offset
             laddie%np1%T( vi)   = laddie%T_amb( vi) + C%laddie_initial_T_offset
             laddie%np1%S( vi)   = laddie%S_amb( vi) + C%laddie_initial_S_offset
         END SELECT
       END IF
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_laddie_model

  SUBROUTINE integrate_euler( mesh, ice, ocean, laddie, tl, dt)
    ! Integrate 1 timestep Euler scheme 

    ! In- and output variables

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_ocean_model),                 INTENT(IN)    :: ocean
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie
    REAL(dp),                               INTENT(INOUT) :: tl
    REAL(dp),                               INTENT(IN)    :: dt

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'integrate_euler'
 
    ! Add routine to path
    CALL init_routine( routine_name)

    ! Integrate H 1 time step
    CALL compute_H_npx( mesh, ice, ocean, laddie, laddie%now, laddie%np1, dt)

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

  SUBROUTINE integrate_fbrk3( mesh, ice, ocean, laddie, tl, dt)
    ! Integrate 1 timestep Forward-Backward Runge Kutta 3 scheme 

    ! Based on Lilly et al (2023, MWR) doi:10.1175/MWR-D-23-0113.1

    ! In- and output variables

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_ocean_model),                 INTENT(IN)    :: ocean
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie
    REAL(dp),                               INTENT(INOUT) :: tl
    REAL(dp),                               INTENT(IN)    :: dt

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'integrate_fbrk3'
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2)                :: Hstar
 
    ! Add routine to path
    CALL init_routine( routine_name)

    ! == Stage 1: explicit 1/3 timestep ==
    ! == RHS terms defined at n ==========
    ! ====================================
 
    ! Integrate H 1/3 time step
    CALL compute_H_npx( mesh, ice, ocean, laddie, laddie%now, laddie%np13, dt/3)

    ! Compute Hstar
    Hstar = C%laddie_fbrk3_beta1 * laddie%np13%H + (1-C%laddie_fbrk3_beta1) * laddie%now%H

    ! Update diffusive terms
    CALL update_diffusive_terms( mesh, ice, laddie, laddie%now)

    ! Integrate U and V 1/3 time step
    CALL compute_UV_npx( mesh, ice, ocean, laddie, laddie%now, laddie%np13, Hstar, dt/3, .false.)

    ! Integrate T and S 1/3 time step
    CALL compute_TS_npx( mesh, ice, laddie, laddie%now, laddie%np13, laddie%now%H, dt/3, .false.)

    ! == Stage 2: explicit 1/2 timestep ==
    ! == RHS terms defined at n + 1/3 ====
    ! ====================================

    ! Integrate H 1/2 time step
    CALL compute_H_npx( mesh, ice, ocean, laddie, laddie%np13, laddie%np12, dt/2)

    ! Compute new Hstar
    Hstar = C%laddie_fbrk3_beta2 * laddie%np12%H + (1-C%laddie_fbrk3_beta2) * laddie%now%H

    ! Update diffusive terms
    !CALL update_diffusive_terms( mesh, ice, laddie, laddie%np13)

    ! Integrate U and V 1/2 time step
    CALL compute_UV_npx( mesh, ice, ocean, laddie, laddie%np13, laddie%np12, Hstar, dt/2, .false.)

    ! Integrate T and S 1/2 time step
    CALL compute_TS_npx( mesh, ice, laddie, laddie%np13, laddie%np12, laddie%np13%H, dt/2, .false.)

    ! == Stage 3: explicit 1 timestep ====
    ! == RHS terms defined at n + 1/2 ====
    ! ====================================

    ! Integrate H 1 time step
    CALL compute_H_npx( mesh, ice, ocean, laddie, laddie%np12, laddie%np1, dt)

    ! Compute new Hstar
    Hstar = C%laddie_fbrk3_beta3 * laddie%np1%H + (1-2*C%laddie_fbrk3_beta3) * laddie%np12%H + C%laddie_fbrk3_beta3 * laddie%now%H

    ! Update diffusive terms
    !CALL update_diffusive_terms( mesh, ice, laddie, laddie%np12)

    ! Integrate U and V 1 time step
    CALL compute_UV_npx( mesh, ice, ocean, laddie, laddie%np12, laddie%np1, Hstar, dt, .true.)

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
      ELSE
        ! Inherit regular masks
        laddie%mask_a( vi)    = ice%mask_floating_ice( vi)
        laddie%mask_gr_a( vi) = ice%mask_grounded_ice( vi) .OR. ice%mask_icefree_land( vi)
        laddie%mask_oc_a( vi) = ice%mask_icefree_ocean( vi)
      END IF
    END DO
      
    ! Mask on b grid
    CALL gather_to_all_logical_1D( laddie%mask_a, mask_a_tot)
    CALL gather_to_all_logical_1D( laddie%mask_gr_a, mask_gr_a_tot)
    CALL gather_to_all_logical_1D( laddie%mask_oc_a, mask_oc_a_tot)

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

END MODULE laddie_main

