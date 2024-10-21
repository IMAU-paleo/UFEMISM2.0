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
                                                                     calc_laddie_flux_divergence_matrix_upwind, &
                                                                     map_laddie_velocities_from_b_to_c_2D
  USE laddie_physics                                         , ONLY: compute_melt_rate, compute_entrainment, &
                                                                     compute_freezing_temperature, compute_buoyancy
  USE laddie_thickness                                       , ONLY: compute_H_npx, compute_divQH
  USE laddie_velocity                                        , ONLY: compute_UV_npx, compute_viscUV, compute_divQUV
  USE laddie_tracers                                         , ONLY: compute_TS_npx, compute_diffTS, compute_divQTS
  USE mesh_operators                                         , ONLY: ddx_a_b_2D, ddy_a_b_2D, map_a_b_2D, map_a_c_2D, map_b_a_2D
  USE petsc_basic                                            , ONLY: multiply_CSR_matrix_with_vector_1D
  USE CSR_sparse_matrix_utilities                            , ONLY: type_sparse_matrix_CSR_dp
  USE mpi_distributed_memory                                 , ONLY: gather_to_all_logical_1D

  IMPLICIT NONE
    
CONTAINS
    
! ===== Main routines =====
! =========================

  SUBROUTINE run_laddie_model( mesh, ice, ocean, laddie, time)
    ! Run the laddie model

    ! In- and output variables

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_ocean_model),                 INTENT(IN)    :: ocean
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'run_laddie_model'
    INTEGER                                               :: vi, ti, nf, i, no
    REAL(dp)                                              :: tl               ! [s] Laddie time
    REAL(dp)                                              :: dt               ! [s] Laddie time step
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2)                :: Hstar            ! [m] Reference thickness in integration
    LOGICAL, DIMENSION(mesh%nV)                           :: mask_a_tot
    LOGICAL, DIMENSION(mesh%nV)                           :: mask_gr_a_tot
    LOGICAL, DIMENSION(mesh%nV)                           :: mask_oc_a_tot
 
    ! Add routine to path
    CALL init_routine( routine_name)

    ! == Preparation ==
    ! =================

    ! Get time step
    tl = 0.0_dp
    dt = C%dt_laddie

    ! == Update masks ==
    ! Mask on a grid
    DO vi = mesh%vi1, mesh%vi2
      ! Check whether vertex on border
      IF (mesh%VBI( vi) > 0) THEN
        laddie%mask_a( vi)    = .false.
        laddie%mask_gr_a( vi) = .true.
      ELSE
        ! Inherit regular masks
        laddie%mask_a( vi)    = ice%mask_floating_ice( vi)
        laddie%mask_gr_a( vi) = ice%mask_grounded_ice( vi)
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

    ! Extrapolate new cells
    ! TODO

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

    ! Update ice shelf draft gradients
    CALL ddx_a_b_2D( mesh, ice%Hib , laddie%dHib_dx_b)
    CALL ddy_a_b_2D( mesh, ice%Hib , laddie%dHib_dy_b)

    ! Update secondary fields
    CALL update_secondary_fields( mesh, ice, ocean, laddie, laddie%now, laddie%now%H, laddie%now%H_b, laddie%now%H_c)

    ! == Main time loop ==
    ! ====================

    DO WHILE (tl <= C%time_duration_laddie * sec_per_day)

      SELECT CASE(C%choice_laddie_integration_scheme)
        CASE DEFAULT
          CALL crash('unknown choice_laddie_integration_scheme "' // TRIM( C%choice_laddie_integration_scheme) // '"')
        CASE ('euler')
          CALL integrate_euler( mesh, ice, ocean, laddie, tl, dt)  
        CASE ('fbrk3')
          CALL integrate_fbrk3( mesh, ice, ocean, laddie, tl, dt)  
      END SELECT

      ! Display or save fields
      ! TODO
      IF (par%master) THEN
        WRITE( *, "(F8.3,A,F8.3,A,F12.7,A,F8.3,A,F8.3)") tl/sec_per_day, '  Dmean ', SUM(laddie%now%H)/SIZE(laddie%now%H), '  Meltmax', MAXVAL(laddie%melt), '   U', MAXVAL(laddie%now%U), '   Tmax', MAXVAL(laddie%now%T)
      END IF     

    END DO !DO WHILE (tl <= C%time_duration_laddie)

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
    INTEGER                                               :: vi, ti
 
    ! Add routine to path
    CALL init_routine( routine_name)

    ! Update diffusive terms based on now time step
    CALL update_diffusive_terms( mesh, ice, laddie)

    ! Integrate H 1 time step
    CALL compute_H_npx( mesh, ice, laddie, laddie%now, laddie%np1, dt)

    ! Map H and Hstar to b grid and c grid
    CALL map_H_a_b( mesh, laddie, laddie%np1%H, laddie%np1%H_b)
    CALL map_H_a_c( mesh, laddie, laddie%np1%H, laddie%np1%H_c)
    CALL map_H_a_b( mesh, laddie, laddie%now%H, laddie%now%H_b)
    CALL map_H_a_c( mesh, laddie, laddie%now%H, laddie%now%H_c)

    ! Integrate U and V 1 time step
    CALL compute_UV_npx( mesh, ice, laddie, laddie%now, laddie%np1, laddie%now%H_b, dt, .true.)

    ! Integrate T and S 1 time step
    CALL compute_TS_npx( mesh, ice, laddie, laddie%now, laddie%np1, dt, .true.)

    ! Update secondary fields
    CALL update_secondary_fields( mesh, ice, ocean, laddie, laddie%np1, laddie%now%H, laddie%now%H_b, laddie%now%H_c)

    ! == Move time ==
    ! Increase laddie time
    tl = tl + C%dt_laddie

    ! Move main variables by 1 time step
    DO vi = mesh%vi1, mesh%vi2
      IF (laddie%mask_a( vi)) THEN
        laddie%now%H( vi) = laddie%np1%H( vi)
        laddie%now%T( vi) = laddie%np1%T( vi)
        laddie%now%S( vi) = laddie%np1%S( vi)
      END IF
    END DO

    ! Move velocities by 1 time step
    DO ti = mesh%ti1, mesh%ti2
      IF (laddie%mask_b( ti)) THEN
        laddie%now%U( ti) = laddie%np1%U( ti)
        laddie%now%V( ti) = laddie%np1%V( ti)
      END IF
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE integrate_euler

  SUBROUTINE integrate_fbrk3( mesh, ice, ocean, laddie, tl, dt)
    ! Integrate 1 timestep Forward-Backward Runge Kutta 3 scheme 

    ! In- and output variables

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_ocean_model),                 INTENT(IN)    :: ocean
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie
    REAL(dp),                               INTENT(INOUT) :: tl
    REAL(dp),                               INTENT(IN)    :: dt

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'integrate_fbrk3'
    INTEGER                                               :: vi, ti
    REAL(dp), PARAMETER                                   :: beta1 = 0.500_dp
    REAL(dp), PARAMETER                                   :: beta2 = 0.500_dp
    REAL(dp), PARAMETER                                   :: beta3 = 0.344_dp
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2)                :: Hstar, Hstarstar, Hstarstarstar
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2)                :: Hstar_b, Hstarstar_b, Hstarstarstar_b
    REAL(dp), DIMENSION(mesh%ei1:mesh%ei2)                :: Hstar_c, Hstarstar_c, Hstarstarstar_c
 
    ! Add routine to path
    CALL init_routine( routine_name)

    ! Update diffusive terms based on now time step
    ! Possibly only add during stage 3
    CALL update_diffusive_terms( mesh, ice, laddie)

    ! == Stage 1: explicit 1/3 timestep ==
    ! == RHS terms defined at n ==========
    ! ====================================
 
    ! Integrate H 1/3 time step
    CALL compute_H_npx( mesh, ice, laddie, laddie%now, laddie%np13, dt/3)

    ! Compute Hstar
    Hstar = 0.0_dp
    DO vi = mesh%vi1, mesh%vi2
      Hstar( vi) = beta1 * laddie%np13%H ( vi) + (1-beta1) * laddie%now%H( vi)
    END DO

    ! Map H and Hstar to b grid and c grid
    CALL map_H_a_b( mesh, laddie, laddie%np13%H, laddie%np13%H_b)
    CALL map_H_a_c( mesh, laddie, laddie%np13%H, laddie%np13%H_c)
    CALL map_H_a_b( mesh, laddie, Hstar, Hstar_b)
    CALL map_H_a_c( mesh, laddie, Hstar, Hstar_c)

    ! Integrate U and V 1/3 time step
    CALL compute_UV_npx( mesh, ice, laddie, laddie%now, laddie%np13, Hstar_b, dt/3, .false.)

    ! Integrate T and S 1/3 time step
    CALL compute_TS_npx( mesh, ice, laddie, laddie%now, laddie%np13, dt/3, .false.)

    ! Update secondary fields
    CALL update_secondary_fields( mesh, ice, ocean, laddie, laddie%np13, Hstar, Hstar_b, Hstar_c)

    ! == Stage 2: explicit 1/2 timestep ==
    ! == RHS terms defined at n + 1/3 ====
    ! ====================================

    ! Integrate H 1/2 time step
    CALL compute_H_npx( mesh, ice, laddie, laddie%np13, laddie%np12, dt/2)

    ! Compute Hstarstar
    Hstarstar = 0.0_dp
    DO vi = mesh%vi1, mesh%vi2
      Hstarstar( vi) = beta2 * laddie%np12%H ( vi) + (1-beta2) * laddie%now%H( vi)
    END DO

    ! Map H and Hstarstar to b and c grid
    CALL map_H_a_b( mesh, laddie, laddie%np12%H, laddie%np12%H_b)
    CALL map_H_a_c( mesh, laddie, laddie%np12%H, laddie%np12%H_c)
    CALL map_H_a_b( mesh, laddie, Hstarstar, Hstarstar_b)
    CALL map_H_a_c( mesh, laddie, Hstarstar, Hstarstar_c)

    ! Integrate U and V 1/2 time step
    CALL compute_UV_npx( mesh, ice, laddie, laddie%np13, laddie%np12, Hstarstar_b, dt/2, .false.)

    ! Integrate T and S 1/2 time step
    CALL compute_TS_npx( mesh, ice, laddie, laddie%np13, laddie%np12, dt/2, .false.)

    ! Update secondary fields
    CALL update_secondary_fields( mesh, ice, ocean, laddie, laddie%np12, Hstarstar, Hstarstar_b, Hstarstar_c)

    ! == Stage 3: explicit 1 timestep ====
    ! == RHS terms defined at n + 1/2 ====
    ! ====================================

    ! Integrate H 1 time step
    CALL compute_H_npx( mesh, ice, laddie, laddie%np12, laddie%np1, dt)

    ! Compute Hstarstarstar
    Hstarstarstar = 0.0_dp
    DO vi = mesh%vi1, mesh%vi2
      Hstarstarstar( vi) = beta3 * laddie%np1%H ( vi) + (1-2*beta3) * laddie%np12%H( vi) + beta3 * laddie%now%H( vi)
    END DO

    ! Map H and Hstarstarstar to b and cgrid
    CALL map_H_a_b( mesh, laddie, laddie%np1%H, laddie%np1%H_b)
    CALL map_H_a_c( mesh, laddie, laddie%np1%H, laddie%np1%H_c)
    CALL map_H_a_b( mesh, laddie, Hstarstarstar, Hstarstarstar_b)
    CALL map_H_a_c( mesh, laddie, Hstarstarstar, Hstarstarstar_c)

    ! Integrate U and V 1 time step
    CALL compute_UV_npx( mesh, ice, laddie, laddie%np12, laddie%np1, Hstarstarstar_b, dt, .true.)

    ! Integrate T and S 1 time step
    CALL compute_TS_npx( mesh, ice, laddie, laddie%np12, laddie%np1, dt, .true.)

    ! Update secondary fields
    CALL update_secondary_fields( mesh, ice, ocean, laddie, laddie%np1, Hstarstarstar, Hstarstarstar_b, Hstarstarstar_c)

    ! =============== 
    ! == Move time ==
    ! Increase laddie time
    tl = tl + C%dt_laddie

    ! Move main variables by 1 time step
    DO vi = mesh%vi1, mesh%vi2
      IF (laddie%mask_a( vi)) THEN
        laddie%now%H( vi) = laddie%np1%H( vi)
        laddie%now%T( vi) = laddie%np1%T( vi)
        laddie%now%S( vi) = laddie%np1%S( vi)
      END IF
    END DO

    ! Move velocities by 1 time step
    DO ti = mesh%ti1, mesh%ti2
      IF (laddie%mask_b( ti)) THEN
        laddie%now%U( ti) = laddie%np1%U( ti)
        laddie%now%V( ti) = laddie%np1%V( ti)
      END IF
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE integrate_fbrk3

  SUBROUTINE update_secondary_fields( mesh, ice, ocean, laddie, npx, Hstar, Hstar_b, Hstar_c)
    ! Update all secondary fields required for next iteration

    ! In- and output variables

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_ocean_model),                 INTENT(IN)    :: ocean
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie
    TYPE(type_laddie_timestep),             INTENT(IN)    :: npx          ! Reference timestep
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)    :: Hstar
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2), INTENT(IN)    :: Hstar_b
    REAL(dp), DIMENSION(mesh%ei1:mesh%ei2), INTENT(IN)    :: Hstar_c

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'update_secondary_fields'
    INTEGER                                               :: vi, ti
    TYPE(type_sparse_matrix_CSR_dp)                       :: M_divQ
    TYPE(type_sparse_matrix_CSR_dp)                       :: M_divQ_b
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2)                :: HstarT
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2)                :: HstarS

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Initialise ambient T and S
    CALL compute_ambient_TS( mesh, ice, ocean, laddie, Hstar)

    ! Compute freezing temperature
    CALL compute_freezing_temperature( mesh, ice, laddie, npx)

    ! Compute buoyancy
    CALL compute_buoyancy( mesh, ice, laddie, npx, Hstar)

    ! Bunch of mappings
    CALL map_a_b_2D( mesh, laddie%Hdrho_amb, laddie%Hdrho_amb_b)
    CALL map_b_a_2D( mesh, npx%U, laddie%U_a)
    CALL map_b_a_2D( mesh, npx%V, laddie%V_a)
    CALL map_laddie_velocities_from_b_to_c_2D( mesh, npx%U, npx%V, laddie%U_c, laddie%V_c)

    ! Bunch of derivatives
    CALL ddx_a_b_2D( mesh, laddie%drho_amb, laddie%ddrho_amb_dx_b)
    CALL ddy_a_b_2D( mesh, laddie%drho_amb, laddie%ddrho_amb_dy_b)
    CALL ddx_a_b_2D( mesh, Hstar, laddie%dH_dx_b)
    CALL ddy_a_b_2D( mesh, Hstar, laddie%dH_dy_b)

    ! Compute melt rate
    CALL compute_melt_rate( mesh, ice, ocean, laddie, npx, Hstar)
    
    ! Compute entrainment
    CALL compute_entrainment( mesh, ice, ocean, laddie, npx, Hstar)

    ! Compute thickness divergence
    CALL compute_divQH( mesh, laddie, npx, laddie%U_c, laddie%V_c, laddie%mask_a, laddie%mask_gr_a, laddie%mask_oc_a)

    ! Compute divergence of momentum
    CALL compute_divQUV( mesh, laddie, npx, laddie%U_c, laddie%V_c, Hstar_c, laddie%mask_b, laddie%mask_gl_b)

    ! Compute divergence of heat and salt
    CALL compute_divQTS( mesh, laddie, npx, laddie%U_c, laddie%V_c, Hstar, laddie%mask_a, laddie%mask_gr_a, laddie%mask_oc_a)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE update_secondary_fields

  SUBROUTINE update_diffusive_terms( mesh, ice, laddie)
    ! Update diffusivity and viscosity. Always based on now timestep for stability

    ! In- and output variables

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'update_diffusive_terms'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Map thickness to b grid
    CALL map_H_a_b( mesh, laddie, laddie%now%H, laddie%now%H_b)

    ! Map thickness to c grid
    CALL map_H_a_c( mesh, laddie, laddie%now%H, laddie%now%H_c)

    ! Compute diffusivities
    CALL compute_diffTS( mesh, ice, laddie)

    ! Compute viscosities
    CALL compute_viscUV( mesh, ice, laddie)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE update_diffusive_terms

END MODULE laddie_main

