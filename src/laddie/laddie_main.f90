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
  USE laddie_model_types                                     , ONLY: type_laddie_model
  USE ocean_model_types                                      , ONLY: type_ocean_model
  USE BMB_model_types                                        , ONLY: type_BMB_model
  USE reallocate_mod                                         , ONLY: reallocate_bounds
  USE laddie_utilities                                       , ONLY: compute_ambient_TS, allocate_laddie_model, &
                                                                     allocate_laddie_timestep, map_H_a_b, &
                                                                     calc_laddie_flux_divergence_matrix_upwind, &
                                                                     map_laddie_velocities_from_b_to_c_2D
  USE laddie_physics                                         , ONLY: compute_melt_rate, compute_entrainment, &
                                                                     compute_freezing_temperature, compute_buoyancy
  USE laddie_thickness                                       , ONLY: compute_H_np1 
  USE laddie_velocity                                        , ONLY: compute_UV_np1, compute_viscUV, compute_divQUV_centered
  USE laddie_tracers                                         , ONLY: compute_TS_np1, compute_diffTS
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
    INTEGER                                               :: vi, ti, nf, i
    REAL(dp)                                              :: tl               ! [s] Laddie time
    REAL(dp)                                              :: dt               ! [s] Laddie time step
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2)                :: Hstar            ! [m] Reference thickness in integration
    LOGICAL, DIMENSION(mesh%nV)                           :: mask_a_tot
    LOGICAL, DIMENSION(mesh%nV)                           :: mask_a_gr_tot
    LOGICAL, DIMENSION(mesh%nV)                           :: mask_a_oc_tot
 
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
      laddie%mask_a( vi)    = ice%mask_floating_ice( vi)
      laddie%mask_gr_a( vi) = ice%mask_grounded_ice( vi)
      laddie%mask_oc_a( vi) = ice%mask_icefree_ocean( vi)
    END DO

    ! Mask on b grid
    CALL gather_to_all_logical_1D( laddie%mask_a, mask_a_tot)
    CALL gather_to_all_logical_1D( laddie%mask_gr_a, mask_a_gr_tot)
    CALL gather_to_all_logical_1D( laddie%mask_oc_a, mask_a_oc_tot)

    DO ti = mesh%ti1, mesh%ti2
      ! Initialise as false to overwrite previous mask
      laddie%mask_b( ti)    = .false.
      laddie%mask_gl_b( ti) = .false.
      laddie%mask_cf_b( ti) = .false.
      ! Loop over connecing vertices and check whether they are floating
      DO i = 1, 3
        vi = mesh%Tri( ti, i)
        IF (mask_a_tot( vi)) THEN
          ! Set true if any of the three vertices is floating
          laddie%mask_b( ti) = .true.
        END IF
      END DO
      ! Loop over connecing vertices 
      DO i = 1, 3
        vi = mesh%Tri( ti, i)
        ! Check if any of them is grounded
        IF (mask_a_gr_tot( vi)) THEN
          ! Omit from mask if any of the three vertices is grounded
          laddie%mask_b( ti) = .false.
          ! Define as grounding line triangle
          laddie%mask_gl_b( ti) = .true.
        END IF
        ! Check if any of them is icefree ocean
        IF (mask_a_oc_tot( vi)) THEN
          ! Define as calving front triangle
          laddie%mask_cf_b( ti) = .true.
        END IF
      END DO
    END DO

    ! Extrapolate new cells
    ! TODO

    ! Set values to zero if outside laddie mask
    DO vi = mesh%vi1, mesh%vi2
      IF (.NOT. laddie%mask_a( vi)) THEN
        laddie%H( vi)     = 0.0_dp
        laddie%T( vi)     = 0.0_dp
        laddie%S( vi)     = 0.0_dp
        laddie%melt( vi)  = 0.0_dp
        laddie%entr( vi)  = 0.0_dp
      END IF
    END DO

    DO ti = mesh%ti1, mesh%ti2
      IF (.NOT. laddie%mask_b( ti)) THEN
        laddie%U( ti)     = 0.0_dp
        laddie%V( ti)     = 0.0_dp
      END IF
    END DO

    ! Update ice shelf draft gradients
    CALL ddx_a_b_2D( mesh, ice%Hib , laddie%dHib_dx_b)
    CALL ddy_a_b_2D( mesh, ice%Hib , laddie%dHib_dy_b)

    ! Update secondary fields
    CALL update_secondary_fields( mesh, ice, ocean, laddie, laddie%H)

    ! == Main time loop ==
    ! ====================

    DO WHILE (tl <= C%time_duration_laddie * sec_per_day)
      CALL integrate_euler( mesh, ice, ocean, laddie, tl, dt)  

      ! Display or save fields
      ! TODO
      IF (par%master) THEN
        WRITE( *, "(A,F8.3,A,F12.7,A,F8.3)") 'Dmax ', MAXVAL(laddie%H), '  Meltmax', MAXVAL(laddie%melt), '   U', MAXVAL(laddie%U)
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

    ! Allocate variables
    CALL allocate_laddie_timestep( mesh, laddie%np1)

    ! Mask on a grid
    DO vi = mesh%vi1, mesh%vi2
      laddie%mask_a( vi)  = ice%mask_floating_ice( vi)
    END DO

    ! Layer thickness 
    DO vi = mesh%vi1, mesh%vi2
       IF (laddie%mask_a( vi)) THEN
         laddie%H( vi)      = C%laddie_initial_thickness
         laddie%np1%H( vi)  = C%laddie_initial_thickness
       END IF
    END DO

    ! Initialise ambient T and S
    CALL compute_ambient_TS( mesh, ice, ocean, laddie, laddie%H)

    ! Initialise main T and S
    DO vi = mesh%vi1, mesh%vi2
       IF (laddie%mask_a( vi)) THEN
         laddie%T( vi)      = laddie%T_amb( vi) + C%laddie_initial_T_offset 
         laddie%np1%T( vi)  = laddie%T_amb( vi) + C%laddie_initial_T_offset
         laddie%S( vi)      = laddie%S_amb( vi) + C%laddie_initial_S_offset
         laddie%np1%S( vi)  = laddie%S_amb( vi) + C%laddie_initial_S_offset
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
 
    ! Integrate H 1 time step
    CALL compute_H_np1( mesh, ice, laddie, dt)

    ! Integrate U and V 1 time step
    CALL map_H_a_b( mesh, laddie, laddie%np1%H, laddie%np1%H_b)
    CALL compute_UV_np1( mesh, ice, laddie, dt)

    ! Integrate T and S 1 time step
    CALL compute_TS_np1( mesh, ice, laddie, dt)

    ! Update secondary fields
    CALL update_secondary_fields( mesh, ice, ocean, laddie, laddie%H)

    ! == Move time ==
    ! Increase laddie time
    tl = tl + C%dt_laddie

    ! Move main variables by 1 time step
    DO vi = mesh%vi1, mesh%vi2
      IF (laddie%mask_a( vi)) THEN
        laddie%H( vi) = laddie%np1%H( vi)
        laddie%T( vi) = laddie%np1%T( vi)
        laddie%S( vi) = laddie%np1%S( vi)
      END IF
    END DO

    ! Move velocities by 1 time step
    DO ti = mesh%ti1, mesh%ti2
      IF (laddie%mask_b( ti)) THEN
        laddie%U( ti) = laddie%np1%U( ti)
        laddie%V( ti) = laddie%np1%V( ti)
      END IF
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE integrate_euler

  SUBROUTINE update_secondary_fields( mesh, ice, ocean, laddie, Hstar)
    ! Update all secondary fields required for next iteration

    ! In- and output variables

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_ocean_model),                 INTENT(IN)    :: ocean
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)    :: Hstar

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'update_secondary_fields'
    INTEGER                                               :: vi, ti
    TYPE(type_sparse_matrix_CSR_dp)                       :: M_divQ
    TYPE(type_sparse_matrix_CSR_dp)                       :: M_divQ_b
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2)                :: HstarT
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2)                :: HstarS
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2)                :: HstarU
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2)                :: HstarV

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Initialise ambient T and S
    CALL compute_ambient_TS( mesh, ice, ocean, laddie, Hstar)

    ! Compute freezing temperature
    CALL compute_freezing_temperature( mesh, ice, laddie)

    ! Compute buoyancy
    CALL compute_buoyancy( mesh, ice, laddie, Hstar)

    ! Map buoyancy to b grid
    CALL map_a_b_2D( mesh, laddie%Hdrho_amb, laddie%Hdrho_amb_b)

    ! Map thickness to b grid
    CALL map_H_a_b( mesh, laddie, Hstar, laddie%H_b)

    ! Map thickness to c grid
    CALL map_a_c_2D( mesh, Hstar, laddie%H_c)

    ! Map next thickness to b grid
    CALL map_a_b_2D( mesh, laddie%np1%H, laddie%np1%H_b)

    ! Map detrainment to b grid
    CALL map_a_b_2D( mesh, laddie%detr, laddie%detr_b)

    ! Map velocities to a grid
    CALL map_b_a_2D( mesh, laddie%U, laddie%U_a)
    CALL map_b_a_2D( mesh, laddie%V, laddie%V_a)

    ! Update buoyancy derivatives
    CALL ddx_a_b_2D( mesh, laddie%drho_amb, laddie%ddrho_amb_dx_b)
    CALL ddy_a_b_2D( mesh, laddie%drho_amb, laddie%ddrho_amb_dy_b)

    ! Update thickness derivatives
    CALL ddx_a_b_2D( mesh, Hstar, laddie%dH_dx_b)
    CALL ddy_a_b_2D( mesh, Hstar, laddie%dH_dy_b)

    ! Compute melt rate
    CALL compute_melt_rate( mesh, ice, ocean, laddie, Hstar)
    
    ! Compute entrainment
    CALL compute_entrainment( mesh, ice, ocean, laddie, Hstar)

    ! Compute diffusivities
    CALL compute_diffTS( mesh, ice, laddie, Hstar)

    ! Compute viscosities
    CALL compute_viscUV( mesh, ice, laddie, Hstar)

    ! Get velocities on c grid
    CALL map_laddie_velocities_from_b_to_c_2D( mesh, laddie%U, laddie%V, laddie%U_c, laddie%V_c)

    ! Compute divergence matrix
    CALL calc_laddie_flux_divergence_matrix_upwind( mesh, laddie%U_c, laddie%V_c, laddie%mask_a, laddie%mask_gr_a, M_divQ)

    ! Compute thickness divergence
    CALL multiply_CSR_matrix_with_vector_1D( M_divQ, laddie%H, laddie%divQ)

    ! Compute Hstar * T
    DO vi = mesh%vi1, mesh%vi2
       IF (laddie%mask_a( vi)) THEN
         HstarT( vi) = Hstar( vi) * laddie%T( vi)
       END IF
    END DO
    ! Compute heat divergence
    CALL multiply_CSR_matrix_with_vector_1D( M_divQ, HstarT, laddie%divQT)

    ! Compute Hstar * S
    DO vi = mesh%vi1, mesh%vi2
       IF (laddie%mask_a( vi)) THEN
         HstarS( vi) = Hstar( vi) * laddie%S( vi)
       END IF
    END DO
    ! Compute salt divergence
    CALL multiply_CSR_matrix_with_vector_1D( M_divQ, HstarS, laddie%divQS)

    ! Compute divergence matrix on b grid
    CALL compute_divQUV_centered( mesh, laddie, laddie%U_c, laddie%V_c, laddie%H_c, laddie%mask_b, laddie%mask_gl_b)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE update_secondary_fields

END MODULE laddie_main

