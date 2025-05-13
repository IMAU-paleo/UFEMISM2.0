MODULE laddie_main

  ! The main laddie model module.

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE control_resources_and_error_messaging                  , ONLY: crash, init_routine, finalise_routine, colour_string, warning
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
                                                                     allocate_laddie_timestep, map_H_a_b, map_H_a_c
  USE laddie_thickness                                       , ONLY: compute_H_npx
  USE laddie_velocity                                        , ONLY: compute_UV_npx, compute_viscUV
  USE laddie_tracers                                         , ONLY: compute_TS_npx, compute_diffTS
  use laddie_operators                                       , only: update_laddie_operators
  USE mesh_utilities                                         , ONLY: extrapolate_Gaussian
  USE mpi_distributed_memory                                 , ONLY: gather_to_all
  use mpi_distributed_shared_memory, only: reallocate_dist_shared, hybrid_to_dist, dist_to_hybrid
  use mesh_halo_exchange, only: exchange_halos
  use laddie_output, only: create_laddie_output_fields_file, create_laddie_output_scalar_file, & 
      write_to_laddie_output_fields_file, write_to_laddie_output_scalar_file, buffer_laddie_scalars

  IMPLICIT NONE

CONTAINS

! ===== Main routines =====
! =========================

  SUBROUTINE run_laddie_model( mesh, ice, ocean, laddie, time, is_initial, region_name)
    ! Run the laddie model

    ! In- and output variables

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_ocean_model),                 INTENT(IN)    :: ocean
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie
    REAL(dp),                               INTENT(IN)    :: time
    logical,                                intent(in   ) :: is_initial
    character(len=3),                       intent(in   ) :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'run_laddie_model'
    INTEGER                                               :: vi, ti
    REAL(dp)                                              :: tl               ! [s] Laddie time
    REAL(dp)                                              :: dt               ! [s] Laddie time step
    REAL(dp)                                              :: duration         ! [days] Duration of run
    REAL(dp)                                              :: ref_time         ! [s] Reference time for writing
    REAL(dp), PARAMETER                                   :: time_relax_laddie = 0.02_dp ! [days]
    REAL(dp), PARAMETER                                   :: fac_dt_relax = 3.0_dp ! Reduction factor of time step
    REAL(dp)                                              :: time_to_write    ! [days] 
    REAL(dp)                                              :: last_write_time  ! [days] 


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
    laddie%now%H_c( mesh%ei1:mesh%ei2) = 0.0_dp

    ! == Update operators ==
    CALL update_laddie_operators( mesh, ice, laddie)

    ! == Main time loop ==
    ! ====================

    ! Determine run duration and apply offset for initial run
    if (is_initial) then
      duration = C%time_duration_laddie_init
      ref_time = time*sec_per_year - duration*sec_per_day
    else
      duration = C%time_duration_laddie
      ref_time = time*sec_per_year
    end if

    tl = 0.0_dp
    last_write_time = 0.0_dp
    time_to_write = C%time_interval_scalar_output

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

      ! Write to output
      if (C%do_write_laddie_output_fields) then
        call write_to_laddie_output_fields_file( mesh, laddie, region_name, ref_time + tl)
      end if

      if (C%do_write_laddie_output_scalar) then
        call buffer_laddie_scalars( mesh, laddie, ref_time + tl)

        ! Write if required
        if (tl > time_to_write * sec_per_day) then
          call write_to_laddie_output_scalar_file( laddie)
          last_write_time = time_to_write
          time_to_write = time_to_write + C%time_interval_scalar_output
        end if
      end if

    END DO !DO WHILE (tl < C%time_duration_laddie)

    ! Write any remaining buffered scalars
    ! if (par%primary .and. laddie%buffer%n > 0) then
    ! if (C%do_write_laddie_output_scalar .and. tl > last_write_time * sec_per_day) then
    !    call write_to_laddie_output_scalar_file( laddie)
    !  end if

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_laddie_model

  SUBROUTINE initialise_laddie_model( mesh, laddie, ocean, ice, region_name)
    ! Initialise the laddie model

    ! In- and output variables

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie
    TYPE(type_ocean_model),                 INTENT(IN)    :: ocean
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    character(len=3),                       intent(in   ) :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'initialise_laddie_model'
    INTEGER                                               :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to terminal
    IF (par%primary)  WRITE(*,"(A)") '   Initialising LADDIE model...'

    ! Allocate variables
    CALL allocate_laddie_model( mesh, laddie)

    ! == Update masks ==
    call update_laddie_masks( mesh, ice, laddie)

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

    ! Create output file
    if (C%do_write_laddie_output_fields) call create_laddie_output_fields_file( mesh, laddie, region_name)
    if (C%do_write_laddie_output_scalar) call create_laddie_output_scalar_file( laddie, region_name)

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
    call exchange_halos( mesh, npx%H)

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
    CALL update_diffusive_terms( mesh, laddie, laddie%now)

    ! Integrate U and V 1 time step
    CALL compute_UV_npx( mesh, ice, ocean, laddie, laddie%now, laddie%np1, laddie%now%H, dt, .true.)

    ! Integrate T and S 1 time step
    CALL compute_TS_npx( mesh, laddie, laddie%now, laddie%np1, laddie%now%H, dt, .true.)

    ! == Move time ==
    CALL move_laddie_timestep( mesh, laddie, tl, dt)

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
    integer                                               :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    ! == Stage 1: explicit 1/3 timestep ==
    ! == RHS terms defined at n ==========
    ! ====================================

    ! Integrate H 1/3 time step
    CALL compute_H_npx( mesh, ice, ocean, laddie, laddie%now, laddie%np13, time, dt/3)

    ! Compute Hstar
    do vi = mesh%vi1, mesh%vi2
      laddie%Hstar( vi) = C%laddie_fbrk3_beta1 * laddie%np13%H( vi) + (1-C%laddie_fbrk3_beta1) * laddie%now%H( vi)
    end do
    call exchange_halos( mesh, laddie%Hstar)

    ! Update diffusive terms
    CALL update_diffusive_terms( mesh, laddie, laddie%now)

    ! Integrate U and V 1/3 time step
    CALL compute_UV_npx( mesh, ice, ocean, laddie, laddie%now, laddie%np13, laddie%Hstar, dt/3, .false.)

    ! Integrate T and S 1/3 time step
    CALL compute_TS_npx( mesh, laddie, laddie%now, laddie%np13, laddie%now%H, dt/3, .false.)

    ! == Stage 2: explicit 1/2 timestep ==
    ! == RHS terms defined at n + 1/3 ====
    ! ====================================

    ! Integrate H 1/2 time step
    CALL compute_H_npx( mesh, ice, ocean, laddie, laddie%np13, laddie%np12, time, dt/2)

    ! Compute new Hstar
    do vi = mesh%vi1, mesh%vi2
      laddie%Hstar( vi) = C%laddie_fbrk3_beta2 * laddie%np12%H( vi) + (1-C%laddie_fbrk3_beta2) * laddie%now%H( vi)
    end do
    call exchange_halos( mesh, laddie%Hstar)

    ! Update diffusive terms
    !CALL update_diffusive_terms( mesh, laddie, laddie%np13)

    ! Integrate U and V 1/2 time step
    CALL compute_UV_npx( mesh, ice, ocean, laddie, laddie%np13, laddie%np12, laddie%Hstar, dt/2, .false.)

    ! Integrate T and S 1/2 time step
    CALL compute_TS_npx( mesh, laddie, laddie%np13, laddie%np12, laddie%np13%H, dt/2, .false.)

    ! == Stage 3: explicit 1 timestep ====
    ! == RHS terms defined at n + 1/2 ====
    ! ====================================

    ! Integrate H 1 time step
    CALL compute_H_npx( mesh, ice, ocean, laddie, laddie%np12, laddie%np1, time, dt)

    ! Compute new Hstar
    do vi = mesh%vi1, mesh%vi2
      laddie%Hstar( vi) = C%laddie_fbrk3_beta3 * laddie%np1%H( vi) + (1-2*C%laddie_fbrk3_beta3) * laddie%np12%H( vi) + C%laddie_fbrk3_beta3 * laddie%now%H( vi)
    end do
    call exchange_halos( mesh, laddie%Hstar)

    ! Update diffusive terms
    !CALL update_diffusive_terms( mesh, laddie, laddie%np12)

    ! Integrate U and V 1 time step
    CALL compute_UV_npx( mesh, ice, ocean, laddie, laddie%np12, laddie%np1, laddie%Hstar, dt, .true.)

    ! Integrate T and S 1 time step
    CALL compute_TS_npx( mesh, laddie, laddie%np12, laddie%np1, laddie%np12%H, dt, .true.)

    ! ===============
    ! == Move time ==
    CALL move_laddie_timestep( mesh, laddie, tl, dt)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE integrate_fbrk3

  SUBROUTINE move_laddie_timestep( mesh, laddie, tl, dt)
    ! Increase laddie time tl by timestep dt and overwrite now timestep

    ! In- and output variables
    type(type_mesh),                        intent(in   ) :: mesh
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
    laddie%now%H  ( mesh%vi1:mesh%vi2) = laddie%np1%H  ( mesh%vi1:mesh%vi2)
    laddie%now%T  ( mesh%vi1:mesh%vi2) = laddie%np1%T  ( mesh%vi1:mesh%vi2)
    laddie%now%S  ( mesh%vi1:mesh%vi2) = laddie%np1%S  ( mesh%vi1:mesh%vi2)
    laddie%now%U  ( mesh%ti1:mesh%ti2) = laddie%np1%U  ( mesh%ti1:mesh%ti2)
    laddie%now%V  ( mesh%ti1:mesh%ti2) = laddie%np1%V  ( mesh%ti1:mesh%ti2)
    laddie%now%H_b( mesh%ti1:mesh%ti2) = laddie%np1%H_b( mesh%ti1:mesh%ti2)
    laddie%now%H_c( mesh%ei1:mesh%ei2) = laddie%np1%H_c( mesh%ei1:mesh%ei2)
    laddie%now%U_a( mesh%vi1:mesh%vi2) = laddie%np1%U_a( mesh%vi1:mesh%vi2)
    laddie%now%U_c( mesh%ei1:mesh%ei2) = laddie%np1%U_c( mesh%ei1:mesh%ei2)
    laddie%now%V_a( mesh%vi1:mesh%vi2) = laddie%np1%V_a( mesh%vi1:mesh%vi2)
    laddie%now%V_c( mesh%ei1:mesh%ei2) = laddie%np1%V_c( mesh%ei1:mesh%ei2)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE move_laddie_timestep

  SUBROUTINE update_diffusive_terms( mesh, laddie, npxref)
    ! Update diffusivity and viscosity. Based on reference timestep npxref

    ! For stability, most studies base diffusive terms on the now timestep

    ! In- and output variables

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie
    TYPE(type_laddie_timestep),             INTENT(IN)    :: npxref

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'update_diffusive_terms'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Compute diffusivities
    CALL compute_diffTS( mesh, laddie, npxref)

    ! Compute viscosities
    CALL compute_viscUV( mesh, laddie, npxref)

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

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Mask on a-grid
    DO vi = mesh%vi1, mesh%vi2
      ! Check whether vertex on border
      IF (mesh%VBI( vi) > 0) THEN
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

    ! Mask on b-grid
    call exchange_halos( mesh, laddie%mask_a)
    call exchange_halos( mesh, laddie%mask_gr_a)
    call exchange_halos( mesh, laddie%mask_oc_a)

    DO ti = mesh%ti1, mesh%ti2
      ! Initialise as false to overwrite previous mask
      laddie%mask_b( ti)    = .false.
      laddie%mask_gl_b( ti) = .false.
      laddie%mask_cf_b( ti) = .false.
      laddie%mask_oc_b( ti) = .false.

      ! Define floating mask if any of the three vertices is floating
      DO i = 1, 3
        vi = mesh%Tri( ti, i)
        IF (laddie%mask_a( vi)) THEN
          ! Set true if any of the three vertices is floating
          laddie%mask_b( ti) = .true.
        END IF
      END DO

      ! Define grounding line triangles
      DO i = 1, 3
        vi = mesh%Tri( ti, i)
        ! Check if any connected vertex is grounded
        IF (laddie%mask_gr_a( vi)) THEN
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
          IF (laddie%mask_oc_a( vi)) THEN
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
        IF (laddie%mask_oc_a( vi)) THEN
          no = no + 1
        END IF
      END DO
      ! Check whether all vertices are icefree ocean
      IF (no == 3) THEN
        ! Define as ocean triangle
        laddie%mask_oc_b( ti) = .true.
      END IF

    END DO !ti = mesh%ti1, mesh%ti2

    call exchange_halos( mesh, laddie%mask_b)
    call exchange_halos( mesh, laddie%mask_gl_b)
    call exchange_halos( mesh, laddie%mask_cf_b)
    call exchange_halos( mesh, laddie%mask_oc_b)

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
    real(dp), dimension(:), pointer :: H_loc, T_loc, S_loc

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Initialise mask
    mask = 0

    ! Determine mask for seed (2: previously floating cells), fill (1: new floating cells), or ignore (0: grounded/ocean)
    DO vi = mesh%vi1, mesh%vi2
      ! Skip if vertex is at border
      IF (mesh%VBI( vi) > 0) CYCLE

      IF (C%choice_calving_law == 'threshold_thickness') THEN
        IF (ice%Hi( vi) < C%calving_threshold_thickness_shelf .and. ice%mask_floating_ice( vi)) CYCLE
      ELSE
        IF (ice%Hi( vi) < 1.0 .and. ice%mask_floating_ice( vi)) CYCLE
      END IF

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
    H_loc => laddie%now%H( mesh%vi1:mesh%vi2)
    T_loc => laddie%now%T( mesh%vi1:mesh%vi2)
    S_loc => laddie%now%S( mesh%vi1:mesh%vi2)
    CALL extrapolate_Gaussian( mesh, mask, H_loc, sigma)
    CALL extrapolate_Gaussian( mesh, mask, T_loc, sigma)
    CALL extrapolate_Gaussian( mesh, mask, S_loc, sigma)
    nullify( H_loc)
    nullify( T_loc)
    nullify( S_loc)

    ! The above should ensure that all (newly) floating vertices have a non-zero thickness
    ! In case the extrapolation did not cover this, apply a backup check to set values
    ! at non-zero initialisation
    DO vi = mesh%vi1, mesh%vi2
      ! Skip if vertex is at border
      IF (mesh%VBI( vi) > 0) CYCLE

      IF (C%choice_calving_law == 'threshold_thickness') THEN
        IF (ice%Hi( vi) < C%calving_threshold_thickness_shelf .and. ice%mask_floating_ice( vi)) CYCLE
      ELSE
        IF (ice%Hi( vi) < 1.0 .and. ice%mask_floating_ice( vi)) CYCLE
      END IF

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

  SUBROUTINE remap_laddie_model( mesh_old, mesh_new, ice, ocean, laddie, time, region_name)
    ! Reallocate and remap laddie variables

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                        INTENT(IN)    :: mesh_new
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_ocean_model),                 INTENT(IN)    :: ocean
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie
    REAL(dp),                               INTENT(IN)    :: time
    character(len=3),                       intent(in   ) :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'remap_laddie_model'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! == Regular variables ==

    ! Thickness
    call reallocate_dist_shared( laddie%dH_dt,          laddie%wdH_dt,          mesh_new%pai_V%n_nih)
    laddie%dH_dt         ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%dH_dt

    ! Temperatures
    call reallocate_dist_shared( laddie%T_amb,          laddie%wT_amb,          mesh_new%pai_V%n_nih)
    call reallocate_dist_shared( laddie%T_base,         laddie%wT_base,         mesh_new%pai_V%n_nih)
    call reallocate_dist_shared( laddie%T_freeze,       laddie%wT_freeze,       mesh_new%pai_V%n_nih)
    laddie%T_amb         ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%T_amb
    laddie%T_base        ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%T_base
    laddie%T_freeze      ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%T_freeze

    ! Salinities
    call reallocate_dist_shared( laddie%S_amb,          laddie%wS_amb,          mesh_new%pai_V%n_nih)
    call reallocate_dist_shared( laddie%S_base,         laddie%wS_base,         mesh_new%pai_V%n_nih)
    laddie%S_amb         ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%S_amb
    laddie%S_base        ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%S_base

    ! Densities and buoyancies
    call reallocate_dist_shared( laddie%rho,            laddie%wrho,            mesh_new%pai_V%n_nih)
    call reallocate_dist_shared( laddie%rho_amb,        laddie%wrho_amb,        mesh_new%pai_V%n_nih)
    call reallocate_dist_shared( laddie%drho_amb,       laddie%wdrho_amb,       mesh_new%pai_V%n_nih)
    call reallocate_dist_shared( laddie%Hdrho_amb,      laddie%wHdrho_amb,      mesh_new%pai_V%n_nih)
    call reallocate_dist_shared( laddie%Hdrho_amb_b,    laddie%wHdrho_amb_b,    mesh_new%pai_Tri%n_nih)
    call reallocate_dist_shared( laddie%drho_base,      laddie%wdrho_base,      mesh_new%pai_V%n_nih)
    laddie%rho           ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%rho
    laddie%rho_amb       ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%rho_amb
    laddie%drho_amb      ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%drho_amb
    laddie%Hdrho_amb     ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%Hdrho_amb
    laddie%Hdrho_amb_b   ( mesh_new%pai_Tri%i1_nih:mesh_new%pai_Tri%i2_nih) => laddie%Hdrho_amb_b
    laddie%drho_base     ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%drho_base

    ! Friction velocity
    call reallocate_dist_shared( laddie%u_star,         laddie%wu_star,         mesh_new%pai_V%n_nih)
    laddie%u_star        ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%u_star

    ! Physical parameter fields
    call reallocate_dist_shared( laddie%gamma_T,        laddie%wgamma_T,        mesh_new%pai_V%n_nih)
    call reallocate_dist_shared( laddie%gamma_S,        laddie%wgamma_S,        mesh_new%pai_V%n_nih)
    call reallocate_dist_shared( laddie%A_h,            laddie%wA_h,            mesh_new%pai_Tri%n_nih)
    call reallocate_dist_shared( laddie%K_h,            laddie%wK_h,            mesh_new%pai_V%n_nih)
    laddie%gamma_T       ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%gamma_T
    laddie%gamma_S       ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%gamma_S
    laddie%A_h           ( mesh_new%pai_Tri%i1_nih:mesh_new%pai_Tri%i2_nih) => laddie%A_h
    laddie%K_h           ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%K_h

    ! Vertical rates
    call reallocate_dist_shared( laddie%melt,           laddie%wmelt,           mesh_new%pai_V%n_nih)
    call reallocate_dist_shared( laddie%entr,           laddie%wentr,           mesh_new%pai_V%n_nih)
    call reallocate_dist_shared( laddie%entr_dmin,      laddie%wentr_dmin,      mesh_new%pai_V%n_nih)
    call reallocate_dist_shared( laddie%detr,           laddie%wdetr,           mesh_new%pai_V%n_nih)
    call reallocate_dist_shared( laddie%entr_tot,       laddie%wentr_tot,       mesh_new%pai_V%n_nih)
    laddie%melt          ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%melt
    laddie%entr          ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%entr
    laddie%entr_dmin     ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%entr_dmin
    laddie%detr          ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%detr
    laddie%entr_tot      ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%entr_tot

    ! Horizontal fluxes
    call reallocate_dist_shared( laddie%divQH,          laddie%wdivQH,          mesh_new%pai_V%n_nih)
    call reallocate_dist_shared( laddie%divQU,          laddie%wdivQU,          mesh_new%pai_Tri%n_nih)
    call reallocate_dist_shared( laddie%divQV,          laddie%wdivQV,          mesh_new%pai_Tri%n_nih)
    call reallocate_dist_shared( laddie%divQT,          laddie%wdivQT,          mesh_new%pai_V%n_nih)
    call reallocate_dist_shared( laddie%divQS,          laddie%wdivQS,          mesh_new%pai_V%n_nih)
    laddie%divQH         ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%divQH
    laddie%divQU         ( mesh_new%pai_Tri%i1_nih:mesh_new%pai_Tri%i2_nih) => laddie%divQU
    laddie%divQV         ( mesh_new%pai_Tri%i1_nih:mesh_new%pai_Tri%i2_nih) => laddie%divQV
    laddie%divQT         ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%divQT
    laddie%divQS         ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%divQS

    ! Viscosities
    call reallocate_dist_shared( laddie%viscU,          laddie%wviscU,          mesh_new%pai_Tri%n_nih)
    call reallocate_dist_shared( laddie%viscV,          laddie%wviscV,          mesh_new%pai_Tri%n_nih)
    laddie%viscU         ( mesh_new%pai_Tri%i1_nih:mesh_new%pai_Tri%i2_nih) => laddie%viscU
    laddie%viscV         ( mesh_new%pai_Tri%i1_nih:mesh_new%pai_Tri%i2_nih) => laddie%viscV

    ! Diffusivities
    call reallocate_dist_shared( laddie%diffT,          laddie%wdiffT,          mesh_new%pai_V%n_nih)
    call reallocate_dist_shared( laddie%diffS,          laddie%wdiffS,          mesh_new%pai_V%n_nih)
    laddie%diffT         ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%diffT
    laddie%diffS         ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%diffS

    ! RHS terms
    call reallocate_dist_shared( laddie%ddrho_amb_dx_b, laddie%wddrho_amb_dx_b, mesh_new%pai_Tri%n_nih)
    call reallocate_dist_shared( laddie%ddrho_amb_dy_b, laddie%wddrho_amb_dy_b, mesh_new%pai_Tri%n_nih)
    call reallocate_dist_shared( laddie%dH_dx_b,        laddie%wdH_dx_b,        mesh_new%pai_Tri%n_nih)
    call reallocate_dist_shared( laddie%dH_dy_b,        laddie%wdH_dy_b,        mesh_new%pai_Tri%n_nih)
    call reallocate_dist_shared( laddie%detr_b,         laddie%wdetr_b,         mesh_new%pai_Tri%n_nih)
    laddie%ddrho_amb_dx_b( mesh_new%pai_Tri%i1_nih:mesh_new%pai_Tri%i2_nih) => laddie%ddrho_amb_dx_b
    laddie%ddrho_amb_dy_b( mesh_new%pai_Tri%i1_nih:mesh_new%pai_Tri%i2_nih) => laddie%ddrho_amb_dy_b
    laddie%dH_dx_b       ( mesh_new%pai_Tri%i1_nih:mesh_new%pai_Tri%i2_nih) => laddie%dH_dx_b
    laddie%dH_dy_b       ( mesh_new%pai_Tri%i1_nih:mesh_new%pai_Tri%i2_nih) => laddie%dH_dy_b
    laddie%detr_b        ( mesh_new%pai_Tri%i1_nih:mesh_new%pai_Tri%i2_nih) => laddie%detr_b

    ! Forward-Backward Runge-Kutta 3 scheme
    call reallocate_dist_shared( laddie%Hstar         , laddie%wHstar         , mesh_new%pai_V%n_nih  )
    laddie%Hstar         ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%Hstar

    ! Mapped variables
    call reallocate_dist_shared( laddie%H_c           , laddie%wH_c           , mesh_new%pai_E%n_nih  )
    call reallocate_dist_shared( laddie%Hstar_b       , laddie%wHstar_b       , mesh_new%pai_Tri%n_nih)
    call reallocate_dist_shared( laddie%Hstar_c       , laddie%wHstar_c       , mesh_new%pai_E%n_nih  )
    laddie%H_c           ( mesh_new%pai_E%i1_nih  :mesh_new%pai_E%i2_nih  ) => laddie%H_c
    laddie%Hstar_b       ( mesh_new%pai_Tri%i1_nih:mesh_new%pai_Tri%i2_nih) => laddie%Hstar_b
    laddie%Hstar_c       ( mesh_new%pai_E%i1_nih  :mesh_new%pai_E%i2_nih  ) => laddie%Hstar_c

    ! Masks
    call reallocate_dist_shared( laddie%mask_a,         laddie%wmask_a,         mesh_new%pai_V%n_nih)
    call reallocate_dist_shared( laddie%mask_gr_a,      laddie%wmask_gr_a,      mesh_new%pai_V%n_nih)
    call reallocate_dist_shared( laddie%mask_oc_a,      laddie%wmask_oc_a,      mesh_new%pai_V%n_nih)
    call reallocate_dist_shared( laddie%mask_b,         laddie%wmask_b,         mesh_new%pai_Tri%n_nih)
    call reallocate_dist_shared( laddie%mask_gl_b,      laddie%wmask_gl_b,      mesh_new%pai_Tri%n_nih)
    call reallocate_dist_shared( laddie%mask_cf_b,      laddie%wmask_cf_b,      mesh_new%pai_Tri%n_nih)
    call reallocate_dist_shared( laddie%mask_oc_b,      laddie%wmask_oc_b,      mesh_new%pai_Tri%n_nih)
    laddie%mask_a        ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%mask_a
    laddie%mask_gr_a     ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%mask_gr_a
    laddie%mask_oc_a     ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%mask_oc_a
    laddie%mask_b        ( mesh_new%pai_Tri%i1_nih:mesh_new%pai_Tri%i2_nih) => laddie%mask_b
    laddie%mask_gl_b     ( mesh_new%pai_Tri%i1_nih:mesh_new%pai_Tri%i2_nih) => laddie%mask_gl_b
    laddie%mask_cf_b     ( mesh_new%pai_Tri%i1_nih:mesh_new%pai_Tri%i2_nih) => laddie%mask_cf_b
    laddie%mask_oc_b     ( mesh_new%pai_Tri%i1_nih:mesh_new%pai_Tri%i2_nih) => laddie%mask_oc_b

    ! == Re-initialise masks ==
    CALL update_laddie_masks( mesh_new, ice, laddie)

    ! == Update operators ==
    CALL update_laddie_operators( mesh_new, ice, laddie)

    ! == Timestep variables ==
    CALL remap_laddie_timestep( mesh_old, mesh_new, laddie, laddie%now)

    SELECT CASE(C%choice_laddie_integration_scheme)
      CASE DEFAULT
        CALL crash('unknown choice_laddie_integration_scheme "' // TRIM( C%choice_laddie_integration_scheme) // '"')
      CASE ('euler')
        CALL remap_laddie_timestep( mesh_old, mesh_new, laddie, laddie%np1)
      CASE ('fbrk3')
        CALL remap_laddie_timestep( mesh_old, mesh_new, laddie, laddie%np13)
        CALL remap_laddie_timestep( mesh_old, mesh_new, laddie, laddie%np12)
        CALL remap_laddie_timestep( mesh_old, mesh_new, laddie, laddie%np1)
    END SELECT

    ! == Re-initialise ==
    CALL run_laddie_model( mesh_new, ice, ocean, laddie, time, .TRUE., region_name)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_laddie_model

  SUBROUTINE remap_laddie_timestep( mesh_old, mesh_new, laddie, npx)
    ! Remap laddie timestep

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                        INTENT(IN)    :: mesh_new
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie
    TYPE(type_laddie_timestep),             INTENT(INOUT) :: npx

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'remap_laddie_timestep'
    real(dp), dimension(:), allocatable                   :: d_loc

    ! Add routine to path
    CALL init_routine( routine_name)

    ! DENK DROM - this should be cleaned up once the remapping code is ported to hybrid memory
    allocate( d_loc( mesh_old%vi1:mesh_old%vi2), source = 0._dp)
    call hybrid_to_dist( mesh_old%pai_V, npx%H, d_loc)
    CALL map_from_mesh_to_mesh_with_reallocation_2D( mesh_old, mesh_new, d_loc, '2nd_order_conservative')
    call reallocate_dist_shared( npx%H, npx%wH, mesh_new%pai_V%n_nih)
    call dist_to_hybrid( mesh_new%pai_V, d_loc, npx%H)
    deallocate( d_loc)

    call reallocate_dist_shared( npx%H_b, npx%wH_b, mesh_new%pai_Tri%n_nih)
    call reallocate_dist_shared( npx%H_c, npx%wH_c, mesh_new%pai_E%n_nih)
    call reallocate_dist_shared( npx%U,   npx%wU,   mesh_new%pai_Tri%n_nih)
    call reallocate_dist_shared( npx%U_a, npx%wU_a, mesh_new%pai_V%n_nih)
    call reallocate_dist_shared( npx%U_c, npx%wU_c, mesh_new%pai_E%n_nih)
    call reallocate_dist_shared( npx%V,   npx%wV,   mesh_new%pai_Tri%n_nih)
    call reallocate_dist_shared( npx%V_a, npx%wV_a, mesh_new%pai_V%n_nih)
    call reallocate_dist_shared( npx%V_c, npx%wV_c, mesh_new%pai_E%n_nih)

    npx%H  ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => npx%H
    npx%H_b( mesh_new%pai_Tri%i1_nih:mesh_new%pai_Tri%i2_nih) => npx%H_b
    npx%H_c( mesh_new%pai_E%i1_nih  :mesh_new%pai_E%i2_nih  ) => npx%H_c
    npx%U  ( mesh_new%pai_Tri%i1_nih:mesh_new%pai_Tri%i2_nih) => npx%U
    npx%U_a( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => npx%U_a
    npx%U_c( mesh_new%pai_E%i1_nih  :mesh_new%pai_E%i2_nih  ) => npx%U_c
    npx%V  ( mesh_new%pai_Tri%i1_nih:mesh_new%pai_Tri%i2_nih) => npx%V
    npx%V_a( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => npx%V_a
    npx%V_c( mesh_new%pai_E%i1_nih  :mesh_new%pai_E%i2_nih  ) => npx%V_c
    npx%T  ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => npx%T
    npx%S  ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => npx%S

    allocate( d_loc( mesh_old%vi1:mesh_old%vi2), source = 0._dp)
    call hybrid_to_dist( mesh_old%pai_V, npx%T, d_loc)
    CALL map_from_mesh_to_mesh_with_reallocation_2D( mesh_old, mesh_new, d_loc, '2nd_order_conservative')
    call reallocate_dist_shared( npx%T, npx%wT, mesh_new%pai_V%n_nih)
    call dist_to_hybrid( mesh_new%pai_V, d_loc, npx%T)
    deallocate( d_loc)

    allocate( d_loc( mesh_old%vi1:mesh_old%vi2), source = 0._dp)
    call hybrid_to_dist( mesh_old%pai_V, npx%S, d_loc)
    CALL map_from_mesh_to_mesh_with_reallocation_2D( mesh_old, mesh_new, d_loc, '2nd_order_conservative')
    call reallocate_dist_shared( npx%S, npx%wS, mesh_new%pai_V%n_nih)
    call dist_to_hybrid( mesh_new%pai_V, d_loc, npx%S)
    deallocate( d_loc)

    ! == Re-initialise ==

    ! Layer thickness on b and c grid
    CALL map_H_a_b( mesh_new, laddie, npx%H, npx%H_b)
    CALL map_H_a_c( mesh_new, laddie, npx%H, npx%H_c)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_laddie_timestep

END MODULE laddie_main

