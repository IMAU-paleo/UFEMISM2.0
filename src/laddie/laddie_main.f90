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
  USE remapping_main                                         , ONLY: map_from_mesh_to_mesh_with_reallocation_2D, &
                                                                     map_from_mesh_tri_to_mesh_tri_with_reallocation_2D
  USE laddie_utilities                                       , ONLY: compute_ambient_TS, allocate_laddie_model, &
                                                                     allocate_laddie_timestep, map_H_a_b, map_H_a_c
  use laddie_operators                                       , only: update_laddie_operators
  use laddie_velocity                                        , only: map_UV_b_c
  use laddie_physics                                         , only: compute_subglacial_discharge
  USE mesh_utilities                                         , ONLY: extrapolate_Gaussian
  use mesh_disc_apply_operators                              , only: map_b_a_2D
  use mesh_integrate_over_domain                             , only: integrate_over_domain, calc_and_print_min_mean_max
  USE mpi_distributed_memory                                 , ONLY: gather_to_all
  use mpi_distributed_shared_memory, only: reallocate_dist_shared, hybrid_to_dist, dist_to_hybrid
  use mesh_halo_exchange, only: exchange_halos
  use laddie_output, only: create_laddie_output_fields_file, create_laddie_output_scalar_file, &
      write_to_laddie_output_fields_file, write_to_laddie_output_scalar_file, buffer_laddie_scalars
  use laddie_integration, only: integrate_euler, integrate_fbrk3, integrate_lfra, move_laddie_timestep
  use mesh_repartitioning, only: repartition_mesh, repartition
  use mesh_memory, only: deallocate_mesh
  use checksum_mod, only: checksum

  IMPLICIT NONE

CONTAINS

! ===== Main routines =====
! =========================

  subroutine run_laddie_model( mesh, ice, ocean, laddie, time, is_initial, region_name)
    ! Run the laddie model

    ! In/output variables
    type(type_mesh),         intent(in   ) :: mesh
    type(type_ice_model),    intent(in   ) :: ice
    type(type_ocean_model),  intent(in   ) :: ocean
    type(type_laddie_model), intent(inout) :: laddie
    real(dp),                intent(in   ) :: time
    logical,                 intent(in   ) :: is_initial
    character(len=3),        intent(in   ) :: region_name

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_laddie_model'
    type(type_mesh)                :: mesh_repartitioned

    ! Add routine to path
    call init_routine( routine_name)

    call update_laddie_forcing( mesh, ice, ocean, laddie)

    ! Only in first time step
    SELECT CASE (C%choice_laddie_SGD)
      CASE DEFAULT
        CALL crash('unknown choice_laddie_SGD "' // trim( C%choice_laddie_SGD) // '"!')
      CASE ('none')
        ! Do nothing
      CASE ('idealised')
        ! Compute SGD
        CALL compute_subglacial_discharge( mesh, laddie)
    END SELECT

    if (C%do_repartition_laddie) then
      ! Repartition the mesh so each process has (approximately)
      ! the same number of ice shelf vertices/triangles
      call repartition_mesh( mesh, mesh_repartitioned, laddie%mask_a, laddie%mask_b)

      ! Repartition laddie
      call repartition_laddie( mesh, mesh_repartitioned, laddie)

      ! Run laddie on the repartitioned mesh
      call run_laddie_model_leg( mesh_repartitioned, laddie, time, is_initial, region_name)

      ! Un-repartition laddie
      call repartition_laddie( mesh_repartitioned, mesh, laddie)
    else
      ! Run laddie on the original mesh
      call run_laddie_model_leg( mesh, laddie, time, is_initial, region_name)
    end if

    ! Clean up after yourself
    call deallocate_mesh( mesh_repartitioned)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_laddie_model

  subroutine run_laddie_model_leg( mesh, laddie, time, is_initial, region_name)
    ! Run one leg of the laddie model

    ! In/output variables
    type(type_mesh),         intent(in   ) :: mesh
    type(type_laddie_model), intent(inout) :: laddie
    real(dp),                intent(in   ) :: time
    logical,                 intent(in   ) :: is_initial
    character(len=3),        intent(in   ) :: region_name

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_laddie_model_leg'
    integer                        :: vi, ti
    real(dp)                       :: tl               ! [s] Laddie time
    real(dp)                       :: dt               ! [s] Laddie time step
    real(dp)                       :: duration         ! [days] Duration of run
    real(dp)                       :: ref_time         ! [s] Reference time for writing
    real(dp), parameter            :: time_relax_laddie = 0.02_dp ! [days]
    real(dp), parameter            :: fac_dt_relax = 3.0_dp ! Reduction factor of time step
    real(dp)                       :: time_to_write    ! [days]
    real(dp)                       :: last_write_time  ! [days]

    ! Add routine to path
    CALL init_routine( routine_name)

    ! == Preparation ==
    ! =================

    ! Extrapolate data into new cells
    CALL extrapolate_laddie_variables( mesh, laddie)

    ! == Update masks ==
    CALL update_laddie_masks( mesh, laddie)

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

    call checksum( laddie%now%H, 'laddie%now%H', mesh%pai_V)
    call checksum( laddie%now%T, 'laddie%now%T', mesh%pai_V)
    call checksum( laddie%now%S, 'laddie%now%S', mesh%pai_V)
    call checksum( laddie%melt , 'laddie%melt' , mesh%pai_V)
    call checksum( laddie%entr , 'laddie%entr' , mesh%pai_V)

    DO ti = mesh%ti1, mesh%ti2
      IF (.NOT. laddie%mask_b( ti)) THEN
        laddie%now%U( ti)     = 0.0_dp
        laddie%now%V( ti)     = 0.0_dp
        laddie%now%H_b( ti)   = 0.0_dp
      END IF
    END DO

    call checksum( laddie%now%U  , 'laddie%now%U'  , mesh%pai_Tri)
    call checksum( laddie%now%V  , 'laddie%now%V'  , mesh%pai_Tri)
    call checksum( laddie%now%H_b, 'laddie%now%H_b', mesh%pai_Tri)

    ! Simply set H_c zero everywhere, will be recomputed through mapping later
    laddie%now%H_c( mesh%ei1:mesh%ei2) = 0.0_dp

    ! == Update operators ==
    CALL update_laddie_operators( mesh, laddie)

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

    ! Perform first integration with half the time step for LFRA scheme
    dt = C%dt_laddie / fac_dt_relax
    if (C%choice_laddie_integration_scheme == 'lfra') then
      call integrate_lfra( mesh, laddie, tl, time, dt)
    end if

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
          CALL integrate_euler( mesh, laddie, tl, time, dt)
        CASE ('fbrk3')
          CALL integrate_fbrk3( mesh, laddie, tl, time, dt)
        CASE ('lfra')
          CALL integrate_lfra( mesh, laddie, tl, time, 2*dt)
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

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  end subroutine run_laddie_model_leg

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

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to terminal
    IF (par%primary)  WRITE(*,"(A)") '   Initialising LADDIE model...'

    ! Allocate variables
    CALL allocate_laddie_model( mesh, laddie)

    call update_laddie_forcing( mesh, ice, ocean, laddie)

    ! == Update masks ==
    call update_laddie_masks( mesh, laddie)

    ! == Update operators ==
    CALL update_laddie_operators( mesh, laddie)

    ! Initialise requested timesteps
    CALL initialise_laddie_model_timestep( mesh, laddie, laddie%now)

    SELECT CASE(C%choice_laddie_integration_scheme)
      CASE DEFAULT
        CALL crash('unknown choice_laddie_integration_scheme "' // TRIM( C%choice_laddie_integration_scheme) // '"')
      CASE ('euler')
        CALL initialise_laddie_model_timestep( mesh, laddie, laddie%np1)
      CASE ('fbrk3')
        CALL initialise_laddie_model_timestep( mesh, laddie, laddie%np13)
        CALL initialise_laddie_model_timestep( mesh, laddie, laddie%np12)
        CALL initialise_laddie_model_timestep( mesh, laddie, laddie%np1)
      CASE ('lfra')
        call crash('LeapFrog RobertAsselin scheme does not work yet, use euler or fbrk3')
        CALL initialise_laddie_model_timestep( mesh, laddie, laddie%nm1)
        CALL initialise_laddie_model_timestep( mesh, laddie, laddie%np1)
    END SELECT

    ! Create output file
    if (C%do_write_laddie_output_fields) call create_laddie_output_fields_file( mesh, laddie, region_name)
    if (C%do_write_laddie_output_scalar) call create_laddie_output_scalar_file( laddie, region_name)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_laddie_model

  SUBROUTINE initialise_laddie_model_timestep( mesh, laddie, npx)
    ! Initialise the laddie model for given timestep

    ! In- and output variables

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie
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
    call checksum( npx%H, 'npx%H', mesh%pai_V)

    ! Layer thickness on b and c grid
    CALL map_H_a_b( mesh, laddie, npx%H, npx%H_b)
    CALL map_H_a_c( mesh, laddie, npx%H, npx%H_c)
    call checksum( npx%H_b, 'npx%H_b', mesh%pai_Tri)
    call checksum( npx%H_c, 'npx%H_c', mesh%pai_E)

    ! Initialise ambient T and S
    CALL compute_ambient_TS( mesh, laddie, npx%H)

    ! Initialise main T and S
    DO vi = mesh%vi1, mesh%vi2
       IF (laddie%mask_a( vi)) THEN
         npx%T( vi)      = laddie%T_amb( vi) + C%laddie_initial_T_offset
         npx%S( vi)      = laddie%S_amb( vi) + C%laddie_initial_S_offset
       END IF
    END DO
    call checksum( npx%T, 'npx%T', mesh%pai_V)
    call checksum( npx%S, 'npx%S', mesh%pai_V)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_laddie_model_timestep

  SUBROUTINE update_laddie_masks( mesh, laddie)
    ! Update bunch of masks for laddie at the start of a new run

    ! In- and output variables

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
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
      ELSE IF (laddie%Hi( vi) < 1.0 .and. laddie%mask_floating_ice( vi)) THEN
        laddie%mask_a( vi)    = .false.
        laddie%mask_oc_a( vi) = .true.
      ELSE
        ! Inherit regular masks
        laddie%mask_a( vi)    = laddie%mask_floating_ice( vi)
        laddie%mask_gr_a( vi) = laddie%mask_grounded_ice( vi) .OR. laddie%mask_icefree_land( vi)
        laddie%mask_oc_a( vi) = laddie%mask_icefree_ocean( vi)
      END IF

      ! Define domain for area integration
      if (laddie%mask_a( vi)) then
        laddie%domain_a( vi) = 1.0_dp
      else
        laddie%domain_a( vi) = 0.0_dp
      end if
    END DO

    call exchange_halos( mesh, laddie%mask_a)
    call exchange_halos( mesh, laddie%mask_gr_a)
    call exchange_halos( mesh, laddie%mask_oc_a)
    call exchange_halos( mesh, laddie%domain_a)
    call checksum( laddie%mask_a   , 'laddie%mask_a'   , mesh%pai_V)
    call checksum( laddie%mask_gr_a, 'laddie%mask_gr_a', mesh%pai_V)
    call checksum( laddie%mask_oc_a, 'laddie%mask_oc_a', mesh%pai_V)
    call checksum( laddie%domain_a , 'laddie%domain_a' , mesh%pai_V)

    call integrate_over_domain( mesh, laddie%domain_a, laddie%area_a)
    call checksum( laddie%area_a, 'laddie%area_a')

    ! Mask on b-grid
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

      ! Define domain for area integration
      if (laddie%mask_b( ti)) then
        laddie%domain_b( ti) = 1.0_dp
      else
        laddie%domain_b( ti) = 0.0_dp
      end if
    END DO !ti = mesh%ti1, mesh%ti2

    call exchange_halos( mesh, laddie%mask_b)
    call exchange_halos( mesh, laddie%mask_gl_b)
    call exchange_halos( mesh, laddie%mask_cf_b)
    call exchange_halos( mesh, laddie%mask_oc_b)
    call exchange_halos( mesh, laddie%domain_b)
    call checksum( laddie%mask_b   , 'laddie%mask_b'   , mesh%pai_Tri)
    call checksum( laddie%mask_gl_b, 'laddie%mask_gl_b', mesh%pai_Tri)
    call checksum( laddie%mask_cf_b, 'laddie%mask_cf_b', mesh%pai_Tri)
    call checksum( laddie%mask_oc_b, 'laddie%mask_oc_b', mesh%pai_Tri)
    call checksum( laddie%domain_b , 'laddie%domain_b' , mesh%pai_Tri)

    call integrate_over_domain( mesh, laddie%domain_b, laddie%area_b)
    call checksum( laddie%area_b, 'laddie%area_b')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE update_laddie_masks

  SUBROUTINE extrapolate_laddie_variables( mesh, laddie)
    ! Update bunch of masks for laddie at the start of a new run

    ! In- and output variables

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
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
        IF (laddie%Hi( vi) < C%calving_threshold_thickness_shelf .and. laddie%mask_floating_ice( vi)) CYCLE
      ELSE
        IF (laddie%Hi( vi) < 1.0 .and. laddie%mask_floating_ice( vi)) CYCLE
      END IF

      ! Currently floating ice, so either seed or fill here
      IF (laddie%mask_floating_ice( vi)) THEN
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
    call checksum( laddie%now%H, 'laddie%now%H', mesh%pai_V)
    call checksum( laddie%now%T, 'laddie%now%T', mesh%pai_V)
    call checksum( laddie%now%S, 'laddie%now%S', mesh%pai_V)

    ! The above should ensure that all (newly) floating vertices have a non-zero thickness
    ! In case the extrapolation did not cover this, apply a backup check to set values
    ! at non-zero initialisation
    DO vi = mesh%vi1, mesh%vi2
      ! Skip if vertex is at border
      IF (mesh%VBI( vi) > 0) CYCLE

      IF (C%choice_calving_law == 'threshold_thickness') THEN
        IF (laddie%Hi( vi) < C%calving_threshold_thickness_shelf .and. laddie%mask_floating_ice( vi)) CYCLE
      ELSE
        IF (laddie%Hi( vi) < 1.0 .and. laddie%mask_floating_ice( vi)) CYCLE
      END IF

      ! Currently floating ice, so either seed or fill here
      IF (laddie%mask_floating_ice( vi)) THEN
        IF (laddie%now%H( vi) == 0.0_dp) THEN
          laddie%now%H( vi) = C%laddie_thickness_minimum
          laddie%now%T( vi) = laddie%T_amb( vi) + C%laddie_initial_T_offset
          laddie%now%S( vi) = laddie%S_amb( vi) + C%laddie_initial_S_offset
        END IF
      END IF
    END DO
    call checksum( laddie%now%H, 'laddie%now%H', mesh%pai_V)
    call checksum( laddie%now%T, 'laddie%now%T', mesh%pai_V)
    call checksum( laddie%now%S, 'laddie%now%S', mesh%pai_V)

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

    ! Forcing
    call reallocate_dist_shared( laddie%Hi                , laddie%wHi                , mesh_new%pai_V%n_nih)
    call reallocate_dist_shared( laddie%Hib               , laddie%wHib               , mesh_new%pai_V%n_nih)
    call reallocate_dist_shared( laddie%dHib_dx_b         , laddie%wdHib_dx_b         , mesh_new%pai_Tri%n_nih)
    call reallocate_dist_shared( laddie%dHib_dy_b         , laddie%wdHib_dy_b         , mesh_new%pai_Tri%n_nih)
    call reallocate_dist_shared( laddie%mask_icefree_land , laddie%wmask_icefree_land , mesh_new%pai_V%n_nih)
    call reallocate_dist_shared( laddie%mask_icefree_ocean, laddie%wmask_icefree_ocean, mesh_new%pai_V%n_nih)
    call reallocate_dist_shared( laddie%mask_grounded_ice , laddie%wmask_grounded_ice , mesh_new%pai_V%n_nih)
    call reallocate_dist_shared( laddie%mask_floating_ice , laddie%wmask_floating_ice , mesh_new%pai_V%n_nih)
    call reallocate_dist_shared( laddie%mask_gl_fl        , laddie%wmask_gl_fl        , mesh_new%pai_V%n_nih)
    call reallocate_dist_shared( laddie%mask_SGD          , laddie%wmask_SGD          , mesh_new%pai_V%n_nih)
    call reallocate_dist_shared( laddie%Ti                , laddie%wTi                , mesh_new%pai_V%n_nih, mesh_new%nz)
    call reallocate_dist_shared( laddie%T_ocean           , laddie%wT_ocean           , mesh_new%pai_V%n_nih, C%nz_ocean)
    call reallocate_dist_shared( laddie%S_ocean           , laddie%wS_ocean           , mesh_new%pai_V%n_nih, C%nz_ocean)
    laddie%Hi                ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih               ) => laddie%Hi
    laddie%Hib               ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih               ) => laddie%Hib
    laddie%dHib_dx_b         ( mesh_new%pai_Tri%i1_nih:mesh_new%pai_Tri%i2_nih             ) => laddie%dHib_dx_b
    laddie%dHib_dy_b         ( mesh_new%pai_Tri%i1_nih:mesh_new%pai_Tri%i2_nih             ) => laddie%dHib_dy_b
    laddie%mask_icefree_land ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih               ) => laddie%mask_icefree_land
    laddie%mask_icefree_ocean( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih               ) => laddie%mask_icefree_ocean
    laddie%mask_grounded_ice ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih               ) => laddie%mask_grounded_ice
    laddie%mask_floating_ice ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih               ) => laddie%mask_floating_ice
    laddie%mask_gl_fl        ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih               ) => laddie%mask_gl_fl
    laddie%mask_SGD          ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih               ) => laddie%mask_SGD
    laddie%Ti                ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih, 1:mesh_new%nz) => laddie%Ti
    laddie%T_ocean           ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih, 1:C%nz_ocean ) => laddie%T_ocean
    laddie%S_ocean           ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih, 1:C%nz_ocean ) => laddie%S_ocean

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
    call reallocate_dist_shared( laddie%SGD,            laddie%wSGD,            mesh_new%pai_V%n_nih)

    laddie%melt          ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%melt
    laddie%entr          ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%entr
    laddie%entr_dmin     ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%entr_dmin
    laddie%detr          ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%detr
    laddie%entr_tot      ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%entr_tot
    laddie%SGD           ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%SGD

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

    ! Domain
    call reallocate_dist_shared( laddie%domain_a      , laddie%wdomain_a      , mesh_new%pai_V%n_nih  )
    call reallocate_dist_shared( laddie%domain_b      , laddie%wdomain_b      , mesh_new%pai_Tri%n_nih)
    laddie%domain_a      ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%domain_a
    laddie%domain_b      ( mesh_new%pai_Tri%i1_nih:mesh_new%pai_Tri%i2_nih) => laddie%domain_b

    call update_laddie_forcing( mesh_new, ice, ocean, laddie)

    ! == Re-initialise masks ==
    CALL update_laddie_masks( mesh_new, laddie)

    ! == Update operators ==
    CALL update_laddie_operators( mesh_new, laddie)

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
      CASE ('lfra')
        CALL remap_laddie_timestep( mesh_old, mesh_new, laddie, laddie%nm1)
        CALL remap_laddie_timestep( mesh_old, mesh_new, laddie, laddie%np1)
    END SELECT

    ! == Re-initialise ==
    CALL run_laddie_model( mesh_new, ice, ocean, laddie, time, .FALSE., region_name)

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
    integer                                               :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    ! DENK DROM - this should be cleaned up once the remapping code is ported to hybrid memory
    allocate( d_loc( mesh_old%vi1:mesh_old%vi2), source = 0._dp)
    call hybrid_to_dist( mesh_old%pai_V, npx%H, d_loc)
    call map_from_mesh_to_mesh_with_reallocation_2D( mesh_old, mesh_new, C%output_dir, d_loc, '2nd_order_conservative')
    call reallocate_dist_shared( npx%H, npx%wH, mesh_new%pai_V%n_nih)
    call dist_to_hybrid( mesh_new%pai_V, d_loc, npx%H)
    deallocate( d_loc)

    allocate( d_loc( mesh_old%vi1:mesh_old%vi2), source = 0._dp)
    call hybrid_to_dist( mesh_old%pai_V, npx%T, d_loc)
    call map_from_mesh_to_mesh_with_reallocation_2D( mesh_old, mesh_new, C%output_dir, d_loc, '2nd_order_conservative')
    call reallocate_dist_shared( npx%T, npx%wT, mesh_new%pai_V%n_nih)
    call dist_to_hybrid( mesh_new%pai_V, d_loc, npx%T)
    deallocate( d_loc)

    allocate( d_loc( mesh_old%vi1:mesh_old%vi2), source = 0._dp)
    call hybrid_to_dist( mesh_old%pai_V, npx%S, d_loc)
    call map_from_mesh_to_mesh_with_reallocation_2D( mesh_old, mesh_new, C%output_dir, d_loc, '2nd_order_conservative')
    call reallocate_dist_shared( npx%S, npx%wS, mesh_new%pai_V%n_nih)
    call dist_to_hybrid( mesh_new%pai_V, d_loc, npx%S)
    deallocate( d_loc)

    allocate( d_loc( mesh_old%ti1:mesh_old%ti2), source = 0._dp)
    call hybrid_to_dist( mesh_old%pai_Tri, npx%U, d_loc)
    call map_from_mesh_tri_to_mesh_tri_with_reallocation_2D( mesh_old, mesh_new, C%output_dir, d_loc, '2nd_order_conservative')
    call reallocate_dist_shared( npx%U, npx%wU, mesh_new%pai_Tri%n_nih)
    call dist_to_hybrid( mesh_new%pai_Tri, d_loc, npx%U)
    deallocate( d_loc)

    allocate( d_loc( mesh_old%ti1:mesh_old%ti2), source = 0._dp)
    call hybrid_to_dist( mesh_old%pai_Tri, npx%V, d_loc)
    call map_from_mesh_tri_to_mesh_tri_with_reallocation_2D( mesh_old, mesh_new, C%output_dir, d_loc, '2nd_order_conservative')
    call reallocate_dist_shared( npx%V, npx%wV, mesh_new%pai_Tri%n_nih)
    call dist_to_hybrid( mesh_new%pai_Tri, d_loc, npx%V)
    deallocate( d_loc)

    call checksum( npx%H, 'npx%H', mesh_new%pai_V)
    call checksum( npx%T, 'npx%T', mesh_new%pai_V)
    call checksum( npx%S, 'npx%S', mesh_new%pai_V)
    call checksum( npx%U, 'npx%U', mesh_new%pai_Tri)
    call checksum( npx%V, 'npx%V', mesh_new%pai_Tri)

    call reallocate_dist_shared( npx%H_b, npx%wH_b, mesh_new%pai_Tri%n_nih)
    call reallocate_dist_shared( npx%H_c, npx%wH_c, mesh_new%pai_E%n_nih)
    call reallocate_dist_shared( npx%U_a, npx%wU_a, mesh_new%pai_V%n_nih)
    call reallocate_dist_shared( npx%U_c, npx%wU_c, mesh_new%pai_E%n_nih)
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

    do vi = mesh_new%vi1, mesh_new%vi2
      npx%H( vi) = max(npx%H( vi), C%laddie_thickness_minimum)
      npx%H( vi) = min(npx%H( vi), C%laddie_thickness_maximum)
    end do

    ! Layer thickness on b and c grid
    call map_H_a_b( mesh_new, laddie, npx%H, npx%H_b)
    call map_H_a_c( mesh_new, laddie, npx%H, npx%H_c)

    ! Velocities to a and c grid
    call map_b_a_2D( mesh_new, npx%U, npx%U_a, d_b_is_hybrid = .true., d_a_is_hybrid = .true.)
    call map_b_a_2D( mesh_new, npx%V, npx%V_a, d_b_is_hybrid = .true., d_a_is_hybrid = .true.)

    call map_UV_b_c( mesh_new, laddie, npx%U, npx%V, npx%U_c, npx%V_c)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_laddie_timestep

  subroutine update_laddie_forcing( mesh, ice, ocean, laddie)

    ! In/output variables
    type(type_mesh),         intent(in   ) :: mesh
    type(type_ice_model),    intent(in   ) :: ice
    type(type_ocean_model),  intent(in   ) :: ocean
    type(type_laddie_model), intent(inout) :: laddie

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'update_laddie_forcing'

    ! Add routine to path
    call init_routine( routine_name)

    laddie%Hi                ( mesh%vi1:mesh%vi2  ) = ice%Hi                ( mesh%vi1:mesh%vi2  )
    laddie%Hib               ( mesh%vi1:mesh%vi2  ) = ice%Hib               ( mesh%vi1:mesh%vi2  )
    laddie%dHib_dx_b         ( mesh%ti1:mesh%ti2  ) = ice%dHib_dx_b         ( mesh%ti1:mesh%ti2  )
    laddie%dHib_dy_b         ( mesh%ti1:mesh%ti2  ) = ice%dHib_dy_b         ( mesh%ti1:mesh%ti2  )
    laddie%mask_icefree_land ( mesh%vi1:mesh%vi2  ) = ice%mask_icefree_land ( mesh%vi1:mesh%vi2  )
    laddie%mask_icefree_ocean( mesh%vi1:mesh%vi2  ) = ice%mask_icefree_ocean( mesh%vi1:mesh%vi2  )
    laddie%mask_grounded_ice ( mesh%vi1:mesh%vi2  ) = ice%mask_grounded_ice ( mesh%vi1:mesh%vi2  )
    laddie%mask_floating_ice ( mesh%vi1:mesh%vi2  ) = ice%mask_floating_ice ( mesh%vi1:mesh%vi2  )

    laddie%mask_gl_fl        ( mesh%vi1:mesh%vi2  ) = ice%mask_gl_fl        ( mesh%vi1:mesh%vi2  )
    laddie%mask_SGD          ( mesh%vi1:mesh%vi2  ) = ice%mask_SGD          ( mesh%vi1:mesh%vi2  )
    
    laddie%Ti                ( mesh%vi1:mesh%vi2,:) = ice%Ti                ( mesh%vi1:mesh%vi2,:) - 273.15 ! [degC]
    laddie%T_ocean           ( mesh%vi1:mesh%vi2,:) = ocean%T               ( mesh%vi1:mesh%vi2,:)
    laddie%S_ocean           ( mesh%vi1:mesh%vi2,:) = ocean%S               ( mesh%vi1:mesh%vi2,:)

    call checksum( laddie%Hi                , 'laddie%Hi'                , mesh%pai_V)
    call checksum( laddie%Hib               , 'laddie%Hib'               , mesh%pai_V)
    call checksum( laddie%dHib_dx_b         , 'laddie%dHib_dx_b'         , mesh%pai_Tri)
    call checksum( laddie%dHib_dy_b         , 'laddie%dHib_dy_b'         , mesh%pai_Tri)
    call checksum( laddie%mask_icefree_land , 'laddie%mask_icefree_land' , mesh%pai_V)
    call checksum( laddie%mask_icefree_ocean, 'laddie%mask_icefree_ocean', mesh%pai_V)
    call checksum( laddie%mask_grounded_ice , 'laddie%mask_grounded_ice' , mesh%pai_V)
    call checksum( laddie%mask_floating_ice , 'laddie%mask_floating_ice' , mesh%pai_V)
    call checksum( laddie%mask_gl_fl        , 'laddie%mask_gl_fl'        , mesh%pai_V)
    call checksum( laddie%mask_SGD          , 'laddie%mask_SGD'          , mesh%pai_V)
    call checksum( laddie%Ti                , 'laddie%Ti'                , mesh%pai_V)
    call checksum( laddie%T_ocean           , 'laddie%T_ocean'           , mesh%pai_V)
    call checksum( laddie%S_ocean           , 'laddie%S_ocean'           , mesh%pai_V)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine update_laddie_forcing

  subroutine repartition_laddie( mesh_old, mesh_new, laddie)

    ! In/output variables
    type(type_mesh),         intent(in   ) :: mesh_old, mesh_new
    type(type_laddie_model), intent(inout) :: laddie

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'repartition_laddie'

    ! Add routine to path
    call init_routine( routine_name)

    ! Thickness
    call repartition( mesh_old, mesh_new, laddie%dH_dt         , laddie%wdH_dt         )    ! [m]             change
    laddie%dH_dt         ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%dH_dt

    ! Forcing
    call repartition( mesh_old, mesh_new, laddie%Hi                , laddie%wHi                )
    call repartition( mesh_old, mesh_new, laddie%Hib               , laddie%wHib               )
    call repartition( mesh_old, mesh_new, laddie%dHib_dx_b         , laddie%wdHib_dx_b         )
    call repartition( mesh_old, mesh_new, laddie%dHib_dy_b         , laddie%wdHib_dy_b         )
    call repartition( mesh_old, mesh_new, laddie%mask_icefree_land , laddie%wmask_icefree_land )
    call repartition( mesh_old, mesh_new, laddie%mask_icefree_ocean, laddie%wmask_icefree_ocean)
    call repartition( mesh_old, mesh_new, laddie%mask_grounded_ice , laddie%wmask_grounded_ice )
    call repartition( mesh_old, mesh_new, laddie%mask_floating_ice , laddie%wmask_floating_ice )
    call repartition( mesh_old, mesh_new, laddie%mask_gl_fl        , laddie%wmask_gl_fl        )
    call repartition( mesh_old, mesh_new, laddie%mask_SGD          , laddie%wmask_SGD          )
    call repartition( mesh_old, mesh_new, laddie%Ti                , laddie%wTi                )
    call repartition( mesh_old, mesh_new, laddie%T_ocean           , laddie%wT_ocean           )
    call repartition( mesh_old, mesh_new, laddie%S_ocean           , laddie%wS_ocean           )

    laddie%Hi                ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih               ) => laddie%Hi
    laddie%Hib               ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih               ) => laddie%Hib
    laddie%dHib_dx_b         ( mesh_new%pai_Tri%i1_nih:mesh_new%pai_Tri%i2_nih             ) => laddie%dHib_dx_b
    laddie%dHib_dy_b         ( mesh_new%pai_Tri%i1_nih:mesh_new%pai_Tri%i2_nih             ) => laddie%dHib_dy_b
    laddie%mask_icefree_land ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih               ) => laddie%mask_icefree_land
    laddie%mask_icefree_ocean( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih               ) => laddie%mask_icefree_ocean
    laddie%mask_grounded_ice ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih               ) => laddie%mask_grounded_ice
    laddie%mask_floating_ice ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih               ) => laddie%mask_floating_ice
    laddie%mask_gl_fl        ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih               ) => laddie%mask_gl_fl
    laddie%mask_SGD          ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih               ) => laddie%mask_SGD
    laddie%Ti                ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih, 1:mesh_new%nz) => laddie%Ti
    laddie%T_ocean           ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih, 1:C%nz_ocean ) => laddie%T_ocean
    laddie%S_ocean           ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih, 1:C%nz_ocean ) => laddie%S_ocean

    call checksum( laddie%Hi                , 'laddie%Hi'                , mesh_new%pai_V)
    call checksum( laddie%Hib               , 'laddie%Hib'               , mesh_new%pai_V)
    call checksum( laddie%dHib_dx_b         , 'laddie%dHib_dx_b'         , mesh_new%pai_Tri)
    call checksum( laddie%dHib_dy_b         , 'laddie%dHib_dy_b'         , mesh_new%pai_Tri)
    call checksum( laddie%mask_icefree_land , 'laddie%mask_icefree_land' , mesh_new%pai_V)
    call checksum( laddie%mask_icefree_ocean, 'laddie%mask_icefree_ocean', mesh_new%pai_V)
    call checksum( laddie%mask_grounded_ice , 'laddie%mask_grounded_ice' , mesh_new%pai_V)
    call checksum( laddie%mask_floating_ice , 'laddie%mask_floating_ice' , mesh_new%pai_V)
    call checksum( laddie%mask_gl_fl        , 'laddie%mask_gl_fl'        , mesh_new%pai_V)
    call checksum( laddie%mask_SGD          , 'laddie%mask_SGD'          , mesh_new%pai_V)
    call checksum( laddie%Ti                , 'laddie%Ti'                , mesh_new%pai_V)
    call checksum( laddie%T_ocean           , 'laddie%T_ocean'           , mesh_new%pai_V)
    call checksum( laddie%S_ocean           , 'laddie%S_ocean'           , mesh_new%pai_V)

    ! Temperatures
    call repartition( mesh_old, mesh_new, laddie%T_amb         , laddie%wT_amb         )    ! [degC]          Temperature layer bottom
    call repartition( mesh_old, mesh_new, laddie%T_base        , laddie%wT_base        )    ! [degC]          Temperature ice shelf base
    call repartition( mesh_old, mesh_new, laddie%T_freeze      , laddie%wT_freeze      )    ! [degC]          Temperature freezing

    laddie%T_amb         ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%T_amb
    laddie%T_base        ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%T_base
    laddie%T_freeze      ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%T_freeze

    call checksum( laddie%T_amb   , 'laddie%T_amb'   , mesh_new%pai_V)
    call checksum( laddie%T_base  , 'laddie%T_base'  , mesh_new%pai_V)
    call checksum( laddie%T_freeze, 'laddie%T_freeze', mesh_new%pai_V)

    ! Salinities
    call repartition( mesh_old, mesh_new, laddie%S_amb         , laddie%wS_amb         )    ! [PSU]           Salinity layer bottom
    call repartition( mesh_old, mesh_new, laddie%S_base        , laddie%wS_base        )    ! [PSU]           Salinity ice shelf base

    laddie%S_amb         ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%S_amb
    laddie%S_base        ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%S_base

    call checksum( laddie%S_amb   , 'laddie%S_amb'   , mesh_new%pai_V)
    call checksum( laddie%S_base  , 'laddie%S_base'  , mesh_new%pai_V)

    ! Densities and buoyancies
    call repartition( mesh_old, mesh_new, laddie%rho           , laddie%wrho           )    ! [kg m^-3]       Layer density
    call repartition( mesh_old, mesh_new, laddie%rho_amb       , laddie%wrho_amb       )    ! [kg m^-3]       Ambient water density
    call repartition( mesh_old, mesh_new, laddie%drho_amb      , laddie%wdrho_amb      )    ! []              Buoyancy at layer bottom
    call repartition( mesh_old, mesh_new, laddie%Hdrho_amb     , laddie%wHdrho_amb     )    ! []              Depth-integrated buoyancy at layer bottom
    call repartition( mesh_old, mesh_new, laddie%Hdrho_amb_b   , laddie%wHdrho_amb_b   )    ! []              Depth-integrated buoyancy at layer bottom
    call repartition( mesh_old, mesh_new, laddie%drho_base     , laddie%wdrho_base     )    ! []              Buoyancy at ice base

    laddie%rho           ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%rho
    laddie%rho_amb       ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%rho_amb
    laddie%drho_amb      ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%drho_amb
    laddie%Hdrho_amb     ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%Hdrho_amb
    laddie%Hdrho_amb_b   ( mesh_new%pai_Tri%i1_nih:mesh_new%pai_Tri%i2_nih) => laddie%Hdrho_amb_b
    laddie%drho_base     ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%drho_base

    call checksum( laddie%rho        , 'laddie%rho        ', mesh_new%pai_V)
    call checksum( laddie%rho_amb    , 'laddie%rho_amb    ', mesh_new%pai_V)
    call checksum( laddie%drho_amb   , 'laddie%drho_amb   ', mesh_new%pai_V)
    call checksum( laddie%Hdrho_amb  , 'laddie%Hdrho_amb  ', mesh_new%pai_V)
    call checksum( laddie%Hdrho_amb_b, 'laddie%Hdrho_amb_b', mesh_new%pai_Tri)
    call checksum( laddie%drho_base  , 'laddie%drho_base  ', mesh_new%pai_V)

    ! Friction velocity
    call repartition( mesh_old, mesh_new, laddie%u_star        , laddie%wu_star        )    ! [m s^-1]        Friction velocity

    laddie%u_star        ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%u_star

    call checksum( laddie%u_star, 'laddie%u_star', mesh_new%pai_V)

    ! Physical parameter fields
    call repartition( mesh_old, mesh_new, laddie%gamma_T       , laddie%wgamma_T       )    ! []              Turbulent heat exchange coefficient
    call repartition( mesh_old, mesh_new, laddie%gamma_S       , laddie%wgamma_S       )    ! []              Turbulent salt exchange coefficient
    call repartition( mesh_old, mesh_new, laddie%A_h           , laddie%wA_h           )    ! [m^2 s^-1]      Horizontal laplacian viscosity
    call repartition( mesh_old, mesh_new, laddie%K_h           , laddie%wK_h           )    ! [m^2 s^-1]      Horizontal diffusivity

    laddie%gamma_T       ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%gamma_T
    laddie%gamma_S       ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%gamma_S
    laddie%A_h           ( mesh_new%pai_Tri%i1_nih:mesh_new%pai_Tri%i2_nih) => laddie%A_h
    laddie%K_h           ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%K_h

    call checksum( laddie%gamma_T, 'laddie%gamma_T', mesh_new%pai_V)
    call checksum( laddie%gamma_S, 'laddie%gamma_S', mesh_new%pai_V)
    call checksum( laddie%A_h    , 'laddie%A_h    ', mesh_new%pai_Tri)
    call checksum( laddie%K_h    , 'laddie%K_h    ', mesh_new%pai_V)

    ! Vertical rates
    call repartition( mesh_old, mesh_new, laddie%melt          , laddie%wmelt          )    ! [m s^-1]        Melting / freezing rate
    call repartition( mesh_old, mesh_new, laddie%entr          , laddie%wentr          )    ! [m s^-1]        Entrainment
    call repartition( mesh_old, mesh_new, laddie%entr_dmin     , laddie%wentr_dmin     )    ! [m s^-1]        Entrainment for D_min
    call repartition( mesh_old, mesh_new, laddie%detr          , laddie%wdetr          )    ! [m s^-1]        Detrainment
    call repartition( mesh_old, mesh_new, laddie%entr_tot      , laddie%wentr_tot      )    ! [m s^-1]        Total (net) entrainment
    call repartition( mesh_old, mesh_new, laddie%SGD           , laddie%wSGD           )    ! [m s^-1]        Subglacial discharge

    laddie%melt          ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%melt
    laddie%entr          ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%entr
    laddie%entr_dmin     ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%entr_dmin
    laddie%detr          ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%detr
    laddie%entr_tot      ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%entr_tot
    laddie%SGD           ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%SGD

    call checksum( laddie%melt     , 'laddie%melt     ', mesh_new%pai_V)
    call checksum( laddie%entr     , 'laddie%entr     ', mesh_new%pai_V)
    call checksum( laddie%entr_dmin, 'laddie%entr_dmin', mesh_new%pai_V)
    call checksum( laddie%detr     , 'laddie%detr     ', mesh_new%pai_V)
    call checksum( laddie%entr_tot , 'laddie%entr_tot ', mesh_new%pai_V)
    call checksum( laddie%SGD      , 'laddie%SGD      ', mesh_new%pai_V)

    ! Horizontal fluxes
    call repartition( mesh_old, mesh_new, laddie%divQH         , laddie%wdivQH         )    ! [m^3 s^-1]      Divergence of layer thickness
    call repartition( mesh_old, mesh_new, laddie%divQU         , laddie%wdivQU         )    ! [m^4 s^-2]      Divergence of momentum
    call repartition( mesh_old, mesh_new, laddie%divQV         , laddie%wdivQV         )    ! [m^4 s^-2]
    call repartition( mesh_old, mesh_new, laddie%divQT         , laddie%wdivQT         )    ! [degC m^3 s^-1] Divergence of heat
    call repartition( mesh_old, mesh_new, laddie%divQS         , laddie%wdivQS         )    ! [PSU m^3 s^-1]  Divergence of salt

    laddie%divQH         ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%divQH
    laddie%divQU         ( mesh_new%pai_Tri%i1_nih:mesh_new%pai_Tri%i2_nih) => laddie%divQU
    laddie%divQV         ( mesh_new%pai_Tri%i1_nih:mesh_new%pai_Tri%i2_nih) => laddie%divQV
    laddie%divQT         ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%divQT
    laddie%divQS         ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%divQS

    call checksum( laddie%divQH, 'laddie%divQH', mesh_new%pai_V)
    call checksum( laddie%divQU, 'laddie%divQU', mesh_new%pai_Tri)
    call checksum( laddie%divQV, 'laddie%divQV', mesh_new%pai_Tri)
    call checksum( laddie%divQT, 'laddie%divQT', mesh_new%pai_V)
    call checksum( laddie%divQS, 'laddie%divQS', mesh_new%pai_V)

    ! Viscosities
    call repartition( mesh_old, mesh_new, laddie%viscU         , laddie%wviscU         )    ! [m^2 s^-2]      Horizontal viscosity term
    call repartition( mesh_old, mesh_new, laddie%viscV         , laddie%wviscV         )    ! [m^2 s^-2]

    laddie%viscU         ( mesh_new%pai_Tri%i1_nih:mesh_new%pai_Tri%i2_nih) => laddie%viscU
    laddie%viscV         ( mesh_new%pai_Tri%i1_nih:mesh_new%pai_Tri%i2_nih) => laddie%viscV

    call checksum( laddie%viscU, 'laddie%viscU', mesh_new%pai_Tri)
    call checksum( laddie%viscV, 'laddie%viscV', mesh_new%pai_Tri)

    ! Diffusivities
    call repartition( mesh_old, mesh_new, laddie%diffT         , laddie%wdiffT         )    ! [degC m s^-1]   Horizontal diffusivity of heat
    call repartition( mesh_old, mesh_new, laddie%diffS         , laddie%wdiffS         )    ! [PSU m s^-1]    Horizontal diffusivity of salt

    laddie%diffT         ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%diffT
    laddie%diffS         ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%diffS

    call checksum( laddie%diffT, 'laddie%diffT', mesh_new%pai_V)
    call checksum( laddie%diffS, 'laddie%diffS', mesh_new%pai_V)

    ! RHS terms
    call repartition( mesh_old, mesh_new, laddie%ddrho_amb_dx_b, laddie%wddrho_amb_dx_b)    ! [m^-1]          Horizontal derivative of buoyancy
    call repartition( mesh_old, mesh_new, laddie%ddrho_amb_dy_b, laddie%wddrho_amb_dy_b)    ! [m^-1]
    call repartition( mesh_old, mesh_new, laddie%dH_dx_b       , laddie%wdH_dx_b       )    ! [m^-2]          Horizontal derivative of thickness
    call repartition( mesh_old, mesh_new, laddie%dH_dy_b       , laddie%wdH_dy_b       )    ! [m^-2]
    call repartition( mesh_old, mesh_new, laddie%detr_b        , laddie%wdetr_b        )    ! [m s^-1]        Detrainment on b grid

    laddie%ddrho_amb_dx_b( mesh_new%pai_Tri%i1_nih:mesh_new%pai_Tri%i2_nih) => laddie%ddrho_amb_dx_b
    laddie%ddrho_amb_dy_b( mesh_new%pai_Tri%i1_nih:mesh_new%pai_Tri%i2_nih) => laddie%ddrho_amb_dy_b
    laddie%dH_dx_b       ( mesh_new%pai_Tri%i1_nih:mesh_new%pai_Tri%i2_nih) => laddie%dH_dx_b
    laddie%dH_dy_b       ( mesh_new%pai_Tri%i1_nih:mesh_new%pai_Tri%i2_nih) => laddie%dH_dy_b
    laddie%detr_b        ( mesh_new%pai_Tri%i1_nih:mesh_new%pai_Tri%i2_nih) => laddie%detr_b

    call checksum( laddie%ddrho_amb_dx_b, 'laddie%ddrho_amb_dx_b', mesh_new%pai_Tri)
    call checksum( laddie%ddrho_amb_dy_b, 'laddie%ddrho_amb_dy_b', mesh_new%pai_Tri)
    call checksum( laddie%dH_dx_b       , 'laddie%dH_dx_b       ', mesh_new%pai_Tri)
    call checksum( laddie%dH_dy_b       , 'laddie%dH_dy_b       ', mesh_new%pai_Tri)
    call checksum( laddie%detr_b        , 'laddie%detr_b        ', mesh_new%pai_Tri)

    ! Forward-Backward Runge-Kutta 3 scheme
    call repartition( mesh_old, mesh_new, laddie%Hstar         , laddie%wHstar         )    ! [m]               Intermediate layer thickness

    laddie%Hstar         ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%Hstar

    call checksum( laddie%Hstar, 'laddie%Hstar', mesh_new%pai_V)

    ! Mapped variables
    call repartition( mesh_old, mesh_new, laddie%H_c           , laddie%wH_c           )
    call repartition( mesh_old, mesh_new, laddie%Hstar_b       , laddie%wHstar_b       )
    call repartition( mesh_old, mesh_new, laddie%Hstar_c       , laddie%wHstar_c       )

    laddie%H_c           ( mesh_new%pai_E%i1_nih  :mesh_new%pai_E%i2_nih  ) => laddie%H_c
    laddie%Hstar_b       ( mesh_new%pai_Tri%i1_nih:mesh_new%pai_Tri%i2_nih) => laddie%Hstar_b
    laddie%Hstar_c       ( mesh_new%pai_E%i1_nih  :mesh_new%pai_E%i2_nih  ) => laddie%Hstar_c

    call checksum( laddie%H_c    , 'laddie%H_c    ', mesh_new%pai_E)
    call checksum( laddie%Hstar_b, 'laddie%Hstar_b', mesh_new%pai_Tri)
    call checksum( laddie%Hstar_c, 'laddie%Hstar_c', mesh_new%pai_E)

    ! Masks
    call repartition( mesh_old, mesh_new, laddie%mask_a        , laddie%wmask_a        )    !                 Mask on a-grid
    call repartition( mesh_old, mesh_new, laddie%mask_gr_a     , laddie%wmask_gr_a     )    !                 Grounded mask on a-grid
    call repartition( mesh_old, mesh_new, laddie%mask_oc_a     , laddie%wmask_oc_a     )    !                 Icefree ocean mask on a-grid
    call repartition( mesh_old, mesh_new, laddie%mask_b        , laddie%wmask_b        )    !                 Mask on b-grid
    call repartition( mesh_old, mesh_new, laddie%mask_gl_b     , laddie%wmask_gl_b     )    !                 Grounding line mask on b-grid
    call repartition( mesh_old, mesh_new, laddie%mask_cf_b     , laddie%wmask_cf_b     )    !                 Calving front mask on b-grid
    call repartition( mesh_old, mesh_new, laddie%mask_oc_b     , laddie%wmask_oc_b     )    !                 Icefree ocean mask on b-grid

    laddie%mask_a        ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%mask_a
    laddie%mask_gr_a     ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%mask_gr_a
    laddie%mask_oc_a     ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%mask_oc_a
    laddie%mask_b        ( mesh_new%pai_Tri%i1_nih:mesh_new%pai_Tri%i2_nih) => laddie%mask_b
    laddie%mask_gl_b     ( mesh_new%pai_Tri%i1_nih:mesh_new%pai_Tri%i2_nih) => laddie%mask_gl_b
    laddie%mask_cf_b     ( mesh_new%pai_Tri%i1_nih:mesh_new%pai_Tri%i2_nih) => laddie%mask_cf_b
    laddie%mask_oc_b     ( mesh_new%pai_Tri%i1_nih:mesh_new%pai_Tri%i2_nih) => laddie%mask_oc_b

    call checksum( laddie%mask_a   , 'laddie%mask_a   ', mesh_new%pai_V)
    call checksum( laddie%mask_gr_a, 'laddie%mask_gr_a', mesh_new%pai_V)
    call checksum( laddie%mask_oc_a, 'laddie%mask_oc_a', mesh_new%pai_V)
    call checksum( laddie%mask_b   , 'laddie%mask_b   ', mesh_new%pai_Tri)
    call checksum( laddie%mask_gl_b, 'laddie%mask_gl_b', mesh_new%pai_Tri)
    call checksum( laddie%mask_cf_b, 'laddie%mask_cf_b', mesh_new%pai_Tri)
    call checksum( laddie%mask_oc_b, 'laddie%mask_oc_b', mesh_new%pai_Tri)

    ! Domains
    call repartition( mesh_old, mesh_new, laddie%domain_a      , laddie%wdomain_a      )    ! []              Floating domain on a-grid
    call repartition( mesh_old, mesh_new, laddie%domain_b      , laddie%wdomain_b      )    ! []              Floating domain on b-grid

    laddie%domain_a      ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih  ) => laddie%domain_a
    laddie%domain_b      ( mesh_new%pai_Tri%i1_nih:mesh_new%pai_Tri%i2_nih) => laddie%domain_b

    call checksum( laddie%domain_a, 'laddie%domain_a', mesh_new%pai_V)
    call checksum( laddie%domain_b, 'laddie%domain_b', mesh_new%pai_Tri)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine repartition_laddie

END MODULE laddie_main

