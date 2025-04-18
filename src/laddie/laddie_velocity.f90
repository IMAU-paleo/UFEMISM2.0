MODULE laddie_velocity

  ! Velocity routines for the laddie model

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
  USE reallocate_mod                                         , ONLY: reallocate_bounds
  USE mpi_distributed_memory                                 , ONLY: gather_to_all
  USE mesh_disc_apply_operators                              , ONLY: ddx_a_b_2D, ddy_a_b_2D, map_a_b_2D, map_b_a_2D
  USE laddie_utilities                                       , ONLY: compute_ambient_TS, map_H_a_b, map_H_a_c
  USE laddie_physics                                         , ONLY: compute_buoyancy
  use CSR_matrix_vector_multiplication, only: multiply_CSR_matrix_with_vector_1D
  use mesh_halo_exchange, only: exchange_halos
  use mesh_integrate_over_domain, only: calc_and_print_min_mean_max

  IMPLICIT NONE

CONTAINS

! ===== Main routines =====
! =========================

  SUBROUTINE compute_UV_npx( mesh, ice, ocean, laddie, npxref, npx, Hstar, dt, include_viscosity_terms)
    ! Integrate U and V by one time step

    ! In- and output variables

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_ocean_model),                 INTENT(IN)    :: ocean
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie
    TYPE(type_laddie_timestep),             INTENT(IN)    :: npxref
    TYPE(type_laddie_timestep),             INTENT(INOUT) :: npx
    REAL(dp),                               INTENT(IN)    :: dt
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)    :: Hstar
    LOGICAL,                                INTENT(IN)    :: include_viscosity_terms

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'compute_UV_npx'
    INTEGER                                               :: ti, ci, nfl, vj
    REAL(dp)                                              :: dHUdt, dHVdt, HU_next, HV_next, PGF_x, PGF_y, Hdrho_fl, Uabs

    ! Add routine to path
    CALL init_routine( routine_name)

    call exchange_halos( mesh, laddie%mask_a)

    ! Initialise ambient T and S
    ! TODO costly, see whether necessary to recompute with Hstar
    CALL compute_ambient_TS( mesh, ice, ocean, laddie, Hstar)

    ! Compute buoyancy
    CALL compute_buoyancy( mesh, ice, laddie, npx, Hstar)

    ! Bunch of mappings - TODO: move b-grid versions to Laddie type, make them hybrid; faster
    call exchange_halos( mesh, laddie%detr)
    call exchange_halos( mesh, laddie%Hdrho_amb)
    ! call exchange_halos( mesh, Hstar) ! Already done in integrate_fbrk3
    CALL map_a_b_2D( mesh, laddie%detr, laddie%detr_b, d_a_is_hybrid = .true., d_b_is_hybrid = .true.)
    CALL map_H_a_b( mesh, laddie, laddie%Hdrho_amb, laddie%Hdrho_amb_b)
    CALL map_H_a_b( mesh, laddie, Hstar, laddie%Hstar_b)
    CALL map_H_a_c( mesh, laddie, Hstar, laddie%Hstar_c)
    call exchange_halos( mesh, laddie%Hstar_b)
    call exchange_halos( mesh, npxref%H_b)

    ! Bunch of derivatives
    call exchange_halos( mesh, laddie%drho_amb)
    CALL ddx_a_b_2D( mesh, laddie%drho_amb, laddie%ddrho_amb_dx_b, d_a_is_hybrid = .true., ddx_b_is_hybrid = .true.)
    CALL ddy_a_b_2D( mesh, laddie%drho_amb, laddie%ddrho_amb_dy_b, d_a_is_hybrid = .true., ddy_b_is_hybrid = .true.)
    CALL ddx_a_b_2D( mesh, Hstar, laddie%dH_dx_b, d_a_is_hybrid = .true., ddx_b_is_hybrid = .true.)
    CALL ddy_a_b_2D( mesh, Hstar, laddie%dH_dy_b, d_a_is_hybrid = .true., ddy_b_is_hybrid = .true.)

    ! Compute divergence of momentum
    SELECT CASE(C%choice_laddie_momentum_advection)
      CASE DEFAULT
        CALL crash('unknown choice_laddie_momentum_advection "' // TRIM( C%choice_laddie_momentum_advection) // '"')
      CASE ('none')
        laddie%divQU( mesh%ti1:mesh%ti2) = 0.0_dp
        laddie%divQV( mesh%ti1:mesh%ti2) = 0.0_dp
      CASE ('upstream')
        ! TODO figure out which of the below lines is best
        !CALL compute_divQUV_upstream( mesh, laddie, npx, laddie%Hstar_b)
        CALL compute_divQUV_upstream( mesh, laddie, npx, npxref%H_b)
    END SELECT

    ! == Integrate U and V ==
    ! =======================

    ! Loop over vertices
    DO ti = mesh%ti1, mesh%ti2
      IF (laddie%mask_b( ti)) THEN

        ! == pressure gradient force ==
        ! =============================

        IF (laddie%mask_cf_b( ti) .OR. laddie%mask_gl_b( ti)) THEN
          ! Assume dH/dx and ddrho/dx = 0

          ! Define PGF at calving front / grounding line
          PGF_x = grav * laddie%Hdrho_amb_b( ti) * ice%dHib_dx_b( ti) &
                  - 0.5*grav * laddie%Hstar_b( ti)**2 * laddie%ddrho_amb_dx_b( ti)

          PGF_y = grav * laddie%Hdrho_amb_b( ti) * ice%dHib_dy_b( ti) &
                  - 0.5*grav * laddie%Hstar_b( ti)**2 * laddie%ddrho_amb_dy_b( ti)
        ELSE
          ! Regular full expression
          PGF_x = - grav * laddie%Hdrho_amb_b( ti) * laddie%dH_dx_b( ti) &
                  + grav * laddie%Hdrho_amb_b( ti) * ice%dHib_dx_b( ti) &
                  - 0.5*grav * laddie%Hstar_b( ti)**2 * laddie%ddrho_amb_dx_b( ti)

          PGF_y = - grav * laddie%Hdrho_amb_b( ti) * laddie%dH_dy_b( ti) &
                  + grav * laddie%Hdrho_amb_b( ti) * ice%dHib_dy_b( ti) &
                  - 0.5*grav * laddie%Hstar_b( ti)**2 * laddie%ddrho_amb_dy_b( ti)
        END IF

        ! == time derivatives ==
        ! ======================

        ! dHU_dt
        dHUdt = - laddie%divQU( ti) &
                + PGF_x &
                + C%uniform_laddie_coriolis_parameter * laddie%Hstar_b( ti) * npxref%V( ti) &
                - C%laddie_drag_coefficient_mom * npxref%U( ti) * (npxref%U( ti)**2 + npxref%V( ti)**2)**.5 &
                - laddie%detr_b( ti) * npxref%U( ti)

        IF (include_viscosity_terms) THEN
          dHUdt = dHUdt + laddie%viscU( ti)
        END IF

        ! dHV_dt
        dHVdt = - laddie%divQV( ti) &
                + PGF_y &
                - C%uniform_laddie_coriolis_parameter * laddie%Hstar_b( ti) * npxref%U( ti) &
                - C%laddie_drag_coefficient_mom * npxref%V( ti) * (npxref%U( ti)**2 + npxref%V( ti)**2)**.5 &
                - laddie%detr_b( ti) * npxref%V( ti)

        IF (include_viscosity_terms) THEN
          dHVdt = dHVdt + laddie%viscV( ti)
        END IF

        ! == next time step ==
        ! ====================

        ! HU_n = HU_n + dHU_dt * dt
        HU_next = laddie%now%U( ti)*laddie%now%H_b( ti) + dHUdt * dt
        HV_next = laddie%now%V( ti)*laddie%now%H_b( ti) + dHVdt * dt

        ! U_n = HU_n / H_n
        npx%U( ti) = HU_next / npx%H_b( ti)
        npx%V( ti) = HV_next / npx%H_b( ti)

      END IF ! (laddie%mask_b( ti))
    END DO !ti = mesh%ti1, mesh%ti2

    ! Cutoff velocities to ensure Uabs <= Uabs_max
    DO ti = mesh%ti1, mesh%ti2
      IF (laddie%mask_b( ti)) THEN
        ! Get absolute velocity
        Uabs = (npx%U( ti)**2 + npx%V( ti)**2)**.5

        ! Scale U and V
        IF (Uabs == 0) CYCLE ! Prevent division by zero
        npx%U( ti) = npx%U( ti) * MIN(1.0_dp, C%laddie_velocity_maximum/Uabs)
        npx%V( ti) = npx%V( ti) * MIN(1.0_dp, C%laddie_velocity_maximum/Uabs)

      END IF ! (laddie%mask_b( ti))
    END DO !ti = mesh%ti1, mesh%ti2

    ! Map velocities to a and c grid
    call exchange_halos( mesh, npx%U)
    call exchange_halos( mesh, npx%V)
    CALL map_UV_b_c( mesh, laddie, npx%U, npx%V, npx%U_c, npx%V_c)
    CALL map_b_a_2D( mesh, npx%U, npx%U_a, d_b_is_hybrid = .true., d_a_is_hybrid = .true.)
    CALL map_b_a_2D( mesh, npx%V, npx%V_a, d_b_is_hybrid = .true., d_a_is_hybrid = .true.)

    call calc_and_print_min_mean_max( mesh, npx%U, 'npx%U')
    call calc_and_print_min_mean_max( mesh, npx%V, 'npx%V')
    call calc_and_print_min_mean_max( mesh, npx%U_c, 'npx%U_c')
    call calc_and_print_min_mean_max( mesh, npx%V_c, 'npx%V_c')
    call calc_and_print_min_mean_max( mesh, npx%U_a, 'npx%U_a')
    call calc_and_print_min_mean_max( mesh, npx%V_a, 'npx%V_a')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE compute_UV_npx

  SUBROUTINE compute_viscUV( mesh, ice, laddie, npxref)
    ! Compute horizontal viscosity of momentum

    ! In- and output variables

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie
    TYPE(type_laddie_timestep),             INTENT(IN)    :: npxref

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'compute_viscUV'
    INTEGER                                               :: ci, ti, tj, ei
    REAL(dp)                                              :: Ah, dUabs

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Exchange halos
    call exchange_halos( mesh, npxref%U)
    call exchange_halos( mesh, npxref%V)
    call exchange_halos( mesh, laddie%mask_oc_b)
    call exchange_halos( mesh, npxref%H_c)

    ! Loop over triangles
    DO ti = mesh%ti1, mesh%ti2
      IF (laddie%mask_b( ti)) THEN

        ! Initialise at 0
        laddie%viscU( ti) = 0.0_dp
        laddie%viscV( ti) = 0.0_dp

        ! Loop over connected triangles
        DO ci = 1, 3
          tj = mesh%TriC( ti, ci)
          ei = mesh%TriE( ti, ci)

          Ah = C%laddie_viscosity ! * 0.5_dp*(SQRT(mesh%TriA( ti)) + SQRT(mesh%TriA( tj))) / 1000.0_dp

          IF (tj==0) THEN
            ! Border or corner. For now, assume no slip. If free slip: CYCLE
            laddie%viscU( ti) = laddie%viscU( ti) - npxref%U( ti) * Ah * npxref%H_b( ti) / mesh%TriA( ti)
            laddie%viscV( ti) = laddie%viscV( ti) - npxref%V( ti) * Ah * npxref%H_b( ti) / mesh%TriA( ti)
          ELSE
            ! Skip calving front - ocean connection: d/dx = d/dy = 0
            IF (laddie%mask_oc_b( tj)) CYCLE

            dUabs = SQRT((npxref%U( tj) - npxref%U( ti))**2 + (npxref%V( tj) - npxref%V( ti))**2)
            Ah = C%laddie_viscosity * dUabs * mesh%triCw( ti, ci) / 100.0_dp

            ! Add viscosity flux based on dU/dx and dV/dy.
            ! Note: for grounded neighbours, npxref%U( tj) = 0, meaning this is a no slip option. Can be expanded
            laddie%viscU( ti) = laddie%viscU( ti) + (npxref%U( tj) - npxref%U( ti)) * Ah * &
              npxref%H_c( ei) / mesh%TriA( ti) * mesh%TriCw( ti, ci) / mesh%TriD( ti, ci)
            laddie%viscV( ti) = laddie%viscV( ti) + (npxref%V( tj) - npxref%V( ti)) * Ah * &
              npxref%H_c( ei) / mesh%TriA( ti) * mesh%TriCw( ti, ci) / mesh%TriD( ti, ci)
          END IF
        END DO

      END IF !(laddie%mask_b( ti)
    END DO !ti = mesh%ti1, mesh%ti2
    call calc_and_print_min_mean_max( mesh, laddie%viscU, 'laddie%viscU')
    call calc_and_print_min_mean_max( mesh, laddie%viscV, 'laddie%viscV')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE compute_viscUV

  subroutine compute_divQUV_upstream( mesh, laddie, npxref, Hstar_b)
    ! Upstream scheme

    ! In/output variables:
    type(type_mesh),                        intent(in)    :: mesh
    type(type_laddie_model),                intent(inout) :: laddie
    type(type_laddie_timestep),             intent(in)    :: npxref
    real(dp), dimension(mesh%pai_Tri%i1_nih:mesh%pai_Tri%i2_nih), intent(in)    :: Hstar_b

    ! Local variables:
    character(len=256), parameter                         :: routine_name = 'compute_divQUV_upstream'
    integer                                               :: ti, tj, ci, ei
    real(dp)                                              :: u_perp

    ! Add routine to path
    call init_routine( routine_name)

    call exchange_halos( mesh, laddie%mask_gl_b)
    call exchange_halos( mesh, laddie%mask_cf_b)
    call exchange_halos( mesh, laddie%mask_b)
    call exchange_halos( mesh, npxref%U)
    call exchange_halos( mesh, npxref%V)
    call exchange_halos( mesh, npxref%U_c)
    call exchange_halos( mesh, npxref%V_c)
    ! call exchange_halos( mesh, Hstar_b, H_b_tot) ! Already done in compute_UV_npx

    ! Initialise with zeros
    laddie%divQU( mesh%ti1:mesh%ti2) = 0.0_dp
    laddie%divQV( mesh%ti1:mesh%ti2) = 0.0_dp

    ! == Loop over triangles ==
    ! =========================

    do ti = mesh%ti1, mesh%ti2

      if (laddie%mask_b( ti)) then

        ! Loop over all connections of triangle ti
        do ci = 1, 3

          tj = mesh%TriC( ti, ci)
          ei = mesh%TriE( ti, ci)

          ! Skip if no connecting triangle on this side
          if (tj == 0) cycle

          ! Skip connection if neighbour is grounded. No flux across grounding line
          if (laddie%mask_gl_b( tj)) cycle

          ! Calculate vertically averaged water velocity component perpendicular to this edge
          u_perp = npxref%U_c( ei) * mesh%TriD_x( ti, ci) / mesh%TriD( ti, ci) &
                 + npxref%V_c( ei) * mesh%TriD_y( ti, ci) / mesh%TriD( ti, ci)

          ! Calculate upstream momentum divergence
          ! =============================
          ! u_perp > 0: flow is exiting this triangle into triangle tj
          if (u_perp > 0) then
            laddie%divQU( ti) = laddie%divQU( ti) + mesh%TriCw( ti, ci) * Hstar_b( ti) * npxref%U( ti)* u_perp / mesh%TriA( ti)
          ! u_perp < 0: flow is entering this triangle into triangle tj
          else
            laddie%divQU( ti) = laddie%divQU( ti) + mesh%TriCw( ti, ci) * Hstar_b( tj) * npxref%U( tj)* u_perp / mesh%TriA( ti)
          end if

          ! V momentum
          if (u_perp > 0) then
            laddie%divQV( ti) = laddie%divQV( ti) + mesh%TriCw( ti, ci) * Hstar_b( ti) * npxref%V( ti)* u_perp / mesh%TriA( ti)
          else
            laddie%divQV( ti) = laddie%divQV( ti) + mesh%TriCw( ti, ci) * Hstar_b( tj) * npxref%V( tj)* u_perp / mesh%TriA( ti)
          end if

        end do ! do ci = 1, 3

      end if ! (laddie%mask_b( ti))

    end do ! do ti = mesh%ti1, mesh%ti2
    call calc_and_print_min_mean_max( mesh, laddie%divQU, 'laddie%divQU')
    call calc_and_print_min_mean_max( mesh, laddie%divQV, 'laddie%divQV')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine compute_divQUV_upstream

  subroutine map_UV_b_c( mesh, laddie, U, V, U_c, V_c)
    ! Calculate velocities on the c-grid for solving the layer thickness equation
    !
    ! Uses a different scheme then the standard mapping operator, as that one is too diffusive

    ! In/output variables:
    type(type_mesh),                        intent(in)    :: mesh
    type(type_laddie_model),                intent(in)    :: laddie
    real(dp), dimension(:),                 intent(in)    :: U
    real(dp), dimension(:),                 intent(in)    :: V
    real(dp), dimension(:),                 intent(out)   :: U_c
    real(dp), dimension(:),                 intent(out)   :: V_c

    ! Local variables:
    character(len=256), parameter                         :: routine_name = 'map_UV_b_c'

    ! Add routine to path
    call init_routine( routine_name)

    call multiply_CSR_matrix_with_vector_1D( laddie%M_map_UV_b_c, &
      mesh%pai_Tri, U, mesh%pai_E, U_c)
    call multiply_CSR_matrix_with_vector_1D( laddie%M_map_UV_b_c, &
      mesh%pai_Tri, V, mesh%pai_E, V_c)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_UV_b_c

END MODULE laddie_velocity


