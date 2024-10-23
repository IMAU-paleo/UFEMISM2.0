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
  USE mpi_distributed_memory                                 , ONLY: gather_to_all_dp_1D, gather_to_all_logical_1D
  USE mesh_operators                                         , ONLY: ddx_a_b_2D, ddy_a_b_2D, map_a_b_2D, map_a_c_2D, map_b_a_2D
  USE math_utilities                                         , ONLY: check_for_NaN_dp_1D
  USE laddie_utilities                                       , ONLY: compute_ambient_TS, map_H_a_b, map_H_a_c
  USE laddie_physics                                         , ONLY: compute_buoyancy

  IMPLICIT NONE
    
CONTAINS
    
! ===== Main routines =====
! =========================

  SUBROUTINE compute_UV_npx( mesh, ice, ocean, laddie, npxref, npx, Hstar, dt, inclvisc)
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
    LOGICAL,                                INTENT(IN)    :: inclvisc


    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'compute_UV_npx'
    INTEGER                                               :: ti, ci, nfl, vj
    REAL(dp)                                              :: dHUdt, dHVdt, HU_next, HV_next, PGF_x, PGF_y, Hdrho_fl, Uabs
    LOGICAL, DIMENSION(mesh%nV)                           :: mask_a_tot
    REAL(dp), DIMENSION(mesh%nV)                          :: Hdrho_amb_tot
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2)                :: detr_b
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2)                :: Hstar_b
    REAL(dp), DIMENSION(mesh%ei1:mesh%ei2)                :: Hstar_c
 
    ! Add routine to path
    CALL init_routine( routine_name)

    CALL gather_to_all_logical_1D( laddie%mask_a, mask_a_tot)
    CALL gather_to_all_dp_1D( laddie%Hdrho_amb, Hdrho_amb_tot)

    ! Initialise ambient T and S                             
    ! TODO costly, see whether necessary to recompute with Hstar
    CALL compute_ambient_TS( mesh, ice, ocean, laddie, Hstar)

    ! Compute buoyancy
    CALL compute_buoyancy( mesh, ice, laddie, npx, Hstar)
 
    ! Bunch of mappings                                      
    CALL map_a_b_2D( mesh, laddie%detr, detr_b)
    CALL map_a_b_2D( mesh, laddie%Hdrho_amb, laddie%Hdrho_amb_b)
    CALL map_H_a_b( mesh, laddie, Hstar, Hstar_b)
    CALL map_H_a_c( mesh, laddie, Hstar, Hstar_c)

    CALL ddx_a_b_2D( mesh, laddie%drho_amb, laddie%ddrho_amb_dx_b)
    CALL ddy_a_b_2D( mesh, laddie%drho_amb, laddie%ddrho_amb_dy_b)
    CALL ddx_a_b_2D( mesh, Hstar, laddie%dH_dx_b)
    CALL ddy_a_b_2D( mesh, Hstar, laddie%dH_dy_b)

    ! Compute divergence of momentum
    CALL compute_divQUV( mesh, laddie, npx, laddie%U_c, laddie%V_c, Hstar_c, laddie%mask_b, laddie%mask_gl_b)

    ! == Integrate U and V ==
    ! =======================

    ! Loop over vertices
    DO ti = mesh%ti1, mesh%ti2
      IF (laddie%mask_b( ti)) THEN

        ! == pressure gradient force ==
        ! =============================

        IF (laddie%mask_cf_b( ti) .OR. laddie%mask_gl_b( ti)) THEN
          ! Assume dH/dx and ddrho/dx = 0
          ! Get nearest neighbour Hdrho from floating vertices
          nfl = 0
          Hdrho_fl = 0
          DO ci = 1, 3
            vj = mesh%Tri( ti, ci)
            IF (vj == 0) CYCLE
            IF (mask_a_tot( vj)) THEN
              nfl = nfl + 1
              Hdrho_fl = Hdrho_fl + Hdrho_amb_tot( vj)
            END IF
          END DO

          ! Define PGF at calving front / grounding line
          PGF_x =   grav * Hdrho_fl/nfl * laddie%dHib_dx_b( ti) &
                  - 0.5*grav * Hstar_b( ti)**2 * laddie%ddrho_amb_dx_b( ti)

          PGF_y =   grav * Hdrho_fl/nfl * laddie%dHib_dy_b( ti) &
                  - 0.5*grav * Hstar_b( ti)**2 * laddie%ddrho_amb_dy_b( ti)
        ELSE
          ! Regular full expression
          PGF_x = - grav * laddie%Hdrho_amb_b( ti) * laddie%dH_dx_b( ti) &
                  + grav * laddie%Hdrho_amb_b( ti) * laddie%dHib_dx_b( ti) &
                  - 0.5*grav * Hstar_b( ti)**2 * laddie%ddrho_amb_dx_b( ti)

          PGF_y = - grav * laddie%Hdrho_amb_b( ti) * laddie%dH_dy_b( ti) &
                  + grav * laddie%Hdrho_amb_b( ti) * laddie%dHib_dy_b( ti) &
                  - 0.5*grav * Hstar_b( ti)**2 * laddie%ddrho_amb_dy_b( ti)
        END IF

        ! == time derivatives ==
        ! ======================

        ! dHU_dt
        dHUdt = - laddie%divQU( ti) &
                + PGF_x &
                + C%uniform_laddie_coriolis_parameter * Hstar_b( ti) * npxref%V( ti) &
                - C%laddie_drag_coefficient_mom * npxref%U( ti) * (npxref%U( ti)**2 + npxref%V( ti)**2)**.5 &
                - detr_b( ti) * npxref%U( ti)

        IF (inclvisc) THEN
          dHUdt = dHUdt + laddie%viscU( ti)
        END IF

        ! dHV_dt
        dHVdt = - laddie%divQV( ti) &
                + PGF_y &
                - C%uniform_laddie_coriolis_parameter * Hstar_b( ti) * npxref%U( ti) &
                - C%laddie_drag_coefficient_mom * npxref%V( ti) * (npxref%U( ti)**2 + npxref%V( ti)**2)**.5 &
                - detr_b( ti) * npxref%V( ti)

        IF (inclvisc) THEN
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
        IF (Uabs > 0) THEN
          npx%U( ti) = npx%U( ti) * MIN(1.0_dp, C%laddie_velocity_maximum/Uabs)
          npx%V( ti) = npx%V( ti) * MIN(1.0_dp, C%laddie_velocity_maximum/Uabs)
        END IF
      END IF ! (laddie%mask_b( ti))
    END DO !ti = mesh%ti1, mesh%ti2

    CALL check_for_NaN_dp_1D( npx%U, 'U_lad')
    CALL check_for_NaN_dp_1D( npx%V, 'V_lad')

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
    REAL(dp)                                              :: D_x, D_y, D, Ah, dUabs
    REAL(dp), DIMENSION(mesh%nTri)                        :: U_tot, V_tot
    LOGICAL, DIMENSION(mesh%nTri)                         :: mask_oc_b_tot
    
    ! Add routine to path
    CALL init_routine( routine_name)        

    ! Gather
    CALL gather_to_all_dp_1D( npxref%U, U_tot)            
    CALL gather_to_all_dp_1D( npxref%V, V_tot)            
    CALL gather_to_all_logical_1D( laddie%mask_oc_b, mask_oc_b_tot)

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

          D_x = mesh%Tri( tj,1) - mesh%Tri( ti,1)
          D_y = mesh%Tri( tj,2) - mesh%Tri( ti,2)
          D   = SQRT( D_x**2 + D_y**2)

          Ah = C%laddie_viscosity ! * 0.5_dp*(SQRT(mesh%TriA( ti)) + SQRT(mesh%TriA( tj))) / 1000.0_dp

          IF (tj==0) THEN
            ! Border or corner. For now, assume no slip. If free slip: CYCLE
            laddie%viscU( ti) = laddie%viscU( ti) - npxref%U( ti) * Ah * npxref%H_b( ti) / mesh%TriA( ti)
            laddie%viscV( ti) = laddie%viscV( ti) - npxref%V( ti) * Ah * npxref%H_b( ti) / mesh%TriA( ti)
          ELSE
            ! Skip calving front - ocean connection: d/dx = d/dy = 0 
            IF (laddie%mask_oc_b( tj)) CYCLE

            dUabs = SQRT((U_tot( tj) - U_tot( ti))**2 + (V_tot( tj) - V_tot( ti))**2)
            Ah = C%laddie_viscosity * dUabs * mesh%triCw( ti, ci) / 100.0_dp
            
            ! Add viscosity flux based on dU/dx and dV/dy. 
            ! Note: for grounded neighbours, U_tot( tj) = 0, meaning this is a no slip option. Can be expanded
            laddie%viscU( ti) = laddie%viscU( ti) + (U_tot( tj) - U_tot( ti)) * Ah * npxref%H_c( ei) / mesh%TriA( ti) * mesh%TriCw( ti, ci) / D
            laddie%viscV( ti) = laddie%viscV( ti) + (V_tot( tj) - V_tot( ti)) * Ah * npxref%H_c( ei) / mesh%TriA( ti) * mesh%TriCw( ti, ci) / D
          END IF
        END DO

      END IF !(laddie%mask_b( ti)
    END DO !ti = mesh%ti1, mesh%ti2
      
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE compute_viscUV

  SUBROUTINE compute_divQUV( mesh, laddie, npxref, U_c, V_c, H_c, mask_b, mask_gl_b)
    ! Calculate the layer flux divergence matrix M_divQ using an upwind scheme
    !
    ! The vertically averaged ice flux divergence represents the net ice volume (which,
    ! assuming constant density, is proportional to the ice mass) entering each Voronoi
    ! cell per unit time. This is found by calculating the ice fluxes through each
    ! shared Voronoi cell boundary, using an upwind scheme: if ice flows from vertex vi
    ! to vertex vj, the flux is found by multiplying the velocity at their shared
    ! boundary u_c with the ice thickness at vi (and, of course, the length L_c of the
    ! shared boundary). If instead it flows from vj to vi, u_c is multiplied with the
    ! ice thickness at vj.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie
    TYPE(type_laddie_timestep),             INTENT(IN)    :: npxref
    REAL(dp), DIMENSION(mesh%ei1:mesh%ei2), INTENT(IN)    :: U_c, V_c, H_c
    LOGICAL, DIMENSION(mesh%ti1:mesh%ti2),  INTENT(IN)    :: mask_b, mask_gl_b

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'compute_divQUV'
    REAL(dp), DIMENSION(mesh%nE)                          :: U_c_tot, V_c_tot, H_c_tot
    REAL(dp), DIMENSION(mesh%nTri)                        :: U_tot, V_tot, H_b_tot
    INTEGER                                               :: ncols, ncols_loc, nrows, nrows_loc, nnz_est_proc
    INTEGER                                               :: ti, tj, ci, ei
    REAL(dp)                                              :: D_x, D_y, D_c, u_perp
    LOGICAL, DIMENSION(mesh%nTri)                         :: mask_gl_b_tot, mask_b_tot
    LOGICAL                                               :: isbound

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate vertically averaged ice velocities on the edges
    CALL gather_to_all_dp_1D( U_c, U_c_tot)
    CALL gather_to_all_dp_1D( V_c, V_c_tot)
    CALL gather_to_all_dp_1D( H_c, H_c_tot)
    CALL gather_to_all_logical_1D( mask_gl_b, mask_gl_b_tot)
    CALL gather_to_all_logical_1D( laddie%mask_b, mask_b_tot)
    CALL gather_to_all_dp_1D( npxref%U, U_tot)
    CALL gather_to_all_dp_1D( npxref%V, V_tot)
    CALL gather_to_all_dp_1D( npxref%H_b, H_b_tot)

    ! Initialise with zeros
    laddie%divQU = 0.0_dp
    laddie%divQV = 0.0_dp

    ! == Loop over triangles ==
    ! =========================

    DO ti = mesh%ti1, mesh%ti2

      IF (mask_b( ti)) THEN

        ! Loop over all connections of triangle ti
        DO ci = 1, 3

          tj = mesh%TriC( ti, ci)

          ! Skip if no connecting triangle on this side
          IF (tj == 0) CYCLE

          ! Skip connection if neighbour is grounded. No flux across grounding line
          ! Can be made more flexible when accounting for partial cells (PMP instead of FCMP)
          ! IF (mask_gl_b_tot( tj)) CYCLE ! Skip grounded
          IF (mask_b_tot( tj)) CYCLE ! Skip grounded and ocean

          ! Get edge index
          ei = mesh%TriE( ti, ci)

          ! The triangle-triangle vector from ti to tj
          D_x = mesh%Tri( tj,1) - mesh%Tri( ti,1)
          D_y = mesh%Tri( tj,2) - mesh%Tri( ti,2)
          D_c = SQRT( D_x**2 + D_y**2)

          ! Calculate vertically averaged ice velocity component perpendicular to this edge
          u_perp = U_c_tot( ei) * D_y/D_c + V_c_tot( ei) * D_x/D_c

          ! Calculate momentum divergence
          ! =============================
          ! Centered:
          laddie%divQU( ti) = laddie%divQU( ti) + mesh%TriCw( ti, ci) * u_perp * U_c_tot( ei) * H_c_tot( ei) / mesh%TriA( ti)
          laddie%divQV( ti) = laddie%divQV( ti) + mesh%TriCw( ti, ci) * u_perp * V_c_tot( ei) * H_c_tot( ei) / mesh%TriA( ti)

          ! Upwind:
          ! u_perp > 0: flow is exiting this vertex into vertex vj
          !IF (u_perp > 0) THEN
          !  laddie%divQU( ti) = laddie%divQU( ti) + mesh%TriCw( ti, ci) * u_perp * H_b_tot( ti) * U_tot( ti) / mesh%TriA( ti)
          !  laddie%divQV( ti) = laddie%divQV( ti) + mesh%TriCw( ti, ci) * u_perp * H_b_tot( ti) * V_tot( ti) / mesh%TriA( ti)
          ! u_perp < 0: flow is entering this vertex from vertex vj
          !ELSE
          !  laddie%divQU( ti) = laddie%divQU( ti) + mesh%TriCw( ti, ci) * u_perp * H_b_tot( tj) * U_tot( tj) / mesh%TriA( ti)
          !  laddie%divQV( ti) = laddie%divQV( ti) + mesh%TriCw( ti, ci) * u_perp * H_b_tot( tj) * V_tot( tj) / mesh%TriA( ti)
          !END IF
        END DO ! DO ci = 1, 3

      END IF ! (mask_b( ti))

    END DO ! DO ti = mesh%ti1, mesh%ti2

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE compute_divQUV

END MODULE laddie_velocity

