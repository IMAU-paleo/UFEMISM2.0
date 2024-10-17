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
  USE mesh_operators                                         , ONLY: map_a_b_2D
  USE math_utilities                                         , ONLY: check_for_NaN_dp_1D

  IMPLICIT NONE
    
CONTAINS
    
! ===== Main routines =====
! =========================

  SUBROUTINE compute_UV_npx( mesh, ice, laddie, npxref, npx, Hstar_b, dt, inclvisc)
    ! Integrate U and V by one time step

    ! In- and output variables

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_laddie_model),                INTENT(IN)    :: laddie
    TYPE(type_laddie_timestep),             INTENT(IN)    :: npxref
    TYPE(type_laddie_timestep),             INTENT(INOUT) :: npx
    REAL(dp),                               INTENT(IN)    :: dt
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2), INTENT(IN)    :: Hstar_b
    LOGICAL,                                INTENT(IN)    :: inclvisc

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'compute_UV_npx'
    INTEGER                                               :: ti, ci, nfl, vj
    REAL(dp)                                              :: dHUdt, dHVdt, HU_next, HV_next, PGF_x, PGF_y, Hdrho_fl, Uabs
    LOGICAL, DIMENSION(mesh%nV)                           :: mask_a_tot
    REAL(dp), DIMENSION(mesh%nV)                          :: Hdrho_amb_tot
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2)                :: detr_b
 
    ! Add routine to path
    CALL init_routine( routine_name)

    CALL gather_to_all_logical_1D( laddie%mask_a, mask_a_tot)
    CALL gather_to_all_dp_1D( laddie%Hdrho_amb, Hdrho_amb_tot)

    ! Map detrainment to b grid
    CALL map_a_b_2D( mesh, laddie%detr, detr_b)

    ! == Integrate U and V ==
    ! =======================

    ! Loop over vertices
    DO ti = mesh%ti1, mesh%ti2
      IF (laddie%mask_b( ti)) THEN

        ! == pressure gradient force ==
        ! =============================

        IF (laddie%mask_cf_b( ti)) THEN
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

          ! Define PGF at calving front
          PGF_x =   grav * Hdrho_fl/nfl * laddie%dHib_dx_b( ti)
          PGF_y =   grav * Hdrho_fl/nfl * laddie%dHib_dy_b( ti)
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
                - C%laddie_drag_coefficient * npxref%U( ti) * (npxref%U( ti)**2 + npxref%V( ti)**2)**.5 &
                - detr_b( ti) * npxref%U( ti)

        IF (inclvisc) THEN
          dHUdt = dHUdt + laddie%viscU( ti)
        END IF

        ! dHV_dt
        dHVdt = - laddie%divQV( ti) &
                + PGF_y &
                - C%uniform_laddie_coriolis_parameter * Hstar_b( ti) * npxref%U( ti) &
                - C%laddie_drag_coefficient * npxref%V( ti) * (npxref%U( ti)**2 + npxref%V( ti)**2)**.5 &
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

  SUBROUTINE compute_viscUV( mesh, ice, laddie)       
    ! Compute horizontal viscosity of momentum          
  
    ! In- and output variables                               

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'compute_viscUV'
    INTEGER                                               :: i, j, vi, vj, ti, tj, nf1, nf2
    REAL(dp), DIMENSION(mesh%nTri)                        :: U_tot, V_tot
    LOGICAL, DIMENSION(mesh%nV)                           :: mask_a_tot
    LOGICAL, DIMENSION(mesh%nTri)                         :: mask_gl_b_tot
    
    ! Add routine to path
    CALL init_routine( routine_name)        

    ! Gather
    CALL gather_to_all_dp_1D( laddie%now%U, U_tot)            
    CALL gather_to_all_dp_1D( laddie%now%V, V_tot)            
    CALL gather_to_all_logical_1D( laddie%mask_a, mask_a_tot)
    CALL gather_to_all_logical_1D( laddie%mask_gl_b, mask_gl_b_tot)

    ! Loop over triangles                                  
    DO ti = mesh%ti1, mesh%ti2
      IF (laddie%mask_b( ti)) THEN
        ! Get viscosity parameter
        laddie%A_h( ti) = C%laddie_viscosity
        ! TODO add scalable options
    
        ! Initialise at 0
        laddie%viscU( ti) = 0.0_dp
        laddie%viscV( ti) = 0.0_dp
      
        ! Loop over connected triangles
        DO i = 1, 3
          tj = mesh%TriC( ti, i)

          IF (tj==0) THEN
            ! Border or corner. For now, assume no slip. If free slip: CYCLE
            laddie%viscU( ti) = laddie%viscU( ti) - laddie%now%U( ti) * laddie%A_h( ti) * laddie%now%H_b( ti) / mesh%TriA( ti)
            laddie%viscV( ti) = laddie%viscV( ti) - laddie%now%V( ti) * laddie%A_h( ti) * laddie%now%H_b( ti) / mesh%TriA( ti)
          ELSE
            ! Skip calving front - ocean connection: d/dx = d/dy = 0 
            IF (laddie%mask_oc_b( tj)) CYCLE
            
            ! Add viscosity flux based on dU/dx and dV/dy. 
            ! Note: for grounded neighbours, U_tot( tj) = 0, meaning this is a no slip option. Can be expanded
            laddie%viscU( ti) = laddie%viscU( ti) + (U_tot( tj)-laddie%now%U( ti)) * laddie%A_h( ti) * laddie%now%H_b( ti) / mesh%TriA( ti)
            laddie%viscV( ti) = laddie%viscV( ti) + (V_tot( tj)-laddie%now%V( ti)) * laddie%A_h( ti) * laddie%now%H_b( ti) / mesh%TriA( ti)
          END IF
        END DO

      END IF !(laddie%mask_b( ti)
    END DO !ti = mesh%ti1, mesh%ti2
      
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE compute_viscUV

  SUBROUTINE compute_divQUV_centered( mesh, laddie, U_c, V_c, H_c, mask_b, mask_gl_b)
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
    REAL(dp), DIMENSION(mesh%ei1:mesh%ei2), INTENT(IN)    :: U_c, V_c, H_c
    LOGICAL, DIMENSION(mesh%ti1:mesh%ti2),  INTENT(IN)    :: mask_b, mask_gl_b

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'compute_divQUV_centered'
    REAL(dp), DIMENSION(mesh%nE)                          :: U_c_tot, V_c_tot, H_c_tot
    INTEGER                                               :: ncols, ncols_loc, nrows, nrows_loc, nnz_est_proc
    INTEGER                                               :: ti, ci, ei, tj, vi, vi1, vi2, i, j, e, k
    REAL(dp)                                              :: A_i, L_c, H_e
    REAL(dp)                                              :: L_x, L_y, u_perp
    LOGICAL, DIMENSION(mesh%nTri)                         :: mask_b_tot
    LOGICAL, DIMENSION(mesh%nTri)                         :: mask_gl_b_tot
    LOGICAL                                               :: isbound

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate vertically averaged ice velocities on the edges
    CALL gather_to_all_dp_1D( U_c, U_c_tot)
    CALL gather_to_all_dp_1D( V_c, V_c_tot)
    CALL gather_to_all_dp_1D( H_c, H_c_tot)
    CALL gather_to_all_logical_1D( mask_b, mask_b_tot)
    CALL gather_to_all_logical_1D( mask_gl_b, mask_gl_b_tot)

    ! Initialise with zeros
    laddie%divQU = 0.0_dp
    laddie%divQV = 0.0_dp

    ! == Loop over triangles ==
    ! =========================

    DO ti = mesh%ti1, mesh%ti2

      IF (mask_b( ti)) THEN

        ! Loop over all connections of triangle ti
        DO ci = 1, 3

          ! TODO move definition of TriE to mesh routine
          ! == Get ei and tj ==
          ei = 0
          tj = 0
          ! Get neighbouring vertex
          vi1 = mesh%Tri( ti, ci)

          ! Get the other neighbouring vertex
          IF (i < 0) THEN
            vi2 = mesh%Tri( ti, ci+1)
          ELSE
            vi2 = mesh%Tri( ti, 1)
          END IF

          ! Loop over edges connected to first vertex
          DO j = 1,mesh%nC_mem
            e = mesh%VE( vi1, j)
            IF ( e == 0) CYCLE
            DO k = 1,2
              IF (mesh%EV( e, k) == vi2) THEN
                ei = e
              END IF
            END DO
          END DO ! j = 1, mesh%nC_mem 
          
          IF ( ei == 0) CYCLE
          ! TODO make sure ei is not 0. Should not be possible
          ! TODO what happens if this is a border?          

          ! Get triangle bordering this shared edge
          IF (mesh%ETri( ei, 1) == ti) THEN
            tj = mesh%Etri( ei, 2)
          ELSEIF (mesh%ETri( ei, 2) == ti) THEN
            tj = mesh%Etri( ei, 1)
          ELSE
            CYCLE
          END IF
          IF (tj == 0) CYCLE
          ! =================

          ! Skip connection if neighbour is grounded. No flux across grounding line
          ! Can be made more flexible when accounting for partial cells (PMP instead of FCMP)
          IF (mask_gl_b_tot( tj)) CYCLE

          ! The Voronoi cell of triangle ti has area A_i
          A_i = mesh%TriA( ti)

          ! Get thickness at edge

          H_e = H_c_tot( ei)
          ! Overwrite with triangle value if at boundary
          isbound = .false.
          DO i = 1, 4
            vi = mesh%EV( ei, i)
            IF (.NOT. laddie%mask_a( vi)) THEN
              isbound = .true.
            END IF
          END DO

          ! Overwrite with triangle (b grid value) if at boundary
          IF (isbound) THEN
            H_e = laddie%now%H_b( ti)
          END IF

          ! The shared edge length of triangles ti and tj has length L_c 
          L_x = mesh%V( vi1,1) - mesh%V( vi2,1)
          L_y = mesh%V( vi1,2) - mesh%V( vi2,2)
          L_c = SQRT( L_x**2 + L_y**2)

          ! Skip if for some reason, L_c = 0
          IF (L_c == 0) CYCLE

          ! Calculate vertically averaged ice velocity component perpendicular to this edge
          u_perp = - U_c_tot( ei) * L_y/L_c + V_c_tot( ei) * L_x/L_c

          ! Calculate momentum divergence
          ! =============================
          ! Upwind:
          ! laddie%divQU( ti) = laddie%divQU( ti) + L_c * u_perp * laddie%U( ti) * laddie%H_b( ti) / A_i
          ! laddie%divQV( ti) = laddie%divQV( ti) + L_c * u_perp * laddie%V( ti) * laddie%H_b( ti) / A_i
          ! Centered:
          laddie%divQU( ti) = laddie%divQU( ti) + L_c * u_perp * U_c_tot( ei) * H_e / A_i
          laddie%divQV( ti) = laddie%divQV( ti) + L_c * u_perp * V_c_tot( ei) * H_e / A_i

        END DO ! DO ci = 1, 3

      END IF ! (mask_b( ti))

    END DO ! DO ti = mesh%ti1, mesh%ti2

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE compute_divQUV_centered

END MODULE laddie_velocity

