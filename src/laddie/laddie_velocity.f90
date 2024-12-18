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
  USE mesh_operators                                         , ONLY: ddx_a_b_2D, ddy_a_b_2D, map_a_b_2D
  USE math_utilities                                         , ONLY: check_for_NaN_dp_1D
  USE laddie_utilities                                       , ONLY: compute_ambient_TS, map_H_a_b, map_H_a_c
  USE laddie_physics                                         , ONLY: compute_buoyancy
  USE mesh_operators                                         , ONLY: map_b_c_2D, map_b_a_2D

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

    ! Initialise ambient T and S                             
    ! TODO costly, see whether necessary to recompute with Hstar
    CALL compute_ambient_TS( mesh, ice, ocean, laddie, Hstar)

    ! Compute buoyancy
    CALL compute_buoyancy( mesh, ice, laddie, npx, Hstar)
 
    ! Bunch of mappings                                      
    CALL map_a_b_2D( mesh, laddie%detr, detr_b)
    CALL map_H_a_b( mesh, laddie, laddie%Hdrho_amb, laddie%Hdrho_amb_b)
    CALL map_H_a_b( mesh, laddie, Hstar, Hstar_b)
    CALL map_H_a_c( mesh, laddie, Hstar, Hstar_c)

    ! Bunch of derivatives
    CALL ddx_a_b_2D( mesh, laddie%drho_amb, laddie%ddrho_amb_dx_b)
    CALL ddy_a_b_2D( mesh, laddie%drho_amb, laddie%ddrho_amb_dy_b)
    CALL ddx_a_b_2D( mesh, Hstar, laddie%dH_dx_b)
    CALL ddy_a_b_2D( mesh, Hstar, laddie%dH_dy_b)

    ! Compute divergence of momentum
    SELECT CASE(C%choice_laddie_momentum_advection)
      CASE DEFAULT
        CALL crash('unknown choice_laddie_momentum_advection "' // TRIM( C%choice_laddie_momentum_advection) // '"')
      CASE ('none')
        laddie%divQU = 0.0_dp
        laddie%divQV = 0.0_dp
      CASE ('upstream')
        !CALL compute_divQUV_upstream( mesh, laddie, npx, Hstar_b)
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
                  - 0.5*grav * Hstar_b( ti)**2 * laddie%ddrho_amb_dx_b( ti)

          PGF_y = grav * laddie%Hdrho_amb_b( ti) * ice%dHib_dy_b( ti) &
                  - 0.5*grav * Hstar_b( ti)**2 * laddie%ddrho_amb_dy_b( ti)
        ELSE
          ! Regular full expression
          PGF_x = - grav * laddie%Hdrho_amb_b( ti) * laddie%dH_dx_b( ti) &
                  + grav * laddie%Hdrho_amb_b( ti) * ice%dHib_dx_b( ti) &
                  - 0.5*grav * Hstar_b( ti)**2 * laddie%ddrho_amb_dx_b( ti)

          PGF_y = - grav * laddie%Hdrho_amb_b( ti) * laddie%dH_dy_b( ti) &
                  + grav * laddie%Hdrho_amb_b( ti) * ice%dHib_dy_b( ti) &
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

    ! Map velocities to a and c grid
    CALL map_laddie_velocities_from_b_to_c_2D( mesh, npx%U, npx%V, npx%U_c, npx%V_c)
    CALL map_b_a_2D( mesh, npx%U, npx%U_a)
    CALL map_b_a_2D( mesh, npx%V, npx%V_a)

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

          D_x = mesh%Tricc( tj,1) - mesh%Tricc( ti,1)
          D_y = mesh%Tricc( tj,2) - mesh%Tricc( ti,2)
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

  SUBROUTINE compute_divQUV_upstream( mesh, laddie, npxref, Hstar_b)
    ! Upstream scheme

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie
    TYPE(type_laddie_timestep),             INTENT(IN)    :: npxref
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2), INTENT(IN)    :: Hstar_b

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'compute_divQUV_upstream'
    REAL(dp), DIMENSION(mesh%nTri)                        :: U_tot, V_tot, H_b_tot
    REAL(dp), DIMENSION(mesh%nE)                          :: U_c_tot, V_c_tot
    INTEGER                                               :: ti, tj, ci, ei
    REAL(dp)                                              :: D_x, D_y, D_c, u_perp_x, u_perp_y
    LOGICAL, DIMENSION(mesh%nTri)                         :: mask_gl_b_tot, mask_cf_b_tot, mask_b_tot

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate vertically averaged ice velocities on the edges
    CALL gather_to_all_logical_1D( laddie%mask_gl_b, mask_gl_b_tot)
    CALL gather_to_all_logical_1D( laddie%mask_cf_b, mask_cf_b_tot)
    CALL gather_to_all_logical_1D( laddie%mask_b, mask_b_tot)
    CALL gather_to_all_dp_1D( npxref%U, U_tot)
    CALL gather_to_all_dp_1D( npxref%V, V_tot)
    CALL gather_to_all_dp_1D( npxref%U_c, U_c_tot)
    CALL gather_to_all_dp_1D( npxref%V_c, V_c_tot)
    CALL gather_to_all_dp_1D( Hstar_b, H_b_tot)

    ! Initialise with zeros
    laddie%divQU = 0.0_dp
    laddie%divQV = 0.0_dp

    ! == Loop over triangles ==
    ! =========================

    DO ti = mesh%ti1, mesh%ti2

      IF (laddie%mask_b( ti)) THEN

        ! Loop over all connections of triangle ti
        DO ci = 1, 3

          tj = mesh%TriC( ti, ci)
          ei = mesh%TriE( ti, ci)

          !IF (par%master) THEN
          !  WRITE( *, *) ti, ci, tj, ei, mask_gl_b_tot( tj)
          !END IF     

          ! Skip if no connecting triangle on this side
          IF (tj == 0) CYCLE

          ! Skip connection if neighbour is grounded. No flux across grounding line
          IF (mask_gl_b_tot( tj)) CYCLE

          ! The triangle-triangle vector from ti to tj
          D_x = mesh%Tricc( tj,1) - mesh%Tricc( ti,1)
          D_y = mesh%Tricc( tj,2) - mesh%Tricc( ti,2)
          D_c = SQRT( D_x**2 + D_y**2)

          ! Calculate vertically averaged ice velocity component perpendicular to this edge
          u_perp_x = U_c_tot( ei) * D_x/D_c
          u_perp_y = V_c_tot( ei) * D_y/D_c

          ! Calculate upstream momentum divergence
          ! =============================
          ! u_perp > 0: flow is exiting this triangle into triangle tj
          IF (u_perp_x > 0) THEN
            laddie%divQU( ti) = laddie%divQU( ti) + mesh%TriCw( ti, ci) * H_b_tot( ti) * U_tot( ti)* u_perp_x / mesh%TriA( ti)
          ! u_perp < 0: flow is entering this triangle into triangle tj
          ELSE
            laddie%divQU( ti) = laddie%divQU( ti) + mesh%TriCw( ti, ci) * H_b_tot( tj) * U_tot( tj)* u_perp_x / mesh%TriA( ti)
          END IF

          ! V momentum
          IF (u_perp_y > 0) THEN
            laddie%divQV( ti) = laddie%divQV( ti) + mesh%TriCw( ti, ci) * H_b_tot( ti) * V_tot( ti)* u_perp_y / mesh%TriA( ti)
          ELSE
            laddie%divQV( ti) = laddie%divQV( ti) + mesh%TriCw( ti, ci) * H_b_tot( tj) * V_tot( tj)* u_perp_y / mesh%TriA( ti)
          END IF

          !IF (par%master) THEN
          !  WRITE( *, "(I3,I3,I3,I3,A,F15.3,A,F12.6,A,F12.6,A,F18.12)") ti, ci, tj, ei, '  D_x: ', D_x, '  H_b_tot(ti)', H_b_tot( ti), '   U_tot( ti)', U_tot( ti), '   divQU(ti)', laddie%divQU( ti)*mesh%TriA( ti)
          !END IF     
        END DO ! DO ci = 1, 3

      END IF ! (laddie%mask_b( ti))

    END DO ! DO ti = mesh%ti1, mesh%ti2

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE compute_divQUV_upstream

  SUBROUTINE map_laddie_velocities_from_b_to_c_2D( mesh, u_b_partial, v_b_partial, u_c, v_c)
    ! Calculate velocities on the c-grid for solving the ice thickness equation
    ! 
    ! Uses a different scheme then the standard mapping operator, as that one is too diffusive
        
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2), INTENT(IN)    :: u_b_partial
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2), INTENT(IN)    :: v_b_partial
    REAL(dp), DIMENSION(mesh%ei1:mesh%ei2), INTENT(OUT)   :: u_c
    REAL(dp), DIMENSION(mesh%ei1:mesh%ei2), INTENT(OUT)   :: v_c
      
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'map_laddie_velocities_from_b_to_c_2D'
    REAL(dp), DIMENSION(:    ), ALLOCATABLE               :: u_b_tot, v_b_tot
    INTEGER                                               :: ei, til, tir
      
    ! Add routine to path
    CALL init_routine( routine_name)
        
    ! Allocate memory
    ALLOCATE( u_b_tot( mesh%nTri))
    ALLOCATE( v_b_tot( mesh%nTri))
        
    ! Gather the full b-grid velocity fields to all processes
    CALL gather_to_all_dp_1D( u_b_partial, u_b_tot)
    CALL gather_to_all_dp_1D( v_b_partial, v_b_tot)

    ! Map velocities from the b-grid (triangles) to the c-grid (edges)
    DO ei = mesh%ei1, mesh%ei2

      til = mesh%ETri( ei,1)
      tir = mesh%ETri( ei,2)

      IF     (til == 0 .AND. tir > 0) THEN
        u_c( ei) = u_b_tot( tir)
        v_c( ei) = v_b_tot( tir)
      ELSEIF (tir == 0 .AND. til > 0) THEN
        u_c( ei) = u_b_tot( til)
        v_c( ei) = v_b_tot( til)
      ELSEIF (til >  0 .AND. tir > 0) THEN
        u_c( ei) = (u_b_tot( til) + u_b_tot( tir)) / 2._dp
        v_c( ei) = (v_b_tot( til) + v_b_tot( tir)) / 2._dp
      ELSE
        CALL crash('something is seriously wrong with the ETri array of this mesh!')
      END IF

    END DO

    ! Clean up after yourself
    DEALLOCATE( u_b_tot)
    DEALLOCATE( v_b_tot)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_laddie_velocities_from_b_to_c_2D


END MODULE laddie_velocity


