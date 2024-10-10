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
  USE laddie_model_types                                     , ONLY: type_laddie_model
  USE ocean_model_types                                      , ONLY: type_ocean_model
  USE reallocate_mod                                         , ONLY: reallocate_bounds
  USE mpi_distributed_memory                                 , ONLY: gather_to_all_dp_1D, gather_to_all_logical_1D
  USE mesh_operators                                         , ONLY: map_a_b_2D

  IMPLICIT NONE
    
CONTAINS
    
! ===== Main routines =====
! =========================

  SUBROUTINE compute_UV_np1( mesh, ice, laddie, dt)
    ! Integrate U and V by one time step

    ! In- and output variables

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie
    REAL(dp),                               INTENT(IN)    :: dt

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'compute_UV_np1'
    INTEGER                                               :: ti, ci, nfl, vj
    REAL(dp)                                              :: dHUdt, dHVdt, HU_next, HV_next, PGF_x, PGF_y, Hdrho_fl, Uabs
    LOGICAL, DIMENSION(mesh%nV)                           :: mask_a_tot
    REAL(dp), DIMENSION(mesh%nV)                          :: Hdrho_amb_tot
 
    ! Add routine to path
    CALL init_routine( routine_name)

    CALL gather_to_all_logical_1D( laddie%mask_a, mask_a_tot)
    CALL gather_to_all_dp_1D( laddie%Hdrho_amb, Hdrho_amb_tot)

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
                  - 0.5*grav * laddie%H_b( ti)**2 * laddie%ddrho_amb_dx_b( ti)

          PGF_y = - grav * laddie%Hdrho_amb_b( ti) * laddie%dH_dy_b( ti) &
                  + grav * laddie%Hdrho_amb_b( ti) * laddie%dHib_dy_b( ti) &
                  - 0.5*grav * laddie%H_b( ti)**2 * laddie%ddrho_amb_dy_b( ti)
        END IF

        ! == time derivatives ==
        ! ======================

        ! dHU_dt
        dHUdt = - laddie%divQU( ti) &
                + PGF_x &
                + C%uniform_laddie_coriolis_parameter * laddie%H_b( ti) * laddie%V( ti) &
                - C%laddie_drag_coefficient * laddie%U( ti) * (laddie%U( ti)**2 + laddie%V( ti)**2)**.5 &
                + laddie%viscU( ti) &
                - laddie%detr_b( ti) * laddie%U( ti)

        ! dHV_dt
        dHVdt = - laddie%divQV( ti) &
                + PGF_y &
                - C%uniform_laddie_coriolis_parameter * laddie%H_b( ti) * laddie%U( ti) &
                - C%laddie_drag_coefficient * laddie%V( ti) * (laddie%U( ti)**2 + laddie%V( ti)**2)**.5 &
                + laddie%viscV( ti) &
                - laddie%detr_b( ti) * laddie%V( ti)

        ! == next time step ==
        ! ====================

        ! HU_n = HU_n + dHU_dt * dt
        HU_next = laddie%U( ti)*laddie%H_b( ti) + dHUdt * dt
        HV_next = laddie%V( ti)*laddie%H_b( ti) + dHVdt * dt

        ! U_n = HU_n / H_n
        laddie%U_next( ti) = HU_next / laddie%H_b_next( ti)
        laddie%V_next( ti) = HV_next / laddie%H_b_next( ti)

      END IF ! (laddie%mask_b( ti))
    END DO !ti = mesh%ti1, mesh%ti2

    ! Cutoff velocities to ensure Uabs <= Uabs_max
    DO ti = mesh%ti1, mesh%ti2
      IF (laddie%mask_b( ti)) THEN
        ! Get absolute velocity
        Uabs = (laddie%U_next( ti)**2 + laddie%V_next( ti)**2)**.5
        
        ! Scale U and V 
        laddie%U_next( ti) = laddie%U_next( ti) * MIN(1.0_dp, C%laddie_velocity_maximum/Uabs)
        laddie%V_next( ti) = laddie%V_next( ti) * MIN(1.0_dp, C%laddie_velocity_maximum/Uabs)
      END IF ! (laddie%mask_b( ti))
    END DO !ti = mesh%ti1, mesh%ti2


    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE compute_UV_np1

  SUBROUTINE compute_viscUV( mesh, ice, laddie, Hstar)       
    ! Compute horizontal viscosity of momentum          
  
    ! In- and output variables                               

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)    :: Hstar

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'compute_viscUV'
    INTEGER                                               :: vi
    INTEGER                                               :: vj
    INTEGER                                               :: i
    INTEGER                                               :: j
    INTEGER                                               :: ti
    INTEGER                                               :: tj
    INTEGER                                               :: nf1
    INTEGER                                               :: nf2
    REAL(dp), DIMENSION(mesh%nTri)                        :: U_tot
    REAL(dp), DIMENSION(mesh%nTri)                        :: V_tot
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2)                :: Hstar_b
    LOGICAL, DIMENSION(mesh%nV)                           :: mask_a_tot
    LOGICAL, DIMENSION(mesh%nTri)                         :: mask_gl_b_tot
    
    ! Add routine to path
    CALL init_routine( routine_name)        

    ! Gather
    CALL gather_to_all_dp_1D( laddie%U, U_tot)            
    CALL gather_to_all_dp_1D( laddie%V, V_tot)            
    CALL gather_to_all_logical_1D( laddie%mask_a, mask_a_tot)
    CALL gather_to_all_logical_1D( laddie%mask_gl_b, mask_gl_b_tot)

    ! Map Hstar
    CALL map_a_b_2D( mesh, Hstar, Hstar_b)

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
            laddie%viscU( ti) = laddie%viscU( ti) - laddie%U( ti) * laddie%A_h( ti) * Hstar_b( ti) / mesh%TriA( ti)
            laddie%viscV( ti) = laddie%viscV( ti) - laddie%V( ti) * laddie%A_h( ti) * Hstar_b( ti) / mesh%TriA( ti)
          ! TODO add CYCLE if neighbour is ocean. Again: d/dx = d/dy = 0
          ELSE
            ! Add viscosity flux based on dU/dx and dV/dy. 
            ! Note: for grounded neighbours, U_tot( tj) = 0, meaning this is a no slip option. Can be expanded
            laddie%viscU( ti) = laddie%viscU( ti) + (U_tot( tj)-laddie%U( ti)) * laddie%A_h( ti) * Hstar_b( ti) / mesh%TriA( ti)
            laddie%viscV( ti) = laddie%viscV( ti) + (V_tot( tj)-laddie%V( ti)) * laddie%A_h( ti) * Hstar_b( ti) / mesh%TriA( ti)
          END IF
        END DO

      END IF !(laddie%mask_b( ti)
    END DO !ti = mesh%ti1, mesh%ti2
      
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE compute_viscUV

END MODULE laddie_velocity

