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
    INTEGER                                               :: ti
    REAL(dp)                                              :: dHUdt
    REAL(dp)                                              :: dHVdt
    REAL(dp)                                              :: HU_next
    REAL(dp)                                              :: HV_next
 
    ! Add routine to path
    CALL init_routine( routine_name)

    ! == Integrate U ==
    ! =================

    ! Loop over vertices
    ! TODO apply only to laddie%mask_b
    DO ti = mesh%ti1, mesh%ti2

      ! Get dHU_dt
      dHUdt = - laddie%divQU( ti) &
              - grav * laddie%Hdrho_amb_b( ti) * laddie%dH_dx_b( ti) &
              + grav * laddie%Hdrho_amb_b( ti) * laddie%dHib_dx_b( ti) &
              - .5*grav * laddie%H_b( ti)**2 * laddie%ddrho_amb_dx_b( ti) &
              + C%uniform_laddie_coriolis_parameter * laddie%H_b( ti) * laddie%V( ti) &
              - C%laddie_drag_coefficient * laddie%U( ti) * (laddie%U( ti)**2 + laddie%V( ti)**2)**.5 &
              + laddie%viscU( ti) &
              - laddie%detr_b( ti) * laddie%U( ti)

      ! HU_n = HU_n + dHU_dt * dt
      HU_next = laddie%U( ti)*laddie%H_b( ti) + dHUdt * dt

      laddie%U_next( ti) = HU_next / laddie%H_b_next( ti)

    END DO !vi = mesh%vi, mesh%v2

    ! == Integrate V ==
    ! =================

    ! Loop over vertices
    ! TODO apply only to laddie%mask_b
    DO ti = mesh%ti1, mesh%ti2

      ! Get dHV_dt
      dHVdt = - laddie%divQV( ti) &
              - grav * laddie%Hdrho_amb_b( ti) * laddie%dH_dy_b( ti) &
              + grav * laddie%Hdrho_amb_b( ti) * laddie%dHib_dy_b( ti) &
              - .5*grav * laddie%H_b( ti)**2 * laddie%ddrho_amb_dy_b( ti) &
              - C%uniform_laddie_coriolis_parameter * laddie%H_b( ti) * laddie%V( ti) &
              - C%laddie_drag_coefficient * laddie%V( ti) * (laddie%U( ti)**2 + laddie%V( ti)**2)**.5 &
              + laddie%viscV( ti) &
              - laddie%detr_b( ti) * laddie%V( ti)

      ! HV_n = HV_n + dHV_dt * dt
      HV_next = laddie%V( ti)*laddie%H_b( ti) + dHVdt * dt

      laddie%V_next( ti) = HV_next / laddie%H_b_next( ti)

    END DO !vi = mesh%vi, mesh%v2

    ! TODO add cutoff by multiplying UVnext by uabsmax/abs(uabs)

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
    
    ! Add routine to path
    CALL init_routine( routine_name)        

    ! Gather
    CALL gather_to_all_dp_1D( laddie%U, U_tot)            
    CALL gather_to_all_dp_1D( laddie%V, V_tot)            
    CALL gather_to_all_logical_1D( laddie%mask_a, mask_a_tot)

    ! Map Hstar
    CALL map_a_b_2D( mesh, Hstar, Hstar_b)

    ! Loop over triangles                                  
    DO ti = mesh%ti1, mesh%ti2
      nf1 = 0

      ! Loop over connecing vertices and check whether any is floating
      DO i = 1, 3
        vi = mesh%Tri( ti, i)
        IF (mask_a_tot( vi)) THEN
          nf1 = nf1 + 1
        END IF
      END DO
      
      IF (nf1 > 0) THEN
        ! Get viscosity parameter
        laddie%A_h( ti) = C%laddie_viscosity
        ! TODO add scalable options
    
        ! Initialise at 0
        laddie%viscU( ti) = 0.0_dp
        laddie%viscV( ti) = 0.0_dp
      
        ! Loop over connected triangles
        DO i = 1, 3
          tj = mesh%TriC( ti, i)

          IF (tj>0) THEN
            ! Check whether at least 2 of 3 vertices are floating to use neighbouring velocity, otherwise d/dx = d/dy = 0
            nf2 = 0
            DO j = 1, 3
              vj = mesh%Tri( tj, j)
              IF (mask_a_tot( vj)) THEN
                nf2 = nf2 + 1
              END IF
            END DO

            ! Can simply skip non-floating vertices to ensure d/dx = d/dy = 0 at boundaries
            IF (nf2 > 1) THEN
              laddie%viscU( ti) = laddie%viscU( ti) + (U_tot( tj)-laddie%U( ti)) * laddie%A_h( ti) * Hstar_b( ti) / mesh%TriA( ti)
              laddie%viscV( ti) = laddie%viscV( ti) + (V_tot( tj)-laddie%V( ti)) * laddie%A_h( ti) * Hstar_b( ti) / mesh%TriA( ti)
            END IF
          END IF
        END DO

      END IF
    END DO
      
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE compute_viscUV

END MODULE laddie_velocity

