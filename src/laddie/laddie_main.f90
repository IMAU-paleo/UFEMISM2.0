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
  USE reallocate_mod                                         , ONLY: reallocate_bounds
  USE ocean_utilities                                        , ONLY: interpolate_ocean_depth

  IMPLICIT NONE
    
CONTAINS
    
! ===== Main routines =====
! =========================
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
 
    ! Allocate memory for main variables
    ALLOCATE( laddie%H( mesh%vi1:mesh%vi2))
    ALLOCATE( laddie%U( mesh%vi1:mesh%vi2))
    ALLOCATE( laddie%V( mesh%vi1:mesh%vi2))
    ALLOCATE( laddie%T( mesh%vi1:mesh%vi2))
    ALLOCATE( laddie%S( mesh%vi1:mesh%vi2))

    laddie%H = C%laddie_thickness_minimum
    laddie%U = 0._dp
    laddie%V = 0._dp
    laddie%T = 0._dp
    laddie%S = 0._dp

    ! Initialise ambient T and S
    ALLOCATE( laddie%T_amb( mesh%vi1:mesh%vi2))
    ALLOCATE( laddie%S_amb( mesh%vi1:mesh%vi2))
    laddie%T_amb = 0._dp
    laddie%S_amb = 0._dp

    ! Get T and S at layer base
    DO vi = mesh%vi1, mesh%vi2
       IF (ice%mask_floating_ice( vi)) THEN
         CALL interpolate_ocean_depth( C%nz_ocean, C%z_ocean, ocean%T( vi,:), laddie%H( vi) - ice%Hib (vi), laddie%T_amb( vi))
         CALL interpolate_ocean_depth( C%nz_ocean, C%z_ocean, ocean%S( vi,:), laddie%H( vi) - ice%Hib (vi), laddie%S_amb( vi))
       END IF
    END DO

    ! Initialise main T and S
    DO vi = mesh%vi1, mesh%vi2
       IF (ice%mask_floating_ice( vi)) THEN
         laddie%T( vi) = laddie%T_amb( vi)
         laddie%S( vi) = laddie%S_amb( vi) -0.1_dp
       END IF
    END DO
 
    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_laddie_model

END MODULE laddie_main

