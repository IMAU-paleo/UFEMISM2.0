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
  USE reallocate_mod                                         , ONLY: reallocate_bounds

  IMPLICIT NONE
    
CONTAINS
    
! ===== Main routines =====
! =========================
  SUBROUTINE initialise_laddie_model( mesh, laddie)
    ! Initialise the laddie model

    ! In- and output variables

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_laddie_model),                INTENT(OUT)   :: laddie

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'initialise_laddie_model'
 
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
 
    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_laddie_model

END MODULE laddie_main

