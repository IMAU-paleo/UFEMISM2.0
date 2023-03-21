MODULE mesh_memory

  ! Routines for allocating, deallocating, extending and cropping the memory for the mesh data.

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, init_routine, finalise_routine
  USE mesh_types                                             , ONLY: type_mesh
  USE reallocate_mod                                         , ONLY: reallocate

  IMPLICIT NONE

CONTAINS

! ===== Subroutines =====
! =======================

  SUBROUTINE allocate_mesh_primary( mesh, name, nV_mem, nTri_mem, nC_mem)
    ! Allocate memory for primary mesh data (everything thats's needed for mesh creation & refinement)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                 INTENT(INOUT)     :: mesh
    CHARACTER(LEN=256),              INTENT(IN)        :: name
    INTEGER,                         INTENT(IN)        :: nV_mem, nTri_mem, nC_mem

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'allocate_mesh_primary'

    ! Add routine to path
    CALL init_routine( routine_name)

    mesh%name     = name
    mesh%nV_mem   = nV_mem
    mesh%nTri_mem = nTri_mem
    mesh%nC_mem   = nC_mem

    ! Safety: check to make sure that no memory is allocated for this mesh yet
    IF (ALLOCATED( mesh%V               ) .OR. &
        ALLOCATED( mesh%nC              ) .OR. &
        ALLOCATED( mesh%C               ) .OR. &
        ALLOCATED( mesh%niTri           ) .OR. &
        ALLOCATED( mesh%iTri            ) .OR. &
        ALLOCATED( mesh%VBI             ) .OR. &
        ALLOCATED( mesh%Tri             ) .OR. &
        ALLOCATED( mesh%Tricc           ) .OR. &
        ALLOCATED( mesh%TriC            ) .OR. &
        ALLOCATED( mesh%Tri_flip_list   ) .OR. &
        ALLOCATED( mesh%refinement_map  ) .OR. &
        ALLOCATED( mesh%refinement_stack) .OR. &
        ALLOCATED( mesh%Tri_li          )) THEN
      CALL crash('memory is already allocated for mesh "' // TRIM( name) // '"!')
    END IF

    ! Allocate memory
    ! ===============

    ! Vertex data
    ALLOCATE( mesh%V                (nV_mem,   2     ), source = 0._dp)
    ALLOCATE( mesh%nC               (nV_mem          ), source = 0    )
    ALLOCATE( mesh%C                (nV_mem,   nC_mem), source = 0    )
    ALLOCATE( mesh%niTri            (nV_mem          ), source = 0    )
    ALLOCATE( mesh%iTri             (nV_mem,   nC_mem), source = 0    )
    ALLOCATE( mesh%VBI              (nV_mem          ), source = 0    )

    ! Triangle data
    ALLOCATE( mesh%Tri              (nTri_mem, 3     ), source = 0    )
    ALLOCATE( mesh%Tricc            (nTri_mem, 2     ), source = 0._dp)
    ALLOCATE( mesh%TriC             (nTri_mem, 3     ), source = 0    )

    ! Mesh generation/refinement data
    ALLOCATE( mesh%Tri_flip_list    (nTri_mem, 2     ), source = 0    )
    ALLOCATE( mesh%refinement_map   (nTri_mem        ), source = 0    )
    ALLOCATE( mesh%refinement_stack (nTri_mem        ), source = 0    )
    ALLOCATE( mesh%Tri_li           (nTri_mem, 2     ), source = 0    )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE allocate_mesh_primary

  SUBROUTINE extend_mesh_primary( mesh, nV_mem_new, nTri_mem_new)
    ! For when we didn't allocate enough. Field by field, copy the data to a temporary array,
    ! deallocate the old field, allocate a new (bigger) one, and copy the data back.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                 INTENT(INOUT)     :: mesh
    INTEGER,                         INTENT(IN)        :: nV_mem_new, nTri_mem_new

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'extend_mesh_primary'

    ! Add routine to path
    CALL init_routine( routine_name)

    mesh%nV_mem   = nV_mem_new
    mesh%nTri_mem = nTri_mem_new

    ! Vertex data
    CALL reallocate( mesh%V               , nV_mem_new  , 2          )
    CALL reallocate( mesh%nC              , nV_mem_new               )
    CALL reallocate( mesh%C               , nV_mem_new  , mesh%nC_mem)
    CALL reallocate( mesh%niTri           , nV_mem_new               )
    CALL reallocate( mesh%iTri            , nV_mem_new  , mesh%nC_mem)
    CALL reallocate( mesh%VBI             , nV_mem_new               )

    ! Triangle data
    CALL reallocate( mesh%Tri             , nTri_mem_new, 3          )
    CALL reallocate( mesh%Tricc           , nTri_mem_new, 2          )
    CALL reallocate( mesh%TriC            , nTri_mem_new, 3          )

    ! Mesh generation/refinement data
    CALL reallocate( mesh%Tri_flip_list   , nTri_mem_new, 2          )
    CALL reallocate( mesh%refinement_map  , nTri_mem_new             )
    CALL reallocate( mesh%refinement_stack, nTri_mem_new             )
    CALL reallocate( mesh%Tri_li          , nTri_mem_new, 2          )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE extend_mesh_primary

  SUBROUTINE crop_mesh_primary( mesh)
    ! For when we allocated too much. Field by field, copy the data to a temporary array,
    ! deallocate the old field, allocate a new (smaller) one, and copy the data back.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                 INTENT(INOUT)     :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'crop_mesh_primary'

    ! Add routine to path
    CALL init_routine( routine_name)

    mesh%nV_mem   = mesh%nV
    mesh%nTri_mem = mesh%nTri

    CALL reallocate( mesh%V               , mesh%nV,   2          )
    CALL reallocate( mesh%nC              , mesh%nV               )
    CALL reallocate( mesh%C               , mesh%nV,   mesh%nC_mem)
    CALL reallocate( mesh%niTri           , mesh%nV               )
    CALL reallocate( mesh%iTri            , mesh%nV,   mesh%nC_mem)
    CALL reallocate( mesh%VBI             , mesh%nV               )

    CALL reallocate( mesh%Tri             , mesh%nTri, 3          )
    CALL reallocate( mesh%Tricc           , mesh%nTri, 2          )
    CALL reallocate( mesh%TriC            , mesh%nTri, 3          )

    ! Mesh generation/refinement data
    CALL reallocate( mesh%Tri_flip_list   , mesh%nTri, 2          )
    CALL reallocate( mesh%refinement_map  , mesh%nTri             )
    CALL reallocate( mesh%refinement_stack, mesh%nTri             )
    CALL reallocate( mesh%Tri_li          , mesh%nTri, 2          )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE crop_mesh_primary

  SUBROUTINE deallocate_mesh( mesh)
    ! Deallocate all mesh data

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                 INTENT(INOUT)     :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'deallocate_mesh'

    ! Add routine to path
    CALL init_routine( routine_name)

  ! Primary mesh data
  ! =================

    ! Vertex data
    IF (ALLOCATED( mesh%V               )) DEALLOCATE( mesh%V               )
    IF (ALLOCATED( mesh%nC              )) DEALLOCATE( mesh%nC              )
    IF (ALLOCATED( mesh%C               )) DEALLOCATE( mesh%C               )
    IF (ALLOCATED( mesh%niTri           )) DEALLOCATE( mesh%niTri           )
    IF (ALLOCATED( mesh%iTri            )) DEALLOCATE( mesh%iTri            )
    IF (ALLOCATED( mesh%VBI             )) DEALLOCATE( mesh%VBI             )

    ! Triangle data
    IF (ALLOCATED( mesh%Tri             )) DEALLOCATE( mesh%Tri             )
    IF (ALLOCATED( mesh%TriC            )) DEALLOCATE( mesh%TriC            )
    IF (ALLOCATED( mesh%Tricc           )) DEALLOCATE( mesh%Tricc           )

  ! Refinement data
  ! ===============

    IF (ALLOCATED( mesh%Tri_flip_list   )) DEALLOCATE( mesh%Tri_flip_list   )
    IF (ALLOCATED( mesh%refinement_map  )) DEALLOCATE( mesh%refinement_map  )
    IF (ALLOCATED( mesh%refinement_stack)) DEALLOCATE( mesh%refinement_stack)
    IF (ALLOCATED( mesh%Tri_li          )) DEALLOCATE( mesh%Tri_li          )

  ! Secondary mesh data (everything that can be calculated after mesh creation is finished)
  ! =======================================================================================

    ! Derived geometry data
    IF (ALLOCATED( mesh%A               )) DEALLOCATE( mesh%A               )
    IF (ALLOCATED( mesh%VorGC           )) DEALLOCATE( mesh%VorGC           )
    IF (ALLOCATED( mesh%R               )) DEALLOCATE( mesh%R               )
    IF (ALLOCATED( mesh%Cw              )) DEALLOCATE( mesh%Cw              )
    IF (ALLOCATED( mesh%TriBI           )) DEALLOCATE( mesh%TriBI           )
    IF (ALLOCATED( mesh%TriGC           )) DEALLOCATE( mesh%TriGC           )
    IF (ALLOCATED( mesh%TriA            )) DEALLOCATE( mesh%TriA            )

    ! lon/lat coordinates
    IF (ALLOCATED( mesh%lat             )) DEALLOCATE( mesh%lat             )
    IF (ALLOCATED( mesh%lon             )) DEALLOCATE( mesh%lon             )

    ! Edges (c-grid)
    IF (ALLOCATED( mesh%E               )) DEALLOCATE( mesh%E               )
    IF (ALLOCATED( mesh%VE              )) DEALLOCATE( mesh%VE              )
    IF (ALLOCATED( mesh%EV              )) DEALLOCATE( mesh%EV              )
    IF (ALLOCATED( mesh%ETri            )) DEALLOCATE( mesh%ETri            )
    IF (ALLOCATED( mesh%EBI             )) DEALLOCATE( mesh%EBI             )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE deallocate_mesh

END MODULE mesh_memory
