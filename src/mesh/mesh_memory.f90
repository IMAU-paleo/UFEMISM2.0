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
    !         (which should never be the case!)

    IF (ALLOCATED( mesh%V         ) .OR. &
        ALLOCATED( mesh%nC        ) .OR. &
        ALLOCATED( mesh%C         ) .OR. &
        ALLOCATED( mesh%niTri     ) .OR. &
        ALLOCATED( mesh%iTri      ) .OR. &
        ALLOCATED( mesh%edge_index) .OR. &
        ALLOCATED( mesh%Tri       ) .OR. &
        ALLOCATED( mesh%Tricc     ) .OR. &
        ALLOCATED( mesh%TriC      )) CALL crash('memory is already allocated for mesh "' // TRIM( name) // '"!')

    ! Allocate memory
    ! ===============

    ! Vertex data
    ALLOCATE( mesh%V              (nV_mem,   2     ), source = 0._dp)
    ALLOCATE( mesh%nC             (nV_mem          ), source = 0    )
    ALLOCATE( mesh%C              (nV_mem,   nC_mem), source = 0    )
    ALLOCATE( mesh%niTri          (nV_mem          ), source = 0    )
    ALLOCATE( mesh%iTri           (nV_mem,   nC_mem), source = 0    )
    ALLOCATE( mesh%edge_index     (nV_mem          ), source = 0    )

    ! Triangle data
    ALLOCATE( mesh%Tri            (nTri_mem, 3     ), source = 0    )
    ALLOCATE( mesh%Tricc          (nTri_mem, 2     ), source = 0._dp)
    ALLOCATE( mesh%TriC           (nTri_mem, 3     ), source = 0    )

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
    CALL reallocate( mesh%V            , nv_mem_new  , 2          )
    CALL reallocate( mesh%V            , nV_mem_new  , 2          )
    CALL reallocate( mesh%nC           , nV_mem_new               )
    CALL reallocate( mesh%C            , nV_mem_new  , mesh%nC_mem)
    CALL reallocate( mesh%niTri        , nV_mem_new               )
    CALL reallocate( mesh%iTri         , nV_mem_new  , mesh%nC_mem)
    CALL reallocate( mesh%edge_index   , nV_mem_new               )

    ! Triangle data
    CALL reallocate( mesh%Tri          , nTri_mem_new, 3          )
    CALL reallocate( mesh%Tricc        , nTri_mem_new, 2          )
    CALL reallocate( mesh%TriC         , nTri_mem_new, 3          )

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

    CALL reallocate( mesh%V             , mesh%nV,   2          )
    CALL reallocate( mesh%nC            , mesh%nV               )
    CALL reallocate( mesh%C             , mesh%nV,   mesh%nC_mem)
    CALL reallocate( mesh%niTri         , mesh%nV               )
    CALL reallocate( mesh%iTri          , mesh%nV,   mesh%nC_mem)
    CALL reallocate( mesh%edge_index    , mesh%nV               )

    CALL reallocate( mesh%Tri           , mesh%nTri, 3          )
    CALL reallocate( mesh%Tricc         , mesh%nTri, 2          )
    CALL reallocate( mesh%TriC          , mesh%nTri, 3          )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE crop_mesh_primary

!  SUBROUTINE allocate_mesh_secondary(     mesh)
!    ! Allocate memory for mesh data
!
!    IMPLICIT NONE
!
!    ! In/output variables:
!    TYPE(type_mesh),                 INTENT(INOUT)     :: mesh
!
!    ! Local variables:
!    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'allocate_mesh_secondary'
!
!    ! Add routine to path
!    CALL init_routine( routine_name)
!
!    if (allocated(mesh% A    )) deallocate(mesh% A    )
!    if (allocated(mesh% VorGC)) deallocate(mesh% VorGC)
!    if (allocated(mesh% R    )) deallocate(mesh% R    )
!    if (allocated(mesh% Cw   )) deallocate(mesh% Cw   )
!    if (allocated(mesh% TriA )) deallocate(mesh% TriA )
!    if (allocated(mesh% TriGC)) deallocate(mesh% TriGC)
!    if (allocated(mesh% lat  )) deallocate(mesh% lat  )
!    if (allocated(mesh% lon  )) deallocate(mesh% lon  )
!
!    allocate(mesh%A    ( 1       :mesh%nV    )) ! Globally used kinda
!    allocate(mesh%VorGC( mesh%vi1:mesh%vi2,2 ))
!    allocate(mesh%R    ( 1       :mesh%nV    )) ! Globally used kinda
!    allocate(mesh%Cw   ( mesh%vi1:mesh%vi2,mesh%nC_mem ))
!
!    allocate(mesh%TriA ( mesh%ti1:mesh%ti2   ))
!    allocate(mesh%TriGC( 1:mesh%nTri,2       )) ! Globally needed
!
!    allocate(mesh%lat  ( mesh%vi1:mesh%vi2   ))
!    allocate(mesh%lon  ( mesh%vi1:mesh%vi2   ))
!
!    ! Finalise routine path
!    CALL finalise_routine( routine_name, n_extra_windows_expected = 9)
!
!  END SUBROUTINE allocate_mesh_secondary
!
!  SUBROUTINE deallocate_mesh_all(         mesh)
!    ! Deallocate memory for mesh data
!
!    IMPLICIT NONE
!
!    ! In/output variables:
!    TYPE(type_mesh),                 INTENT(INOUT)     :: mesh
!
!    ! Local variables:
!    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'deallocate_mesh_all'
!
!    ! Add routine to path
!    CALL init_routine( routine_name)
!
!    ! Basic meta properties
!    ! =====================
!
!!   CALL deallocate_shared( mesh%wlambda_M)
!!   CALL deallocate_shared( mesh%wphi_M)
!!   CALL deallocate_shared( mesh%walpha_stereo)
!!   CALL deallocate_shared( mesh%wxmin)
!!   CALL deallocate_shared( mesh%wxmax)
!!   CALL deallocate_shared( mesh%wymin)
!!   CALL deallocate_shared( mesh%wymax)
!!   CALL deallocate_shared( mesh%wtol_dist)
!!   CALL deallocate_shared( mesh%wnV_mem)
!!   CALL deallocate_shared( mesh%wnTri_mem)
!!   CALL deallocate_shared( mesh%wnC_mem)
!!   CALL deallocate_shared( mesh%wnV)
!!   CALL deallocate_shared( mesh%wnTri)
!!   CALL deallocate_shared( mesh%wperturb_dir)
!!   CALL deallocate_shared( mesh%walpha_min)
!!   CALL deallocate_shared( mesh%wdz_max_ice)
!!   CALL deallocate_shared( mesh%wres_max)
!!   CALL deallocate_shared( mesh%wres_max_margin)
!!   CALL deallocate_shared( mesh%wres_max_gl)
!!   CALL deallocate_shared( mesh%wres_max_cf)
!!   CALL deallocate_shared( mesh%wres_max_mountain)
!!   CALL deallocate_shared( mesh%wres_max_coast)
!!   CALL deallocate_shared( mesh%wres_min)
!!   CALL deallocate_shared( mesh%wresolution_min)
!!   CALL deallocate_shared( mesh%wresolution_max)
!
!!   ! Primary mesh data (needed for mesh creation & refinement)
!!   ! =========================================================
!!
!!   CALL deallocate_shared( mesh%wV)
!!   CALL deallocate_shared( mesh%wnC)
!!   CALL deallocate_shared( mesh%wC)
!!   CALL deallocate_shared( mesh%wniTri)
!!   CALL deallocate_shared( mesh%wiTri)
!!   CALL deallocate_shared( mesh%wedge_index)
!!   CALL deallocate_shared( mesh%wmesh_old_ti_in)
!!
!!   CALL deallocate_shared( mesh%wTri)
!!   CALL deallocate_shared( mesh%wTricc)
!!   CALL deallocate_shared( mesh%wTriC)
!!   CALL deallocate_shared( mesh%wTri_edge_index)
!!
!!   CALL deallocate_shared( mesh%wTriflip)
!!   CALL deallocate_shared( mesh%wRefMap)
!!   CALL deallocate_shared( mesh%wRefStack)
!!   CALL deallocate_shared( mesh%wRefStackN)
!!
!!   CALL deallocate_shared( mesh%wVMap)
!!   CALL deallocate_shared( mesh%wVStack1)
!!   CALL deallocate_shared( mesh%wVStack2)
!!   CALL deallocate_shared( mesh%wTriMap)
!!   CALL deallocate_shared( mesh%wTriStack1)
!!   CALL deallocate_shared( mesh%wTriStack2)
!!
!!   CALL deallocate_shared( mesh%wnPOI               )
!!   CALL deallocate_shared( mesh%wPOI_coordinates    )
!!   CALL deallocate_shared( mesh%wPOI_XY_coordinates )
!!   CALL deallocate_shared( mesh%wPOI_resolutions    )
!!   CALL deallocate_shared( mesh%wPOI_vi             )
!!   CALL deallocate_shared( mesh%wPOI_w              )
!!
!!   ! Secondary mesh data
!!   ! ===================
!
!!   CALL deallocate_shared( mesh%wA)
!!   CALL deallocate_shared( mesh%wVorGC)
!!   CALL deallocate_shared( mesh%wR)
!!   CALL deallocate_shared( mesh%wCw)
!!   CALL deallocate_shared( mesh%wTriGC)
!!   CALL deallocate_shared( mesh%wTriA)
!!
!!   CALL deallocate_shared( mesh%wlat)
!!   CALL deallocate_shared( mesh%wlon)
!!
!!   CALL deallocate_shared( mesh%wnAc)
!!   CALL deallocate_shared( mesh%wVAc)
!!   CALL deallocate_shared( mesh%wAci)
!!   CALL deallocate_shared( mesh%wiAci)
!!   CALL deallocate_shared( mesh%wedge_index_Ac)
!
!    CALL MatDestroy( mesh%M_map_a_b             , perr)
!    CALL MatDestroy( mesh%M_map_a_c             , perr)
!    CALL MatDestroy( mesh%M_map_b_a             , perr)
!    CALL MatDestroy( mesh%M_map_b_c             , perr)
!    CALL MatDestroy( mesh%M_map_c_a             , perr)
!    CALL MatDestroy( mesh%M_map_c_b             , perr)
!
!    CALL MatDestroy( mesh%M_ddx_a_a             , perr)
!    CALL MatDestroy( mesh%M_ddx_a_b             , perr)
!    CALL MatDestroy( mesh%M_ddx_a_c             , perr)
!    CALL MatDestroy( mesh%M_ddx_b_a             , perr)
!    CALL MatDestroy( mesh%M_ddx_b_b             , perr)
!    CALL MatDestroy( mesh%M_ddx_b_c             , perr)
!    CALL MatDestroy( mesh%M_ddx_c_a             , perr)
!    CALL MatDestroy( mesh%M_ddx_c_b             , perr)
!    CALL MatDestroy( mesh%M_ddx_c_c             , perr)
!
!    CALL MatDestroy( mesh%M_ddy_a_a             , perr)
!    CALL MatDestroy( mesh%M_ddy_a_b             , perr)
!    CALL MatDestroy( mesh%M_ddy_a_c             , perr)
!    CALL MatDestroy( mesh%M_ddy_b_a             , perr)
!    CALL MatDestroy( mesh%M_ddy_b_b             , perr)
!    CALL MatDestroy( mesh%M_ddy_b_c             , perr)
!    CALL MatDestroy( mesh%M_ddy_c_a             , perr)
!    CALL MatDestroy( mesh%M_ddy_c_b             , perr)
!    CALL MatDestroy( mesh%M_ddy_c_c             , perr)
!
!    CALL MatDestroy( mesh%M2_ddx_b_b            , perr)
!    CALL MatDestroy( mesh%M2_ddy_b_b            , perr)
!    CALL MatDestroy( mesh%M2_d2dx2_b_b          , perr)
!    CALL MatDestroy( mesh%M2_d2dxdy_b_b         , perr)
!    CALL MatDestroy( mesh%M2_d2dy2_b_b          , perr)
!
!    CALL deallocate_matrix_CSR( mesh%M2_ddx_b_b_CSR    )
!    CALL deallocate_matrix_CSR( mesh%M2_ddy_b_b_CSR    )
!    CALL deallocate_matrix_CSR( mesh%M2_d2dx2_b_b_CSR  )
!    CALL deallocate_matrix_CSR( mesh%M2_d2dxdy_b_b_CSR )
!    CALL deallocate_matrix_CSR( mesh%M2_d2dy2_b_b_CSR  )
!
!    CALL MatDestroy( mesh%M_Neumann_BC_b        , perr)
!    CALL deallocate_matrix_CSR( mesh%M_Neumann_BC_b_CSR)
!
!!    CALL deallocate_shared( mesh%wnV_transect)
!!    CALL deallocate_shared( mesh%wvi_transect)
!!    CALL deallocate_shared( mesh%ww_transect )
!
!    ! Finalise routine path
!    CALL finalise_routine( routine_name)
!
!   END SUBROUTINE deallocate_mesh_all

END MODULE mesh_memory
