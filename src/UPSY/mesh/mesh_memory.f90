MODULE mesh_memory

  ! Routines for allocating, deallocating, extending and cropping the memory for the mesh data.

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
  USE mesh_types                                             , ONLY: type_mesh
  USE reallocate_mod                                         , ONLY: reallocate
  USE CSR_matrix_basics                            , ONLY: deallocate_matrix_CSR_dist
  use mpi_distributed_shared_memory, only: deallocate_dist_shared

  IMPLICIT NONE

CONTAINS

! ===== Subroutines =====
! =======================

  SUBROUTINE allocate_mesh_primary( mesh, name, nV_mem, nTri_mem)
    ! Allocate memory for primary mesh data (everything thats's needed for mesh creation & refinement)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                 INTENT(INOUT)     :: mesh
    CHARACTER(LEN=*),                INTENT(IN)        :: name
    INTEGER,                         INTENT(IN)        :: nV_mem, nTri_mem

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'allocate_mesh_primary'

    ! Add routine to path
    CALL init_routine( routine_name)

    mesh%name     = trim(name)
    mesh%nV_mem   = nV_mem
    mesh%nTri_mem = nTri_mem

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
    ALLOCATE( mesh%V                (nV_mem,   2          ), source = 0._dp)
    ALLOCATE( mesh%nC               (nV_mem               ), source = 0    )
    ALLOCATE( mesh%C                (nV_mem,   mesh%nC_mem), source = 0    )
    ALLOCATE( mesh%niTri            (nV_mem               ), source = 0    )
    ALLOCATE( mesh%iTri             (nV_mem,   mesh%nC_mem), source = 0    )
    ALLOCATE( mesh%VBI              (nV_mem               ), source = 0    )

    ! Triangle data
    ALLOCATE( mesh%Tri              (nTri_mem, 3     ), source = 0    )
    ALLOCATE( mesh%Tricc            (nTri_mem, 2     ), source = 0._dp)
    ALLOCATE( mesh%TriC             (nTri_mem, 3     ), source = 0    )

    ! Mesh generation/refinement data
    ALLOCATE( mesh%Tri_flip_list    (nTri_mem*2, 2   ), source = 0    )
    ALLOCATE( mesh%refinement_map   (nTri_mem        ), source = 0    )
    ALLOCATE( mesh%refinement_stack (nTri_mem        ), source = 0    )
    mesh%refinement_stackN = 0
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
    CALL reallocate( mesh%Tri_flip_list   , nTri_mem_new*2, 2        )
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

    IF (mesh%nV_mem > mesh%nV) THEN

      mesh%nV_mem   = mesh%nV

      ! Vertex data
      CALL reallocate( mesh%V               , mesh%nV,   2          )
      CALL reallocate( mesh%nC              , mesh%nV               )
      CALL reallocate( mesh%C               , mesh%nV,   mesh%nC_mem)
      CALL reallocate( mesh%niTri           , mesh%nV               )
      CALL reallocate( mesh%iTri            , mesh%nV,   mesh%nC_mem)
      CALL reallocate( mesh%VBI             , mesh%nV               )

    END IF ! IF (mesh%nV_mem > mesh%nV) THEN

    IF (mesh%nTri_mem > mesh%nTri) THEN

      mesh%nTri_mem = mesh%nTri

      ! Triangle data
      CALL reallocate( mesh%Tri             , mesh%nTri, 3          )
      CALL reallocate( mesh%Tricc           , mesh%nTri, 2          )
      CALL reallocate( mesh%TriC            , mesh%nTri, 3          )

      ! Mesh generation/refinement data
      CALL reallocate( mesh%Tri_flip_list   , mesh%nTri*2, 2        )
      CALL reallocate( mesh%refinement_map  , mesh%nTri             )
      CALL reallocate( mesh%refinement_stack, mesh%nTri             )
      CALL reallocate( mesh%Tri_li          , mesh%nTri, 2          )

    END IF ! IF (mesh%nTri_mem > mesh%nTri) THEN

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

    ! Basic meta properties
    ! =====================

    IF (ALLOCATED( mesh%zeta            )) DEALLOCATE( mesh%zeta            )
    IF (ALLOCATED( mesh%zeta_stag       )) DEALLOCATE( mesh%zeta_stag       )

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
    IF (ALLOCATED( mesh%D_x             )) DEALLOCATE( mesh%D_x             )
    IF (ALLOCATED( mesh%D_y             )) DEALLOCATE( mesh%D_y             )
    IF (ALLOCATED( mesh%D               )) DEALLOCATE( mesh%D               )
    IF (ALLOCATED( mesh%TriCw           )) DEALLOCATE( mesh%TriCw           )
    IF (ALLOCATED( mesh%TriBI           )) DEALLOCATE( mesh%TriBI           )
    IF (ALLOCATED( mesh%TriGC           )) DEALLOCATE( mesh%TriGC           )
    IF (ALLOCATED( mesh%TriA            )) DEALLOCATE( mesh%TriA            )
    IF (ALLOCATED( mesh%TriD_x          )) DEALLOCATE( mesh%TriD_x          )
    IF (ALLOCATED( mesh%TriD_y          )) DEALLOCATE( mesh%TriD_y          )
    IF (ALLOCATED( mesh%TriD            )) DEALLOCATE( mesh%TriD            )

    ! lon/lat coordinates
    IF (ALLOCATED( mesh%lat             )) DEALLOCATE( mesh%lat             )
    IF (ALLOCATED( mesh%lon             )) DEALLOCATE( mesh%lon             )

    ! Edges (c-grid)
    IF (ALLOCATED( mesh%E               )) DEALLOCATE( mesh%E               )
    IF (ALLOCATED( mesh%VE              )) DEALLOCATE( mesh%VE              )
    IF (ALLOCATED( mesh%EV              )) DEALLOCATE( mesh%EV              )
    IF (ALLOCATED( mesh%ETri            )) DEALLOCATE( mesh%ETri            )
    IF (ALLOCATED( mesh%TriE            )) DEALLOCATE( mesh%TriE            )
    IF (ALLOCATED( mesh%EBI             )) DEALLOCATE( mesh%EBI             )

    ! Deallocate buffer shared memory
    if (associated( mesh%buffer1_d_a_nih )) call deallocate_dist_shared( mesh%buffer1_d_a_nih , mesh%wbuffer1_d_a_nih )
    if (associated( mesh%buffer2_d_a_nih )) call deallocate_dist_shared( mesh%buffer2_d_a_nih , mesh%wbuffer2_d_a_nih )
    if (associated( mesh%buffer1_d_ak_nih)) call deallocate_dist_shared( mesh%buffer1_d_ak_nih, mesh%wbuffer1_d_ak_nih)
    if (associated( mesh%buffer2_d_ak_nih)) call deallocate_dist_shared( mesh%buffer2_d_ak_nih, mesh%wbuffer2_d_ak_nih)
    if (associated( mesh%buffer1_d_b_nih )) call deallocate_dist_shared( mesh%buffer1_d_b_nih , mesh%wbuffer1_d_b_nih )
    if (associated( mesh%buffer2_d_b_nih )) call deallocate_dist_shared( mesh%buffer2_d_b_nih , mesh%wbuffer2_d_b_nih )
    if (associated( mesh%buffer1_d_bk_nih)) call deallocate_dist_shared( mesh%buffer1_d_bk_nih, mesh%wbuffer1_d_bk_nih)
    if (associated( mesh%buffer2_d_bk_nih)) call deallocate_dist_shared( mesh%buffer2_d_bk_nih, mesh%wbuffer2_d_bk_nih)
    if (associated( mesh%buffer1_d_c_nih )) call deallocate_dist_shared( mesh%buffer1_d_c_nih , mesh%wbuffer1_d_c_nih )
    if (associated( mesh%buffer2_d_c_nih )) call deallocate_dist_shared( mesh%buffer2_d_c_nih , mesh%wbuffer2_d_c_nih )
    if (associated( mesh%buffer1_d_ck_nih)) call deallocate_dist_shared( mesh%buffer1_d_ck_nih, mesh%wbuffer1_d_ck_nih)
    if (associated( mesh%buffer2_d_ck_nih)) call deallocate_dist_shared( mesh%buffer2_d_ck_nih, mesh%wbuffer2_d_ck_nih)

    ! Matrix operators
    ! ================

    ! Grid-cell-to-matrix-row translation tables

    ! a-grid (vertices)
    IF (ALLOCATED( mesh%n2vi            )) DEALLOCATE( mesh%n2vi            )
    IF (ALLOCATED( mesh%n2viuv          )) DEALLOCATE( mesh%n2viuv          )
    IF (ALLOCATED( mesh%n2vik           )) DEALLOCATE( mesh%n2vik           )
    IF (ALLOCATED( mesh%n2vikuv         )) DEALLOCATE( mesh%n2vikuv         )
    IF (ALLOCATED( mesh%n2viks          )) DEALLOCATE( mesh%n2viks          )
    IF (ALLOCATED( mesh%n2viksuv        )) DEALLOCATE( mesh%n2viksuv        )
    IF (ALLOCATED( mesh%vi2n            )) DEALLOCATE( mesh%vi2n            )
    IF (ALLOCATED( mesh%viuv2n          )) DEALLOCATE( mesh%viuv2n          )
    IF (ALLOCATED( mesh%vik2n           )) DEALLOCATE( mesh%vik2n           )
    IF (ALLOCATED( mesh%vikuv2n         )) DEALLOCATE( mesh%vikuv2n         )
    IF (ALLOCATED( mesh%viks2n          )) DEALLOCATE( mesh%viks2n          )
    IF (ALLOCATED( mesh%viksuv2n        )) DEALLOCATE( mesh%viksuv2n        )

    ! b-grid (triangles)
    IF (ALLOCATED( mesh%n2ti            )) DEALLOCATE( mesh%n2ti            )
    IF (ALLOCATED( mesh%n2tiuv          )) DEALLOCATE( mesh%n2tiuv          )
    IF (ALLOCATED( mesh%n2tik           )) DEALLOCATE( mesh%n2tik           )
    IF (ALLOCATED( mesh%n2tikuv         )) DEALLOCATE( mesh%n2tikuv         )
    IF (ALLOCATED( mesh%n2tiks          )) DEALLOCATE( mesh%n2tiks          )
    IF (ALLOCATED( mesh%n2tiksuv        )) DEALLOCATE( mesh%n2tiksuv        )
    IF (ALLOCATED( mesh%ti2n            )) DEALLOCATE( mesh%ti2n            )
    IF (ALLOCATED( mesh%tiuv2n          )) DEALLOCATE( mesh%tiuv2n          )
    IF (ALLOCATED( mesh%tik2n           )) DEALLOCATE( mesh%tik2n           )
    IF (ALLOCATED( mesh%tikuv2n         )) DEALLOCATE( mesh%tikuv2n         )
    IF (ALLOCATED( mesh%tiks2n          )) DEALLOCATE( mesh%tiks2n          )
    IF (ALLOCATED( mesh%tiksuv2n        )) DEALLOCATE( mesh%tiksuv2n        )

    ! c-grid (edges)
    IF (ALLOCATED( mesh%n2ei            )) DEALLOCATE( mesh%n2ei            )
    IF (ALLOCATED( mesh%n2eiuv          )) DEALLOCATE( mesh%n2eiuv          )
    IF (ALLOCATED( mesh%n2eik           )) DEALLOCATE( mesh%n2eik           )
    IF (ALLOCATED( mesh%n2eikuv         )) DEALLOCATE( mesh%n2eikuv         )
    IF (ALLOCATED( mesh%n2eiks          )) DEALLOCATE( mesh%n2eiks          )
    IF (ALLOCATED( mesh%n2eiksuv        )) DEALLOCATE( mesh%n2eiksuv        )
    IF (ALLOCATED( mesh%ei2n            )) DEALLOCATE( mesh%ei2n            )
    IF (ALLOCATED( mesh%eiuv2n          )) DEALLOCATE( mesh%eiuv2n          )
    IF (ALLOCATED( mesh%eik2n           )) DEALLOCATE( mesh%eik2n           )
    IF (ALLOCATED( mesh%eikuv2n         )) DEALLOCATE( mesh%eikuv2n         )
    IF (ALLOCATED( mesh%eiks2n          )) DEALLOCATE( mesh%eiks2n          )
    IF (ALLOCATED( mesh%eiksuv2n        )) DEALLOCATE( mesh%eiksuv2n        )

    ! Basic mapping and gradient operators

    ! a-grid (vertices) to a-grid (vertices)
    CALL deallocate_matrix_CSR_dist( mesh%M_ddx_a_a)
    CALL deallocate_matrix_CSR_dist( mesh%M_ddy_a_a)
    ! a-grid (vertices) to b-grid (triangles)
    CALL deallocate_matrix_CSR_dist( mesh%M_map_a_b)
    CALL deallocate_matrix_CSR_dist( mesh%M_ddx_a_b)
    CALL deallocate_matrix_CSR_dist( mesh%M_ddy_a_b)
    ! b-grid (triangles) to a-grid (vertices)
    CALL deallocate_matrix_CSR_dist( mesh%M_map_b_a)
    CALL deallocate_matrix_CSR_dist( mesh%M_ddx_b_a)
    CALL deallocate_matrix_CSR_dist( mesh%M_ddy_b_a)
    ! b-grid (triangles) to b-grid (triangles)
    CALL deallocate_matrix_CSR_dist( mesh%M_ddx_b_b)
    CALL deallocate_matrix_CSR_dist( mesh%M_ddy_b_b)

    ! b-grid (triangles) to b-grid (triangles), 2nd-order accurate
    CALL deallocate_matrix_CSR_dist( mesh%M2_ddx_b_b   )
    CALL deallocate_matrix_CSR_dist( mesh%M2_ddy_b_b   )
    CALL deallocate_matrix_CSR_dist( mesh%M2_d2dx2_b_b )
    CALL deallocate_matrix_CSR_dist( mesh%M2_d2dxdy_b_b)
    CALL deallocate_matrix_CSR_dist( mesh%M2_d2dy2_b_b )

    ! Operators on the zeta grids
    CALL deallocate_matrix_CSR_dist( mesh%M_ddzeta_k_k_1D)
    CALL deallocate_matrix_CSR_dist( mesh%M_d2dzeta2_k_k_1D)
    CALL deallocate_matrix_CSR_dist( mesh%M_map_k_ks_1D)
    CALL deallocate_matrix_CSR_dist( mesh%M_ddzeta_k_ks_1D)
    CALL deallocate_matrix_CSR_dist( mesh%M_map_ks_k_1D)
    CALL deallocate_matrix_CSR_dist( mesh%M_ddzeta_ks_k_1D)

    ! Zeta operators in tridiagonal form for efficient use in thermodynamics
    IF (ALLOCATED( mesh%M_ddzeta_k_k_ldiag  )) DEALLOCATE( mesh%M_ddzeta_k_k_ldiag  )
    IF (ALLOCATED( mesh%M_ddzeta_k_k_diag   )) DEALLOCATE( mesh%M_ddzeta_k_k_diag   )
    IF (ALLOCATED( mesh%M_ddzeta_k_k_udiag  )) DEALLOCATE( mesh%M_ddzeta_k_k_udiag  )
    IF (ALLOCATED( mesh%M_d2dzeta2_k_k_ldiag)) DEALLOCATE( mesh%M_d2dzeta2_k_k_ldiag)
    IF (ALLOCATED( mesh%M_d2dzeta2_k_k_diag )) DEALLOCATE( mesh%M_d2dzeta2_k_k_diag )
    IF (ALLOCATED( mesh%M_d2dzeta2_k_k_udiag)) DEALLOCATE( mesh%M_d2dzeta2_k_k_udiag)

    ! 3-D gradient operators

    ! bk to ak (for calculating the horizontal stretch/shear strain rates in the BPA)
    CALL deallocate_matrix_CSR_dist( mesh%M_ddx_bk_ak)
    CALL deallocate_matrix_CSR_dist( mesh%M_ddy_bk_ak)

    ! ak to bk (for calculating the horizontal gradients of the effective viscosity in the BPA)
    CALL deallocate_matrix_CSR_dist( mesh%M_ddx_ak_bk)
    CALL deallocate_matrix_CSR_dist( mesh%M_ddy_ak_bk)

    ! bk to bks (for calculating the vertical shear strain rates in the BPA)
    CALL deallocate_matrix_CSR_dist( mesh%M_ddz_bk_bks)

    ! bks to bk (for calculating (the vertical gradient of) the effective viscosity in the BPA)
    CALL deallocate_matrix_CSR_dist( mesh%M_map_bks_bk)
    CALL deallocate_matrix_CSR_dist( mesh%M_ddz_bks_bk)

    ! Map between the bks-grid and the ak-grid (for calculating strain rates in the BPA)
    CALL deallocate_matrix_CSR_dist( mesh%M_map_bks_ak)
    CALL deallocate_matrix_CSR_dist( mesh%M_map_ak_bks)

    ! bk to bk (for constructing the BPA stiffness matrix)
    CALL deallocate_matrix_CSR_dist( mesh%M2_ddx_bk_bk)
    CALL deallocate_matrix_CSR_dist( mesh%M2_ddy_bk_bk)
    CALL deallocate_matrix_CSR_dist( mesh%M2_d2dx2_bk_bk)
    CALL deallocate_matrix_CSR_dist( mesh%M2_d2dxdy_bk_bk)
    CALL deallocate_matrix_CSR_dist( mesh%M2_d2dy2_bk_bk)
    CALL deallocate_matrix_CSR_dist( mesh%M2_ddz_bk_bk)
    CALL deallocate_matrix_CSR_dist( mesh%M2_d2dz2_bk_bk)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE deallocate_mesh

END MODULE mesh_memory
