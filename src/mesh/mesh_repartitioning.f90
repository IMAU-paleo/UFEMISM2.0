module mesh_repartitioning

  use precisions, only: dp
  use mesh_types, only: type_mesh
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, warning, crash
  use mpi_basic, only: par
  use mesh_memory, only: allocate_mesh_primary
  use mesh_secondary, only: calc_all_secondary_mesh_data
  use mesh_disc_calc_matrix_operators_2D, only: calc_all_matrix_operators_mesh
  use mpi_distributed_shared_memory, only: allocate_dist_shared, gather_dist_shared_to_all, &
    deallocate_dist_shared
  use mpi_f08, only: MPI_WIN

  implicit none

  private

  public :: repartition_mesh

contains

  subroutine repartition_mesh( mesh, mesh_new, mask_active_a_nih, mask_active_b_nih)
    !< Create a copy of the mesh with the vertices and triangles
    !< partitioned over the processes to balance the number of
    !< active vertices/triangles specified in the masks.

    ! In/output variables:
    type(type_mesh),                                             intent(in   ) :: mesh
    type(type_mesh),                                             intent(inout) :: mesh_new
    logical, dimension(mesh%pai_V%i1_nih:mesh%pai_V%i2_nih),     intent(in   ) :: mask_active_a_nih
    logical, dimension(mesh%pai_Tri%i1_nih:mesh%pai_Tri%i2_nih), intent(in   ) :: mask_active_b_nih

    ! Local variables:
    character(len=1024), parameter             :: routine_name = 'repartition_mesh'
    logical, dimension(:), contiguous, pointer :: mask_active_a_tot => null()
    logical, dimension(:), contiguous, pointer :: mask_active_b_tot => null()
    type(MPI_WIN)                              :: wmask_active_a_tot, wmask_active_b_tot

    ! Add routine to path
    call init_routine( routine_name)

    ! Copy mesh metadata
    mesh_new%nV       = mesh%nV
    mesh_new%nTri     = mesh%nTri

    mesh_new%xmin     = mesh%xmin
    mesh_new%xmax     = mesh%xmax
    mesh_new%ymin     = mesh%ymin
    mesh_new%ymax     = mesh%ymax

    mesh_new%tol_dist = mesh%tol_dist

    ! Copy primary mesh data

    call allocate_mesh_primary( mesh_new, trim(mesh%name) // '_repartitioned', &
      mesh%nV, mesh%nTri, mesh%nC_mem)

    mesh_new%V(                1:mesh%nV  ,:) = mesh%V(                1:mesh%nV  ,:)
    mesh_new%nC(               1:mesh%nV    ) = mesh%nC(               1:mesh%nV    )
    mesh_new%C(                1:mesh%nV  ,:) = mesh%C(                1:mesh%nV  ,:)
    mesh_new%niTri(            1:mesh%nV    ) = mesh%niTri(            1:mesh%nV    )
    mesh_new%iTri(             1:mesh%nV  ,:) = mesh%iTri(             1:mesh%nV  ,:)
    mesh_new%VBI(              1:mesh%nV    ) = mesh%VBI(              1:mesh%nV    )

    mesh_new%Tri(              1:mesh%nTri,:) = mesh%Tri(              1:mesh%nTri,:)
    mesh_new%Tricc(            1:mesh%nTri,:) = mesh%Tricc(            1:mesh%nTri,:)
    mesh_new%TriC(             1:mesh%nTri,:) = mesh%TriC(             1:mesh%nTri,:)

    mesh_new%Tri_flip_list(    1:mesh%nTri,:) = mesh%Tri_flip_list(    1:mesh%nTri,:)
    mesh_new%refinement_map(   1:mesh%nTri  ) = mesh%refinement_map(   1:mesh%nTri  )
    mesh_new%refinement_stack( 1:mesh%nTri  ) = mesh%refinement_stack( 1:mesh%nTri  )
    mesh_new%Tri_li(           1:mesh%nTri,:) = mesh%Tri_li(           1:mesh%nTri,:)

    ! Gather masks
    call allocate_dist_shared( mask_active_a_tot, wmask_active_a_tot, mesh%nV)
    call allocate_dist_shared( mask_active_b_tot, wmask_active_b_tot, mesh%nTri)
    call gather_dist_shared_to_all( mesh%pai_V,   mask_active_a_nih, mask_active_a_tot)
    call gather_dist_shared_to_all( mesh%pai_Tri, mask_active_b_nih, mask_active_b_tot)

    ! Simnly recalculate secondary mesh data and mesh operators
    ! (fast enough not to be a problem, and prevents code duplication)
    call calc_all_secondary_mesh_data( mesh_new, mesh%lambda_M, mesh%phi_M, mesh%beta_stereo, &
      mask_active_a_tot, mask_active_b_tot)
    call calc_all_matrix_operators_mesh( mesh_new)

    call deallocate_dist_shared( mask_active_a_tot, wmask_active_a_tot)
    call deallocate_dist_shared( mask_active_b_tot, wmask_active_b_tot)

    call crash('whoopsiedaisy')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine repartition_mesh

end module mesh_repartitioning
