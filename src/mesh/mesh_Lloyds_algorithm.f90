module mesh_Lloyds_algorithm

  ! Lloyd's algorithm for "smoothing" a mesh

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use mesh_types, only: type_mesh
  use plane_geometry, only: cross2
  use move_vertices, only: move_vertex
  use mesh_refinement_basic, only: refine_mesh_split_encroaching_triangles_all

  implicit none

contains

subroutine Lloyds_algorithm_single_iteration( mesh, alpha_min)
  ! Lloyd's algorithm: move all vertices to the geometric centers of their Voronoi cells, and update the triangulation.
  ! This "smooths" the mesh, reducing resolution gradients and widening internal angles, thus making it more
  ! suitable for numerical methods (particularly the SSA).

  ! In/output variables:
  type(type_mesh),            intent(inout)     :: mesh
  real(dp),                   intent(in)        :: alpha_min     ! minimum allowed internal triangle angle

  ! Local variables:
  character(len=256), parameter                 :: routine_name = 'Lloyds_algorithm_single_iteration'
  integer                                       :: vi, ci, cip1
  real(dp)                                      :: VorTriA, sumVorTriA
  real(dp), dimension(2)                        :: pa, pb, pc, VorGC

  ! Add routine to path
  call init_routine( routine_name)

  ! Move all non-border vertices to their Voronoi cell geometric centre
  do vi = 1, mesh%nV

    ! Leave border vertices where they are
    if (mesh%VBI( vi) > 0) cycle

    ! Find the geometric centre of this vertex' Voronoi cell
    VorGC      = 0._dp
    sumVorTriA = 0._dp

    do ci = 1, mesh%nC( vi)

      cip1 = ci + 1
      if (cip1 > mesh%nC( vi)) cip1 = 1

      pa = mesh%V( vi,:)
      pb = mesh%V( mesh%C( vi,ci  ),:)
      pc = mesh%V( mesh%C( vi,cip1),:)

      VorTriA = cross2( pb - pa, pc - pa)

      VorGC = VorGC + VorTriA * (pa + pb + pc) / 3._dp
      sumVorTriA = sumVorTriA + VorTriA

    end do ! do ci = 1, mesh%nC( vi)

    VorGC = VorGC / sumVorTriA

    ! Move the vertex
    call move_vertex( mesh, vi, VorGC)

  end do

  ! Final step to ensure a nice clean mesh
  call refine_mesh_split_encroaching_triangles_all( mesh, alpha_min)

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine Lloyds_algorithm_single_iteration

end module mesh_Lloyds_algorithm
