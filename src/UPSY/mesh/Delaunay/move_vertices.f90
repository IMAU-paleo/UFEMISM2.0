module move_vertices

  ! Move a vertex to a new position and update the Delaunay triangulation accordingly.

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use mesh_types, only: type_mesh
  use flip_triangles, only: initialise_Delaunay_check_stack, &
    add_triangle_pairs_around_vertex_to_Delaunay_check_stack, flip_triangles_until_Delaunay
  use mesh_utilities, only: update_triangle_circumcenter

  implicit none

  private

  public :: move_vertex

contains

  subroutine move_vertex( mesh, vi, p)
    ! Move vertex vi of the mesh to point p

    ! In/output variables:
    type(type_mesh),            intent(inout)     :: mesh
    integer,                    intent(in)        :: vi
    real(dp), dimension(2),     intent(in)        :: p

    ! Local variables:
    character(len=256), parameter                 :: routine_name = 'move_vertex'
    integer                                       :: iti, ti

    ! Add routine to path
    call init_routine( routine_name)

    ! Move the vertex
    mesh%V( vi,:) = p

    ! Update surrounding triangle circumcentres
    do iti = 1, mesh%niTri( vi)
      ti = mesh%iTri( vi,iti)
      call update_triangle_circumcenter( mesh, ti)
    end do

    ! Update triangulation
    call initialise_Delaunay_check_stack( mesh)
    call add_triangle_pairs_around_vertex_to_Delaunay_check_stack( mesh, vi)
    call flip_triangles_until_Delaunay( mesh)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine move_vertex

end module move_vertices
