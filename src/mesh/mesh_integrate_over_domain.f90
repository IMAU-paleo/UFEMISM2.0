module mesh_integrate_over_domain

  ! Integrate a function defined on the mesh vertices/triangles/edges over the domain

  use precisions, only: dp
  use mesh_types, only: type_mesh
  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use mpi_f08, only: MPI_ALLREDUCE, MPI_IN_PLACE, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD

  implicit none

  private

  public :: integrate_over_domain_a, average_over_domain_a

contains

  subroutine integrate_over_domain_a( mesh, d, int_d)
    !< Integrate a function defined on the mesh vertices over the domain

    ! In/output variables:
    type(type_mesh),                         intent(in   ) :: mesh
    real(dp), dimension( mesh%vi1:mesh%vi2), intent(in   ) :: d
    real(dp),                                intent(  out) :: int_d

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'integrate_over_domain_a'
    integer                        :: vi, ierr

    ! Add routine to path
    call init_routine( routine_name)

    ! Integrate over process domain
    int_d = 0._dp
    do vi = mesh%vi1, mesh%vi2
      int_d = int_d + d( vi) * mesh%A( vi)
    end do

    ! Reduce over processes to find total domain integral
    call MPI_ALLREDUCE( MPI_IN_PLACE, int_d, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine integrate_over_domain_a

  subroutine average_over_domain_a( mesh, d, av_d)
    !< Average a function defined on the mesh vertices over the domain

    ! In/output variables:
    type(type_mesh),                         intent(in   ) :: mesh
    real(dp), dimension( mesh%vi1:mesh%vi2), intent(in   ) :: d
    real(dp),                                intent(  out) :: av_d

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'average_over_domain_a'
    real(dp)                       :: int_d

    ! Add routine to path
    call init_routine( routine_name)

    call integrate_over_domain_a( mesh, d, int_d)
    av_d = int_d / ((mesh%xmax - mesh%xmin) * (mesh%ymax - mesh%ymin))

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine average_over_domain_a

end module mesh_integrate_over_domain
