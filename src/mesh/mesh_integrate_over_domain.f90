module mesh_integrate_over_domain

  ! Integrate a function defined on the mesh vertices/triangles/edges over the domain

  use precisions, only: dp
  use mesh_types, only: type_mesh
  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use mpi_f08, only: MPI_ALLREDUCE, MPI_IN_PLACE, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD

  implicit none

  private

  public :: integrate_over_domain_a, integrate_over_domain_b, integrate_over_domain_c, &
    average_over_domain_a, average_over_domain_b, average_over_domain_c

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

    int_d = 0._dp
    do vi = mesh%vi1, mesh%vi2
      int_d = int_d + d( vi) * mesh%A( vi)
    end do
    call MPI_ALLREDUCE( MPI_IN_PLACE, int_d, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine integrate_over_domain_a

  subroutine integrate_over_domain_b( mesh, d, int_d)
    !< Integrate a function defined on the mesh vertices over the domain

    ! In/output variables:
    type(type_mesh),                         intent(in   ) :: mesh
    real(dp), dimension( mesh%ti1:mesh%ti2), intent(in   ) :: d
    real(dp),                                intent(  out) :: int_d

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'integrate_over_domain_b'
    integer                        :: ti, ierr

    ! Add routine to path
    call init_routine( routine_name)

    int_d = 0._dp
    do ti = mesh%ti1, mesh%ti2
      int_d = int_d + d( ti) * mesh%TriA( ti)
    end do
    call MPI_ALLREDUCE( MPI_IN_PLACE, int_d, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine integrate_over_domain_b

  subroutine integrate_over_domain_c( mesh, d, int_d)
    !< Integrate a function defined on the mesh edges over the domain

    ! In/output variables:
    type(type_mesh),                         intent(in   ) :: mesh
    real(dp), dimension( mesh%ei1:mesh%ei2), intent(in   ) :: d
    real(dp),                                intent(  out) :: int_d

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'integrate_over_domain_c'
    integer                        :: ei, ierr

    ! Add routine to path
    call init_routine( routine_name)

    int_d = 0._dp
    do ei = mesh%ei1, mesh%ei2
      int_d = int_d + d( ei) * mesh%EA( ei)
    end do
    call MPI_ALLREDUCE( MPI_IN_PLACE, int_d, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine integrate_over_domain_c

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

  subroutine average_over_domain_b( mesh, d, av_d)
    !< Average a function defined on the mesh triangles over the domain

    ! In/output variables:
    type(type_mesh),                         intent(in   ) :: mesh
    real(dp), dimension( mesh%ti1:mesh%ti2), intent(in   ) :: d
    real(dp),                                intent(  out) :: av_d

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'average_over_domain_b'
    real(dp)                       :: int_d

    ! Add routine to path
    call init_routine( routine_name)

    call integrate_over_domain_b( mesh, d, int_d)
    av_d = int_d / ((mesh%xmax - mesh%xmin) * (mesh%ymax - mesh%ymin))

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine average_over_domain_b

  subroutine average_over_domain_c( mesh, d, av_d)
    !< Average a function defined on the mesh edges over the domain

    ! In/output variables:
    type(type_mesh),                         intent(in   ) :: mesh
    real(dp), dimension( mesh%ei1:mesh%ei2), intent(in   ) :: d
    real(dp),                                intent(  out) :: av_d

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'average_over_domain_c'
    real(dp)                       :: int_d

    ! Add routine to path
    call init_routine( routine_name)

    call integrate_over_domain_c( mesh, d, int_d)
    av_d = int_d / ((mesh%xmax - mesh%xmin) * (mesh%ymax - mesh%ymin))

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine average_over_domain_c

end module mesh_integrate_over_domain
