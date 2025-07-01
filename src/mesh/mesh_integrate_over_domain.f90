module mesh_integrate_over_domain

  ! Integrate a function defined on the mesh vertices/triangles/edges over the domain

  use precisions, only: dp
  use mesh_types, only: type_mesh
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash, warning
  use mpi_f08, only: MPI_ALLREDUCE, MPI_IN_PLACE, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, &
    MPI_MIN, MPI_MAX
  use mpi_basic, only: par, sync

  implicit none

  private

  public :: integrate_over_domain, average_over_domain, calc_and_print_min_mean_max

contains

  subroutine integrate_over_domain( mesh, d, int_d, max_d, min_d)
    !< Integrate a function defined on the mesh over the domain

    ! In/output variables:
    type(type_mesh),                intent(in   ) :: mesh
    real(dp), dimension(:), target, intent(in   ) :: d
    real(dp),                       intent(  out) :: int_d
    real(dp),             optional, intent(  out) :: max_d, min_d

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'integrate_over_domain'
    real(dp), dimension(:), pointer :: d_nih, d_loc

    ! Add routine to path
    call init_routine( routine_name)

    if     (size( d,1) == mesh%pai_V%n_loc) then
      call integrate_over_domain_a( mesh, d, int_d, max_d, min_d)
    elseif (size( d,1) == mesh%pai_V%n_nih) then
      d_nih( mesh%pai_V%i1_nih:mesh%pai_V%i2_nih) => d
      d_loc => d_nih( mesh%pai_V%i1:mesh%pai_V%i2)
      call integrate_over_domain_a( mesh, d_loc, int_d, max_d, min_d)
    elseif (size( d,1) == mesh%pai_Tri%n_loc) then
      call integrate_over_domain_b( mesh, d, int_d, max_d, min_d)
    elseif (size( d,1) == mesh%pai_Tri%n_nih) then
      d_nih( mesh%pai_Tri%i1_nih:mesh%pai_Tri%i2_nih) => d
      d_loc => d_nih( mesh%pai_Tri%i1:mesh%pai_Tri%i2)
      call integrate_over_domain_b( mesh, d_loc, int_d, max_d, min_d)
    elseif (size( d,1) == mesh%pai_E%n_loc) then
      call integrate_over_domain_c( mesh, d, int_d, max_d, min_d)
    elseif (size( d,1) == mesh%pai_E%n_nih) then
      d_nih( mesh%pai_E%i1_nih:mesh%pai_E%i2_nih) => d
      d_loc => d_nih( mesh%pai_E%i1:mesh%pai_E%i2)
      call integrate_over_domain_c( mesh, d_loc, int_d, max_d, min_d)
    else
      call crash('invalid vector size')
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine integrate_over_domain

  subroutine average_over_domain( mesh, d, av_d)
    !< Average a function defined on the mesh vertices over the domain

    ! In/output variables:
    type(type_mesh),        intent(in   ) :: mesh
    real(dp), dimension(:), intent(in   ) :: d
    real(dp),               intent(  out) :: av_d

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'average_over_domain'
    real(dp)                       :: int_d

    ! Add routine to path
    call init_routine( routine_name)

    call integrate_over_domain( mesh, d, int_d)
    av_d = int_d / ((mesh%xmax - mesh%xmin) * (mesh%ymax - mesh%ymin))

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine average_over_domain

  subroutine calc_and_print_min_mean_max( mesh, d, name)

    ! In/output variables:
    type(type_mesh),                intent(in   ) :: mesh
    real(dp), dimension(:), target, intent(in   ) :: d
    character(len=*),               intent(in   ) :: name

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'calc_and_print_min_mean_max'
    real(dp), dimension(:), pointer :: d_nih, d_loc
    real(dp)                        :: d_min, d_max, d_av
    integer                         :: ierr
    character(len=20)               :: name_

    ! Add routine to path
    call init_routine( routine_name)

    call sync

    if     (size( d,1) == mesh%pai_V%n_loc) then
      d_loc => d
    elseif (size( d,1) == mesh%pai_V%n_nih) then
      d_nih( mesh%pai_V%i1_nih:mesh%pai_V%i2_nih) => d
      d_loc => d_nih( mesh%pai_V%i1:mesh%pai_V%i2)
    elseif (size( d,1) == mesh%pai_Tri%n_loc) then
      d_loc => d
    elseif (size( d,1) == mesh%pai_Tri%n_nih) then
      d_nih( mesh%pai_Tri%i1_nih:mesh%pai_Tri%i2_nih) => d
      d_loc => d_nih( mesh%pai_Tri%i1:mesh%pai_Tri%i2)
    elseif (size( d,1) == mesh%pai_E%n_loc) then
      d_loc => d
    elseif (size( d,1) == mesh%pai_E%n_nih) then
      d_nih( mesh%pai_E%i1_nih:mesh%pai_E%i2_nih) => d
      d_loc => d_nih( mesh%pai_E%i1:mesh%pai_E%i2)
    else
      call crash('invalid vector size')
    end if

    d_min = minval( d_loc)
    d_max = maxval( d_loc)
    call MPI_ALLREDUCE( MPI_IN_PLACE, d_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, d_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

    call average_over_domain( mesh, d_loc, d_av)

    name_ = ' '
    name_( len( name_)-len_trim( name)+1:len( name_)) = name( 1:len_trim( name))

    if (par%primary) call warning( name_ // ': [{dp_01} - {dp_02} - {dp_03}]', &
      dp_01 = d_min, dp_02 = d_av, dp_03 = d_max)

    nullify( d_nih)
    nullify( d_loc)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_and_print_min_mean_max

  subroutine integrate_over_domain_a( mesh, d, int_d, max_d, min_d)
    !< Integrate a function defined on the mesh vertices over the domain

    ! In/output variables:
    type(type_mesh),                         intent(in   ) :: mesh
    real(dp), dimension( mesh%vi1:mesh%vi2), intent(in   ) :: d
    real(dp),                                intent(  out) :: int_d
    real(dp),                      optional, intent(  out) :: max_d, min_d

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

    if (present( max_d)) then
      max_d = maxval( d)
      call MPI_ALLREDUCE( MPI_IN_PLACE, max_d, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    end if

    if (present( min_d)) then
      min_d = minval( d)
      call MPI_ALLREDUCE( MPI_IN_PLACE, min_d, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine integrate_over_domain_a

  subroutine integrate_over_domain_b( mesh, d, int_d, max_d, min_d)
    !< Integrate a function defined on the mesh vertices over the domain

    ! In/output variables:
    type(type_mesh),                         intent(in   ) :: mesh
    real(dp), dimension( mesh%ti1:mesh%ti2), intent(in   ) :: d
    real(dp),                                intent(  out) :: int_d
    real(dp),                      optional, intent(  out) :: max_d, min_d

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

    if (present( max_d)) then
      max_d = maxval( d)
      call MPI_ALLREDUCE( MPI_IN_PLACE, max_d, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    end if

    if (present( min_d)) then
      min_d = minval( d)
      call MPI_ALLREDUCE( MPI_IN_PLACE, min_d, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine integrate_over_domain_b

  subroutine integrate_over_domain_c( mesh, d, int_d, max_d, min_d)
    !< Integrate a function defined on the mesh edges over the domain

    ! In/output variables:
    type(type_mesh),                         intent(in   ) :: mesh
    real(dp), dimension( mesh%ei1:mesh%ei2), intent(in   ) :: d
    real(dp),                                intent(  out) :: int_d
    real(dp),                      optional, intent(  out) :: max_d, min_d

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

    if (present( max_d)) then
      max_d = maxval( d)
      call MPI_ALLREDUCE( MPI_IN_PLACE, max_d, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    end if

    if (present( min_d)) then
      min_d = minval( d)
      call MPI_ALLREDUCE( MPI_IN_PLACE, min_d, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine integrate_over_domain_c

end module mesh_integrate_over_domain
