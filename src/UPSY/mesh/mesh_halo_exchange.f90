module mesh_halo_exchange

  ! See mesh_parallelisation for an explanation of the halos.

  use precisions, only: dp
  use mesh_types, only: type_mesh
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, warning, crash
  use mpi_basic, only: par, sync
  use halo_exchange_mod, only: basic_halo_exchange

  implicit none

  private

  public :: exchange_halos

  interface exchange_halos
    procedure :: exchange_halos_logical
    procedure :: exchange_halos_int
    procedure :: exchange_halos_int_3D
    procedure :: exchange_halos_dp
    procedure :: exchange_halos_dp_3D
  end interface exchange_halos

contains

  subroutine exchange_halos_logical( mesh, d_nih)

    ! In/output variables:
    type(type_mesh),       intent(in   ) :: mesh
    logical, dimension(:), intent(inout) :: d_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'exchange_halos_logical'

    ! Add routine to path
    call init_routine( routine_name)

    ! If running on one node, do nothing
    if (par%n_nodes == 1) then
      call finalise_routine( routine_name)
      return
    end if

    if (size( d_nih,1) == mesh%pai_V%n_nih) then
      call basic_halo_exchange( mesh%pai_V, d_nih)
    elseif (size( d_nih,1) == mesh%pai_Tri%n_nih) then
      call basic_halo_exchange( mesh%pai_Tri, d_nih)
    elseif (size( d_nih,1) == mesh%pai_E%n_nih) then
      call basic_halo_exchange( mesh%pai_E, d_nih)
    else
      call crash('unexpected array size')
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine exchange_halos_logical

  subroutine exchange_halos_int( mesh, d_nih)

    ! In/output variables:
    type(type_mesh),       intent(in   ) :: mesh
    integer, dimension(:), intent(inout) :: d_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'exchange_halos_int'

    ! Add routine to path
    call init_routine( routine_name)

    ! If running on one node, do nothing
    if (par%n_nodes == 1) then
      call finalise_routine( routine_name)
      return
    end if

    if (size( d_nih,1) == mesh%pai_V%n_nih) then
      call basic_halo_exchange( mesh%pai_V, d_nih)
    elseif (size( d_nih,1) == mesh%pai_Tri%n_nih) then
      call basic_halo_exchange( mesh%pai_Tri, d_nih)
    elseif (size( d_nih,1) == mesh%pai_E%n_nih) then
      call basic_halo_exchange( mesh%pai_E, d_nih)
    else
      call crash('unexpected array size')
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine exchange_halos_int

  subroutine exchange_halos_int_3D( mesh, d_nih)

    ! In/output variables:
    type(type_mesh),         intent(in   ) :: mesh
    integer, dimension(:,:), intent(inout) :: d_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'exchange_halos_int_3D'

    ! Add routine to path
    call init_routine( routine_name)

    ! If running on one node, do nothing
    if (par%n_nodes == 1) then
      call finalise_routine( routine_name)
      return
    end if

    if (size( d_nih,1) == mesh%pai_V%n_nih) then
      call basic_halo_exchange( mesh%pai_V, size( d_nih,2), d_nih)
    elseif (size( d_nih,1) == mesh%pai_Tri%n_nih) then
      call basic_halo_exchange( mesh%pai_Tri, size( d_nih,2), d_nih)
    elseif (size( d_nih,1) == mesh%pai_E%n_nih) then
      call basic_halo_exchange( mesh%pai_E, size( d_nih,2), d_nih)
    else
      call crash('unexpected array size')
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine exchange_halos_int_3D

  subroutine exchange_halos_dp( mesh, d_nih)

    ! In/output variables:
    type(type_mesh),        intent(in   ) :: mesh
    real(dp), dimension(:), intent(inout) :: d_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'exchange_halos_dp'

    ! Add routine to path
    call init_routine( routine_name)

    ! If running on one node, do nothing
    if (par%n_nodes == 1) then
      call finalise_routine( routine_name)
      return
    end if

    if (size( d_nih,1) == mesh%pai_V%n_nih) then
      call basic_halo_exchange( mesh%pai_V, d_nih)
    elseif (size( d_nih,1) == mesh%pai_Tri%n_nih) then
      call basic_halo_exchange( mesh%pai_Tri, d_nih)
    elseif (size( d_nih,1) == mesh%pai_E%n_nih) then
      call basic_halo_exchange( mesh%pai_E, d_nih)
    else
      call crash('unexpected array size')
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine exchange_halos_dp

  subroutine exchange_halos_dp_3D( mesh, d_nih)

    ! In/output variables:
    type(type_mesh),          intent(in   ) :: mesh
    real(dp), dimension(:,:), intent(inout) :: d_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'exchange_halos_dp_3D'

    ! Add routine to path
    call init_routine( routine_name)

    ! If running on one node, do nothing
    if (par%n_nodes == 1) then
      call finalise_routine( routine_name)
      return
    end if

    if (size( d_nih,1) == mesh%pai_V%n_nih) then
      call basic_halo_exchange( mesh%pai_V, size( d_nih,2), d_nih)
    elseif (size( d_nih,1) == mesh%pai_Tri%n_nih) then
      call basic_halo_exchange( mesh%pai_Tri, size( d_nih,2), d_nih)
    elseif (size( d_nih,1) == mesh%pai_E%n_nih) then
      call basic_halo_exchange( mesh%pai_E, size( d_nih,2), d_nih)
    else
      call crash('unexpected array size')
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine exchange_halos_dp_3D

end module mesh_halo_exchange
