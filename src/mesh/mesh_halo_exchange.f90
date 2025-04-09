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

    if (size( d_nih,1) == mesh%nV_nih) then
      call exchange_halos_a_logical( mesh, d_nih)
    elseif (size( d_nih,1) == mesh%nTri_nih) then
      call exchange_halos_b_logical( mesh, d_nih)
    elseif (size( d_nih,1) == mesh%nE_nih) then
      call exchange_halos_c_logical( mesh, d_nih)
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

    if (size( d_nih,1) == mesh%nV_nih) then
      call exchange_halos_a_int( mesh, d_nih)
    elseif (size( d_nih,1) == mesh%nTri_nih) then
      call exchange_halos_b_int( mesh, d_nih)
    elseif (size( d_nih,1) == mesh%nE_nih) then
      call exchange_halos_c_int( mesh, d_nih)
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

    if (size( d_nih,1) == mesh%nV_nih) then
      call exchange_halos_a_int_3D( mesh, d_nih)
    elseif (size( d_nih,1) == mesh%nTri_nih) then
      call exchange_halos_b_int_3D( mesh, d_nih)
    elseif (size( d_nih,1) == mesh%nE_nih) then
      call exchange_halos_c_int_3D( mesh, d_nih)
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

    if (size( d_nih,1) == mesh%nV_nih) then
      call exchange_halos_a_dp( mesh, d_nih)
    elseif (size( d_nih,1) == mesh%nTri_nih) then
      call exchange_halos_b_dp( mesh, d_nih)
    elseif (size( d_nih,1) == mesh%nE_nih) then
      call exchange_halos_c_dp( mesh, d_nih)
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

    if (size( d_nih,1) == mesh%nV_nih) then
      call exchange_halos_a_dp_3D( mesh, d_nih)
    elseif (size( d_nih,1) == mesh%nTri_nih) then
      call exchange_halos_b_dp_3D( mesh, d_nih)
    elseif (size( d_nih,1) == mesh%nE_nih) then
      call exchange_halos_c_dp_3D( mesh, d_nih)
    else
      call crash('unexpected array size')
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine exchange_halos_dp_3D

  ! == a-grid

  subroutine exchange_halos_a_logical( mesh, d_a_nih)
    !< Exhchange halos for a field defined on the vertices

    ! In/output variables:
    type(type_mesh),       intent(in   ) :: mesh
    logical, dimension(:), intent(inout) :: d_a_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'exchange_halos_a_logical'

    ! Add routine to path
    call init_routine( routine_name)

    call basic_halo_exchange( d_a_nih, mesh%vi1_nih, mesh%vi2_nih, &
      mesh%vi1_hle, mesh%vi2_hle, mesh%vi1_hli, mesh%vi2_hli, &
      mesh%vi1_hre, mesh%vi2_hre, mesh%vi1_hri, mesh%vi2_hri)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine exchange_halos_a_logical

  subroutine exchange_halos_a_int( mesh, d_a_nih)
    !< Exhchange halos for a field defined on the vertices

    ! In/output variables:
    type(type_mesh),       intent(in   ) :: mesh
    integer, dimension(:), intent(inout) :: d_a_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'exchange_halos_a_int'

    ! Add routine to path
    call init_routine( routine_name)

    call basic_halo_exchange( d_a_nih, mesh%vi1_nih, mesh%vi2_nih, &
      mesh%vi1_hle, mesh%vi2_hle, mesh%vi1_hli, mesh%vi2_hli, &
      mesh%vi1_hre, mesh%vi2_hre, mesh%vi1_hri, mesh%vi2_hri)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine exchange_halos_a_int

  subroutine exchange_halos_a_int_3D( mesh, d_a_nih)
    !< Exhchange halos for a field defined on the vertices

    ! In/output variables:
    type(type_mesh),         intent(in   ) :: mesh
    integer, dimension(:,:), intent(inout) :: d_a_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'exchange_halos_a_int_3D'

    ! Add routine to path
    call init_routine( routine_name)

    call basic_halo_exchange( d_a_nih, mesh%vi1_nih, mesh%vi2_nih, &
      mesh%vi1_hle, mesh%vi2_hle, mesh%vi1_hli, mesh%vi2_hli, &
      mesh%vi1_hre, mesh%vi2_hre, mesh%vi1_hri, mesh%vi2_hri, size( d_a_nih,2))

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine exchange_halos_a_int_3D

  subroutine exchange_halos_a_dp( mesh, d_a_nih)
    !< Exhchange halos for a field defined on the vertices

    ! In/output variables:
    type(type_mesh),        intent(in   ) :: mesh
    real(dp), dimension(:), intent(inout) :: d_a_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'exchange_halos_a_dp'

    ! Add routine to path
    call init_routine( routine_name)

    call basic_halo_exchange( d_a_nih, mesh%vi1_nih, mesh%vi2_nih, &
      mesh%vi1_hle, mesh%vi2_hle, mesh%vi1_hli, mesh%vi2_hli, &
      mesh%vi1_hre, mesh%vi2_hre, mesh%vi1_hri, mesh%vi2_hri)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine exchange_halos_a_dp

  subroutine exchange_halos_a_dp_3D( mesh, d_a_nih)
    !< Exhchange halos for a field defined on the vertices

    ! In/output variables:
    type(type_mesh),          intent(in   ) :: mesh
    real(dp), dimension(:,:), intent(inout) :: d_a_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'exchange_halos_a_dp_3D'

    ! Add routine to path
    call init_routine( routine_name)

    call basic_halo_exchange( d_a_nih, mesh%vi1_nih, mesh%vi2_nih, &
      mesh%vi1_hle, mesh%vi2_hle, mesh%vi1_hli, mesh%vi2_hli, &
      mesh%vi1_hre, mesh%vi2_hre, mesh%vi1_hri, mesh%vi2_hri, size( d_a_nih,2))

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine exchange_halos_a_dp_3D

  ! == b-grid

  subroutine exchange_halos_b_logical( mesh, d_b_nih)
    !< Exhchange halos for a field defined on the triangles

    ! In/output variables:
    type(type_mesh),       intent(in   ) :: mesh
    logical, dimension(:), intent(inout) :: d_b_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'exchange_halos_b_logical'

    ! Add routine to path
    call init_routine( routine_name)

    call basic_halo_exchange( d_b_nih, mesh%ti1_nih, mesh%ti2_nih, &
      mesh%ti1_hle, mesh%ti2_hle, mesh%ti1_hli, mesh%ti2_hli, &
      mesh%ti1_hre, mesh%ti2_hre, mesh%ti1_hri, mesh%ti2_hri)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine exchange_halos_b_logical

  subroutine exchange_halos_b_int( mesh, d_b_nih)
    !< Exhchange halos for a field defined on the triangles

    ! In/output variables:
    type(type_mesh),       intent(in   ) :: mesh
    integer, dimension(:), intent(inout) :: d_b_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'exchange_halos_b_int'

    ! Add routine to path
    call init_routine( routine_name)

    call basic_halo_exchange( d_b_nih, mesh%ti1_nih, mesh%ti2_nih, &
      mesh%ti1_hle, mesh%ti2_hle, mesh%ti1_hli, mesh%ti2_hli, &
      mesh%ti1_hre, mesh%ti2_hre, mesh%ti1_hri, mesh%ti2_hri)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine exchange_halos_b_int

  subroutine exchange_halos_b_int_3D( mesh, d_b_nih)
    !< Exhchange halos for a field defined on the triangles

    ! In/output variables:
    type(type_mesh),         intent(in   ) :: mesh
    integer, dimension(:,:), intent(inout) :: d_b_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'exchange_halos_b_int_3D'

    ! Add routine to path
    call init_routine( routine_name)

    call basic_halo_exchange( d_b_nih, mesh%ti1_nih, mesh%ti2_nih, &
      mesh%ti1_hle, mesh%ti2_hle, mesh%ti1_hli, mesh%ti2_hli, &
      mesh%ti1_hre, mesh%ti2_hre, mesh%ti1_hri, mesh%ti2_hri, size( d_b_nih,2))

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine exchange_halos_b_int_3D

  subroutine exchange_halos_b_dp( mesh, d_b_nih)
    !< Exhchange halos for a field defined on the triangles

    ! In/output variables:
    type(type_mesh),        intent(in   ) :: mesh
    real(dp), dimension(:), intent(inout) :: d_b_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'exchange_halos_b_dp'

    ! Add routine to path
    call init_routine( routine_name)

    call basic_halo_exchange( d_b_nih, mesh%ti1_nih, mesh%ti2_nih, &
      mesh%ti1_hle, mesh%ti2_hle, mesh%ti1_hli, mesh%ti2_hli, &
      mesh%ti1_hre, mesh%ti2_hre, mesh%ti1_hri, mesh%ti2_hri)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine exchange_halos_b_dp

  subroutine exchange_halos_b_dp_3D( mesh, d_b_nih)
    !< Exhchange halos for a field defined on the triangles

    ! In/output variables:
    type(type_mesh),          intent(in   ) :: mesh
    real(dp), dimension(:,:), intent(inout) :: d_b_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'exchange_halos_b_dp_3D'

    ! Add routine to path
    call init_routine( routine_name)

    call basic_halo_exchange( d_b_nih, mesh%ti1_nih, mesh%ti2_nih, &
      mesh%ti1_hle, mesh%ti2_hle, mesh%ti1_hli, mesh%ti2_hli, &
      mesh%ti1_hre, mesh%ti2_hre, mesh%ti1_hri, mesh%ti2_hri, size( d_b_nih,2))

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine exchange_halos_b_dp_3D

  ! == c-grid

  subroutine exchange_halos_c_logical( mesh, d_c_nih)
    !< Exhchange halos for a field defined on the edges

    ! In/output variables:
    type(type_mesh),       intent(in   ) :: mesh
    logical, dimension(:), intent(inout) :: d_c_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'exchange_halos_c_logical'

    ! Add routine to path
    call init_routine( routine_name)

    call basic_halo_exchange( d_c_nih, mesh%ei1_nih, mesh%ei2_nih, &
      mesh%ei1_hle, mesh%ei2_hle, mesh%ei1_hli, mesh%ei2_hli, &
      mesh%ei1_hre, mesh%ei2_hre, mesh%ei1_hri, mesh%ei2_hri)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine exchange_halos_c_logical

  subroutine exchange_halos_c_int( mesh, d_c_nih)
    !< Exhchange halos for a field defined on the edges

    ! In/output variables:
    type(type_mesh),       intent(in   ) :: mesh
    integer, dimension(:), intent(inout) :: d_c_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'exchange_halos_c_int'

    ! Add routine to path
    call init_routine( routine_name)

    call basic_halo_exchange( d_c_nih, mesh%ei1_nih, mesh%ei2_nih, &
      mesh%ei1_hle, mesh%ei2_hle, mesh%ei1_hli, mesh%ei2_hli, &
      mesh%ei1_hre, mesh%ei2_hre, mesh%ei1_hri, mesh%ei2_hri)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine exchange_halos_c_int

  subroutine exchange_halos_c_int_3D( mesh, d_c_nih)
    !< Exhchange halos for a field defined on the edges

    ! In/output variables:
    type(type_mesh),         intent(in   ) :: mesh
    integer, dimension(:,:), intent(inout) :: d_c_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'exchange_halos_c_int_3D'

    ! Add routine to path
    call init_routine( routine_name)

    call basic_halo_exchange( d_c_nih, mesh%ei1_nih, mesh%ei2_nih, &
      mesh%ei1_hle, mesh%ei2_hle, mesh%ei1_hli, mesh%ei2_hli, &
      mesh%ei1_hre, mesh%ei2_hre, mesh%ei1_hri, mesh%ei2_hri, size( d_c_nih,2))

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine exchange_halos_c_int_3D

  subroutine exchange_halos_c_dp( mesh, d_c_nih)
    !< Exhchange halos for a field defined on the edges

    ! In/output variables:
    type(type_mesh),        intent(in   ) :: mesh
    real(dp), dimension(:), intent(inout) :: d_c_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'exchange_halos_c_dp'

    ! Add routine to path
    call init_routine( routine_name)

    call basic_halo_exchange( d_c_nih, mesh%ei1_nih, mesh%ei2_nih, &
      mesh%ei1_hle, mesh%ei2_hle, mesh%ei1_hli, mesh%ei2_hli, &
      mesh%ei1_hre, mesh%ei2_hre, mesh%ei1_hri, mesh%ei2_hri)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine exchange_halos_c_dp

  subroutine exchange_halos_c_dp_3D( mesh, d_c_nih)
    !< Exhchange halos for a field defined on the edges

    ! In/output variables:
    type(type_mesh),          intent(in   ) :: mesh
    real(dp), dimension(:,:), intent(inout) :: d_c_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'exchange_halos_c_dp_3D'

    ! Add routine to path
    call init_routine( routine_name)

    call basic_halo_exchange( d_c_nih, mesh%ei1_nih, mesh%ei2_nih, &
      mesh%ei1_hle, mesh%ei2_hle, mesh%ei1_hli, mesh%ei2_hli, &
      mesh%ei1_hre, mesh%ei2_hre, mesh%ei1_hri, mesh%ei2_hri, size( d_c_nih,2))

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine exchange_halos_c_dp_3D

end module mesh_halo_exchange
