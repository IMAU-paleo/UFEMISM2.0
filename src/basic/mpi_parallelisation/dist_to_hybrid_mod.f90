module dist_to_hybrid_mod

  use precisions, only: dp
  use mpi_basic, only: par, sync_node
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use parallel_array_info_type, only: type_par_arr_info

  implicit none

  private

  public :: dist_to_hybrid, hybrid_to_dist

  interface dist_to_hybrid
    !< Convert distributed to hybrid distributed/shared memory
    procedure :: dist_to_hybrid_logical_1D
    procedure :: dist_to_hybrid_logical_2D
    procedure :: dist_to_hybrid_logical_3D
    procedure :: dist_to_hybrid_int_1D
    procedure :: dist_to_hybrid_int_2D
    procedure :: dist_to_hybrid_int_3D
    procedure :: dist_to_hybrid_dp_1D
    procedure :: dist_to_hybrid_dp_2D
    procedure :: dist_to_hybrid_dp_3D
    procedure :: dist_to_hybrid_complex_1D
    procedure :: dist_to_hybrid_complex_2D
    procedure :: dist_to_hybrid_complex_3D
  end interface dist_to_hybrid

  interface hybrid_to_dist
    !< Convert hybrid distributed/shared to distributed memory
    procedure :: hybrid_to_dist_logical_1D
    procedure :: hybrid_to_dist_logical_2D
    procedure :: hybrid_to_dist_logical_3D
    procedure :: hybrid_to_dist_int_1D
    procedure :: hybrid_to_dist_int_2D
    procedure :: hybrid_to_dist_int_3D
    procedure :: hybrid_to_dist_dp_1D
    procedure :: hybrid_to_dist_dp_2D
    procedure :: hybrid_to_dist_dp_3D
    procedure :: hybrid_to_dist_complex_1D
    procedure :: hybrid_to_dist_complex_2D
    procedure :: hybrid_to_dist_complex_3D
  end interface hybrid_to_dist

contains

  ! == dist to hybrid

  subroutine dist_to_hybrid_logical_1D( pai, d, d_nih)

    ! In/output variables:
    type(type_par_arr_info),                   intent(in   ) :: pai
    logical, dimension(pai%i1:pai%i2),         intent(in   ) :: d
    logical, dimension(pai%i1_nih:pai%i2_nih), intent(  out) :: d_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'dist_to_hybrid_logical_1D'

    ! Add routine to path
    call init_routine( routine_name)

    d_nih( pai%i1:pai%i2) = d
    call sync_node

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine dist_to_hybrid_logical_1D

  subroutine dist_to_hybrid_logical_2D( pai, nz, d, d_nih)

    ! In/output variables:
    type(type_par_arr_info),                        intent(in   ) :: pai
    integer,                                        intent(in   ) :: nz
    logical, dimension(pai%i1:pai%i2,1:nz),         intent(in   ) :: d
    logical, dimension(pai%i1_nih:pai%i2_nih,1:nz), intent(  out) :: d_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'dist_to_hybrid_logical_2D'

    ! Add routine to path
    call init_routine( routine_name)

    d_nih( pai%i1:pai%i2,:) = d
    call sync_node

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine dist_to_hybrid_logical_2D

  subroutine dist_to_hybrid_logical_3D( pai, nz, nl, d, d_nih)

    ! In/output variables:
    type(type_par_arr_info),                             intent(in   ) :: pai
    integer,                                             intent(in   ) :: nz, nl
    logical, dimension(pai%i1:pai%i2,1:nz,1:nl),         intent(in   ) :: d
    logical, dimension(pai%i1_nih:pai%i2_nih,1:nz,1:nl), intent(  out) :: d_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'dist_to_hybrid_logical_3D'

    ! Add routine to path
    call init_routine( routine_name)

    d_nih( pai%i1:pai%i2,:,:) = d
    call sync_node

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine dist_to_hybrid_logical_3D

  subroutine dist_to_hybrid_int_1D( pai, d, d_nih)

    ! In/output variables:
    type(type_par_arr_info),                   intent(in   ) :: pai
    integer, dimension(pai%i1:pai%i2),         intent(in   ) :: d
    integer, dimension(pai%i1_nih:pai%i2_nih), intent(  out) :: d_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'dist_to_hybrid_int_1D'

    ! Add routine to path
    call init_routine( routine_name)

    d_nih( pai%i1:pai%i2) = d
    call sync_node

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine dist_to_hybrid_int_1D

  subroutine dist_to_hybrid_int_2D( pai, nz, d, d_nih)

    ! In/output variables:
    type(type_par_arr_info),                        intent(in   ) :: pai
    integer,                                        intent(in   ) :: nz
    integer, dimension(pai%i1:pai%i2,1:nz),         intent(in   ) :: d
    integer, dimension(pai%i1_nih:pai%i2_nih,1:nz), intent(  out) :: d_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'dist_to_hybrid_int_2D'

    ! Add routine to path
    call init_routine( routine_name)

    d_nih( pai%i1:pai%i2,:) = d
    call sync_node

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine dist_to_hybrid_int_2D

  subroutine dist_to_hybrid_int_3D( pai, nz, nl, d, d_nih)

    ! In/output variables:
    type(type_par_arr_info),                             intent(in   ) :: pai
    integer,                                             intent(in   ) :: nz, nl
    integer, dimension(pai%i1:pai%i2,1:nz,1:nl),         intent(in   ) :: d
    integer, dimension(pai%i1_nih:pai%i2_nih,1:nz,1:nl), intent(  out) :: d_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'dist_to_hybrid_int_3D'

    ! Add routine to path
    call init_routine( routine_name)

    d_nih( pai%i1:pai%i2,:,:) = d
    call sync_node

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine dist_to_hybrid_int_3D

  subroutine dist_to_hybrid_dp_1D( pai, d, d_nih)

    ! In/output variables:
    type(type_par_arr_info),                    intent(in   ) :: pai
    real(dp), dimension(pai%i1:pai%i2),         intent(in   ) :: d
    real(dp), dimension(pai%i1_nih:pai%i2_nih), intent(  out) :: d_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'dist_to_hybrid_dp_1D'

    ! Add routine to path
    call init_routine( routine_name)

    d_nih( pai%i1:pai%i2) = d
    call sync_node

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine dist_to_hybrid_dp_1D

  subroutine dist_to_hybrid_dp_2D( pai, nz, d, d_nih)

    ! In/output variables:
    type(type_par_arr_info),                         intent(in   ) :: pai
    integer,                                         intent(in   ) :: nz
    real(dp), dimension(pai%i1:pai%i2,1:nz),         intent(in   ) :: d
    real(dp), dimension(pai%i1_nih:pai%i2_nih,1:nz), intent(  out) :: d_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'dist_to_hybrid_dp_2D'

    ! Add routine to path
    call init_routine( routine_name)

    d_nih( pai%i1:pai%i2,:) = d
    call sync_node

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine dist_to_hybrid_dp_2D

  subroutine dist_to_hybrid_dp_3D( pai, nz, nl, d, d_nih)

    ! In/output variables:
    type(type_par_arr_info),                              intent(in   ) :: pai
    integer,                                              intent(in   ) :: nz, nl
    real(dp), dimension(pai%i1:pai%i2,1:nz,1:nl),         intent(in   ) :: d
    real(dp), dimension(pai%i1_nih:pai%i2_nih,1:nz,1:nl), intent(  out) :: d_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'dist_to_hybrid_dp_3D'

    ! Add routine to path
    call init_routine( routine_name)

    d_nih( pai%i1:pai%i2,:,:) = d
    call sync_node

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine dist_to_hybrid_dp_3D

  subroutine dist_to_hybrid_complex_1D( pai, d, d_nih)

    ! In/output variables:
    type(type_par_arr_info),                      intent(in   ) :: pai
    complex*16, dimension(pai%i1:pai%i2),         intent(in   ) :: d
    complex*16, dimension(pai%i1_nih:pai%i2_nih), intent(  out) :: d_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'dist_to_hybrid_complex_1D'

    ! Add routine to path
    call init_routine( routine_name)

    d_nih( pai%i1:pai%i2) = d
    call sync_node

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine dist_to_hybrid_complex_1D

  subroutine dist_to_hybrid_complex_2D( pai, nz, d, d_nih)

    ! In/output variables:
    type(type_par_arr_info),                           intent(in   ) :: pai
    integer,                                           intent(in   ) :: nz
    complex*16, dimension(pai%i1:pai%i2,1:nz),         intent(in   ) :: d
    complex*16, dimension(pai%i1_nih:pai%i2_nih,1:nz), intent(  out) :: d_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'dist_to_hybrid_complex_2D'

    ! Add routine to path
    call init_routine( routine_name)

    d_nih( pai%i1:pai%i2,:) = d
    call sync_node

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine dist_to_hybrid_complex_2D

  subroutine dist_to_hybrid_complex_3D( pai, nz, nl, d, d_nih)

    ! In/output variables:
    type(type_par_arr_info),                                intent(in   ) :: pai
    integer,                                                intent(in   ) :: nz, nl
    complex*16, dimension(pai%i1:pai%i2,1:nz,1:nl),         intent(in   ) :: d
    complex*16, dimension(pai%i1_nih:pai%i2_nih,1:nz,1:nl), intent(  out) :: d_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'dist_to_hybrid_complex_3D'

    ! Add routine to path
    call init_routine( routine_name)

    d_nih( pai%i1:pai%i2,:,:) = d
    call sync_node

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine dist_to_hybrid_complex_3D

  ! == hybrid to dist

  subroutine hybrid_to_dist_logical_1D( pai, d_nih, d)

    ! In/output variables:
    type(type_par_arr_info),                   intent(in   ) :: pai
    logical, dimension(pai%i1_nih:pai%i2_nih), intent(in   ) :: d_nih
    logical, dimension(pai%i1:pai%i2),         intent(  out) :: d

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'hybrid_to_dist_logical_1D'

    ! Add routine to path
    call init_routine( routine_name)

    d = d_nih( pai%i1:pai%i2)
    call sync_node

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine hybrid_to_dist_logical_1D

  subroutine hybrid_to_dist_logical_2D( pai, nz, d_nih, d)

    ! In/output variables:
    type(type_par_arr_info),                        intent(in   ) :: pai
    integer,                                        intent(in   ) :: nz
    logical, dimension(pai%i1_nih:pai%i2_nih,1:nz), intent(in   ) :: d_nih
    logical, dimension(pai%i1:pai%i2,1:nz),         intent(  out) :: d

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'hybrid_to_dist_logical_2D'

    ! Add routine to path
    call init_routine( routine_name)

    d = d_nih( pai%i1:pai%i2,:)
    call sync_node

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine hybrid_to_dist_logical_2D

  subroutine hybrid_to_dist_logical_3D( pai, nz, nl, d_nih, d)

    ! In/output variables:
    type(type_par_arr_info),                             intent(in   ) :: pai
    integer,                                             intent(in   ) :: nz, nl
    logical, dimension(pai%i1_nih:pai%i2_nih,1:nz,1:nl), intent(in   ) :: d_nih
    logical, dimension(pai%i1:pai%i2,1:nz,1:nl),         intent(  out) :: d

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'hybrid_to_dist_logical_3D'

    ! Add routine to path
    call init_routine( routine_name)

    d = d_nih( pai%i1:pai%i2,:,:)
    call sync_node

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine hybrid_to_dist_logical_3D

  subroutine hybrid_to_dist_int_1D( pai, d_nih, d)

    ! In/output variables:
    type(type_par_arr_info),                   intent(in   ) :: pai
    integer, dimension(pai%i1_nih:pai%i2_nih), intent(in   ) :: d_nih
    integer, dimension(pai%i1:pai%i2),         intent(  out) :: d

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'hybrid_to_dist_int_1D'

    ! Add routine to path
    call init_routine( routine_name)

    d = d_nih( pai%i1:pai%i2)
    call sync_node

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine hybrid_to_dist_int_1D

  subroutine hybrid_to_dist_int_2D( pai, nz, d_nih, d)

    ! In/output variables:
    type(type_par_arr_info),                        intent(in   ) :: pai
    integer,                                        intent(in   ) :: nz
    integer, dimension(pai%i1_nih:pai%i2_nih,1:nz), intent(in   ) :: d_nih
    integer, dimension(pai%i1:pai%i2,1:nz),         intent(  out) :: d

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'hybrid_to_dist_int_2D'

    ! Add routine to path
    call init_routine( routine_name)

    d = d_nih( pai%i1:pai%i2,:)
    call sync_node

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine hybrid_to_dist_int_2D

  subroutine hybrid_to_dist_int_3D( pai, nz, nl, d_nih, d)

    ! In/output variables:
    type(type_par_arr_info),                             intent(in   ) :: pai
    integer,                                             intent(in   ) :: nz, nl
    integer, dimension(pai%i1_nih:pai%i2_nih,1:nz,1:nl), intent(in   ) :: d_nih
    integer, dimension(pai%i1:pai%i2,1:nz,1:nl),         intent(  out) :: d

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'hybrid_to_dist_int_3D'

    ! Add routine to path
    call init_routine( routine_name)

    d = d_nih( pai%i1:pai%i2,:,:)
    call sync_node

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine hybrid_to_dist_int_3D

  subroutine hybrid_to_dist_dp_1D( pai, d_nih, d)

    ! In/output variables:
    type(type_par_arr_info),                    intent(in   ) :: pai
    real(dp), dimension(pai%i1_nih:pai%i2_nih), intent(in   ) :: d_nih
    real(dp), dimension(pai%i1:pai%i2),         intent(  out) :: d

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'hybrid_to_dist_dp_1D'

    ! Add routine to path
    call init_routine( routine_name)

    d = d_nih( pai%i1:pai%i2)
    call sync_node

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine hybrid_to_dist_dp_1D

  subroutine hybrid_to_dist_dp_2D( pai, nz, d_nih, d)

    ! In/output variables:
    type(type_par_arr_info),                         intent(in   ) :: pai
    integer,                                         intent(in   ) :: nz
    real(dp), dimension(pai%i1_nih:pai%i2_nih,1:nz), intent(in   ) :: d_nih
    real(dp), dimension(pai%i1:pai%i2,1:nz),         intent(  out) :: d

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'hybrid_to_dist_dp_2D'

    ! Add routine to path
    call init_routine( routine_name)

    d = d_nih( pai%i1:pai%i2,:)
    call sync_node

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine hybrid_to_dist_dp_2D

  subroutine hybrid_to_dist_dp_3D( pai, nz, nl, d_nih, d)

    ! In/output variables:
    type(type_par_arr_info),                              intent(in   ) :: pai
    integer,                                              intent(in   ) :: nz, nl
    real(dp), dimension(pai%i1_nih:pai%i2_nih,1:nz,1:nl), intent(in   ) :: d_nih
    real(dp), dimension(pai%i1:pai%i2,1:nz,1:nl),         intent(  out) :: d

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'hybrid_to_dist_dp_3D'

    ! Add routine to path
    call init_routine( routine_name)

    d = d_nih( pai%i1:pai%i2,:,:)
    call sync_node

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine hybrid_to_dist_dp_3D

  subroutine hybrid_to_dist_complex_1D( pai, d_nih, d)

    ! In/output variables:
    type(type_par_arr_info),                      intent(in   ) :: pai
    complex*16, dimension(pai%i1_nih:pai%i2_nih), intent(in   ) :: d_nih
    complex*16, dimension(pai%i1:pai%i2),         intent(  out) :: d

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'hybrid_to_dist_complex_1D'

    ! Add routine to path
    call init_routine( routine_name)

    d = d_nih( pai%i1:pai%i2)
    call sync_node

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine hybrid_to_dist_complex_1D

  subroutine hybrid_to_dist_complex_2D( pai, nz, d_nih, d)

    ! In/output variables:
    type(type_par_arr_info),                           intent(in   ) :: pai
    integer,                                           intent(in   ) :: nz
    complex*16, dimension(pai%i1_nih:pai%i2_nih,1:nz), intent(in   ) :: d_nih
    complex*16, dimension(pai%i1:pai%i2,1:nz),         intent(  out) :: d

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'hybrid_to_dist_complex_2D'

    ! Add routine to path
    call init_routine( routine_name)

    d = d_nih( pai%i1:pai%i2,:)
    call sync_node

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine hybrid_to_dist_complex_2D

  subroutine hybrid_to_dist_complex_3D( pai, nz, nl, d_nih, d)

    ! In/output variables:
    type(type_par_arr_info),                                intent(in   ) :: pai
    integer,                                                intent(in   ) :: nz, nl
    complex*16, dimension(pai%i1_nih:pai%i2_nih,1:nz,1:nl), intent(in   ) :: d_nih
    complex*16, dimension(pai%i1:pai%i2,1:nz,1:nl),         intent(  out) :: d

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'hybrid_to_dist_complex_3D'

    ! Add routine to path
    call init_routine( routine_name)

    d = d_nih( pai%i1:pai%i2,:,:)
    call sync_node

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine hybrid_to_dist_complex_3D

end module dist_to_hybrid_mod
