module ct_mass_conservation

  use precisions, only: dp
  use mpi_f08, only: MPI_COMM_WORLD, MPI_BCAST, MPI_CHAR
  use model_configuration, only: C
  use mpi_basic, only: par
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, colour_string
  use mesh_types, only: type_mesh
  use netcdf_io_main, only: open_existing_netcdf_file_for_reading, setup_mesh_from_file, &
    close_netcdf_file, create_new_netcdf_file_for_writing, setup_mesh_in_netcdf_file, &
    add_field_mesh_dp_2D_notime, write_to_field_multopt_mesh_dp_2D_notime
  use tests_main, only: test_mesh_is_self_consistent
  use assertions_basic, only: assert
  use Halfar_SIA_solution, only: Halfar
  use conservation_of_mass_main, only: calc_dHi_dt

  implicit none

  private

  public :: run_all_mass_cons_component_tests

contains

  !> Run all mass conservation component tests.
  subroutine run_all_mass_cons_component_tests( test_mesh_filenames)

    ! In/output variables:
    character(len=*), dimension(:), intent(in) :: test_mesh_filenames

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_all_mass_cons_component_tests'
    character(len=1024)            :: foldername_mass_cons
    integer                        :: i

    ! Add routine to call stack
    call init_routine( routine_name)

    if (par%primary) write(0,*) '  Running mass_cons component tests...'
    if (par%primary) write(0,*) ''

    call create_mass_cons_component_tests_output_folder( foldername_mass_cons)

    do i = 1, size( test_mesh_filenames,1)
      call run_mass_cons_test_on_mesh( foldername_mass_cons, test_mesh_filenames( i))
    end do

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_all_mass_cons_component_tests

  !> Create the output folder for the mass conservation component tests
  subroutine create_mass_cons_component_tests_output_folder( foldername_mass_cons)

    ! In/output variables:
    character(len=*), intent(out) :: foldername_mass_cons

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_mass_cons_component_tests_output_folder'
    logical                        :: ex
    integer                        :: ierr

    ! Add routine to path
    call init_routine( routine_name)

    foldername_mass_cons = trim(C%output_dir) // '/mass_conservation'

    if (par%primary) then

      ! Remove existing folder if necessary
      inquire( file = trim( foldername_mass_cons) // '/.', exist = ex)
      if (ex) then
        call system('rm -rf ' // trim( foldername_mass_cons))
      end if

      ! Create the directory
      call system('mkdir ' // trim( foldername_mass_cons))

    end if
    call MPI_BCAST( foldername_mass_cons, len(foldername_mass_cons), MPI_CHAR, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_mass_cons_component_tests_output_folder

  !> Run the mass conservation test on a particular mesh
  subroutine run_mass_cons_test_on_mesh( foldername_mass_cons, test_mesh_filename)

    ! In/output variables:
    character(len=*), intent(in) :: foldername_mass_cons
    character(len=*), intent(in) :: test_mesh_filename

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'run_mass_cons_test_on_mesh'
    type(type_mesh)                     :: mesh
    integer                             :: ncid
    real(dp)                            :: A, n, H0, R0, x, y, t
    integer                             :: vi, ti
    real(dp), dimension(:), allocatable :: Hi, Hb, Hs, SL
    real(dp), dimension(:), allocatable :: u_vav_b, v_vav_b, dHi_dt_ex
    real(dp), dimension(:), allocatable :: SMB, BMB, LMB, AMB
    real(dp), dimension(:), allocatable :: fraction_margin
    logical,  dimension(:), allocatable :: mask_noice
    real(dp)                            :: dt
    real(dp), dimension(:), allocatable :: Hi_tplusdt
    real(dp), dimension(:), allocatable :: divQ, dHi_dt_target
    real(dp), dimension(:), allocatable :: dHi_dt_expl, dHi_dt_semiimpl, dHi_dt_impl, dHi_dt_overimpl

    ! Add routine to call stack
    call init_routine( routine_name)

    if (par%primary) write(0,*) '      Running mass conservation test on mesh ', &
      colour_string(trim(test_mesh_filename( index( test_mesh_filename,'/',back=.true.)+1:&
      len_trim( test_mesh_filename))),'light blue'), '...'

    ! Set up the mesh from the file (includes calculating secondary geometry data and matrix operators)
    call open_existing_netcdf_file_for_reading( trim(test_mesh_filename), ncid)
    call setup_mesh_from_file( test_mesh_filename, ncid, mesh)
    call close_netcdf_file( ncid)

    ! Check mesh self-consistency
    call assert( test_mesh_is_self_consistent( mesh), 'mesh is not self-consistent')

    ! Allocate memory
    allocate( Hi             ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( Hb             ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( Hs             ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( SL             ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( u_vav_b        ( mesh%ti1:mesh%ti2), source = 0._dp)
    allocate( v_vav_b        ( mesh%ti1:mesh%ti2), source = 0._dp)
    allocate( dHi_dt_ex      ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( SMB            ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( BMB            ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( LMB            ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( AMB            ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( fraction_margin( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( mask_noice     ( mesh%vi1:mesh%vi2), source = .false.)
    allocate( Hi_tplusdt     ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( divQ           ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( dHi_dt_target  ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( dHi_dt_expl    ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( dHi_dt_semiimpl( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( dHi_dt_impl    ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( dHi_dt_overimpl( mesh%vi1:mesh%vi2), source = 0._dp)

    ! Calculate ice thickness, vertically averaged velocities and thinning rates from
    ! the extended Halfar solution

    A  = 1e-16_dp
    n  = 3._dp
    H0 = 6000._dp
    R0 = 1500e3_dp
    t  = 0._dp

    do vi = mesh%vi1, mesh%vi2
      x = mesh%V( vi,1)
      y = mesh%V( vi,2)

      Hi       ( vi) = Halfar%H    ( A, n, H0, R0, x, y, t)
      dHi_dt_ex( vi) = Halfar%dH_dt( A, n, H0, R0, x, y, t)

      Hb( vi) = 0._dp
      Hs( vi) = Hi( vi)
      SL( vi) = -100._dp
    end do

    do ti = mesh%ti1, mesh%ti2
      u_vav_b( ti) = Halfar%u_vav( A, n, H0, R0, mesh%Tricc( ti,1), mesh%Tricc( ti,2), t)
      v_vav_b( ti) = Halfar%v_vav( A, n, H0, R0, mesh%Tricc( ti,1), mesh%Tricc( ti,2), t)
    end do

    ! Calculate modelled thinning rates using different solvers
    SMB             = 0._dp
    BMB             = 0._dp
    LMB             = 0._dp
    AMB             = 0._dp
    fraction_margin = 1._dp
    mask_noice      = .false.
    dHi_dt_target   = 0._dp
    dt              = 0.1_dp

    ! Explicit
    C%choice_ice_integration_method = 'explicit'
    call calc_dHi_dt( mesh, Hi, Hb, SL, u_vav_b, v_vav_b, SMB, BMB, LMB, AMB, &
      fraction_margin, mask_noice, dt, dHi_dt_expl, Hi_tplusdt, divQ, dHi_dt_target)

    ! Semi-implicit
    C%choice_ice_integration_method = 'semi-implicit'
    C%dHi_semiimplicit_fs = 0.5_dp
    call calc_dHi_dt( mesh, Hi, Hb, SL, u_vav_b, v_vav_b, SMB, BMB, LMB, AMB, &
      fraction_margin, mask_noice, dt, dHi_dt_semiimpl, Hi_tplusdt, divQ, dHi_dt_target)

    ! Implicit
    C%choice_ice_integration_method = 'semi-implicit'
    C%dHi_semiimplicit_fs = 1._dp
    call calc_dHi_dt( mesh, Hi, Hb, SL, u_vav_b, v_vav_b, SMB, BMB, LMB, AMB, &
      fraction_margin, mask_noice, dt, dHi_dt_impl, Hi_tplusdt, divQ, dHi_dt_target)

    ! Over-implicit
    C%choice_ice_integration_method = 'semi-implicit'
    C%dHi_semiimplicit_fs = 1.5_dp
    call calc_dHi_dt( mesh, Hi, Hb, SL, u_vav_b, v_vav_b, SMB, BMB, LMB, AMB, &
      fraction_margin, mask_noice, dt, dHi_dt_overimpl, Hi_tplusdt, divQ, dHi_dt_target)

    ! Write results to output
    call write_mass_cons_test_results_to_file( &
      foldername_mass_cons, test_mesh_filename, mesh, dHi_dt_ex, &
      dHi_dt_expl, dHi_dt_semiimpl, dHi_dt_impl, dHi_dt_overimpl)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_mass_cons_test_on_mesh

  ! !> Write the results of mass conservation test= for a particular mesh to a file.
  subroutine write_mass_cons_test_results_to_file( &
    foldername_mass_cons, test_mesh_filename, mesh, dHi_dt_ex, &
    dHi_dt_expl, dHi_dt_semiimpl, dHi_dt_impl, dHi_dt_overimpl)

    ! In/output variables:
    character(len=*),                       intent(in) :: foldername_mass_cons
    character(len=*),                       intent(in) :: test_mesh_filename
    type(type_mesh),                        intent(in) :: mesh
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in) :: dHi_dt_ex
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in) :: dHi_dt_expl, dHi_dt_semiimpl, dHi_dt_impl, dHi_dt_overimpl

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'write_mass_cons_test_results_to_file'
    integer                        :: ncid
    character(len=1024)            :: mesh_name
    character(len=1024)            :: filename
    integer                        :: i

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Create a file and write the mesh to it
    mesh_name = test_mesh_filename( 1:len_trim( test_mesh_filename)-3)
    i = index( mesh_name, '/', back = .true.)
    mesh_name = mesh_name( i+1:len_trim( mesh_name))
    filename = trim( foldername_mass_cons) // '/res_' // &
      trim( mesh_name) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)

    ! ! Add all the variables
    call add_field_mesh_dp_2D_notime( filename, ncid, 'dHi_dt_ex')
    call add_field_mesh_dp_2D_notime( filename, ncid, 'dHi_dt_explicit')
    call add_field_mesh_dp_2D_notime( filename, ncid, 'dHi_dt_semiimplicit')
    call add_field_mesh_dp_2D_notime( filename, ncid, 'dHi_dt_implicit')
    call add_field_mesh_dp_2D_notime( filename, ncid, 'dHi_dt_overimplicit')

    ! Write all the variables
    call write_to_field_multopt_mesh_dp_2D_notime( mesh, filename, ncid, 'dHi_dt_ex'          , dHi_dt_ex)
    call write_to_field_multopt_mesh_dp_2D_notime( mesh, filename, ncid, 'dHi_dt_explicit'    , dHi_dt_expl)
    call write_to_field_multopt_mesh_dp_2D_notime( mesh, filename, ncid, 'dHi_dt_semiimplicit', dHi_dt_semiimpl)
    call write_to_field_multopt_mesh_dp_2D_notime( mesh, filename, ncid, 'dHi_dt_implicit'    , dHi_dt_impl)
    call write_to_field_multopt_mesh_dp_2D_notime( mesh, filename, ncid, 'dHi_dt_overimplicit', dHi_dt_overimpl)

    ! Close the file
    call close_netcdf_file( ncid)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine write_mass_cons_test_results_to_file


end module ct_mass_conservation
