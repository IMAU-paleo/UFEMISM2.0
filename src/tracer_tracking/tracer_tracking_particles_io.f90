module tracer_tracking_model_particles_io

  use tests_main
  use assertions_basic
  use precisions, only: dp
  use mpi_basic, only: par
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash, warning
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model
  use tracer_tracking_model_types, only: type_tracer_tracking_model_particles
  use model_configuration, only: C
  use mesh_utilities, only: find_containing_triangle, find_containing_vertex, &
    interpolate_to_point_dp_2D, interpolate_to_point_dp_3D
  use reallocate_mod, only: reallocate
  use netcdf, only: NF90_UNLIMITED, NF90_INT, NF90_DOUBLE
  use netcdf_io_main
  use mpi_distributed_memory, only: gather_to_primary

  implicit none

  private

  public :: create_particles_netcdf_file, write_to_particles_netcdf_file

contains

  subroutine write_to_particles_netcdf_file( particles, time)
    !< Write particle data to NetCDF

    ! In/output variables:
    type(type_tracer_tracking_model_particles), intent(in) :: particles
    real(dp),                                   intent(in) :: time

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'write_to_particles_netcdf_file'
    character(len=1024)                     :: filename
    integer                                 :: ncid, id_dim_time, ti
    integer,  dimension(:    ), allocatable :: id_tot
    real(dp), dimension(:,:  ), allocatable :: r_tot
    real(dp), dimension(:    ), allocatable :: t_origin_tot
    integer,  dimension(:,:  ), allocatable :: id_tot_with_time
    real(dp), dimension(:,:,:), allocatable :: r_tot_with_time
    real(dp), dimension(:,:  ), allocatable :: t_tot_origin_with_time

    ! Add routine to path
    call init_routine( routine_name)

    filename = particles%nc%filename

    call open_existing_netcdf_file_for_writing( filename, ncid)

    ! Write time and find timeframe
    call write_time_to_file( filename, ncid, time)
    call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Gather all particles to the primary
    if (par%primary) then
      allocate( id_tot      ( particles%n_max * par%n   ))
      allocate( r_tot       ( particles%n_max * par%n, 3))
      allocate( t_origin_tot( particles%n_max * par%n   ))
    else
      allocate( id_tot      ( 0   ))
      allocate( r_tot       ( 0, 0))
      allocate( t_origin_tot( 0   ))
    end if

    call gather_to_primary( particles%id      , id_tot)
    call gather_to_primary( particles%r       , r_tot)
    call gather_to_primary( particles%t_origin, t_origin_tot)

    ! Add "pretend" time dimension
    if (par%primary) then

      allocate( id_tot_with_time      ( particles%n_max * par%n,    1))
      allocate( r_tot_with_time       ( particles%n_max * par%n, 3, 1))
      allocate( t_tot_origin_with_time( particles%n_max * par%n,    1))

      id_tot_with_time      ( :  ,1) = id_tot
      r_tot_with_time       ( :,:,1) = r_tot
      t_tot_origin_with_time( :  ,1) = t_origin_tot

    end if

    ! Write data
    call write_var_primary( filename, ncid, particles%nc%id_var_id       , id_tot_with_time      , &
      start = [1,ti], count = [particles%n_max * par%n,1])
    call write_var_primary( filename, ncid, particles%nc%id_var_r        , r_tot_with_time       , &
      start = [1,1,ti], count = [particles%n_max * par%n,3,1])
    call write_var_primary( filename, ncid, particles%nc%id_var_t_origin , t_tot_origin_with_time, &
      start = [1,ti], count = [particles%n_max * par%n,1])

    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_particles_netcdf_file

  subroutine create_particles_netcdf_file( filename, particles)
    !< Create a NetCDF output file for the raw particle data

    ! In/output variables:
    character(len=*),                           intent(in   ) :: filename
    type(type_tracer_tracking_model_particles), intent(inout) :: particles

    ! Local variables
    character(len=1024), parameter :: routine_name = 'create_particles_netcdf_file'
    integer                        :: ncid, n, three, time

    ! Add routine to path
    call init_routine( routine_name)

    ! Create and open new file
    call create_new_netcdf_file_for_writing( filename, ncid)

    particles%nc%filename = trim(filename)

    ! Define dimensions
    call create_dimension( filename, ncid, 'n',     particles%n_max * par%n, particles%nc%id_dim_n    )
    call create_dimension( filename, ncid, 'three', 3                      , particles%nc%id_dim_three)
    call create_dimension( filename, ncid, 'time',  NF90_UNLIMITED         , particles%nc%id_dim_time )

    ! Abbreviations for shorter code
    n     = particles%nc%id_dim_n
    three = particles%nc%id_dim_three
    time  = particles%nc%id_dim_time

    ! Define variables
    call create_variable( filename, ncid, 'time', &
      NF90_DOUBLE, [time], particles%nc%id_var_time)
    call create_variable( filename, ncid, 'id', &
      NF90_INT, [n, time], particles%nc%id_var_id)
    call create_variable( filename, ncid, 'r', &
      NF90_DOUBLE, [n, three, time], particles%nc%id_var_r)
    call create_variable( filename, ncid, 't_origin', &
      NF90_DOUBLE, [n, time], particles%nc%id_var_t_origin)

    ! Close file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_particles_netcdf_file

end module tracer_tracking_model_particles_io