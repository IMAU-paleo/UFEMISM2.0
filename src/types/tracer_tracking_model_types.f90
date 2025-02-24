module tracer_tracking_model_types

  use precisions, only: dp

  implicit none

  private

  public :: type_tracer_tracking_model, type_tracer_tracking_model_particles, &
    type_map_particles_to_mesh

  type type_map_particles_to_mesh
    !< Indices and interpolation weights to map tracers from the particles to the model mesh

    integer                                 :: n    !< Number of nearest particles to interpolate between
    integer,  dimension(:,:,:), allocatable :: ip   !< Indices of   the n nearest particles
    real(dp), dimension(:,:,:), allocatable :: d    !< Distances to the n nearest particles

  end type type_map_particles_to_mesh

  type type_tracer_tracking_model_particles_netcdf
    !< NetCDF IDs for the particle-tracing output file

    character(len=1024) :: filename

    integer :: id_dim_n
    integer :: id_dim_three
    integer :: id_dim_time

    integer :: id_var_time
    integer :: id_var_id
    integer :: id_var_r
    integer :: id_var_t_origin

  end type type_tracer_tracking_model_particles_netcdf

  type type_tracer_tracking_model_particles
    !< The data structure for the particle-based tracer tracking model

    integer                               :: n_max          !< Size of allocated memory
    logical,  dimension(:  ), allocatable :: is_in_use      !< Whether or not memory slot i is in use
    integer,  dimension(:  ), allocatable :: id             !< Unique ID for each particle
    integer                               :: id_max         !< Largest ID used so far (i.e. how many particles have existed in this simulation)
    real(dp), dimension(:,:), allocatable :: r              !< Particle coordinates (x,y,z)
    real(dp), dimension(:  ), allocatable :: zeta           !< Particle zeta coordinate
    integer,  dimension(:  ), allocatable :: vi_in          !< Which mesh vertex' Voronoi cell each particle is located in
    integer,  dimension(:  ), allocatable :: ti_in          !< Which mesh triangle             each particle is located in
    real(dp), dimension(:,:), allocatable :: u              !< Particle velocity (u,v,w)
    real(dp), dimension(:,:), allocatable :: r_origin       !< Coordinates of particle origin (x,y,z)
    real(dp), dimension(:  ), allocatable :: t_origin       !< Time of particle origin (i.e. time of snow deposition)
    real(dp), dimension(:,:), allocatable :: tracers        !< Values of different tracers

    real(dp), dimension(:  ), allocatable :: t0             !< Previous timeframe
    real(dp), dimension(:  ), allocatable :: t1             !< Next timeframe
    real(dp), dimension(:,:), allocatable :: r_t0           !< Particle position at time t0
    real(dp), dimension(:,:), allocatable :: r_t1           !< Particle position at time t1

    type(type_map_particles_to_mesh)      :: map            !< Map from the particles to the mesh

    type(type_tracer_tracking_model_particles_netcdf) :: nc   !< NetCDF IDs for the output file

  end type type_tracer_tracking_model_particles

  type type_tracer_tracking_model
    !< The main tracer tracking model data structure.

    real(dp), dimension(:,:  ), allocatable :: age           !< Age of ice [nV, zeta]
    real(dp), dimension(:,:,:), allocatable :: tracers       !< Age of ice [nV, zeta, i_tracer]

    type(type_tracer_tracking_model_particles) :: particles

  end type type_tracer_tracking_model

contains

end module tracer_tracking_model_types