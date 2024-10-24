module tracer_tracking_model_types

  use precisions, only: dp, int8

  implicit none

  private

  public :: type_tracer_tracking_model, type_tracer_tracking_model_particles, &
    type_map_particles_to_mesh

  type type_map_particles_to_mesh
    integer                                 :: n    !
    integer,  dimension(:,:,:), allocatable :: ip   ! Indices of   the n nearest particles
    real(dp), dimension(:,:,:), allocatable :: d    ! Distances to the n nearest particles
  end type type_map_particles_to_mesh

  type type_tracer_tracking_model_particles
    ! The data structure for the particle-based tracer tracking model

    integer                                    :: n           ! Size of allocated memory
    logical,       dimension(:  ), allocatable :: is_in_use   ! Whether or not memory slot i is in use
    integer(int8), dimension(:  ), allocatable :: id          ! Unique ID for each particle
    integer(int8)                              :: id_max      ! Largest ID used so far (i.e. how many particles have existed in this simulation)
    real(dp),      dimension(:,:), allocatable :: r           ! Particle coordinates (x,y,z)
    real(dp),      dimension(:  ), allocatable :: zeta        ! Particle zeta coordinate
    integer,       dimension(:  ), allocatable :: vi_in       ! Which mesh vertex' Voronoi cell each particle is located in
    integer,       dimension(:  ), allocatable :: ti_in       ! Which mesh triangle             each particle is located in
    real(dp),      dimension(:,:), allocatable :: u           ! Particle velocity (u,v,w)
    real(dp),      dimension(:,:), allocatable :: r_origin    ! Coordinates of particle origin (x,y,z)
    real(dp),      dimension(:  ), allocatable :: t_origin    ! Time of particle origin (i.e. time of snow deposition)
    real(dp),      dimension(:,:), allocatable :: tracers     ! Values of different tracers
    type(type_map_particles_to_mesh)           :: map         ! Map from the particles to the mesh

  end type type_tracer_tracking_model_particles

  type type_tracer_tracking_model
    ! The main tracer tracking model data structure.

    real(dp), dimension(:,:  ), allocatable :: age           ! Age of ice [nV, zeta]
    real(dp), dimension(:,:,:), allocatable :: tracers       ! Age of ice [nV, zeta, i_tracer]

    type(type_tracer_tracking_model_particles) :: particles

  end type type_tracer_tracking_model

contains

end module tracer_tracking_model_types