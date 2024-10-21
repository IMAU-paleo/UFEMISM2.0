module tracer_tracking_model_types

  use precisions, only: dp

  implicit none

  private

  public :: type_tracer_tracking_model, type_tracer_tracking_model_particles

  type type_tracer_tracking_model_particles
    ! The data structure for the particle-based tracer tracking model

    logical,  dimension(:  ), allocatable :: is_in_use     ! Whether or not memory slot i is in use
    real(dp), dimension(:,:), allocatable :: r             ! Particle coordinates (x,y,z)
    integer,  dimension(:  ), allocatable :: vi_in         ! Which mesh vertex' Voronoi cell each particle is located in
    integer,  dimension(:  ), allocatable :: ti_in         ! Which mesh triangle             each particle is located in
    real(dp), dimension(:,:), allocatable :: u             ! Particle velocity    (u,v,w)
    real(dp), dimension(:,:), allocatable :: r_origin      ! Coordinates of particle origin (x,y,z)
    real(dp), dimension(:  ), allocatable :: t_origin      ! Time of particle origin (i.e. time of snow deposition)
    real(dp), dimension(:,:), allocatable :: tracers       ! Values of different tracers

  end type type_tracer_tracking_model_particles

  type type_tracer_tracking_model
    ! The main tracer tracking model data structure.

    real(dp), dimension(:,:  ), allocatable :: age           ! Age of ice [nV, zeta]
    real(dp), dimension(:,:,:), allocatable :: tracers       ! Age of ice [nV, zeta, i_tracer]

    type(type_tracer_tracking_model_particles) :: particles

  end type type_tracer_tracking_model

contains

end module tracer_tracking_model_types