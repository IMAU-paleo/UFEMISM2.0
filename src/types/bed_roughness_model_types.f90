module bed_roughness_model_types

  use precisions, only: dp

  implicit none

  type type_bed_roughness_model

    ! Main data fields
    real(dp), dimension(:), allocatable :: generic_bed_roughness_1
    real(dp), dimension(:), allocatable :: generic_bed_roughness_2

    ! Timestepping
    real(dp), dimension(:), allocatable :: generic_bed_roughness_1_prev
    real(dp), dimension(:), allocatable :: generic_bed_roughness_2_prev
    real(dp), dimension(:), allocatable :: generic_bed_roughness_1_next
    real(dp), dimension(:), allocatable :: generic_bed_roughness_2_next
    real(dp)                            :: t_prev, t_next

  end type type_bed_roughness_model

contains

end module bed_roughness_model_types