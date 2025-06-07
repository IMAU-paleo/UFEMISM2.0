module bed_roughness_model_types

  use precisions, only: dp

  implicit none

  type type_bed_roughness_model

    ! Bed roughness as described in different sliding laws
    real(dp), dimension(:), allocatable :: till_friction_angle         ! [degrees]          Till friction angle
    real(dp), dimension(:), allocatable :: alpha_sq                    ! [-]                Coulomb-law friction coefficient (used when choice_sliding_law = "Tsai2015", or "Schoof2005")
    real(dp), dimension(:), allocatable :: beta_sq                     ! [Pa m^−1/m yr^1/m] Power-law friction coefficient (used when choice_sliding_law = "Weertman", "Tsai2015", or "Schoof2005")

    ! Main data fields
    real(dp), dimension(:), allocatable :: generic_bed_roughness

    ! Timestepping
    real(dp), dimension(:), allocatable :: generic_bed_roughness_prev
    real(dp), dimension(:), allocatable :: generic_bed_roughness_next
    real(dp)                            :: t_prev, t_next

  end type type_bed_roughness_model

contains

end module bed_roughness_model_types