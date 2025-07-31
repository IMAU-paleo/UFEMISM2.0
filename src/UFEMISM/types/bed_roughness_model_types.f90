module bed_roughness_model_types

  use precisions, only: dp

  implicit none

  type type_bed_roughness_nudging_model_H_dHdt_flowline

    ! Nudging masks
    logical,  dimension(:), allocatable :: mask_calc_dCdt_from_nudging
    logical,  dimension(:), allocatable :: mask_calc_dCdt_from_extrapolation
    integer,  dimension(:), allocatable :: mask_extrapolation

    ! Half-flowline-averaged deltaHs and dHs/dt
    real(dp), dimension(:), allocatable :: deltaHs_av_up
    real(dp), dimension(:), allocatable :: deltaHs_av_down
    real(dp), dimension(:), allocatable :: dHs_dt_av_up
    real(dp), dimension(:), allocatable :: dHs_dt_av_down

    ! Intermediate terms
    real(dp), dimension(:), allocatable :: R
    real(dp), dimension(:), allocatable :: I_tot
    real(dp), dimension(:), allocatable :: dC_dt

  end type type_bed_roughness_nudging_model_H_dHdt_flowline

  type type_bed_roughness_nudging_model_H_dHdt_local

    ! Nudging masks
    logical,  dimension(:), allocatable :: mask_calc_dCdt_from_nudging
    logical,  dimension(:), allocatable :: mask_calc_dCdt_from_extrapolation
    integer,  dimension(:), allocatable :: mask_extrapolation

    ! Intermediate terms
    real(dp), dimension(:), allocatable :: C
    real(dp), dimension(:), allocatable :: Laplac_C
    real(dp), dimension(:), allocatable :: dC_dt

  end type type_bed_roughness_nudging_model_H_dHdt_local

  type type_bed_roughness_nudging_model_H_u_flowline

    ! Nudging masks
    logical,  dimension(:), allocatable :: mask_calc_dCdt_from_nudging
    logical,  dimension(:), allocatable :: mask_calc_dCdt_from_extrapolation
    integer,  dimension(:), allocatable :: mask_extrapolation

    ! Target velocity field
    real(dp), dimension(:), allocatable :: uabs_surf_target_b

    ! Half-flowline-averaged deltaHs and deltau
    real(dp), dimension(:), allocatable :: deltaHs_av_up
    real(dp), dimension(:), allocatable :: deltaHs_av_down
    real(dp), dimension(:), allocatable :: deltau_av_up
    real(dp), dimension(:), allocatable :: deltau_av_down

    ! Intermediate terms
    real(dp), dimension(:), allocatable :: R
    real(dp), dimension(:), allocatable :: I_tot
    real(dp), dimension(:), allocatable :: Laplac_C
    real(dp), dimension(:), allocatable :: dC_dt

  end type type_bed_roughness_nudging_model_H_u_flowline

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

    ! Different nudging models
    type(type_bed_roughness_nudging_model_H_dHdt_flowline) :: nudging_H_dHdt_flowline
    type(type_bed_roughness_nudging_model_H_dHdt_local)    :: nudging_H_dHdt_local
    type(type_bed_roughness_nudging_model_H_u_flowline   ) :: nudging_H_u_flowline

  end type type_bed_roughness_model

contains

end module bed_roughness_model_types