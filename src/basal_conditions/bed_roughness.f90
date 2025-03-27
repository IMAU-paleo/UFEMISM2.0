MODULE bed_roughness

  ! Contains all the routines for calculating the parameters
  ! that determine bed roughness.
  ! NOTE: the variable ice%bed_roughness itself is computed
  ! from these parameters within each sliding law in their
  ! respective module

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE grid_basic                                             , ONLY: type_grid
  USE ice_model_types                                        , ONLY: type_ice_model
  USE reference_geometry_types                               , ONLY: type_reference_geometry
  USE basal_inversion_types                                  , ONLY: type_basal_inversion
  USE analytical_solutions                                   , ONLY: Schoof2006_icestream
  use netcdf_io_main

  IMPLICIT NONE

CONTAINS

! ===== Main routines =====
! =========================

  SUBROUTINE run_bed_roughness_model( mesh, grid_smooth, ice, refgeo, BIV, time)
    ! Run the chosen bed roughness model

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_grid),                     INTENT(IN)    :: grid_smooth
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo
    TYPE(type_basal_inversion),          INTENT(IN)    :: BIV
    REAL(dp),                            INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_bed_roughness_model'
    INTEGER                                            :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate bed roughness using the chosen model
    ! ==============================================

    IF (C%choice_bed_roughness == 'uniform') THEN
      ! Apply a uniform bed roughness

      ! No need to do anything

    ELSEIF (C%choice_bed_roughness == 'parameterised') THEN
      ! Apply the chosen parameterisation of bed roughness

      CALL calc_bed_roughness_parameterised( mesh, ice)

    ELSEIF (C%choice_bed_roughness == 'read_from_file') THEN
      ! Initialise bed roughness from a NetCDF file

      ! No need to do anything

    ELSE
      CALL crash('unknown choice_bed_roughness "' // TRIM( C%choice_bed_roughness) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_bed_roughness_model

  SUBROUTINE initialise_bed_roughness( mesh, ice, region_name)
    ! Initialise the bed roughness

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_bed_roughness'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Initialisation
    ! ==============

    IF (C%choice_bed_roughness == 'uniform') THEN
      ! Apply a uniform bed roughness

      CALL initialise_bed_roughness_uniform( mesh, ice)

    ELSEIF (C%choice_bed_roughness == 'parameterised') THEN
      ! Apply the chosen parameterisation of bed roughness

      CALL initialise_bed_roughness_parameterised( mesh, ice)

    ELSEIF (C%choice_bed_roughness == 'read_from_file') THEN
      ! Initialise bed roughness from a NetCDF file

      CALL initialise_bed_roughness_from_file( mesh, ice, region_name)

    ELSE
      CALL crash('unknown choice_bed_roughness "' // TRIM( C%choice_bed_roughness) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_bed_roughness

  SUBROUTINE remap_bed_roughness( mesh_old, mesh_new, ice)
    ! Remap or reallocate all the data fields

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_new
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'remap_bed_roughness'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! DENK DROM
    CALL crash('fixme!')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_bed_roughness

! ===== Different bed roughness models ====
! =========================================

  ! == Uniform bed roughness
  SUBROUTINE initialise_bed_roughness_uniform( mesh, ice)
    ! Initialise bed roughness
    !
    ! Use a uniform value over the whole domain

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_bed_roughness_uniform'
    INTEGER                                            :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Initialise bed roughness
    ! ========================

    ! Initialise field based on chosen sliding law
    SELECT CASE (C%choice_sliding_law)

      CASE ('no_sliding')
        ! No need to do anything

      CASE ('idealised')
        ! No need to do anything

      CASE ('Weertman')
        ! Weertman sliding law; bed roughness is described by slid_beta_sq
        ice%slid_beta_sq = C%slid_Weertman_beta_sq_uniform

      CASE ('Coulomb')
        ! Coulomb sliding law; bed roughness is described by till_friction_angle
        ice%till_friction_angle = C%slid_Coulomb_phi_fric_uniform

      CASE ('Budd')
        ! Budd-type sliding law; bed roughness is described by till_friction_angle
        ice%till_friction_angle = C%slid_Coulomb_phi_fric_uniform

      CASE ('Tsai2015')
        ! Tsai2015 sliding law; bed roughness is described by slid_alpha_sq for the Coulomb part, and slid_beta_sq for the Weertman part
        ice%slid_alpha_sq = C%slid_Tsai2015_alpha_sq_uniform
        ice%slid_beta_sq  = C%slid_Tsai2015_beta_sq_uniform

      CASE ('Schoof2005')
        ! Schoof2005 sliding law; bed roughness is described by slid_alpha_sq for the Coulomb part, and slid_beta_sq for the Weertman part
        ice%slid_alpha_sq = C%slid_Schoof2005_alpha_sq_uniform
        ice%slid_beta_sq  = C%slid_Schoof2005_beta_sq_uniform

      CASE ('Zoet-Iverson')
        ! Zoet-Iverson sliding law; bed roughness is described by till_friction_angle
        ice%till_friction_angle = C%slid_ZI_phi_fric_uniform

      CASE DEFAULT
        CALL crash('unknown choice_sliding_law "' // TRIM( C%choice_sliding_law) // '"!')

    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_bed_roughness_uniform

  ! == Parameterised bed roughness
  SUBROUTINE calc_bed_roughness_parameterised( mesh, ice)
    ! Compute bed roughness
    ! Use a simple parameterisation to calculate bed roughness

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_bed_roughness_parameterised'

    ! Add routine to path
    CALL init_routine( routine_name)

    SELECT CASE (C%choice_bed_roughness_parameterised)
      CASE ('Martin2011')
        CALL calc_bed_roughness_Martin2011( mesh, ice)
      CASE ('SSA_icestream')
        ! No need to do anything
      CASE ('MISMIPplus')
        ! No need to do anything
      CASE ('MISMIP+')
        ! No need to do anything
      CASE DEFAULT
        CALL crash('unknown choice_bed_roughness_parameterised "' // TRIM( C%choice_bed_roughness_parameterised) // '"!')
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_bed_roughness_parameterised

  SUBROUTINE initialise_bed_roughness_parameterised( mesh, ice)
    ! Initialise the bed roughness
    ! Use a simple parameterisation to calculate bed roughness

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_bed_roughness_parameterised'

    ! Add routine to path
    CALL init_routine( routine_name)

    SELECT CASE (C%choice_bed_roughness_parameterised)
      CASE ('Martin2011')
        CALL initialise_bed_roughness_Martin2011( mesh, ice)
      CASE ('SSA_icestream')
        CALL initialise_bed_roughness_SSA_icestream( mesh, ice)
      CASE ('MISMIPplus')
        CALL initialise_bed_roughness_MISMIPplus( mesh, ice)
      CASE ('MISMIP+')
        CALL initialise_bed_roughness_MISMIPplus( mesh, ice)
      CASE DEFAULT
        CALL crash('unknown choice_bed_roughness_parameterised "' // TRIM( C%choice_bed_roughness_parameterised) // '"!')
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_bed_roughness_parameterised

  ! The Martin et al. (2011) till parameterisation
  SUBROUTINE calc_bed_roughness_Martin2011( mesh, ice)
    ! Calculate the till friction angle using the till model by Martin et al. (2011).
    !
    ! Only applicable when choice_sliding_law = "Coulomb", "Budd", or "Zoet-Iverson"

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_bed_roughness_Martin2011'
    INTEGER                                            :: vi
    REAL(dp)                                           :: weight_Hb

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (.NOT. (C%choice_sliding_law == 'Coulomb' .OR. C%choice_sliding_law == 'Budd' .OR. C%choice_sliding_law == 'Zoet-Iverson')) THEN
      CALL crash('only applicable when choice_sliding_law = "Coulomb", "Budd", or "Zoet-Iverson"!')
    END IF

    DO vi = mesh%vi1, mesh%vi2

      ! Compute till friction angle based on Martin et al. (2011) Eq. 10
      weight_Hb = MIN( 1._dp, MAX( 0._dp, (ice%Hb( vi) - C%Martin2011till_phi_Hb_min) / (C%Martin2011till_phi_Hb_max - C%Martin2011till_phi_Hb_min) ))
      ice%till_friction_angle( vi) = (1._dp - weight_Hb) * C%Martin2011till_phi_min + weight_Hb * C%Martin2011till_phi_max

    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_bed_roughness_Martin2011

  SUBROUTINE initialise_bed_roughness_Martin2011( mesh, ice)
    ! Calculate the till friction angle using the till model by Martin et al. (2011).
    !
    ! Only applicable when choice_sliding_law = "Coulomb", "Budd", or "Zoet-Iverson"

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_bed_roughness_Martin2011'
    INTEGER                                            :: vi
    REAL(dp)                                           :: weight_Hb

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (.NOT. (C%choice_sliding_law == 'Coulomb' .OR. C%choice_sliding_law == 'Budd' .OR. C%choice_sliding_law == 'Zoet-Iverson')) THEN
      CALL crash('only applicable when choice_sliding_law = "Coulomb", "Budd", or "Zoet-Iverson"!')
    END IF

    ! Compute till friction angle based on Martin et al. (2011) Eq. 10
    DO vi = mesh%vi1, mesh%vi2

      weight_Hb = MIN( 1._dp, MAX( 0._dp, (ice%Hb( vi) - C%Martin2011till_phi_Hb_min) / (C%Martin2011till_phi_Hb_max - C%Martin2011till_phi_Hb_min) ))

      ice%till_friction_angle( vi) = (1._dp - weight_Hb) * C%Martin2011till_phi_min + weight_Hb * C%Martin2011till_phi_max

    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_bed_roughness_Martin2011

  ! Idealised cases
  SUBROUTINE initialise_bed_roughness_SSA_icestream( mesh, ice)
    ! Determine the basal conditions underneath the ice
    !
    ! Idealised case: SSA_icestream (i.e. the Schoof 2006 analytical solution)

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_bed_roughness_SSA_icestream'
    INTEGER                                            :: vi
    REAL(dp)                                           :: y, u

    ! Add routine to path
    CALL init_routine( routine_name)

    DO vi = mesh%vi1, mesh%vi2
      y = mesh%V( vi,2)
      CALL Schoof2006_icestream( C%uniform_Glens_flow_factor, C%Glens_flow_law_exponent, C%refgeo_idealised_SSA_icestream_Hi, &
        C%refgeo_idealised_SSA_icestream_dhdx, C%refgeo_idealised_SSA_icestream_L, C%refgeo_idealised_SSA_icestream_m, &
        y, u, ice%till_yield_stress( vi))
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_bed_roughness_SSA_icestream

  SUBROUTINE initialise_bed_roughness_MISMIPplus( mesh, ice)
    ! Determine the basal conditions underneath the ice
    !
    ! Idealised case: MISMIP+ (i.e. just a uniform value)

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_bed_roughness_MISMIPplus'
    REAL(dp), PARAMETER                                :: MISMIPplus_alpha_sq = 0.5_dp   ! Coulomb-law friction coefficient [unitless];         see Asay-Davis et al., 2016
    REAL(dp), PARAMETER                                :: MISMIPplus_beta_sq  = 1.0E4_dp ! Power-law friction coefficient   [Pa m^−1/3 yr^1/3]; idem dito

    ! Add routine to path
    CALL init_routine( routine_name)

    SELECT CASE (C%choice_sliding_law)
      CASE ('Weertman')
        ! Uniform sliding factor for the MISMIP+ configuration, using the first (Weertman) sliding law option

        ice%slid_beta_sq  = MISMIPplus_beta_sq

      CASE ('Tsai2015')
        ! Uniform sliding factor for the MISMIP+ configuration, using the second (Tsai et al., 2015) sliding law option

        ice%slid_alpha_sq = MISMIPplus_alpha_sq
        ice%slid_beta_sq  = MISMIPplus_beta_sq

      CASE ('Schoof2005')
        ! Uniform sliding factor for the MISMIP+ configuration, using the third (Schoof, 2005) sliding law option

        ice%slid_alpha_sq = MISMIPplus_alpha_sq
        ice%slid_beta_sq  = MISMIPplus_beta_sq

      CASE DEFAULT
        CALL crash('only defined when choice_sliding_law = "Weertman", "Tsai2015", or "Schoof2005"!')
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_bed_roughness_MISMIPplus

  ! == Bed roughness from an external file
  SUBROUTINE initialise_bed_roughness_from_file( mesh, ice, region_name)
    ! Initialise bed roughness with data from an external NetCDF file

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_bed_roughness_from_file'
    CHARACTER(LEN=256)                                 :: filename_bed_roughness
    REAL(dp)                                           :: timeframe_bed_roughness

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine filename and timeframe for this model region
    IF     (region_name == 'NAM') THEN
      filename_bed_roughness  = C%filename_bed_roughness_NAM
      timeframe_bed_roughness = C%timeframe_bed_roughness_NAM
    ELSEIF (region_name == 'EAS') THEN
      filename_bed_roughness  = C%filename_bed_roughness_EAS
      timeframe_bed_roughness = C%timeframe_bed_roughness_EAS
    ELSEIF (region_name == 'GRL') THEN
      filename_bed_roughness  = C%filename_bed_roughness_GRL
      timeframe_bed_roughness = C%timeframe_bed_roughness_GRL
    ELSEIF (region_name == 'ANT') THEN
      filename_bed_roughness  = C%filename_bed_roughness_ANT
      timeframe_bed_roughness = C%timeframe_bed_roughness_ANT
    ELSE
      CALL crash('unknown region_name "' // TRIM( region_name) // '"!')
    END IF

    IF     (C%choice_sliding_law == 'no_sliding') THEN
      ! No need to do anything
    ELSEIF (C%choice_sliding_law == 'idealised') THEN
      ! No need to do anything
    ELSEIF (C%choice_sliding_law == 'Weertman') THEN
      ! Weertman sliding law; bed roughness is described by slid_beta_sq

      IF (timeframe_bed_roughness == 1E9_dp) THEN
        CALL read_field_from_file_2D( filename_bed_roughness, 'slid_beta_sq', mesh, ice%slid_beta_sq)
      ELSE
        CALL read_field_from_file_2D( filename_bed_roughness, 'slid_beta_sq', mesh, ice%slid_beta_sq, timeframe_bed_roughness)
      END IF

    ELSEIF (C%choice_sliding_law == 'Coulomb') THEN
      ! Coulomb sliding law; bed roughness is described by till_friction_angle

      IF (timeframe_bed_roughness == 1E9_dp) THEN
        CALL read_field_from_file_2D( filename_bed_roughness, 'till_friction_angle||phi_fric', mesh, ice%till_friction_angle)
      ELSE
        CALL read_field_from_file_2D( filename_bed_roughness, 'till_friction_angle||phi_fric', mesh, ice%till_friction_angle, timeframe_bed_roughness)
      END IF

    ELSEIF (C%choice_sliding_law == 'Budd') THEN
      ! Budd-type sliding law; bed roughness is described by till_friction_angle

      IF (timeframe_bed_roughness == 1E9_dp) THEN
        CALL read_field_from_file_2D( filename_bed_roughness, 'till_friction_angle||phi_fric', mesh, ice%till_friction_angle)
      ELSE
        CALL read_field_from_file_2D( filename_bed_roughness, 'till_friction_angle||phi_fric', mesh, ice%till_friction_angle, timeframe_bed_roughness)
      END IF

    ELSEIF (C%choice_sliding_law == 'Tsai2015') THEN
      ! Tsai2015 sliding law; bed roughness is described by slid_alpha_sq for the Coulomb part, and slid_beta_sq for the Weertman part

      IF (timeframe_bed_roughness == 1E9_dp) THEN
        CALL read_field_from_file_2D( filename_bed_roughness, 'slid_alpha_sq', mesh, ice%slid_alpha_sq)
        CALL read_field_from_file_2D( filename_bed_roughness, 'slid_beta_sq' , mesh, ice%slid_beta_sq)
      ELSE
        CALL read_field_from_file_2D( filename_bed_roughness, 'slid_alpha_sq', mesh, ice%slid_alpha_sq, timeframe_bed_roughness)
        CALL read_field_from_file_2D( filename_bed_roughness, 'slid_beta_sq' , mesh, ice%slid_beta_sq, timeframe_bed_roughness)
      END IF

    ELSEIF (C%choice_sliding_law == 'Schoof2005') THEN
      ! Schoof2005 sliding law; bed roughness is described by slid_alpha_sq for the Coulomb part, and slid_beta_sq for the Weertman part

      IF (timeframe_bed_roughness == 1E9_dp) THEN
        CALL read_field_from_file_2D( filename_bed_roughness, 'slid_alpha_sq', mesh, ice%slid_alpha_sq)
        CALL read_field_from_file_2D( filename_bed_roughness, 'slid_beta_sq' , mesh, ice%slid_beta_sq)
      ELSE
        CALL read_field_from_file_2D( filename_bed_roughness, 'slid_alpha_sq', mesh, ice%slid_alpha_sq, timeframe_bed_roughness)
        CALL read_field_from_file_2D( filename_bed_roughness, 'slid_beta_sq' , mesh, ice%slid_beta_sq, timeframe_bed_roughness)
      END IF

    ELSEIF (C%choice_sliding_law == 'Zoet-Iverson') THEN
      ! Zoet-Iverson sliding law; bed roughness is described by till_friction_angle

      IF (timeframe_bed_roughness == 1E9_dp) THEN
        CALL read_field_from_file_2D( filename_bed_roughness, 'till_friction_angle||phi_fric', mesh, ice%till_friction_angle)
      ELSE
        CALL read_field_from_file_2D( filename_bed_roughness, 'till_friction_angle||phi_fric', mesh, ice%till_friction_angle, timeframe_bed_roughness)
      END IF

    ELSE
      CALL crash('unknown choice_sliding_law "' // TRIM( C%choice_sliding_law) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_bed_roughness_from_file

END MODULE bed_roughness
