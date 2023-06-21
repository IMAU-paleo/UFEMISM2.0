MODULE bed_roughness

  ! Contains all the routines for calculating the bed roughness.

! ===== Preamble =====
! ====================

#include <petsc/finclude/petscksp.h>
  USE petscksp
  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  USE reallocate_mod                                         , ONLY: reallocate_clean_dp_1D
  USE analytical_solutions                                   , ONLY: Schoof2006_icestream
  USE netcdf_input                                           , ONLY: read_field_from_file_2D

  IMPLICIT NONE

CONTAINS

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

    ! Allocate memory for bed roughness fields
    ALLOCATE( ice%phi_fric( mesh%vi1:mesh%vi2))
    ALLOCATE( ice%tau_c(    mesh%vi1:mesh%vi2))
    ALLOCATE( ice%alpha_sq( mesh%vi1:mesh%vi2))
    ALLOCATE( ice%beta_sq(  mesh%vi1:mesh%vi2))

    ice%phi_fric = 0._dp
    ice%tau_c    = 0._dp
    ice%alpha_sq = 0._dp
    ice%beta_sq  = 0._dp

!    ! Inversion-stuff allocation
!    ! ==========================
!
!    IF (C%do_BIVgeo) THEN
!      IF (C%choice_sliding_law == 'Weertman' .OR. &
!          C%choice_sliding_law == 'Tsai2015' .OR. &
!          C%choice_sliding_law == 'Schoof2005') THEN
!
!        CALL allocate_shared_dp_1D( mesh%nV, ice%beta_sq_inv_a , ice%wbeta_sq_inv_a )
!
!      ELSEIF (C%choice_sliding_law == 'Coulomb' .OR. &
!              C%choice_sliding_law == 'Coulomb_regularised' .OR. &
!              C%choice_sliding_law == 'Zoet-Iverson') THEN
!
!        CALL allocate_shared_dp_1D( mesh%nV, ice%phi_fric_inv_a, ice%wphi_fric_inv_a)
!
!      ELSE
!        CALL crash('choice_sliding_law "' // TRIM( C%choice_sliding_law) // '" not compatible with basal sliding inversion!')
!      END IF
!    END IF

    ! Initialisation
    ! ==============

    IF (C%choice_bed_roughness == 'uniform') THEN
      ! Apply a uniform bed roughness

      IF     (C%choice_sliding_law == 'no_sliding') THEN
        ! No need to do anything
      ELSEIF (C%choice_sliding_law == 'idealised') THEN
        ! No need to do anything
      ELSEIF (C%choice_sliding_law == 'Weertman') THEN
        ! Weertman sliding law; bed roughness is described by beta_sq
        ice%beta_sq = C%slid_Weertman_beta_sq_uniform
      ELSEIF (C%choice_sliding_law == 'Coulomb') THEN
        ! Coulomb sliding law; bed roughness is described by phi_fric
        ice%phi_fric = C%slid_Coulomb_phi_fric_uniform
      ELSEIF (C%choice_sliding_law == 'Budd') THEN
        ! Budd-type sliding law; bed roughness is described by phi_fric
        ice%phi_fric = C%slid_Coulomb_phi_fric_uniform
      ELSEIF (C%choice_sliding_law == 'Tsai2015') THEN
        ! Tsai2015 sliding law; bed roughness is described by alpha_sq for the Coulomb part, and beta_sq for the Weertman part
        ice%alpha_sq = C%slid_Tsai2015_alpha_sq_uniform
        ice%beta_sq  = C%slid_Tsai2015_beta_sq_uniform
      ELSEIF (C%choice_sliding_law == 'Schoof2005') THEN
        ! Schoof2005 sliding law; bed roughness is described by alpha_sq for the Coulomb part, and beta_sq for the Weertman part
        ice%alpha_sq = C%slid_Schoof2005_alpha_sq_uniform
        ice%beta_sq  = C%slid_Schoof2005_beta_sq_uniform
      ELSEIF (C%choice_sliding_law == 'Zoet-Iverson') THEN
        ! Zoet-Iverson sliding law; bed roughness is described by phi_fric
        ice%phi_fric = C%slid_ZI_phi_fric_uniform
      ELSE
        CALL crash('unknown choice_sliding_law "' // TRIM( C%choice_sliding_law) // '"!')
      END IF

    ELSEIF (C%choice_bed_roughness == 'parameterised') THEN
      ! Apply the chosen parameterisation of bed roughness

      CALL initialise_bed_roughness_parameterised( mesh, ice)

    ELSEIF (C%choice_bed_roughness == 'read_from_file') THEN
      ! Initialise bed roughness from a NetCDF file

      CALL initialise_bed_roughness_from_file( mesh, ice, region_name)

    ELSE
      CALL crash('unknown choice_bed_roughness "' // TRIM( C%choice_bed_roughness) // '"!')
    END IF

!    ! Inversion initialisation
!    ! ========================
!
!    IF (C%do_BIVgeo) THEN
!      IF (C%choice_sliding_law == 'Weertman' .OR. &
!          C%choice_sliding_law == 'Tsai2015' .OR. &
!          C%choice_sliding_law == 'Schoof2005') THEN
!
!        ice%beta_sq_inv_a(  mesh%vi1:mesh%vi2) = ice%beta_sq_a(  mesh%vi1:mesh%vi2)
!
!      ELSEIF (C%choice_sliding_law == 'Coulomb' .OR. &
!              C%choice_sliding_law == 'Coulomb_regularised' .OR. &
!              C%choice_sliding_law == 'Zoet-Iverson') THEN
!
!        ice%phi_fric_inv_a( mesh%vi1:mesh%vi2) = ice%phi_fric_a( mesh%vi1:mesh%vi2)
!
!      ELSE
!        CALL crash('choice_sliding_law "' // TRIM( C%choice_sliding_law) // '" not compatible with basal sliding inversion!')
!      END IF
!    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_bed_roughness

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

    IF     (C%choice_bed_roughness_parameterised == 'Martin2011') THEN
      CALL calc_bed_roughness_Martin2011( mesh, ice)
    ELSEIF (C%choice_bed_roughness_parameterised == 'SSA_icestream') THEN
      CALL calc_bed_roughness_SSA_icestream( mesh, ice)
    ELSEIF (C%choice_bed_roughness_parameterised == 'MISMIPplus') THEN
      CALL calc_bed_roughness_MISMIPplus( mesh, ice)
    ELSE
      CALL crash('unknown choice_bed_roughness_parameterised "' // TRIM( C%choice_bed_roughness_parameterised) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_bed_roughness_parameterised

  ! The Martin et al. (2011) till parameterisation
  SUBROUTINE calc_bed_roughness_Martin2011( mesh, ice)
    ! Calculate the till friction angle phi_fric and till yield stress tauc,
    ! using the till model by Martin et al. (2011).
    !
    ! Only applicable when choice_sliding_law = "Coulomb", "Budd", or "Zoet-Iverson"

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_bed_roughness_Martin2011'
    INTEGER                                            :: vi
    REAL(dp)                                           :: w_Hb

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (.NOT. (C%choice_sliding_law == 'Coulomb' .OR. C%choice_sliding_law == 'Budd' .OR. C%choice_sliding_law == 'Zoet-Iverson')) THEN
      CALL crash('only applicable when choice_sliding_law = "Coulomb", "Budd", or "Zoet-Iverson"!')
    END IF

    DO vi = mesh%vi1, mesh%vi2

      ! Martin et al. (2011) Eq. 10
      w_Hb = MIN( 1._dp, MAX( 0._dp, (ice%Hb( vi) - C%Martin2011till_phi_Hb_min) / (C%Martin2011till_phi_Hb_max - C%Martin2011till_phi_Hb_min) ))
      ice%phi_fric( vi) = (1._dp - w_Hb) * C%Martin2011till_phi_min + w_Hb * C%Martin2011till_phi_max

    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_bed_roughness_Martin2011

  ! Idealised cases
  SUBROUTINE calc_bed_roughness_SSA_icestream( mesh, ice)
    ! Determine the basal conditions underneath the ice
    !
    ! Idealised case: SSA_icestream (i.e. the Schoof 2006 analytical solution)

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_bed_roughness_SSA_icestream'
    INTEGER                                            :: vii, vi
    REAL(dp)                                           :: y, u

    ! Add routine to path
    CALL init_routine( routine_name)

    DO vii = 1, mesh%nV_loc
      vi = vii + mesh%vi1 - 1
      y = mesh%V( vi,2)
      CALL Schoof2006_icestream( C%uniform_flow_factor, C%Glens_flow_law_exponent, C%refgeo_idealised_SSA_icestream_Hi, &
        C%refgeo_idealised_SSA_icestream_dhdx, C%refgeo_idealised_SSA_icestream_L, C%refgeo_idealised_SSA_icestream_m, &
        y, u, ice%tau_c( vii))
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_bed_roughness_SSA_icestream

  SUBROUTINE calc_bed_roughness_MISMIPplus( mesh, ice)
    ! Determine the basal conditions underneath the ice
    !
    ! Idealised case: MISMIP+ (see Asay-Davis et al., 2016)

    IMPLICIT NONE

    ! Local variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_bed_roughness_MISMIPplus'
    REAL(dp), PARAMETER                                :: MISMIPplus_alpha_sq = 0.5_dp   ! Coulomb-law friction coefficient [unitless];         see Asay-Davis et al., 2016
    REAL(dp), PARAMETER                                :: MISMIPplus_beta_sq  = 1.0E4_dp ! Power-law friction coefficient   [Pa m^−1/3 yr^1/3]; idem dito
    INTEGER                                            :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    IF     (C%choice_sliding_law == 'Weertman') THEN
      ! Uniform sliding factor for the MISMIP+ configuration, using the first (Weertman) sliding law option

      ice%beta_sq = MISMIPplus_beta_sq

    ELSEIF (C%choice_sliding_law == 'Tsai2015') THEN
      ! Uniform sliding factor for the MISMIP+ configuration, using the second (Tsai et al., 2015) sliding law option

      ice%alpha_sq = MISMIPplus_alpha_sq
      ice%beta_sq  = MISMIPplus_beta_sq

    ELSEIF (C%choice_sliding_law == 'Schoof2005') THEN
      ! Uniform sliding factor for the MISMIP+ configuration, using the third (Schoof, 2005) sliding law option

      ice%alpha_sq = MISMIPplus_alpha_sq
      ice%beta_sq  = MISMIPplus_beta_sq

    ELSE
      CALL crash('only defined when choice_sliding_law = "Weertman", "Tsai2015", or "Schoof2005"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_bed_roughness_MISMIPplus

  ! Initialise bed roughness from a file
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
      ! Weertman sliding law; bed roughness is described by beta_sq

      IF (timeframe_bed_roughness == 1E9_dp) THEN
        CALL read_field_from_file_2D( filename_bed_roughness, 'beta_sq', mesh, ice%beta_sq)
      ELSE
        CALL read_field_from_file_2D( filename_bed_roughness, 'beta_sq', mesh, ice%beta_sq, timeframe_bed_roughness)
      END IF

    ELSEIF (C%choice_sliding_law == 'Coulomb') THEN
      ! Coulomb sliding law; bed roughness is described by phi_fric

      IF (timeframe_bed_roughness == 1E9_dp) THEN
        CALL read_field_from_file_2D( filename_bed_roughness, 'phi_fric', mesh, ice%phi_fric)
      ELSE
        CALL read_field_from_file_2D( filename_bed_roughness, 'phi_fric', mesh, ice%phi_fric, timeframe_bed_roughness)
      END IF

    ELSEIF (C%choice_sliding_law == 'Budd') THEN
      ! Budd-type sliding law; bed roughness is described by phi_fric

      IF (timeframe_bed_roughness == 1E9_dp) THEN
        CALL read_field_from_file_2D( filename_bed_roughness, 'phi_fric', mesh, ice%phi_fric)
      ELSE
        CALL read_field_from_file_2D( filename_bed_roughness, 'phi_fric', mesh, ice%phi_fric, timeframe_bed_roughness)
      END IF

    ELSEIF (C%choice_sliding_law == 'Tsai2015') THEN
      ! Tsai2015 sliding law; bed roughness is described by alpha_sq for the Coulomb part, and beta_sq for the Weertman part

      IF (timeframe_bed_roughness == 1E9_dp) THEN
        CALL read_field_from_file_2D( filename_bed_roughness, 'alpha_sq', mesh, ice%alpha_sq)
        CALL read_field_from_file_2D( filename_bed_roughness, 'beta_sq' , mesh, ice%beta_sq)
      ELSE
        CALL read_field_from_file_2D( filename_bed_roughness, 'alpha_sq', mesh, ice%alpha_sq, timeframe_bed_roughness)
        CALL read_field_from_file_2D( filename_bed_roughness, 'beta_sq' , mesh, ice%beta_sq, timeframe_bed_roughness)
      END IF

    ELSEIF (C%choice_sliding_law == 'Schoof2005') THEN
      ! Schoof2005 sliding law; bed roughness is described by alpha_sq for the Coulomb part, and beta_sq for the Weertman part

      IF (timeframe_bed_roughness == 1E9_dp) THEN
        CALL read_field_from_file_2D( filename_bed_roughness, 'alpha_sq', mesh, ice%alpha_sq)
        CALL read_field_from_file_2D( filename_bed_roughness, 'beta_sq' , mesh, ice%beta_sq)
      ELSE
        CALL read_field_from_file_2D( filename_bed_roughness, 'alpha_sq', mesh, ice%alpha_sq, timeframe_bed_roughness)
        CALL read_field_from_file_2D( filename_bed_roughness, 'beta_sq' , mesh, ice%beta_sq, timeframe_bed_roughness)
      END IF

    ELSEIF (C%choice_sliding_law == 'Zoet-Iverson') THEN
      ! Zoet-Iverson sliding law; bed roughness is described by phi_fric

      IF (timeframe_bed_roughness == 1E9_dp) THEN
        CALL read_field_from_file_2D( filename_bed_roughness, 'phi_fric', mesh, ice%phi_fric)
      ELSE
        CALL read_field_from_file_2D( filename_bed_roughness, 'phi_fric', mesh, ice%phi_fric, timeframe_bed_roughness)
      END IF

    ELSE
      CALL crash('unknown choice_sliding_law "' // TRIM( C%choice_sliding_law) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_bed_roughness_from_file

  SUBROUTINE remap_bed_roughness( mesh_old, mesh_new, ice)
    ! Remap or reallocate all the data fields

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_new
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'remap_bed_roughness'
    INTEGER                                            :: int_dummy

    ! Add routine to path
    CALL init_routine( routine_name)

    ! DENK DROM
    CALL crash('fixme!')

!    ! To prevent compiler warnings for unused variables
!    int_dummy = mesh_old%nV
!    int_dummy = mesh_new%nV
!
!    ! == Reallocate/Remap shared memory
!    ! =================================
!
!    IF     (C%choice_sliding_law == 'no_sliding') THEN
!      ! No sliding allowed
!
!    ELSEIF (C%choice_sliding_law == 'idealised') THEN
!      ! Sliding laws for some idealised experiments
!
!    ELSEIF (C%choice_sliding_law == 'Weertman') THEN
!      ! Power-law sliding law
!      IF (C%do_BIVgeo .OR. C%choice_basal_roughness == 'restart') THEN
!        CALL remap_field_dp_2D( mesh_old, mesh_new, ice%beta_sq_a, ice%wbeta_sq_a, 'cons_2nd_order')
!        IF (C%do_BIVgeo) THEN
!          CALL remap_field_dp_2D( mesh_old, mesh_new, ice%beta_sq_inv_a, ice%wbeta_sq_inv_a, 'cons_2nd_order')
!        END IF
!      ELSE
!        CALL reallocate_shared_dp_1D( mesh_new%nV, ice%beta_sq_a , ice%wbeta_sq_a )
!      END IF
!
!    ELSEIF (C%choice_sliding_law == 'Tsai2015' .OR. &
!            C%choice_sliding_law == 'Schoof2005') THEN
!      ! Modified power-law relation
!      IF (C%do_BIVgeo .OR. C%choice_basal_roughness == 'restart') THEN
!        CALL remap_field_dp_2D( mesh_old, mesh_new, ice%alpha_sq_a, ice%walpha_sq_a, 'cons_2nd_order')
!        CALL remap_field_dp_2D( mesh_old, mesh_new, ice%beta_sq_a,  ice%wbeta_sq_a,  'cons_2nd_order')
!        IF (C%do_BIVgeo) THEN
!          CALL remap_field_dp_2D( mesh_old, mesh_new, ice%beta_sq_inv_a, ice%wbeta_sq_inv_a, 'cons_2nd_order')
!        END IF
!      ELSE
!        CALL reallocate_shared_dp_1D( mesh_new%nV, ice%alpha_sq_a, ice%walpha_sq_a)
!        CALL reallocate_shared_dp_1D( mesh_new%nV, ice%beta_sq_a , ice%wbeta_sq_a )
!      END IF
!
!    ELSEIF (C%choice_sliding_law == 'Coulomb' .OR. &
!            C%choice_sliding_law == 'Coulomb_regularised' .OR. &
!            C%choice_sliding_law == 'Zoet-Iverson') THEN
!      ! Yield-stress sliding law
!      IF (C%do_BIVgeo .OR. C%choice_basal_roughness == 'restart') THEN
!        CALL remap_field_dp_2D( mesh_old, mesh_new, ice%phi_fric_a, ice%wphi_fric_a, 'cons_2nd_order')
!        CALL remap_field_dp_2D( mesh_old, mesh_new, ice%tauc_a,     ice%wtauc_a,     'cons_2nd_order')
!        IF (C%do_BIVgeo) THEN
!          CALL remap_field_dp_2D( mesh_old, mesh_new, ice%phi_fric_inv_a, ice%wphi_fric_inv_a, 'cons_2nd_order')
!        END IF
!      ELSE
!        CALL reallocate_shared_dp_1D( mesh_new%nV, ice%phi_fric_a, ice%wphi_fric_a)
!        CALL reallocate_shared_dp_1D( mesh_new%nV, ice%tauc_a    , ice%wtauc_a    )
!      END IF
!
!    ELSE
!      CALL crash('unknown choice_sliding_law "' // TRIM( C%choice_sliding_law) // '"!')
!    END IF
!
!    ! == Reinitialise values
!    ! ======================
!
!    IF (.NOT. C%do_BIVgeo) THEN
!      ! Do not reset the values if we are doing an inversion (values from previous
!      ! time-step and mesh remapped in the previous step above)
!
!      IF (C%choice_basal_roughness == 'uniform') THEN
!        ! Apply a uniform bed roughness
!
!        IF     (C%choice_sliding_law == 'Weertman') THEN
!          ice%beta_sq_a( mesh_new%vi1:mesh_new%vi2) = C%slid_Weertman_beta_sq_uniform
!        ELSEIF (C%choice_sliding_law == 'Coulomb' .OR. &
!                C%choice_sliding_law == 'Coulomb_regularised') THEN
!          ice%phi_fric_a( mesh_new%vi1:mesh_new%vi2) = C%slid_Coulomb_phi_fric_uniform
!        ELSEIF (C%choice_sliding_law == 'Tsai2015') THEN
!          ice%alpha_sq_a( mesh_new%vi1:mesh_new%vi2) = C%slid_Tsai2015_alpha_sq_uniform
!          ice%beta_sq_a(  mesh_new%vi1:mesh_new%vi2) = C%slid_Tsai2015_beta_sq_uniform
!        ELSEIF (C%choice_sliding_law == 'Schoof2005') THEN
!          ice%alpha_sq_a( mesh_new%vi1:mesh_new%vi2) = C%slid_Schoof2005_alpha_sq_uniform
!          ice%beta_sq_a(  mesh_new%vi1:mesh_new%vi2) = C%slid_Schoof2005_beta_sq_uniform
!        ELSEIF (C%choice_sliding_law == 'Zoet-Iverson') THEN
!          ice%phi_fric_a( mesh_new%vi1:mesh_new%vi2) = C%slid_ZI_phi_fric_uniform
!        ELSE
!          CALL crash('unknown choice_sliding_law "' // TRIM( C%choice_sliding_law) // '"!')
!        END IF ! (C%choice_sliding_law)
!
!      ELSEIF (C%choice_basal_roughness == 'parameterised') THEN
!        ! Apply the chosen parameterisation of bed roughness
!        CALL calc_bed_roughness( mesh_new, ice)
!
!      ELSEIF (C%choice_basal_roughness == 'prescribed') THEN
!        ! If bed roughness is prescribed, read it from the provided NetCDF file
!        CALL initialise_bed_roughness_from_file( mesh_new, ice)
!
!      ELSEIF (C%choice_basal_roughness == 'restart') THEN
!        ! Do nothing, as these values were already remapped in the previous step above
!
!      ELSE
!        CALL crash('unknown choice_basal_roughness "' // TRIM( C%choice_basal_roughness) // '"!')
!
!      END IF ! (C%choice_basal_roughness)
!      CALL sync
!
!    END IF ! (.NOT. C%do_BIVgeo)
!
!    ! Basal inversion target velocity
!    ! ===============================
!
!    IF (C%do_BIVgeo) THEN
!      IF (C%choice_BIVgeo_method == 'Berends2022') THEN
!        CALL deallocate_shared( ice%wBIV_uabs_surf_target)
!        CALL initialise_basal_inversion_target_velocity( mesh_new, ice)
!      END IF
!    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_bed_roughness

END MODULE bed_roughness
