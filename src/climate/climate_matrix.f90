module climate_matrix

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE control_resources_and_error_messaging                  , ONLY: crash, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  USE climate_model_types                                    , ONLY: type_climate_model
  USE climate_idealised                                      , ONLY: initialise_climate_model_idealised, run_climate_model_idealised
  USE climate_realistic                                      , ONLY: initialise_climate_model_realistic, run_climate_model_realistic
  USE reallocate_mod                                         , ONLY: reallocate_bounds
  use netcdf_io_main
  use mesh_data_smoothing, only: smooth_Gaussian
  
  ! check the previous calleds
  use forcing_module, only: get_insolation_at_time, update_CO2_at_model_time
! added in climate_realistic the subroutines initialise_global_forcing, get_insolation_at_time, update_insolation_timeframes_from_file
  IMPLICIT NONE

CONTAINS

! == Climate matrix
! ===========================

  ! Climate matrix with warm + cold snapshots, forced with CO2 (from record or from inverse routine) from Berends et al., 2018
  ! Generalised for different timeframes, L.B. Stap (2021)
  SUBROUTINE run_climate_model_matrix( mesh, grid, ice, SMB, climate, region_name, time)
    ! Use CO2 (either prescribed or inversely modelled) to force the 2-snapshot (PI-LGM) climate matrix (Berends et al., 2018)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_SMB_model),                INTENT(IN)    :: SMB
    TYPE(type_climate_model),            INTENT(INOUT) :: climate
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    REAL(dp),                            INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_climate_model_matrix'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Update forcing at model time
    CALL get_insolation_at_time( grid, time, climate%Q_TOA)
    CALL update_CO2_at_model_time( time) ! forcing module == NOT YET IN climate_realistic!

    ! Use the (CO2 + absorbed insolation)-based interpolation scheme for temperature
    CALL run_climate_model_matrix_temperature( mesh, grid, ice, SMB, climate, region_name)

    ! Use the (CO2 + ice-sheet geometry)-based interpolation scheme for precipitation
    CALL run_climate_model_matrix_precipitation( mesh, grid, ice, climate, region_name)

    ! == Safety checks
    ! ================

! copied and pasted from UFE1.x Check the names of the variables!
    DO vi = mesh%vi1, mesh%vi2
    DO m = 1, 12
      IF (climate_matrix%applied%T2m( vi,m) < 150._dp) THEN
        CALL warning('excessively low temperatures (<150K) detected!')
      ELSEIF (climate_matrix%applied%T2m( vi,m) < 0._dp) THEN
        CALL crash('negative temperatures (<0K) detected!')
      ELSEIF (climate_matrix%applied%T2m( vi,m) /= climate_matrix%applied%T2m( vi,m)) THEN
        CALL crash('NaN temperatures  detected!')
      ELSEIF (climate_matrix%applied%Precip( vi,m) <= 0._dp) THEN
        CALL crash('zero/negative precipitation detected!')
      ELSEIF (climate_matrix%applied%Precip( vi,m) /= climate_matrix%applied%Precip( vi,m)) THEN
        CALL crash('NaN precipitation detected!')
      END IF
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_climate_model_matrix
  SUBROUTINE run_climate_model_matrix_temperature( mesh, grid, ice, SMB, climate, region_name)
    ! The (CO2 + absorbed insolation)-based matrix interpolation for temperature, from Berends et al. (2018)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_SMB_model),                INTENT(IN)    :: SMB
    TYPE(type_climate_model),            INTENT(INOUT) :: climate
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_climate_model_matrix_temperature'
    INTEGER                                            :: vi ,m
    REAL(dp)                                           :: CO2, w_CO2
    REAL(dp), DIMENSION(:    ), POINTER                ::  w_ins,  w_ins_smooth,  w_ice,  w_tot
    INTEGER                                            :: ww_ins, ww_ins_smooth, ww_ice, ww_tot
    REAL(dp)                                           :: w_ins_av
    REAL(dp), DIMENSION(:,:  ), POINTER                :: T_ref_GCM
    REAL(dp), DIMENSION(:    ), POINTER                :: Hs_GCM, lambda_GCM
    INTEGER                                            :: wT_ref_GCM, wHs_GCM, wlambda_GCM

    REAL(dp), PARAMETER                                :: w_cutoff = 0.5_dp        ! Crop weights to [-w_cutoff, 1 + w_cutoff]
    REAL(dp), PARAMETER                                :: P_offset = 0.008_dp       ! Normalisation term in precipitation anomaly to avoid divide-by-nearly-zero

    ! Add routine to path
    CALL init_routine( routine_name)
!!!
!!! UFEMISM1.x has the same type of mesh so is better to copy and paste from there, it is straightforward
!!! 
    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nV,     w_ins,        ww_ins        )
    CALL allocate_shared_dp_1D( mesh%nV,     w_ins_smooth, ww_ins_smooth )
    CALL allocate_shared_dp_1D( mesh%nV,     w_ice,        ww_ice        )
    CALL allocate_shared_dp_1D( mesh%nV,     w_tot,        ww_tot        )
    CALL allocate_shared_dp_2D( mesh%nV, 12, T_ref_GCM,    wT_ref_GCM    )
    CALL allocate_shared_dp_1D( mesh%nV,     Hs_GCM,       wHs_GCM       )
    CALL allocate_shared_dp_1D( mesh%nV,     lambda_GCM,   wlambda_GCM   )
    
    ! Find CO2 interpolation weight (use either prescribed or modelled CO2)
    ! =====================================================================

    IF     (C%choice_forcing_method == 'CO2_direct') THEN
      CO2 = forcing%CO2_obs
    ELSEIF (C%choice_forcing_method == 'd18O_inverse_CO2') THEN
      CO2 = forcing%CO2_mod
    ELSEIF (C%choice_forcing_method == 'd18O_inverse_dT_glob') THEN
      CO2 = 0._dp
      CALL crash('must only be called with the correct forcing method, check your code!')
    ELSE
      CO2 = 0._dp
      CALL crash('unknown choice_forcing_method"' // TRIM(C%choice_forcing_method) // '"!')
    END IF

    ! If CO2 ~= warm snap -> weight is 1. If ~= cold snap -> weight is 0.
    ! Otherwise interpolate. Berends et al., 2018 - Eq. 1
    w_CO2 = MAX( -w_cutoff, MIN( 1._dp + w_cutoff, (CO2 - C%matrix_low_CO2_level) / &
                               (C%matrix_high_CO2_level - C%matrix_low_CO2_level) ))

    ! Find the interpolation weights based on absorbed insolation
    ! ===========================================================

    ! Calculate modelled absorbed insolation
    ! CHECK where I_Abs will be stored it will be called as climate%matrix? 
    climate%matrix%I_abs( mesh%vi1:mesh%vi2) = 0._dp
    DO m = 1, 12
    DO vi = mesh%vi1, mesh%vi2
      ! Calculate modelled absorbed insolation. Berends et al., 2018 - Eq. 2
      climate%matrix%I_abs( vi) = climate%matrix%I_abs( vi) + & 
      climate%Q_TOA( vi,m) * (1._dp - SMB%Albedo( vi, m))  
    END DO
    END DO

    ! Calculate "direct" weighting field
    ! Berends et al., 2018 - Eq. 3
    DO vi = mesh%vi1, mesh%vi2
      ! If absorbed insolation ~= warm snap -> weight is 1.
      ! If ~= cold snap -> weight is 0. Otherwise interpolate
      w_ins( vi) = MAX( -w_cutoff, MIN( 1._dp + w_cutoff, &
                      ( climate%matrix%I_abs( vi) - climate%matrix%GCM_cold%I_abs( vi)) / &
                      ( climate%matrix%GCM_warm%I_abs( vi) - climate%matrix%GCM_cold%I_abs( vi)) ))
    END DO
    w_ins_av      = MAX( -w_cutoff, MIN( 1._dp + w_cutoff, (SUM( climate%matrix%I_abs         )      - SUM( climate%matrix%GCM_cold%I_abs)     ) / &
                                                           (SUM( climate%matrix%GCM_warm%I_abs)      - SUM( climate%matrix%GCM_cold%I_abs)     ) ))

    ! Smooth the weighting field
    w_ins_smooth( mesh%vi1:mesh%vi2) = w_ins( mesh%vi1:mesh%vi2)
    CALL smooth_Gaussian_2D( mesh, grid, w_ins_smooth, 200000._dp)

    ! Combine unsmoothed, smoothed, and regional average weighting fields (Berends et al., 2018, Eq. 4)
    IF (region_name == 'NAM' .OR. region_name == 'EAS') THEN
      w_ice( mesh%vi1:mesh%vi2) = (1._dp * w_ins(        mesh%vi1:mesh%vi2) + &
                                   3._dp * w_ins_smooth( mesh%vi1:mesh%vi2) + &
                                   3._dp * w_ins_av) / 7._dp
    ELSEIF (region_name == 'GRL' .OR. region_name == 'ANT') THEN
      ! Use only the regional (averaged) and smoothed weights
      w_ice( mesh%vi1:mesh%vi2) = (1._dp * w_ins_smooth( mesh%vi1:mesh%vi2) + &
                                   6._dp * w_ins_av) / 7._dp
    END IF


    ! Combine interpolation weights from absorbed insolation and CO2 into the final weights fields
    ! Berends et al., 2018 - Eqs. 5, 9 with weights 0.5 for NAM & EAS, and 0.75 for ANT
    ! Generalised: "switch" between matrix method and glacial index method by altering C%climate_matrix_CO2vsice_<region>
    IF         (region_name == 'NAM') THEN
      w_tot( mesh%vi1:mesh%vi2) = (C%climate_matrix_CO2vsice_NAM * w_CO2) + ((1._dp - C%climate_matrix_CO2vsice_NAM) * w_ice( mesh%vi1:mesh%vi2))
    ELSEIF     (region_name == 'EAS') THEN
      w_tot( mesh%vi1:mesh%vi2) = (C%climate_matrix_CO2vsice_EAS * w_CO2) + ((1._dp - C%climate_matrix_CO2vsice_EAS) * w_ice( mesh%vi1:mesh%vi2))
    ELSEIF     (region_name == 'GRL') THEN
      w_tot( mesh%vi1:mesh%vi2) = (C%climate_matrix_CO2vsice_GRL * w_CO2) + ((1._dp - C%climate_matrix_CO2vsice_GRL) * w_ice( mesh%vi1:mesh%vi2))
    ELSEIF     (region_name == 'ANT') THEN
      w_tot( mesh%vi1:mesh%vi2) = (C%climate_matrix_CO2vsice_ANT * w_CO2) + ((1._dp - C%climate_matrix_CO2vsice_ANT) * w_ice( mesh%vi1:mesh%vi2))
    END IF

!! ==============================================================================================================
!! In UFE1.x here are two more options, already asked Jorge about it, glacial matrix and glacial index
!! lines 1050 - 1080 in climate_module.f90
!! ==============================================================================================================

    ! Interpolate between the GCM snapshots
    ! =====================================

    DO vi = mesh%vi1, mesh%vi2

      ! Find matrix-interpolated orography, lapse rate, and temperature
      Hs_GCM( vi     ) = (w_tot( vi) * climate%matrix%GCM_warm%Hs( vi    )) + &
                         ((1._dp - w_tot( vi)) * climate%matrix%GCM_cold%Hs( vi    ))  ! Berends et al., 2018 - Eq. 8
      lambda_GCM( vi ) = (w_tot( vi) * climate%matrix%GCM_warm%lambda( vi)) + &
                         ((1._dp - w_tot( vi)) * climate%matrix%GCM_cold%lambda( vi))  ! Not listed in the article, shame on me!
      T_ref_GCM( vi,:) = (w_tot( vi) * climate%matrix%GCM_warm%T2m( vi,: )) + &
                         ((1._dp - w_tot( vi)) * climate%matrix%GCM_cold%T2m( vi,: ))  ! Berends et al., 2018 - Eq. 6

      ! Adapt temperature to model orography using matrix-derived lapse-rate
      DO m = 1, 12
      !! check name variable Hs_a in UFE2 should be just Hs !
        climate%T2m( vi,:) = T_ref_GCM( vi, m) - lambda_GCM( vi) * (ice%Hs_a( vi) - Hs_GCM( vi))  ! Berends et al., 2018 - Eq. 11
      END DO

    END DO

    ! Clean up after yourself
    CALL deallocate_shared( ww_ins)
    CALL deallocate_shared( ww_ins_smooth)
    CALL deallocate_shared( ww_ice)
    CALL deallocate_shared( ww_tot)
    CALL deallocate_shared( wT_ref_GCM)
    CALL deallocate_shared( wHs_GCM)
    CALL deallocate_shared( wlambda_GCM)

    ! Safety checks
    CALL check_safety_temperature( climate%T2m)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_climate_model_matrix_temperature
  SUBROUTINE run_climate_model_matrix_precipitation( mesh, grid, ice, climate, region_name)
    ! The (CO2 + ice geometry)-based matrix interpolation for precipitation, from Berends et al. (2018)
    ! For NAM and EAS, this is based on local ice geometry and uses the Roe&Lindzen precipitation model for downscaling.
    ! For GRL and ANT, this is based on total ice volume,  and uses the simple CC   precipitation model for downscaling.
    ! The rationale for this difference is that glacial-interglacial differences in ice geometry are much more
    ! dramatic in NAM and EAS than they are in GRL and ANT.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_climate_model),            INTENT(INOUT) :: climate
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_climate_model_matrix_precipitation'
    INTEGER                                            :: vi, m
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  w_warm,  w_cold
    INTEGER                                            :: ww_warm, ww_cold
    REAL(dp)                                           :: w_tot
    REAL(dp), DIMENSION(:,:,:), POINTER                :: T_ref_GCM, P_ref_GCM
    REAL(dp), DIMENSION(:,:  ), POINTER                :: Hs_GCM
    INTEGER                                            :: wT_ref_GCM, wP_ref_GCM, wHs_GCM

    REAL(dp), PARAMETER                                :: w_cutoff = 0.25_dp        ! Crop weights to [-w_cutoff, 1 + w_cutoff]

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_dp_1D(     mesh%nV, w_warm,         ww_warm        )
    CALL allocate_shared_dp_1D(     mesh%nV, w_cold,         ww_cold        )
    CALL allocate_shared_dp_2D( 12, mesh%nV, T_ref_GCM,      wT_ref_GCM     )
    CALL allocate_shared_dp_2D( 12, mesh%nV, P_ref_GCM,      wP_ref_GCM     )
    CALL allocate_shared_dp_1D(     mesh%nV, Hs_GCM,         wHs_GCM        )

    ! Calculate interpolation weights based on ice geometry
    ! =====================================================
!! check ice%Hs_a should be ice%Hs now
    ! First calculate the total ice volume term (second term in the equation)
    w_tot = MAX(-w_cutoff, MIN(1._dp + w_cutoff, &
      (SUM( ice%Hs_a) - SUM( climate%matrix%GCM_warm%Hs)) / (SUM( climate%matrix%GCM_cold%Hs) - SUM( climate%matrix%GCM_warm%Hs)) ))

    IF (region_name == 'NAM' .OR. region_name == 'EAS') THEN
      ! Combine total + local ice thicness; Berends et al., 2018, Eq. 12

      ! Then the local ice thickness term
      DO vi = mesh%vi1, mesh%vi2

        IF (climate%matrix%GCM_warm%Hs( vi) < climate%matrix%GCM_PI%Hs( vi) + 50._dp) THEN
          IF (climate%matrix%GCM_cold%Hs( vi) < climate%matrix%GCM_PI%Hs( vi) + 50._dp) THEN
            ! No ice in any GCM state. Use only total ice volume.
            w_cold( vi) = MAX(-w_cutoff, MIN(1._dp + w_cutoff, w_tot ))
            w_warm( vi) = 1._dp - w_cold( vi)
          ELSE
            ! No ice in warm climate, ice in cold climate. Linear inter- / extrapolation.
            w_cold( vi) = MAX(-w_cutoff, MIN(1._dp + w_cutoff, ((ice%Hs_a( vi) - climate%matrix%GCM_PI%Hs( vi)) / (climate%matrix%GCM_cold%Hs( vi) - climate%matrix%GCM_PI%Hs( vi))) * w_tot ))
            w_warm( vi)  = 1._dp - w_cold( vi)
          END IF
        ELSE
          ! Ice in both GCM states.  Linear inter- / extrapolation
          w_cold( vi) = MAX(-w_cutoff, MIN(1._dp + w_cutoff, ((ice%Hs_a( vi) - climate%matrix%GCM_PI%Hs( vi)) / (climate%matrix%GCM_cold%Hs( vi) - climate%matrix%GCM_PI%Hs( vi))) * w_tot ))
          w_warm( vi)  = 1._dp - w_cold( vi)
        END IF

      END DO

      w_cold( mesh%vi1:mesh%vi2) = w_cold( mesh%vi1:mesh%vi2) * w_tot 
      
      ! Smooth the weighting field
      CALL smooth_Gaussian_2D( grid, w_cold, 200000._dp)

      w_warm( mesh%vi1:mesh%vi2) = 1._dp - w_cold( mesh%vi1:mesh%vi2)

    ELSEIF (region_name == 'GRL' .OR. region_name == 'ANT') THEN
      ! Use only total ice volume and CO2; Berends et al., 2018, Eq. 13

      w_cold( mesh%vi1:mesh%vi2) = w_tot 
      w_warm( mesh%vi1:mesh%vi2) = 1._dp - w_cold( mesh%vi1:mesh%vi2)

    END IF
!
!!! CHECK THE ChANGES DONE ABOvE WITH IMAU-ICE, MAKE SENSE AS BEFORE?
!
    IF (C%switch_glacial_index_precip) THEN ! If a glacial index is used for the precipitation forcing, it will only depend on CO2
      w_tot = 1._dp - (MAX( -w_cutoff, MIN( 1._dp + w_cutoff, (forcing%CO2_obs - C%matrix_low_CO2_level) / (C%matrix_high_CO2_level - C%matrix_low_CO2_level) )) )
      w_cold( mesh%vi1:mesh%vi2) = w_tot
      w_warm( mesh%vi1:mesh%vi2) = 1._dp - w_cold( mesh%vi1:mesh%vi2)
    END IF

    ! Interpolate the GCM snapshots
    ! =============================

    DO vi = mesh%vi1, mesh%vi2

      T_ref_GCM( vi,:) =      (w_warm( vi) *     climate%matrix%GCM_warm%T2m(    vi,: )) + & 
                              (w_cold( vi) *     climate%matrix%GCM_cold%T2m(    vi,: ))   ! Berends et al., 2018 - Eq. 6
                              
      P_ref_GCM( vi,:) = EXP( (w_warm( vi) * LOG(climate%matrix%GCM_warm%Precip( vi,:))) + &
                              (w_cold( vi) * LOG(climate%matrix%GCM_cold%Precip( vi,:)))) ! Berends et al., 2018 - Eq. 7
                              
      Hs_GCM(    vi  ) =      (w_warm( vi) *     climate%matrix%GCM_warm%Hs(     vi   )) + &
                              (w_cold( vi) *     climate%matrix%GCM_cold%Hs(     vi   ))   ! Berends et al., 2018 - Eq. 8

    END DO

    ! Downscale precipitation from the coarse-resolution reference
    ! GCM orography to the fine-resolution ice-model orography
    ! ========================================================

    IF (region_name == 'NAM' .OR. region_name == 'EAS') THEN
      ! Use the Roe&Lindzen precipitation model to do this; Berends et al., 2018, Eqs. A3-A7
      CALL adapt_precip_Roe( mesh, Hs_GCM,   T_ref_GCM  , climate%matrix%PD_obs%Wind_LR, climate%matrix%PD_obs%Wind_DU, P_ref_GCM, &
                                   ice%Hs_a, climate%T2m, climate%matrix%PD_obs%Wind_LR, climate%matrix%PD_obs%Wind_DU, climate%Precip)
    ELSEIF (region_name == 'GRL' .OR. region_name == 'ANT') THEN
      ! Use a simpler temperature-based correction; Berends et al., 2018, Eq. 14
      CALL adapt_precip_CC( mesh, ice%Hs_a, Hs_GCM, T_ref_GCM, P_ref_GCM, climate%Precip, region_name)
    END IF

    ! Clean up after yourself
    CALL deallocate_shared( ww_warm)
    CALL deallocate_shared( ww_cold)
    CALL deallocate_shared( wT_ref_GCM)
    CALL deallocate_shared( wP_ref_GCM)
    CALL deallocate_shared( wHs_GCM)

    ! Safety checks
    CALL check_safety_precipitation( climate%Precip)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_climate_model_matrix_precipitation
  SUBROUTINE initialise_climate_matrix( mesh, grid, climate, region_name, mask_noice)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    type(type_mesh),                     intent(in)    :: grid !used to smooth later on, check if grid is called during initialise
    TYPE(type_climate_model),            INTENT(INOUT) :: climate
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)    :: mask_noice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_climate_matrix'
    INTEGER                                            :: vi, m
    LOGICAL                                            :: found_winds_PD_obs, found_winds_PI, found_winds_warm, found_winds_cold

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    allocate( climate%matrix%I_abs(           mesh%vi1:mesh%vi2))
    allocate( climate%matrix%GCM_bias_T2m(    mesh%vi1:mesh%vi2), 12)
    allocate( climate%matrix%GCM_bias_Precip( mesh%vi1:mesh%vi2), 12)
    !CALL allocate_shared_dp_2D(     grid%ny, grid%nx, climate%matrix%I_abs          , climate%matrix%wI_abs          )
    !CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, climate%matrix%GCM_bias_T2m   , climate%matrix%wGCM_bias_T2m   )
    !CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, climate%matrix%GCM_bias_Precip, climate%matrix%wGCM_bias_Precip)

    ! Allocate memory for the regional ERA40 climate and the final applied climate
    CALL allocate_climate_snapshot( mesh, climate%matrix%PD_obs,   name = 'PD_obs'  )
    CALL allocate_climate_snapshot( mesh, climate%matrix%GCM_PI,   name = 'GCM_PI'  )
    CALL allocate_climate_snapshot( mesh, climate%matrix%GCM_warm, name = 'GCM_warm')
    CALL allocate_climate_snapshot( mesh, climate%matrix%GCM_cold, name = 'GCM_cold')

    call read_climate_snapshot( C%filename_PD_obs_climate       , mesh, climate%matrix%PD_obs  )
    call read_climate_snapshot( C%filename_climate_snapshot_PI  , mesh, climate%matrix%GCM_PI  )
    call read_climate_snapshot( C%filename_climate_snapshot_warm, mesh, climate%matrix%GCM_warm)
    call read_climate_snapshot( C%filename_climate_snapshot_cold, mesh, climate%matrix%GCM_cold)

! here appear some checks if wind is included or not in the data
! after that cames the snapshot ocean and ice mask in dev branch, not in UFE1.x
! need to check what it is? maybe is good to implement it. LINE 1215

    ! Calculate spatially variable lapse rate

    ! Use a uniform value for the warm snapshot [this assumes "warm" is actually identical to PI!]
    climate%matrix%GCM_warm%lambda( mesh%vi1:mesh%vi2) = C%constant_lapserate

    IF     (region_name == 'NAM' .OR. region_name == 'EAS') THEN
      CALL initialise_matrix_calc_spatially_variable_lapserate( mesh, grid, climate%matrix%GCM_PI, climate%matrix%GCM_cold)
    ELSEIF (region_name == 'GLR' .OR. region_name == 'ANT') THEN
      climate%matrix%GCM_cold%lambda( mesh%vi1:mesh%vi2) = C%constant_lapserate
      CALL sync
    END IF
!
! the GCM bias is also applied to Hs in Ufe1.x for consistency, should I add it??
!
    ! Calculate GCM bias
    CALL initialise_matrix_calc_GCM_bias( mesh, climate%matrix%GCM_PI, climate%matrix%PD_obs, &
      climate%matrix%GCM_bias_T2m, climate%matrix%GCM_bias_Precip)

    ! Apply bias correction
    IF (C%climate_matrix_biascorrect_warm) CALL initialise_matrix_apply_bias_correction( mesh, climate%matrix%GCM_warm, &
      climate%matrix%GCM_bias_T2m, climate%matrix%GCM_bias_Precip)
    IF (C%climate_matrix_biascorrect_warm) CALL initialise_matrix_apply_bias_correction( mesh, climate%matrix%GCM_cold, &
      climate%matrix%GCM_bias_T2m, climate%matrix%GCM_bias_Precip)

! working on this now
!
    ! Get reference absorbed insolation for the GCM snapshots
    CALL initialise_matrix_calc_absorbed_insolation( grid, climate%matrix%GCM_warm, region_name, mask_noice)
    CALL initialise_matrix_calc_absorbed_insolation( grid, climate%matrix%GCM_cold, region_name, mask_noice)

    ! Initialise applied climate with present-day observations
    DO i = grid%i2, grid%i2
    DO j = 1, grid%ny
    DO m = 1, 12
      climate%T2m(     m,j,i) = climate%matrix%PD_obs%T2m(     m,j,i)
      climate%Precip(  m,j,i) = climate%matrix%PD_obs%Precip(  m,j,i)
      climate%Wind_LR( m,j,i) = climate%matrix%PD_obs%Wind_LR( m,j,i)
      climate%Wind_DU( m,j,i) = climate%matrix%PD_obs%Wind_DU( m,j,i)
    END DO
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_climate_matrix
  
! == Some generally useful tools
! ==============================

  ! Allocate memory for a single climate snapshot
  SUBROUTINE allocate_climate_snapshot( mesh, snapshot, name)
    ! Allocate shared memory for a single climate snapshot

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_climate_snapshot),         INTENT(INOUT) :: snapshot
    CHARACTER(LEN=*),                    INTENT(IN)    :: name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'allocate_climate_snapshot'

    ! Add routine to path
    CALL init_routine( routine_name)

    snapshot%name = name

    allocate( snapshot%Hs( mesh%vi1:mesh%vi2))
    allocate( snapshot%T2m(     mesh%vi1:mesh%vi2, 12))
    allocate( snapshot%Precip(  mesh%vi1:mesh%vi2, 12))
    allocate( snapshot%Wind_WE( mesh%vi1:mesh%vi2, 12))
    allocate( snapshot%Wind_SN( mesh%vi1:mesh%vi2, 12))
    allocate( snapshot%Wind_LR( mesh%vi1:mesh%vi2, 12))
    allocate( snapshot%Wind_DU( mesh%vi1:mesh%vi2, 12))
    !CALL allocate_shared_dp_2D(     grid%ny, grid%nx, snapshot%Hs,             snapshot%wHs            )
    !CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, snapshot%T2m,            snapshot%wT2m           )
    !CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, snapshot%Precip,         snapshot%wPrecip        )
    !CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, snapshot%Wind_WE,        snapshot%wWind_WE       )
    !CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, snapshot%Wind_SN,        snapshot%wWind_SN       )
    !CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, snapshot%Wind_LR,        snapshot%wWind_LR       )
    !CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, snapshot%Wind_DU,        snapshot%wWind_DU       )

! this is not needed anymore, check!they need to be defined in types!
!    CALL allocate_shared_dp_0D(                       snapshot%CO2,            snapshot%wCO2           )
!    CALL allocate_shared_dp_0D(                       snapshot%orbit_time,     snapshot%worbit_time    )
!    CALL allocate_shared_dp_0D(                       snapshot%orbit_ecc,      snapshot%worbit_ecc     )
!    CALL allocate_shared_dp_0D(                       snapshot%orbit_obl,      snapshot%worbit_obl     )
!    CALL allocate_shared_dp_0D(                       snapshot%orbit_pre,      snapshot%worbit_pre     )
!    CALL allocate_shared_dp_0D(                       snapshot%sealevel,       snapshot%wsealevel      )

    allocate( snapshot%lambda( mesh%vi1:mesh%vi2))
    !CALL allocate_shared_dp_2D(     grid%ny, grid%nx, snapshot%lambda,         snapshot%wlambda        )
    
    allocate( snapshot%Q_TOA(  mesh%vi1:mesh%vi2, 12))
    allocate( snapshot%Albedo( mesh%vi1:mesh%vi2, 12))
    allocate( snapshot%I_abs(  mesh%vi1:mesh%vi2))
    !CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, snapshot%Q_TOA,          snapshot%wQ_TOA         )
    !CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, snapshot%Albedo,         snapshot%wAlbedo        )
    !CALL allocate_shared_dp_2D(     grid%ny, grid%nx, snapshot%I_abs,          snapshot%wI_abs         )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE allocate_climate_snapshot
    SUBROUTINE read_climate_snapshot( filename, mesh, snapshot)
    ! Read a climate snapshot from a NetCDF file. Works both for global lon/lat files and regional x/y files

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=256),                 INTENT(IN)    :: filename
    TYPE(type_mesh),                    INTENT(IN)    :: mesh
    TYPE(type_climate_snapshot),        INTENT(INOUT) :: snapshot

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                     :: routine_name = 'read_climate_snapshot'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Write message to screen
    IF (par%master) WRITE(0,*) '  Reading climate for snapshot "' // TRIM( snapshot%name) // '" from file ' // TRIM( filename)

    ! here in IMAU-ICE it check the name variable of the wind. In UFEMISM2.0 the subroutine is called as inquire_var_multopt
    ! to use it would need to have ncid this means open the netcdf file and it will be a lot of code writing that do not look nice
    ! For simplicity now, we will assume that wind is WE, SN as UFE1.x
    ! call inquire_var_multopt( filename, ncid, var_name_options = 'x||X||x1||X1||x-coordinate||X-coordinate||easting', id_var)

! 5th input now is: time_to_read = timeframe_SMB_prescribed. Is optional but think about it
    CALL read_field_from_file_2D(         filename, field_name_options_Hs , mesh, snapshot%Hs     )
    CALL read_field_from_file_2D_monthly( filename, 'T2m'                 , mesh, snapshot%T2m    )
    CALL read_field_from_file_2D_monthly( filename, 'Precip'              , mesh, snapshot%Precip )
    call read_field_from_file_2D_monthly( filename, 'Wind_WE||uas||'      , mesh, snapshot%Wind_WE) ! is needed the last ||? I copy it from SMB_realistic
    call read_field_from_file_2D_monthly( filename, 'Wind_SN||vas||'      , mesh, snapshot%Wind_SN)

    call rotate_wind_to_model_mesh( mesh, snapshot%Wind_WE, snapshot%Wind_SN, snapshot%Wind_LR, snapshot%Wind_DU)
    
    ! option to do safety check here, I moved it after run, so it will check also with the forcings! functions not needed
    !CALL check_safety_temperature(   snapshot%T2m   )
    !CALL check_safety_precipitation( snapshot%Precip)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_climate_snapshot
  
  SUBROUTINE initialise_matrix_calc_GCM_bias( mesh, GCM_PI, PD_obs, GCM_bias_T2m, GCM_bias_Precip)
    ! Calculate the GCM bias in temperature and precipitation
    !
    ! Account for the fact that the GCM PI snapshot has a lower resolution, and therefore
    ! a different surface elevation than the PD observed climatology!

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_climate_snapshot),         INTENT(IN)    :: GCM_PI, PD_obs
    REAL(dp), DIMENSION(:,:),          INTENT(OUT)   :: GCM_bias_T2m
    REAL(dp), DIMENSION(:,:),          INTENT(OUT)   :: GCM_bias_Precip

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_matrix_calc_GCM_bias'
    INTEGER                                            :: vi,m
    REAL(dp)                                           :: T2m_SL_GCM, T2m_SL_obs

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate bias
    DO vi = mesh%vi1, mesh%vi2
    DO m = 1, 12

      ! Scale modelled and observed temperature to sea level using a constant lapse rate
      T2m_SL_GCM = GCM_PI%T2m( vi,m) + GCM_PI%Hs( vi) * C%constant_lapserate
      T2m_SL_obs = PD_obs%T2m( vi,m) + PD_obs%Hs( vi) * C%constant_lapserate

      ! Calculate bias
      GCM_bias_T2m(    vi,m) = T2m_SL_GCM           - T2m_SL_obs
      GCM_bias_Precip( vi,m) = GCM_PI%Precip( vi,m) / PD_obs%Precip( vi,m)

    END DO
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_matrix_calc_GCM_bias
  SUBROUTINE initialise_matrix_apply_bias_correction( mesh, snapshot, bias_T2m, bias_Precip)
    ! Apply a bias correction to this (GCM) snapshot

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                      INTENT(IN)    :: mesh
    TYPE(type_climate_snapshot),          INTENT(INOUT) :: snapshot
    REAL(dp), DIMENSION(:,:),           INTENT(IN)    :: bias_T2m
    REAL(dp), DIMENSION(:,:),           INTENT(IN)    :: bias_Precip

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                       :: routine_name = 'initialise_matrix_apply_bias_correction'
    INTEGER                                             :: vi,m

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Apply bias correction
    DO vi = mesh%vi1, mesh%vi2
    DO m = 1, 12
      snapshot%T2m(    vi,m) = snapshot%T2m(    vi,m) - bias_T2m(    vi,m)
      snapshot%Precip( vi,m) = snapshot%Precip( vi,m) / bias_Precip( vi,m)
    END DO
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_matrix_apply_bias_correction
  SUBROUTINE initialise_matrix_calc_spatially_variable_lapserate( mesh, grid_smooth, snapshot_PI, snapshot)
    ! Calculate the spatially variable lapse-rate (for non-PI GCM climates; see Berends et al., 2018)
    ! Only meaningful for climates where there is ice (LGM, M2_Medium, M2_Large),
    ! and only intended for North America and Eurasia

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                      INTENT(IN)    :: mesh
    TYPE(type_grid),                      INTENT(IN)    :: grid_smooth
    TYPE(type_climate_snapshot),          INTENT(IN)    :: snapshot_PI
    TYPE(type_climate_snapshot),          INTENT(INOUT) :: snapshot

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                       :: routine_name = 'initialise_matrix_calc_spatially_variable_lapserate'
    INTEGER                                             :: i,j,m
    INTEGER,  DIMENSION(:    ), POINTER                 ::  mask_calc_lambda
    REAL(dp)                                            :: dT_mean_nonice
    INTEGER                                             :: n_nonice, n_ice
    REAL(dp)                                            :: lambda_mean_ice

    REAL(dp), PARAMETER                                 :: lambda_min = 0.002_dp
    REAL(dp), PARAMETER                                 :: lambda_max = 0.05_dp

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    allocate( mask_calc_lambda( mesh%vi1:mesh%vi2))
    !CALL allocate_shared_int_2D( grid%ny, grid%nx, mask_calc_lambda, wmask_calc_lambda)

    ! Determine where the variable lapse rate should be calculated
    ! (i.e. where has the surface elevation increased substantially)
    ! ==============================================================

    DO vi = mesh%vi1, mesh%vi2
! this is done in a different way in Ufe1.x
      IF (snapshot%Hs( vi) > snapshot_PI%Hs( vi) + 100._dp) THEN
        mask_calc_lambda( vi) = 1
      ELSE
        mask_calc_lambda( vi) = 0
      END IF

    END DO

    ! Calculate the regional average temperature change outside of the ice sheet
    ! ==========================================================================

    dT_mean_nonice = 0._dp
    n_nonice       = 0
    DO vi = mesh%vi1, mesh%vi2
    DO m = 1, 12
      IF (mask_calc_lambda( vi) == 0) THEN
        dT_mean_nonice = dT_mean_nonice + snapshot%T2m( vi,m) - snapshot_PI%T2m( vi,m)
        n_nonice = n_nonice + 1
      END IF
    END DO
    END DO
    
!! this call still work? or is different now? check
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, dT_mean_nonice, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, n_nonice,       1, MPI_INTEGER,          MPI_SUM, MPI_COMM_WORLD, ierr)

    dT_mean_nonice = dT_mean_nonice / REAL(n_nonice,dp)

    ! Calculate the lapse rate over the ice itself
    ! ============================================

    lambda_mean_ice = 0._dp
    n_ice           = 0

    DO vi = mesh%vi1, mesh%vi2

      IF (mask_calc_lambda( vi) == 1) THEN

        DO m = 1, 12
          ! Berends et al., 2018 - Eq. 10
          snapshot%lambda( vi) = snapshot%lambda( vi) + 1/12._dp * MAX(lambda_min, MIN(lambda_max, &
            -(snapshot%T2m( vi,m) - (snapshot_PI%T2m( vi,m) + dT_mean_nonice)) / (snapshot%Hs( vi) - snapshot_PI%Hs( vi))))
        END DO

        lambda_mean_ice = lambda_mean_ice + snapshot%lambda( vi)
        n_ice = n_ice + 1

      END IF

    END DO
! again check this call later
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, lambda_mean_ice, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, n_ice,           1, MPI_INTEGER,          MPI_SUM, MPI_COMM_WORLD, ierr)

    lambda_mean_ice = lambda_mean_ice / n_ice

    ! Apply mean lapse-rate over ice to the rest of the region
    ! ========================================================

    DO vi = mesh%vi1, mesh%vi2
      IF (mask_calc_lambda( vi) == 0) snapshot%lambda( vi) = lambda_mean_ice
    END DO

    ! Smooth the lapse rate field with a 160 km Gaussian filter
    CALL smooth_Gaussian_2D( mesh, grid_smooth, snapshot%lambda, 160000._dp)

    ! Normalise the entire region to a mean lapse rate of 8 K /km
    snapshot%lambda( mesh%vi1:mesh%vi2) = snapshot%lambda( mesh%vi1:mesh%vi2) * (C%constant_lapserate / lambda_mean_ice)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_matrix_calc_spatially_variable_lapserate
  SUBROUTINE initialise_matrix_calc_absorbed_insolation( mesh, snapshot, region_name, mask_noice)
    ! Calculate the yearly absorbed insolation for this (regional) GCM snapshot, to be used in the matrix interpolation

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                      INTENT(IN)    :: mesh
    TYPE(type_climate_snapshot),          INTENT(INOUT) :: snapshot
    CHARACTER(LEN=3),                     INTENT(IN)    :: region_name
    INTEGER,  DIMENSION(:,:  ),           INTENT(IN)    :: mask_noice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                       :: routine_name = 'initialise_matrix_calc_absorbed_insolation'
    INTEGER                                             :: i,j,m
    TYPE(type_ice_model)                                :: ice_dummy
    TYPE(type_climate_model)                            :: climate_dummy
    TYPE(type_SMB_model)                                :: SMB_dummy

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Get insolation at the desired time from the insolation NetCDF file
    ! ==================================================================

    CALL get_insolation_at_time( mesh, snapshot%orbit_time, snapshot%Q_TOA)

    ! Create temporary "dummy" climate, ice & SMB data structures,
    ! so we can run the SMB model and determine the reference albedo field
    ! ====================================================================

    ! Climate
    ! =======

    ! Allocate shared memory
    allocate( climate_dummy%T2m(    mesh%vi1:mesh%vi2, 12))
    allocate( climate_dummy%Precip( mesh%vi1:mesh%vi2, 12))
    allocate( climate_dummy%Q_TOA(  mesh%vi1:mesh%vi2, 12))
    !CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, climate_dummy%T2m,    climate_dummy%wT2m)
    !CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, climate_dummy%Precip, climate_dummy%wPrecip)
    !CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, climate_dummy%Q_TOA,  climate_dummy%wQ_TOA)

    ! Copy climate fields
    climate_dummy%T2m(    mesh%vi1:mesh%vi2,:) = snapshot%T2m(    mesh%vi1:mesh%vi2,:)
    climate_dummy%Precip( mesh%vi1:mesh%vi2,:) = snapshot%Precip( mesh%vi1:mesh%vi2,:)
    climate_dummy%Q_TOA(  mesh%vi1:mesh%vi2,:) = snapshot%Q_TOA(  mesh%vi1:mesh%vi2,:)

    ! Ice
    ! ===
    allocate( ice_dummy%mask_ocean_a( mesh%vi1:mesh%vi2))
    allocate( ice_dummy%mask_ice_a(   mesh%vi1:mesh%vi2))
    allocate( ice_dummy%mask_shelf_a( mesh%vi1:mesh%vi2))

    ! Fill in masks for the SMB model
    DO vi = mesh%vi1, mesh%vi2

      IF (snapshot%Hs( vi) == MINVAL(snapshot%Hs)) THEN
        ice_dummy%mask_ocean_a( vi) = 1
      ELSE
        ice_dummy%mask_ocean_a( vi) = 0
      END IF

      ! this IF is like (climate%Mask_ice( vi) > .3_dp) in Ufe1.x
      IF (snapshot%Hs( vi) > 100._dp .AND. SUM(snapshot%T2m( vi,:)) / 12._dp < 0._dp) THEN
        ice_dummy%mask_ice_a(   vi) = 1
      ELSE
        ice_dummy%mask_ice_a(   vi) = 0
      END IF

      ! mask_shelf is used in the SMB model only to find open ocean; since mask_ocean
      ! in this case already marks only open ocean, no need to look for shelves
      ice_dummy%mask_shelf_a( vi) = 0

    END DO

    ! SMB
    ! ===
    allocate( SMB_dummy%AlbedoSurf(       mesh%vi1:mesh%vi2)
    allocate( SMB_dummy%MeltPreviousYear( mesh%vi1:mesh%vi2))
    allocate( SMB_dummy%FirnDepth(        mesh%vi1:mesh%vi2, 12))
    allocate( SMB_dummy%Rainfall(         mesh%vi1:mesh%vi2, 12))
    allocate( SMB_dummy%Snowfall(         mesh%vi1:mesh%vi2, 12))
    allocate( SMB_dummy%AddedFirn(        mesh%vi1:mesh%vi2, 12))
    allocate( SMB_dummy%Melt(             mesh%vi1:mesh%vi2, 12))
    allocate( SMB_dummy%Refreezing(       mesh%vi1:mesh%vi2, 12))
    allocate( SMB_dummy%Refreezing_year(  mesh%vi1:mesh%vi2))
    allocate( SMB_dummy%Runoff(           mesh%vi1:mesh%vi2, 12))
    allocate( SMB_dummy%Albedo(           mesh%vi1:mesh%vi2, 12))
    allocate( SMB_dummy%Albedo_year(      mesh%vi1:mesh%vi2))
    allocate( SMB_dummy%SMB(              mesh%vi1:mesh%vi2, 12))
    allocate( SMB_dummy%SMB_year(         mesh%vi1:mesh%vi2))

! not needed anymore.. commment for now
!    CALL allocate_shared_dp_0D( SMB_dummy%C_abl_constant, SMB_dummy%wC_abl_constant)
!    CALL allocate_shared_dp_0D( SMB_dummy%C_abl_Ts,       SMB_dummy%wC_abl_Ts      )
!    CALL allocate_shared_dp_0D( SMB_dummy%C_abl_Q,        SMB_dummy%wC_abl_Q       )
!    CALL allocate_shared_dp_0D( SMB_dummy%C_refr,         SMB_dummy%wC_refr        )

    IF (par%master) THEN
      IF     (region_name == 'NAM') THEN
        SMB_dummy%C_abl_constant = C%SMB_IMAUITM_C_abl_constant_NAM
        SMB_dummy%C_abl_Ts       = C%SMB_IMAUITM_C_abl_Ts_NAM
        SMB_dummy%C_abl_Q        = C%SMB_IMAUITM_C_abl_Q_NAM
        SMB_dummy%C_refr         = C%SMB_IMAUITM_C_refr_NAM
      ELSEIF (region_name == 'EAS') THEN
        SMB_dummy%C_abl_constant = C%SMB_IMAUITM_C_abl_constant_EAS
        SMB_dummy%C_abl_Ts       = C%SMB_IMAUITM_C_abl_Ts_EAS
        SMB_dummy%C_abl_Q        = C%SMB_IMAUITM_C_abl_Q_EAS
        SMB_dummy%C_refr         = C%SMB_IMAUITM_C_refr_EAS
      ELSEIF (region_name == 'GRL') THEN
        SMB_dummy%C_abl_constant = C%SMB_IMAUITM_C_abl_constant_GRL
        SMB_dummy%C_abl_Ts       = C%SMB_IMAUITM_C_abl_Ts_GRL
        SMB_dummy%C_abl_Q        = C%SMB_IMAUITM_C_abl_Q_GRL
        SMB_dummy%C_refr         = C%SMB_IMAUITM_C_refr_GRL
      ELSEIF (region_name == 'ANT') THEN
        SMB_dummy%C_abl_constant = C%SMB_IMAUITM_C_abl_constant_ANT
        SMB_dummy%C_abl_Ts       = C%SMB_IMAUITM_C_abl_Ts_ANT
        SMB_dummy%C_abl_Q        = C%SMB_IMAUITM_C_abl_Q_ANT
        SMB_dummy%C_refr         = C%SMB_IMAUITM_C_refr_ANT
      END IF
    END IF ! IF (par%master) THEN

    ! Run the SMB model for 10 years for this particular climate
    ! (experimentally determined to be long enough to converge)
    DO i = 1, 10
      CALL run_SMB_model( mesh, ice_dummy, climate_dummy, 0._dp, SMB_dummy, mask_noice)
    END DO

    ! Calculate yearly total absorbed insolation
    snapshot%I_abs( mesh%vi1:mesh%vi2,:) = 0._dp
    DO vi = mesh%vi1, mesh%vi2
    DO m = 1, 12
      snapshot%I_abs( vi) = snapshot%I_abs( vi) + snapshot%Q_TOA( vi,m) * (1._dp - SMB_dummy%Albedo( vi,m))
    END DO
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_matrix_calc_absorbed_insolation


  ! Two different parameterised precipitation models:
  ! - a simply Clausius-Clapeyron-based method            (used for GRL and ANT)
  ! - the Roe & Lindzen temperature/orography-based model (used for NAM and EAS)
  SUBROUTINE adapt_precip_CC( mesh, Hs, Hs_GCM, T_ref_GCM, P_ref_GCM, Precip_GCM, region_name)

    USE parameters_module, ONLY: T0

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: Hs              ! Model orography (m)
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: Hs_GCM          ! Reference orography (m)           - total ice-weighted
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: T_ref_GCM       ! Reference temperature (K)         - total ice-weighted
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: P_ref_GCM       ! Reference precipitation (m/month) - total ice-weighted
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Output variables:
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: Precip_GCM      ! Climate matrix precipitation

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'adapt_precip_CC'
    INTEGER                                            :: vi,m
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  T_inv,  T_inv_ref
    INTEGER                                            :: wT_inv, wT_inv_ref

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_dp_2D( mesh%nV, 12, T_inv,     wT_inv    )
    CALL allocate_shared_dp_2D( mesh%nV, 12, T_inv_ref, wT_inv_ref)

    ! Calculate inversion layer temperatures
    DO m = 1, 12
    DO vi = mesh%vi1, mesh%vi2
      T_inv_ref( vi,m) = 88.9_dp + 0.67_dp *  T_ref_GCM( vi,m)
      T_inv(     vi,m) = 88.9_dp + 0.67_dp * (T_ref_GCM( vi,m) - C%constant_lapserate * (Hs( vi) - Hs_GCM( vi)))
    END DO
    END DO
    CALL sync

    IF     (region_name == 'GRL') THEN
      ! Method of Jouzel and Merlivat (1984), see equation (4.82) in Huybrechts (1992)

      DO m = 1, 12
      DO vi = mesh%vi1, mesh%vi2
        Precip_GCM( vi,m) = P_ref_GCM( vi,m) * 1.04**(T_inv( vi,m) - T_inv_ref( vi,m))
      END DO
      END DO
      CALL sync

    ELSEIF (region_name == 'ANT') THEN
      ! As with Lorius/Jouzel method (also Huybrechts, 2002

      DO m = 1, 12
      DO vi = mesh%vi1, mesh%vi2
        Precip_GCM( vi,m) = P_ref_GCM( vi,m) * (T_inv_ref( vi,m) / T_inv( vi,m))**2 * EXP(22.47_dp * (T0 / T_inv_ref( vi,m) - T0 / T_inv( vi,m)))
      END DO
      END DO
      CALL sync

    ELSE
      IF (par%master) WRITE(0,*) '  ERROR - adapt_precip_CC should only be used for Greenland and Antarctica!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF

    ! Clean up after yourself
    CALL deallocate_shared( wT_inv)
    CALL deallocate_shared( wT_inv_ref)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE adapt_precip_CC
  SUBROUTINE adapt_precip_Roe( mesh, Hs1, T2m1, Wind_LR1, Wind_DU1, Precip1, &
                                     Hs2, T2m2, Wind_LR2, Wind_DU2, Precip2)
    ! Adapt precipitation from reference state 1 to model state 2, using the Roe&Lindzen precipitation model

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: Hs1,      Hs2
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: T2m1,     T2m2
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: Wind_LR1, Wind_LR2
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: Wind_DU1, Wind_DU2
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: Precip1
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: Precip2

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'adapt_precip_Roe'
    INTEGER                                            :: vi,m
    REAL(dp), DIMENSION(:    ), POINTER                ::  dHs_dx1,  dHs_dx2
    REAL(dp), DIMENSION(:    ), POINTER                ::  dHs_dy1,  dHs_dy2
    INTEGER                                            :: wdHs_dx1, wdHs_dx2
    INTEGER                                            :: wdHs_dy1, wdHs_dy2
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  Precip_RL1,  Precip_RL2,  dPrecip_RL
    INTEGER                                            :: wPrecip_RL1, wPrecip_RL2, wdPrecip_RL

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nv,     dHs_dx1,     wdHs_dx1   )
    CALL allocate_shared_dp_1D( mesh%nv,     dHs_dx2,     wdHs_dx2   )
    CALL allocate_shared_dp_1D( mesh%nv,     dHs_dy1,     wdHs_dy1   )
    CALL allocate_shared_dp_1D( mesh%nv,     dHs_dy2,     wdHs_dy2   )
    CALL allocate_shared_dp_2D( mesh%nv, 12, Precip_RL1,  wPrecip_RL1)
    CALL allocate_shared_dp_2D( mesh%nv, 12, Precip_RL2,  wPrecip_RL2)
    CALL allocate_shared_dp_2D( mesh%nv, 12, dPrecip_RL,  wdPrecip_RL)

    ! Calculate surface slopes for both states
    CALL ddx_a_to_a_2D( mesh, Hs1, dHs_dx1)
    CALL ddx_a_to_a_2D( mesh, Hs2, dHs_dx2)
    CALL ddy_a_to_a_2D( mesh, Hs1, dHs_dy1)
    CALL ddy_a_to_a_2D( mesh, Hs2, dHs_dy2)

    DO vi = mesh%vi1, mesh%vi2
    DO m = 1, 12

      ! Calculate precipitation with the Roe&Lindzen model for both states
      CALL precipitation_model_Roe( T2m1( vi,m), dHs_dx1( vi), dHs_dy1( vi), Wind_LR1( vi,m), Wind_DU1( vi,m), Precip_RL1( vi,m))
      CALL precipitation_model_Roe( T2m2( vi,m), dHs_dx2( vi), dHs_dy2( vi), Wind_LR2( vi,m), Wind_DU2( vi,m), Precip_RL2( vi,m))

      ! Calculate the ratio between those two precipitation rates
      dPrecip_RL( vi,m) = MAX(0.01_dp, MIN( 2._dp, Precip_RL2( vi,m) / Precip_RL1( vi,m) ))

      ! Applied model precipitation = (matrix-interpolated GCM reference precipitation) * RL ratio
      Precip2( vi,m) = Precip1( vi,m) * dPrecip_RL( vi,m)

    END DO
    END DO
    CALL sync

    ! Clean up after yourself
    CALL deallocate_shared( wdHs_dx1)
    CALL deallocate_shared( wdHs_dx2)
    CALL deallocate_shared( wdHs_dy1)
    CALL deallocate_shared( wdHs_dy2)
    CALL deallocate_shared( wPrecip_RL1)
    CALL deallocate_shared( wPrecip_RL2)
    CALL deallocate_shared( wdPrecip_RL)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE adapt_precip_Roe
  SUBROUTINE precipitation_model_Roe( T2m, dHs_dx, dHs_dy, Wind_LR, Wind_DU, Precip)
    ! Precipitation model of Roe (J. Glac, 2002), integration from Roe and Lindzen (J. Clim. 2001)

    USE parameters_module, ONLY: T0, pi, sec_per_year

    ! In/output variables:
    REAL(dp),                            INTENT(IN)    :: T2m                  ! 2-m air temperature [K]
    REAL(dp),                            INTENT(IN)    :: dHs_dx               ! Surface slope in the x-direction [m/m]
    REAL(dp),                            INTENT(IN)    :: dHs_dy               ! Surface slope in the y-direction [m/m]
    REAL(dp),                            INTENT(IN)    :: Wind_LR              ! Wind speed    in the x-direction [m/s]
    REAL(dp),                            INTENT(IN)    :: Wind_DU              ! Wind speed    in the y-direction [m/s]
    REAL(dp),                            INTENT(OUT)   :: Precip               ! Modelled precipitation

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'precipitation_model_Roe'
    REAL(dp)                                           :: upwind_slope         ! Upwind slope
    REAL(dp)                                           :: E_sat                ! Saturation vapour pressure as function of temperature [Pa]
    REAL(dp)                                           :: x0                   ! Integration parameter x0 [m s-1]
    REAL(dp)                                           :: err_in,err_out

    REAL(dp), PARAMETER                                :: e_sat0  = 611.2_dp   ! Saturation vapour pressure at 273.15 K [Pa]
    REAL(dp), PARAMETER                                :: c_one   = 17.67_dp   ! Constant c1 []
    REAL(dp), PARAMETER                                :: c_two   = 243.5_dp   ! Constant c2 [Celcius]

    REAL(dp), PARAMETER                                :: a_par   = 2.5E-11_dp ! Constant a [m2 s  kg-1] (from Roe et al., J. Clim. 2001)
    REAL(dp), PARAMETER                                :: b_par   = 5.9E-09_dp ! Constant b [m  s2 kg-1] (from Roe et al., J. Clim. 2001)
    REAL(dp), PARAMETER                                :: alpha   = 100.0_dp   ! Constant alpha [s m-1]

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate the upwind slope
    upwind_slope = MAX(0._dp, Wind_LR * dHs_dx + Wind_DU * dHs_dy)

    ! Calculate the saturation vapour pressure E_sat:
    E_sat = e_sat0 * EXP( c_one * (T2m - T0) / (c_two + T2m - T0) )

    ! Calculate integration parameter x0 = a/b + w (with w = wind times slope)
    x0 = a_par / b_par + upwind_slope

    ! Calculate the error function (2nd term on the r.h.s.)
    err_in = alpha * ABS(x0)
    CALL error_function(err_in,err_out)

    ! Calculate precipitation rate as in Appendix of Roe et al. (J. Clim, 2001)
    Precip = ( b_par * E_sat ) * ( x0 / 2._dp + x0**2 * err_out / (2._dp * ABS(x0)) + &
                                         EXP (-alpha**2 * x0**2) / (2._dp * SQRT(pi) * alpha) ) * sec_per_year

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE precipitation_model_Roe
  ! Rotate wind_WE, wind_SN to wind_LR, wind_DU
  SUBROUTINE rotate_wind_to_model_mesh( mesh, wind_WE, wind_SN, wind_LR, wind_DU)
    ! Code copied from ANICE.

    USE parameters, ONLY: pi

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: wind_WE
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: wind_SN
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: wind_LR
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: wind_DU

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'rotate_wind_to_model_mesh'
    INTEGER                                            :: vi,m
    REAL(dp)                                           :: longitude_start, Uwind_x, Uwind_y, Vwind_x, Vwind_y

    ! Add routine to path
    CALL init_routine( routine_name)

    ! First find the first longitude which defines the start of quadrant I:
    longitude_start = mesh%lambda_M - 90._dp

    DO vi = mesh%vi1, mesh%vi2
    DO m = 1, 12

      ! calculate x and y from the zonal wind
      Uwind_x =   wind_WE( vi,m) * SIN((pi/180._dp) * (mesh%lon( vi) - longitude_start))
      Uwind_y = - wind_WE( vi,m) * COS((pi/180._dp) * (mesh%lon( vi) - longitude_start))

      ! calculate x and y from the meridional winds
      Vwind_x =   wind_SN( vi,m) * COS((pi/180._dp) * (mesh%lon( vi) - longitude_start))
      Vwind_y =   wind_SN( vi,m) * SIN((pi/180._dp) * (mesh%lon( vi) - longitude_start))

      ! Sum up wind components
      wind_LR( vi,m) = Uwind_x + Vwind_x   ! winds left to right
      wind_DU( vi,m) = Uwind_y + Vwind_y   ! winds bottom to top

    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE rotate_wind_to_model_mesh
  ! Safety checks
  SUBROUTINE check_safety_temperature( T2m)
    ! Safety checks on a monthly temperature field

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp), DIMENSION(:,:,:),          INTENT(IN)    :: T2m

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_safety_temperature'
    INTEGER                                            :: n1,n2,n3,i1,i2,i,j,k

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Get field size
    n1 = SIZE( T2m,1)
    n2 = SIZE( T2m,2)
    n3 = SIZE( T2m,3)

    ! Parallelisation
    CALL partition_list( n3, par%i, par%n, i1, i2)

    ! Perform safety check
    DO i = i1, i2
    DO j = 1, n2
    DO k = 1, n1

      ! Temperature errors
      IF     (T2m( k,j,i) < 0._dp) THEN
        CALL crash('negative temperature detected at [{int_01},{int_02},{int_03}]', int_01 = k, int_02 = j, int_03 = i)
      ELSEIF (ISNAN( T2m( k,j,i))) THEN
        CALL crash('NaN temperature detected at [{int_01},{int_02},{int_03}]', int_01 = k, int_02 = j, int_03 = i)
      END IF

      ! Temperature warnings
      IF     (T2m( k,j,i) < 150._dp) THEN
        CALL warning('excessively low temperature (< 150 K) detected at [{int_01},{int_02},{int_03}]', int_01 = k, int_02 = j, int_03 = i)
      ELSEIF (T2m( k,j,i) > 350._dp) THEN
        CALL warning('excessively high temperature (> 350 K) detected at [{int_01},{int_02},{int_03}]', int_01 = k, int_02 = j, int_03 = i)
      END IF

    END DO
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_safety_temperature
  SUBROUTINE check_safety_precipitation( Precip)
    ! Safety checks on a monthly precipitation field

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp), DIMENSION(:,:,:),          INTENT(IN)    :: Precip

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_safety_precipitation'
    INTEGER                                            :: n1,n2,n3,i1,i2,i,j,k

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Get field size
    n1 = SIZE( Precip,1)
    n2 = SIZE( Precip,2)
    n3 = SIZE( Precip,3)

    ! Parallelisation
    CALL partition_list( n3, par%i, par%n, i1, i2)

    ! Perform safety check
    DO i = i1, i2
    DO j = 1, n2
    DO k = 1, n1

      ! Precipitation errors
      IF     (Precip( k,j,i) <= 0._dp) THEN
        CALL crash('zero/negative precipitation detected at [{int_01},{int_02},{int_03}]', int_01 = k, int_02 = j, int_03 = i)
      ELSEIF (ISNAN(Precip( k,j,i))) THEN
        CALL crash('NaN precipitation detected at [{int_01},{int_02},{int_03}]', int_01 = k, int_02 = j, int_03 = i)
      END IF

      ! Precipitation warnings
      IF     (Precip( k,j,i) > 10._dp) THEN
        CALL warning('excessively high precipitation (> 10 m/month) detected at [{int_01},{int_02},{int_03}]', int_01 = k, int_02 = j, int_03 = i)
      END IF

    END DO
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_safety_precipitation

end module climate_matrix