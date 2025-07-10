module climate_matrix

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE control_resources_and_error_messaging                  , ONLY: crash, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE parameters
  use mpi_f08, only: MPI_ALLREDUCE, MPI_IN_PLACE, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, MPI_INTEGER
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  USE grid_types                                             , ONLY: type_grid
  USE climate_model_types                                    , ONLY: type_climate_model, type_climate_model_matrix, type_climate_model_snapshot
  USE global_forcing_types                                   , ONLY: type_global_forcing
  use SMB_model_types, only: type_SMB_model
  USE climate_realistic                                      , ONLY: initialise_climate_model_realistic, initialise_insolation_forcing, initialise_CO2_record, get_insolation_at_time, update_CO2_at_model_time
  USE reallocate_mod                                         , ONLY: reallocate_bounds
  use netcdf_io_main
  use mesh_data_smoothing, only: smooth_Gaussian
  use mesh_disc_apply_operators, only: ddx_a_a_2D, ddy_a_a_2D
  use erf_mod, only: error_function
  use SMB_IMAU_ITM, only: run_SMB_model_IMAUITM

 implicit none

  private

  public :: run_climate_model_matrix
  public :: initialise_climate_matrix
  public :: remap_climate_matrix_model

contains

! == Climate matrix
! ===========================

  ! Climate matrix with warm + cold snapshots, forced with CO2 (from record or from inverse routine) from Berends et al., 2018
  ! Generalised for different timeframes, L.B. Stap (2021)
  SUBROUTINE run_climate_model_matrix( mesh, grid, ice, SMB, climate, region_name, time, forcing)
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
    TYPE(type_global_forcing),           INTENT(INOUT) :: forcing

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_climate_model_matrix'
    INTEGER                                            :: vi ,m

    ! Add routine to path
    CALL init_routine( routine_name)

    !IF (par%primary)  WRITE(*,"(A)") '      Running climate matrix model...'

    ! Update insolation forcing at model time
    CALL get_insolation_at_time( mesh, time, climate%snapshot)

    ! Update CO2 forcing at model time
    CALL update_CO2_at_model_time( time, forcing)

    ! Use the (CO2 + absorbed insolation)-based interpolation scheme for temperature
    CALL run_climate_model_matrix_temperature( mesh, grid, ice, SMB, climate, region_name, forcing)

    ! Use the (CO2 + ice-sheet geometry)-based interpolation scheme for precipitation
    CALL run_climate_model_matrix_precipitation( mesh, grid, ice, climate, region_name, forcing)

    ! == Safety checks from UFE1.x
    ! ================
    DO vi = mesh%vi1, mesh%vi2
    DO m = 1, 12
      IF (climate%T2m( vi,m) < 150._dp) THEN
        CALL crash('excessively low temperatures (<150K) detected!')
      ELSEIF (climate%T2m( vi,m) < 0._dp) THEN
        CALL crash('negative temperatures (<0K) detected!')
      ELSEIF (climate%T2m( vi,m) /= climate%T2m( vi,m)) THEN
        CALL crash('NaN temperatures  detected!')
      ELSEIF (climate%Precip( vi,m) <= 0._dp) THEN
        CALL crash('zero/negative precipitation detected!')
      ELSEIF (climate%Precip( vi,m) /= climate%Precip( vi,m)) THEN
        CALL crash('NaN precipitation detected!')
      END IF
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_climate_model_matrix
  SUBROUTINE run_climate_model_matrix_temperature( mesh, grid, ice, SMB, climate, region_name, forcing)
    ! The (CO2 + absorbed insolation)-based matrix interpolation for temperature, from Berends et al. (2018)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_SMB_model),                INTENT(IN)    :: SMB
    TYPE(type_climate_model),            INTENT(INOUT) :: climate
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    type(type_global_forcing), intent(inout) :: forcing

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_climate_model_matrix_temperature'
    INTEGER                                            :: vi ,m
    REAL(dp)                                           :: CO2, w_CO2
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            ::  w_ins,  w_ins_smooth,  w_ice,  w_tot
    REAL(dp)                                           :: w_ins_av, w_ins_denominator
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: T_ref_GCM
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: Hs_GCM, lambda_GCM
    REAL(dp), PARAMETER                                :: w_cutoff = 0.5_dp        ! Crop weights to [-w_cutoff, 1 + w_cutoff]
    REAL(dp), PARAMETER                                :: P_offset = 0.008_dp       ! Normalisation term in precipitation anomaly to avoid divide-by-nearly-zero

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    allocate( w_ins(        mesh%vi1:mesh%vi2))
    allocate( w_ins_smooth( mesh%vi1:mesh%vi2))
    allocate( w_ice(        mesh%vi1:mesh%vi2))
    allocate( w_tot(        mesh%vi1:mesh%vi2))
    allocate( T_ref_GCM(    mesh%vi1:mesh%vi2, 12))
    allocate( Hs_GCM(       mesh%vi1:mesh%vi2))
    allocate( lambda_GCM(   mesh%vi1:mesh%vi2))
    
    !IF (par%primary)  WRITE(*,"(A)") '   Running climate matrix temperature model...'
    
    ! Find CO2 interpolation weight (use either prescribed or modelled CO2)
    ! =====================================================================

    select case (C%choice_matrix_forcing)

    case default
      call crash('unknown choice_forcing_method"' // TRIM(C%choice_matrix_forcing) // '"!')

    case ('CO2_direct')
      CO2 = forcing%CO2_obs

    case ('d18O_inverse_CO2')
      call crash('d18O_inverse_CO2 not implemented yet!')

    case ('d18O_inverse_dT_glob')
      call crash('must only be called with the correct forcing method, check your code!')

    end select
    call sync
    
    ! If CO2 ~= warm snap -> weight is 1. If ~= cold snap -> weight is 0.
    ! Otherwise interpolate. Berends et al., 2018 - Eq. 1
    w_CO2 = MAX( -w_cutoff, MIN( 1._dp + w_cutoff, (CO2 - C%matrix_low_CO2_level) / &
                               (C%matrix_high_CO2_level - C%matrix_low_CO2_level) ))

    ! Find the interpolation weights based on absorbed insolation
    ! ===========================================================

    ! Calculate modelled absorbed insolation

    climate%matrix%I_abs( mesh%vi1:mesh%vi2) = 0._dp
    call sync
    DO vi = mesh%vi1, mesh%vi2
    DO m = 1, 12
      ! Calculate modelled absorbed insolation. Berends et al., 2018 - Eq. 2
      climate%matrix%I_abs( vi) = climate%matrix%I_abs( vi) + & 
                                  climate%snapshot%Q_TOA( vi,m) * (1._dp - SMB%IMAUITM%Albedo( vi, m))  
    END DO
    END DO
    call sync

    ! Calculate "direct" weighting field
    ! Berends et al., 2018 - Eq. 3

    DO vi = mesh%vi1, mesh%vi2
      ! If absorbed insolation ~= warm snap -> weight is 1.
      ! If ~= cold snap -> weight is 0. Otherwise interpolate
      w_ins_denominator = climate%matrix%GCM_warm%I_abs( vi) - climate%matrix%GCM_cold%I_abs( vi)
      if (abs(w_ins_denominator) > 1.0e-10_dp) then
        w_ins( vi) = MAX( -w_cutoff, MIN( 1._dp + w_cutoff, &
                        ( climate%matrix%I_abs( vi) - climate%matrix%GCM_cold%I_abs( vi)) / &
                        ( climate%matrix%GCM_warm%I_abs( vi) - climate%matrix%GCM_cold%I_abs( vi)) ))
      else
      ! If all the inputs are working, this should not be triggered
        call crash('absorbed insolation in warm and cold snap are equal - causing a division by zero')
      end if

    END DO
    call sync
    
    w_ins_av      = MAX( -w_cutoff, MIN( 1._dp + w_cutoff, (SUM( climate%matrix%I_abs         )      - SUM( climate%matrix%GCM_cold%I_abs)     ) / &
                                                           (SUM( climate%matrix%GCM_warm%I_abs)      - SUM( climate%matrix%GCM_cold%I_abs)     ) ))
    ! Smooth the weighting field
    w_ins_smooth( mesh%vi1:mesh%vi2) = w_ins( mesh%vi1:mesh%vi2)

    CALL smooth_Gaussian( mesh, grid, w_ins_smooth, 200000._dp)

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
!! In UFE1.x here are two more options, glacial matrix and glacial index
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
        climate%T2m( vi,:) = T_ref_GCM( vi, m) - lambda_GCM( vi) * (ice%Hs( vi) - Hs_GCM( vi))  ! Berends et al., 2018 - Eq. 11
      END DO

    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_climate_model_matrix_temperature
  SUBROUTINE run_climate_model_matrix_precipitation( mesh, grid, ice, climate, region_name, forcing)
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
    type(type_global_forcing), intent(in) :: forcing
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_climate_model_matrix_precipitation'
    INTEGER                                            :: vi, m
    REAL(dp), DIMENSION(:), ALLOCATABLE                ::  w_warm,  w_cold
    REAL(dp)                                           :: w_tot
    REAL(dp), DIMENSION(:,:), ALLOCATABLE                :: T_ref_GCM, P_ref_GCM
    REAL(dp), DIMENSION(:),   ALLOCATABLE                :: Hs_GCM

    REAL(dp), PARAMETER                                :: w_cutoff = 0.25_dp        ! Crop weights to [-w_cutoff, 1 + w_cutoff]

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    allocate( w_warm(    mesh%vi1:mesh%vi2))
    allocate( w_cold(    mesh%vi1:mesh%vi2))
    allocate( T_ref_GCM( mesh%vi1:mesh%vi2, 12))
    allocate( P_ref_GCM( mesh%vi1:mesh%vi2, 12))
    allocate( Hs_GCM(    mesh%vi1:mesh%vi2))

    !IF (par%primary)  WRITE(*,"(A)") '   Running climate matrix precipitation model...'

    ! Calculate interpolation weights based on ice geometry
    ! =====================================================
    ! First calculate the total ice volume term (second term in the equation)
    w_tot = MAX(-w_cutoff, MIN(1._dp + w_cutoff, &
      (SUM( ice%Hs) - SUM( climate%matrix%GCM_warm%Hs)) / (SUM( climate%matrix%GCM_cold%Hs) - SUM( climate%matrix%GCM_warm%Hs)) ))

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
            w_cold( vi) = MAX(-w_cutoff, MIN(1._dp + w_cutoff, ((ice%Hs( vi) - climate%matrix%GCM_PI%Hs( vi)) / (climate%matrix%GCM_cold%Hs( vi) - climate%matrix%GCM_PI%Hs( vi))) * w_tot ))
            w_warm( vi)  = 1._dp - w_cold( vi)
          END IF
        ELSE
          ! Ice in both GCM states.  Linear inter- / extrapolation
          w_cold( vi) = MAX(-w_cutoff, MIN(1._dp + w_cutoff, ((ice%Hs( vi) - climate%matrix%GCM_PI%Hs( vi)) / (climate%matrix%GCM_cold%Hs( vi) - climate%matrix%GCM_PI%Hs( vi))) * w_tot ))
          w_warm( vi)  = 1._dp - w_cold( vi)
        END IF

      END DO

      w_cold( mesh%vi1:mesh%vi2) = w_cold( mesh%vi1:mesh%vi2) * w_tot 
      
      ! Smooth the weighting field
      CALL smooth_Gaussian( mesh, grid, w_cold, 200000._dp)

      w_warm( mesh%vi1:mesh%vi2) = 1._dp - w_cold( mesh%vi1:mesh%vi2)

    ELSEIF (region_name == 'GRL' .OR. region_name == 'ANT') THEN
      ! Use only total ice volume and CO2; Berends et al., 2018, Eq. 13

      w_cold( mesh%vi1:mesh%vi2) = w_tot 
      w_warm( mesh%vi1:mesh%vi2) = 1._dp - w_cold( mesh%vi1:mesh%vi2)

    END IF

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
                                   ice%Hs, climate%T2m, climate%matrix%PD_obs%Wind_LR, climate%matrix%PD_obs%Wind_DU, climate%Precip)
    ELSEIF (region_name == 'GRL' .OR. region_name == 'ANT') THEN
      ! Use a simpler temperature-based correction; Berends et al., 2018, Eq. 14
      CALL adapt_precip_CC( mesh, ice%Hs, Hs_GCM, T_ref_GCM, P_ref_GCM, climate%Precip, region_name)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_climate_model_matrix_precipitation
  SUBROUTINE initialise_climate_matrix( mesh, grid, ice, climate, region_name, forcing)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    type(type_grid),                     intent(in)    :: grid !used to smooth later on, check if grid is called during initialise
    type(type_ice_model),                intent(in)    :: ice
    TYPE(type_climate_model),            INTENT(INOUT) :: climate
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    TYPE(type_global_forcing),              INTENT(INOUT) :: forcing    

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_climate_matrix'
    INTEGER                                            :: vi, m

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    allocate( climate%matrix%I_abs(           mesh%vi1:mesh%vi2))
    allocate( climate%matrix%GCM_bias_T2m(    mesh%vi1:mesh%vi2, 12))
    allocate( climate%matrix%GCM_bias_Precip( mesh%vi1:mesh%vi2, 12))
    !print *, "size of GCM_bias_T2m called just before allocate climate snapshot", size(climate%matrix%GCM_bias_T2m, dim=1)

    ! Allocate memory for the regional ERA40 climate and the final applied climate
    CALL allocate_climate_snapshot( mesh, climate%matrix%PD_obs,   name = 'PD_obs'  )
    CALL allocate_climate_snapshot( mesh, climate%matrix%GCM_PI,   name = 'GCM_PI'  )
    CALL allocate_climate_snapshot( mesh, climate%matrix%GCM_warm, name = 'GCM_warm')
    CALL allocate_climate_snapshot( mesh, climate%matrix%GCM_cold, name = 'GCM_cold')

    call read_climate_snapshot( C%filename_PD_obs_climate       , mesh, climate%matrix%PD_obs  )
    call read_climate_snapshot( C%filename_climate_snapshot_PI  , mesh, climate%matrix%GCM_PI  )
    call read_climate_snapshot( C%filename_climate_snapshot_warm, mesh, climate%matrix%GCM_warm)
    call read_climate_snapshot( C%filename_climate_snapshot_cold, mesh, climate%matrix%GCM_cold)
    
    ! Get the orbit time
    climate%matrix%GCM_PI%orbit_time   = 0._dp
    climate%matrix%GCM_warm%orbit_time = C%matrix_warm_orbit_time
    climate%matrix%GCM_cold%orbit_time = C%matrix_cold_orbit_time

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

    ! Get reference absorbed insolation for the GCM snapshots
    CALL initialise_matrix_calc_absorbed_insolation( mesh, climate%matrix%GCM_warm, region_name, forcing, ice)
    CALL initialise_matrix_calc_absorbed_insolation( mesh, climate%matrix%GCM_cold, region_name, forcing, ice)

   ! initialise the insolation forcing
   call initialise_insolation_forcing( climate%snapshot, mesh) ! this will initialise climate%snapshot%Q_TOA

   ! initialise CO2 forcing
   call initialise_CO2_record( forcing)

    ! Initialise applied climate with present-day observations

    DO vi = mesh%vi1, mesh%vi2
    DO m = 1, 12
      climate%T2m(     vi,m) = climate%matrix%PD_obs%T2m(     vi,m)
      climate%Precip(  vi,m) = climate%matrix%PD_obs%Precip(  vi,m)
      climate%Wind_LR( vi,m) = climate%matrix%PD_obs%Wind_LR( vi,m)
      climate%Wind_DU( vi,m) = climate%matrix%PD_obs%Wind_DU( vi,m)
    END DO
    END DO

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
    TYPE(type_climate_model_snapshot),   INTENT(INOUT) :: snapshot
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

    allocate( snapshot%lambda( mesh%vi1:mesh%vi2))
    ! these are now allocated in initalise_insolation_forcing, commented for now?
    !allocate( snapshot%Q_TOA(  mesh%vi1:mesh%vi2, 12))
    !allocate( snapshot%Albedo( mesh%vi1:mesh%vi2, 12))
    !allocate( snapshot%I_abs(  mesh%vi1:mesh%vi2))
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE allocate_climate_snapshot
    SUBROUTINE read_climate_snapshot( filename, mesh, snapshot)
    ! Read a climate snapshot from a NetCDF file. Works both for global lon/lat files and regional x/y files

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=256),                 INTENT(IN)    :: filename
    TYPE(type_mesh),                    INTENT(IN)    :: mesh
    TYPE(type_climate_model_snapshot),        INTENT(INOUT) :: snapshot

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                     :: routine_name = 'read_climate_snapshot'
    integer                                            :: vi, m
    !REAL(dp)                                           :: longitude_start, Uwind_x, Uwind_y, Vwind_x, Vwind_y
    ! Add routine to path
    CALL init_routine( routine_name)

    ! Write message to screen
    IF (par%primary) WRITE(0,*) '  Reading climate for snapshot "' // TRIM( snapshot%name) // '" from file ' // TRIM( filename)

    CALL read_field_from_file_2D(         filename, field_name_options_Hs , mesh, snapshot%Hs     )
    CALL read_field_from_file_2D_monthly( filename, 'T2m'                 , mesh, snapshot%T2m    )
    CALL read_field_from_file_2D_monthly( filename, 'Precip'              , mesh, snapshot%Precip )
    call read_field_from_file_2D_monthly( filename, 'Wind_WE||uas||'      , mesh, snapshot%Wind_WE) ! is needed the last ||? I copy it from SMB_realistic
    call read_field_from_file_2D_monthly( filename, 'Wind_SN||vas||'      , mesh, snapshot%Wind_SN)
    
    call rotate_wind_to_model_mesh( mesh, snapshot%Wind_WE, snapshot%Wind_SN, snapshot%Wind_LR, snapshot%Wind_DU)

! Check if the snapshot have negative values in the Precip after mapping to mesh
  do m = 1, 12
  do vi = mesh%vi1, mesh%vi2
    if (snapshot%Precip(vi, m) < 0) then
      call crash('snapshot%Precip(vi,m) with negative values')
    end if
  end do
  end do

    !print *, "no error?"
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
    TYPE(type_climate_model_snapshot),   INTENT(IN)    :: GCM_PI, PD_obs 
    REAL(dp), dimension(mesh%vi1:mesh%vi2, 12),          INTENT(OUT)   :: GCM_bias_T2m
    REAL(dp), dimension(mesh%vi1:mesh%vi2, 12),          INTENT(OUT)   :: GCM_bias_Precip

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
    TYPE(type_climate_model_snapshot),          INTENT(INOUT) :: snapshot
    REAL(dp), dimension(mesh%vi1:mesh%vi2, 12),           INTENT(IN)    :: bias_T2m
    REAL(dp), dimension(mesh%vi1:mesh%vi2, 12),           INTENT(IN)    :: bias_Precip

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
    TYPE(type_climate_model_snapshot),          INTENT(IN)    :: snapshot_PI
    TYPE(type_climate_model_snapshot),          INTENT(INOUT) :: snapshot

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                       :: routine_name = 'initialise_matrix_calc_spatially_variable_lapserate'
    INTEGER                                             :: vi,m
    INTEGER,  DIMENSION(:    ), ALLOCATABLE             ::  mask_calc_lambda
    REAL(dp)                                            :: dT_mean_nonice
    INTEGER                                             :: n_nonice, n_ice
    REAL(dp)                                            :: lambda_mean_ice
    integer                               :: ierr

    REAL(dp), PARAMETER                                 :: lambda_min = 0.002_dp
    REAL(dp), PARAMETER                                 :: lambda_max = 0.05_dp

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    allocate( mask_calc_lambda( mesh%vi1:mesh%vi2))

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

    CALL MPI_ALLREDUCE( MPI_IN_PLACE, lambda_mean_ice, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, n_ice,           1, MPI_INTEGER,          MPI_SUM, MPI_COMM_WORLD, ierr)

    lambda_mean_ice = lambda_mean_ice / n_ice

    ! Apply mean lapse-rate over ice to the rest of the region
    ! ========================================================

    DO vi = mesh%vi1, mesh%vi2
      IF (mask_calc_lambda( vi) == 0) snapshot%lambda( vi) = lambda_mean_ice
    END DO

    ! Smooth the lapse rate field with a 160 km Gaussian filter
    CALL smooth_Gaussian( mesh, grid_smooth, snapshot%lambda, 160000._dp)

    ! Normalise the entire region to a mean lapse rate of 8 K /km
    snapshot%lambda( mesh%vi1:mesh%vi2) = snapshot%lambda( mesh%vi1:mesh%vi2) * (C%constant_lapserate / lambda_mean_ice)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_matrix_calc_spatially_variable_lapserate
  SUBROUTINE initialise_matrix_calc_absorbed_insolation( mesh, snapshot, region_name, forcing, ice)
    ! Calculate the yearly absorbed insolation for this (regional) GCM snapshot, to be used in the matrix interpolation

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                      INTENT(IN)    :: mesh
    TYPE(type_climate_model_snapshot),    INTENT(INOUT) :: snapshot
    CHARACTER(LEN=3),                     INTENT(IN)    :: region_name
    TYPE(type_global_forcing),              INTENT(INOUT) :: forcing
    type(type_ice_model),                   intent(in)    :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                       :: routine_name = 'initialise_matrix_calc_absorbed_insolation'
    INTEGER                                             :: vi,m,i
    TYPE(type_ice_model)                                :: ice_dummy
    TYPE(type_climate_model)                            :: climate_dummy
    TYPE(type_SMB_model)                                :: SMB_dummy
    CHARACTER(LEN=256)                                  :: choice_SMB_IMAUITM_init_firn_dummy

    ! Add routine to path
    CALL init_routine( routine_name)
    ! Initialise the insolation variables inside snapshot
    call initialise_insolation_forcing(snapshot, mesh)
    ! Get insolation at the desired time from the insolation NetCDF file
    ! ==================================================================

    CALL get_insolation_at_time( mesh, snapshot%orbit_time, snapshot)

    ! Create temporary "dummy" climate, ice & SMB data structures,
    ! so we can run the SMB model and determine the reference albedo field
    ! ====================================================================

    ! Climate
    ! =======

    ! Allocate shared memory
    allocate( climate_dummy%T2m(    mesh%vi1:mesh%vi2, 12))
    allocate( climate_dummy%Precip( mesh%vi1:mesh%vi2, 12))
    allocate( climate_dummy%snapshot%Q_TOA(  mesh%vi1:mesh%vi2, 12))
    !CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, climate_dummy%T2m,    climate_dummy%wT2m)
    !CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, climate_dummy%Precip, climate_dummy%wPrecip)
    !CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, climate_dummy%Q_TOA,  climate_dummy%wQ_TOA)

    ! Copy climate fields
    climate_dummy%T2m(    mesh%vi1:mesh%vi2,:) = snapshot%T2m(    mesh%vi1:mesh%vi2,:)
    climate_dummy%Precip( mesh%vi1:mesh%vi2,:) = snapshot%Precip( mesh%vi1:mesh%vi2,:)
    ! is needed to allocate it as climate%snapshot because is used in that way later on SMB-ITM
    climate_dummy%snapshot%Q_TOA(  mesh%vi1:mesh%vi2,:) = snapshot%Q_TOA(  mesh%vi1:mesh%vi2,:)

    ! Ice
    ! ===
    allocate( ice_dummy%mask_icefree_ocean( mesh%vi1:mesh%vi2))
    allocate( ice_dummy%mask_grounded_ice(   mesh%vi1:mesh%vi2))
    allocate( ice_dummy%mask_floating_ice( mesh%vi1:mesh%vi2))
    allocate( ice_dummy%mask_noice(        mesh%vi1:mesh%vi2))

    ! Fill in masks for the SMB model
    DO vi = mesh%vi1, mesh%vi2
    
   ! In IMAU-ICE SMB it uses region%mask_noice in UFE2 is ice%mask_noice, I will keep the masks from above for ice_dummy
   ! and make ice_dummy%mask_noice = ice%mask_noice to run the SMB using the dummy, following IMAU-ICE code..
      ice_dummy%mask_noice( vi) = ice%mask_noice( vi) 

      IF (snapshot%Hs( vi) == MINVAL(snapshot%Hs)) THEN
        ice_dummy%mask_icefree_ocean( vi) = .true.
      ELSE
        ice_dummy%mask_icefree_ocean( vi) = .false.
      END IF

      ! this IF is like (climate%Mask_ice( vi) > .3_dp) in Ufe1.x
      IF (snapshot%Hs( vi) > 100._dp .AND. SUM(snapshot%T2m( vi,:)) / 12._dp < 0._dp) THEN
        ice_dummy%mask_grounded_ice(   vi) = .true.
      ELSE
        ice_dummy%mask_grounded_ice(   vi) = .false.
      END IF

      ! mask_shelf is used in the SMB model only to find open ocean; since mask_ocean
      ! in this case already marks only open ocean, no need to look for shelves
      ice_dummy%mask_floating_ice( vi) = .false.

    END DO

    ! SMB
    ! ===
    allocate( SMB_dummy%IMAUITM%AlbedoSurf      (mesh%vi1:mesh%vi2))
    allocate( SMB_dummy%IMAUITM%Rainfall        (mesh%vi1:mesh%vi2, 12))
    allocate( SMB_dummy%IMAUITM%Snowfall        (mesh%vi1:mesh%vi2, 12))
    allocate( SMB_dummy%IMAUITM%AddedFirn       (mesh%vi1:mesh%vi2, 12))
    allocate( SMB_dummy%IMAUITM%Melt            (mesh%vi1:mesh%vi2, 12))
    allocate( SMB_dummy%IMAUITM%Refreezing      (mesh%vi1:mesh%vi2, 12))
    allocate( SMB_dummy%IMAUITM%Refreezing_year (mesh%vi1:mesh%vi2))
    allocate( SMB_dummy%IMAUITM%Runoff          (mesh%vi1:mesh%vi2, 12))
    allocate( SMB_dummy%IMAUITM%Albedo          (mesh%vi1:mesh%vi2, 12))
    allocate( SMB_dummy%IMAUITM%Albedo_year     (mesh%vi1:mesh%vi2))
    allocate( SMB_dummy%IMAUITM%SMB_monthly     (mesh%vi1:mesh%vi2,12))
    allocate( SMB_dummy%IMAUITM%FirnDepth       (mesh%vi1:mesh%vi2,12))
    allocate( SMB_dummy%IMAUITM%MeltPreviousYear(mesh%vi1:mesh%vi2))
    allocate( SMB_dummy%SMB             (mesh%vi1:mesh%vi2))
    SMB_dummy%SMB = 0._dp

    IF (par%primary) THEN
      IF     (region_name == 'NAM') THEN
        SMB_dummy%IMAUITM%C_abl_constant = C%SMB_IMAUITM_C_abl_constant_NAM
        SMB_dummy%IMAUITM%C_abl_Ts       = C%SMB_IMAUITM_C_abl_Ts_NAM
        SMB_dummy%IMAUITM%C_abl_Q        = C%SMB_IMAUITM_C_abl_Q_NAM
        SMB_dummy%IMAUITM%C_refr         = C%SMB_IMAUITM_C_refr_NAM
        choice_SMB_IMAUITM_init_firn_dummy = C%choice_SMB_IMAUITM_init_firn_NAM
      ELSEIF (region_name == 'EAS') THEN
        SMB_dummy%IMAUITM%C_abl_constant = C%SMB_IMAUITM_C_abl_constant_EAS
        SMB_dummy%IMAUITM%C_abl_Ts       = C%SMB_IMAUITM_C_abl_Ts_EAS
        SMB_dummy%IMAUITM%C_abl_Q        = C%SMB_IMAUITM_C_abl_Q_EAS
        SMB_dummy%IMAUITM%C_refr         = C%SMB_IMAUITM_C_refr_EAS
        choice_SMB_IMAUITM_init_firn_dummy = C%choice_SMB_IMAUITM_init_firn_EAS
      ELSEIF (region_name == 'GRL') THEN
        SMB_dummy%IMAUITM%C_abl_constant = C%SMB_IMAUITM_C_abl_constant_GRL
        SMB_dummy%IMAUITM%C_abl_Ts       = C%SMB_IMAUITM_C_abl_Ts_GRL
        SMB_dummy%IMAUITM%C_abl_Q        = C%SMB_IMAUITM_C_abl_Q_GRL
        SMB_dummy%IMAUITM%C_refr         = C%SMB_IMAUITM_C_refr_GRL
        choice_SMB_IMAUITM_init_firn_dummy = C%choice_SMB_IMAUITM_init_firn_GRL
      ELSEIF (region_name == 'ANT') THEN
        SMB_dummy%IMAUITM%C_abl_constant = C%SMB_IMAUITM_C_abl_constant_ANT
        SMB_dummy%IMAUITM%C_abl_Ts       = C%SMB_IMAUITM_C_abl_Ts_ANT
        SMB_dummy%IMAUITM%C_abl_Q        = C%SMB_IMAUITM_C_abl_Q_ANT
        SMB_dummy%IMAUITM%C_refr         = C%SMB_IMAUITM_C_refr_ANT
        choice_SMB_IMAUITM_init_firn_dummy = C%choice_SMB_IMAUITM_init_firn_ANT
      END IF
    END IF ! IF (par%primary) THEN

    IF (par%primary) THEN
      IF     (choice_SMB_IMAUITM_init_firn_dummy == 'uniform') THEN
        ! do nothing
      ELSE
        CALL crash('climate matrix only implemented with uniform init firn"' // TRIM( choice_SMB_IMAUITM_init_firn_dummy) // '"!')
      END IF
    END IF ! IF (par%primary) THEN

    ! Initialise with a uniform firn layer over the ice sheet
    DO vi = mesh%vi1, mesh%vi2
      IF (ice%Hi( vi) > 0._dp) THEN
        SMB_dummy%IMAUITM%FirnDepth(        vi,:) = C%SMB_IMAUITM_initial_firn_thickness
        SMB_dummy%IMAUITM%MeltPreviousYear(   vi) = 0._dp
      ELSE
        SMB_dummy%IMAUITM%FirnDepth(        vi,:) = 0._dp
        SMB_dummy%IMAUITM%MeltPreviousYear(   vi) = 0._dp
      END IF
    END DO

    ! Initialise albedo
    DO vi = mesh%vi1, mesh%vi2
      ! Background albedo
      IF (ice%Hb( vi) < 0._dp) THEN
        SMB_dummy%IMAUITM%AlbedoSurf( vi) = 0.1_dp ! albedo_water
      ELSE
        SMB_dummy%IMAUITM%AlbedoSurf( vi) = 0.2_dp ! albedo_soil
      END IF
      IF (ice%Hi( vi) > 0._dp) THEN
        SMB_dummy%IMAUITM%AlbedoSurf(  vi) = 0.85_dp ! albedo_snow
      END IF
      SMB_dummy%IMAUITM%Albedo( vi,:) = SMB_dummy%IMAUITM%AlbedoSurf( vi)
    END DO

    ! Run the SMB model for 10 years for this particular climate
    ! (experimentally determined to be long enough to converge)
    DO i = 1, 10
      CALL run_SMB_model_IMAUITM( mesh, ice_dummy, SMB_dummy, climate_dummy)
    END DO

    ! Calculate yearly total absorbed insolation
    snapshot%I_abs( mesh%vi1:mesh%vi2) = 0._dp
    DO vi = mesh%vi1, mesh%vi2
    DO m = 1, 12
      snapshot%I_abs( vi) = snapshot%I_abs( vi) + snapshot%Q_TOA( vi,m) * (1._dp - SMB_dummy%IMAUITM%Albedo( vi,m))
    END DO
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_matrix_calc_absorbed_insolation


  ! Two different parameterised precipitation models:
  ! - a simply Clausius-Clapeyron-based method            (used for GRL and ANT)
  ! - the Roe & Lindzen temperature/orography-based model (used for NAM and EAS)
  SUBROUTINE adapt_precip_CC( mesh, Hs, Hs_GCM, T_ref_GCM, P_ref_GCM, Precip_GCM, region_name)

    USE parameters, ONLY: T0

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2),          INTENT(IN)    :: Hs              ! Model orography (m)
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2),          INTENT(IN)    :: Hs_GCM          ! Reference orography (m)           - total ice-weighted
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2,12),          INTENT(IN)    :: T_ref_GCM       ! Reference temperature (K)         - total ice-weighted
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2,12),          INTENT(IN)    :: P_ref_GCM       ! Reference precipitation (m/month) - total ice-weighted
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Output variables:
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2,12),          INTENT(INOUT)   :: Precip_GCM      ! Climate matrix precipitation

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'adapt_precip_CC'
    INTEGER                                            :: vi,m
    REAL(dp), DIMENSION(:, :), ALLOCATABLE            ::  T_inv,  T_inv_ref

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    allocate( T_inv(     mesh%vi1:mesh%vi2, 12))
    allocate( T_inv_ref( mesh%vi1:mesh%vi2, 12))

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
      IF (par%primary) THEN
        CALL crash('ERROR - adapt_precip_CC should only be used for Greenland and Antarctica!')
      END IF
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE adapt_precip_CC
  SUBROUTINE adapt_precip_Roe( mesh, Hs1, T2m1, Wind_LR1, Wind_DU1, Precip1, &
                                     Hs2, T2m2, Wind_LR2, Wind_DU2, Precip2)
    ! Adapt precipitation from reference state 1 to model state 2, using the Roe&Lindzen precipitation model

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2),          INTENT(IN)    :: Hs1,      Hs2
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2, 12),          INTENT(IN)    :: T2m1,     T2m2
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2, 12),          INTENT(IN)    :: Wind_LR1, Wind_LR2
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2, 12),          INTENT(IN)    :: Wind_DU1, Wind_DU2
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2, 12),          INTENT(IN)    :: Precip1
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2, 12),          INTENT(OUT)   :: Precip2

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'adapt_precip_Roe'
    INTEGER                                            :: vi,m
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            ::  dHs_dx1,  dHs_dx2
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            ::  dHs_dy1,  dHs_dy2
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            ::  Precip_RL1,  Precip_RL2,  dPrecip_RL

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    allocate( dHs_dx1( mesh%vi1:mesh%vi2))
    allocate( dHs_dx2( mesh%vi1:mesh%vi2))
    allocate( dHs_dy1( mesh%vi1:mesh%vi2))
    allocate( dHs_dy2( mesh%vi1:mesh%vi2))
    allocate( Precip_RL1( mesh%vi1:mesh%vi2, 12))
    allocate( Precip_RL2( mesh%vi1:mesh%vi2, 12))
    allocate( dPrecip_RL( mesh%vi1:mesh%vi2, 12))

    ! Calculate surface slopes for both states
    CALL ddx_a_a_2D( mesh, Hs1, dHs_dx1)
    CALL ddx_a_a_2D( mesh, Hs2, dHs_dx2)
    CALL ddy_a_a_2D( mesh, Hs1, dHs_dy1)
    CALL ddy_a_a_2D( mesh, Hs2, dHs_dy2)

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

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE adapt_precip_Roe
  SUBROUTINE precipitation_model_Roe( T2m, dHs_dx, dHs_dy, Wind_LR, Wind_DU, Precip)
    ! Precipitation model of Roe (J. Glac, 2002), integration from Roe and Lindzen (J. Clim. 2001)

    USE parameters, ONLY: T0, pi, sec_per_year

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
    REAL(dp)                                           :: err

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
    err = alpha * ABS(x0)
    !CALL error_function(err)

    ! Calculate precipitation rate as in Appendix of Roe et al. (J. Clim, 2001)
    Precip = ( b_par * E_sat ) * ( x0 / 2._dp + x0**2 * error_function(err) / (2._dp * ABS(x0)) + &
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
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2,12),          INTENT(IN)    :: wind_WE
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2,12),          INTENT(IN)    :: wind_SN
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2,12),          INTENT(OUT)   :: wind_LR
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2,12),          INTENT(OUT)   :: wind_DU

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

    call sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE rotate_wind_to_model_mesh

  SUBROUTINE remap_climate_matrix_model( mesh_new, climate, region_name, grid, ice, forcing)

  implicit none

  ! In- and output variables
    !TYPE(type_mesh),                        INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                        INTENT(IN)    :: mesh_new
    TYPE(type_climate_model),               INTENT(INOUT) :: climate
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name
    type(type_grid),                    intent(in)    :: grid
    type(type_ice_model),               intent(in)    :: ice
    type(type_global_forcing),          intent(inout) :: forcing
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'remap_climate_matrix_model'

  IF (par%primary)  WRITE(*,"(A)") '      Remapping climate matrix model data to the new mesh...'
  
  ! reallocate main variables of GCM snapshots
  call reallocate_bounds(climate%matrix%PD_obs%Precip, mesh_new%vi1, mesh_new%vi2, 12)
  call reallocate_bounds(climate%matrix%PD_obs%T2m, mesh_new%vi1, mesh_new%vi2, 12)
  call reallocate_bounds(climate%matrix%PD_obs%Wind_WE, mesh_new%vi1, mesh_new%vi2, 12)
  call reallocate_bounds(climate%matrix%PD_obs%Wind_SN, mesh_new%vi1, mesh_new%vi2, 12)
  call reallocate_bounds(climate%matrix%PD_obs%Wind_LR, mesh_new%vi1, mesh_new%vi2, 12)
  call reallocate_bounds(climate%matrix%PD_obs%Wind_DU, mesh_new%vi1, mesh_new%vi2, 12)
  call reallocate_bounds(climate%matrix%PD_obs%Hs, mesh_new%vi1, mesh_new%vi2)
  call reallocate_bounds(climate%matrix%PD_obs%lambda, mesh_new%vi1, mesh_new%vi2)

  call reallocate_bounds(climate%matrix%GCM_PI%Precip, mesh_new%vi1, mesh_new%vi2, 12)
  call reallocate_bounds(climate%matrix%GCM_PI%T2m, mesh_new%vi1, mesh_new%vi2, 12)
  call reallocate_bounds(climate%matrix%GCM_PI%Wind_WE, mesh_new%vi1, mesh_new%vi2, 12)
  call reallocate_bounds(climate%matrix%GCM_PI%Wind_SN, mesh_new%vi1, mesh_new%vi2, 12)
  call reallocate_bounds(climate%matrix%GCM_PI%Wind_LR, mesh_new%vi1, mesh_new%vi2, 12)
  call reallocate_bounds(climate%matrix%GCM_PI%Wind_DU, mesh_new%vi1, mesh_new%vi2, 12)
  call reallocate_bounds(climate%matrix%GCM_PI%Hs, mesh_new%vi1, mesh_new%vi2)
  call reallocate_bounds(climate%matrix%GCM_PI%lambda, mesh_new%vi1, mesh_new%vi2)

  call reallocate_bounds(climate%matrix%GCM_warm%Precip, mesh_new%vi1, mesh_new%vi2, 12)
  call reallocate_bounds(climate%matrix%GCM_warm%T2m, mesh_new%vi1, mesh_new%vi2, 12)
  call reallocate_bounds(climate%matrix%GCM_warm%Wind_WE, mesh_new%vi1, mesh_new%vi2, 12)
  call reallocate_bounds(climate%matrix%GCM_warm%Wind_SN, mesh_new%vi1, mesh_new%vi2, 12)
  call reallocate_bounds(climate%matrix%GCM_warm%Wind_LR, mesh_new%vi1, mesh_new%vi2, 12)
  call reallocate_bounds(climate%matrix%GCM_warm%Wind_DU, mesh_new%vi1, mesh_new%vi2, 12)
  call reallocate_bounds(climate%matrix%GCM_warm%Hs, mesh_new%vi1, mesh_new%vi2)
  call reallocate_bounds(climate%matrix%GCM_warm%lambda, mesh_new%vi1, mesh_new%vi2)

  call reallocate_bounds(climate%matrix%GCM_cold%Precip, mesh_new%vi1, mesh_new%vi2, 12)
  call reallocate_bounds(climate%matrix%GCM_cold%T2m, mesh_new%vi1, mesh_new%vi2, 12)
  call reallocate_bounds(climate%matrix%GCM_cold%Wind_WE, mesh_new%vi1, mesh_new%vi2, 12)
  call reallocate_bounds(climate%matrix%GCM_cold%Wind_SN, mesh_new%vi1, mesh_new%vi2, 12)
  call reallocate_bounds(climate%matrix%GCM_cold%Wind_LR, mesh_new%vi1, mesh_new%vi2, 12)
  call reallocate_bounds(climate%matrix%GCM_cold%Wind_DU, mesh_new%vi1, mesh_new%vi2, 12)
  call reallocate_bounds(climate%matrix%GCM_cold%Hs, mesh_new%vi1, mesh_new%vi2)
  call reallocate_bounds(climate%matrix%GCM_cold%lambda, mesh_new%vi1, mesh_new%vi2) 

  ! reallocate main variables of climate%snapshot for insolation
  call reallocate_bounds( climate%snapshot%ins_Q_TOA0, mesh_new%vi1, mesh_new%vi2,12)
  call reallocate_bounds( climate%snapshot%ins_Q_TOA1, mesh_new%vi1, mesh_new%vi2,12)
  call reallocate_bounds( climate%snapshot%Q_TOA, mesh_new%vi1, mesh_new%vi2,12)
  
  call reallocate_bounds(climate%matrix%I_abs, mesh_new%vi1, mesh_new%vi2)
  call reallocate_bounds(climate%matrix%GCM_bias_T2m, mesh_new%vi1, mesh_new%vi2, 12)
  call reallocate_bounds(climate%matrix%GCM_bias_Precip, mesh_new%vi1, mesh_new%vi2, 12)
  
  ! read the snapshots for the new mesh
    call read_climate_snapshot( C%filename_PD_obs_climate       , mesh_new, climate%matrix%PD_obs  )
    call read_climate_snapshot( C%filename_climate_snapshot_PI  , mesh_new, climate%matrix%GCM_PI  )
    call read_climate_snapshot( C%filename_climate_snapshot_warm, mesh_new, climate%matrix%GCM_warm)
    call read_climate_snapshot( C%filename_climate_snapshot_cold, mesh_new, climate%matrix%GCM_cold)

    ! Use a uniform value for the warm snapshot [this assumes "warm" is actually identical to PI!]
    climate%matrix%GCM_warm%lambda( mesh_new%vi1:mesh_new%vi2) = C%constant_lapserate

    IF     (region_name == 'NAM' .OR. region_name == 'EAS') THEN
      CALL initialise_matrix_calc_spatially_variable_lapserate( mesh_new, grid, climate%matrix%GCM_PI, climate%matrix%GCM_cold)
    ELSEIF (region_name == 'GLR' .OR. region_name == 'ANT') THEN
      climate%matrix%GCM_cold%lambda( mesh_new%vi1:mesh_new%vi2) = C%constant_lapserate
      CALL sync
    END IF

    ! Calculate GCM bias
    CALL initialise_matrix_calc_GCM_bias( mesh_new, climate%matrix%GCM_PI, climate%matrix%PD_obs, &
      climate%matrix%GCM_bias_T2m, climate%matrix%GCM_bias_Precip)

    ! Apply bias correction
    IF (C%climate_matrix_biascorrect_warm) CALL initialise_matrix_apply_bias_correction( mesh_new, climate%matrix%GCM_warm, &
      climate%matrix%GCM_bias_T2m, climate%matrix%GCM_bias_Precip)
    IF (C%climate_matrix_biascorrect_warm) CALL initialise_matrix_apply_bias_correction( mesh_new, climate%matrix%GCM_cold, &
      climate%matrix%GCM_bias_T2m, climate%matrix%GCM_bias_Precip)

    ! deallocate variables of the insolation from the warm and cold snapshots
    deallocate(climate%matrix%GCM_warm%ins_t0)
    deallocate(climate%matrix%GCM_warm%ins_t1)
    deallocate(climate%matrix%GCM_warm%ins_ti0)
    deallocate(climate%matrix%GCM_warm%ins_ti1)
    deallocate(climate%matrix%GCM_warm%ins_nlat)
    deallocate(climate%matrix%GCM_warm%ins_nlon)
    deallocate(climate%matrix%GCM_warm%ins_lat)
    deallocate(climate%matrix%GCM_warm%ins_Q_TOA0)
    deallocate(climate%matrix%GCM_warm%ins_Q_TOA1)
    deallocate(climate%matrix%GCM_warm%Q_TOA)
    deallocate(climate%matrix%GCM_warm%Albedo)
    deallocate(climate%matrix%GCM_warm%I_abs)

    deallocate(climate%matrix%GCM_cold%ins_t0)
    deallocate(climate%matrix%GCM_cold%ins_t1)
    deallocate(climate%matrix%GCM_cold%ins_ti0)
    deallocate(climate%matrix%GCM_cold%ins_ti1)
    deallocate(climate%matrix%GCM_cold%ins_nlat)
    deallocate(climate%matrix%GCM_cold%ins_nlon)
    deallocate(climate%matrix%GCM_cold%ins_lat)
    deallocate(climate%matrix%GCM_cold%ins_Q_TOA0)
    deallocate(climate%matrix%GCM_cold%ins_Q_TOA1)
    deallocate(climate%matrix%GCM_cold%Q_TOA)
    deallocate(climate%matrix%GCM_cold%Albedo)
    deallocate(climate%matrix%GCM_cold%I_abs)

    ! Get reference absorbed insolation for the GCM snapshots
    CALL initialise_matrix_calc_absorbed_insolation( mesh_new, climate%matrix%GCM_warm, region_name, forcing, ice)
    CALL initialise_matrix_calc_absorbed_insolation( mesh_new, climate%matrix%GCM_cold, region_name, forcing, ice)

  END SUBROUTINE remap_climate_matrix_model

end module climate_matrix