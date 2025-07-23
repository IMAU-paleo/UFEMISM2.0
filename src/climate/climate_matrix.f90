module climate_matrix

! ===== Preamble =====
! ====================

  use precisions                                             , only: dp
  use mpi_basic                                              , only: par, sync
  use control_resources_and_error_messaging                  , only: crash, init_routine, finalise_routine, colour_string
  use model_configuration                                    , only: C
  use parameters
  use mpi_f08, only: MPI_ALLREDUCE, MPI_IN_PLACE, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, MPI_INTEGER
  use mesh_types                                             , only: type_mesh
  use ice_model_types                                        , only: type_ice_model
  use grid_types                                             , only: type_grid
  use climate_model_types                                    , only: type_climate_model, type_climate_model_matrix, type_climate_model_snapshot
  use global_forcing_types                                   , only: type_global_forcing
  use SMB_model_types, only: type_SMB_model
  use climate_realistic                                      , only: initialise_climate_model_realistic, initialise_insolation_forcing
  use reallocate_mod                                         , only: reallocate_bounds
  use netcdf_io_main
  use mesh_data_smoothing, only: smooth_Gaussian
  use SMB_IMAU_ITM, only: run_SMB_model_IMAUITM, initialise_SMB_model_IMAUITM
  use climate_matrix_utilities, only: allocate_climate_snapshot, read_climate_snapshot, adapt_precip_CC, adapt_precip_Roe, get_insolation_at_time
  use assertions_basic, only: assert

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
  subroutine run_climate_model_matrix( mesh, grid, ice, SMB, climate, region_name, time, forcing)
    ! Use CO2 (either prescribed or inversely modelled) to force the 2-snapshot (PI-LGM) climate matrix (Berends et al., 2018)

    ! In/output variables:
    type(type_mesh),                     intent(in)    :: mesh
    type(type_grid),                     intent(in)    :: grid
    type(type_ice_model),                intent(in)    :: ice
    type(type_SMB_model),                intent(in)    :: SMB
    type(type_climate_model),            intent(inout) :: climate
    character(LEN=3),                    intent(in)    :: region_name
    real(dp),                            intent(in)    :: time
    type(type_global_forcing),           intent(in)    :: forcing

    ! Local variables:
    character(LEN=256), parameter                      :: routine_name = 'run_climate_model_matrix'
    integer                                            :: vi ,m

    ! Add routine to path
    call init_routine( routine_name)

    !IF (par%primary)  WRITE(*,"(A)") '      Running climate matrix model...'

    ! Update insolation forcing at model time
    call get_insolation_at_time( mesh, time, climate%snapshot)

    ! moved to global_forcings_main
    !CALL update_CO2_at_model_time( time, forcing)

    ! Use the (CO2 + absorbed insolation)-based interpolation scheme for temperature
    call run_climate_model_matrix_temperature( mesh, grid, ice, SMB, climate, region_name, forcing)

    ! Use the (CO2 + ice-sheet geometry)-based interpolation scheme for precipitation
    call run_climate_model_matrix_precipitation( mesh, grid, ice, climate, region_name, forcing)

    ! == Safety checks from UFE1.x
    ! ================
#if (DO_ASSERTIONS)
    do vi = mesh%vi1, mesh%vi2
    do m = 1, 12
      call assert( climate%T2m( vi,m) >= 150._dp, 'excessively low temperatures (<150K) detected!')
      call assert( climate%T2m( vi,m) >= 0._dp, 'negative temperatures (<0K) detected!')
      call assert( climate%T2m( vi,m) == climate%T2m( vi,m), 'NaN temperatures  detected!')
      call assert( climate%Precip( vi,m) > 0._dp, 'zero/negative precipitation detected!')
      call assert( climate%Precip( vi,m) == climate%Precip( vi,m), 'NaN precipitation detected!')
    end do
    end do
    call sync
#endif

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_climate_model_matrix
  subroutine run_climate_model_matrix_temperature( mesh, grid, ice, SMB, climate, region_name, forcing)
    ! The (CO2 + absorbed insolation)-based matrix interpolation for temperature, from Berends et al. (2018)

    ! In/output variables:
    type(type_mesh),                     intent(in)    :: mesh
    type(type_grid),                     intent(in)    :: grid
    type(type_ice_model),                intent(in)    :: ice
    type(type_SMB_model),                intent(in)    :: SMB
    type(type_climate_model),            intent(inout) :: climate
    character(LEN=3),                    intent(in)    :: region_name
    type(type_global_forcing),           intent(in)    :: forcing

    ! Local variables:
    character(LEN=256), parameter                      :: routine_name = 'run_climate_model_matrix_temperature'
    integer                                            :: vi ,m
    real(dp)                                           :: CO2, w_CO2
    real(dp), DIMENSION(:    ), allocatable            :: w_tot
    real(dp), DIMENSION(:,:  ), allocatable            :: T_ref_GCM
    real(dp), DIMENSION(:    ), allocatable            :: Hs_GCM, lambda_GCM
    real(dp), parameter                                :: w_cutoff = 0.5_dp        ! Crop weights to [-w_cutoff, 1 + w_cutoff]

    ! Add routine to path
    call init_routine( routine_name)

    ! Allocate shared memory
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
    w_CO2 = MAX( -w_cutoff, MIN( 1._dp + w_cutoff, (CO2 - C%climate_matrix_low_CO2_level) / &
                               (C%climate_matrix_high_CO2_level - C%climate_matrix_low_CO2_level) ))

    ! Find the interpolation weights based on absorbed insolation
    ! ===========================================================

    ! Calculate modelled absorbed insolation

    climate%matrix%I_abs( mesh%vi1:mesh%vi2) = 0._dp
    call sync
    do vi = mesh%vi1, mesh%vi2
    do m = 1, 12
      ! Calculate modelled absorbed insolation. Berends et al., 2018 - Eq. 2
      climate%matrix%I_abs( vi) = climate%matrix%I_abs( vi) + & 
                                  climate%snapshot%Q_TOA( vi,m) * (1._dp - SMB%IMAUITM%Albedo( vi, m))  
    end do
    end do
    call sync

    call weighting_fields_matrix_temperature( climate, mesh, grid, region_name, w_CO2, w_tot)

!! ==============================================================================================================
!! In UFE1.x here are two more options, glacial matrix and glacial index
!! lines 1050 - 1080 in climate_module.f90
!! ==============================================================================================================

    ! Interpolate between the GCM snapshots
    ! =====================================

    do vi = mesh%vi1, mesh%vi2

      ! Find matrix-interpolated orography, lapse rate, and temperature
      Hs_GCM( vi     ) = (w_tot( vi) * climate%matrix%GCM_warm%Hs( vi    )) + &
                         ((1._dp - w_tot( vi)) * climate%matrix%GCM_cold%Hs( vi    ))  ! Berends et al., 2018 - Eq. 8
      lambda_GCM( vi ) = (w_tot( vi) * climate%matrix%GCM_warm%lambda( vi)) + &
                         ((1._dp - w_tot( vi)) * climate%matrix%GCM_cold%lambda( vi))  ! Not listed in the article, shame on me!
      T_ref_GCM( vi,:) = (w_tot( vi) * climate%matrix%GCM_warm%T2m( vi,: )) + &
                         ((1._dp - w_tot( vi)) * climate%matrix%GCM_cold%T2m( vi,: ))  ! Berends et al., 2018 - Eq. 6

      ! Adapt temperature to model orography using matrix-derived lapse-rate
      do m = 1, 12
        climate%T2m( vi,:) = T_ref_GCM( vi, m) - lambda_GCM( vi) * (ice%Hs( vi) - Hs_GCM( vi))  ! Berends et al., 2018 - Eq. 11
      end do

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_climate_model_matrix_temperature
  subroutine weighting_fields_matrix_temperature( climate, mesh, grid, region_name, w_CO2, w_tot)
    ! In/output variables:
    type(type_mesh),                          intent(in)    :: mesh
    type(type_grid),                          intent(in)    :: grid
    type(type_climate_model),                 intent(in)    :: climate
    character(LEN=3),                         intent(in)    :: region_name
    real(dp),                                 intent(in)    :: w_CO2
    real(dp), dimension(mesh%vi1:mesh%vi2),   intent(out)   :: w_tot

    ! Local variables:
    character(LEN=256), parameter                      :: routine_name = 'weighting_fields_matrix_temperature'
    integer                                            :: vi ,m
    real(dp), DIMENSION(:    ), allocatable            :: w_ins,  w_ins_smooth,  w_ice
    real(dp)                                           :: w_ins_av, w_ins_denominator
    real(dp), parameter                                :: w_cutoff = 0.5_dp        ! Crop weights to [-w_cutoff, 1 + w_cutoff]

    ! Add routine to path
    call init_routine( routine_name)

    ! Allocate shared memory
    allocate( w_ins(        mesh%vi1:mesh%vi2))
    allocate( w_ins_smooth( mesh%vi1:mesh%vi2))
    allocate( w_ice(        mesh%vi1:mesh%vi2))

    ! Calculate "direct" weighting field
    ! Berends et al., 2018 - Eq. 3

    do vi = mesh%vi1, mesh%vi2
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

    end do
    call sync
    
    w_ins_av      = MAX( -w_cutoff, MIN( 1._dp + w_cutoff, (SUM( climate%matrix%I_abs         )      - SUM( climate%matrix%GCM_cold%I_abs)     ) / &
                                                           (SUM( climate%matrix%GCM_warm%I_abs)      - SUM( climate%matrix%GCM_cold%I_abs)     ) ))
    ! Smooth the weighting field
    w_ins_smooth( mesh%vi1:mesh%vi2) = w_ins( mesh%vi1:mesh%vi2)

    call smooth_Gaussian( mesh, grid, C%output_dir, w_ins_smooth, 200000._dp)

    ! Combine unsmoothed, smoothed, and regional average weighting fields (Berends et al., 2018, Eq. 4)
    if (region_name == 'NAM' .OR. region_name == 'EAS') then
      w_ice( mesh%vi1:mesh%vi2) = (1._dp * w_ins(        mesh%vi1:mesh%vi2) + &
                                   3._dp * w_ins_smooth( mesh%vi1:mesh%vi2) + &
                                   3._dp * w_ins_av) / 7._dp
    elseif (region_name == 'GRL' .OR. region_name == 'ANT') then
      ! Use only the regional (averaged) and smoothed weights
      w_ice( mesh%vi1:mesh%vi2) = (1._dp * w_ins_smooth( mesh%vi1:mesh%vi2) + &
                                   6._dp * w_ins_av) / 7._dp
    end if

    ! Combine interpolation weights from absorbed insolation and CO2 into the final weights fields
    ! Berends et al., 2018 - Eqs. 5, 9 with weights 0.5 for NAM & EAS, and 0.75 for ANT
    ! Generalised: "switch" between matrix method and glacial index method by altering C%climate_matrix_CO2vsice_<region>
    if         (region_name == 'NAM') then
      w_tot( mesh%vi1:mesh%vi2) = (C%climate_matrix_CO2vsice_NAM * w_CO2) + ((1._dp - C%climate_matrix_CO2vsice_NAM) * w_ice( mesh%vi1:mesh%vi2))
    elseif     (region_name == 'EAS') then
      w_tot( mesh%vi1:mesh%vi2) = (C%climate_matrix_CO2vsice_EAS * w_CO2) + ((1._dp - C%climate_matrix_CO2vsice_EAS) * w_ice( mesh%vi1:mesh%vi2))
    elseif     (region_name == 'GRL') then
      w_tot( mesh%vi1:mesh%vi2) = (C%climate_matrix_CO2vsice_GRL * w_CO2) + ((1._dp - C%climate_matrix_CO2vsice_GRL) * w_ice( mesh%vi1:mesh%vi2))
    elseif     (region_name == 'ANT') then
      w_tot( mesh%vi1:mesh%vi2) = (C%climate_matrix_CO2vsice_ANT * w_CO2) + ((1._dp - C%climate_matrix_CO2vsice_ANT) * w_ice( mesh%vi1:mesh%vi2))
    end if

    ! If a glacial index is used, weight will depend only on CO2
    !IF (C%climate_matrix_switch_glacial_index) THEN
    !  w_tot = w_CO2
    !END IF

    ! Finalise routine path
    call finalise_routine( routine_name)
    
  end subroutine weighting_fields_matrix_temperature
  subroutine run_climate_model_matrix_precipitation( mesh, grid, ice, climate, region_name, forcing)
    ! The (CO2 + ice geometry)-based matrix interpolation for precipitation, from Berends et al. (2018)
    ! For NAM and EAS, this is based on local ice geometry and uses the Roe&Lindzen precipitation model for downscaling.
    ! For GRL and ANT, this is based on total ice volume,  and uses the simple CC   precipitation model for downscaling.
    ! The rationale for this difference is that glacial-interglacial differences in ice geometry are much more
    ! dramatic in NAM and EAS than they are in GRL and ANT.

    ! In/output variables:
    type(type_mesh),                     intent(in)    :: mesh
    type(type_grid),                     intent(in)    :: grid
    type(type_ice_model),                intent(in)    :: ice
    type(type_climate_model),            intent(inout) :: climate
    type(type_global_forcing),           intent(in)    :: forcing
    character(LEN=3),                    intent(in)    :: region_name

    ! Local variables:
    character(LEN=256), parameter                      :: routine_name = 'run_climate_model_matrix_precipitation'
    integer                                            :: vi, m
    real(dp), DIMENSION(:), allocatable                ::  w_warm,  w_cold
    real(dp)                                           :: w_tot
    real(dp), DIMENSION(:,:), allocatable              :: T_ref_GCM, P_ref_GCM
    real(dp), DIMENSION(:),   allocatable              :: Hs_GCM

    real(dp), parameter                                :: w_cutoff = 0.25_dp        ! Crop weights to [-w_cutoff, 1 + w_cutoff]

    ! Add routine to path
    call init_routine( routine_name)

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

    call weighting_fields_matrix_precipitation( climate, mesh, grid, ice, region_name, forcing, w_tot, w_warm, w_cold)

    ! Interpolate the GCM snapshots
    ! =============================

    do vi = mesh%vi1, mesh%vi2

      T_ref_GCM( vi,:) =      (w_warm( vi) *     climate%matrix%GCM_warm%T2m(    vi,: )) + & 
                              (w_cold( vi) *     climate%matrix%GCM_cold%T2m(    vi,: ))   ! Berends et al., 2018 - Eq. 6
                              
      P_ref_GCM( vi,:) = EXP( (w_warm( vi) * LOG(climate%matrix%GCM_warm%Precip( vi,:))) + &
                              (w_cold( vi) * LOG(climate%matrix%GCM_cold%Precip( vi,:)))) ! Berends et al., 2018 - Eq. 7
                              
      Hs_GCM(    vi  ) =      (w_warm( vi) *     climate%matrix%GCM_warm%Hs(     vi   )) + &
                              (w_cold( vi) *     climate%matrix%GCM_cold%Hs(     vi   ))   ! Berends et al., 2018 - Eq. 8                  
    end do

    ! Downscale precipitation from the coarse-resolution reference
    ! GCM orography to the fine-resolution ice-model orography
    ! ========================================================

    if (region_name == 'NAM' .OR. region_name == 'EAS') then
      ! Use the Roe&Lindzen precipitation model to do this; Berends et al., 2018, Eqs. A3-A7
      call adapt_precip_Roe( mesh, Hs_GCM,   T_ref_GCM  , climate%matrix%PD_obs%Wind_LR, climate%matrix%PD_obs%Wind_DU, P_ref_GCM, &
                                   ice%Hs, climate%T2m, climate%matrix%PD_obs%Wind_LR, climate%matrix%PD_obs%Wind_DU, climate%Precip)
    elseif (region_name == 'GRL' .OR. region_name == 'ANT') then
      ! Use a simpler temperature-based correction; Berends et al., 2018, Eq. 14
      call adapt_precip_CC( mesh, ice%Hs, Hs_GCM, T_ref_GCM, P_ref_GCM, climate%Precip, region_name)
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_climate_model_matrix_precipitation

  subroutine weighting_fields_matrix_precipitation( climate, mesh, grid, ice, region_name, forcing, w_tot, w_warm, w_cold)
    ! In/output variables
    type(type_mesh),                         intent(in)    :: mesh
    type(type_grid),                         intent(in)    :: grid
    type(type_ice_model),                    intent(in)    :: ice
    type(type_climate_model),                intent(inout) :: climate
    character(LEN=3),                        intent(in)    :: region_name
    type(type_global_forcing),               intent(in)    :: forcing
    real(dp),                                intent(inout) :: w_tot
    real(dp), dimension(mesh%vi1:mesh%vi2),  intent(out)   :: w_warm, w_cold

    ! Local variables:
    character(LEN=256), parameter                          :: routine_name = 'weighting_fields_matrix_precipitation'
    integer                                                :: vi, m
    real(dp), parameter                                    :: w_cutoff = 0.25_dp        ! Crop weights to [-w_cutoff, 1 + w_cutoff]

    ! Add routine to path
    call init_routine( routine_name)

    select case (region_name)

    case ('NAM' , 'EAS')
      ! Combine total + local ice thicness; Berends et al., 2018, Eq. 12

      ! Then the local ice thickness term
      do vi = mesh%vi1, mesh%vi2

        if (climate%matrix%GCM_warm%Hs( vi) < climate%matrix%GCM_PI%Hs( vi) + 50._dp) then
          if (climate%matrix%GCM_cold%Hs( vi) < climate%matrix%GCM_PI%Hs( vi) + 50._dp) then
            ! No ice in any GCM state. Use only total ice volume.
            w_cold( vi) = MAX(-w_cutoff, MIN(1._dp + w_cutoff, w_tot ))
            w_warm( vi) = 1._dp - w_cold( vi)
          else
            ! No ice in warm climate, ice in cold climate. Linear inter- / extrapolation.
            w_cold( vi) = MAX(-w_cutoff, MIN(1._dp + w_cutoff, ((ice%Hs( vi) - climate%matrix%GCM_PI%Hs( vi)) / (climate%matrix%GCM_cold%Hs( vi) - climate%matrix%GCM_PI%Hs( vi))) * w_tot ))
            w_warm( vi)  = 1._dp - w_cold( vi)
          end if
        else
          ! Ice in both GCM states.  Linear inter- / extrapolation
          w_cold( vi) = MAX(-w_cutoff, MIN(1._dp + w_cutoff, ((ice%Hs( vi) - climate%matrix%GCM_PI%Hs( vi)) / (climate%matrix%GCM_cold%Hs( vi) - climate%matrix%GCM_PI%Hs( vi))) * w_tot ))
          w_warm( vi)  = 1._dp - w_cold( vi)
        end if

      end do

      w_cold( mesh%vi1:mesh%vi2) = w_cold( mesh%vi1:mesh%vi2) * w_tot 
      
      ! Smooth the weighting field
      call smooth_Gaussian( mesh, grid, C%output_dir, w_cold, 200000._dp)

      w_warm( mesh%vi1:mesh%vi2) = 1._dp - w_cold( mesh%vi1:mesh%vi2)

    case ('GRL' , 'ANT')
      ! Use only total ice volume and CO2; Berends et al., 2018, Eq. 13

      w_cold( mesh%vi1:mesh%vi2) = w_tot 
      w_warm( mesh%vi1:mesh%vi2) = 1._dp - w_cold( mesh%vi1:mesh%vi2)

    end select

    if (C%climate_matrix_switch_glacial_index_precip) then ! If a glacial index is used for the precipitation forcing, it will only depend on CO2
      w_tot = 1._dp - (MAX( -w_cutoff, MIN( 1._dp + w_cutoff, (forcing%CO2_obs - C%climate_matrix_low_CO2_level) & 
              / (C%climate_matrix_high_CO2_level - C%climate_matrix_low_CO2_level) )) )
      w_cold( mesh%vi1:mesh%vi2) = w_tot
      w_warm( mesh%vi1:mesh%vi2) = 1._dp - w_cold( mesh%vi1:mesh%vi2)
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine weighting_fields_matrix_precipitation

  subroutine initialise_climate_matrix( mesh, grid, ice, climate, region_name, forcing)

    ! In/output variables:
    type(type_mesh),                     intent(in)    :: mesh
    type(type_grid),                     intent(in)    :: grid !used to smooth later on, check if grid is called during initialise
    type(type_ice_model),                intent(in)    :: ice
    type(type_climate_model),            intent(inout) :: climate
    character(LEN=3),                    intent(in)    :: region_name
    type(type_global_forcing),           intent(in)    :: forcing    

    ! Local variables:
    character(LEN=256), parameter                      :: routine_name = 'initialise_climate_matrix'
    integer                                            :: vi, m

    ! Add routine to path
    call init_routine( routine_name)

    ! Allocate shared memory
    allocate( climate%matrix%I_abs(           mesh%vi1:mesh%vi2))
    allocate( climate%matrix%GCM_bias_T2m(    mesh%vi1:mesh%vi2, 12))
    allocate( climate%matrix%GCM_bias_Precip( mesh%vi1:mesh%vi2, 12))

    ! Allocate memory for the regional ERA40 climate and the final applied climate
    call allocate_climate_snapshot( mesh, climate%matrix%PD_obs,   name = 'PD_obs'  )
    call allocate_climate_snapshot( mesh, climate%matrix%GCM_PI,   name = 'GCM_PI'  )
    call allocate_climate_snapshot( mesh, climate%matrix%GCM_warm, name = 'GCM_warm')
    call allocate_climate_snapshot( mesh, climate%matrix%GCM_cold, name = 'GCM_cold')

    call read_climate_snapshot( C%climate_matrix_filename_PD_obs_climate       , mesh, climate%matrix%PD_obs  )
    call read_climate_snapshot( C%climate_matrix_filename_climate_snapshot_PI  , mesh, climate%matrix%GCM_PI  )
    call read_climate_snapshot( C%climate_matrix_filename_climate_snapshot_warm, mesh, climate%matrix%GCM_warm)
    call read_climate_snapshot( C%climate_matrix_filename_climate_snapshot_cold, mesh, climate%matrix%GCM_cold)
    
    ! Get the orbit time
    climate%matrix%GCM_PI%orbit_time   = 0._dp
    climate%matrix%GCM_warm%orbit_time = C%climate_matrix_warm_orbit_time
    climate%matrix%GCM_cold%orbit_time = C%climate_matrix_cold_orbit_time

    ! Calculate spatially variable lapse rate

    ! Use a uniform value for the warm snapshot [this assumes "warm" is actually identical to PI!]
    climate%matrix%GCM_warm%lambda( mesh%vi1:mesh%vi2) = C%climate_matrix_constant_lapserate

    if     (region_name == 'NAM' .OR. region_name == 'EAS') then
      call initialise_matrix_calc_spatially_variable_lapserate( mesh, grid, climate%matrix%GCM_PI, climate%matrix%GCM_cold)
    elseif (region_name == 'GLR' .OR. region_name == 'ANT') then
      climate%matrix%GCM_cold%lambda( mesh%vi1:mesh%vi2) = C%climate_matrix_constant_lapserate
      call sync
    end if

    ! Calculate GCM bias
    call initialise_matrix_calc_GCM_bias( mesh, climate%matrix%GCM_PI, climate%matrix%PD_obs, &
      climate%matrix%GCM_bias_T2m, climate%matrix%GCM_bias_Precip)

    ! Apply bias correction
    if (C%climate_matrix_biascorrect_warm) call initialise_matrix_apply_bias_correction( mesh, climate%matrix%GCM_warm, &
      climate%matrix%GCM_bias_T2m, climate%matrix%GCM_bias_Precip)
    if (C%climate_matrix_biascorrect_warm) call initialise_matrix_apply_bias_correction( mesh, climate%matrix%GCM_cold, &
      climate%matrix%GCM_bias_T2m, climate%matrix%GCM_bias_Precip)

    ! Get reference absorbed insolation for the GCM snapshots
    call initialise_matrix_calc_absorbed_insolation( mesh, climate%matrix%GCM_warm, region_name, forcing, ice)
    call initialise_matrix_calc_absorbed_insolation( mesh, climate%matrix%GCM_cold, region_name, forcing, ice)

    ! initialise the insolation forcing
    call initialise_insolation_forcing( climate%snapshot, mesh) ! this will initialise climate%snapshot%Q_TOA

    ! Initialise applied climate with present-day observations

    do vi = mesh%vi1, mesh%vi2
    do m = 1, 12
      climate%T2m(     vi,m) = climate%matrix%PD_obs%T2m(     vi,m)
      climate%Precip(  vi,m) = climate%matrix%PD_obs%Precip(  vi,m)
      climate%Wind_LR( vi,m) = climate%matrix%PD_obs%Wind_LR( vi,m)
      climate%Wind_DU( vi,m) = climate%matrix%PD_obs%Wind_DU( vi,m)
    end do
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_climate_matrix
  
  subroutine initialise_matrix_calc_GCM_bias( mesh, GCM_PI, PD_obs, GCM_bias_T2m, GCM_bias_Precip)
    ! Calculate the GCM bias in temperature and precipitation
    !
    ! Account for the fact that the GCM PI snapshot has a lower resolution, and therefore
    ! a different surface elevation than the PD observed climatology!

    ! In/output variables:
    type(type_mesh),                            intent(in)    :: mesh
    type(type_climate_model_snapshot),          intent(in)    :: GCM_PI, PD_obs 
    real(dp), dimension(mesh%vi1:mesh%vi2, 12), intent(out)   :: GCM_bias_T2m
    real(dp), dimension(mesh%vi1:mesh%vi2, 12), intent(out)   :: GCM_bias_Precip

    ! Local variables:
    character(LEN=256), parameter                             :: routine_name = 'initialise_matrix_calc_GCM_bias'
    integer                                                   :: vi,m
    real(dp)                                                  :: T2m_SL_GCM, T2m_SL_obs

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate bias
    do vi = mesh%vi1, mesh%vi2
    do m = 1, 12

      ! Scale modelled and observed temperature to sea level using a constant lapse rate
      T2m_SL_GCM = GCM_PI%T2m( vi,m) + GCM_PI%Hs( vi) * C%climate_matrix_constant_lapserate
      T2m_SL_obs = PD_obs%T2m( vi,m) + PD_obs%Hs( vi) * C%climate_matrix_constant_lapserate

      ! Calculate bias
      GCM_bias_T2m(    vi,m) = T2m_SL_GCM           - T2m_SL_obs
      GCM_bias_Precip( vi,m) = GCM_PI%Precip( vi,m) / PD_obs%Precip( vi,m)

    end do
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_matrix_calc_GCM_bias
  subroutine initialise_matrix_apply_bias_correction( mesh, snapshot, bias_T2m, bias_Precip)
    ! Apply a bias correction to this (GCM) snapshot

    ! In/output variables:
    type(type_mesh),                            intent(in)    :: mesh
    type(type_climate_model_snapshot),          intent(inout) :: snapshot
    real(dp), dimension(mesh%vi1:mesh%vi2, 12), intent(in)    :: bias_T2m
    real(dp), dimension(mesh%vi1:mesh%vi2, 12), intent(in)    :: bias_Precip

    ! Local variables:
    character(LEN=256), parameter                       :: routine_name = 'initialise_matrix_apply_bias_correction'
    integer                                             :: vi,m

    ! Add routine to path
    call init_routine( routine_name)

    ! Apply bias correction
    do vi = mesh%vi1, mesh%vi2
    do m = 1, 12
      snapshot%T2m(    vi,m) = snapshot%T2m(    vi,m) - bias_T2m(    vi,m)
      snapshot%Precip( vi,m) = snapshot%Precip( vi,m) / bias_Precip( vi,m)
    end do
    end do
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_matrix_apply_bias_correction
  subroutine initialise_matrix_calc_spatially_variable_lapserate( mesh, grid_smooth, snapshot_PI, snapshot)
    ! Calculate the spatially variable lapse-rate (for non-PI GCM climates; see Berends et al., 2018)
    ! Only meaningful for climates where there is ice (LGM, M2_Medium, M2_Large),
    ! and only intended for North America and Eurasia

    ! In/output variables:
    type(type_mesh),                      intent(in)    :: mesh
    type(type_grid),                      intent(in)    :: grid_smooth
    type(type_climate_model_snapshot),    intent(in)    :: snapshot_PI
    type(type_climate_model_snapshot),    intent(inout) :: snapshot

    ! Local variables:
    character(LEN=256), parameter                       :: routine_name = 'initialise_matrix_calc_spatially_variable_lapserate'
    integer                                             :: vi,m
    integer,  DIMENSION(:    ), allocatable             :: mask_calc_lambda
    real(dp)                                            :: dT_mean_nonice
    real(dp)                                            :: lambda_mean_ice

    ! Add routine to path
    call init_routine( routine_name)

    ! Allocate shared memory
    allocate( mask_calc_lambda( mesh%vi1:mesh%vi2))

    ! Determine where the variable lapse rate should be calculated
    ! (i.e. where has the surface elevation increased substantially)
    ! ==============================================================

    do vi = mesh%vi1, mesh%vi2

      if (snapshot%Hs( vi) > snapshot_PI%Hs( vi) + 100._dp) then
        mask_calc_lambda( vi) = 1
      else
        mask_calc_lambda( vi) = 0
      end if

    end do

    ! Calculate the regional average temperature change outside of the ice sheet
    call calc_regional_avrg_temperature_outside_icesheet(mesh, snapshot, snapshot_PI, mask_calc_lambda, dT_mean_nonice)
    ! Calculate the lapse rate over the ice itself
    call calc_lapse_rate_over_ice(mesh, snapshot, snapshot_PI, dT_mean_nonice, mask_calc_lambda, lambda_mean_ice)

    ! Apply mean lapse-rate over ice to the rest of the region
    ! ========================================================

    do vi = mesh%vi1, mesh%vi2
      if (mask_calc_lambda( vi) == 0) snapshot%lambda( vi) = lambda_mean_ice
    end do

    ! Smooth the lapse rate field with a 160 km Gaussian filter
    call smooth_Gaussian( mesh, grid_smooth, C%output_dir, snapshot%lambda, 160000._dp)

    ! Normalise the entire region to a mean lapse rate of 8 K /km
    snapshot%lambda( mesh%vi1:mesh%vi2) = snapshot%lambda( mesh%vi1:mesh%vi2) * (C%climate_matrix_constant_lapserate / lambda_mean_ice)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_matrix_calc_spatially_variable_lapserate

  subroutine calc_regional_avrg_temperature_outside_icesheet(mesh, snapshot, snapshot_PI, mask_calc_lambda, dT_mean_nonice)
    ! In/output variables:
    type(type_mesh),                       intent(in)    :: mesh
    type(type_climate_model_snapshot),     intent(in)    :: snapshot_PI
    type(type_climate_model_snapshot),     intent(in)    :: snapshot
    integer, dimension(mesh%vi1:mesh%vi2), intent(in)    :: mask_calc_lambda
    real(dp),                              intent(out)   :: dT_mean_nonice

    ! Local variables:
    character(LEN=256), parameter                        :: routine_name = 'calc_regional_avrg_temperature_outside_icesheet'
    integer                                              :: vi,m
    integer                                              :: n_nonice
    integer                                              :: ierr

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate the regional average temperature change outside of the ice sheet
    dT_mean_nonice = 0._dp
    n_nonice       = 0
    do vi = mesh%vi1, mesh%vi2
    do m = 1, 12
      if (mask_calc_lambda( vi) == 0) then
        dT_mean_nonice = dT_mean_nonice + snapshot%T2m( vi,m) - snapshot_PI%T2m( vi,m)
        n_nonice = n_nonice + 1
      end if
    end do
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, dT_mean_nonice, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, n_nonice,       1, MPI_INTEGER,          MPI_SUM, MPI_COMM_WORLD, ierr)

    dT_mean_nonice = dT_mean_nonice / real(n_nonice,dp)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_regional_avrg_temperature_outside_icesheet
  subroutine calc_lapse_rate_over_ice(mesh, snapshot, snapshot_PI, dT_mean_nonice, mask_calc_lambda, lambda_mean_ice)
    ! In/output variables:
    type(type_mesh),                       intent(in)    :: mesh
    type(type_climate_model_snapshot),     intent(in)    :: snapshot_PI
    type(type_climate_model_snapshot),     intent(inout) :: snapshot
    real(dp),                              intent(in)    :: dT_mean_nonice
    integer, dimension(mesh%vi1:mesh%vi2), intent(in)    :: mask_calc_lambda
    real(dp),                              intent(out)   :: lambda_mean_ice

    ! Local variables:
    character(LEN=256), parameter                        :: routine_name = 'calc_lapse_rate_over_ice'
    integer                                              :: vi,m
    integer                                              :: n_ice
    integer                                              :: ierr
    real(dp), parameter                                  :: lambda_min = 0.002_dp
    real(dp), parameter                                  :: lambda_max = 0.05_dp

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate the lapse rate over the ice itself

    lambda_mean_ice = 0._dp
    n_ice           = 0

    do vi = mesh%vi1, mesh%vi2

      if (mask_calc_lambda( vi) == 1) then

        do m = 1, 12
          ! Berends et al., 2018 - Eq. 10
          snapshot%lambda( vi) = snapshot%lambda( vi) + 1/12._dp * MAX(lambda_min, MIN(lambda_max, &
            -(snapshot%T2m( vi,m) - (snapshot_PI%T2m( vi,m) + dT_mean_nonice)) / (snapshot%Hs( vi) - snapshot_PI%Hs( vi))))
        end do

        lambda_mean_ice = lambda_mean_ice + snapshot%lambda( vi)
        n_ice = n_ice + 1

      end if

    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, lambda_mean_ice, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, n_ice,           1, MPI_INTEGER,          MPI_SUM, MPI_COMM_WORLD, ierr)

    lambda_mean_ice = lambda_mean_ice / n_ice

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_lapse_rate_over_ice

  subroutine initialise_matrix_calc_absorbed_insolation( mesh, snapshot, region_name, forcing, ice)
    ! Calculate the yearly absorbed insolation for this (regional) GCM snapshot, to be used in the matrix interpolation

    ! In/output variables:
    type(type_mesh),                      intent(in)    :: mesh
    type(type_climate_model_snapshot),    intent(inout) :: snapshot
    character(LEN=3),                     intent(in)    :: region_name
    type(type_global_forcing),            intent(in)    :: forcing
    type(type_ice_model),                 intent(in)    :: ice

    ! Local variables:
    character(LEN=256), parameter                       :: routine_name = 'initialise_matrix_calc_absorbed_insolation'
    integer                                             :: vi,m,i
    type(type_ice_model)                                :: ice_dummy
    type(type_climate_model)                            :: climate_dummy
    type(type_SMB_model)                                :: SMB_dummy
    character(LEN=256)                                  :: choice_SMB_IMAUITM_init_firn_dummy

    ! Add routine to path
    call init_routine( routine_name)
    ! Initialise the insolation variables inside snapshot
    call initialise_insolation_forcing(snapshot, mesh)
    ! Get insolation at the desired time from the insolation NetCDF file
    ! ==================================================================

    call get_insolation_at_time( mesh, snapshot%orbit_time, snapshot)

    ! Create temporary "dummy" climate, ice & SMB data structures,
    ! so we can run the SMB model and determine the reference albedo field
    ! ====================================================================

    ! Climate
    ! =======

    ! Allocate shared memory
    allocate( climate_dummy%T2m(    mesh%vi1:mesh%vi2, 12))
    allocate( climate_dummy%Precip( mesh%vi1:mesh%vi2, 12))
    allocate( climate_dummy%snapshot%Q_TOA(  mesh%vi1:mesh%vi2, 12))

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
    do vi = mesh%vi1, mesh%vi2
    
   ! In IMAU-ICE SMB it uses region%mask_noice in UFE2 is ice%mask_noice, I will keep the masks from above for ice_dummy
   ! and make ice_dummy%mask_noice = ice%mask_noice to run the SMB using the dummy, following IMAU-ICE code..
      ice_dummy%mask_noice( vi) = ice%mask_noice( vi) 

      if (snapshot%Hs( vi) == MINVAL(snapshot%Hs)) then
        ice_dummy%mask_icefree_ocean( vi) = .true.
      else
        ice_dummy%mask_icefree_ocean( vi) = .false.
      end if

      ! this IF is like (climate%Mask_ice( vi) > .3_dp) in Ufe1.x
      if (snapshot%Hs( vi) > 100._dp .AND. SUM(snapshot%T2m( vi,:)) / 12._dp < 0._dp) then
        ice_dummy%mask_grounded_ice(   vi) = .true.
      else
        ice_dummy%mask_grounded_ice(   vi) = .false.
      end if

      ! mask_shelf is used in the SMB model only to find open ocean; since mask_ocean
      ! in this case already marks only open ocean, no need to look for shelves
      ice_dummy%mask_floating_ice( vi) = .false.

    end do

    ! SMB
    ! ===
    call initialise_SMB_model_IMAUITM( mesh, ice, SMB_dummy%IMAUITM, region_name)
    allocate( SMB_dummy%SMB             (mesh%vi1:mesh%vi2))
    SMB_dummy%SMB = 0._dp

    ! Initialisation choice
    if     (region_name == 'NAM') then
      choice_SMB_IMAUITM_init_firn_dummy = C%choice_SMB_IMAUITM_init_firn_NAM
    elseif (region_name == 'EAS') then
      choice_SMB_IMAUITM_init_firn_dummy = C%choice_SMB_IMAUITM_init_firn_EAS
    elseif (region_name == 'GRL') then
      choice_SMB_IMAUITM_init_firn_dummy = C%choice_SMB_IMAUITM_init_firn_GRL
    elseif (region_name == 'ANT') then
      choice_SMB_IMAUITM_init_firn_dummy = C%choice_SMB_IMAUITM_init_firn_ANT
    end if

    if     (choice_SMB_IMAUITM_init_firn_dummy == 'uniform') then
      ! do nothing
    else
      call crash('climate matrix only implemented with uniform init firn"' // TRIM( choice_SMB_IMAUITM_init_firn_dummy) // '"!')
    end if

    ! Run the SMB model for 10 years for this particular climate
    ! (experimentally determined to be long enough to converge)
    do i = 1, 10
      call run_SMB_model_IMAUITM( mesh, ice_dummy, SMB_dummy, climate_dummy)
    end do

    ! Calculate yearly total absorbed insolation
    snapshot%I_abs( mesh%vi1:mesh%vi2) = 0._dp
    do vi = mesh%vi1, mesh%vi2
    do m = 1, 12
      snapshot%I_abs( vi) = snapshot%I_abs( vi) + snapshot%Q_TOA( vi,m) * (1._dp - SMB_dummy%IMAUITM%Albedo( vi,m))
    end do
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_matrix_calc_absorbed_insolation

  subroutine remap_climate_matrix_model( mesh_new, climate, region_name, grid, ice, forcing)

    ! In- and output variables
    !TYPE(type_mesh),                        INTENT(IN)    :: mesh_old
    type(type_mesh),                        intent(in)    :: mesh_new
    type(type_climate_model),               intent(inout) :: climate
    character(LEN=3),                       intent(in)    :: region_name
    type(type_grid),                        intent(in)    :: grid
    type(type_ice_model),                   intent(in)    :: ice
    type(type_global_forcing),              intent(in)    :: forcing
    ! Local variables:
    character(LEN=256), parameter                         :: routine_name = 'remap_climate_matrix_model'

    ! Add routine to path
    call init_routine( routine_name)

    if (par%primary)  WRITE(*,"(A)") '      Remapping climate matrix model data to the new mesh...'
  
    ! reallocate main variables of GCM snapshots
    call remap_climate_matrix_snapshot( mesh_new, climate%matrix%PD_obs)
    call remap_climate_matrix_snapshot( mesh_new, climate%matrix%GCM_PI)
    call remap_climate_matrix_snapshot( mesh_new, climate%matrix%GCM_warm)
    call remap_climate_matrix_snapshot( mesh_new, climate%matrix%GCM_cold)

    ! reallocate main variables of climate%snapshot for insolation
    call reallocate_bounds( climate%snapshot%ins_Q_TOA0, mesh_new%vi1, mesh_new%vi2,12)
    call reallocate_bounds( climate%snapshot%ins_Q_TOA1, mesh_new%vi1, mesh_new%vi2,12)
    call reallocate_bounds( climate%snapshot%Q_TOA, mesh_new%vi1, mesh_new%vi2,12)
  
    call reallocate_bounds(climate%matrix%I_abs, mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds(climate%matrix%GCM_bias_T2m, mesh_new%vi1, mesh_new%vi2, 12)
    call reallocate_bounds(climate%matrix%GCM_bias_Precip, mesh_new%vi1, mesh_new%vi2, 12)
  
    ! read the snapshots for the new mesh
    call read_climate_snapshot( C%climate_matrix_filename_PD_obs_climate       , mesh_new, climate%matrix%PD_obs  )
    call read_climate_snapshot( C%climate_matrix_filename_climate_snapshot_PI  , mesh_new, climate%matrix%GCM_PI  )
    call read_climate_snapshot( C%climate_matrix_filename_climate_snapshot_warm, mesh_new, climate%matrix%GCM_warm)
    call read_climate_snapshot( C%climate_matrix_filename_climate_snapshot_cold, mesh_new, climate%matrix%GCM_cold)

    ! Use a uniform value for the warm snapshot [this assumes "warm" is actually identical to PI!]
    climate%matrix%GCM_warm%lambda = C%climate_matrix_constant_lapserate

    if     (region_name == 'NAM' .OR. region_name == 'EAS') then
      call initialise_matrix_calc_spatially_variable_lapserate( mesh_new, grid, climate%matrix%GCM_PI, climate%matrix%GCM_cold)
    elseif (region_name == 'GLR' .OR. region_name == 'ANT') then
      climate%matrix%GCM_cold%lambda = C%climate_matrix_constant_lapserate
      call sync
    end if

    ! Calculate GCM bias
    call initialise_matrix_calc_GCM_bias( mesh_new, climate%matrix%GCM_PI, climate%matrix%PD_obs, &
      climate%matrix%GCM_bias_T2m, climate%matrix%GCM_bias_Precip)

    ! Apply bias correction
    if (C%climate_matrix_biascorrect_warm) call initialise_matrix_apply_bias_correction( mesh_new, climate%matrix%GCM_warm, &
      climate%matrix%GCM_bias_T2m, climate%matrix%GCM_bias_Precip)
    if (C%climate_matrix_biascorrect_warm) call initialise_matrix_apply_bias_correction( mesh_new, climate%matrix%GCM_cold, &
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
    call initialise_matrix_calc_absorbed_insolation( mesh_new, climate%matrix%GCM_warm, region_name, forcing, ice)
    call initialise_matrix_calc_absorbed_insolation( mesh_new, climate%matrix%GCM_cold, region_name, forcing, ice)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remap_climate_matrix_model

  subroutine remap_climate_matrix_snapshot (mesh_new, snapshot)

    type(type_mesh),                        intent(in)    :: mesh_new
    type(type_climate_model_snapshot),      intent(inout) :: snapshot

    ! local variables
    character(LEN=256), parameter                         :: routine_name = 'remap_climate_matrix_snapshot'

    ! Add routine to path
    call init_routine( routine_name)
  
    call reallocate_bounds(snapshot%Precip, mesh_new%vi1, mesh_new%vi2, 12)
    call reallocate_bounds(snapshot%T2m, mesh_new%vi1, mesh_new%vi2, 12)
    call reallocate_bounds(snapshot%Wind_WE, mesh_new%vi1, mesh_new%vi2, 12)
    call reallocate_bounds(snapshot%Wind_SN, mesh_new%vi1, mesh_new%vi2, 12)
    call reallocate_bounds(snapshot%Wind_LR, mesh_new%vi1, mesh_new%vi2, 12)
    call reallocate_bounds(snapshot%Wind_DU, mesh_new%vi1, mesh_new%vi2, 12)
    call reallocate_bounds(snapshot%Hs, mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds(snapshot%lambda, mesh_new%vi1, mesh_new%vi2)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine


end module climate_matrix