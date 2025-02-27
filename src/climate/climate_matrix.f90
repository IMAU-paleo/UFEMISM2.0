module climate_matrix

  SUBROUTINE initialise_climate_matrix( grid, climate, region_name, mask_noice)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_climate_model),            INTENT(INOUT) :: climate
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)    :: mask_noice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_climate_matrix'
    INTEGER                                            :: i,j,m
    LOGICAL                                            :: found_winds_PD_obs, found_winds_PI, found_winds_warm, found_winds_cold

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, climate%matrix%I_abs          , climate%matrix%wI_abs          )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, climate%matrix%GCM_bias_T2m   , climate%matrix%wGCM_bias_T2m   )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, climate%matrix%GCM_bias_Precip, climate%matrix%wGCM_bias_Precip)

    ! Allocate memory for the regional ERA40 climate and the final applied climate
    CALL allocate_climate_snapshot( grid, climate%matrix%PD_obs,   name = 'PD_obs'  )
    CALL allocate_climate_snapshot( grid, climate%matrix%GCM_PI,   name = 'GCM_PI'  )
    CALL allocate_climate_snapshot( grid, climate%matrix%GCM_warm, name = 'GCM_warm')
    CALL allocate_climate_snapshot( grid, climate%matrix%GCM_cold, name = 'GCM_cold')

    ! Read climate data from files
    CALL read_climate_snapshot( C%filename_PD_obs_climate       , grid, climate%matrix%PD_obs  , found_winds_PD_obs)
    CALL read_climate_snapshot( C%filename_climate_snapshot_PI  , grid, climate%matrix%GCM_PI  , found_winds_PI    )
    CALL read_climate_snapshot( C%filename_climate_snapshot_warm, grid, climate%matrix%GCM_warm, found_winds_warm  )
    CALL read_climate_snapshot( C%filename_climate_snapshot_cold, grid, climate%matrix%GCM_cold, found_winds_cold  )

    ! Safety
    IF (.NOT. found_winds_PD_obs) CALL crash('couldnt find wind fields for PD climate in file ' // TRIM( C%filename_PD_obs_climate))

    ! If no wind fields are provided in the GCM snapshots (as is usually the case), use the present-day observed winds instead
    IF (.NOT. found_winds_PI) THEN
      IF (par%master) CALL warning('no wind fields found for PI   climate snapshot; using PD observed winds instead')
      climate%matrix%GCM_PI%Wind_LR(   :,:,grid%i1:grid%i2) = climate%matrix%PD_obs%wind_LR( :,:,grid%i1:grid%i2)
      climate%matrix%GCM_PI%Wind_DU(   :,:,grid%i1:grid%i2) = climate%matrix%PD_obs%Wind_DU( :,:,grid%i1:grid%i2)
    END IF
    IF (.NOT. found_winds_warm) THEN
      IF (par%master) CALL warning('no wind fields found for warm climate snapshot; using PD observed winds instead')
      climate%matrix%GCM_warm%Wind_LR( :,:,grid%i1:grid%i2) = climate%matrix%PD_obs%wind_LR( :,:,grid%i1:grid%i2)
      climate%matrix%GCM_warm%Wind_DU( :,:,grid%i1:grid%i2) = climate%matrix%PD_obs%Wind_DU( :,:,grid%i1:grid%i2)
    END IF
    IF (.NOT. found_winds_cold) THEN
      IF (par%master) CALL warning('no wind fields found for cold climate snapshot; using PD observed winds instead')
      climate%matrix%GCM_cold%Wind_LR( :,:,grid%i1:grid%i2) = climate%matrix%PD_obs%wind_LR( :,:,grid%i1:grid%i2)
      climate%matrix%GCM_cold%Wind_DU( :,:,grid%i1:grid%i2) = climate%matrix%PD_obs%Wind_DU( :,:,grid%i1:grid%i2)
    END IF
    CALL sync

    ! Calculate spatially variable lapse rate

    ! Use a uniform value for the warm snapshot [this assumes "warm" is actually identical to PI!]
    climate%matrix%GCM_warm%lambda( :,grid%i1:grid%i2) = C%constant_lapserate
    CALL sync

    IF     (region_name == 'NAM' .OR. region_name == 'EAS') THEN
      CALL initialise_matrix_calc_spatially_variable_lapserate( grid, climate%matrix%GCM_PI, climate%matrix%GCM_cold)
    ELSEIF (region_name == 'GLR' .OR. region_name == 'ANT') THEN
      climate%matrix%GCM_cold%lambda( :,grid%i1:grid%i2) = C%constant_lapserate
      CALL sync
    END IF

    ! Calculate GCM bias
    CALL initialise_matrix_calc_GCM_bias( grid, climate%matrix%GCM_PI, climate%matrix%PD_obs, &
      climate%matrix%GCM_bias_T2m, climate%matrix%GCM_bias_Precip)

    ! Apply bias correction
    IF (C%climate_matrix_biascorrect_warm) CALL initialise_matrix_apply_bias_correction( grid, climate%matrix%GCM_warm, &
      climate%matrix%GCM_bias_T2m, climate%matrix%GCM_bias_Precip)
    IF (C%climate_matrix_biascorrect_warm) CALL initialise_matrix_apply_bias_correction( grid, climate%matrix%GCM_cold, &
      climate%matrix%GCM_bias_T2m, climate%matrix%GCM_bias_Precip)

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
  
end module climate_matrix