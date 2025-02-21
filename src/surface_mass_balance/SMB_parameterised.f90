module SMB_Parameterised

  ! Parameterised SMB models

  use precisions, only: dp
  use mpi_basic, only: par, sync
  use control_resources_and_error_messaging, only: crash, init_routine, finalise_routine, colour_string
  use model_configuration, only: C
  use parameters
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model
  use SMB_model_types, only: type_SMB_model
  USE parameters,      only: T0, L_fusion, sec_per_year, pi, ice_density

  implicit none

  REAL(dp), PARAMETER :: albedo_water        = 0.1_dp
  REAL(dp), PARAMETER :: albedo_soil         = 0.2_dp
  REAL(dp), PARAMETER :: albedo_ice          = 0.5_dp
  REAL(dp), PARAMETER :: albedo_snow         = 0.85_dp

contains

  subroutine run_SMB_model_parameterised( mesh, ice, SMB, climate, time)
    ! Calculate the surface mass balance
    !
    ! use an idealised SMB scheme

    ! In/output variables:
    type(type_mesh),      intent(in)    :: mesh
    type(type_ice_model), intent(in)    :: ice
    type(type_SMB_model), intent(inout) :: SMB
    type(climate),        intent(in)    :: climate
    real(dp),             intent(in)    :: time

    ! Local variables:
    character(len=256), parameter :: routine_name = 'run_SMB_model_parameterised'
    integer                       :: vi

    ! Run the chosen parameterised SMB model
    SELECT CASE (C%choice_SMB_model_parameterised)
      CASE ('IMAU-ITM')
        CALL run_SMB_model_parameterised_IMAUITM( mesh, ice, climate, SMB)
      CASE DEFAULT
        CALL crash('unknown choice_SMB_model_parameterised "' // TRIM( C%choice_SMB_model_parameterised) // '"')
    END SELECT

    ! Add routine to path
    call init_routine( routine_name)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_SMB_model_parameterised

  subroutine run_SMB_model_parameterised_IMAUITM( mesh, ice, SMB, climate, SMB)
    ! Initialise the SMB model
    !
    ! use a parameterised SMB scheme

    ! In- and output variables
    type(type_mesh),          intent(in)    :: mesh
    type(type_ice_model),     intent(in)    :: ice
    type(type_climate_model), intent(in)    :: climate
    type(type_SMB_model),     intent(inout) :: SMB

    ! Local variables:
    character(len=256), parameter :: routine_name = 'run_SMB_model_parameterised_IMAUITM'
    integer                       :: vi
    integer                       :: m,mprev
    real(dp)                      :: snowfrac, liquid_water, sup_imp_wat
    real(dp)                      :: dummy_dp

    ! Add routine to path
    call init_routine( routine_name)

    ! Print to terminal
    if (par%master) write(*,"(a)") '   Initialising parameterised SMB model "' // &
      colour_string( trim( C%choice_SMB_model_parameterised),'light blue') // '"...'

    ! Run the chosen parameterised SMB model

    ! Initialise
    ! To prevent compiler warnings
    dummy_dp = mesh%xmin
    dummy_dp = SMB%SMB( mesh%vi1)

    DO vi = mesh%vi1, mesh%vi2
      ! Background albedo
      SMB%AlbedoSurf( vi) = albedo_soil
      IF ((ice%mask_ocean( vi) == 1 .AND. ice%mask_shelf( vi) == 0) .OR. ice%mask_noice( vi) == 1) SMB%AlbedoSurf( vi) = albedo_water
      IF (ice%mask_ice(    vi) == 1) SMB%AlbedoSurf( vi) = albedo_ice

      DO m = 1, 12  ! Month loop

        mprev = m - 1
        IF (mprev==0) mprev = 12

        SMB%Albedo( m,vi) = MIN(albedo_snow, MAX( SMB%AlbedoSurf( vi), albedo_snow - (albedo_snow - SMB%AlbedoSurf( vi))  * &
                             EXP(-15._dp * SMB%FirnDepth( mprev,vi)) - 0.015_dp * SMB%MeltPreviousYear( vi)))
        IF ((ice%mask_ocean( vi) == 1 .AND. ice%mask_shelf( vi) == 0) .OR. ice%mask_noice( vi) == 1) SMB%Albedo( m,vi) = albedo_water

        ! Determine ablation as a function of surface temperature and albedo/insolation according to Bintanja et al. (2002)
        SMB%Melt( m,vi) = MAX(0._dp, ( SMB%C_abl_Ts         * (climate%T2m( m,vi) - T0) + &
                                       SMB%C_abl_Q          * (1.0_dp - SMB%Albedo( m,vi)) * climate%Q_TOA( m,vi) - &
                                       SMB%C_abl_constant)  * sec_per_year / (L_fusion * 1000._dp * 12._dp))

        ! Determine accumulation with snow/rain fraction from Ohmura et al. (1999), liquid water content (rain and melt water) and snow depth

        ! NOTE: commented version is the old ANICE version, supposedly based on "physics" (which we cant check), but
        !       the new version was tuned to RACMO output and produced significantly better snow fractions...
        !snowfrac = MAX(0._dp, MIN(1._dp, 0.5_dp   * (1 - ATAN((climate%T2m(vi,m) - T0) / 3.5_dp)  / 1.25664_dp)))
        snowfrac = MAX(0._dp, MIN(1._dp, 0.725_dp * (1 - ATAN((climate%T2m( m,vi) - T0) / 5.95_dp) / 1.8566_dp)))

        SMB%Snowfall( m,vi) = climate%Precip( m,vi) *          snowfrac
        SMB%Rainfall( m,vi) = climate%Precip( m,vi) * (1._dp - snowfrac)

        ! Refreezing according to Janssens & Huybrechts (2000)
        ! The refreezing (=effective retention) is the minimum value of the amount of super imposed
        ! water and the available liquid water, with a maximum value of the total precipitation.
        ! (see also Huybrechts & de Wolde, 1999)

        ! Add this month's snow accumulation to next month's initial snow depth.
        SMB%AddedFirn( m,vi) = SMB%Snowfall( m,vi) - SMB%Melt( m,vi)
        SMB%FirnDepth( m,vi) = MIN(10._dp, MAX(0._dp, SMB%FirnDepth( mprev,vi) + SMB%AddedFirn( m,vi) ))

      END DO ! DO m = 1, 12

      ! Calculate refreezing for the whole year, divide equally over the 12 months, then calculate resulting runoff and SMB.
      ! This resolves the problem with refreezing where liquid water is mostly available in summer
      ! but "refreezing potential" mostly in winter, and there is no proper meltwater retention.
      sup_imp_wat  = SMB%C_refr * MAX(0._dp, T0 - SUM(climate%T2m( :,vi))/12._dp)
      liquid_water = SUM(SMB%Rainfall( :,vi)) + SUM(SMB%Melt( :,vi))

      SMB%Refreezing_year( vi) = MIN( MIN( sup_imp_wat, liquid_water), SUM(climate%Precip( :,vi)))
      IF (ice%mask_ice( vi)==0) SMB%Refreezing_year( vi) = 0._dp

      DO m = 1, 12
        SMB%Refreezing(  m,vi) = SMB%Refreezing_year( vi) / 12._dp
        SMB%Runoff(      m,vi) = SMB%Melt( m,vi) + SMB%Rainfall( m,vi) - SMB%Refreezing( m,vi)
        SMB%SMB_monthly( m,vi) = SMB%Snowfall( m,vi) + SMB%Refreezing( m,vi) - SMB%Melt( m,vi)
      END DO

      SMB%SMB( vi) = SUM(SMB%SMB_monthly( :,vi))

      ! Calculate total melt over this year, to be used for determining next year's albedo
      SMB%MeltPreviousYear( vi) = SUM(SMB%Melt( :,vi))
    END DO
    

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_SMB_model_parameterised_IMAUITM

SUBROUTINE initialise_SMB_model_parameterised( mesh, SMB, climate, region_name)
    ! Initialise the SMB model
    !
    ! Use a parameterised SMB scheme

    ! In- and output variables
    TYPE(type_mesh),          INTENT(IN)    :: mesh
    TYPE(type_climate_model), INTENT(IN)    :: climate
    TYPE(type_SMB_model),     INTENT(INOUT) :: SMB
    CHARACTER(LEN=3),         INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'initialise_SMB_model_parameterised'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to terminal
    IF (par%master)  WRITE(*,"(A)") '   Initialising parameterised SMB model "' // &
      colour_string( TRIM( C%choice_SMB_model_parameterised),'light blue') // '"...'

    ! Initialise the chosen parameterised SMB model
    SELECT CASE (C%choice_SMB_model_parameterised)
      CASE ('IMAU-ITM')
        CALL initialise_SMB_model_parameterised_IMAUITM( mesh, ice, SMB, climate, region_name)
      CASE DEFAULT
        CALL crash('unknown choice_SMB_model_parameterised "' // TRIM( C%choice_SMB_model_parameterised) // '"')
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_SMB_model_parameterised

  SUBROUTINE initialise_SMB_model_parameterised_IMAUITM( mesh, ice, SMB, climate, region_name)
    ! Allocate memory for the data fields of the SMB model.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),          INTENT(IN)    :: mesh
    TYPE(type_ice_model),     INTENT(IN)    :: ice
    TYPE(type_climate_model), INTENT(IN)    :: climate
    CHARACTER(LEN=3),         INTENT(IN)    :: region_name
    TYPE(type_SMB_model),     INTENT(INOUT) :: SMB

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_SMB_model_parameterised_IMAUITM'
    INTEGER                                            :: vi
    CHARACTER(LEN=256)                                 :: choice_SMB_IMAUITM_init_firn

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Data fields
    ALLOCATE(SMB%AlbedoSurf      (   mesh%vi1:mesh%vi2))
    ALLOCATE(SMB%MeltPreviousYear(   mesh%vi1:mesh%vi2))
    ALLOCATE(SMB%FirnDepth       (12,mesh%vi1:mesh%vi2))
    ALLOCATE(SMB%Rainfall        (12,mesh%vi1:mesh%vi2))
    ALLOCATE(SMB%Snowfall        (12,mesh%vi1:mesh%vi2))
    ALLOCATE(SMB%AddedFirn       (12,mesh%vi1:mesh%vi2))
    ALLOCATE(SMB%Melt            (12,mesh%vi1:mesh%vi2))
    ALLOCATE(SMB%Refreezing      (12,mesh%vi1:mesh%vi2))
    ALLOCATE(SMB%Refreezing_year (   mesh%vi1:mesh%vi2))
    ALLOCATE(SMB%Runoff          (12,mesh%vi1:mesh%vi2))
    ALLOCATE(SMB%Albedo          (12,mesh%vi1:mesh%vi2))
    ALLOCATE(SMB%Albedo_year     (   mesh%vi1:mesh%vi2))
    ALLOCATE(SMB%SMB_monthly     (12,mesh%vi1:mesh%vi2))
    ALLOCATE(SMB%SMB             (   mesh%vi1:mesh%vi2))

    ! Tuning parameters
    ALLOCATE( SMB%C_abl_constant)
    ALLOCATE( SMB%C_abl_Ts)
    ALLOCATE( SMB%C_abl_Q)
    ALLOCATE( SMB%C_refr)

    ! Determine which constants to use for this region
    IF     (region_name == 'NAM') THEN
      SMB%C_abl_constant           = C%SMB_IMAUITM_C_abl_constant_NAM
      SMB%C_abl_Ts                 = C%SMB_IMAUITM_C_abl_Ts_NAM
      SMB%C_abl_Q                  = C%SMB_IMAUITM_C_abl_Q_NAM
      SMB%C_refr                   = C%SMB_IMAUITM_C_refr_NAM
    ELSEIF (region_name == 'EAS') THEN
      SMB%C_abl_constant           = C%SMB_IMAUITM_C_abl_constant_EAS
      SMB%C_abl_Ts                 = C%SMB_IMAUITM_C_abl_Ts_EAS
      SMB%C_abl_Q                  = C%SMB_IMAUITM_C_abl_Q_EAS
      SMB%C_refr                   = C%SMB_IMAUITM_C_refr_EAS
    ELSEIF (region_name == 'GRL') THEN
      SMB%C_abl_constant           = C%SMB_IMAUITM_C_abl_constant_GRL
      SMB%C_abl_Ts                 = C%SMB_IMAUITM_C_abl_Ts_GRL
      SMB%C_abl_Q                  = C%SMB_IMAUITM_C_abl_Q_GRL
      SMB%C_refr                   = C%SMB_IMAUITM_C_refr_GRL
    ELSEIF (region_name == 'ANT') THEN
      SMB%C_abl_constant           = C%SMB_IMAUITM_C_abl_constant_ANT
      SMB%C_abl_Ts                 = C%SMB_IMAUITM_C_abl_Ts_ANT
      SMB%C_abl_Q                  = C%SMB_IMAUITM_C_abl_Q_ANT
      SMB%C_refr                   = C%SMB_IMAUITM_C_refr_ANT
    END IF


    ! Initialisation choice
    IF     (region_name == 'NAM') THEN
      choice_SMB_IMAUITM_init_firn = C%choice_SMB_IMAUITM_init_firn_NAM
    ELSEIF (region_name == 'EAS') THEN
      choice_SMB_IMAUITM_init_firn = C%choice_SMB_IMAUITM_init_firn_EAS
    ELSEIF (region_name == 'GRL') THEN
      choice_SMB_IMAUITM_init_firn = C%choice_SMB_IMAUITM_init_firn_GRL
    ELSEIF (region_name == 'ANT') THEN
      choice_SMB_IMAUITM_init_firn = C%choice_SMB_IMAUITM_init_firn_ANT
    END IF

    ! Initialise the firn layer
    IF     (choice_SMB_IMAUITM_init_firn == 'uniform') THEN
      ! Initialise with a uniform firn layer over the ice sheet

      DO vi = mesh%vi1, mesh%vi2
        IF (ice%Hi( vi) > 0._dp) THEN
          SMB%FirnDepth(        :,vi) = C%SMB_IMAUITM_initial_firn_thickness
          SMB%MeltPreviousYear(   vi) = 0._dp
        ELSE
          SMB%FirnDepth(        :,vi) = 0._dp
          SMB%MeltPreviousYear(   vi) = 0._dp
        END IF
      END DO

    ELSEIF (choice_SMB_IMAUITM_init_firn == 'read_from_file') THEN
      ! Initialise with the firn layer of a previous run
      CALL initialise_IMAUITM_firn_from_file( mesh, SMB, region_name)
    ELSE
      CALL crash('unknown choice_SMB_IMAUITM_init_firn "' // TRIM( SMB_IMAUITM_choice_init_firn) // '"!')
    END IF

    ! Initialise albedo
    DO vi = mesh%vi1, mesh%vi2
      ! Background albedo
      IF (ice%Hb_a( vi) < 0._dp) THEN
        SMB%AlbedoSurf( vi) = albedo_water
      ELSE
        SMB%AlbedoSurf( vi) = albedo_soil
      END IF
      IF (ice%Hi( vi) > 0._dp) THEN
        SMB%AlbedoSurf(  vi) = albedo_snow
      END IF
      SMB%Albedo( :,vi) = SMB%AlbedoSurf( vi)
    END DO
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)
  END SUBROUTINE initialise_climate_model_parameterised_IMAUITM

  SUBROUTINE initialise_IMAUITM_firn_from_file( mesh, SMB, region_name)
    ! If this is a restarted run, read the firn depth and meltpreviousyear data from the restart file

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_SMB_model),                INTENT(INOUT) :: SMB
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_IMAUITM_firn_restart'
    CHARACTER(LEN=256)                                 :: filename_restart_firn
    REAL(dp)                                           :: timeframe_restart_firn
    TYPE(type_restart_data)                            :: restart

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Assume that SMB and geometry are read from the same restart file
    SELECT CASE (region_name)
    CASE('NAM') 
      filename_restart_firn = C%filename_firn_IMAUITM_NAM
      timeframe_restart_firn = C%timeframe_firn_IMAUITM_NAM
    CASE('EAS') 
      filename_restart_firn = C%filename_firn_IMAUITM_EAS
      timeframe_restart_firn = C%timeframe_firn_IMAUITM_EAS
    CASE('GRL') 
      filename_restart_firn = C%filename_firn_IMAUITM_GRL
      timeframe_restart_firn = C%timeframe_firn_IMAUITM_GRL
    CASE('ANT') 
      filename_restart_firn = C%filename_firn_IMAUITM_ANT
      timeframe_restart_firn = C%timeframe_firn_IMAUITM_ANT
    CASE DEFAULT
        CALL crash('unknown region_name "' // TRIM( region_name) // '"!')
    END SELECT

     ! Print to terminal
    IF (par%master)  WRITE(*,"(A)") '   Initialising SMB-model firn layer from file "' // colour_string( TRIM(filename_restart_SMB),'light blue') // '"...'


    ! Read firn layer from file
    IF (timeframe_restart_firn == 1E9_dp) THEN
      ! Assume the file has no time dimension
      CALL read_field_from_file_2D_monthly( filename_restart_firn, 'FirnDepth', mesh, SMB%FirnDepth)
      CALL read_field_from_file_2D( filename_restart_firn, 'MeltPreviousYear', mesh, SMB%MeltPreviousYear)
    ELSE
      ! Assume the file has a time dimension, and read the specified timeframe
      CALL read_field_from_file_2D_monthly( filename_restart_firn, 'FirnDepth', mesh, SMB%FirnDepth, time_to_read = timeframe_restart_firn)
      CALL read_field_from_file_2D( filename_restart_firn, 'MeltPreviousYear', mesh, SMB%MeltPreviousYear, time_to_read = timeframe_restart_firn)
    END IF


    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_IMAUITM_firn_from_file

end module SMB_parameterised
