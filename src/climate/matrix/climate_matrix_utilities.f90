module climate_matrix_utilities
  ! check which of these is needed!
  use precisions                                             , only: dp
  use mpi_basic                                              , only: par, sync
  use control_resources_and_error_messaging                  , only: crash, init_routine, finalise_routine, colour_string, warning
  use model_configuration                                    , only: C
  use parameters
  use mpi_f08, only: MPI_ALLREDUCE, MPI_IN_PLACE, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, MPI_INTEGER
  use mesh_types                                             , only: type_mesh
  use ice_model_types                                        , only: type_ice_model
  use grid_types                                             , only: type_grid
  use climate_model_types                                    , only: type_climate_model, type_climate_model_matrix, type_climate_model_snapshot
  use global_forcing_types                                   , only: type_global_forcing
  use SMB_model_types, only: type_SMB_model
  use netcdf_io_main
  use mesh_disc_apply_operators, only: ddx_a_a_2D, ddy_a_a_2D
  use erf_mod, only: error_function
  use assertions_basic, only: assert

 implicit none

  private

  public :: allocate_climate_snapshot
  public :: read_climate_snapshot
  public :: adapt_precip_CC
  public :: adapt_precip_Roe
  !public :: rotate_wind_to_model_mesh ! I don't think is needed...
  public :: get_insolation_at_time
  public :: update_insolation_timeframes_from_file

  contains

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
    TYPE(type_climate_model_snapshot),  INTENT(INOUT) :: snapshot

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
#if (DO_ASSERTIONS)
    do m = 1, 12
    do vi = mesh%vi1, mesh%vi2
      call assert( snapshot%Precip(vi, m) >= 0, 'snapshot%Precip with negative values')
    end do
    end do
#endif

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_climate_snapshot

  ! Two different parameterised precipitation models:
  ! - a simply Clausius-Clapeyron-based method            (used for GRL and ANT)
  ! - the Roe & Lindzen temperature/orography-based model (used for NAM and EAS)
  SUBROUTINE adapt_precip_CC( mesh, Hs, Hs_GCM, T_ref_GCM, P_ref_GCM, Precip_GCM, region_name)

    USE parameters, ONLY: T0

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                                 INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2),          INTENT(IN)    :: Hs              ! Model orography (m)
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2),          INTENT(IN)    :: Hs_GCM          ! Reference orography (m)           - total ice-weighted
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2,12),       INTENT(IN)    :: T_ref_GCM       ! Reference temperature (K)         - total ice-weighted
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2,12),       INTENT(IN)    :: P_ref_GCM       ! Reference precipitation (m/month) - total ice-weighted
    CHARACTER(LEN=3),                                INTENT(IN)    :: region_name

    ! Output variables:
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2,12),       INTENT(INOUT)   :: Precip_GCM      ! Climate matrix precipitation

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
      T_inv(     vi,m) = 88.9_dp + 0.67_dp * (T_ref_GCM( vi,m) - C%climate_matrix_constant_lapserate * (Hs( vi) - Hs_GCM( vi)))
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

  SUBROUTINE get_insolation_at_time( mesh, time, snapshot)
    ! Get monthly insolation at time t on the regional grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_climate_model_snapshot),      INTENT(INOUT) :: snapshot
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                    :: routine_name = 'get_insolation_at_time'
    REAL(dp)                                         :: time_applied
    INTEGER                                          :: vi,m 
    REAL(dp)                                         :: wt0, wt1

    ! Add routine to path
    CALL init_routine( routine_name)

    time_applied = 0._dp

    ! Safety
    IF     (C%choice_insolation_forcing == 'none') THEN
      CALL crash('insolation should not be used when choice_insolation_forcing = "none"!')
    ELSEIF (C%choice_insolation_forcing == 'static') THEN
      time_applied = C%static_insolation_time
    ELSEIF (C%choice_insolation_forcing == 'realistic') THEN
      time_applied = time
    ELSE
      CALL crash('unknown choice_insolation_forcing "' // TRIM( C%choice_insolation_forcing) // '"!')
    END IF

    ! Check if the requested time is enveloped by the two timeframes;
    ! if not, read the two relevant timeframes from the NetCDF file
    IF (time_applied < snapshot%ins_t0 .OR. time_applied > snapshot%ins_t1) THEN
      IF (par%primary)  WRITE(0,*) '   Model time is out of the current insolation timeframes. Updating timeframes...'
      CALL update_insolation_timeframes_from_file( snapshot, time_applied, mesh)
    END IF

    ! Calculate timeframe interpolation weights (plus safety checks for when the extend beyond the record)
    if (snapshot%ins_t1 == snapshot%ins_t0) then
      wt0 = 0._dp
      wt1 = 1._dp
    else
      if (time_applied > snapshot%ins_t1) then
        wt0 = 0._dp
      elseif (time_applied < snapshot%ins_t0) then
        wt0 = 1._dp
      else
        wt0 = (snapshot%ins_t1 - time_applied) / (snapshot%ins_t1 - snapshot%ins_t0)
      end if
      wt1 = 1._dp - wt0
    end if

    ! Interpolate the two timeframes
    do vi = mesh%vi1, mesh%vi2
      do m = 1, 12
      !print *, "value of Q_TOA0 ", snapshot%ins_Q_TOA0(vi, m), "and Q_TOA1 ", snapshot%ins_Q_TOA1( vi, m)
        snapshot%Q_TOA(vi, m) = wt0 * snapshot%ins_Q_TOA0(vi, m) + wt1 * snapshot%ins_Q_TOA1(vi, m)
      end do
    end do

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE get_insolation_at_time

  SUBROUTINE update_insolation_timeframes_from_file( snapshot, time, mesh)
    ! Read the NetCDF file containing the insolation forcing data. Only read the time frames enveloping the current
    ! coupling timestep to save on memory usage. Only done by master.

    ! NOTE: assumes time in forcing file is in kyr

    IMPLICIT NONE

    TYPE(type_mesh),                  INTENT(IN)     :: mesh
    TYPE(type_climate_model_snapshot), INTENT(INOUT)  :: snapshot
    REAL(dp),                         INTENT(IN)     :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'update_insolation_timeframes_from_file'
    INTEGER                                            :: ti0, ti1, ncid
    CHARACTER(LEN=256)                                 :: str

    ! Add routine to path
    CALL init_routine( routine_name)

    IF     (C%choice_insolation_forcing == 'none') THEN
      CALL crash('insolation should not be used when choice_insolation_forcing = "none"!')
    ELSEIF (C%choice_insolation_forcing == 'static' .OR. &
            C%choice_insolation_forcing == 'realistic') THEN

      ! Update insolation
      ! Find time indices to be read
      !IF (par%primary) THEN

        call read_field_from_file_0D( C%filename_insolation, field_name_options_time, snapshot%ins_t0, time_to_read = time)
        
        ! if the desired time is after t0, we read one record after for t1
        if (time >= snapshot%ins_t0) then
          call read_field_from_file_0D( C%filename_insolation, field_name_options_time, snapshot%ins_t1, time_to_read = time+1000._dp)
        else
        ! otherwise we read one record before for t0, and that record becomes t1
          call read_field_from_file_0D( C%filename_insolation, field_name_options_time, snapshot%ins_t1, time_to_read = time)
          call read_field_from_file_0D( C%filename_insolation, field_name_options_time, snapshot%ins_t0, time_to_read = time-1000._dp)
        end if

      !END IF ! IF (par%primary) THEN
      call read_field_from_file_2D_monthly( C%filename_insolation, field_name_options_insolation, mesh, snapshot%ins_Q_TOA0, time_to_read = snapshot%ins_t0)
      call read_field_from_file_2D_monthly( C%filename_insolation, field_name_options_insolation, mesh, snapshot%ins_Q_TOA1, time_to_read = snapshot%ins_t1)

      call warning('insolation timeframes at t = {dp_01} are ins_t0={dp_02} and ins_t1={dp_03}', dp_01 =  time, dp_02 = snapshot%ins_t0, dp_03 = snapshot%ins_t1)

    ELSE
      CALL crash('unknown choice_insolation_forcing "' // TRIM( C%choice_insolation_forcing) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE update_insolation_timeframes_from_file

end module