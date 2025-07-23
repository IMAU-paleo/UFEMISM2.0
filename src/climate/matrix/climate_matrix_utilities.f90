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
  subroutine allocate_climate_snapshot( mesh, snapshot, name)
    ! Allocate shared memory for a single climate snapshot

    ! In/output variables:
    type(type_mesh),                     intent(in)    :: mesh
    type(type_climate_model_snapshot),   intent(inout) :: snapshot
    character(len=*),                    intent(in)    :: name

    ! Local variables:
    character(len=256), PARAMETER                      :: routine_name = 'allocate_climate_snapshot'

    ! Add routine to path
    call init_routine( routine_name)

    snapshot%name = name

    allocate( snapshot%Hs( mesh%vi1:mesh%vi2))
    allocate( snapshot%T2m(     mesh%vi1:mesh%vi2, 12))
    allocate( snapshot%Precip(  mesh%vi1:mesh%vi2, 12))
    allocate( snapshot%Wind_WE( mesh%vi1:mesh%vi2, 12))
    allocate( snapshot%Wind_SN( mesh%vi1:mesh%vi2, 12))
    allocate( snapshot%Wind_LR( mesh%vi1:mesh%vi2, 12))
    allocate( snapshot%Wind_DU( mesh%vi1:mesh%vi2, 12))
    allocate( snapshot%lambda( mesh%vi1:mesh%vi2))
    
    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine allocate_climate_snapshot
  subroutine read_climate_snapshot( filename, mesh, snapshot)
    ! Read a climate snapshot from a NetCDF file. Works both for global lon/lat files and regional x/y files

    ! In/output variables:
    character(len=256),                 intent(in)    :: filename
    type(type_mesh),                    intent(in)    :: mesh
    type(type_climate_model_snapshot),  intent(inout) :: snapshot

    ! Local variables:
    character(len=256), parameter                     :: routine_name = 'read_climate_snapshot'
    integer                                           :: vi, m

    ! Add routine to path
    call init_routine( routine_name)

    ! Write message to screen
    if (par%primary) write(0,*) '  Reading climate for snapshot "' // trim( snapshot%name) // '" from file ' // trim( filename)

    call read_field_from_file_2D(         filename, field_name_options_Hs , mesh, C%output_dir, snapshot%Hs     )
    call read_field_from_file_2D_monthly( filename, 'T2m'                 , mesh, C%output_dir, snapshot%T2m    )
    call read_field_from_file_2D_monthly( filename, 'Precip'              , mesh, C%output_dir, snapshot%Precip )
    call read_field_from_file_2D_monthly( filename, 'Wind_WE||uas||'      , mesh, C%output_dir, snapshot%Wind_WE) ! is needed the last ||? I copy it from SMB_realistic
    call read_field_from_file_2D_monthly( filename, 'Wind_SN||vas||'      , mesh, C%output_dir, snapshot%Wind_SN)
    
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
    call finalise_routine( routine_name)

  end subroutine read_climate_snapshot

  ! Two different parameterised precipitation models:
  ! - a simply Clausius-Clapeyron-based method            (used for GRL and ANT)
  ! - the Roe & Lindzen temperature/orography-based model (used for NAM and EAS)
  subroutine adapt_precip_CC( mesh, Hs, Hs_GCM, T_ref_GCM, P_ref_GCM, Precip_GCM, region_name)

    use parameters, only: T0

    ! Input variables:
    type(type_mesh),                                 intent(in)    :: mesh
    real(dp), dimension(mesh%vi1:mesh%vi2),          intent(in)    :: Hs              ! Model orography (m)
    real(dp), dimension(mesh%vi1:mesh%vi2),          intent(in)    :: Hs_GCM          ! Reference orography (m)           - total ice-weighted
    real(dp), dimension(mesh%vi1:mesh%vi2,12),       intent(in)    :: T_ref_GCM       ! Reference temperature (K)         - total ice-weighted
    real(dp), dimension(mesh%vi1:mesh%vi2,12),       intent(in)    :: P_ref_GCM       ! Reference precipitation (m/month) - total ice-weighted
    character(len=3),                                intent(in)    :: region_name

    ! Output variables:
    real(dp), dimension(mesh%vi1:mesh%vi2,12),       intent(inout)   :: Precip_GCM      ! Climate matrix precipitation

    ! Local variables:
    character(len=256), parameter                      :: routine_name = 'adapt_precip_CC'
    integer                                            :: vi,m
    real(dp), dimension(:, :), allocatable            ::  T_inv,  T_inv_ref

    ! Add routine to path
    call init_routine( routine_name)

    ! Allocate shared memory
    allocate( T_inv(     mesh%vi1:mesh%vi2, 12))
    allocate( T_inv_ref( mesh%vi1:mesh%vi2, 12))

    ! Calculate inversion layer temperatures
    do m = 1, 12
    do vi = mesh%vi1, mesh%vi2
      T_inv_ref( vi,m) = 88.9_dp + 0.67_dp *  T_ref_GCM( vi,m)
      T_inv(     vi,m) = 88.9_dp + 0.67_dp * (T_ref_GCM( vi,m) - C%climate_matrix_constant_lapserate * (Hs( vi) - Hs_GCM( vi)))
    end do
    end do
    call sync

    if     (region_name == 'GRL') then
      ! Method of Jouzel and Merlivat (1984), see equation (4.82) in Huybrechts (1992)

      do m = 1, 12
      do vi = mesh%vi1, mesh%vi2
        Precip_GCM( vi,m) = P_ref_GCM( vi,m) * 1.04**(T_inv( vi,m) - T_inv_ref( vi,m))
      end do
      end do
      call sync

    elseif (region_name == 'ANT') then
      ! As with Lorius/Jouzel method (also Huybrechts, 2002

      do m = 1, 12
      do vi = mesh%vi1, mesh%vi2

        Precip_GCM( vi,m) = P_ref_GCM( vi,m) * (T_inv_ref( vi,m) / T_inv( vi,m))**2 * EXP(22.47_dp * (T0 / T_inv_ref( vi,m) - T0 / T_inv( vi,m)))

      end do
      end do
      call sync

    else
      if (par%primary) then
        call crash('ERROR - adapt_precip_CC should only be used for Greenland and Antarctica!')
      end if
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine adapt_precip_CC
  subroutine adapt_precip_Roe( mesh, Hs1, T2m1, Wind_LR1, Wind_DU1, Precip1, &
                                      Hs2, T2m2, Wind_LR2, Wind_DU2, Precip2)
    ! Adapt precipitation from reference state 1 to model state 2, using the Roe&Lindzen precipitation model

    ! In/output variables:
    type(type_mesh),                                 intent(in)    :: mesh
    real(dp), dimension(mesh%vi1:mesh%vi2),          intent(in)    :: Hs1,      Hs2
    real(dp), dimension(mesh%vi1:mesh%vi2, 12),      intent(in)    :: T2m1,     T2m2
    real(dp), dimension(mesh%vi1:mesh%vi2, 12),      intent(in)    :: Wind_LR1, Wind_LR2
    real(dp), dimension(mesh%vi1:mesh%vi2, 12),      intent(in)    :: Wind_DU1, Wind_DU2
    real(dp), dimension(mesh%vi1:mesh%vi2, 12),      intent(in)    :: Precip1
    real(dp), dimension(mesh%vi1:mesh%vi2, 12),      intent(out)   :: Precip2

    ! Local variables:
    character(len=256), parameter                      :: routine_name = 'adapt_precip_Roe'
    integer                                            :: vi,m
    real(dp), dimension(:    ), allocatable            ::  dHs_dx1,  dHs_dx2
    real(dp), dimension(:    ), allocatable            ::  dHs_dy1,  dHs_dy2
    real(dp), dimension(:,:  ), allocatable            ::  Precip_RL1,  Precip_RL2,  dPrecip_RL

    ! Add routine to path
    call init_routine( routine_name)

    ! Allocate shared memory
    allocate( dHs_dx1( mesh%vi1:mesh%vi2))
    allocate( dHs_dx2( mesh%vi1:mesh%vi2))
    allocate( dHs_dy1( mesh%vi1:mesh%vi2))
    allocate( dHs_dy2( mesh%vi1:mesh%vi2))
    allocate( Precip_RL1( mesh%vi1:mesh%vi2, 12))
    allocate( Precip_RL2( mesh%vi1:mesh%vi2, 12))
    allocate( dPrecip_RL( mesh%vi1:mesh%vi2, 12))

    ! Calculate surface slopes for both states
    call ddx_a_a_2D( mesh, Hs1, dHs_dx1)
    call ddx_a_a_2D( mesh, Hs2, dHs_dx2)
    call ddy_a_a_2D( mesh, Hs1, dHs_dy1)
    call ddy_a_a_2D( mesh, Hs2, dHs_dy2)

    do vi = mesh%vi1, mesh%vi2
    do m = 1, 12

      ! Calculate precipitation with the Roe&Lindzen model for both states
      call precipitation_model_Roe( T2m1( vi,m), dHs_dx1( vi), dHs_dy1( vi), Wind_LR1( vi,m), Wind_DU1( vi,m), Precip_RL1( vi,m))
      call precipitation_model_Roe( T2m2( vi,m), dHs_dx2( vi), dHs_dy2( vi), Wind_LR2( vi,m), Wind_DU2( vi,m), Precip_RL2( vi,m))

      ! Calculate the ratio between those two precipitation rates
      dPrecip_RL( vi,m) = max(0.01_dp, min( 2._dp, Precip_RL2( vi,m) / Precip_RL1( vi,m) ))

      ! Applied model precipitation = (matrix-interpolated GCM reference precipitation) * RL ratio
      Precip2( vi,m) = Precip1( vi,m) * dPrecip_RL( vi,m)

    end do
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine adapt_precip_Roe


  subroutine precipitation_model_Roe( T2m, dHs_dx, dHs_dy, Wind_LR, Wind_DU, Precip)
    ! Precipitation model of Roe (J. Glac, 2002), integration from Roe and Lindzen (J. Clim. 2001)

    use parameters, only: T0, pi, sec_per_year

    ! In/output variables:
    real(dp),                            intent(in)    :: T2m
    real(dp),                            intent(in)    :: dHs_dx
    real(dp),                            intent(in)    :: dHs_dy
    real(dp),                            intent(in)    :: Wind_LR
    real(dp),                            intent(in)    :: Wind_DU
    real(dp),                            intent(out)   :: Precip

    ! Local variables:
    character(len=256), parameter                      :: routine_name = 'precipitation_model_Roe'
    real(dp)                                           :: upwind_slope, E_sat, x0, err

    real(dp), parameter                                :: e_sat0  = 611.2_dp
    real(dp), parameter                                :: c_one   = 17.67_dp
    real(dp), parameter                                :: c_two   = 243.5_dp

    real(dp), parameter                                :: a_par   = 2.5E-11_dp
    real(dp), parameter                                :: b_par   = 5.9E-09_dp
    real(dp), parameter                                :: alpha   = 100.0_dp

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate the upwind slope
    upwind_slope = max(0._dp, Wind_LR * dHs_dx + Wind_DU * dHs_dy)

    ! Calculate the saturation vapour pressure E_sat:
    E_sat = e_sat0 * exp( c_one * (T2m - T0) / (c_two + T2m - T0) )

    ! Calculate integration parameter x0 = a/b + w
    x0 = a_par / b_par + upwind_slope

    ! Calculate the error function argument
    err = alpha * abs(x0)

    ! Calculate precipitation rate
    Precip = ( b_par * E_sat ) * ( x0 / 2._dp + x0**2 * error_function(err) / (2._dp * abs(x0)) + &
                                         exp (-alpha**2 * x0**2) / (2._dp * sqrt(pi) * alpha) ) * sec_per_year

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine precipitation_model_Roe
  ! Rotate wind_WE, wind_SN to wind_LR, wind_DU
  subroutine rotate_wind_to_model_mesh( mesh, wind_WE, wind_SN, wind_LR, wind_DU)
    ! Code copied from ANICE.

    use parameters, only: pi

    ! In/output variables:
    type(type_mesh),                                   intent(in)    :: mesh
    real(dp), dimension(mesh%vi1:mesh%vi2,12),         intent(in)    :: wind_WE
    real(dp), dimension(mesh%vi1:mesh%vi2,12),         intent(in)    :: wind_SN
    real(dp), dimension(mesh%vi1:mesh%vi2,12),         intent(out)   :: wind_LR
    real(dp), dimension(mesh%vi1:mesh%vi2,12),         intent(out)   :: wind_DU

    ! Local variables:
    character(len=256), parameter                      :: routine_name = 'rotate_wind_to_model_mesh'
    integer                                            :: vi,m
    real(dp)                                           :: longitude_start, Uwind_x, Uwind_y, Vwind_x, Vwind_y

    ! Add routine to path
    call init_routine( routine_name)

    ! First find the first longitude which defines the start of quadrant I:
    longitude_start = mesh%lambda_M - 90._dp

    do vi = mesh%vi1, mesh%vi2
    do m = 1, 12
      ! calculate x and y from the zonal wind
      Uwind_x =   wind_WE( vi,m) * sin((pi/180._dp) * (mesh%lon( vi) - longitude_start))
      Uwind_y = - wind_WE( vi,m) * cos((pi/180._dp) * (mesh%lon( vi) - longitude_start))

      ! calculate x and y from the meridional winds
      Vwind_x =   wind_SN( vi,m) * cos((pi/180._dp) * (mesh%lon( vi) - longitude_start))
      Vwind_y =   wind_SN( vi,m) * sin((pi/180._dp) * (mesh%lon( vi) - longitude_start))

      ! Sum up wind components
      wind_LR( vi,m) = Uwind_x + Vwind_x   ! winds left to right
      wind_DU( vi,m) = Uwind_y + Vwind_y   ! winds bottom to top

    end do
    end do

    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine rotate_wind_to_model_mesh

  subroutine get_insolation_at_time( mesh, time, snapshot)
    ! Get monthly insolation at time t on the regional grid

    ! In/output variables:
    type(type_mesh),                        intent(in)    :: mesh
    type(type_climate_model_snapshot),      intent(inout) :: snapshot
    real(dp),                               intent(in)    :: time

    ! Local variables:
    character(len=256), parameter                    :: routine_name = 'get_insolation_at_time'
    real(dp)                                         :: time_applied
    integer                                          :: vi,m 
    real(dp)                                         :: wt0, wt1

    ! Add routine to path
    call init_routine( routine_name)

    time_applied = 0._dp

    ! Safety
    if     (C%choice_insolation_forcing == 'none') then
      call crash('insolation should not be used when choice_insolation_forcing = "none"!')
    elseif (C%choice_insolation_forcing == 'static') then
      time_applied = C%static_insolation_time
    elseif (C%choice_insolation_forcing == 'realistic') then
      time_applied = time
    else
      call crash('unknown choice_insolation_forcing "' // trim( C%choice_insolation_forcing) // '"!')
    end if

    ! Check if the requested time is enveloped by the two timeframes;
    ! if not, read the two relevant timeframes from the NetCDF file
    if (time_applied < snapshot%ins_t0 .or. time_applied > snapshot%ins_t1) then
      if (par%primary)  write(0,*) '   Model time is out of the current insolation timeframes. Updating timeframes...'
      call update_insolation_timeframes_from_file( snapshot, time_applied, mesh)
    end if

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
        snapshot%Q_TOA(vi, m) = wt0 * snapshot%ins_Q_TOA0(vi, m) + wt1 * snapshot%ins_Q_TOA1(vi, m)
      end do
    end do

  ! Finalise routine path
  call finalise_routine( routine_name)

  end subroutine get_insolation_at_time

  subroutine update_insolation_timeframes_from_file( snapshot, time, mesh)
    ! Read the NetCDF file containing the insolation forcing data. Only read the time frames enveloping the current
    ! coupling timestep to save on memory usage. Only done by master.

    ! NOTE: assumes time in forcing file is in kyr

    type(type_mesh),                   intent(in)     :: mesh
    type(type_climate_model_snapshot), intent(inout)  :: snapshot
    real(dp),                          intent(in)     :: time

    ! Local variables:
    character(len=256), parameter                      :: routine_name = 'update_insolation_timeframes_from_file'
    integer                                            :: ti0, ti1, ncid
    character(len=256)                                 :: str

    ! Add routine to path
    call init_routine( routine_name)

    if     (C%choice_insolation_forcing == 'none') then
      call crash('insolation should not be used when choice_insolation_forcing = "none"!')
    elseif (C%choice_insolation_forcing == 'static' .or. &
            C%choice_insolation_forcing == 'realistic') then

      ! Update insolation
      call read_field_from_file_0D( C%filename_insolation, field_name_options_time, snapshot%ins_t0, time_to_read = time)

      if (time >= snapshot%ins_t0) then
        call read_field_from_file_0D( C%filename_insolation, field_name_options_time, snapshot%ins_t1, time_to_read = time+1000._dp)
      else
        call read_field_from_file_0D( C%filename_insolation, field_name_options_time, snapshot%ins_t1, time_to_read = time)
        call read_field_from_file_0D( C%filename_insolation, field_name_options_time, snapshot%ins_t0, time_to_read = time-1000._dp)
      end if

      call read_field_from_file_2D_monthly( C%filename_insolation, field_name_options_insolation, mesh, C%output_dir, snapshot%ins_Q_TOA0, time_to_read = snapshot%ins_t0)
      call read_field_from_file_2D_monthly( C%filename_insolation, field_name_options_insolation, mesh, C%output_dir, snapshot%ins_Q_TOA1, time_to_read = snapshot%ins_t1)

      call warning('insolation timeframes at t = {dp_01} are ins_t0={dp_02} and ins_t1={dp_03}', dp_01 =  time, dp_02 = snapshot%ins_t0, dp_03 = snapshot%ins_t1)

    else
      call crash('unknown choice_insolation_forcing "' // trim( C%choice_insolation_forcing) // '"!')
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine update_insolation_timeframes_from_file

end module