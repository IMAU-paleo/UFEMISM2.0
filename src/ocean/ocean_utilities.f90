MODULE ocean_utilities

  ! Realistic ocean models

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE control_resources_and_error_messaging                  , ONLY: crash, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  USE ocean_model_types                                      , ONLY: type_ocean_model

  IMPLICIT NONE

CONTAINS

! ===== Mixed layer =====
! =======================

  subroutine calc_ocean_temperature_at_shelf_base( mesh, ice, ocean)
    ! Calculate ocean temperature at the base of the shelf by interpolating
    ! the 3-D ocean temperature field in the vertical column

    implicit none

    ! In/output variables
    type(type_mesh),                    intent(in)    :: mesh
    type(type_ice_model),               intent(in)    :: ice
    type(type_ocean_model),             intent(inout) :: ocean

    ! Local variables:
    character(len=256), parameter                     :: routine_name = 'calc_ocean_temperature_at_shelf_base'
    integer                                           :: vi
    real(dp)                                          :: depth

    ! Add routine to path
    call init_routine( routine_name)

    do vi = mesh%vi1, mesh%vi2

      ! Initialise at zero
      ocean%T_draft( vi) = 0._dp

      ! Depth (positive downwards)
      depth = ice%SL( vi) - ice%Hib( vi)

      if (depth > 0._dp) then

        ! Find ocean temperature at this depth
        call interpolate_ocean_depth( C%nz_ocean, C%z_ocean, ocean%T( vi,:), depth, ocean%T_draft( vi))

      end if

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_ocean_temperature_at_shelf_base

  subroutine calc_ocean_freezing_point_at_shelf_base( mesh, ice, ocean)
    ! Calculate the ocean freezing point at the base of the shelf, needed to calculate
    ! basal melt in the different parameterisations from Favier et al. (2019)

    implicit none

    ! In/output variables
    type(type_mesh),                    intent(in)    :: mesh
    type(type_ice_model),               intent(in)    :: ice
    type(type_ocean_model),             intent(inout) :: ocean

    ! Local variables:
    character(len=256), parameter                     :: routine_name = 'calc_ocean_freezing_point_at_shelf_base'
    integer                                           :: vi
    real(dp)                                          :: depth
    real(dp)                                          :: S0                   ! Practical salinity [PSU]
    real(dp), parameter                               :: lambda1 = -0.0575_dp ! Liquidus slope                [degC PSU^-1] (Favier et al. (2019), Table 2)
    real(dp), parameter                               :: lambda2 = 0.0832_dp  ! Liquidus intercept            [degC]        (Favier et al. (2019), Table 2)
    real(dp), parameter                               :: lambda3 = 7.59E-4_dp ! Liquidus pressure coefficient [degC m^-1]   (Favier et al. (2019), Table 2)

    ! Add routine to path
    call init_routine( routine_name)

    do vi = mesh%vi1, mesh%vi2

      ! Initialise at zero
      ocean%T_freezing_point( vi) = 0._dp

      ! Depth (positive downwards)
      depth = ice%SL( vi) - ice%Hib( vi)

      if (depth > 0._dp) then

        ! Find salinity at this depth
        call interpolate_ocean_depth( C%nz_ocean, C%z_ocean, ocean%S( vi,:), depth, S0)

        ! Calculate ocean freezing temperature (Favier et al. (2019), Eq. 3) in degrees Celsius
        ocean%T_freezing_point( vi) = lambda1 * S0 + lambda2 - lambda3 * depth

      end if

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_ocean_freezing_point_at_shelf_base

  SUBROUTINE interpolate_ocean_depth( nz_ocean, z_ocean, f_ocean, z_query, f_query)
    ! Interpolate ocean column data to a queried depth using a simple bisection method.

    IMPLICIT NONE

    ! In/output variables:
    INTEGER,                             INTENT(IN)    :: nz_ocean    ! Number of vertical layers
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: z_ocean     ! Depth of layers (assumed to be monotonically increasing, does not need to be regularly spaced)
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: f_ocean     ! Value of whatever function we want to interpolate
    REAL(dp),                            INTENT(IN)    :: z_query     ! Depth at which we want to know the function
    REAL(dp),                            INTENT(OUT)   :: f_query     ! Interpolated function value

    ! Local variables:
    INTEGER                                            :: k_lo,k_hi,k_mid
    LOGICAL                                            :: foundit
    REAL(dp)                                           :: w

    ! Safety
    IF (z_query < 0._dp) THEN
      CALL crash('  interpolate_ocean_depth - ERROR: z_query = {dp_01} < 0; cannot extrapolate above the sea surface, obviously!', dp_01 = z_query)
    ELSEIF (z_query > 12000._dp) THEN
      CALL crash('  interpolate_ocean_depth - ERROR: z_query = {dp_01} > 12 km; the ocean is not that deep!', dp_01 = z_query)
    ELSEIF (SIZE(z_ocean,1) /= nz_ocean) THEN
      CALL crash('  interpolate_ocean_depth - ERROR: SIZE(z_ocean,1) = {int_01} /= nz_ocean = {int_02}!', int_01 = SIZE(z_ocean,1), int_02 = nz_ocean)
    ELSEIF (SIZE(f_ocean,1) /= nz_ocean) THEN
      CALL crash('  interpolate_ocean_depth - ERROR: SIZE(f_ocean,1) = {int_01} /= nz_ocean = {int_02}!', int_01 = SIZE(f_ocean,1), int_02 = nz_ocean)
    ELSEIF (z_query > MAXVAL(z_ocean)) THEN
      ! Nearest-neighbour extrapolation when querying data beneath the end of the ocean data column
      f_query = f_ocean( nz_ocean)
      RETURN
    END IF

    ! Exception for when z_query = 0 (the World Ocean Atlas depth starts at 1.25...)
    IF (z_query < MINVAL(z_ocean)) THEN
      f_query = f_ocean(1)
      RETURN
    END IF

    ! Bisection method
    k_lo  = 1
    k_hi  = nz_ocean
    k_mid = INT( REAL(k_lo + k_hi,dp) / 2._dp)

    ! Exceptions
    IF     (ABS(z_query - z_ocean( k_lo )) < 1E-4_dp) THEN
      f_query = f_ocean( k_lo)
      RETURN
    ELSEIF (ABS(z_query - z_ocean( k_hi )) < 1E-4_dp) THEN
      f_query = f_ocean( k_hi)
      RETURN
    ELSEIF (ABS(z_query - z_ocean( k_mid)) < 1E-4_dp) THEN
      f_query = f_ocean( k_mid)
      RETURN
    END IF

    ! Bisection method
    foundit = .FALSE.
    DO WHILE (.NOT. foundit)

      IF (ABS(z_query - z_ocean( k_mid)) < 1E-4_dp) THEN
        ! Exception for when the queried depth is exactly at the midpoint index depth
        f_query = f_ocean( k_mid)
        RETURN
      ELSEIF (z_query > z_ocean( k_mid)) THEN
        ! Queried depth lies to the right of the midpoint
        k_lo = k_mid
        k_mid = INT( REAL(k_lo + k_hi,dp) / 2._dp)
      ELSE
        ! Queried depth lies to the left of the midpoint
        k_hi = k_mid
        k_mid = INT( REAL(k_lo + k_hi,dp) / 2._dp)
      END IF

      ! Stop iterating when endpoints lie next to each other; then just do linear interpolation between those two.
      IF (k_hi == k_lo+1) foundit = .TRUE.

    END DO ! DO WHILE (.NOT. foundit)

    ! Get interpolated value
    if (isnan(f_ocean( k_hi)) .and. isnan(f_ocean( k_lo))) then
      ! Both NaNs, so output NaN
      f_query = f_ocean( k_hi)
    elseif (isnan(f_ocean( k_hi))) then
      ! Only deeper value is NaN, output non-NaN value
      f_query = f_ocean( k_lo)
    elseif (isnan(f_ocean( k_lo))) then
      ! Only shallower value is NaN, output non-NaN value
      f_query = f_ocean( k_hi)
    else
      ! Both values non-NaN
      ! Linear interpolation between nearest layers
      w = (z_query - z_ocean( k_lo)) / (z_ocean( k_hi) - z_ocean( k_lo))
      f_query = w * f_ocean( k_hi) + (1._dp - w) * f_ocean( k_lo)
    end if

  END SUBROUTINE interpolate_ocean_depth

! ===== Regridding ======
! =======================

  subroutine initialise_ocean_vertical_grid
    ! Set up the vertical grid used for ocean data - regular grid

    implicit none

    ! Local variables:
    character(len=256), parameter :: routine_name = 'initialise_ocean_vertical_grid'
    INTEGER                       :: k

    ! Add routine to path
    call init_routine( routine_name)

    ! Determine the number of vertical layers to be used
    C%nz_ocean = 1 + floor( C%ocean_vertical_grid_max_depth / C%ocean_vertical_grid_dz)

    ! Allocate memory
    allocate( C%z_ocean( C%nz_ocean))

    ! Fill in the values
    do k = 1, C%nz_ocean
      C%z_ocean( k) = real(k-1,dp) * C%ocean_vertical_grid_dz
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_ocean_vertical_grid

END MODULE ocean_utilities
