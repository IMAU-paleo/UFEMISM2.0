MODULE thermodynamics_3D_heat_equation

  ! All the routines to solve the 3-D heat equation in the ice

! ===== Preamble =====
! ====================

  use mpi_f08, only: MPI_COMM_WORLD, MPI_ALLREDUCE, MPI_IN_PLACE, MPI_INTEGER, MPI_SUM
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  USE climate_model_types                                    , ONLY: type_climate_model
  USE SMB_model_types                                        , ONLY: type_SMB_model
  USE BMB_model_types                                        , ONLY: type_BMB_model
  use zeta_gradients, only: calc_zeta_gradients
  USE thermodynamics_utilities                               , ONLY: calc_heat_capacity, calc_thermal_conductivity, calc_pressure_melting_point, &
                                                                     calc_upwind_heat_flux_derivatives, calc_strain_heating, calc_frictional_heating, &
                                                                     replace_Ti_with_robin_solution
  use tridiagonal_solver, only: solve_tridiagonal_matrix_equation
  use netcdf_io_main

  IMPLICIT NONE

CONTAINS

! ===== Main routines =====
! =========================

  SUBROUTINE solve_3D_heat_equation( mesh, ice, climate, SMB, dt)
    ! Solve the three-dimensional heat equation
    !
    ! (See solve_1D_heat_equation for the derivation)
    !
    ! Uses time-step reduction on a per-vertex basis; for each vertex, it
    ! checks for instability in T(t+dt). If that is detected, it tries
    ! again with dt' = dt/2, repeatedly halving the time step until the
    ! instability disappears. This is done separately for every vertex.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                      INTENT(INOUT) :: mesh
    TYPE(type_ice_model),                 INTENT(INOUT) :: ice
    TYPE(type_climate_model),             INTENT(IN)    :: climate
    TYPE(type_SMB_model),                 INTENT(IN)    :: SMB
    REAL(dp),                             INTENT(IN)    :: dt

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'solve_3D_heat_equation'
    integer                                            :: ierr
    INTEGER                                            :: vi, k
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: u_times_dTdxp_upwind, v_times_dTdyp_upwind
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: T_surf_annual, Q_base_grnd, T_base_float
    REAL(dp)                                           :: dt_applied
    INTEGER                                            :: it_dt, it_it_dt
    LOGICAL                                            :: found_stable_solution
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: Ti_tplusdt
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: is_unstable
    INTEGER                                            :: n_unstable

    REAL(dp), DIMENSION( C%nz)                         :: icecol_Ti                   ! Vertical profile of ice temperature at time t                     [K]
    REAL(dp), DIMENSION( C%nz)                         :: icecol_u                    !   "         "    "  horizontal ice velocity in the x-direction    [m yr^-1]
    REAL(dp), DIMENSION( C%nz)                         :: icecol_v                    !   "         "    "  horizontal ice velocity in the y-direction    [m yr^-1]
    REAL(dp), DIMENSION( C%nz)                         :: icecol_w                    !   "         "    "  vertical   ice velocity                       [m yr^-1]
    REAL(dp), DIMENSION( C%nz)                         :: icecol_u_times_dTdxp_upwind !   "         "    "  u * dT/dxp in the upwind direction            [K yr^-1]
    REAL(dp), DIMENSION( C%nz)                         :: icecol_v_times_dTdyp_upwind !   "         "    "  u * dT/dxp in the upwind direction            [K yr^-1]
    REAL(dp), DIMENSION( C%nz)                         :: icecol_Ti_pmp               !   "         "    "  pressure melting point temperature            [K]
    REAL(dp), DIMENSION( C%nz)                         :: icecol_Ki                   !   "         "    "  thermal conductivity of ice                   [J yr^-1 m^-1 K^-1]
    REAL(dp), DIMENSION( C%nz)                         :: icecol_Cpi                  !   "         "    "  specific heat of ice                          [J kg^-1 K^-1]
    REAL(dp), DIMENSION( C%nz)                         :: icecol_dzeta_dx             !   "         "    "  x-derivative of the scaled coordinate zeta    [m^-1]
    REAL(dp), DIMENSION( C%nz)                         :: icecol_dzeta_dy             !   "         "    "  y-derivative of the scaled coordinate zeta    [m^-1]
    REAL(dp), DIMENSION( C%nz)                         :: icecol_dzeta_dz             !   "         "    "  z-derivative of the scaled coordinate zeta    [m^-1]
    REAL(dp), DIMENSION( C%nz)                         :: icecol_dzeta_dt             !   "         "    "  time-derivative of the scaled coordinate zeta [yr^-1]
    REAL(dp), DIMENSION( C%nz)                         :: icecol_Phi                  !   "         "    "  internal heat production                      [J kg^1 yr^1]
    REAL(dp), DIMENSION( C%nz)                         :: icecol_Ti_tplusdt           ! Vertical profile of ice temperature at time t + dt                [K]
    REAL(dp), DIMENSION( C%nz)                         :: icecol_Ti_tplusdt_gl_fl     ! Vertical profile of ice temperature at time t + dt (floating GL)  [K]

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate memory
    ALLOCATE( u_times_dTdxp_upwind( mesh%vi1:mesh%vi2,mesh%nz), source = 0._dp)
    ALLOCATE( v_times_dTdyp_upwind( mesh%vi1:mesh%vi2,mesh%nz), source = 0._dp)
    ALLOCATE( T_surf_annual       ( mesh%vi1:mesh%vi2        ), source = 0._dp)
    ALLOCATE( Q_base_grnd         ( mesh%vi1:mesh%vi2        ), source = 0._dp)
    ALLOCATE( T_base_float        ( mesh%vi1:mesh%vi2        ), source = 0._dp)
    ALLOCATE( Ti_tplusdt          ( mesh%vi1:mesh%vi2,mesh%nz), source = 0._dp)
    ALLOCATE( is_unstable         ( mesh%vi1:mesh%vi2        ), source = 0    )

    ! Calculate zeta gradients
    CALL calc_zeta_gradients( mesh, ice)

    ! Calculate temperature-dependent heat capacity
    CALL calc_heat_capacity( mesh, ice)

    ! Calculate temperature-dependent thermal conductivity
    CALL calc_thermal_conductivity( mesh, ice)

    ! Calculate pressure melting point
    CALL calc_pressure_melting_point( mesh, ice)

    ! Calculate upwind velocity times temperature gradients
    CALL calc_upwind_heat_flux_derivatives( mesh, ice, u_times_dTdxp_upwind, v_times_dTdyp_upwind)

    ! Calculate heating terms
    CALL calc_strain_heating(     mesh, ice)
    CALL calc_frictional_heating( mesh, ice)

    ! Calculate annual mean surface temperature
    DO vi = mesh%vi1, mesh%vi2
      T_surf_annual( vi) = SUM( climate%T2m( vi,:)) / REAL( SIZE( climate%T2m,2),dp)
    END DO

    ! For floating ice, basal temperatures are assumed to be always at
    ! the pressure melting point (since ocean water cannot be colder than
    ! this, and the ice itself cannot be warmer than this).
    DO vi = mesh%vi1, mesh%vi2
     T_base_float( vi) = ice%Ti_pmp( vi,C%nz)
    END DO

    ! Calculate heat flux at the base of the grounded ice
    DO vi = mesh%vi1, mesh%vi2
      Q_base_grnd( vi) = ice%frictional_heating( vi) + ice%geothermal_heat_flux( vi)
    END DO

    ! Solve the heat equation for all vertices
    is_unstable = 0
    n_unstable  = 0
    DO vi = mesh%vi1, mesh%vi2

      ! For very thin ice, just let the profile equal the surface temperature
      IF (ice%Hi_eff( vi) < C%Hi_min_thermo) THEN
        is_unstable( vi) = 0
        Ti_tplusdt( vi,:) = T_surf_annual( vi)
        CYCLE
      END IF

      ! Gather all the data needed to solve the 1-D heat equation
      icecol_Ti                   = ice%Ti(               vi,:)
      icecol_u                    = ice%u_3D(             vi,:)
      icecol_v                    = ice%v_3D(             vi,:)
      icecol_w                    = ice%w_3D(             vi,:)
      icecol_u_times_dTdxp_upwind = u_times_dTdxp_upwind( vi,:)
      icecol_v_times_dTdyp_upwind = v_times_dTdyp_upwind( vi,:)
      icecol_Ti_pmp               = ice%Ti_pmp(           vi,:)
      icecol_Ki                   = ice%Ki(               vi,:)
      icecol_Cpi                  = ice%Cpi(              vi,:)
      icecol_dzeta_dx             = ice%dzeta_dx_ak(      vi,:)
      icecol_dzeta_dy             = ice%dzeta_dy_ak(      vi,:)
      icecol_dzeta_dz             = ice%dzeta_dz_ak(      vi,:)
      icecol_dzeta_dt             = ice%dzeta_dt_ak(      vi,:)
      icecol_Phi                  = ice%internal_heating( vi,:)

      ! Solve the 1-D heat equation in the vertical column;
      ! if the solution is unstable, try again with a smaller timestep

      found_stable_solution = .FALSE.
      it_dt = 0
      DO WHILE ((.NOT. found_stable_solution) .AND. it_dt < 10)

        it_dt = it_dt + 1
        dt_applied = dt * (0.5_dp**(REAL( it_dt-1,dp)))   ! When it_dt = 0, dt_applied = dt; when it_dt = 1, dt_applied = dt/2, etc.

        ! If dt_applied = dt, solve it once; if dt_applied = dt/2, solve it twice; etc.
        DO it_it_dt = 1, 2**(it_dt-1)

          ! Solve the heat equation in the vertical column
          IF (ice%mask_gl_gr( vi)) THEN
            ! Grounding line: use some combination of the solutions using Q_base_grnd and T_base_float as boundary conditions

            ! Fully grounded solution (default: immediately assigned to final solution)
            CALL solve_1D_heat_equation( mesh, icecol_Ti, icecol_u, icecol_v, icecol_w, &
              icecol_u_times_dTdxp_upwind, icecol_v_times_dTdyp_upwind, T_surf_annual( vi), &
              icecol_Ti_pmp, icecol_Ki, icecol_Cpi, icecol_dzeta_dx, icecol_dzeta_dy, icecol_dzeta_dz, icecol_dzeta_dt, &
              icecol_Phi, dt_applied, icecol_Ti_tplusdt, Q_base_grnd = Q_base_grnd( vi))

            IF (C%choice_GL_temperature_BC == 'subgrid' .OR. C%choice_GL_temperature_BC == 'pmp') THEN
              ! Fully floating solution: assumes base is at the pressure melting point
              CALL solve_1D_heat_equation( mesh, icecol_Ti, icecol_u, icecol_v, icecol_w, &
                icecol_u_times_dTdxp_upwind, icecol_v_times_dTdyp_upwind, T_surf_annual( vi), &
                icecol_Ti_pmp, icecol_Ki, icecol_Cpi, icecol_dzeta_dx, icecol_dzeta_dy, icecol_dzeta_dz, icecol_dzeta_dt, &
                icecol_Phi, dt_applied, icecol_Ti_tplusdt_gl_fl, T_base_float = T_base_float( vi))
            ELSE
              ! Safety: just copy the profile from the default, fully grounded solution
              icecol_Ti_tplusdt_gl_fl = icecol_Ti_tplusdt
            END IF

              ! Combine both solutions
              SELECT CASE (C%choice_GL_temperature_BC)
                CASE ('grounded')
                  ! Use solution that assumes a fully grounded ice column (already assigned)
                CASE ('subgrid')
                  ! Interpolate grounded and floating solutions based on grounded fraction
                  icecol_Ti_tplusdt = ice%fraction_gr( vi) * icecol_Ti_tplusdt + (1._dp - ice%fraction_gr( vi)) * icecol_Ti_tplusdt_gl_fl
                CASE ('pmp')
                  ! Use solution that assumes ice base is at pressure melting point
                  icecol_Ti_tplusdt = icecol_Ti_tplusdt_gl_fl
                CASE DEFAULT
                  CALL crash('unknown choice_GL_temperature_BC "' // TRIM( C%choice_GL_temperature_BC) // '"!')
              END SELECT

          ELSEIF (ice%mask_grounded_ice( vi)) THEN
            ! Grounded ice: use Q_base_grnd as boundary condition

            CALL solve_1D_heat_equation( mesh, icecol_Ti, icecol_u, icecol_v, icecol_w, &
              icecol_u_times_dTdxp_upwind, icecol_v_times_dTdyp_upwind, T_surf_annual( vi), &
              icecol_Ti_pmp, icecol_Ki, icecol_Cpi, icecol_dzeta_dx, icecol_dzeta_dy, icecol_dzeta_dz, icecol_dzeta_dt, &
              icecol_Phi, dt_applied, icecol_Ti_tplusdt, Q_base_grnd = Q_base_grnd( vi))

          ELSEIF (ice%mask_floating_ice( vi)) THEN
            ! Floating ice: use T_base_float as boundary condition

            CALL solve_1D_heat_equation( mesh, icecol_Ti, icecol_u, icecol_v, icecol_w, &
              icecol_u_times_dTdxp_upwind, icecol_v_times_dTdyp_upwind, T_surf_annual( vi), &
              icecol_Ti_pmp, icecol_Ki, icecol_Cpi, icecol_dzeta_dx, icecol_dzeta_dy, icecol_dzeta_dz, icecol_dzeta_dt, &
              icecol_Phi, dt_applied, icecol_Ti_tplusdt, T_base_float = T_base_float( vi))

          ELSE
            CALL crash('Hi_eff > Hi_min_thermo, but mask_grounded_ice and mask_floating_ice are both .false.')
          END IF

          ! Update temperature solution for next semi-time-step
          icecol_Ti = icecol_Ti_tplusdt

        END DO ! DO it_it_dt = 1, 2**(it_dt-1)

        ! Check if we found a stable solution
        found_stable_solution = .TRUE.
        DO k = 1, mesh%nz
          IF (icecol_Ti( k) /= icecol_Ti( k)) found_stable_solution = .FALSE.   ! If we found NaN, the solution is unstable
          IF (icecol_Ti( k) < 180._dp)        found_stable_solution = .FALSE.   ! If we found temperatures below 180 K, the solution is unstable
          IF (icecol_Ti( k) > T0)             found_stable_solution = .FALSE.   ! If we found temperatures above freezing point, the solution is unstable
        END DO

      END DO ! DO WHILE ((.NOT. found_stable_solution) .AND. it_dt < 10)

      ! If the solution is still unstable, set the vertical temperature profile
      ! here to the Robin solution. Not the most accurate, but at least it will
      ! keep the model running. If too many grid cells need this, crash the model.
      IF (.NOT. found_stable_solution) THEN
        is_unstable( vi) = 1
        n_unstable = n_unstable + 1
      END IF

      ! Copy temperature solution
      Ti_tplusdt( vi,:) = icecol_Ti

    END DO

    ! Cope with instability
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, n_unstable, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    IF (n_unstable < CEILING( REAL( mesh%nV) / 100._dp)) THEN
      ! Instability is limited to an acceptably small number (< 1%) of grid cells;
      ! replace the temperature profile in those cells with the Robin solution

      DO vi = mesh%vi1, mesh%vi2
        IF (is_unstable( vi) == 1) THEN
          CALL replace_Ti_with_robin_solution( mesh, ice, climate, SMB, Ti_tplusdt, vi)
        END IF
      END DO

    ELSE
      ! An unacceptably large number of grid cells was unstable; throw an error.

      ! CALL save_variable_as_netcdf_dp_1D(  C%output_dir, ice%Hi                  , 'Hi'                  )
      ! CALL save_variable_as_netcdf_dp_1D(  C%output_dir, ice%Hi_eff              , 'Hi_eff'              )
      ! CALL save_variable_as_netcdf_dp_1D(  C%output_dir, ice%Hb                  , 'Hb'                  )
      ! CALL save_variable_as_netcdf_dp_1D(  C%output_dir, ice%SL                  , 'SL'                  )
      ! CALL save_variable_as_netcdf_dp_2D(  C%output_dir, ice%dzeta_dx_ak         , 'dzeta_dx_ak'         )
      ! CALL save_variable_as_netcdf_dp_2D(  C%output_dir, ice%dzeta_dy_ak         , 'dzeta_dy_ak'         )
      ! CALL save_variable_as_netcdf_dp_2D(  C%output_dir, ice%dzeta_dz_ak         , 'dzeta_dz_ak'         )
      ! CALL save_variable_as_netcdf_dp_2D(  C%output_dir, ice%dzeta_dt_ak         , 'dzeta_dt_ak'         )
      ! CALL save_variable_as_netcdf_dp_2D(  C%output_dir, ice%u_3D                , 'u_3D'                )
      ! CALL save_variable_as_netcdf_dp_2D(  C%output_dir, ice%v_3D                , 'v_3D'                )
      ! CALL save_variable_as_netcdf_dp_2D(  C%output_dir, ice%w_3D                , 'w_3D'                )
      ! CALL save_variable_as_netcdf_dp_1D(  C%output_dir, ice%frictional_heating  , 'frictional_heating'  )
      ! CALL save_variable_as_netcdf_dp_2D(  C%output_dir, ice%internal_heating    , 'internal_heating'    )
      ! CALL save_variable_as_netcdf_dp_1D(  C%output_dir, T_surf_annual           , 'T_surf_annual'       )
      ! CALL save_variable_as_netcdf_dp_1D(  C%output_dir, Q_base_grnd             , 'Q_base_grnd'         )
      ! CALL save_variable_as_netcdf_dp_1D(  C%output_dir, T_base_float            , 'T_base_float'        )
      ! CALL save_variable_as_netcdf_dp_2D(  C%output_dir, ice%Ti                  , 'Ti'                  )
      ! CALL save_variable_as_netcdf_dp_2D(  C%output_dir, Ti_tplusdt              , 'Ti_tplusdt'          )
      ! CALL save_variable_as_netcdf_int_1D( C%output_dir, is_unstable             , 'is_unstable'         )
      CALL crash('heat equation solver unstable for more than 1% of vertices!')

    END IF

    ! Update modelled ice temperature for time stepping
    ice%Ti_next = Ti_tplusdt

    ! Clean up after yourself
    DEALLOCATE( u_times_dTdxp_upwind)
    DEALLOCATE( v_times_dTdyp_upwind)
    DEALLOCATE( T_surf_annual       )
    DEALLOCATE( Q_base_grnd         )
    DEALLOCATE( T_base_float        )
    DEALLOCATE( Ti_tplusdt          )
    DEALLOCATE( is_unstable         )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE solve_3D_heat_equation

  SUBROUTINE solve_1D_heat_equation( mesh, Ti, u, v, w, u_times_dTdxp_upwind, v_times_dTdyp_upwind, T_surf, &
    Ti_pmp, Ki, Cpi, dzeta_dx, dzeta_dy, dzeta_dz, dzeta_dt, Phi, dt, Ti_tplusdt, Q_base_grnd, T_base_float)
    ! Solve the heat equation in the vertical column, i.e. using implicit discretisation of dT/dz, but explicit for dT/dx, dT/dy
    !
    ! The general heat equation (i.e. conservation of energy) inside the ice reads:
    !
    !   dT/dt = k / (rho cp) grad^2 T - u dT/dx - v dT/dy - w dT/dz + Phi / (rho cp)
    !
    ! With the following quantities:
    !
    !   T     - ice temperature
    !   k     - thermal conductivity of the ice
    !   rho   - density of the ice
    !   cp    - specific heat of the ice
    !   u,v,w - velocity of the ice
    !   Phi   - internal heating of the ice (due to strain rates)
    !
    ! By neglecting horizontal diffusion of heat (which can be done because of the flat
    ! geometry of the ice sheet), the grad^2 operator reduces to d2/dz2:
    !
    !   dT/dt = k / (rho cp) d2T/dz2 - u dT/dx - v dT/dy - w dT/dz + Phi / (rho cp)
    !
    ! Transforming into [xp,yp,zeta]-coordinates (xp = x, yp = y, zeta = (Hs(x,y) - z) / Hi(x,y)) yields:
    !
    !   dT/dt + dT/dzeta dzeta/dt = k / (rho cp) d2T/dzeta2 (dzeta/dz)^2 ...
    !                                - u (dT/dxp + dT/dzeta dzeta/dx) ...
    !                                - v (dT/dyp + dT/dzeta dzeta/dy) ...
    !                                - w (         dT/dzeta dzeta/dz) ...
    !                                + Phi / (rho cp)
    !
    ! The horizontal temperature gradients dT/dxp, dT/dyp are discretised explicitly in time,
    ! whereas the vertical gradients dT/dzeta, d2T/dzeta2 are discretised implicitly, yielding:
    !
    !   (T( t+dt) - T( t)) / dt + dzeta/dt d/dzeta( T( t+dt)) = ...
    !       k / (rho cp) (dzeta/dz)^2 d2/dzeta2( T( t+dt)) ...
    !     - u (dT( t)/dxp + dzeta/dx  d /dzeta ( T( t+dt))) ...
    !     - v (dT( t)/dyp + dzeta/dy  d /dzeta ( T( t+dt))) ...
    !     - w (             dzeta/dz  d /dzeta ( T( t+dt))) ...
    !     + Phi / (rho cp)
    !
    ! Moving all terms involving the unknown T( t+dt) to the left-hand side,
    ! and all other terms to the right-hand side, yields:
    !
    !   [ 1 / dt + dzeta/dt d/dzeta ...
    !     - k / (rho cp) (dzeta/dz)^2 d2/dzeta2 ...
    !     + u dzeta/dx d/dzeta ...
    !     + v dzeta/dy d/dzeta ...
    !     + w dzeta/dz d/dzeta ] T( t+dt) = ...
    !      T( t) / dt - u dT( t)/dxp - v dT( t)/dyp + Phi / (rho cp)
    !
    ! This can be further rearranged to read:
    !
    !   [ 1/dt + (dzeta/dt + u dzeta/dx + v dzeta/dy + w dzeta/dz) d/dzeta - k / (rho cp) (dzeta/dz)^2 d2/dzeta2 ] T( t+dt) = ...
    !      T( t) / dt - u dT( t)/dxp - v dT( t)/dyp + Phi / (rho cp)
    !
    ! Discretising this expression in space, the d/dzeta and d2/dzeta2 operators on the
    ! left-hand side become matrices, while all other terms become vectors. The equation
    ! thus describes a system of linear equations, with T( t+dt) the solution.
    !
    ! In order to solve this, boundary conditions must be applied at the surface and
    ! base of the ice. At the surface, the ice temperature is assumed to equal the
    ! annual mean surface temperature. At the base of grounded ice, the temperature
    ! gradient dT/dz must follow the heat flux at the base: dT/dz = -Q_base / k, and
    ! the temperature must also not exceed the pressure melting point temperature.
    ! At the base of floating ice, the temperature is assumed to equal the temperature
    ! of the ocean, limited by the pressure melting point temperature.
    !
    ! NOTE: with the current vertical discretisation scheme, d/dzeta and d2/dzeta2 are
    !       tridiagonal matrices, so that the entire matrix is tridiagonal. The LAPACK
    !       solver depends on this fact; changing the vertical discretisation to a
    !       higher-order scheme will require a different matrix solver.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh                 ! Contains the d/dzeta, d2/dzeta2 operators
    REAL(dp), DIMENSION( mesh%nz),          INTENT(IN)    :: Ti                   ! Vertical profile of ice temperature at time t                     [K]
    REAL(dp), DIMENSION( mesh%nz),          INTENT(IN)    :: u                    !   "         "    "  horizontal ice velocity in the x-direction    [m yr^-1]
    REAL(dp), DIMENSION( mesh%nz),          INTENT(IN)    :: v                    !   "         "    "  horizontal ice velocity in the y-direction    [m yr^-1]
    REAL(dp), DIMENSION( mesh%nz),          INTENT(IN)    :: w                    !   "         "    "  vertical   ice velocity                       [m yr^-1]
    REAL(dp), DIMENSION( mesh%nz),          INTENT(IN)    :: u_times_dTdxp_upwind !   "         "    "  u * dT/dxp in the upwind direction            [K yr^-1]
    REAL(dp), DIMENSION( mesh%nz),          INTENT(IN)    :: v_times_dTdyp_upwind !   "         "    "  u * dT/dxp in the upwind direction            [K yr^-1]
    REAL(dp),                               INTENT(IN)    :: T_surf               ! Annual mean surface temperature                                   [K]
    REAL(dp), DIMENSION( mesh%nz),          INTENT(IN)    :: Ti_pmp               !   "         "    "  pressure melting point temperature            [K]
    REAL(dp), DIMENSION( mesh%nz),          INTENT(IN)    :: Ki                   !   "         "    "  thermal conductivity of ice                   [J yr^-1 m^-1 K^-1]
    REAL(dp), DIMENSION( mesh%nz),          INTENT(IN)    :: Cpi                  !   "         "    "  specific heat of ice                          [J kg^-1 K^-1]
    REAL(dp), DIMENSION( mesh%nz),          INTENT(IN)    :: dzeta_dx             !   "         "    "  x-derivative of the scaled coordinate zeta    [m^-1]
    REAL(dp), DIMENSION( mesh%nz),          INTENT(IN)    :: dzeta_dy             !   "         "    "  y-derivative of the scaled coordinate zeta    [m^-1]
    REAL(dp), DIMENSION( mesh%nz),          INTENT(IN)    :: dzeta_dz             !   "         "    "  z-derivative of the scaled coordinate zeta    [m^-1]
    REAL(dp), DIMENSION( mesh%nz),          INTENT(IN)    :: dzeta_dt             !   "         "    "  time-derivative of the scaled coordinate zeta [yr^-1]
    REAL(dp), DIMENSION( mesh%nz),          INTENT(IN)    :: Phi                  !   "         "    "  internal heat production                      [J kg^1 yr^1]
    REAL(dp),                               INTENT(IN)    :: dt                   ! Time step                                                         [yr]
    REAL(dp), DIMENSION( mesh%nz),          INTENT(OUT)   :: Ti_tplusdt           ! Vertical profile of ice temperature at time t + dt                [K]
    REAL(dp),                     OPTIONAL, INTENT(IN)    :: Q_base_grnd   ! Heat flux at the base of grounded ice                             [J m^-2 yr^-1]
    REAL(dp),                     OPTIONAL, INTENT(IN)    :: T_base_float  ! Heat flux at the ice base                                         [J m^-2 yr^-1]

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'solve_1D_heat_equation'
    REAL(dp), DIMENSION( mesh%nz-1)                       :: AA_ldiag           ! Lower diagonal of A
    REAL(dp), DIMENSION( mesh%nz  )                       :: AA_diag            !       Diagonal of A
    REAL(dp), DIMENSION( mesh%nz-1)                       :: AA_udiag           ! Upper diagonal of A
    REAL(dp), DIMENSION( mesh%nz)                         :: bb                 ! Right-hand side b
    INTEGER                                               :: k
    REAL(dp)                                              :: c_ddzeta, c_d2dzeta2

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Initialise
    AA_ldiag = 0._dp
    AA_diag  = 0._dp
    AA_udiag = 0._dp
    bb       = 0._dp

    ! Ice surface boundary conditions: T = MIN( T_surf, T0)
    AA_diag(  1) = 1._dp
    bb(       1) = MIN( T_surf, T0)

    ! Ice base: either dT( t+dt)/dz = -Q_base / k, T( t+dt) <= T_pmp  for grounded ice, or:
    !                   T           = MIN( T_base_float, T_pmp)       for floating ice

    IF (PRESENT( Q_base_grnd)) THEN
      ! For grounded ice, let dT/dz = -Q_base / k

      ! Safety
      IF (PRESENT( T_base_float)) CALL crash('must provide either Q_base_grnd or T_base_float, but not both!')

      AA_diag(  mesh%nz) = 1._dp
      bb(       mesh%nz) = MIN( Ti_pmp( mesh%nz), Ti( mesh%nz-1) - (mesh%zeta( mesh%nz) - mesh%zeta( mesh%nz-1)) * Q_base_grnd / (dzeta_dz( mesh%nz) * Ki( mesh%nz)))

    ELSEIF (PRESENT( T_base_float)) THEN
      ! For floating ice, let T = MIN( T_base_float, T_pmp)

      ! Safety
      IF (PRESENT( Q_base_grnd)) CALL crash('must provide either Q_base_grnd or T_base_float, but not both!')

      AA_diag(  mesh%nz) = 1._dp
      bb(       mesh%nz) = MIN( T_base_float, Ti_pmp( mesh%nz))

    ELSE
      CALL crash('must provide either Q_base_grnd or T_base_float, but not both!')
    END IF ! IF (PRESENT( Q_base_grnd)) THEN

    ! Ice column
    DO k = 2, mesh%nz-1

      ! Calculate matrix coefficients
      c_ddzeta   = dzeta_dt( k) + (u( k) * dzeta_dx( k)) + (v( k) * dzeta_dy( k)) + (w( k) * dzeta_dz( k))
      c_d2dzeta2 = -Ki( k) / (ice_density * Cpi( k)) * dzeta_dz( k)**2

      AA_ldiag( k-1) =              (c_ddzeta * mesh%M_ddzeta_k_k_ldiag( k-1)) + (c_d2dzeta2 * mesh%M_d2dzeta2_k_k_ldiag( k-1))
      AA_diag(  k  ) = 1._dp / dt + (c_ddzeta * mesh%M_ddzeta_k_k_diag(  k  )) + (c_d2dzeta2 * mesh%M_d2dzeta2_k_k_diag(  k  ))
      AA_udiag( k  ) =              (c_ddzeta * mesh%M_ddzeta_k_k_udiag( k  )) + (c_d2dzeta2 * mesh%M_d2dzeta2_k_k_udiag( k  ))

      ! Calculate right-hand side
      bb(       k) = Ti( k) / dt - u_times_dTdxp_upwind( k) - v_times_dTdyp_upwind( k) + Phi( k) / (ice_density * Cpi( k))

    END DO ! DO k = 1, mesh%nz

    ! Solve the tridiagonal matrix equation representing the heat equation for this grid cell
    call solve_tridiagonal_matrix_equation( C%nz, AA_ldiag, AA_diag, AA_udiag, bb, Ti_tplusdt)

    ! Make sure ice temperature doesn't exceed pressure melting point
    DO k = 1, mesh%nz
      Ti_tplusdt( k) = MIN( Ti_tplusdt( k), Ti_pmp( k))
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE solve_1D_heat_equation

  SUBROUTINE write_to_restart_file_thermo_3D_heat_equation( mesh, ice, time)
    ! Write to the restart NetCDF file for the thermodynamics

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_model),                INTENT(IN)              :: ice
    REAL(dp),                            INTENT(IN)              :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'write_to_restart_file_thermo_3D_heat_equation'
    INTEGER                                                      :: ncid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If no NetCDF output should be created, do nothing
    IF (.NOT. C%do_create_netcdf_output) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Print to terminal
    IF (par%primary) WRITE(0,'(A)') '   Writing to thermodynamics restart file "' // &
      colour_string( TRIM( ice%thermo_restart_filename), 'light blue') // '"...'

    ! Open the NetCDF file
    CALL open_existing_netcdf_file_for_writing( ice%thermo_restart_filename, ncid)

    ! Write the time to the file
    CALL write_time_to_file( ice%thermo_restart_filename, ncid, time)

    ! Write the velocity fields to the file
    CALL write_to_field_multopt_mesh_dp_3D( mesh, ice%thermo_restart_filename, ncid, 'Ti', ice%Ti)

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_restart_file_thermo_3D_heat_equation

  SUBROUTINE create_restart_file_thermo_3D_heat_equation( mesh, ice)
    ! Create a restart NetCDF file for the thermodynamics
    ! Includes generation of the procedural filename (e.g. "restart_thermo_00001.nc")

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_model),                INTENT(INOUT)           :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'create_restart_file_thermo_3D_heat_equation'
    CHARACTER(LEN=256)                                           :: filename_base
    INTEGER                                                      :: ncid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If no NetCDF output should be created, do nothing
    IF (.NOT. C%do_create_netcdf_output) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Set the filename
    filename_base = TRIM( C%output_dir) // 'restart_thermodynamics'
    CALL generate_filename_XXXXXdotnc( filename_base, ice%thermo_restart_filename)

    ! Print to terminal
    IF (par%primary) WRITE(0,'(A)') '   Creating thermodynamics restart file "' // &
      colour_string( TRIM( ice%thermo_restart_filename), 'light blue') // '"...'

    ! Create the NetCDF file
    CALL create_new_netcdf_file_for_writing( ice%thermo_restart_filename, ncid)

    ! Set up the mesh in the file
    CALL setup_mesh_in_netcdf_file( ice%thermo_restart_filename, ncid, mesh)

    ! Add a time dimension to the file
    CALL add_time_dimension_to_file( ice%thermo_restart_filename, ncid)

    ! Add a zeta dimension to the file
    CALL add_zeta_dimension_to_file( ice%thermo_restart_filename, ncid, mesh%zeta)

    ! Add the temperature field to the file
    CALL add_field_mesh_dp_3D( ice%thermo_restart_filename, ncid, 'Ti', long_name = 'Englacial temperature', units = 'K')

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_restart_file_thermo_3D_heat_equation

END MODULE thermodynamics_3D_heat_equation
