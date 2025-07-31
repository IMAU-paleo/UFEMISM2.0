MODULE thermodynamics_utilities

  ! General physical terms required for the thermodynamics module

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  USE climate_model_types                                    , ONLY: type_climate_model
  USE SMB_model_types                                        , ONLY: type_SMB_model
  use mesh_disc_apply_operators, only: ddx_a_b_3D, ddy_a_b_3D
  use plane_geometry, only: cross2
  use mpi_distributed_memory, only: gather_to_all

  IMPLICIT NONE

CONTAINS

  SUBROUTINE calc_strain_heating( mesh, ice)
    ! Calculate internal heating due to strain rates
    !
    ! Bueler and Brown (2009), Eq. 8 (though they use Sigma instead of Phi):
    !
    !   Phi = 2 B(T*) D^(1/n + 1) = 2 A(T*)^(-1/n) D^(1/n + 1)
    !
    ! From the text just after their Eq. 6:
    !
    !   2D^2 = Dij Dij, so: D = SQRT( Dij Dij / 2)
    !
    ! Here, Dij is the strain rate tensor: Dij = 1/2 (dui/dxj + duj/dxi)
    !
    !         |          du/dx          1/2 (du/dy + dv/dx)     1/2 (du/dz + dw/dx) |
    !         |                                                                     |
    !   Dij = | 1/2 (du/dy + dv/dx)              dv/dy          1/2 (dv/dz + dw/dy) |
    !         |                                                                     |
    !         | 1/2 (du/dz + dw/dx)     1/2 (dv/dz + dw/dy)              dw/dz      |
    !
    ! So:
    !
    !   D = SQRT( 1/2 [ ...
    !                   (du/dx)^2     + 1/4 (du/dy + dv/dx)^2 + 1/4 (du/dz + dw/dx)^2 + ...
    !           1/4 (du/dy + dv/dx)^2 +         (dv/dy)^2     + 1/4 (dv/dz + dw/dy)^2 + ...
    !           1/4 (du/dz + dw/dx)^2 + 1/4 (dv/dx + dw/dy)^2 +         (dw/dz)^2 ])
    !
    !     = SQRT( 1/2 [ (du/dx)^2 + (dv/dy)^2 + (dv/dz)^2 + ...
    !                    1/2 (du/dy + dv/dx)^2 + 1/2 (du/dz + dw/dx)^2 + 1/2 (dv/dz + dw/dy)^2])

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_strain_heating'
    INTEGER                                            :: vi, k
    REAL(dp)                                           :: D

    ! Add routine to path
    CALL init_routine( routine_name)

    DO vi = mesh%vi1, mesh%vi2
    DO k = 1, C%nz

      ! No ice means no heating
      IF (ice%Hi( vi) < 0.1_dp) THEN
        ice%internal_heating( vi,k) = 0._dp
        CYCLE
      END IF

      ! Calculate the total strain rate D
      D = SQRT( 0.5_dp * (ice%du_dx_3D( vi,k)**2 + ice%dv_dy_3D( vi,k)**2 + ice%dw_dz_3D( vi,k)**2 + &
                0.5_dp * (ice%du_dy_3D( vi,k)    + ice%dv_dx_3D( vi,k))**2 + &
                0.5_dp * (ice%du_dz_3D( vi,k)    + ice%dw_dx_3D( vi,k))**2 + &
                0.5_dp * (ice%dv_dz_3D( vi,k)    + ice%dw_dy_3D( vi,k))**2 ))

      ! Calculate the strain heating rate Phi
      ice%internal_heating( vi,k) = 2._dp * ice%A_flow( vi,k)**(-1._dp / C%Glens_flow_law_exponent) * D**(1._dp / C%Glens_flow_law_exponent + 1._dp)

    END DO
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_strain_heating

  SUBROUTINE calc_frictional_heating( mesh, ice)
    ! Calculate frictional heating at the base due to sliding

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_frictional_heating'
    INTEGER                                            :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    ! No sliding means no friction
    IF (C%choice_sliding_law == 'no_sliding') THEN
      ice%frictional_heating = 0._dp
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Calculate frictional heating
    DO vi = mesh%vi1, mesh%vi2
      IF (ice%mask_grounded_ice( vi)) THEN
        ice%frictional_heating( vi) = ice%basal_friction_coefficient( vi) * ice%uabs_base( vi)
      ELSE
        ice%frictional_heating( vi) = 0._dp
      END IF
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_frictional_heating

  SUBROUTINE calc_heat_capacity( mesh, ice)
    ! Calculate the heat capacity of the ice

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_heat_capacity'
    INTEGER                                            :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    IF     (C%choice_ice_heat_capacity == 'uniform') THEN
      ! Apply a uniform value for the heat capacity

      ice%Cpi = C%uniform_ice_heat_capacity

    ELSEIF (C%choice_ice_heat_capacity == 'Pounder1965') THEN
      ! Calculate the heat capacity of ice according to Pounder: The Physics of Ice (1965)

      DO vi = mesh%vi1, mesh%vi2
        ice%Cpi( vi,:) = 2115.3_dp + 7.79293_dp * (ice%Ti( vi,:) - T0)
      END DO ! DO vi = mesh%vi1, mesh%vi2

    ELSE
      CALL crash('unknown choice_ice_heat_capacity "' // TRIM( C%choice_ice_heat_capacity) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_heat_capacity

  SUBROUTINE calc_thermal_conductivity( mesh, ice)
    ! Calculate the thermal conductivity of the ice

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_thermal_conductivity'
    INTEGER                                            :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    IF     (C%choice_ice_thermal_conductivity == 'uniform') THEN
      ! Apply a uniform value for the thermal conductivity

      ice%Ki = C%uniform_ice_thermal_conductivity

    ELSEIF (C%choice_ice_thermal_conductivity == 'Ritz1987') THEN
      ! Calculate the thermal conductivity of ice according to Ritz (1987)

      DO vi = mesh%vi1, mesh%vi2
        ice%Ki( vi,:) = 3.101E08_dp * EXP(-0.0057_dp * ice%Ti( vi,:))
      END DO ! DO vi = mesh%vi1, mesh%vi2

    ELSE
      CALL crash('unknown choice_ice_thermal_conductivity "' // TRIM( C%choice_ice_thermal_conductivity) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_thermal_conductivity

  SUBROUTINE calc_pressure_melting_point( mesh, ice)
    ! Calculate the pressure melting point of the ice according to Huybrechts (1992)

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_pressure_melting_point'
    INTEGER                                            :: vi,k

    ! Add routine to path
    CALL init_routine( routine_name)

    DO vi = mesh%vi1, mesh%vi2
    DO k = 1, mesh%nz
      ice%Ti_pmp( vi,k) = T0 - Clausius_Clapeyron_gradient * ice%Hi_eff( vi) * mesh%zeta( k)
    END DO
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_pressure_melting_point

  SUBROUTINE calc_homologous_temperature( mesh, ice)
    ! Calculate the pressure melting point of the ice according to Huybrechts (1992)

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_homologous_temperature'
    INTEGER                                            :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    DO vi = mesh%vi1, mesh%vi2
      ice%Ti_hom( vi) = MIN( 0._dp, ice%Ti( vi,C%nz) - ice%Ti_pmp( vi,C%nz))
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_homologous_temperature

  SUBROUTINE replace_Ti_with_robin_solution( mesh, ice, climate, SMB, Ti, vi)
    ! This function calculates for one horizontal grid point the temperature profiles
    ! using the surface temperature and the geothermal heat flux as boundary conditions.
    ! See Robin solution in: Cuffey & Paterson 2010, 4th ed, chapter 9, eq. (9.13) - (9.22).

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                                INTENT(IN)    :: mesh
    TYPE(type_ice_model),                           INTENT(INOUT) :: ice
    TYPE(type_climate_model),                       INTENT(IN)    :: climate
    TYPE(type_SMB_model),                           INTENT(IN)    :: SMB
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2,mesh%nz), INTENT(INOUT) :: Ti
    INTEGER,                                        INTENT(IN)    :: vi

    ! Local variables:
    INTEGER                                                       :: k
    REAL(dp)                                                      :: Ts
    REAL(dp)                                                      :: thermal_length_scale
    REAL(dp)                                                      :: distance_above_bed
    REAL(dp)                                                      :: erf1
    REAL(dp)                                                      :: erf2

    REAL(dp)                                                      :: thermal_conductivity_Robin
    REAL(dp)                                                      :: thermal_diffusivity_Robin
    REAL(dp)                                                      :: bottom_temperature_gradient_Robin

    REAL(dp), PARAMETER                                           :: kappa_0_ice_conductivity     = 9.828_dp          ! The linear constant in the thermal conductivity of ice [J m^-1 K^-1 s^-1], see equation (12.6), Ritz (1987), Cuffey & Paterson (2010, p. 400), Zwinger (2007)
    REAL(dp), PARAMETER                                           :: kappa_e_ice_conductivity     = 0.0057_dp         ! The exponent constant in the thermal conductivity of ice [K^-1], see equation (12.6), Ritz (1987), Cuffey & Paterson (2010, p. 400), Zwinger (2007)
    REAL(dp), PARAMETER                                           :: c_0_specific_heat            = 2127.5_dp         ! The constant in the specific heat capacity of ice [J kg^-1 K^-1], see equation (12.5), Zwinger (2007), Cuffey & Paterson (2010, p. 400)

    thermal_conductivity_Robin        = kappa_0_ice_conductivity * sec_per_year * EXP(-kappa_e_ice_conductivity * T0) ! Thermal conductivity            [J m^-1 K^-1 y^-1]
    thermal_diffusivity_Robin         = thermal_conductivity_Robin / (ice_density * c_0_specific_heat)                ! Thermal diffusivity             [m^2 y^-1]
    bottom_temperature_gradient_Robin = -1._dp * ice%geothermal_heat_flux( vi) / thermal_conductivity_Robin           ! Temperature gradient at bedrock

    Ts = MIN( T0, SUM( climate%T2m( vi,:)) / REAL( SIZE( climate%T2m,2),dp))

    IF (ice%Hi_eff( vi) > C%Hi_min_thermo) THEN
      ! This vertex has enough ice to have a noticeable temperature profile

      IF (ice%mask_grounded_ice( vi)) THEN
        ! This vertex has more than 1m of grounded ice

        IF (SMB%SMB( vi) > 0._dp) THEN
          ! The Robin solution can be used to estimate the subsurface temperature profile in an accumulation area

          thermal_length_scale = SQRT( 2._dp * thermal_diffusivity_Robin * ice%Hi_eff( vi) / SMB%SMB( vi))
          DO k = 1, C%nz
            distance_above_bed = (1._dp - mesh%zeta( k)) * ice%Hi_eff( vi)
            erf1 = erf( distance_above_bed / thermal_length_scale)
            erf2 = erf( ice%Hi_eff( vi) / thermal_length_scale)
            Ti( vi,k) = Ts + SQRT(pi) / 2._dp * thermal_length_scale * bottom_temperature_gradient_Robin * (erf1 - erf2)
          END DO

        ELSE ! IF (SMB%SMB( vi) > 0._dp) THEN

          ! Ablation area: use linear temperature profile from Ts to (offset below) T_pmp
          Ti( vi,:) = Ts + ((T0 - Clausius_Clapeyron_gradient * ice%Hi_eff( vi)) - Ts) * mesh%zeta

        END IF ! IF (SMB%SMB( vi) > 0._dp) THEN

      ELSEIF( ice%mask_floating_ice( vi)) THEN
        ! This vertex has more than 1m of floating ice
        ! Set a linear profile between T_surf and Ti_pmp_base

        Ti( vi,:) = Ts + mesh%zeta * (ice%Ti_pmp( vi,mesh%nz) - Ts)

      END IF ! IF (ice%mask_grounded_ice( vi)) THEN

    ELSE ! IF (ice%Hi_eff( vi) > C%Hi_min_thermo) THEN
      ! No (significant) ice present; set temperature to annual mean surface temperature

      Ti( vi,:) = Ts

    END IF ! IF (ice%Hi_eff( vi) > C%Hi_min_thermo) THEN

    ! Safety: limit temperatures to the pressure melting point
    DO k = 1, mesh%nz
      Ti( vi,k) = MIN( Ti( vi,k), ice%Ti_pmp( vi,k))
    END DO

  END SUBROUTINE replace_Ti_with_robin_solution

  SUBROUTINE calc_upwind_heat_flux_derivatives( mesh, ice, u_times_dTdxp_upwind, v_times_dTdyp_upwind)
    ! Calculate upwind heat flux derivatives at vertex vi, vertical layer k

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                                INTENT(IN)  :: mesh
    TYPE(type_ice_model),                           INTENT(IN)  :: ice
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2,mesh%nz), INTENT(OUT) :: u_times_dTdxp_upwind, v_times_dTdyp_upwind

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                               :: routine_name = 'calc_upwind_heat_flux_derivatives'
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE                     :: dTi_dxp_3D_b, dTi_dyp_3D_b
    INTEGER                                                     :: vi, k, vti, ti, n1, n2, n3, vib, vic, ti_upwind
    REAL(dp), DIMENSION(2)                                      :: u_upwind, ab, ac
    REAL(dp), DIMENSION(mesh%nTri,mesh%nz)                      :: u_3D_b_tot, v_3D_b_tot
    REAL(dp), DIMENSION(mesh%nTri,mesh%nz)                      :: dTi_dxp_3D_b_tot, dTi_dyp_3D_b_tot

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate memory
    ALLOCATE( dTi_dxp_3D_b( mesh%ti1:mesh%ti2,mesh%nz))
    ALLOCATE( dTi_dyp_3D_b( mesh%ti1:mesh%ti2,mesh%nz))

    ! Calculate dT/dxp, dT/dyp on the b-grid
    CALL ddx_a_b_3D( mesh, ice%Ti, dTi_dxp_3D_b)
    CALL ddy_a_b_3D( mesh, ice%Ti, dTi_dyp_3D_b)

    ! Gather full velocity fields
    CALL gather_to_all( ice%u_3D_b  , u_3D_b_tot      )
    CALL gather_to_all( ice%v_3D_b  , v_3D_b_tot      )
    CALL gather_to_all( dTi_dxp_3D_b, dTi_dxp_3D_b_tot)
    CALL gather_to_all( dTi_dyp_3D_b, dTi_dyp_3D_b_tot)

    DO vi = mesh%vi1, mesh%vi2

      ! Exception for the trivial case of no ice
      IF (ice%Hi( vi) < 1._dp) THEN
        u_times_dTdxp_upwind( vi,:) = 0._dp
        v_times_dTdyp_upwind( vi,:) = 0._dp
        CYCLE
      END IF

      ! The upwind velocity vector
      u_upwind = [-ice%u_vav( vi), -ice%v_vav( vi)]

      ! Find the upwind triangle
      ti_upwind = 0
      DO vti = 1, mesh%niTri( vi)

        ! Triangle ti is spanned counter-clockwise by vertices [vi,vib,vic]
        ti  = mesh%iTri( vi,vti)
        vib = 0
        vic = 0
        DO n1 = 1, 3
          n2 = n1 + 1
          IF (n2 == 4) n2 = 1
          n3 = n2 + 1
          IF (n3 == 4) n3 = 1

          IF (mesh%Tri( ti,n1) == vi) THEN
            vib = mesh%Tri( ti,n2)
            vic = mesh%Tri( ti,n3)
            EXIT
          END IF
        END DO

        ! Check if the upwind velocity vector points into this triangle
        ab = mesh%V( vib,:) - mesh%V( vi,:)
        ac = mesh%V( vic,:) - mesh%V( vi,:)

        IF (cross2( ab, u_upwind) >= 0._dp .AND. cross2( u_upwind, ac) >= 0._dp) THEN
          ti_upwind = ti
          EXIT
        END IF

      END DO ! DO iti = 1, mesh%niTri( vi)

      ! Safety
      IF (ti_upwind == 0) THEN
        CALL crash('could not find upwind triangle!')
      END IF

      ! Calculate u * dT/dx, v * dT/dy
      DO k = 1, C%nz
        u_times_dTdxp_upwind( vi,k) = u_3D_b_tot( ti_upwind,k) * dTi_dxp_3D_b_tot( ti_upwind,k)
        v_times_dTdyp_upwind( vi,k) = v_3D_b_tot( ti_upwind,k) * dTi_dyp_3D_b_tot( ti_upwind,k)
      END DO

    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_upwind_heat_flux_derivatives

END MODULE thermodynamics_utilities
