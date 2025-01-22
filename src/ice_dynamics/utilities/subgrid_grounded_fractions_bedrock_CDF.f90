module subgrid_grounded_fractions_bedrock_CDF
  !< Routines for calculating sub-grid grounded fractions using
  !< the sub-grid bedrock cumulative density functions

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use model_configuration, only: C
  use parameters, only: ice_density, seawater_density
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model
  use mesh_disc_apply_operators, only: map_a_b_2D
  use mpi_distributed_memory, only: gather_to_all

  implicit none

  private

  public :: calc_grounded_fractions_bedrock_CDF_a, calc_grounded_fractions_bedrock_CDF_b

contains

  subroutine calc_grounded_fractions_bedrock_CDF_a( mesh, ice, fraction_gr)
    !< Calculate the sub-grid grounded fractions of the vertices,
    !< using the sub-grid bedrock cumulative density functions (CDFs)

    ! In- and output variables
    type(type_mesh),                        intent(in   ) :: mesh
    type(type_ice_model),                   intent(in   ) :: ice
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(  out) :: fraction_gr

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_grounded_fractions_bedrock_CDF_a'
    integer                        :: vi, il, iu
    real(dp)                       :: Hb_float, wl, wu

    ! Add routine to path
    call init_routine( routine_name)

    do vi = mesh%vi1, mesh%vi2

      ! Compute the bedrock depth at which the current ice thickness and sea level
      ! will make this point afloat. Account for GIA here so we don't have to do it in
      ! the computation of the cumulative density function (CDF).

      Hb_float = ice%SL( vi) - ice%Hi( vi) * ice_density/seawater_density - ice%dHb( vi)

      ! Get the fraction of bedrock within vertex coverage that is below
      ! Hb_float as a linear interpolation of the numbers in the CDF.

      if     (Hb_float <= ice%bedrock_cdf( vi,1)) then
        ! All sub-grid points are above the floating bedrock elevation

        fraction_gr( vi) = 1._dp

      elseif (Hb_float >= ice%bedrock_cdf( vi,C%subgrid_bedrock_cdf_nbins)) then
        ! All sub-grid points are below the floating bedrock elevation

        fraction_gr( vi) = 0._dp

      else

        ! Find the 2 elements in the CDF surrounding Hb_float
        iu = 1
        do WHILE (ice%bedrock_cdf( vi,iu) < Hb_float)
          iu = iu+1
        end do
        il = iu-1

        ! Interpolate the two enveloping bedrock bins to find the grounded fraction
        wl = (ice%bedrock_cdf( vi,iu) - Hb_float) / (ice%bedrock_cdf( vi,iu) - ice%bedrock_cdf( vi,il))
        wu = 1._dp - wl
        fraction_gr( vi) = 1._dp - (real( il-1,dp) * wl + real( iu-1) * wu) / real( C%subgrid_bedrock_cdf_nbins-1,dp)

        ! Safety
        fraction_gr( vi) = min( 1._dp, max( 0._dp, fraction_gr( vi)))

      end if

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_grounded_fractions_bedrock_CDF_a

  subroutine calc_grounded_fractions_bedrock_CDF_b( mesh, ice, fraction_gr_b)
    !< Calculate the sub-grid grounded fractions of the triangles,
    !< using the sub-grid bedrock cumulative density functions (CDFs)

    ! In- and output variables
    type(type_mesh),                        intent(in   ) :: mesh
    type(type_ice_model),                   intent(in   ) :: ice
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(  out) :: fraction_gr_b

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'calc_grounded_fractions_bedrock_CDF_b'
    real(dp), dimension(mesh%nV)           :: TAF_tot
    real(dp), dimension(mesh%ti1:mesh%ti2) :: Hi_b, SL_b, dHb_b
    integer                                :: ti, via, vib, vic, il, iu
    real(dp)                               :: Hb_float, wl, wu

    ! Add routine to path
    call init_routine( routine_name)

    ! Map ice thickness, sea level, and bedrock deformation to the b-grid (triangles)
    call map_a_b_2D( mesh, ice%Hi , Hi_b )
    call map_a_b_2D( mesh, ice%SL , SL_b )
    call map_a_b_2D( mesh, ice%dHb, dHb_b)

    ! Gather global thickness above floatation
    call gather_to_all( ice%TAF, TAF_tot)

    do ti = mesh%ti1, mesh%ti2

      ! On the domain border, remapping issues make this answer unreliable
      ! (NOTE: only relevant when there's ice at the domain border, which in
      !        realistic experiments should never be the case; only happens
      !        in idealised geometries (e.g. MISMIP+))
      if (mesh%TriBI( ti) > 0) then
        ! if any of the three vertices spanning this triangle are grounded, treat it as grounded
        via = mesh%Tri( ti,1)
        vib = mesh%Tri( ti,2)
        vic = mesh%Tri( ti,3)
        if (TAF_tot( via) > 0._dp .OR. TAF_tot( vib) > 0._dp .OR. TAF_tot( vic) > 0._dp) then
          fraction_gr_b( ti) = 1._dp
        else
          fraction_gr_b( ti) = 0._dp
        end if
        cycle
      end if

      ! Compute the bedrock depth at which the current ice thickness and sea level
      ! will make this point afloat. Account for GIA here so we don't have to do it in
      ! the computation of the cumulative density function (CDF).

      Hb_float = SL_b( ti) - Hi_b( ti) * ice_density/seawater_density - dHb_b( ti)

      ! Get the fraction of bedrock within vertex coverage that is below
      ! Hb_float as a linear interpolation of the numbers in the CDF.

      if     (Hb_float <= ice%bedrock_cdf_b( ti,1)) then
        ! All sub-grid points are above the floating bedrock elevation

        fraction_gr_b( ti) = 1._dp

      elseif (Hb_float >= ice%bedrock_cdf_b( ti,C%subgrid_bedrock_cdf_nbins)) then
        ! All sub-grid points are below the floating bedrock elevation

        fraction_gr_b( ti) = 0._dp

      else

        ! Find the 2 elements in the CDF surrounding Hb_float
        iu = 1
        do WHILE (ice%bedrock_cdf_b( ti,iu) < Hb_float)
          iu = iu+1
        end do
        il = iu-1

        ! Interpolate the two enveloping bedrock bins to find the grounded fraction
        wl = (ice%bedrock_cdf_b( ti,iu) - Hb_float) / (ice%bedrock_cdf_b( ti,iu) - ice%bedrock_cdf_b( ti,il))
        wu = 1._dp - wl
        fraction_gr_b( ti) = 1._dp - (real( il-1,dp) * wl + real( iu-1) * wu) / real( C%subgrid_bedrock_cdf_nbins-1,dp)

        ! Safety
        fraction_gr_b( ti) = min( 1._dp, max( 0._dp, fraction_gr_b( ti)))

      end if

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_grounded_fractions_bedrock_CDF_b

end module subgrid_grounded_fractions_bedrock_CDF