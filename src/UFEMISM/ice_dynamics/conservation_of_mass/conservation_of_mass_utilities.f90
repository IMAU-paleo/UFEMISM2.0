module conservation_of_mass_utilities

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use CSR_sparse_matrix_type, only: type_sparse_matrix_CSR_dp
  use CSR_matrix_basics, only: allocate_matrix_CSR_dist, add_entry_CSR_dist, finalise_matrix_CSR_dist
  use map_velocities_to_c_grid, only: map_velocities_from_b_to_c_2D
  use mpi_distributed_memory, only: gather_to_all
  use mpi_f08, only: MPI_ALLREDUCE, MPI_IN_PLACE, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD

  implicit none

  private

  public :: calc_ice_flux_divergence_matrix_upwind, apply_mask_noice_direct, calc_flux_limited_timestep, &
    calc_n_interior_neighbours

contains

  subroutine calc_ice_flux_divergence_matrix_upwind( mesh, u_vav_b, v_vav_b, fraction_margin, M_divQ)
    !< Calculate the ice flux divergence matrix M_divQ using an upwind scheme

    ! The vertically averaged ice flux divergence represents the net ice volume (which,
    ! assuming constant density, is proportional to the ice mass) entering each Voronoi
    ! cell per unit time. This is found by calculating the ice fluxes through each
    ! shared Voronoi cell boundary, using an upwind scheme: if ice flows from vertex vi
    ! to vertex vj, the flux is found by multiplying the velocity at their shared
    ! boundary u_c with the ice thickness at vi (and, of course, the length L_c of the
    ! shared boundary). If instead it flows from vj to vi, u_c is multiplied with the
    ! ice thickness at vj.

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    real(dp), dimension(mesh%ti1:mesh%ti1), intent(in   ) :: u_vav_b
    real(dp), dimension(mesh%ti1:mesh%ti1), intent(in   ) :: v_vav_b
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   ) :: fraction_margin
    type(type_sparse_matrix_CSR_dp),        intent(  out) :: M_divQ

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'calc_ice_flux_divergence_matrix_upwind'
    real(dp), dimension(mesh%ei1:mesh%ei2) :: u_vav_c, v_vav_c
    real(dp), dimension(mesh%nE)           :: u_vav_c_tot, v_vav_c_tot
    real(dp), dimension(mesh%nV)           :: fraction_margin_tot
    integer                                :: ncols, ncols_loc, nrows, nrows_loc, nnz_est_proc
    integer                                :: vi, ci, ei, vj
    real(dp)                               :: A_i, L_c
    real(dp)                               :: u_perp
    real(dp), dimension(0:mesh%nC_mem)     :: cM_divQ

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate vertically averaged ice velocities on the edges
    call map_velocities_from_b_to_c_2D( mesh, u_vav_b, v_vav_b, u_vav_c, v_vav_c)
    call gather_to_all( u_vav_c, u_vav_c_tot)
    call gather_to_all( v_vav_c, v_vav_c_tot)
    call gather_to_all( fraction_margin, fraction_margin_tot)

    ! == Initialise the matrix using the native UFEMISM CSR-matrix format
    ! ===================================================================

    ! Matrix size
    ncols           = mesh%nV      ! from
    ncols_loc       = mesh%nV_loc
    nrows           = mesh%nV      ! to
    nrows_loc       = mesh%nV_loc
    nnz_est_proc    = mesh%nV_loc + SUM( mesh%nC( mesh%vi1:mesh%vi2))

    call allocate_matrix_CSR_dist( M_divQ, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc, &
      pai_x = mesh%pai_V, pai_y = mesh%pai_V)

    ! == Calculate coefficients
    ! =========================

    do vi = mesh%vi1, mesh%vi2

      ! Initialise
      cM_divQ = 0._dp

      ! Loop over all connections of vertex vi
      do ci = 1, mesh%nC( vi)

        ! Connection ci from vertex vi leads through edge ei to vertex vj
        ei = mesh%VE( vi,ci)
        vj = mesh%C(  vi,ci)

        ! The Voronoi cell of vertex vi has area A_i
        A_i = mesh%A( vi)

        ! The shared Voronoi cell boundary section between the Voronoi cells
        ! of vertices vi and vj has length L_c
        L_c = mesh%Cw( vi,ci)

        ! Calculate vertically averaged ice velocity component perpendicular to this shared Voronoi cell boundary section
        u_perp = u_vav_c_tot( ei) * mesh%D_x( vi, ci)/mesh%D( vi, ci) + v_vav_c_tot( ei) * mesh%D_y( vi, ci)/mesh%D( vi, ci)

        ! Calculate matrix coefficients
        ! =============================

        ! u_perp > 0: flow is exiting this vertex into vertex vj
        if (fraction_margin_tot( vi) >= 1._dp) then
          cM_divQ( 0) = cM_divQ( 0) + L_c * max( 0._dp, u_perp) / A_i
        else
          ! if this vertex is not completely covering its assigned area, then don't let ice out of it yet.
        end if

        ! u_perp < 0: flow is entering this vertex from vertex vj
        if (fraction_margin_tot( vj) >= 1._dp) then
          cM_divQ( ci) = L_c * MIN( 0._dp, u_perp) / A_i
        else
          ! if that vertex is not completely covering its assigned area, then don't let ice out of it yet.
        end if

      end do ! do ci = 1, mesh%nC( vi)

      ! Add coefficients to matrix
      call add_entry_CSR_dist( M_divQ, vi, vi, cM_divQ( 0))
      do ci = 1, mesh%nC( vi)
        vj = mesh%C(  vi,ci)
        call add_entry_CSR_dist( M_divQ, vi, vj, cM_divQ( ci))
      end do ! do ci = 1, mesh%nC( vi)

    end do ! do vi = mesh%vi1, mesh%vi2

    call finalise_matrix_CSR_dist( M_divQ)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_ice_flux_divergence_matrix_upwind

  subroutine apply_mask_noice_direct( mesh, mask_noice, Hi)
    !< Enforce Hi = 0 where told to do so

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    logical,  dimension(mesh%vi1:mesh%vi2), intent(in   ) :: mask_noice
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(inout) :: Hi

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'apply_mask_noice_direct'
    integer                        :: vi

    ! Add routine to path
    call init_routine( routine_name)

    do vi = mesh%vi1, mesh%vi2
      if (mask_noice( vi)) Hi( vi) = 0._dp
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine apply_mask_noice_direct

  subroutine calc_flux_limited_timestep( mesh, Hi, dHi_dt, dt_max)
    !< Calculate the largest time step that does not result in more
    !< ice flowing out of a cell than is contained within it.

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   ) :: Hi
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   ) :: dHi_dt
    real(dp),                               intent(  out) :: dt_max

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'calc_flux_limited_timestep'
    integer                                :: vi
    real(dp), dimension(mesh%vi1:mesh%vi2) :: dt_lim
    integer                                :: ierr

    ! Add routine to path
    call init_routine( routine_name)

    ! Initialise
    dt_lim = C%dt_ice_max

    ! Loop over each mesh vertex within this process
    do vi = mesh%vi1, mesh%vi2
      ! if there is [non-negligible] ice, and there is mass loss
      if (Hi( vi) > C%Hi_min .and. dHi_dt( vi) < 0._dp) then

        ! Compute time step limit (in yr) based on
        ! available ice thickness and flux divergence
        dt_lim( vi) = Hi( vi) / max( dHi_dt( vi), 1E-9_dp)

      end if
    end do

    ! Get most strict time step limit for this process
    dt_max = minval( dt_lim)

    ! Get most strict time step limit among all processes
    call MPI_ALLREDUCE( MPI_IN_PLACE, dt_max, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)

    ! Limit to minimum ice model time step
    dt_max = max( C%dt_ice_min, dt_max)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_flux_limited_timestep

  subroutine calc_n_interior_neighbours( mesh, mask_noice, n_interior_neighbours)

    ! In/output variables:
    type(type_mesh),                       intent(in   ) :: mesh
    logical, dimension(mesh%vi1:mesh%vi2), intent(in   ) :: mask_noice
    integer, dimension(mesh%vi1:mesh%vi2), intent(  out) :: n_interior_neighbours

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_flux_limited_timestep'
    logical,  dimension(mesh%nV)   :: mask_noice_tot
    integer                        :: vi, ci, vj

    ! Gather global data fields
    call gather_to_all( mask_noice, mask_noice_tot)

    do vi = mesh%vi1, mesh%vi2

      n_interior_neighbours( vi) = 0

      do ci = 1, mesh%nC( vi)
        vj = mesh%C( vi,ci)
        if (mesh%VBI( vj) == 0 .and. .not. mask_noice_tot( vj)) then
          n_interior_neighbours( vi) = n_interior_neighbours( vi) + 1
        end if
      end do

    end do

  end subroutine calc_n_interior_neighbours

end module conservation_of_mass_utilities