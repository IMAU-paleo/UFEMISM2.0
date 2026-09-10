module momentum_balance_solver_SSADIVA

  ! Intermediate class containing data and functions that are
  ! shared between the SSA and DIVA momentum balance solvers

  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use crash_mod, only: crash
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use momentum_balance_solver_basic, only: atype_momentum_balance_solver
  use mesh_disc_apply_operators, only: map_a_b_2D, map_b_a_2D, ddx_a_b_2D, ddy_a_b_2D, ddx_b_a_2D, ddy_b_a_2D
  use remapping_main, only: map_from_mesh_to_mesh_with_reallocation_2D
  use reallocate_mod, only: reallocate_bounds, reallocate_clean
  use ice_geometry_model_data, only: atype_ice_geometry_model_data
  use mpi_f08, only: MPI_COMM_WORLD, MPI_ALLREDUCE, MPI_DOUBLE_PRECISION, MPI_IN_PLACE, MPI_SUM
  use parameters, only: ice_density, grav
  use checksum_mod, only: checksum
  use CSR_matrix_mod, only: type_CSR_matrix_dp
  use mesh_utilities, only: find_ti_copy_ISMIP_HOM_periodic, find_ti_copy_SSA_icestream_infinite
  use mpi_distributed_memory, only: gather_to_all
  use mpi_distributed_shared_memory, only: gather_dist_shared_to_all
  use petsc_basic, only: solve_matrix_equation_CSR_PETSc
  use mpi_f08, only: MPI_WIN
  use Arakawa_grid_mod, only: Arakawa_grid
  use fields_dimensions, only: third_dimension

  implicit none

  private

  public :: atype_momentum_balance_solver_SSADIVA

  type, abstract, extends(atype_momentum_balance_solver) :: atype_momentum_balance_solver_SSADIVA

      ! Solution
      real(dp), dimension(:), contiguous, pointer :: u_vav_b                       => null()   ! [m yr^-1] 2-D horizontal ice velocity
      real(dp), dimension(:), contiguous, pointer :: v_vav_b                       => null()
      type(MPI_WIN) :: wu_vav_b, wv_vav_b

      ! Intermediate data fields
      real(dp), dimension(:), contiguous, pointer :: du_dx_a                       => null()   ! [yr^-1] 2-D horizontal strain rates
      real(dp), dimension(:), contiguous, pointer :: du_dy_a                       => null()
      real(dp), dimension(:), contiguous, pointer :: dv_dx_a                       => null()
      real(dp), dimension(:), contiguous, pointer :: dv_dy_a                       => null()
      real(dp), dimension(:), allocatable :: eta_vav_a                   ! Effective viscosity
      real(dp), dimension(:), allocatable :: N_a                         ! Product term N = eta * H
      real(dp), dimension(:), allocatable :: N_b
      real(dp), dimension(:), allocatable :: dN_dx_b                     ! Gradients of N
      real(dp), dimension(:), allocatable :: dN_dy_b
      real(dp), dimension(:), allocatable :: basal_friction_coefficient_b! Basal friction coefficient (tau_b = u * beta_b)
      real(dp), dimension(:), allocatable :: tau_dx_b                    ! Driving stress
      real(dp), dimension(:), allocatable :: tau_dy_b
      real(dp), dimension(:), allocatable :: u_vav_b_prev                ! Velocity solution from previous viscosity iteration
      real(dp), dimension(:), allocatable :: v_vav_b_prev
      type(MPI_WIN) :: wdu_dx_a, wdu_dy_a, wdv_dx_a, wdv_dy_a

      ! Restart file
      character(len=256)                  :: restart_filename

    contains

      procedure, public :: allocate_shared_SSA_DIVA_variables
      procedure, public :: deallocate_shared_SSA_DIVA_variables
      procedure, public :: remap_shared_SSA_DIVA_variables

      procedure, public :: calc_driving_stress
      procedure, public :: calc_horizontal_strain_rates
      procedure, public :: relax_viscosity_iterations
      procedure, public :: apply_velocity_limits
      procedure, public :: calc_L2_norm_uv

      procedure, public :: solve_SSA_DIVA_linearised
      procedure, public :: assemble_SSA_DIVA_linearised_matrix_eq
      procedure, public :: calc_SSA_DIVA_stiffness_matrix_row_free
      procedure, public :: calc_SSA_DIVA_sans_stiffness_matrix_row_free
      procedure, public :: calc_SSA_DIVA_stiffness_matrix_row_BC

  end type atype_momentum_balance_solver_SSADIVA

  ! Interfaces for procedures defined in submodules
  interface

    module subroutine solve_SSA_DIVA_linearised( self, u_ii_term, n_Axb_its, BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b)
      class(atype_momentum_balance_solver_SSADIVA), intent(inout) :: self
      real(dp), dimension(self%mesh%ti1:self%mesh%ti2),   intent(in   ) :: u_ii_term             ! Term to add to the diagonal; either the basal friction coefficient in the SSA, or beta_eff in the DIVA
      integer,                                            intent(  out) :: n_Axb_its             ! Number of iterations used in the iterative solver
      integer,  dimension(self%mesh%ti1:self%mesh%ti2),   intent(in   ) :: BC_prescr_mask_b      ! Mask of triangles where velocity is prescribed
      real(dp), dimension(self%mesh%ti1:self%mesh%ti2),   intent(in   ) :: BC_prescr_u_b         ! Prescribed velocities in the x-direction
      real(dp), dimension(self%mesh%ti1:self%mesh%ti2),   intent(in   ) :: BC_prescr_v_b         ! Prescribed velocities in the y-direction
    end subroutine solve_SSA_DIVA_linearised

    module subroutine assemble_SSA_DIVA_linearised_matrix_eq( self, u_ii_term, &
      BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b, &
      A_CSR, bb, uv_buv)
      class(atype_momentum_balance_solver_SSADIVA),     intent(in   ) :: self
      real(dp), dimension(self%mesh%ti1:self%mesh%ti2), intent(in   ) :: u_ii_term             ! Term to add to the diagonal; either the basal friction coefficient in the SSA, or beta_eff in the DIVA
      integer,  dimension(self%mesh%ti1:self%mesh%ti2), intent(in   ) :: BC_prescr_mask_b      ! Mask of triangles where velocity is prescribed
      real(dp), dimension(self%mesh%ti1:self%mesh%ti2), intent(in   ) :: BC_prescr_u_b         ! Prescribed velocities in the x-direction
      real(dp), dimension(self%mesh%ti1:self%mesh%ti2), intent(in   ) :: BC_prescr_v_b         ! Prescribed velocities in the y-direction
      type(type_CSR_matrix_dp),                         intent(inout) :: A_CSR
      real(dp), dimension(:), allocatable,              intent(inout) :: bb
      real(dp), dimension(:), allocatable,              intent(inout) :: uv_buv

    end subroutine assemble_SSA_DIVA_linearised_matrix_eq

    module subroutine calc_SSA_DIVA_stiffness_matrix_row_free( self, u_ii_term, A_CSR, bb, row_tiuv)
      class(atype_momentum_balance_solver_SSADIVA), intent(in   ) :: self
      real(dp), dimension(self%mesh%ti1:self%mesh%ti2),   intent(in   ) :: u_ii_term             ! Term to add to the diagonal; either the basal friction coefficient in the SSA, or beta_eff in the DIVA
      type(type_CSR_matrix_dp),                           intent(inout) :: A_CSR
      real(dp), dimension(A_CSR%i1:A_CSR%i2),             intent(inout) :: bb
      integer,                                            intent(in   ) :: row_tiuv
    end subroutine calc_SSA_DIVA_stiffness_matrix_row_free

    module subroutine calc_SSA_DIVA_sans_stiffness_matrix_row_free( self, u_ii_term, A_CSR, bb, row_tiuv)
      class(atype_momentum_balance_solver_SSADIVA), intent(in   ) :: self
      real(dp), dimension(self%mesh%ti1:self%mesh%ti2),   intent(in   ) :: u_ii_term             ! Term to add to the diagonal; either the basal friction coefficient in the SSA, or beta_eff in the DIVA
      type(type_CSR_matrix_dp),                           intent(inout) :: A_CSR
      real(dp), dimension(A_CSR%i1:A_CSR%i2),             intent(inout) :: bb
      integer,                                            intent(in   ) :: row_tiuv
    end subroutine calc_SSA_DIVA_sans_stiffness_matrix_row_free

    module subroutine calc_SSA_DIVA_stiffness_matrix_row_BC( self, A_CSR, bb, row_tiuv, choice_BC_u, choice_BC_v)
      class(atype_momentum_balance_solver_SSADIVA), intent(in   ) :: self
      type(type_CSR_matrix_dp),                           intent(inout) :: A_CSR
      real(dp), dimension(A_CSR%i1:A_CSR%i2),             intent(inout) :: bb
      integer,                                            intent(in   ) :: row_tiuv
      character(len=*),                                   intent(in   ) :: choice_BC_u, choice_BC_v
    end subroutine calc_SSA_DIVA_stiffness_matrix_row_BC

  end interface

contains

  subroutine allocate_shared_SSA_DIVA_variables( self)

    ! In/output variables:
    class(atype_momentum_balance_solver_SSADIVA), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'allocate_shared_SSA_DIVA_variables'

    ! Add routine to path
    call init_routine( routine_name)

    ! Solution
    call self%create_field( self%u_vav_b, self%wu_vav_b, &
      self%mesh, Arakawa_grid%b(), &
      name      = 'u_vav_b', &
      long_name = 'Vertically averaged ice velocity in the x-direction on the triangles', &
      units     = 'm yr^-1', &
      remap_method = 'reallocate')
    call self%create_field( self%v_vav_b, self%wv_vav_b, &
      self%mesh, Arakawa_grid%b(), &
      name      = 'v_vav_b', &
      long_name = 'Vertically averaged ice velocity in the y-direction on the triangles', &
      units     = 'm yr^-1', &
      remap_method = 'reallocate')

    ! Intermediate data fields
    call self%create_field( self%du_dx_a, self%wdu_dx_a, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'du_dx_a', &
      long_name = 'Vertically averaged xx strain rate', &
      units     = 'yr^-1', &
      remap_method = 'reallocate')
    call self%create_field( self%du_dy_a, self%wdu_dy_a, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'du_dy_a', &
      long_name = 'Vertically averaged xy strain rate', &
      units     = 'yr^-1', &
      remap_method = 'reallocate')
    call self%create_field( self%dv_dx_a, self%wdv_dx_a, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'dv_dx_a', &
      long_name = 'Vertically averaged yx strain rate', &
      units     = 'yr^-1', &
      remap_method = 'reallocate')
    call self%create_field( self%dv_dy_a, self%wdv_dy_a, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'dv_dy_a', &
      long_name = 'Vertically averaged yy strain rate', &
      units     = 'yr^-1', &
      remap_method = 'reallocate')
    allocate( self%eta_vav_a(                    self%mesh%vi1:self%mesh%vi2), source = 0._dp)
    allocate( self%N_a(                          self%mesh%vi1:self%mesh%vi2), source = 0._dp)
    allocate( self%N_b(                          self%mesh%ti1:self%mesh%ti2), source = 0._dp)
    allocate( self%dN_dx_b(                      self%mesh%ti1:self%mesh%ti2), source = 0._dp)
    allocate( self%dN_dy_b(                      self%mesh%ti1:self%mesh%ti2), source = 0._dp)
    allocate( self%basal_friction_coefficient_b( self%mesh%ti1:self%mesh%ti2), source = 0._dp)
    allocate( self%tau_dx_b(                     self%mesh%ti1:self%mesh%ti2), source = 0._dp)
    allocate( self%tau_dy_b(                     self%mesh%ti1:self%mesh%ti2), source = 0._dp)
    allocate( self%u_vav_b_prev(                 self%mesh%nTri             ), source = 0._dp)
    allocate( self%v_vav_b_prev(                 self%mesh%nTri             ), source = 0._dp)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine allocate_shared_SSA_DIVA_variables

  subroutine deallocate_shared_SSA_DIVA_variables( self)

    ! In/output variables:
    class(atype_momentum_balance_solver_SSADIVA), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'deallocate_shared_SSA_DIVA_variables'

    ! Add routine to path
    call init_routine( routine_name)

    ! Solution
    nullify( self%u_vav_b)
    nullify( self%v_vav_b)

    ! Intermediate data fields
    nullify( self%du_dx_a)
    nullify( self%du_dy_a)
    nullify( self%dv_dx_a)
    nullify( self%dv_dy_a)
    deallocate( self%eta_vav_a)
    deallocate( self%N_a)
    deallocate( self%N_b)
    deallocate( self%dN_dx_b)
    deallocate( self%dN_dy_b)
    deallocate( self%basal_friction_coefficient_b)
    deallocate( self%tau_dx_b)
    deallocate( self%tau_dy_b)
    deallocate( self%u_vav_b_prev)
    deallocate( self%v_vav_b_prev)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine deallocate_shared_SSA_DIVA_variables

  subroutine remap_shared_SSA_DIVA_variables( self, mesh_old, mesh_new)

    ! In/output variables:
    class(atype_momentum_balance_solver_SSADIVA), intent(inout) :: self
    type(type_mesh),                                    intent(in   ) :: mesh_old
    type(type_mesh), target,                            intent(in   ) :: mesh_new

    ! Local variables:
    character(len=*), parameter         :: routine_name = 'remap_shared_SSA_DIVA_variables'
    real(dp), dimension(:), allocatable :: u_vav_a
    real(dp), dimension(:), allocatable :: v_vav_a

    ! Add routine to path
    call init_routine( routine_name)

    ! Remap the fields that are re-used during the viscosity iteration
    ! ================================================================

    ! allocate memory for velocities on the a-grid (vertices)
    allocate( u_vav_a ( mesh_old%vi1: mesh_old%vi2             ))
    allocate( v_vav_a ( mesh_old%vi1: mesh_old%vi2             ))

    ! Map data from the triangles of the old mesh to the vertices of the old mesh
    call map_b_a_2D( mesh_old, self%u_vav_b , u_vav_a )
    call map_b_a_2D( mesh_old, self%v_vav_b , v_vav_a )

    ! Remap data from the vertices of the old mesh to the vertices of the new mesh
    call map_from_mesh_to_mesh_with_reallocation_2D( mesh_old, mesh_new, C%output_dir, u_vav_a , '2nd_order_conservative')
    call map_from_mesh_to_mesh_with_reallocation_2D( mesh_old, mesh_new, C%output_dir, v_vav_a , '2nd_order_conservative')

    ! reallocate memory for the data on the triangles
    call self%remap_field( mesh_new, 'u_vav_b', self%u_vav_b)
    call self%remap_field( mesh_new, 'v_vav_b', self%v_vav_b)

    ! Map data from the vertices of the new mesh to the triangles of the new mesh
    call map_a_b_2D( mesh_new, u_vav_a , self%u_vav_b )
    call map_a_b_2D( mesh_new, v_vav_a , self%v_vav_b )

    ! reallocate everything else
    ! ==========================

   !call reallocate_bounds( DIVA%u_vav_b                     , mesh_new%ti1, mesh_new%ti2)
   !call reallocate_bounds( DIVA%v_vav_b                     , mesh_new%ti1, mesh_new%ti2)
    call self%remap_field( mesh_new, 'du_dx_a'                     , self%du_dx_a                     )
    call self%remap_field( mesh_new, 'du_dy_a'                     , self%du_dy_a                     )
    call self%remap_field( mesh_new, 'dv_dx_a'                     , self%dv_dx_a                     )
    call self%remap_field( mesh_new, 'dv_dy_a'                     , self%dv_dy_a                     )
    call reallocate_bounds( self%eta_vav_a                   , mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( self%N_a                         , mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( self%N_b                         , mesh_new%ti1, mesh_new%ti2)
    call reallocate_bounds( self%dN_dx_b                     , mesh_new%ti1, mesh_new%ti2)
    call reallocate_bounds( self%dN_dy_b                     , mesh_new%ti1, mesh_new%ti2)
    call reallocate_bounds( self%basal_friction_coefficient_b, mesh_new%ti1, mesh_new%ti2)
    call reallocate_bounds( self%tau_dx_b                    , mesh_new%ti1, mesh_new%ti2)
    call reallocate_bounds( self%tau_dy_b                    , mesh_new%ti1, mesh_new%ti2)
    call reallocate_clean ( self%u_vav_b_prev                , mesh_new%nTri             )
    call reallocate_clean ( self%v_vav_b_prev                , mesh_new%nTri             )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remap_shared_SSA_DIVA_variables

  subroutine calc_driving_stress( self, geom)

    ! In/output variables:
    class(atype_momentum_balance_solver_SSADIVA), intent(inout) :: self
    class(atype_ice_geometry_model_data),               intent(in   ) :: geom

    ! Local variables:
    character(len=*), parameter         :: routine_name = 'calc_driving_stress'
    real(dp), dimension(:), allocatable :: Hi_b
    real(dp), dimension(:), allocatable :: dHs_dx_b
    real(dp), dimension(:), allocatable :: dHs_dy_b
    integer                             :: ti

    ! Add routine to path
    call init_routine( routine_name)

    ! allocate shared memory
    allocate( Hi_b(     self%mesh%ti1:self%mesh%ti2))
    allocate( dHs_dx_b( self%mesh%ti1:self%mesh%ti2))
    allocate( dHs_dy_b( self%mesh%ti1:self%mesh%ti2))

    ! Calculate Hi, dHs/dx, and dHs/dy on the b-grid
    call map_a_b_2D( self%mesh, geom%Hi, Hi_b    )
    call ddx_a_b_2D( self%mesh, geom%Hs, dHs_dx_b)
    call ddy_a_b_2D( self%mesh, geom%Hs, dHs_dy_b)

    ! Calculate the driving stress
    do ti = self%mesh%ti1, self%mesh%ti2
      self%tau_dx_b( ti) = -ice_density * grav * Hi_b( ti) * dHs_dx_b( ti)
      self%tau_dy_b( ti) = -ice_density * grav * Hi_b( ti) * dHs_dy_b( ti)
    end do

    call checksum( self%mesh%pai_Tri, self%tau_dx_b, 'tau_dx_b')
    call checksum( self%mesh%pai_Tri, self%tau_dy_b, 'tau_dy_b')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_driving_stress

  subroutine calc_horizontal_strain_rates( self)
    !< Calculate the vertically averaged horizontal strain rates

    ! In/output variables:
    class(atype_momentum_balance_solver_SSADIVA), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'calc_horizontal_strain_rates'

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate the strain rates
    call ddx_b_a_2D( self%mesh, self%u_vav_b, self%du_dx_a)
    call ddy_b_a_2D( self%mesh, self%u_vav_b, self%du_dy_a)
    call ddx_b_a_2D( self%mesh, self%v_vav_b, self%dv_dx_a)
    call ddy_b_a_2D( self%mesh, self%v_vav_b, self%dv_dy_a)

    call checksum( self%mesh%pai_V, self%du_dx_a, 'du_dx_a')
    call checksum( self%mesh%pai_V, self%du_dy_a, 'du_dy_a')
    call checksum( self%mesh%pai_V, self%dv_dx_a, 'dv_dx_a')
    call checksum( self%mesh%pai_V, self%dv_dy_a, 'dv_dy_a')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_horizontal_strain_rates

  subroutine relax_viscosity_iterations( self, visc_it_relax)
    !< Reduce the change between velocity solutions

    ! In/output variables:
    class(atype_momentum_balance_solver_SSADIVA), intent(inout) :: self
    real(dp),                                           intent(in   ) :: visc_it_relax

    ! Local variables:
    character(len=*), parameter :: routine_name = 'relax_viscosity_iterations'
    integer                     :: ti

    ! Add routine to path
    call init_routine( routine_name)

    do ti = self%mesh%ti1, self%mesh%ti2
      self%u_vav_b( ti) = (visc_it_relax * self%u_vav_b( ti)) + ((1._dp - visc_it_relax) * self%u_vav_b_prev( ti))
      self%v_vav_b( ti) = (visc_it_relax * self%v_vav_b( ti)) + ((1._dp - visc_it_relax) * self%v_vav_b_prev( ti))
    end do

    call checksum( self%mesh%pai_Tri, self%u_vav_b, 'u_vav_b')
    call checksum( self%mesh%pai_Tri, self%v_vav_b, 'v_vav_b')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine relax_viscosity_iterations

  subroutine apply_velocity_limits( self)
    !< Limit velocities for improved stability

    ! In/output variables:
    class(atype_momentum_balance_solver_SSADIVA), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'apply_velocity_limits'
    integer                     :: ti
    real(dp)                    :: uabs

    ! Add routine to path
    call init_routine( routine_name)

    do ti = self%mesh%ti1, self%mesh%ti2

      ! Calculate absolute speed
      uabs = sqrt( self%u_vav_b( ti)**2 + self%v_vav_b( ti)**2)

      ! Reduce velocities if necessary
      if (uabs > C%vel_max) then
        self%u_vav_b( ti) = self%u_vav_b( ti) * C%vel_max / uabs
        self%v_vav_b( ti) = self%v_vav_b( ti) * C%vel_max / uabs
      end if

    end do

    call checksum( self%mesh%pai_Tri, self%u_vav_b, 'u_vav_b')
    call checksum( self%mesh%pai_Tri, self%v_vav_b, 'v_vav_b')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine apply_velocity_limits

  subroutine calc_L2_norm_uv( self, L2_uv)
    !< Calculate the L2-norm of the two consecutive velocity solutions

    ! In/output variables:
    class(atype_momentum_balance_solver_SSADIVA), intent(in   ) :: self
    real(dp),                                           intent(  out) :: L2_uv

    ! Local variables:
    character(len=*), parameter :: routine_name = 'calc_visc_iter_UV_resid'
    integer                     :: ierr
    integer                     :: ti
    real(dp)                    :: res1, res2

    ! Add routine to path
    call init_routine( routine_name)

    res1 = 0._dp
    res2 = 0._dp

    do ti = self%mesh%ti1, self%mesh%ti2

      if (.not. isnan( self%u_vav_b( ti)) .and. .not. isnan( self%v_vav_b( ti))) then

        res1 = res1 + (self%u_vav_b( ti) - self%u_vav_b_prev( ti))**2
        res1 = res1 + (self%v_vav_b( ti) - self%v_vav_b_prev( ti))**2

        res2 = res2 + (self%u_vav_b( ti) + self%u_vav_b_prev( ti))**2
        res2 = res2 + (self%v_vav_b( ti) + self%v_vav_b_prev( ti))**2

      end if

    end do

    ! Combine results from all processes
    call MPI_ALLREDUCE( MPI_IN_PLACE, res1, 1, MPI_doUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, res2, 1, MPI_doUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! Calculate L2-norm
    L2_uv = 2._dp * res1 / max( res2, 1E-8_dp)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_L2_norm_uv

end module momentum_balance_solver_SSADIVA
