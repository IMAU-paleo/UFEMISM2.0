module momentum_balance_solver_SSA_FD_SNES

  ! Routines for calculating ice velocities using the Shallow Shelf Approximation (SSA),
  ! discretised with the in-house finite-difference operators (exactly as in
  ! momentum_balance_solver_SSA) but solved with PETSc's SNES non-linear solver
  ! instead of the hand-rolled viscosity (Picard) iteration.
  !
  ! This is "Tier 1" of SSA_FD_SNES_implementation_plan.md (next to this file):
  ! defect-correction Newton. The non-linear residual is F(u) = A(u) u - b(u), where
  ! A(u) and b(u) are the linear system assembled by
  ! atype_momentum_balance_solver_SSADIVA%assemble_SSA_DIVA_linearised_matrix_eq -
  ! the same system the SSA solver builds once per viscosity iteration. The Jacobian
  ! handed to SNES is that same A(u) (the "Picard operator"), so no analytic tangent
  ! is derived. SNES's line search and convergence control replace the manual
  ! relaxation / divergence back-off logic.
  !
  ! It is selected with
  !   choice_stress_balance_approximation = 'SSA_FD_SNES'
  ! and is added alongside the existing solvers without affecting any of them. The
  ! type extends type_momentum_balance_solver_SSA, so all of the finite-difference
  ! machinery (driving stress, strain rates, effective viscosity, applied basal
  ! friction, field allocation, remap, restart files) is inherited unchanged; only
  ! the non-linear solve in run_momentum_balance_solver is new.
  !
  ! Velocities are defined on the b-grid (triangles), exactly like the SSA solver,
  ! so ice-thickness evolution and all downstream code are unaffected.
  !
  ! Implementation status: STEPS 1-2 of the plan (skeleton + registration, plus the
  ! "recompute coefficients from velocity" helper). The SNES driver, residual and
  ! Jacobian callbacks are added in steps 3-6; run_momentum_balance_solver is a stub
  ! until then.

  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash
  use model_configuration, only: C
  use ice_model_data, only: atype_ice_model_data
  use ice_geometry_model_data, only: atype_ice_geometry_model_data
  use bed_roughness_model_types, only: type_bed_roughness_model
  use momentum_balance_solver_SSA, only: type_momentum_balance_solver_SSA

  implicit none

  private

  public :: type_momentum_balance_solver_SSA_FD_SNES

  type, extends(type_momentum_balance_solver_SSA) :: type_momentum_balance_solver_SSA_FD_SNES

    contains

      ! Only the solver name and the non-linear solve differ from the SSA solver;
      ! everything else is inherited from type_momentum_balance_solver_SSA.
      procedure, public :: get_momentum_balance_solver_name => get_momentum_balance_solver_name_SSA_FD_SNES
      procedure, public :: run_momentum_balance_solver      => momentum_balance_solver_SSA_FD_SNES_run

  end type type_momentum_balance_solver_SSA_FD_SNES

contains

  function get_momentum_balance_solver_name_SSA_FD_SNES( self) result( model_name)
    class(type_momentum_balance_solver_SSA_FD_SNES), intent(in) :: self
    character(len=:), allocatable                               :: model_name
    model_name = 'SSA_FD_SNES'
  end function get_momentum_balance_solver_name_SSA_FD_SNES

  subroutine momentum_balance_solver_SSA_FD_SNES_run( self, ice, geom, bed_roughness, &
    BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b, BC_prescr_mask_bk, BC_prescr_u_bk, BC_prescr_v_bk)
    !< Calculate ice velocities by solving the non-linear SSA with PETSc's SNES,
    !< using the in-house finite-difference discretisation for the residual and
    !< the Picard operator as the (approximate) Jacobian.

    ! In/output variables:
    class(type_momentum_balance_solver_SSA_FD_SNES), intent(inout) :: self
    class(atype_ice_model_data),                     intent(inout) :: ice
    class(atype_ice_geometry_model_data),            intent(in   ) :: geom
    type(type_bed_roughness_model),                  intent(in   ) :: bed_roughness
    integer,  dimension(:  ), optional,             intent(in   ) :: BC_prescr_mask_b   ! Mask of triangles where velocity is prescribed
    real(dp), dimension(:  ), optional,             intent(in   ) :: BC_prescr_u_b      ! Prescribed velocities in the x-direction
    real(dp), dimension(:  ), optional,             intent(in   ) :: BC_prescr_v_b      ! Prescribed velocities in the y-direction
    integer,  dimension(:,:), optional,             intent(in   ) :: BC_prescr_mask_bk  ! Mask of triangles where velocity is prescribed
    real(dp), dimension(:,:), optional,             intent(in   ) :: BC_prescr_u_bk     ! Prescribed velocities in the x-direction
    real(dp), dimension(:,:), optional,             intent(in   ) :: BC_prescr_v_bk     ! Prescribed velocities in the y-direction

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_SSA_FD_SNES_run'

    ! Add routine to path
    call init_routine( routine_name)

    call crash('the SSA_FD_SNES solver is not implemented yet - only steps 1-2 ' // &
      'of SSA_FD_SNES_implementation_plan.md are done; the SNES driver, residual ' // &
      'and Jacobian callbacks are steps 3-6')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_SSA_FD_SNES_run

end module momentum_balance_solver_SSA_FD_SNES
