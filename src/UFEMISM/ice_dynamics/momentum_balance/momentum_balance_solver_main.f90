module momentum_balance_solver_main

  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use crash_mod, only: crash
  use momentum_balance_solver_basic, only: atype_momentum_balance_solver
  use momentum_balance_solver_none, only: type_momentum_balance_solver_none
  use momentum_balance_solver_SIA, only: type_momentum_balance_solver_SIA
  use momentum_balance_solver_SSA, only: type_momentum_balance_solver_SSA
  use momentum_balance_solver_SSA_FD_SNES, only: type_momentum_balance_solver_SSA_FD_SNES
  use momentum_balance_solver_SSA_FEM_PETSc, only: type_momentum_balance_solver_SSA_FEM_PETSc
  use momentum_balance_solver_SIASSA, only: type_momentum_balance_solver_SIASSA
  use momentum_balance_solver_DIVA, only: type_momentum_balance_solver_DIVA
  use momentum_balance_solver_BPA, only: type_momentum_balance_solver_BPA
  use momentum_balance_solver_hybrid_DIVA_BPA, only: type_momentum_balance_solver_hybrid_DIVA_BPA

  implicit none

  private

  public :: atype_momentum_balance_solver, create_momentum_balance_solver

contains

  subroutine create_momentum_balance_solver( momentum_balance_solver, choice_stress_balance_approximation)
    !< Allocate a concrete implementation of a momentum balance solver

    ! In/output variables:
    class(atype_momentum_balance_solver), allocatable, intent(inout) :: momentum_balance_solver
    character(len=*),                                  intent(in   ) :: choice_stress_balance_approximation

    ! Local variables:
    character(len=*), parameter :: routine_name = 'create_momentum_balance_solver'

    ! Add routine to call stack
    call init_routine( routine_name)

    select case (choice_stress_balance_approximation)
    case default
      call crash('unknown choice_stress_balance_approximation "' // trim( choice_stress_balance_approximation) // '"!')
    case ('none')
      allocate( type_momentum_balance_solver_none :: momentum_balance_solver)
    case ('SIA')
      allocate( type_momentum_balance_solver_SIA :: momentum_balance_solver)
    case ('SSA')
      allocate( type_momentum_balance_solver_SSA :: momentum_balance_solver)
    case ('SSA_FD_SNES')
      allocate( type_momentum_balance_solver_SSA_FD_SNES :: momentum_balance_solver)
    case ('SSA_FEM_PETSc')
      allocate( type_momentum_balance_solver_SSA_FEM_PETSc :: momentum_balance_solver)
    case ('SIA/SSA')
      allocate( type_momentum_balance_solver_SIASSA :: momentum_balance_solver)
    case ('DIVA')
      allocate( type_momentum_balance_solver_DIVA :: momentum_balance_solver)
    case ('BPA')
      allocate( type_momentum_balance_solver_BPA :: momentum_balance_solver)
    case ('hybrid DIVA/BPA')
      allocate( type_momentum_balance_solver_hybrid_DIVA_BPA :: momentum_balance_solver)
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_momentum_balance_solver

end module momentum_balance_solver_main