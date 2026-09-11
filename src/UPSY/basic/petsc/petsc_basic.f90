module petsc_basic

#include <petsc/finclude/petscsys.h>
  use precisions, only: dp
  use CSR_matrix_mod, only: type_CSR_matrix_dp
  use petsc, only: PetscErrorF, tVec, VecCreate, VecSetSizes, VecSetFromOptions, VecSetValues, &
    VecAssemblyBegin, VecAssemblyEnd, VecGetLocalSize, VecGetSize, VecGetOwnershipRange, &
    VecGetValues, INSERT_VALUES, PETSC_COMM_WORLD, &
    tMat, MatGetSize, MatGetLocalSize, MatGetOwnershipRange, &
    MatGetRow, MatRestoreRow, MatCreateMPIAIJWithArrays, &
    MatAssemblyBegin, MatAssemblyEnd, MAT_FINAL_ASSEMBLY, MatMult, VecDestroy, &
    tKSP, tPC, KSPcreate, KSPSetOperators, KSPSetInitialGuessNonzero, &
    PETSC_TRUE, KSPSetTolerances, PETSC_DEFAULT_REAL, KSPSetFromOptions, &
    KSPGetIterationNumber, KSPDestroy, KSPGetType, PCGetType, &
    KSPSolve, KSPGMRES, KSPLGMRES, KSPFGMRES, KSPPGMRES, KSPBICG, KSPBCGS, KSPIBCGS, &
    PCBJACOBI, PCASM, PCGASM, PCGAMG, PCNONE, KSPSetType, PCSetType, KSPGetPC, MatDestroy
  use assertions_basic, only: assert
  use mpi_basic, only: par
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use crash_mod, only: crash
  use mpi_f08, only: MPI_ALLGATHER, MPI_INTEGER, MPI_COMM_WORLD
  use string_module, only: colour_string

  implicit none

  private

  public :: vec_double2petsc, vec_petsc2double
  public :: mat_petsc2CSR, mat_CSR2petsc
  public :: multiply_PETSc_matrix_with_vector_1D, multiply_PETSc_matrix_with_vector_2D
  public :: solve_matrix_equation_CSR_PETSc, solve_matrix_equation_PETSc

  ! Interfaces for procedures defined in submodules
  interface

    module subroutine vec_double2petsc( xx, x)
      real(dp), dimension(:), intent(in)  :: xx
      type(tVec),             intent(out) :: x
    end subroutine vec_double2petsc

    module subroutine vec_petsc2double( x, xx)
      type(tVec),             intent(in)  :: x
      real(dp), dimension(:), intent(out) :: xx
    end subroutine vec_petsc2double

    module subroutine mat_petsc2CSR( A, AA)
      type(tMat),               intent(in   ) :: A
      type(type_CSR_matrix_dp), intent(  out) :: AA
    end subroutine mat_petsc2CSR

    module subroutine mat_CSR2petsc( AA, A)
      type(type_CSR_matrix_dp), intent(in   ) :: AA
      type(tMat),               intent(  out) :: A
    end subroutine mat_CSR2petsc

    module subroutine multiply_PETSc_matrix_with_vector_1D( A, xx, yy)
      type(tMat),                     intent(in   ) :: A
      real(dp), dimension(:), target, intent(in   ) :: xx
      real(dp), dimension(:), target, intent(  out) :: yy
    end subroutine multiply_PETSc_matrix_with_vector_1D

    module subroutine multiply_PETSc_matrix_with_vector_2D( A, xx, yy)
      type(tMat),               intent(in   ) :: A
      real(dp), dimension(:,:), intent(in   ) :: xx
      real(dp), dimension(:,:), intent(  out) :: yy
    end subroutine multiply_PETSc_matrix_with_vector_2D

    module subroutine solve_matrix_equation_CSR_PETSc( AA, bb, xx, rtol, abstol, &
      n_Axb_its,  PETSc_KSPtype, PETSc_PCtype)
      type(type_CSR_matrix_dp),            intent(in   ) :: AA
      real(dp), dimension(:),              intent(in   ) :: bb
      real(dp), dimension(:),              intent(inout) :: xx
      real(dp),                            intent(in   ) :: rtol, abstol
      integer,                             intent(out)   :: n_Axb_its
      character(len=*), optional,          intent(in   ) :: PETSc_KSPtype
      character(len=*), optional,          intent(in   ) :: PETSc_PCtype
    end subroutine solve_matrix_equation_CSR_PETSc

    module subroutine solve_matrix_equation_PETSc( A, bb, xx, rtol, abstol, &
      n_Axb_its,  PETSc_KSPtype, PETSc_PCtype)
      type(tMat),                          intent(in   ) :: A
      real(dp), dimension(:),              intent(in   ) :: bb
      real(dp), dimension(:),              intent(inout) :: xx
      real(dp),                            intent(in   ) :: rtol, abstol
      integer,                             intent(out)   :: n_Axb_its
      character(len=*), optional,          intent(in   ) :: PETSc_KSPtype
      character(len=*), optional,          intent(in   ) :: PETSc_PCtype
    end subroutine solve_matrix_equation_PETSc

  end interface

end module petsc_basic
