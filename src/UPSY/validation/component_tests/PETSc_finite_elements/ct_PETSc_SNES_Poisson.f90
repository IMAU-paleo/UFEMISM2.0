module ct_PETSc_SNES_Poisson

#include <petsc/finclude/petscsys.h>
  use precisions, only: dp
  use iso_c_binding, only: c_bool, c_char, c_double, c_funloc, c_funptr, c_int, c_intptr_t, c_loc, c_null_char, &
    c_null_funptr, c_null_ptr, c_ptr, c_f_pointer
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use crash_mod, only: crash, warning
  use mpi_f08, only: MPI_ALLTOALL, MPI_ALLTOALLV, MPI_COMM_WORLD, MPI_DOUBLE_PRECISION, MPI_INTEGER
  use mpi_basic, only: par
  use mesh_types, only: type_mesh
  use netcdf_io_main, only: open_existing_netcdf_file_for_reading, setup_mesh_from_file, &
    close_netcdf_file, save_variable_as_netcdf_dp_1D
  use petsc, only: PetscErrorF, PETSC_COMM_SELF, PETSC_COMM_WORLD, PETSC_TRUE, PETSC_FALSE, &
    PETSC_NULL_DMLABEL, PETSC_NULL_VEC, DM_BC_ESSENTIAL, tDM, tVec, tMat, tSNES, tKSP, tPC, &
    tPetscObject, tPetscFE, tPetscDS, tDMLabel, tPetscSection, DMDestroy, MatDestroy, &
    VecDestroy, SNESDestroy, PetscFEDestroy, PetscFECreateLagrange, PetscObjectSetName, DMSetField, &
    DMCreateDS, DMGetDS, PetscDSSetConstants, DMCreateLabel, DMGetLabel, DMPlexMarkBoundaryFaces, &
    SNESCreate, SNESSetDM, SNESSetType, SNESSetTolerances, SNESGetKSP, SNESNEWTONLS, &
    KSPSetType, KSPGetPC, KSPPREONLY, PCSetType, PCLU, DMCreateGlobalVector, DMCreateLocalVector, DMCreateMatrix, &
    VecDuplicate, &
    DMGlobalToLocalBegin, DMGlobalToLocalEnd, DMGetLocalSection, PetscSectionGetOffset, VecGetArrayRead, &
    VecRestoreArrayRead, &
    DMPlexGetDepthStratum, DMLabelGetValue, &
    INSERT_VALUES, VecSet, VecNorm, MatNorm, NORM_INFINITY, SNESComputeFunction, SNESComputeJacobian, MatGetDiagonal, &
    SNESSolve, SNESGetIterationNumber
  use petsc_dmplex, only: mesh_to_dmplex, dmplex_upsy_vertex_id_label_name
  use string_module, only: strrep

  implicit none

  private

  public :: ct_solve_Poisson_eq_with_PETSc_SNES

  real(dp), dimension(4) :: poisson_domain
  integer                :: poisson_g3_call_count = 0
  integer(c_int), parameter :: dm_bc_essential_value = 1_c_int
  real(dp), parameter :: poisson_snes_relative_tolerance = 1.0e-10_dp
  real(dp), parameter :: poisson_snes_absolute_tolerance = 1.0e-10_dp
  real(dp), parameter :: poisson_snes_step_tolerance     = 1.0e-12_dp
  integer,  parameter :: poisson_snes_max_iterations     = 20
  integer,  parameter :: poisson_snes_max_function_evaluations = 1000
  character(kind=c_char), dimension(9), parameter :: boundary_label_name = [character(kind=c_char) :: &
    'b', 'o', 'u', 'n', 'd', 'a', 'r', 'y', c_null_char]

  interface
    subroutine DMPlexSetSNESLocalFEM( dm, has_boundary, ctx, ierr)
      import :: c_bool, c_intptr_t, tDM
      type(tDM),              intent(inout) :: dm
      logical(kind=c_bool),   intent(in)    :: has_boundary
      integer(c_intptr_t),    intent(in)    :: ctx
      integer,                intent(out)   :: ierr
    end subroutine DMPlexSetSNESLocalFEM

    integer(c_int) function petsc_ds_set_residual( ds, field, f0, f1) bind(C, name='PetscDSSetResidual')
      import :: c_funptr, c_int, c_intptr_t
      integer(c_intptr_t), value :: ds
      integer(c_intptr_t), value :: field
      type(c_funptr),      value :: f0, f1
    end function petsc_ds_set_residual

    integer(c_int) function petsc_ds_set_jacobian( ds, field_test, field_trial, g0, g1, g2, g3) &
      bind(C, name='PetscDSSetJacobian')
      import :: c_funptr, c_int, c_intptr_t
      integer(c_intptr_t), value :: ds
      integer(c_intptr_t), value :: field_test, field_trial
      type(c_funptr),      value :: g0, g1, g2, g3
    end function petsc_ds_set_jacobian

    integer(c_int) function petsc_ds_has_jacobian( ds, has_jacobian) bind(C, name='PetscDSHasJacobian')
      import :: c_bool, c_int, c_intptr_t
      integer(c_intptr_t), value       :: ds
      logical(kind=c_bool), intent(out) :: has_jacobian
    end function petsc_ds_has_jacobian

    integer(c_int) function petsc_ds_set_exact_solution( ds, field, function, ctx) &
      bind(C, name='PetscDSSetExactSolution')
      import :: c_funptr, c_int, c_intptr_t, c_ptr
      integer(c_intptr_t), value :: ds
      integer(c_intptr_t), value :: field
      type(c_funptr),      value :: function
      type(c_ptr),          value :: ctx
    end function petsc_ds_set_exact_solution

    integer(c_int) function dm_add_boundary( dm, boundary_type, name, label, nvalues, values, field, ncomponents, &
      components, boundary_function, boundary_time_derivative, ctx, boundary_index) bind(C, name='DMAddBoundary')
      import :: c_char, c_funptr, c_int, c_intptr_t, c_ptr
      integer(c_intptr_t),                value :: dm, label
      integer(c_int),                     value :: boundary_type
      integer(c_intptr_t),                value :: nvalues, field, ncomponents
      character(kind=c_char), dimension(*), intent(in) :: name
      type(c_ptr),                        value :: values, components, ctx
      type(c_funptr),                     value :: boundary_function, boundary_time_derivative
      integer(c_intptr_t),             intent(out) :: boundary_index
    end function dm_add_boundary

    integer(c_int) function snes_get_converged_reason( snes, reason) bind(C, name='SNESGetConvergedReason')
      import :: c_int, c_intptr_t
      integer(c_intptr_t), value :: snes
      integer(c_int),      intent(out) :: reason
    end function snes_get_converged_reason

    integer(c_int) function snes_set_jacobian( snes, jacobian, preconditioner, function, ctx) &
      bind(C, name='SNESSetJacobian')
      import :: c_funptr, c_int, c_intptr_t, c_ptr
      integer(c_intptr_t), value :: snes, jacobian, preconditioner
      type(c_funptr),      value :: function
      type(c_ptr),         value :: ctx
    end function snes_set_jacobian

    integer(c_int) function petsc_section_has_constraints( section, has_constraints) &
      bind(C, name='PetscSectionHasConstraints')
      import :: c_bool, c_int, c_intptr_t
      integer(c_intptr_t),      value       :: section
      logical(kind=c_bool), intent(out) :: has_constraints
    end function petsc_section_has_constraints

    integer(c_int) function petsc_section_get_storage_size( section, storage_size) &
      bind(C, name='PetscSectionGetStorageSize')
      import :: c_int, c_intptr_t
      integer(c_intptr_t),   value       :: section
      integer(c_intptr_t), intent(out) :: storage_size
    end function petsc_section_get_storage_size

    integer(c_int) function petsc_section_get_constrained_storage_size( section, constrained_storage_size) &
      bind(C, name='PetscSectionGetConstrainedStorageSize')
      import :: c_int, c_intptr_t
      integer(c_intptr_t),   value       :: section
      integer(c_intptr_t), intent(out) :: constrained_storage_size
    end function petsc_section_get_constrained_storage_size
  end interface

contains

  subroutine ct_solve_Poisson_eq_with_PETSc_SNES( foldername_output, test_mesh_filenames, test_grid_filenames)

    ! In/output variables:
    character(len=*),               intent(in) :: foldername_output
    character(len=*), dimension(:), intent(in) :: test_mesh_filenames
    character(len=*), dimension(:), intent(in) :: test_grid_filenames

    ! Local variables:
    character(len=*), parameter   :: routine_name = 'ct_solve_Poisson_eq_with_PETSc_SNES'
    integer                       :: i_mesh
    type(type_mesh), allocatable  :: mesh
    character(len=:), allocatable :: filename_mesh
    integer                       :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    if (par%primary) write(0,*) '   Running PETSc SNES test: solve a simple Poisson equation'
    if (par%primary) write(0,*) ''

    ! Meshes read from files
    do i_mesh = 1, size( test_mesh_filenames)
      filename_mesh = test_mesh_filenames( i_mesh)
      allocate( mesh)
      call open_existing_netcdf_file_for_reading( filename_mesh, ncid)
      call setup_mesh_from_file( filename_mesh, ncid, mesh)
      call close_netcdf_file( ncid)
      call ct_solve_Poisson_eq_with_PETSc_SNES_on_mesh( foldername_output, mesh)
      deallocate( mesh)
    end do

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine ct_solve_Poisson_eq_with_PETSc_SNES

  subroutine ct_solve_Poisson_eq_with_PETSc_SNES_on_mesh( foldername_output, mesh)

    ! In/output variables:
    character(len=*), intent(in) :: foldername_output
    type(type_mesh),  intent(in) :: mesh

    ! Local variables:
    character(len=*), parameter   :: routine_name = 'ct_solve_Poisson_eq_with_PETSc_SNES_on_mesh'
    type(tDM)                     :: dm
    type(tVec)                    :: solution, residual
    type(tMat)                    :: jacobian
    type(tSNES)                   :: snes
    type(tPetscFE)                :: fe
    type(tPetscObject)            :: fe_object
    type(tPetscDS)                :: ds
    type(tDMLabel)                :: boundary_label
    type(tPetscSection)           :: local_section
    integer                       :: ierr, snes_iterations
    logical(kind=c_bool)          :: ds_has_jacobian, local_section_has_constraints
    integer(c_int)                :: snes_reason
    integer(c_intptr_t)           :: boundary_index
    integer(c_intptr_t)           :: local_section_storage_size, local_section_constrained_storage_size
    integer(c_intptr_t), target, dimension(1) :: boundary_ids, unused_components
    integer(c_intptr_t)           :: no_context
    type(type_poisson_context), target :: poisson_context
    real(dp)                      :: initial_residual_norm, jacobian_diagonal_norm, jacobian_norm, solution_norm
    real(dp), dimension(:), allocatable, target :: solution_on_vertices
    character(len=:), allocatable :: mesh_name_cleaned
    character(len=:), allocatable :: filename

    ! Add routine to call stack
    call init_routine( routine_name)

    if (par%primary) write(0,*) '    Using PETSc SNES to solve the Poisson equation on mesh ', trim( mesh%name), '...'

    call mesh_to_dmplex( mesh, dm)

    ! The analytical solution is defined on the mesh bounding box. For a rectangular
    ! mesh its value is zero on every boundary face.
    poisson_domain( 1) = minval( mesh%V(:,1))
    poisson_domain( 2) = minval( mesh%V(:,2))
    poisson_domain( 3) = maxval( mesh%V(:,1)) - poisson_domain( 1)
    poisson_domain( 4) = maxval( mesh%V(:,2)) - poisson_domain( 2)
    no_context = 0_c_intptr_t

    ! Attach one scalar P1 field, u, to the DMPlex.
    PetscCall( PetscFECreateLagrange( PETSC_COMM_SELF, 2, 1, PETSC_TRUE, 1, -1, fe, ierr))
    PetscCall( PetscObjectSetName( fe, 'u', ierr))
    PetscObjectSpecificCast( fe_object, fe)
    PetscCall( DMSetField( dm, 0, PETSC_NULL_DMLABEL, fe_object, ierr))
    PetscCall( DMCreateDS( dm, ierr))
    PetscCall( DMGetDS( dm, ds, ierr))
    PetscCall( PetscDSSetConstants( ds, 4, poisson_domain, ierr))

    ! Define Poisson weak form
    !     Strong form: -grad^2 u = f within the domain, u = 0 on the border
    !            with: f = pi^2 (1 / w_x^2 + 1 / w_y^2) sin( pi x / w_x) sin( pi y / w_y)
    !            where [w_x,w_y] is the size of the rectangular domain.
    ierr = petsc_ds_set_residual( ds%v, 0_c_intptr_t, c_funloc( poisson_f0), c_funloc( poisson_f1))
    CHKERRQ( ierr)
    ierr = petsc_ds_set_jacobian( ds%v, 0_c_intptr_t, 0_c_intptr_t, c_null_funptr, c_null_funptr, c_null_funptr, &
      c_funloc( poisson_g3))
    CHKERRQ( ierr)
    ierr = petsc_ds_has_jacobian( ds%v, ds_has_jacobian)
    CHKERRQ( ierr)
    ierr = petsc_ds_set_exact_solution( ds%v, 0_c_intptr_t, c_funloc( poisson_exact_solution), c_null_ptr)
    CHKERRQ( ierr)

    ! Mark all exterior faces and impose the manufactured solution on them. The
    ! prescribed value is identically zero for the rectangular-domain problem.
    PetscCall( DMCreateLabel( dm, 'boundary', ierr))
    PetscCall( DMGetLabel( dm, 'boundary', boundary_label, ierr))
    PetscCall( DMPlexMarkBoundaryFaces( dm, 1, boundary_label, ierr))
    boundary_ids = [1_c_intptr_t]
    unused_components = 0_c_intptr_t
    ierr = dm_add_boundary( dm%v, dm_bc_essential_value, boundary_label_name, boundary_label%v, 1_c_intptr_t, &
      c_loc( boundary_ids), 0_c_intptr_t, 0_c_intptr_t, c_loc( unused_components), c_funloc( poisson_dirichlet_boundary), &
      c_null_funptr, c_null_ptr, boundary_index)
    CHKERRQ( ierr)

    ! Let DMPlex assemble the residual and Jacobian from the PetscDS weak form.
    PetscCall( SNESCreate( PETSC_COMM_WORLD, snes, ierr))
    PetscCall( SNESSetDM( snes, dm, ierr))
    PetscCall( DMPlexSetSNESLocalFEM( dm, PETSC_FALSE, no_context, ierr))
    PetscCall( DMCreateMatrix( dm, jacobian, ierr))
    PetscCall( DMGetLocalSection( dm, local_section, ierr))
    ierr = petsc_section_has_constraints( local_section%v, local_section_has_constraints)
    CHKERRQ( ierr)
    ierr = petsc_section_get_storage_size( local_section%v, local_section_storage_size)
    CHKERRQ( ierr)
    ierr = petsc_section_get_constrained_storage_size( local_section%v, local_section_constrained_storage_size)
    CHKERRQ( ierr)
    ierr = snes_set_jacobian( snes%v, jacobian%v, jacobian%v, c_null_funptr, c_null_ptr)
    CHKERRQ( ierr)
    call configure_PETSc_SNES_for_Poisson( snes)

    PetscCall( DMCreateGlobalVector( dm, solution, ierr))
    PetscCall( VecSet( solution, 0._dp, ierr))
    PetscCall( VecDuplicate( solution, residual, ierr))
    PetscCall( SNESComputeFunction( snes, solution, residual, ierr))
    PetscCall( VecNorm( residual, NORM_INFINITY, initial_residual_norm, ierr))
    poisson_g3_call_count = 0
    PetscCall( SNESComputeJacobian( snes, solution, jacobian, jacobian, ierr))
    PetscCall( MatNorm( jacobian, NORM_INFINITY, jacobian_norm, ierr))
    PetscCall( MatGetDiagonal( jacobian, residual, ierr))
    PetscCall( VecNorm( residual, NORM_INFINITY, jacobian_diagonal_norm, ierr))
    PetscCall( SNESSolve( snes, PETSC_NULL_VEC, solution, ierr))
    PetscCall( VecNorm( solution, NORM_INFINITY, solution_norm, ierr))
    PetscCall( SNESGetIterationNumber( snes, snes_iterations, ierr))
    ierr = snes_get_converged_reason( snes%v, snes_reason)
    CHKERRQ( ierr)
    if (par%primary) write(0,*) '      PetscDS has Jacobian  = ', ds_has_jacobian
    if (par%primary) write(0,*) '      poisson_g3 calls      = ', poisson_g3_call_count
    if (par%primary) write(0,*) '      Local section constrained = ', local_section_has_constraints
    if (par%primary) write(0,*) '      Local section storage size = ', local_section_storage_size
    if (par%primary) write(0,*) '      Local section free dofs    = ', local_section_constrained_storage_size
    if (par%primary) write(0,*) '      Initial residual norm = ', initial_residual_norm
    if (par%primary) write(0,*) '      Jacobian matrix norm   = ', jacobian_norm
    if (par%primary) write(0,*) '      Jacobian diagonal norm = ', jacobian_diagonal_norm
    if (par%primary) write(0,*) '      PETSc solution norm   = ', solution_norm
    if (par%primary) write(0,*) '      SNES iterations       = ', snes_iterations
    if (par%primary) write(0,*) '      SNES convergence code = ', snes_reason
    call copy_PETSc_solution_to_mesh_vertices( dm, solution, mesh, solution_on_vertices)
    write(0,*) 'max = ', maxval( solution_on_vertices)

    mesh_name_cleaned = trim( mesh%name)
    mesh_name_cleaned = strrep( mesh_name_cleaned, '"', '')
    mesh_name_cleaned = strrep( mesh_name_cleaned, '.', '_')
    mesh_name_cleaned = strrep( mesh_name_cleaned, '/', '_')
    filename = trim( mesh_name_cleaned) // '_solution'
    call save_variable_as_netcdf_dp_1D( foldername_output, solution_on_vertices, filename)

    ! ========================================================================
    ! NEXT STEPS
    ! ========================================================================
    ! 1. Run this component test on one small mesh and inspect PETSc's solver
    !    behaviour. Solver choices and tolerances are set in
    !    configure_PETSc_SNES_for_Poisson below; move those named parameters to
    !    the UPSY configuration type when this becomes model functionality.
    !
    ! 2. Add an L2-error calculation against poisson_exact_solution. PETSc's
    !    Fortran module in the current environment does not expose
    !    DMComputeL2Diff, but the function is exported by libpetsc. Follow the
    !    direct bind(C) pattern used above for the PetscDS registrations, pass
    !    dm%v and solution%v, and use poisson_exact_solution via c_funloc.
    !
    ! 3. Check convergence under mesh refinement. The default PetscFE is a P1
    !    simplex element, so the L2 error should decrease at approximately
    !    second order for this smooth manufactured solution.
    !
    ! 4. The manufactured solution and the zero Dirichlet condition assume an
    !    axis-aligned rectangular boundary. Replace them before using meshes
    !    with another outer boundary, or prescribe the exact value on that
    !    boundary instead of zero.
    !
    ! 5. DMPlex is distributed before the solve. mesh_to_dmplex preserves the
    !    original UPSY vertex ID in a DMLabel, and the copy-back routine
    !    combines the local vertex values collectively in UPSY vertex order.
    !    Recheck this mapping after changes to the DMPlex construction.
    !
    ! 6. The direct C bindings in this module work around PETSc routines whose
    !    legacy Fortran wrappers are absent from the installed PETSc 3.25.5
    !    package. Prefer generated PETSc Fortran bindings if a later PETSc
    !    installation supplies them.


    ! Clean up after yourself
    PetscCall( VecDestroy( residual, ierr))
    PetscCall( VecDestroy( solution, ierr))
    PetscCall( MatDestroy( jacobian, ierr))
    PetscCall( SNESDestroy( snes, ierr))
    PetscCall( PetscFEDestroy( fe, ierr))
    PetscCall( DMDestroy( dm, ierr))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine ct_solve_Poisson_eq_with_PETSc_SNES_on_mesh

  subroutine configure_PETSc_SNES_for_Poisson( snes)

    type(tSNES), intent(inout) :: snes

    type(tKSP) :: ksp
    type(tPC)  :: preconditioner
    integer    :: ierr

    PetscCall( SNESSetType( snes, SNESNEWTONLS, ierr))
    call SNESSetTolerances( snes, poisson_snes_absolute_tolerance, poisson_snes_relative_tolerance, &
      poisson_snes_step_tolerance, poisson_snes_max_iterations, poisson_snes_max_function_evaluations, ierr)
    CHKERRQ( ierr)

    PetscCall( SNESGetKSP( snes, ksp, ierr))
    PetscCall( KSPSetType( ksp, KSPPREONLY, ierr))
    PetscCall( KSPGetPC( ksp, preconditioner, ierr))
    PetscCall( PCSetType( preconditioner, PCLU, ierr))

  end subroutine configure_PETSc_SNES_for_Poisson

  subroutine copy_PETSc_solution_to_mesh_vertices( dm, solution, mesh, solution_on_vertices)

    type(tDM),                                  intent(in)  :: dm
    type(tVec),                                 intent(in)  :: solution
    type(type_mesh),                            intent(in)  :: mesh
    real(dp), dimension(:), allocatable, target, intent(out) :: solution_on_vertices

    type(tVec)                 :: local_solution
    type(tPetscSection)        :: local_section
    type(tDMLabel)             :: upsy_vertex_id_label
    real(dp), dimension(:), pointer :: local_solution_values
    integer, dimension(:), allocatable :: copies_per_vertex, send_counts, receive_counts, send_displacements, &
      receive_displacements, send_positions, send_vertex_ids, receive_vertex_ids
    real(dp), dimension(:), allocatable :: send_values, receive_values
    integer                    :: ierr, vi, point, local_offset, vertex_start, vertex_end, ip, destination, &
      send_count, receive_count, send_index, receive_index

    allocate( solution_on_vertices( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( copies_per_vertex( mesh%vi1:mesh%vi2), source = 0)
    allocate( send_counts( 0:par%n-1), receive_counts( 0:par%n-1), source = 0)
    allocate( send_displacements( 0:par%n-1), receive_displacements( 0:par%n-1), send_positions( 0:par%n-1))

    PetscCall( DMCreateLocalVector( dm, local_solution, ierr))
    PetscCall( DMGlobalToLocalBegin( dm, solution, INSERT_VALUES, local_solution, ierr))
    PetscCall( DMGlobalToLocalEnd(   dm, solution, INSERT_VALUES, local_solution, ierr))
    PetscCall( DMGetLocalSection( dm, local_section, ierr))
    PetscCall( DMGetLabel( dm, dmplex_upsy_vertex_id_label_name, upsy_vertex_id_label, ierr))
    PetscCall( DMPlexGetDepthStratum( dm, 0, vertex_start, vertex_end, ierr))

    do point = vertex_start, vertex_end - 1
      PetscCall( DMLabelGetValue( upsy_vertex_id_label, point, vi, ierr))
      if (vi < 1 .or. vi > mesh%nV) call crash('DMPlex vertex lacks a valid UPSY vertex ID')
      destination = mesh%V_owning_process( vi)
      if (destination < 0 .or. destination >= par%n) call crash('UPSY vertex has an invalid owning process')
      send_counts( destination) = send_counts( destination) + 1
    end do

    send_displacements( 0) = 0
    do ip = 1, par%n-1
      send_displacements( ip) = send_displacements( ip-1) + send_counts( ip-1)
    end do
    send_count = sum( send_counts)
    send_positions = send_displacements
    allocate( send_vertex_ids( max( 1, send_count)), send_values( max( 1, send_count)))

    PetscCall( VecGetArrayRead( local_solution, local_solution_values, ierr))
    do point = vertex_start, vertex_end - 1
      PetscCall( DMLabelGetValue( upsy_vertex_id_label, point, vi, ierr))
      destination = mesh%V_owning_process( vi)
      send_index = send_positions( destination) + 1
      PetscCall( PetscSectionGetOffset( local_section, point, local_offset, ierr))
      send_vertex_ids( send_index) = vi
      send_values( send_index) = local_solution_values( local_offset + 1)
      send_positions( destination) = send_positions( destination) + 1
    end do

    PetscCall( VecRestoreArrayRead( local_solution, local_solution_values, ierr))
    PetscCall( VecDestroy( local_solution, ierr))

    call MPI_ALLTOALL( send_counts, 1, MPI_INTEGER, receive_counts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    receive_displacements( 0) = 0
    do ip = 1, par%n-1
      receive_displacements( ip) = receive_displacements( ip-1) + receive_counts( ip-1)
    end do
    receive_count = sum( receive_counts)
    allocate( receive_vertex_ids( max( 1, receive_count)), receive_values( max( 1, receive_count)))

    call MPI_ALLTOALLV( send_vertex_ids, send_counts, send_displacements, MPI_INTEGER, receive_vertex_ids, receive_counts, &
      receive_displacements, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLTOALLV( send_values, send_counts, send_displacements, MPI_DOUBLE_PRECISION, receive_values, receive_counts, &
      receive_displacements, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)

    do receive_index = 1, receive_count
      vi = receive_vertex_ids( receive_index)
      if (vi < mesh%vi1 .or. vi > mesh%vi2) call crash('DMPlex solution was sent to the wrong UPSY process')
      solution_on_vertices( vi) = solution_on_vertices( vi) + receive_values( receive_index)
      copies_per_vertex( vi) = copies_per_vertex( vi) + 1
    end do

    if (any( copies_per_vertex == 0)) call crash('DMPlex distribution omitted an UPSY vertex')
    solution_on_vertices = solution_on_vertices / real( copies_per_vertex, dp)

  end subroutine copy_PETSc_solution_to_mesh_vertices

  integer(c_int) function poisson_exact_solution( dim, time, x, ncomp, u, ctx) bind(C)

    integer(c_intptr_t),   value :: dim, ncomp
    real(c_double),        value :: time
    type(c_ptr),           value :: x, u, ctx
    real(c_double), pointer      :: x_values(:), u_values(:)

    call c_f_pointer( x, x_values, [int( dim)])
    call c_f_pointer( u, u_values, [int( ncomp)])
    u_values( 1) = sin( acos( -1._c_double) * (x_values( 1) - poisson_domain( 1)) / poisson_domain( 3)) * &
                    sin( acos( -1._c_double) * (x_values( 2) - poisson_domain( 2)) / poisson_domain( 4))
    poisson_exact_solution = 0_c_int

  end function poisson_exact_solution

  integer(c_int) function poisson_dirichlet_boundary( dim, time, x, ncomp, u, ctx) bind(C)

    integer(c_intptr_t),   value :: dim, ncomp
    real(c_double),        value :: time
    type(c_ptr),           value :: x, u, ctx
    real(c_double), pointer      :: u_values(:)

    call c_f_pointer( u, u_values, [int( ncomp)])
    u_values( 1) = 0._c_double
    poisson_dirichlet_boundary = 0_c_int

  end function poisson_dirichlet_boundary

  subroutine poisson_f0( dim, nf, nfaux, uoff, uoff_x, u, u_t, u_x, aoff, aoff_x, a, a_t, a_x, &
    time, x, nconstants, constants, f0) bind(C)

    integer(c_intptr_t), value :: dim, nf, nfaux, nconstants
    type(c_ptr),    value :: uoff, uoff_x, u, u_t, u_x, aoff, aoff_x, a, a_t, a_x, x, constants, f0
    real(c_double), value :: time
    real(c_double), pointer :: x_values(:), constants_values(:), f0_values(:)
    real(c_double)          :: pi, solution

    call c_f_pointer( x, x_values, [int( dim)])
    call c_f_pointer( constants, constants_values, [int( nconstants)])
    call c_f_pointer( f0, f0_values, [1])
    pi = acos( -1._c_double)
    solution = sin( pi * (x_values( 1) - constants_values( 1)) / constants_values( 3)) * &
               sin( pi * (x_values( 2) - constants_values( 2)) / constants_values( 4))

    ! PETSc assembles f0 * v + f1 . grad(v), hence f0 is minus the RHS.
    f0_values( 1) = -pi**2 * (1._c_double / constants_values( 3)**2 + &
      1._c_double / constants_values( 4)**2) * solution

  end subroutine poisson_f0

  subroutine poisson_f1( dim, nf, nfaux, uoff, uoff_x, u, u_t, u_x, aoff, aoff_x, a, a_t, a_x, &
    time, x, nconstants, constants, f1) bind(C)

    integer(c_intptr_t), value :: dim, nf, nfaux, nconstants
    type(c_ptr),    value :: uoff, uoff_x, u, u_t, u_x, aoff, aoff_x, a, a_t, a_x, x, constants, f1
    real(c_double), value :: time
    real(c_double), pointer :: u_x_values(:), f1_values(:)

    call c_f_pointer( u_x, u_x_values, [int( dim)])
    call c_f_pointer( f1, f1_values, [int( dim)])
    f1_values( 1:dim) = u_x_values( 1:dim)

  end subroutine poisson_f1

  subroutine poisson_g3( dim, nf, nfaux, uoff, uoff_x, u, u_t, u_x, aoff, aoff_x, a, a_t, a_x, &
    time, u_tshift, x, nconstants, constants, g3) bind(C)

    integer(c_intptr_t), value :: dim, nf, nfaux, nconstants
    type(c_ptr),    value :: uoff, uoff_x, u, u_t, u_x, aoff, aoff_x, a, a_t, a_x, x, constants
    real(c_double), value :: time, u_tshift
    real(c_double), intent(out) :: g3(*)
    integer                 :: d

    poisson_g3_call_count = poisson_g3_call_count + 1
    g3( 1:dim*dim) = 0._c_double
    do d = 1, dim
      g3( (d - 1) * dim + d) = 1._c_double
    end do

  end subroutine poisson_g3

end module ct_PETSc_SNES_Poisson