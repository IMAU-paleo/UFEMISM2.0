# SSA FD + PETSc SNES solver — Tier 1 implementation plan

## Goal

Add a **new, standalone** momentum-balance solver that keeps the existing
finite-difference SSA discretisation but replaces the hand-rolled viscosity
(Picard) iteration in
[`momentum_balance_solver_SSA_run`](src/UFEMISM/ice_dynamics/momentum_balance/SSA_DIVA/momentum_balance_solver_SSA.f90#L212)
with a PETSc `SNES` non-linear solve.

"Tier 1" means **defect-correction Newton** (a.k.a. Picard-as-Newton): the
non-linear residual is `F(u) = A(u) u - b(u)`, where `A(u)` and `b(u)` are
*exactly* the linear system the current code already assembles per viscosity
iteration. The Jacobian handed to `SNES` is that same `A(u)` — no analytic
tangent, no new derivation. The payoff is SNES's line search and convergence
control replacing the manual relaxation / divergence back-off logic, plus
optional inexact-Newton (Eisenstat–Walker) tolerancing on the inner solve.

Expected outcome: same solution as the `SSA` solver (identical discretisation
and fixed point) to solver tolerance, with fewer / more robust outer iterations
and no `visc_it_relax` / `Glens_flow_law_epsilon_sq_0` babysitting.

## Progress

| Step | Status | Notes |
| --- | --- | --- |
| Prereq 1 — extract assembly | done | `assemble_SSA_DIVA_linearised_matrix_eq` on `atype_momentum_balance_solver_SSADIVA` (commit "Extract SSA/DIVA matrix assembly") |
| Prereq 2 — expose SSA helpers | done | `calc_effective_viscosity`, `calc_applied_basal_friction_coefficient`, `calc_vertically_averaged_flow_parameter`, `initialise_SSA_velocities_from_file` are now `public` |
| Step 0 — config plumbing | done | `SSA_FD_SNES_snes_rtol` (1E-6), `_snes_abstol` (1E-4), `_snes_maxits` (50), `_use_EW` (`.true.`) added to `model_configuration_type_and_namelist.f90` (declaration + type + namelist group + assignment). Inner KSP/PC reuse `stress_balance_PETSc_*`. Not yet consumed (Step 5). |
| Step 1 — new module + registration | done | `momentum_balance_solver_SSA_FD_SNES.f90` (type `extends type_momentum_balance_solver_SSA`); `run` is a `crash` stub; registered in `create_momentum_balance_solver`; 46 output-file `select case` lists gained `'SSA_FD_SNES'` in the b-grid group. `allocate`/`deallocate`/`initialise`/`remap` are inherited unchanged (no PETSc state yet — deferred to Step 6). |
| Step 2 — `update_SSA_coefficients_from_velocity` | done | private method on the new type: `calc_horizontal_strain_rates` + `calc_effective_viscosity` (fixed `eps0`) + `calc_applied_basal_friction_coefficient`; no `apply_velocity_limits` / `relax_viscosity_iterations`. Not called yet (the `run` stub still crashes). |
| Step 3 — residual callback | done (compiles; runtime-untested) | `SSA_FD_SNES_form_function`, `bind(C)`: unpack `x` → `u_vav_b`/`v_vav_b`, `update_SSA_coefficients_from_velocity`, `assemble_SSA_DIVA_linearised_matrix_eq`, `f = A·x − b` via `multiply_PETSc_matrix_with_vector_1D` + `VecCopy`. |
| Step 4 — Jacobian callback | done (compiles; runtime-untested) | `SSA_FD_SNES_form_jacobian`, `bind(C)`: same assembly, then `mat_CSR2petsc` → `MatCopy(…, DIFFERENT_NONZERO_PATTERN)` into the registered `A_petsc`. `pmat` == `amat`. The "assemble once per SNES step + cache" optimisation is **not** done. |
| Step 5 — SNES driver | done (compiles; runtime-untested) | `run` is now thin (early-out + BC prep into `self` components) and calls `solve_SSA_FD_SNES` (separate routine so its dummies can be `target`). Prime assembly → `A_petsc`/`sol`, `SNESCreate`, register callbacks, `SNESNEWTONLS`, KSP/PC from `stress_balance_PETSc_*`, EW, `SNESSolve`, disentangle, `apply_velocity_limits`, destroy. `crash` on `SNESConvergedReason < 0`. |
| Step 6 — persistent state & lifecycle | simplified — not needed | PETSc objects (`snes`, `A_petsc`, `sol`, `res_vec`) are created and destroyed **within each solve**, so `allocate`/`deallocate`/`initialise`/`remap` are **not** overridden. Only new non-PETSc state: the `p_ice`/`p_geom`/`p_bed_roughness` context pointers and `BC_prescr_*_applied` arrays. Making the objects persistent (rebuilt on remap) is a later optimisation. |
| Step 7 — PETSc `bind(C)` interfaces | done | `bind(C)` interfaces for `SNESSetFunction`, `SNESSetJacobian`, `SNESKSPSetUseEW`, `SNESGetConvergedReason` (link-resolved against libpetsc). Callbacks are `bind(C)` functions taking raw `c_intptr_t` handles, wrapped into `tVec`/`tMat` via `%v`. `self` reaches the callbacks through the module pointer `SSA_FD_SNES_active_solver` (a polymorphic `self` can't go through `c_loc`), set around `SNESSolve`. |

## Run findings (MISMIP_mod, `integrated_test_SSA_notime_MISMIP_mod_full`, 2 ranks)

1. **Interop is correct.** `||F(u0)||` computed directly from Fortran arrays and
   via PETSc calling the `bind(C)` residual callback match to all digits
   (1.96025E+06). Raw-handle marshalling, `c_funloc` callbacks and the
   vector/matrix partitioning all work.
2. **Tier 1 (Picard operator as the Jacobian) does not converge.** It is
   mathematically undamped Picard: with `NEWTONLS` + `bt` it fails the line
   search after ~30 iterations; with `basic` line search + damping 0.2 (i.e.
   relaxed Picard) the residual falls ~5 orders over ~370 iterations and then
   **limit-cycles**. This is the same behaviour the hand-rolled viscosity
   iteration shows without its adaptive relaxation / `eps0` inflation.
3. **Non-dimensionalisation is required** and is now implemented (the SNES
   unknown is `u_hat = u / velocity_scale`, the residual is `f_hat = (A u - b) /
   stress_scale`, `velocity_scale = 1e3`, `stress_scale = 1e5`, same as the FEM
   solver). This brings `||F_hat(u0)||` to ~19.6 and makes the residual decrease
   monotonically.
4. **Tier 2 (JFNK) converges only at Picard rate and then stalls.** The
   assembled Picard operator is an excellent *preconditioner* (the Krylov solve
   needs ~1 iteration when the KSP tolerance is loose), so each SNES step is
   essentially a Picard step: `||F_hat||` falls 19.6 -> ~1 over ~15 iterations,
   then the line search fails (~1.3 orders of reduction). Trying to resolve the
   matrix-free true Jacobian accurately instead (Eisenstat-Walker off, tight
   `ksp_rtol`, GMRES restart 200, exact block-LU preconditioner) **fails**: GMRES
   cannot drive the linear residual below ~1e-4 in 200+ iterations. The
   matrix-free `J*v` is too noisy - the residual re-assembles the viscosity /
   friction coefficients and applies mesh operators each evaluation, and the
   `eta` and `max(0.1,H)` clamps make it non-smooth, so finite-difference
   Jacobian-vector products have only a few correct digits.

**Conclusion.** JFNK is not viable here because the FD residual is not smooth /
clean enough for matrix-free differencing. Real Newton convergence needs an
**analytic Jacobian** (Tier 3) - at minimum the viscosity shear-thinning term
`d(N)/d(u)` - which is exactly what the FEM solver assembles and why it
converges. **Tier 3 is the chosen direction.**

## Tier 3 status — WORKING

Implemented and validated on MISMIP_mod (`integrated_test_SSA_notime_MISMIP_mod_full`,
solver `SSA_FD_SNES`, 2 ranks).

- **Analytic Jacobian** `dF_hat/du_hat` = the frozen-coefficient ("Picard")
  operator `A(u)` **plus** the coefficient-derivative term `d/du[A(u)] u`:
  - the Glen shear-thinning term `dN/du` (chain `u_b -> strain rates on a -> eta
    on a -> N on a and b`);
  - the sliding-law term `d beta_b/du` (Zoet-Iverson only; the assembly `crash`es
    for other laws).
  Assembled in `assemble_SSA_coeff_jacobian_CSR`, added to the Picard operator
  with `MatAXPY` in `build_SSA_FD_SNES_jacobian_petsc`, then scaled by
  `velocity_scale/stress_scale`. Full derivation + assembly notes in
  `SSA_FD_SNES_jacobian_derivation.tex` (a standalone LaTeX document; build with
  `pdflatex SSA_FD_SNES_jacobian_derivation.tex`, run twice for the
  cross-references and table of contents).
- **Jacobian verified**: central finite-difference check `J v` vs
  `(F(u+hv)-F(u-hv))/2h` at a non-zero base velocity gives rel. error ~8e-6.
- **Convergence**: Newton is quadratic - e.g. `||F_hat||` 12.5 -> 8.6 -> 4.8 ->
  2.0 -> 0.57 -> 0.095 -> 0.005, then it bottoms out at the ~1e-3 noise floor of
  the non-smooth `eta`/friction/`max(0.1,H)` clamps, where the line search
  reports failure (reason -6). The solve converges (`SNESConvergedReason = 3`) in
  ~7-14 Newton iterations; where it stops on -6 instead, the solution is accepted
  if `||F_hat|| < SSA_FD_SNES_resid_floor` (1e-1).
- **Solution check**: agrees with the FD Picard solver run to convergence
  (`visc_it_nit = 2000`; it does *not* converge in the default 50) to ~4-5% RMS
  on the velocity, with larger localised differences near the grounding line -
  consistent with both hitting the clamp noise floor. The `SSA_FEM_PETSc` solver
  is not a usable reference on this config (natural BCs only -> a different
  problem, ~120% different).

### Supporting machinery added

- **Direct LU inner solve** (`KSPPREONLY` + `PCLU`): the analytic Jacobian is too
  stiff for the FD solver's gmres+bjacobi, which returns poor Newton directions.
  Matches the FEM solver's default.
- **Picard warm start**: 5 heavily-relaxed Picard iterations before the SNES
  solve, to move off `u = 0` where the velocity-weakening sliding law is nearly
  non-differentiable and Newton cannot start.
- **`gather_CSR_to_all`**: `M_ddx_b_a` / `M_ddy_b_a` / `M_map_b_a` are broadcast
  to full local copies once per solve so the 2-hop Jacobian stencil can read
  operator rows for vertices this rank does not own.
- **Config** `SSA_FD_SNES_snes_rtol` / `_abstol` defaults changed to `1e-3`
  (the reachable range given the noise floor); `_use_EW` is now unused (direct
  solve).

### Follow-ups (not blocking)

- Other sliding laws (`Weertman`, `Coulomb`, `Budd`, ...) need their own
  `d beta_a/d|u|` in `assemble_SSA_coeff_jacobian_CSR`.
- Investigate the ~5% grounding-line discrepancy vs converged Picard (likely the
  clamp noise floor, but worth confirming with a smoother test case).
- Persistent PETSc objects / operator gathers across solves (rebuilt every solve
  now); a config knob for the warm-start iteration count and the LU vs iterative
  choice.
- DIVA (deferred - see the scope note).

**Scope note.** SSA only for now. The residual / Jacobian / driver could be
hoisted to `atype_momentum_balance_solver_SSADIVA` and parameterised by two
deferred hooks (`recompute_nonlinear_coefficients`,
`reconstruct_after_nonlinear_solve`) to serve the DIVA as well — the linear
system, the unknown vector and the boundary conditions are already shared. This
is deliberately deferred to a later iteration; until then the SNES machinery
lives on the concrete `type_momentum_balance_solver_SSA_FD_SNES`.

## Design decisions

- **Separate solver class.** New type
  `type_momentum_balance_solver_SSA_FD_SNES`, selected with
  `choice_stress_balance_approximation = 'SSA_FD_SNES'`. The existing `SSA`,
  `DIVA`, `SSA_FEM_PETSc` solvers are untouched. (Name is provisional — rename
  freely.)
- **Inherit the FD machinery.** Extend
  [`type_momentum_balance_solver_SSA`](src/UFEMISM/ice_dynamics/momentum_balance/SSA_DIVA/momentum_balance_solver_SSA.f90#L34)
  so we reuse `calc_driving_stress`, `calc_horizontal_strain_rates`,
  `calc_effective_viscosity`, `calc_applied_basal_friction_coefficient`,
  `calc_vertically_averaged_flow_parameter`, the shared field allocation, and
  remap. Only `run_momentum_balance_solver` and the solver name are genuinely
  new; `allocate` / `deallocate` / `initialise` / `remap` are thin overrides
  that call the parent and then manage the persistent PETSc objects.
- **Reuse the existing assembly.** The per-row assembly in
  [`solve_linearised_SSA_DIVA_infinite_slab.f90`](src/UFEMISM/ice_dynamics/momentum_balance/SSA_DIVA/solve_linearised_SSA_DIVA_infinite_slab.f90)
  (`calc_SSA_DIVA_stiffness_matrix_row_free`,
  `calc_SSA_DIVA_sans_stiffness_matrix_row_free`,
  `calc_SSA_DIVA_stiffness_matrix_row_BC`) is factored so both the legacy
  `solve_SSA_DIVA_linearised` and the new residual/Jacobian callbacks build the
  identical `A_CSR` / `bb`.
- **Reuse existing PETSc plumbing.**
  [`petsc_basic`](src/UPSY/basic/petsc/petsc_basic.f90) already provides
  `mat_CSR2petsc`, `vec_double2petsc`, `vec_petsc2double`,
  `multiply_PETSc_matrix_with_vector_1D`. The inner KSP/PC reuse the existing
  `C%stress_balance_PETSc_KSPtype` / `_PCtype` / `_rtol` / `_abstol` config.

## Background maths

The FD SSA fixed point is: assemble `A(u_k)` and `b(u_k)` from the current
velocity (via `N = eta(u_k) H` on the a-grid and `beta_b(u_k)` on the b-grid),
solve `A(u_k) u_{k+1} = b(u_k)`, relax, repeat until `u` stops changing. The
converged `u*` satisfies

```
F(u*) := A(u*) u* - b(u*) = 0
```

so `F` is a ready-made non-linear residual. Feeding `SNES`:

- `FormFunction(u)  -> F(u) = A(u) u - b(u)`
- `FormJacobian(u)  -> J = A(u)`   (the "Picard operator"; not the true tangent)

gives `SNESNEWTONLS` doing damped Picard. Because `J` omits `dN/du` and
`dbeta_b/du`, convergence is linear, not quadratic — but the line search makes
it robust without tuning, and each SNES step costs exactly one linear solve,
i.e. the same as one of today's viscosity iterations.

Nice property: the Dirichlet / `zero` / `infinite` BC rows the current code
writes (identity or Laplacian-style row + matching `bb`) make `F_i` the correct
residual for those rows automatically — `F_i = (A u - b)_i` already encodes
`u_i - u_prescribed` or `du/dn = 0`.

## Prerequisite refactor (no behaviour change) — DONE

1. **Extract the linear-system assembly.** Implemented as a public method on
   `atype_momentum_balance_solver_SSADIVA`:

   ```fortran
   subroutine assemble_SSA_DIVA_linearised_matrix_eq( self, u_ii_term, &
       BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b, A_CSR, bb, uv_buv)
   ```

   It holds the body of `solve_SSA_DIVA_linearised` from "Initialise the
   stiffness matrix" through `call A_CSR%finalise()` (the `do row_tiuv` loop
   and its BC / free / sans branches), and also returns `uv_buv` — the current
   `[u_vav_b, v_vav_b]` packed into the `2*nTri` vector layout, i.e. the SNES
   initial guess. `A_CSR` / `bb` / `uv_buv` are `allocatable, intent(inout)` and
   are allocated inside. `solve_SSA_DIVA_linearised` is now `gather prev` →
   `assemble_SSA_DIVA_linearised_matrix_eq` → `solve_matrix_equation_CSR_PETSc`
   → disentangle.

2. **Expose parent helpers.** `calc_effective_viscosity`,
   `calc_applied_basal_friction_coefficient`,
   `calc_vertically_averaged_flow_parameter` (and
   `initialise_SSA_velocities_from_file`) on `type_momentum_balance_solver_SSA`
   are now `procedure, public`. Pure visibility change.

## Implementation steps

### Step 0 — Config plumbing — DONE

In
[`model_configuration_type_and_namelist.f90`](src/UPSY/basic/model_configuration/model_configuration_type_and_namelist.f90),
directly after the `SSA_FEM_PETSc_*` block:

- `choice_stress_balance_approximation` doc comment lists `'SSA_FD_SNES'` (done
  in Step 1).
- New `_config` fields + type fields + namelist entries + assignment for:
  - `SSA_FD_SNES_snes_rtol`   (default `1E-6_dp`) — relative residual reduction
  - `SSA_FD_SNES_snes_abstol` (default `1E-4_dp`) — on the raw (dimensional)
    residual norm, so rtol drives convergence in practice
  - `SSA_FD_SNES_snes_maxits` (default `50`)
  - `SSA_FD_SNES_use_EW`      (logical, default `.true.`) — Eisenstat–Walker
    inexact-Newton on the inner KSP
- The inner KSP/PC reuse `stress_balance_PETSc_KSPtype` / `_PCtype` /
  `_rtol` / `_abstol` — no new fields there.
- Existing `.cfg` files need no change: unspecified namelist variables take the
  `_config` defaults.

### Step 1 — New module + registration — DONE

- `src/UFEMISM/ice_dynamics/momentum_balance/SSA_DIVA/momentum_balance_solver_SSA_FD_SNES.f90`
  (keep it next to the other SSA_DIVA solvers).
- `type_momentum_balance_solver_SSA_FD_SNES`, `extends(type_momentum_balance_solver_SSA)`.
- New/overridden bindings:
  - `get_momentum_balance_solver_name` → `'SSA_FD_SNES'`
  - `run_momentum_balance_solver` → `..._SSA_FD_SNES_run` (Step 5)
  - `allocate_momentum_balance_solver` / `deallocate...` / `initialise...` /
    `remap...` → call the parent method, then create / destroy / rebuild the
    persistent PETSc objects and the `n2tiuv` packing (Step 6)
  - `create_restart_file_old` / `write_to_restart_file_old` → inherit the
    parent's (`u_vav_b` / `v_vav_b` only), or override the filename string.
- Register in
  [`create_momentum_balance_solver`](src/UFEMISM/ice_dynamics/momentum_balance/momentum_balance_solver_main.f90#L36):
  `case ('SSA_FD_SNES') ; allocate( type_momentum_balance_solver_SSA_FD_SNES :: momentum_balance_solver)`
  plus the `use` line.
- The build uses `GLOB_RECURSE` (`src/UFEMISM/CMakeLists.txt`), so the new file
  is picked up on the next `cmake` configure — no source-list edit needed.
- Also updated: the 46 `select case (C%choice_stress_balance_approximation)`
  branches in `ice_velocity_model_basic.f90`, `grid_output_files.f90` and
  `mesh_output_files.f90` — `'SSA_FD_SNES'` joins the b-grid group
  (`'none','SIA','SSA', ...`), since the solver exposes velocities on the
  triangles exactly like `SSA`.
- `allocate` / `deallocate` / `initialise` / `remap` are **not** overridden yet;
  the type inherits them from `type_momentum_balance_solver_SSA`. They become
  thin overrides in Step 6 once there is persistent PETSc state to manage.
- `run_momentum_balance_solver` is a `crash` stub until Step 5.

### Step 2 — "recompute coefficients" helper — DONE

Private method `update_SSA_coefficients_from_velocity( self, ice, geom, bed_roughness)`
on `type_momentum_balance_solver_SSA_FD_SNES`, running one pass of the
coefficient chain at the current `self%u_vav_b` / `self%v_vav_b`:

```
call self%calc_horizontal_strain_rates()
call self%calc_effective_viscosity( ice, geom, C%Glens_flow_law_epsilon_sq_0)   ! FIXED eps0, no adaptive inflation
call self%calc_applied_basal_friction_coefficient( ice, geom, bed_roughness)
```

It will be called at the top of both the residual and the Jacobian callback
(Steps 3–4); it is currently unreferenced because the `run` stub still crashes.
Note the eta `min/max` clamp inside `calc_effective_viscosity` stays — it is a
coefficient bound, and the Picard Jacobian ignores its derivative anyway, so it
is consistent. `apply_velocity_limits` and `relax_viscosity_iterations` are
**not** called here — SNES owns the iterate.

### Step 3 — Residual callback `FormFunction`

C-interoperable module procedure with the PETSc SNES residual signature
`(SNES, Vec x, Vec f, void *ctx, PetscErrorCode *ierr)`:

1. Recover `self` from `ctx` (`c_f_pointer`; pass `c_loc(self)` at
   `SNESSetFunction` time — same pattern as the FEM solver's context handling).
2. `vec_petsc2double(x, uv)` then unpack into `self%u_vav_b` / `self%v_vav_b`
   using `self%mesh%n2tiuv` (identical ordering to `solve_SSA_DIVA_linearised`
   lines 135–145).
3. `call self%update_SSA_coefficients_from_velocity( ice, geom, bed_roughness)`.
4. `call self%assemble_SSA_DIVA_linearised_matrix_eq( self%basal_friction_coefficient_b, mask, u_bc, v_bc, A_CSR, bb, uv_buv)`.
5. `call mat_CSR2petsc(A_CSR, A_petsc)` (or reuse the cached Mat from the
   Jacobian call — see Step 4 optimisation).
6. `F = A_petsc * x - b`: `multiply_PETSc_matrix_with_vector_1D(A_petsc, x, f)`
   then `VecAXPY(f, -1, b_petsc)` (with `b_petsc = vec_double2petsc(bb)`).
7. Return `ierr = 0`.

`ice`, `geom`, `bed_roughness` are not SNES arguments — stash pointers to them
on `self` for the duration of the solve (set in the `run` routine, nullify
after).

### Step 4 — Jacobian callback `FormJacobian`

Signature `(SNES, Vec x, Mat Amat, Mat Pmat, void *ctx, PetscErrorCode *ierr)`:

1. Recover `self`, unpack `x`, `update_SSA_coefficients_from_velocity`,
   `assemble_SSA_DIVA_linearised_matrix_eq` — same as Step 3 items 1–4.
2. Fill `Amat` (and `Pmat`, same matrix) from `A_CSR`. Either `mat_CSR2petsc`
   into a fresh Mat and `MatCopy`, or preallocate a persistent `self%A_petsc`
   once (fixed sparsity from the mesh stencil) and overwrite its values each
   call via `MatSetValues` / `MatZeroEntries` + refill.
3. `MatAssemblyBegin/End`.

**Optimisation (recommended):** SNES calls `FormJacobian` then `FormFunction`
at the same `x` within a step. Assemble `A_CSR` + `bb` **once** per SNES
iteration — cache them (and the PETSc Mat/Vec) on `self`, keyed by a solve-local
counter — and have both callbacks reuse the cache. This halves the assembly
cost and guarantees `F` and `J` are consistent.

### Step 5 — `run` routine (SNES driver)

Replaces the viscosity-iteration `do while` loop. Structure:

```
call init_routine(...)

! early-out: no grounded ice / no sliding  -> u = v = 0, return   (copy from SSA_run lines 247-255)
! handle optional BC_prescr_* args                                (copy from SSA_run lines 257-273)
call self%calc_driving_stress( geom)

! stash context pointers for the callbacks
self%p_ice => ice ; self%p_geom => geom ; self%p_bed_roughness => bed_roughness
self%BC_prescr_mask_b_applied = ... etc.

! pack current velocity guess into the PETSc solution vec
do ti ... ; sol(n2tiuv) = u_vav_b / v_vav_b ; end do ; vec_double2petsc(sol, self%sol)

! (re)build SNES if mesh changed
SNESSetFunction( snes, self%res_vec, FormFunction, c_loc(self))
SNESSetJacobian( snes, self%A_petsc, self%A_petsc, FormJacobian, c_loc(self))
SNESSetType( snes, SNESNEWTONLS)
SNESSetTolerances( snes, abstol, rtol, PETSC_DEFAULT, maxits, PETSC_DEFAULT)
SNESGetKSP( snes, ksp)
  KSPSetType( ksp, C%stress_balance_PETSc_KSPtype)
  KSPGetPC( ksp, pc) ; PCSetType( pc, C%stress_balance_PETSc_PCtype)
  KSPSetTolerances( ksp, C%stress_balance_PETSc_rtol, C%stress_balance_PETSc_abstol, ...)
if (C%SSA_FD_SNES_use_EW) SNESKSPSetUseEW( snes, PETSC_TRUE)

SNESSolve( snes, PETSC_NULL_VEC, self%sol)
SNESGetConvergedReason( snes, reason) ; if (reason < 0) call crash/warning
SNESGetIterationNumber( snes, self%n_visc_its)
SNESGetLinearSolveIterations( snes, self%n_Axb_its)

! unpack self%sol -> u_vav_b / v_vav_b
call self%apply_velocity_limits()          ! post-solve clamp only
self%p_ice => null() ; ...

call finalise_routine(...)
```

Reuse `SNESNEWTONLS`, `SNESGetIterationNumber`, `SNESGetLinearSolveIterations`,
`SNESSetTolerances`, `SNESGetKSP` — already imported by the FEM solver from the
`petsc` module, so the import list is known-good.

### Step 6 — Persistent state & lifecycle

Add to the type:

- `type(tSNES) :: snes`
- `type(tMat)  :: A_petsc`  (persistent, fixed sparsity)
- `type(tVec)  :: sol, res_vec, b_petsc`
- context pointers: `p_ice`, `p_geom`, `p_bed_roughness` (class pointers),
  `BC_prescr_mask_b_applied` / `_u_` / `_v_` arrays
- cached `A_CSR` (`type_CSR_matrix_dp`) + `bb` + the solve-local assembly
  counter for the Step 4 optimisation

Lifecycle:

- `allocate` / `initialise`: call parent, then `SNESCreate` and create the
  Vecs (sizes `2*nTri_loc`). Set `self%PETSc_rtol` / `_abstol` as the parent
  does (still used by the fallback path, harmless otherwise).
- `deallocate`: `SNESDestroy`, `MatDestroy`, `VecDestroy`, then parent.
- `remap`: call parent `remap`, then destroy and recreate all PETSc objects
  (sizes and sparsity change with the mesh). Simplest is a
  `destroy_petsc_objects` / `create_petsc_objects` pair called from
  `initialise`, `remap`, `deallocate`.

### Step 7 — PETSc `bind(C)` interfaces (only if missing)

The FEM solver imports `SNES*` names from the `petsc` module but reaches
`SNESSetJacobian` through a manual `bind(C)` interface
([`snes_set_jacobian`, line 268](src/UFEMISM/ice_dynamics/momentum_balance/SSA_FEM_PETSc/momentum_balance_solver_SSA_FEM_PETSc.f90#L268)).
Check whether the installed PETSc Fortran module exposes usable
`SNESSetFunction` / `SNESSetJacobian` / `SNESKSPSetUseEW` wrappers that accept
a Fortran callback + context. If not, add `bind(C)` interfaces modelled on the
existing `snes_set_jacobian` one, and make the callbacks `bind(C)` with
`type(c_ptr), value :: ctx`.

## Testing & validation

1. **Residual sanity check (unit).** Run the `SSA` solver to convergence on a
   small mesh, copy its `u_vav_b` / `v_vav_b` into the new solver, call
   `FormFunction` once — `‖F‖` should be at solver-tolerance zero.
2. **Solution equivalence (integrated).** Add
   `automated_testing/UFEMISM/integrated_test_SSA_notime_MISMIP_mod_full/config_SSA_FD_SNES.cfg`
   (copy `config_SSA_FEM_PETSc.cfg`, set `choice_stress_balance_approximation = 'SSA_FD_SNES'`).
   The velocity field must match the `SSA` solver's to within KSP tolerance.
3. **MISMIP full.** Run the
   `integrated_test_MISMIP_mod_full` configs with the new solver; compare
   grounding-line trajectory against the `SSA` baseline.
4. **Iteration counts.** Log `n_visc_its` (SNES iterations) and `n_Axb_its`
   (total KSP iterations) vs. the `SSA` solver on the same runs — this is the
   headline result: fewer / steadier outer iterations, and whether EW cuts
   total KSP work.
5. **Regression.** Confirm `SSA` and `DIVA` integrated tests are unchanged. The
   only shared code touched is the prerequisite assembly extraction and the
   additive `'SSA_FD_SNES'` `case` labels in the three output-file dispatchers —
   none of which alters the `SSA` / `DIVA` code paths.

## Tuning levers (once it runs)

- `SNESKSPSetUseEW` — inexact Newton; biggest expected KSP-work saving.
- Line-search type: default `bt` (cubic backtracking); try `l2` / `cp` if it
  stalls.
- `do_include_SSADIVA_crossterms` — the `sans` assembly (divide through by `N`)
  is better conditioned; still valid as a residual.
- Inner KSP tolerance floor — with EW on, `stress_balance_PETSc_rtol` becomes
  the loosest allowed, not the target.
- Non-dimensional scaling of `F`: if `snes_rtol` behaves oddly, row-scale the
  residual by `1/|diag(A)|`. The FEM solver non-dimensionalises heavily; for
  Tier 1 the exact linear preconditioner usually makes this unnecessary, but
  keep it in reserve.

## Risks & non-goals

- **Idealised BCs.** `periodic_ISMIP-HOM` and `infinite_SSA_icestream` rows
  reference `self%u_vav_b_prev` (a lagged copy). Inside SNES that term must
  reference the current iterate instead. For the target cases (MISMIP,
  realistic domains) only `zero` / `infinite` are used and are fine — restrict
  the new solver to those initially, or refactor those two BC branches to read
  the current `x`.
- **Not Tier 2/3.** No matrix-free Jacobian, no analytic `dN/du` /
  `dbeta_b/du`. Convergence stays linear. If profiling later shows the Picard
  operator is the bottleneck in outer iterations, that is the Tier 2/3 follow-up
  (`-snes_mf_operator` with this `A(u)` as preconditioner, then the analytic
  tangent).
- **Non-smooth clamps.** Keep the eta `min/max` and `vel_max` clamps out of the
  Newton loop's residual derivative (they already are — eta clamp is a
  coefficient, `apply_velocity_limits` moves to post-solve only).
- **Assembly cost.** A fresh `A_CSR` per SNES iteration is the same allocation
  pattern as today; the Step 4 cache removes the duplicate assembly within a
  step. Persistent-sparsity `MatSetValues` refill is a later optimisation.
- **Parallel.** `A_CSR` is already row-distributed (`i1:i2`); `mat_CSR2petsc`
  and SNES run on `PETSC_COMM_WORLD` over the full mesh — no sub-communicator /
  sub-mesh machinery like the FEM solver needs.

## File checklist

| File | Change |
| --- | --- |
| `src/UPSY/basic/model_configuration/model_configuration_type_and_namelist.f90` | new `SSA_FD_SNES_*` config fields + namelist + assignment; extend `choice_stress_balance_approximation` doc — **done** |
| `src/UFEMISM/ice_dynamics/momentum_balance/SSA_DIVA/momentum_balance_solver_SSADIVA.f90` | new public `assemble_SSA_DIVA_linearised_matrix_eq` (extracted) — **done** |
| `src/UFEMISM/ice_dynamics/momentum_balance/SSA_DIVA/solve_linearised_SSA_DIVA_infinite_slab.f90` | split assembly out of `solve_SSA_DIVA_linearised` |
| `src/UFEMISM/ice_dynamics/momentum_balance/SSA_DIVA/momentum_balance_solver_SSA.f90` | make three coefficient helpers `public` |
| `src/UFEMISM/ice_dynamics/momentum_balance/SSA_DIVA/momentum_balance_solver_SSA_FD_SNES.f90` | **new** solver module (type, run, residual/Jacobian callbacks, lifecycle) |
| `src/UFEMISM/ice_dynamics/momentum_balance/momentum_balance_solver_main.f90` | register `'SSA_FD_SNES'` |
| build files (CMake/Make) for the momentum-balance dir | add the new source |
| `automated_testing/UFEMISM/integrated_test_SSA_notime_MISMIP_mod_full/config_SSA_FD_SNES.cfg` | **new** test config |
