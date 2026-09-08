# Implementation plan: SSA via PetscFE + PetscSNES

A new momentum-balance solver that discretises and solves the Shallow Shelf
Approximation entirely with PETSc (DMPlex + PetscFE for the discretisation,
PetscSNES for the nonlinear solve), added as one more choice for
`choice_stress_balance_approximation` **without touching any existing solver**.

The solver is selected with `choice_stress_balance_approximation = 'SSA_FEM_PETSc'`
(the `FEM` distinguishes it from the finite-difference-based SSA/DIVA solvers).

**The mathematical derivation of the PETSc residual/Jacobian callbacks
(`f0`, `f1`, `g0`, `g3`) from the SSA PDEs is in
[`SSA_FEM_PETSc_weak_form_derivation.md`](SSA_FEM_PETSc_weak_form_derivation.md).**

## Progress

| Phase | Status | Notes |
| --- | --- | --- |
| 0 - Scaffolding | **done** | `type_momentum_balance_solver_SSA_FEM_PETSc` in `src/UFEMISM/ice_dynamics/momentum_balance/SSA_FEM_PETSc/`; dispatch case + config comment added; integrated test `integrated_test_SSA_notime_MISMIP_mod_full` now has `config_SSA.cfg` + `config_SSA_FEM_PETSc.cfg` and its `test_script.csh` runs both solvers. |
| 1 - DMPlex + PetscFE + SNES skeleton, constant-coefficient linear SSA | **done** | Full `DMPlex -> PetscFE (P1, 2-comp) -> PetscDS (f0/f1 + analytic g0/g3) -> SNES` pipeline. On the integrated-test mesh, 2 MPI ranks: SNES converges in 1 iteration and the nodal field matches the closed-form uniform solution `-tau/beta` to `~1e-16`. `overlap = 0` assembly confirmed correct in parallel for this case. |
| 2 - Auxiliary fields (real spatially varying coefficients) | **done** | 4-component P1 aux field on a `DMClone`d DM, attached with `DMSetAuxiliaryVec`; `f0/f1/g0/g3` read `a[]`. Per-vertex coefficients from the existing UFEMISM machinery, incl. the sub-grid grounded-fraction scaling of `beta` (friction vanishes under floating ice). Superseded by Phase 3's residual (viscosity is no longer frozen). |
| 3 - Fully non-linear Newton residual (viscosity **and** friction), analytic Jacobian, no outer loop | **done** | Glen `eta(grad u)` **and** the Zoet-Iverson `beta(|u|)` both evaluated pointwise in the residual, with analytic Jacobians (`g3` = frozen membrane + rank-1 `d eta / d grad u`; `g0` = `beta I` + rank-1 `d beta / d|u|`). **One Newton solve, no Picard loop.** Aux field `[Abar, H, tauc_eff, tau_dx, tau_dy]`; `eps0`, `n` and the ZI params via `PetscDSSetConstants`. Warm-started persistent solution vector. Sign fix: `f0 = beta*u - tau_d`. On the integrated test, 2 ranks, cold start: Newton converges in **10 iterations**; signed `u_vav`/`v_vav` correlate **+0.98** with the FD `SSA` solver. Magnitudes are ~1.4x the FD solver's, which itself does not converge here (still climbing toward the FE result as its Picard count is raised 50 -> 500); the residual gap is FD under-convergence + different strain-rate discretisation + grounding-line `beta` representation (Phase 5). Pointwise sliding relation is `SSA_FEM_PETSc_sliding_beta` (ZI only; TODO to merge into `sliding_laws`). |
| 6 - Scaling, solver, preconditioner (pulled forward) | **done** | Nondimensionalisation (`u_hat = u/velocity_scale`, residual/Jacobian / `stress_scale`) - transparent (identical answer). Config-selectable `pc_type` (`lu`/`gamg`/`bjacobi`) + rigid-body near-null-space; all three converge to the identical answer on the integrated test, with `gamg` needing under half `bjacobi`'s KSP iterations (362 vs 857). New config knobs `SSA_FEM_PETSc_{pc_type,snes_rtol,snes_abstol,snes_maxits}_config`. `lu` stays the default. Done ahead of Phases 4/5 at the repo owner's request. FD-solver nondimensionalisation deferred to a separate PR (Section 9). |
| 5 - Boundary conditions (ice-front back-pressure) | **in progress** | `mesh_to_dmplex_masked` (ice-covered sub-mesh topology) done; the *entire* PETSc build (`self%dm` and everything downstream) now runs on that sub-mesh instead of the whole mesh, rebuilt from scratch every `run()` call (ice mask/margin can move every timestep) including a re-seeded (not just zeroed) `self%sol` warm start. Verified on the integrated test, 2 ranks: Newton still converges in 10 iterations (10 KSP its, `lu`); max speed shifts from 2.120e3 to 1.193e3 m/yr, the expected result of excluding ice-free area from the discretisation, not a regression. **Still to do**: the actual back-pressure term (`DMAddBoundary(DM_BC_NATURAL)` + `PetscDSSetBdResidual`, `f0_bd,i = -tau_o,i`) on the sub-mesh's new exterior boundary - until that lands, the margin is a plain natural (zero-traction) boundary, not yet the physical condition. |
| 4, 7, 8 | not started (Phase 4 analytic Jacobian folded into Phase 3; outer-loop / Picard-option is now moot for the pointwise-friction path) |

Implementation notes that deviate from the original plan:

- **Step-plan order changed 2025-09-08**: Phase 6's nondimensionalisation is being
  done now, before Phases 4 (BCs) and 5 (Picard/robustness options), at the repo
  owner's request.
- **Single module for now.** The `bind(C)` PETSc interfaces, the pointwise
  residual/Jacobian functions and the orchestration all live in
  `momentum_balance_solver_SSA_FEM_PETSc.f90`, following the proven structure of
  `ct_PETSc_SNES_Poisson.f90`. Split into `SSA_FEM_PETSc_weak_form.f90` etc. once
  it grows unwieldy.
- **Integrated-test wiring** was done by the repo owner: two config files and a
  loop in `test_script.csh`, rather than a separate sibling test directory.

## 1. Guiding constraints

- **Additive only.** No changes to `momentum_balance_solver_SSA`,
  `_DIVA`, `_SSADIVA`, `_BPA`, `_SIA`, etc. The only edits to existing files are:
  one `case` in `create_momentum_balance_solver`, one config option, and
  (later) test-harness wiring.
- **Same public contract.** The new class implements the same deferred
  procedures as every other solver (`allocate`/`deallocate`/`initialise`/`run`/
  `set_velocities_to_solver_results`/`remap` + name + restart hooks, see
  `src/UFEMISM/ice_dynamics/momentum_balance/basic/momentum_balance_solver_basic.f90`).
- **Downstream stays on the b-grid.** The solver's *internal* unknown is a
  nodal (vertex, P1) velocity field, but it exposes the result as
  `u_vav_b`/`v_vav_b` on the triangles, exactly like the existing SSA solver.
  Mass continuity, CFL, `calc_secondary_velocities`, output and restart are
  therefore unaffected (this is the low-risk "project nodal velocity to the
  b-grid" option; a native FE transport scheme is out of scope here).
- **Reuse the existing PETSc/DMPlex layer.** `mesh_to_dmplex` and the
  `bind(C)` PetscFE/PetscDS/SNES pattern already exist and work; see
  `src/UPSY/basic/petsc/petsc_dmplex.f90` and
  `src/UPSY/validation/component_tests/PETSc_finite_elements/ct_PETSc_SNES_Poisson.f90`.

## 2. Where it slots into the architecture

| Concern | Existing mechanism | New solver |
| --- | --- | --- |
| Dispatch | `create_momentum_balance_solver` `select case` | add `case ('SSA_FEM_PETSc')` |
| Config | `choice_stress_balance_approximation_config` in `model_configuration_type_and_namelist.f90` | extend the comment list of valid values; add new PETSc-FE knobs |
| Base class | `atype_momentum_balance_solver` (via `atype_momentum_balance_solver_data` / `atype_model`) | extend **`atype_momentum_balance_solver` directly** (not `_SSADIVA`, since we do not reuse the CSR stiffness assembly) |
| Result hand-off | `set_velocities_to_solver_results` writes `vel%u_3D_b`, strain rates, etc. | same, filled from the projected nodal solution |

### Files

```
src/UFEMISM/ice_dynamics/momentum_balance/SSA_FEM_PETSc/
  momentum_balance_solver_SSA_FEM_PETSc.f90     ! DONE (Phase 0); grows through the phases
  SSA_FEM_PETSc_fields.f90                       ! (Phase 2) aux-field DM: build + fill + DMSetAuxiliaryVec
src/UPSY/basic/petsc/petsc_fe.f90              ! (later) shared bind(C) bindings missing from PETSc's Fortran module
```

The `bind(C)` interfaces and pointwise functions currently live inside
`momentum_balance_solver_SSA_FEM_PETSc.f90`; promote them to `petsc_fe.f90` when a
second caller appears. CMake picks up new `.f90` files under `src/` automatically
(`GLOB_RECURSE`), no `CMakeLists.txt` edit needed.

## 3. Discretisation (recap of the design decision)

- **Velocity**: continuous Lagrange **P1, 2 components**, on mesh vertices.
  `PetscFECreateLagrange(comm, dim=2, Nc=2, isSimplex=PETSC_TRUE, k=1, qorder=-1, fe, ierr)`.
- **Auxiliary (data) fields**, P1 on vertices, on a cloned DM, read-only in the
  residual: ice thickness `H`, surface elevation `s` (for the driving stress
  `-rho g H grad s`), vertically averaged flow factor `A`, basal friction
  coefficient `beta` (already includes the sub-grid grounded fraction), and any
  masks needed for the friction regularisation.
- **Weak form** (single vector field, index 0):
  - `f1[i][j]` = membrane stress tensor `2 eta H (2 eps_dot + tr(eps_dot) I)_{ij}`,
    with `eta = 1/2 A^{-1/n} eps_eff^{(1-n)/n}` computed pointwise from `u_x`.
  - `f0[i]` = basal drag `beta(|u|) u_i` + driving stress `rho g H (grad s)_i`
    (driving stress taken from the aux surface-elevation gradient `a_x`).
  - **Ice-front / calving-front** back-pressure `1/2 rho g H^2 (...) n_i`: boundary
    residual via `PetscDSSetBdResidual` + `DMAddBoundary(DM_BC_NATURAL, ...)`.
- **Regularisation**: keep the existing `Glens_flow_law_epsilon_sq_0` on the
  effective strain rate; regularise `|u|` in the sliding law so Newton stays
  differentiable near stagnation.

## 4. Phased steps

### Phase 0 - Scaffolding (compiles, selectable, does nothing) - **done**

- `type_momentum_balance_solver_SSA_FEM_PETSc`
  (`src/UFEMISM/ice_dynamics/momentum_balance/SSA_FEM_PETSc/momentum_balance_solver_SSA_FEM_PETSc.f90`)
  extends `atype_momentum_balance_solver`; all deferred procedures implemented.
- `run` is a no-op (does **not** crash - the repo owner asked for a runnable
  placeholder); `set_velocities_to_solver_results` writes zeros to the b-grid
  velocities and vertex strain rates.
- Dispatch `case ('SSA_FEM_PETSc')` and the config valid-values comment added.
- Verified: `dev/changed` build is clean and the integrated test runs the
  `SSA_FEM_PETSc` config end to end (zero velocity, no crash).

### Phase 1 - DMPlex + PetscFE field + SNES skeleton, linear constant-coefficient SSA - **done**

Stood up the full DMPlex -> PetscFE -> PetscDS -> SNES pipeline and solved a
**constant-coefficient linear** SSA, verified against a closed-form result.
All of it lives in `momentum_balance_solver_SSA_FEM_PETSc.f90`:

- `build_petsc_objects` / `destroy_petsc_objects` (called from
  `initialise`/`deallocate`, and both from `remap`).
- `bind(C)` interfaces for `PetscDSSetResidual`, `PetscDSSetJacobian`,
  `SNESSetJacobian`, `SNESGetConvergedReason`, `DMPlexSetSNESLocalFEM`.
- pointwise `SSA_FEM_PETSc_f0/_f1/_g0/_g3` (constant `N`, `beta`, `tau_d` via
  `PetscDSSetConstants`).
- `copy_PETSc_solution_to_mesh_vertices_vec2`: `DMGlobalToLocal` + local
  `PetscSection` read + `MPI_Alltoallv` on `upsy_vertex_id` /
  `mesh%V_owning_process`, then `map_a_b_2D` to the b-grid.
- KSP `preonly` + PC `lu` (MUMPS on >1 rank); `run` checks the result against
  `u = -tau_dx/beta`, `v = -tau_dy/beta` and `warning`s on mismatch.

Result on `integrated_test_SSA_notime_MISMIP_mod_full` (5511-vertex mesh, 2 ranks):
`SNES its = 1`, `max|u - u_exact| = 7e-16`. The `overlap = 0` open check is
resolved for this case (parallel result is exact); revisit once coefficients vary
in space.

Original target for reference:

**Weak form** (PETSc convention `residual = integral( f0.phi + f1:grad(phi) ) = 0`),
single 2-component vector field, constant `N = eta*H`, constant `beta`, constant
driving stress `(tau_dx, tau_dy)`:

- `f1[u,x] = 2N(2 du/dx + dv/dy)`, `f1[u,y] = N(du/dy + dv/dx)`,
  `f1[v,x] = N(du/dy + dv/dx)`, `f1[v,y] = 2N(2 dv/dy + du/dx)`
- `f0[u] = beta*u + tau_dx`, `f0[v] = beta*v + tau_dy`
- analytic Jacobian: `g0 = beta*I` (2x2), `g3 = d f1 / d grad(u)` (8 non-zeros of 16).

**Well-posedness / BCs.** The `beta*I` (basal drag) term makes the operator
positive-definite with *natural* (do-nothing) boundary conditions, so Phase 1
imposes **no** essential BCs. Real BCs (prescribed velocity, ice front, periodic)
are Phase 5.

**Closed-form check.** With constant coefficients, constant forcing and natural
BCs the exact solution is the spatially uniform field
`u = -tau_dx/beta`, `v = -tau_dy/beta`. Phase 1 solves on the integrated-test
mesh and checks the returned nodal field is uniform and equals `-tau/beta` to
solver tolerance.

Steps:

1. `initialise` -> `build_petsc_objects`: `mesh_to_dmplex` (stored on `self%dm`);
   `PetscFECreateLagrange(comm, dim=2, Nc=2, simplex, k=1)`; `DMSetField`;
   `DMCreateDS`; `PetscDSSetConstants([N, beta, tau_dx, tau_dy])`;
   `PetscDSSetResidual`/`PetscDSSetJacobian` via `bind(C)` (pattern from
   `ct_PETSc_SNES_Poisson.f90`); `SNESCreate` -> `SNESSetDM` ->
   `DMPlexSetSNESLocalFEM` -> `DMCreateMatrix` -> `SNESSetJacobian(snes,J,J,NULL)`;
   KSP `preonly` + PC `lu` (auto-MUMPS on >1 rank, as in the Poisson test).
2. `run`: `DMCreateGlobalVector`, `VecSet(0)`, `SNESSolve`; check
   `SNESGetConvergedReason` (`bind(C)` wrapper); copy back.
3. Copy back: `DMGlobalToLocal` -> read the 2 dofs per vertex via the local
   `PetscSection` -> `MPI_Alltoallv` on `upsy_vertex_id` + `mesh%V_owning_process`
   into `u_vav_a`/`v_vav_a` (`vi1:vi2`) -> `map_a_b_2D` to `u_vav_b`/`v_vav_b`.
4. `set_velocities_to_solver_results`: `vel%u_3D_b(ti,:) = u_vav_b(ti)` etc.
   (strain rates still zeroed until Phase 3).
5. `deallocate` -> `destroy_petsc_objects` (`MatDestroy`/`SNESDestroy`/
   `PetscFEDestroy`/`DMDestroy`); `remap` -> destroy + rebuild for the new mesh.

**Open check (deferred within Phase 1):** FE assembly overlap. `mesh_to_dmplex`
distributes with `overlap = 0`. If the 2-rank result disagrees with the 1-rank
result, add an `overlap = 1` path to `mesh_to_dmplex` (new optional argument,
default unchanged).

### Phase 2 - Auxiliary fields (real spatially varying coefficients) - **done**

What was built (all in `momentum_balance_solver_SSA_FEM_PETSc.f90`):

1. **Aux DM.** `bind(C)` `DMClone` of the primary DM, one P1 `PetscFE` with
   `Nc = 4`, `DMCreateDS`, a persistent local `Vec` (`aux_vec`). `bind(C)`
   `DMSetAuxiliaryVec(dm, NULL, 0, 0, aux_vec)`.
2. **`calc_auxiliary_fields`** computes the per-vertex coefficients
   `a = [N, beta, tau_dx, tau_dy]` from the current velocity solution:
   `eta` from `calc_effective_viscosity_Glen_2D` (strain rates via `ddx_a_a_2D`
   on the nodal field, `A_vav` from `calc_ice_rheology_Glen` + `vertical_average`),
   `N = eta * max(0.1, H)`, driving stress `-rho g H grad(Hs)` via `ddx_a_a_2D`,
   `beta` from `calc_basal_friction_coefficient`, then scaled by the sub-grid
   grounded fraction: `beta *= geom%fraction_gr**C%subgrid_friction_exponent_on_B_grid`
   when `C%do_GL_subgrid_friction` (the a-grid analogue of what
   `calc_applied_basal_friction_coefficient` does on the b-grid). This is
   essential, not optional: without it, full basal friction under the floating
   shelf makes it a different physical problem - it was the main reason the
   Phase-2 solution first looked nothing like the finite-difference `SSA` result.
3. **`fill_PETSc_aux_from_mesh_vertices`** - inverse of the solution copy-back:
   each rank requests, for its local DMPlex vertices, the 4 coefficients from the
   UFEMISM process that owns that vertex (two `MPI_Alltoallv`), then
   `VecSetValues` at the local section offsets.
4. **Weak form** `f0/f1/g0/g3` now read `a[]` instead of `PetscDSSetConstants`.
   `f1` needs only `N` (undifferentiated - the FE weak form has `integral(N ... :
   grad phi)`, so no `grad N` term, unlike the finite-difference assembly).
5. **Outer loop.** A Picard viscosity iteration in `run` (freeze coefficients ->
   linear SNES solve -> relax -> velocity limit -> L2 stop), reusing
   `C%visc_it_nit`, `C%visc_it_relax`, `C%visc_it_norm_dUV_tol`, `C%vel_max`.

Result on `integrated_test_SSA_notime_MISMIP_mod_full` (2 ranks): runs clean, no
NaN, velocities grow from rest to a physically scaled, spatially structured field.

Vertexwise correlation of `uabs_vav` with the finite-difference `SSA` solver is
**0.95** (same spatial structure; FE mean speed ~316 m/yr vs FD ~551 m/yr).

**Convergence caveat.** With fixed relaxation the Picard loop still hits the
50-iteration cap (L2 ~4e-6, decreasing; velocities still ramping) - and **the
finite-difference `SSA` solver fails to converge on this same config too**, so it
is inherent to the setup/settings, not the FE discretisation. The remaining
FD-vs-FE gap is dominated by both solvers being under-converged with different
Picard damping and different strain-rate discretisations. The proper fix is
Phase 3: replace the hand-rolled Picard with a single SNES solve of the true
nonlinear residual (`eta = eta(grad u)` pointwise), optionally with the
adaptive-relaxation Picard from `momentum_balance_solver_SSA` as a fallback
(Phase 4).

### Phase 3 - Non-linear Newton residual + analytic Jacobian - **done**

(Phase 3 and the analytic-Jacobian half of Phase 4 were done together - once the
`d eta / d grad u` term was written for the residual it was trivially also the
Jacobian, so the FD-coloured-Jacobian intermediate step was skipped.)

- **`f1`** computes `eta` pointwise from `grad u`
  (`eps2 = ux^2 + vy^2 + ux*vy + (uy+vx)^2/4 + eps0`,
  `eta = 1/2 Abar^(-1/n) eps2^((1-n)/2n)`), then `f1[c,d] = 2 (eta*H) D[c,d]`
  with `D` the SSA strain tensor. `Abar`, `H` come from the aux field; `eps0`,
  `n` from `PetscDSSetConstants`.
- **`f0` = `beta*u - tau_d`** (basal drag minus the driving stress). The minus
  sign matters - an earlier `+` produced a velocity field anti-correlated with
  the FD SSA solver (right speed, flipped direction).
- **`g3`** (analytic) = `2 N dD/d(grad u)` (the old frozen term, `N` now
  pointwise) `+ coef * D[m] D[k]` with
  `coef = 2 H eta ((1-n)/2n) / eps2` (rank-1 shear-thinning term). `g0 = beta*I`.
- **`f0` = `beta(|u|)*u - tau_d`** with `beta` also evaluated pointwise
  (`SSA_FEM_PETSc_sliding_beta`, Zoet-Iverson: `beta = tauc_eff |u|^(1/p-1)
  (|u|+u_t)^(-1/p)`, `|u|` regularised by `delta_v`). `g0` = `beta*I` + rank-1
  `(dbeta/d|u| / |u|) u_c u_c'`.
- **No outer loop.** Both non-linearities are in the residual, so `run` does a
  single `SNESSolve` on a persistent, warm-started solution vector. The Picard
  viscosity/friction iteration is gone.
- **Aux field** is `[Abar, H, tauc_eff, tau_dx, tau_dy]` (`tauc_eff` = till yield
  stress * `fraction_gr**exp`); `fill_PETSc_aux_from_mesh_vertices` generalised to
  `n_aux_comp`. `eps0`, `n`, ZI `p`, `u_t`, `delta_v`, `beta_max` via
  `PetscDSSetConstants`.
- `momentum_balance_solver_SSA_FEM_PETSc_initialise` guards
  `C%choice_sliding_law` (only `Zoet-Iverson` / `no_sliding` supported so far;
  TODO in `SSA_FEM_PETSc_sliding_beta` to move the pointwise kernel into
  `sliding_laws` and cover the other laws).

Result (integrated test, 2 ranks, **cold start**): the single Newton solve
converges in **10 iterations**. Signed `u_vav`/`v_vav` correlate **+0.98** with
the FD `SSA` solver (`uabs_vav` +0.96). FE speeds are ~1.4x the FD solver's,
because the FD solver's own Picard does **not** converge on this near-plastic-
friction config (mean speed climbs 550 -> 617 m/yr and correlation with the FE
result rises 0.983 -> 0.987 as its cap is raised 50 -> 500), while the single
Newton solve is properly converged. The rest of the gap is the different strain-
rate discretisation and the grounding-line `beta` representation (nodal `tauc`
vs FD's triangle-averaged `beta`), which Phase 5 (BCs / GL) should narrow.

Verification still worth doing: a manufactured / Schoof case where both solvers
converge, to confirm the FE solution is the correct one (Phase 7).

**Time-evolving runs needed a much smaller timestep than with the FD `SSA`
solver - root-caused and fixed.** Diagnosis: `calc_critical_timestep_adv`
(the hard advective CFL) takes a global minimum over mesh edges of
`dist/(|u_c|+|v_c|)`, with zero smoothing, so a single anomalous edge speed sets
the timestep for the whole domain. The edge (c-grid) velocity `u_c` was produced
by remapping the solver's b-grid output (itself already a vertex -> triangle
remap of the FE solution) triangle -> edge - i.e. **two consecutive averaging
remaps**, vertex -> triangle -> edge, each edge value an average of 4 vertex
values - which added real numerical diffusion on top of whatever the momentum
solve itself produced. Fixed with a new `map_velocities_from_a_to_c_2D`
(`map_velocities_to_c_grid.f90`) that maps vertex -> edge **directly**, and
critically, picks the single **upwind** vertex value per edge instead of
averaging the two endpoints - consistent with the general finding (see
`SSA_FEM_PETSc_weak_form_derivation.md`'s companion discussion earlier in this
project) that mass-continuity stability comes from upwinding the advective flux,
not from where the velocity unknown formally lives. This resolves the concern
raised when Phase 0 started, that moving the velocity unknown off the staggered
b-grid might by itself destabilise mass continuity - it does not, provided the
consumer-facing remap is direct and upwinded. No DMPlex/mesh-topology change
(e.g. building the DMPlex from the dual/Voronoi mesh) was needed or is planned;
`PetscFECreateLagrange`'s continuous Lagrange elements only support simplices/
tensor cells, not the general polygons of a Voronoi dual, so that path would have
required abandoning the current PetscFE approach rather than adapting it.

### Phase 4 - Analytic pointwise Jacobian + Picard option

1. Implement `g0` (d f0 / d u: basal-drag linearisation), `g3` (d f1 / d u_x:
   the frozen-viscosity term **plus** the `d eta / d eps_dot` shear-thinning
   terms), `g1`/`g2` if any cross terms remain.
2. Add a config switch `SSA_PETSc_nonlinear_solver = 'Newton' | 'Picard'`:
   - `Newton`: analytic `g*` as above.
   - `Picard`: assemble only the frozen-viscosity part as the preconditioning
     matrix while the residual stays fully nonlinear (`SNESSetPicard`, or a
     Newton-LS with the Picard operator as PC). Expected to be the robust
     default far from the solution; Newton for polish.
3. Cross-check the analytic Jacobian against `-snes_test_jacobian`.

### Phase 5 - Boundary conditions from UFEMISM

Starting with the ice-front ocean back-pressure (the repo owner's priority - the
finite-difference solver doesn't have this at all); `zero`/`periodic`/etc.
domain-edge BCs come later.

#### Ice-front ocean back-pressure - design (in progress)

**The problem.** UFEMISM's mesh spans the *entire* fixed config rectangle
(confirmed: `xmin/xmax/ymin/ymax`, open ocean explicitly included, mesh extent
asserted unchanged across remeshing) - the calving front is an *interior*
curve of that mesh, not its outer edge. `DMAddBoundary`/`PetscDSSetBdResidual`
are built around the DMPlex's *topological* exterior boundary (faces with
support size 1); checked directly against installed-PETSc source references
that the boundary-residual assembly reads a face's `support[0]` only, so
feeding it a label of genuinely interior faces (two neighbouring cells, ice and
open ocean) would silently use whichever cell happens to be `support[0]` -
arbitrary, unsafe, not something to build physics on. This is also why the
existing FD solver's `choice_BC_u/v_*` only ever fire on `mesh%TriBI/VBI`
(border-of-the-config-rectangle) triangles/vertices, never at the margin, and
why there is currently no calving-front treatment anywhere in the momentum
balance (confirmed by reading `solve_linearised_SSA_DIVA_infinite_slab.f90`
and `momentum_balance_solver_SSADIVA.f90`).

**Rejected alternatives:**
- Feeding `DMAddBoundary` an interior-face label directly - unsafe (`support[0]`
  ambiguity above).
- A hand-rolled nodal force added post-hoc into the assembled residual - works
  in principle (the term has zero Jacobian, being a function of `H`/`Ho` only)
  but sidesteps PETSc's normal machinery entirely and needs its own geometry
  (margin edges, ad hoc normal/length weighting) invented from scratch.
- Reusing UFEMISM's existing `graph`/`is_border`/`border_nhat` abstraction
  (`src/UPSY/mesh/graph/`) - **rejected by the repo owner**: that graph *is*
  UFEMISM's own in-house implementation of "a mesh built from only the
  ice-covered vertices", i.e. it duplicates exactly the PETSc-native mechanism
  below, in a non-DMPlex data structure we'd rather not depend on.

**Chosen approach: a genuine ice-covered sub-DMPlex, built directly (not via
`DMPlexFilter`).** Restrict the solver's DMPlex to ice-covered cells so the
calving front becomes the sub-mesh's *real* topological exterior boundary,
where the standard `DMPlexMarkBoundaryFaces` + `DMAddBoundary(DM_BC_NATURAL)` +
`PetscDSSetBdResidual` pipeline applies exactly as designed (outward normals
supplied automatically by PETSc's boundary pointwise-function arguments - no
manual normal/length geometry needed at all). `DMPlexFilter` (confirmed
available with a Fortran binding) was the first idea, but the repo owner
redirected to something cleaner: **`mesh_to_dmplex_masked`, a copy of
`mesh_to_dmplex` that takes a triangle mask and builds the restricted topology
directly**, rather than building the full DMPlex first and filtering it down.
This reuses the exact same, already-proven cone/chart construction as
`mesh_to_dmplex` (just skipping non-masked triangles/their unused
vertices/edges) instead of depending on `DMPlexFilter`'s less-well-documented
behaviour (halo/ownership-transfer handling, coordinate transfer).

Steps:
1. **`mesh_to_dmplex_masked( mesh, mask_tri, dm)` - done**, in the shared
   `src/UPSY/basic/petsc/petsc_dmplex.f90` (alongside `mesh_to_dmplex`, publicly
   exported). Given a triangle mask, it: marks which vertices/edges are touched
   by at least one masked triangle; builds point-translation tables and the
   chart/cone topology using *only* masked triangles and the edges/vertices
   they touch (so an edge bordering exactly one masked triangle - its other
   neighbour excluded, or on the parent mesh's own outer border - becomes a
   genuine exterior face of the sub-mesh); sets the same `upsy_vertex_id` label
   (restricted to included vertices) and coordinates (restricted likewise); then
   distributes exactly as before. `mask_tri` must be identical on every
   process, since mesh connectivity (`mesh%V`, `mesh%Tri`, `mesh%TriE`,
   `mesh%EV`) is itself fully replicated on every rank.
2. **Verified end-to-end on the integrated test**, first with a temporary
   diagnostic helper (`verify_ice_covered_dmplex`, since superseded - built the
   ice mask, called `mesh_to_dmplex_masked`, marked boundary faces, reported
   counts, then discarded the sub-mesh without using it in the actual solve).
   The ice mask: a triangle counts as ice-covered if all 3 vertices have
   `mask_grounded_ice .or. mask_floating_ice` (matches the existing graph
   abstraction's `Hi > 0` rule); gathered from the "dist-shared" mask fields via
   `gather_dist_shared_to_all( mesh%pai_V, ...)` - **not** `gather_to_all`,
   which is for plainly-distributed (non-shared-memory) fields and errors on
   these (`combined sizes of d_partial dont match size of d_tot`).
   Result on the MISMIP_mod test, 2 ranks: **8581 / 10805 triangles
   ice-covered; 179 margin (exterior) faces per rank** - a plausible calving-
   front perimeter for this ice sheet, and no crash. That mask logic is now the
   permanent `calc_ice_covered_triangle_mask` helper, and the sub-mesh it
   builds is what the actual solve runs on (item 4 below) - it's no longer
   just a diagnostic.
3. **Done**: `self%dm` (and everything downstream of it - `self%fe`, the
   primary `PetscDS`, the residual/Jacobian callbacks, `self%dm_aux`/
   `self%fe_aux`/`self%aux_vec`, `self%jac`, the near-null-space, `self%snes`)
   is now built from `mesh_to_dmplex_masked` (the ice-covered sub-mesh) instead
   of the whole-mesh `mesh_to_dmplex`. **Still to do**: the back-pressure BC
   term itself - `DMPlexMarkBoundaryFaces` + `DMAddBoundary(DM_BC_NATURAL,
   ...)` + `PetscDSSetBdResidual` on the sub-mesh's new exterior (the ice
   margin), with `f0_bd,i = -tau_o,i`,
   `tau_o,i = (1/2 rho_i g H^2 - 1/2 rho_sw g Ho^2) n_i` (`n` supplied by
   PETSc; no boundary Jacobian needed, since `tau_o` depends only on the aux
   field, not on `u`). `Ho` matches `geom%Ho` exactly (`height_of_water_column`,
   already computed live every step - reuse it, don't recompute):
   `Ho = min(max(SL - Hb, 0), (rho_i/rho_sw) H)`. Until this lands, the
   sub-mesh's margin is a plain natural (zero-traction) boundary - i.e. the
   weak form's implicit "no boundary term" default - not yet the physical
   back-pressure condition.
4. **Done**: the full rebuild moves from "once at `initialise`, again on
   `remap`" to every `run()` call, since the grounding line / calving front can
   migrate every timestep even without a full remesh - and it really is the
   *full* rebuild, not just `self%dm`. Implemented exactly as the checklist
   below describes; `initialise` no longer builds anything (it doesn't have
   `geom`), `remap` only destroys + reallocates the plain arrays, and `run`
   unconditionally does `destroy_petsc_objects` (if built) + `build_petsc_objects(
   geom)` right after the "no grounded ice" early return (which itself now
   guards its `VecSet(self%sol, ...)` behind `self%petsc_is_built`, since
   `self%sol` may not exist yet on a first call with no ice anywhere). Checklist,
   all now implemented inside `build_petsc_objects( self, geom)`:
     1. ice mask -> `self%dm` (`mesh_to_dmplex_masked`, via the extracted
        `calc_ice_covered_triangle_mask` helper);
     2. `self%fe` + `DMSetField` + `DMCreateDS` -> a brand new `PetscDS`;
     3. `PetscDSSetConstants` re-set on the new DS;
     4. residual/Jacobian callbacks (`petsc_ds_set_residual`/`_jacobian`)
        re-registered on the new DS;
     5. the ice-margin boundary face label + `DMAddBoundary`/
        `PetscDSSetBdResidual` - deferred along with the BC term itself (item 3
        above);
     6. `self%dm_aux` (clone of the *new* `self%dm`) + `self%fe_aux` +
        `DMSetField` + `DMCreateDS` + `self%aux_vec`;
     7. `self%jac` - new sparsity pattern (`DMCreateMatrix` on the new `dm`);
     8. the rigid-body near-null-space - tied to the new `jac` and the new
        `dm`'s coordinates;
     9. `self%snes` itself - `SNESSetDM`, `DMPlexSetSNESLocalFEM`,
        `SNESSetJacobian`, KSP/PC type, KSP/SNES tolerances all reapplied to
        fresh objects (the existing `destroy_petsc_objects` ->
        `build_petsc_objects` cycle, just triggered every `run()` instead of
        only at `initialise`/`remap`);
     10. `self%sol` - a new `DMCreateGlobalVector` on the new `dm`, **re-seeded,
         not just zeroed**: DOF numbering from `mesh_to_dmplex_masked`/
         `DMPlexDistribute` is not guaranteed stable across independent
         rebuilds, so the old `self%sol` would be the wrong layout for the new
         DM. `run()` scatters the previous *physical* solution
         (`self%u_vav_a`/`v_vav_a`, converted back to dimensionless units) into
         the new `self%sol` right after the rebuild, reusing
         `fill_PETSc_aux_from_mesh_vertices` (now generalised to an arbitrary
         component count via `size(coeffs, 2)`, rather than the hard-coded
         `n_aux_comp`) with `dm_topo = dm_aux = self%dm` and a temporary local
         vector scattered in via `DMLocalToGlobalBegin`/`End`.

   Also updated as part of this: `copy_PETSc_solution_to_mesh_vertices_vec2` no
   longer crashes when `ncopies == 0` (a vertex outside the ice-covered
   sub-mesh) - such vertices are simply left at their `u = v = 0` default
   instead, matching the "no grounded ice" convention already used elsewhere.

   Verified on the integrated test, 2 ranks: builds and runs cleanly, Newton
   converges in the same **10 iterations** as before (10 cumulative KSP its,
   `lu`); max speed is now **1.193e3 m/yr** (was 2.120e3 on the whole-mesh
   solve) - expected to shift, since the domain, and therefore the discrete
   problem being solved, has genuinely changed (ice-free area, and its
   zero-thickness/zero-stress contribution, is now excluded rather than padded
   in with the `max(0.1, Hi)` floor, and the ice margin is a real, if still
   physically-incomplete, boundary rather than an artefact of the whole-mesh
   discretisation) - not a regression, but the expected result of switching
   discretisations ahead of adding the actual back-pressure term (item 3).

   This is a real performance cost (rebuilding DMPlex+FE+DS+SNES every
   timestep) to revisit in Phase 6 tuning once correctness is established -
   e.g. detecting an unchanged ice mask and skipping the rebuild, or keeping
   `self%fe`/`self%fe_aux` across rebuilds (they are reference-element
   objects, not mesh-sized, so in principle they don't need recreating, only
   re-attaching via `DMSetField` to each new `dm`/`dm_aux`) - not attempted
   in the first, correctness-first pass.

   Not yet done: the `max(0.1, Hi)` thickness floor is still applied inside
   `calc_auxiliary_fields` even though every vertex in the sub-mesh now has
   real ice - harmless (the floor never binds there) but could be dropped as a
   cleanup once the BC work (item 3) lands. `BC_prescr_mask_b`/
   `BC_prescr_u_b`/`BC_prescr_v_b` (prescribed velocity on triangles) and the
   domain-edge `choice_BC_u/v_*` options are unaffected by this change and
   remain future work (item 6 below).

Formula, config option and reference already exist in the repo, unused -
match them exactly rather than reinventing:
- `C%BC_ice_front` (`'infinite_slab'` | `'ocean_pressure'`,
  `model_configuration_type_and_namelist.f90:326`) is declared and copied to
  `C%BC_ice_front` but **never read anywhere** - wire into this option, don't
  add a new one.
- The exact formula, sign convention and citation (**Robinson et al., 2020, Eq.
  19**) are fully written out in the dead code
  `solve_linearised_SSA_DIVA_ocean_pressure.f90:443-466` and
  `DIVA_solver_ocean_pressure.f90:616-675` (both commented out, for the
  graph-based DIVA path) - use these as the reference derivation, and note in
  code comments that this FE implementation is the first *live* one.

6. Domain-edge choices (`choice_BC_u/v`: `zero`, `infinite slab`,
   `periodic_ISMIP-HOM`) and prescribed-velocity masks - unchanged from the
   original plan, still future work:
   - zero / prescribed: `DM_BC_ESSENTIAL`.
   - periodic: build the DMPlex with periodicity, or add the periodic face pairs
     as a constraint (`find_ti_copy_*` in `mesh_utilities` gives the partner);
     simplest first target is the non-periodic benchmarks.

### Phase 6 - Scaling, solver, preconditioner

1. **Nondimensionalisation** - **done** (pulled forward ahead of Phases 4/5, see
   the Progress table). Implemented as a pure change of variables around the
   existing pointwise functions, no re-derivation needed:
   - The PetscFE field (and `self%sol`) holds `u_hat = u / velocity_scale`
     (`velocity_scale = 1e3 m/yr`), not physical velocity.
   - `f0`/`f1` convert `u_hat`/`grad(u_hat)` to physical units, run the physics
     unchanged, then divide by `stress_scale` (`= 1e5 Pa`).
   - `g0`/`g3` get the same physical-unit treatment, then an extra factor
     `velocity_scale/stress_scale` (from the chain rule through
     `u = velocity_scale * u_hat`).
   - Only the copy-back to `u_vav_a`/`v_vav_a` (in `solve_SSA_Newton`) multiplies
     back by `velocity_scale`; the SNES tolerances (item 4 below) now apply to
     this dimensionless residual, which is more sensible than raw SI units.
   - **Result on the integrated test, 2 ranks**: bit-for-bit reproduces the
     pre-nondim answer with direct LU (`Newton its = 10`, `max speed =
     2.120e3 m/yr`) - confirms it's a transparent change of variables, not a
     physics change.
   - **A/B experiment**: swapped `KSPPREONLY`+`PCLU` for `KSPGMRES`+`PCBJACOBI`
     (the finite-difference `SSA` solver's own default, see
     `solve_matrix_equation_PETSc`). Converges to the identical answer in the
     same number of Newton iterations either way; non-dimensionalising reduces
     the cumulative KSP iteration count from **625 to 583** (~7%) on this test.
     A real but modest effect here - as expected, since a single global
     constant rescaling cannot equalise the intrinsic ~1e8 ratio between the
     membrane-stress coefficient `N` and the basal-friction coefficient `beta`
     (that ratio is physical - shelf/fast-stream vs. slow interior - not a
     units artefact). LU stays the default in the committed code; switching it
     is left to item 2 below, a separate decision.
2. **Solver config** - **done**. `SNESNEWTONLS` (unchanged) + a config-selected
   KSP/PC (`build_petsc_objects`, `select case (C%SSA_FEM_PETSc_pc_type)`):
   `'lu'` (`KSPPREONLY`+`PCLU`, the validated default), `'gamg'`
   (`KSPGMRES`+`PCGAMG`), `'bjacobi'` (`KSPGMRES`+`PCBJACOBI`, matching the
   finite-difference solver's own default). Explicit `KSPSetTolerances`
   (`rtol=1e-8`, `abstol=1e-12`, `maxits=10000`) for the two iterative cases.
3. **Near-null-space** - **done**. Rigid-body modes (2 translations + 1
   rotation in 2D) built from the DM's own coordinates and attached
   unconditionally in `build_petsc_objects`:
   `MatSetBlockSize(jac,2)` -> `DMGetCoordinates` -> `MatNullSpaceCreateRigidBody`
   -> `MatSetNearNullSpace` -> `MatNullSpaceDestroy`. Harmless for `'lu'`/`'bjacobi'`
   (near-null-space is only consulted by AMG); required for `'gamg'` to coarsen
   this vector-valued (elasticity-like) operator well.

   **Result on the integrated test, 2 ranks, all three `pc_type` choices** (with
   the explicit KSP tolerances above): all converge to the *identical* answer
   (`Newton its = 10`, `max speed = 2.120e3 m/yr`) - only the cumulative KSP
   iteration count differs:

   | `pc_type` | cumulative KSP its |
   | --- | --- |
   | `lu` | 10 (one direct solve per Newton step) |
   | `gamg` | 362 |
   | `bjacobi` | 857 |

   `gamg` needs less than half `bjacobi`'s iterations, confirming the
   near-null-space is doing its job. `lu` stays the committed default (small
   test meshes, zero solver-tuning risk); `gamg` is the path to scaling up
   later, per the very first architecture discussion in this project.
4. **Config knobs** - **done** (added to `model_configuration_type_and_namelist.f90`
   and `config_SSA_FEM_PETSc.cfg`): `SSA_FEM_PETSc_pc_type_config` (`'lu'` |
   `'gamg'` | `'bjacobi'`, default `'lu'`), `SSA_FEM_PETSc_snes_rtol_config`
   (default `1e-8`), `SSA_FEM_PETSc_snes_abstol_config` (default `1e-10`),
   `SSA_FEM_PETSc_snes_maxits_config` (default `50`). These are solver-specific,
   distinct from the generic `stress_balance_PETSc_rtol/abstol` the
   finite-difference solvers use, both because ours is a non-linear (SNES, not
   plain KSP) solve and because non-dimensionalisation changed what these
   tolerances mean for us.
   `SSA_PETSc_nonlinear_solver` (Newton vs Picard) was already moot (Phase 4:
   no outer loop). `SSA_PETSc_velocity_element_order` (P1 vs P2) is deferred -
   no current need, and a bigger change (a second `PetscFECreateLagrange`
   degree, revisit only if accuracy demands it in Phase 7).

   Also added a permanent diagnostic: `run` now prints `Newton its`,
   cumulative `KSP its` (`self%n_Axb_its`, via `SNESGetLinearSolveIterations`)
   and `max speed` every solve, and `self%n_visc_its`/`n_Axb_its` now hold their
   intended meanings (nonlinear vs. linear iteration counts) instead of both
   being set from the same number.

### Phase 7 - Verification

1. **Analytic**: Schoof SSA ice stream — `src/UPSY/basic/analytical_solutions/Schoof_SSA_solution.f90`;
   compare L2 error and its order under uniform refinement (expect ~2 for P1).
2. **Cross-solver**: run the existing SSA integrated/component tests with
   `SSA_PETSc` and diff velocity fields against `SSA` (Halfar/dome, ISMIP-HOM,
   MISMIP+).
3. Add an L2-error helper via `bind(C)` to `DMComputeL2FieldDiff`
   (not in the 3.25.5 Fortran module - see the note in `ct_PETSc_SNES_Poisson.f90`).

### Phase 8 - Remap, restart, cleanup

1. `remap_momentum_balance_solver`: destroy DM/DS/FE/SNES, rebuild from
   `mesh_new`, re-init aux DM; remap `u_vav_b`/`v_vav_b` as the SSA solver does
   (via the a-grid).
2. Restart: reuse the SSA restart file layout (`u_vav_b`, `v_vav_b`).
3. `deallocate`: `SNESDestroy`, `PetscFEDestroy`, `DMDestroy` (main + aux),
   `MatDestroy`.
4. Remove the Phase 0 `crash`; document the solver in the config-file docs and
   the wiki page for stress-balance approximations.

## 5. New PETSc Fortran bindings likely required

PETSc 3.25.5's Fortran module is missing several symbols; add `bind(C)`
interfaces (pattern: `ct_PETSc_SNES_Poisson.f90` lines 49-116) in a shared
`src/UPSY/basic/petsc/petsc_fe.f90`:

- `PetscDSSetResidual`, `PetscDSSetJacobian` (already prototyped in the test - promote to the shared module)
- `PetscDSSetBdResidual`, `PetscDSSetBdJacobian`
- `DMAddBoundary`, `SNESGetConvergedReason`, `SNESSetJacobian` (promote from the test)
- `DMSetAuxiliaryVec`
- `DMComputeL2FieldDiff` (verification)
- `MatNullSpaceCreateRigidBody`, `MatSetNearNullSpace`
- `DMClone` if not exposed
- `DMProjectFunctionLocal` / `DMProjectFieldLocal` if used to fill aux fields

## 6. Edits to existing files (the complete list)

- `momentum_balance_solver_main.f90` - one `use`, one `case`.
- `model_configuration_type_and_namelist.f90` - valid-values comment on
  `choice_stress_balance_approximation_config`; new `_config` fields + their
  second declaration + namelist block + assignment (4 spots each, per existing
  convention).
- Source list / `CMakeLists.txt` - register new modules.
- `automated_testing/` - add an integrated/component test entry for the
  `SSA_PETSc` choice (new config, reference data); optionally extend
  `UPSY_component_test_program_PETSc_DMPLEX` with a standalone SSA-FE check.

## 7. Risks / open questions

- **Parallel FE assembly overlap** - confirm `overlap = 0` from
  `mesh_to_dmplex` is sufficient for `DMPlexSetSNESLocalFEM`; add an
  `overlap = 1` path if not.
- **b-grid <-> vertex BC translation** - prescribed-velocity masks live on
  triangles; the vertex Dirichlet mapping is approximate near the mask edge.
- **Periodic benchmarks** (ISMIP-HOM) need periodic DMPlex or explicit
  constraints; defer past first validation.
- **Grounding line** - non-smooth `beta` and grounded fraction hurt Newton;
  rely on Picard there, and on the existing sub-grid `fraction_gr_b` already
  folded into the `beta` aux field.
- **Fortran <-> C pointwise callbacks** - long signatures that must match PETSc
  exactly; keep them all in `SSA_PETSc_weak_form.f90` and unit-test the
  viscosity/Jacobian algebra against `constitutive_equation` in isolation.

## 8. Milestone checklist

- [x] Phase 0: `SSA_FEM_PETSc` selectable, runs as a no-op placeholder.
- [x] Phase 1: linear constant-coefficient SSA solves via SNES, uniform-field
      closed-form check passes (`~1e-16`, 2 ranks), result on the b-grid.
- [x] Phase 2: 4-component aux field `[N, beta, tau_dx, tau_dy]` (beta incl.
      sub-grid grounded-fraction scaling) drives the weak form; per-vertex
      coefficients + Picard loop; runs on 2 ranks; `uabs_vav` correlates 0.95
      with the FD SSA solver (tight Picard convergence deferred to Phase 3, as
      the FD SSA also stalls on this config).
- [x] Phase 3: fully non-linear residual - pointwise Glen `eta(grad u)` AND
      Zoet-Iverson `beta(|u|)` - with analytic Jacobians; **one Newton solve, no
      outer loop**; converges in 10 its (cold start); signed velocity correlates
      +0.98 with FD `SSA`. (`f0` sign fixed.)
- [x] Phase 4: analytic Jacobian done; the Picard-option item is moot for the
      pointwise-friction path (no outer loop). Adaptive relaxation only matters if
      a non-pointwise sliding law is added later.
- [ ] Phase 5: UFEMISM BCs (prescribed, ice front; periodic later). Ice-covered
      sub-mesh (`mesh_to_dmplex_masked`) now underlies the entire solve, rebuilt
      every `run()` call (incl. re-seeded `self%sol` warm start); the actual
      back-pressure BC term on its exterior (ice margin) is still to do.
- [x] Phase 6: nondimensionalisation (transparent, verified), config-selectable
      `pc_type` (`lu`/`gamg`/`bjacobi`) + rigid-body near-null-space (`gamg`
      under half `bjacobi`'s KSP its), new SNES config knobs. `lu` stays the
      default.
- [ ] Phase 7: Schoof convergence order + benchmark cross-checks.
- [ ] Phase 8: remap, restart, cleanup, docs.

## 9. Deferred follow-up work (separate PRs)

- **Nondimensionalise the finite-difference `SSA`/`DIVA` solver's linear solve.**
  `solve_matrix_equation_PETSc` (called from `solve_SSA_DIVA_linearised` via
  `solve_matrix_equation_CSR_PETSc`) already defaults to an **iterative** KSP
  (`gmres` + `bjacobi`), unlike `SSA_FEM_PETSc`'s direct LU - so it has been
  exposed to the same raw-SI-unit conditioning problem (coefficients spanning
  `~1e-8` to `~1e13`) the whole time, with no direct-solve fallback cushioning
  it. Plausibly a bigger win there than the modest 625->583 KSP-iteration
  reduction (~7%) measured on `SSA_FEM_PETSc` (see Phase 6), precisely because
  it has no such cushion.
  - **Proposed retrofit** (small, contained, no change to the row-by-row
    assembly in `calc_SSA_DIVA_stiffness_matrix_row_free`): wrap the already-
    assembled system in `solve_SSA_DIVA_linearised`, right around the existing
    `solve_matrix_equation_CSR_PETSc` call - scale the RHS by `1/stress_scale`
    before the solve, scale the returned solution by `velocity_scale` after it.
    Same `u_hat = u/velocity_scale`, `stress_scale` idea as `SSA_FEM_PETSc`,
    applied as a diagonal row/column rescaling around the existing black-box
    solve rather than inside the weak form.
  - **Verification plan**: A/B on the same integrated test as Phase 6's
    experiment - compare `n_Axb_its` (KSP iteration count, already reported by
    this solver) and the solution, before/after, to confirm it's transparent
    (identical answer) and measure the iteration-count effect.
  - **Why deferred**: touches a shared code path used by every existing
    production config and benchmark (SSA and DIVA both go through
    `solve_SSA_DIVA_linearised`), unlike everything else in this plan, which is
    purely additive. Repo owner wants this as its own, separate PR - explicitly
    requested 2025-09-08, not started.
