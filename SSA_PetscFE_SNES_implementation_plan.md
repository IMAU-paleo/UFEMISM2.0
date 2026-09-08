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
| 4, 5, 7, 8 | not started (Phase 4 analytic Jacobian folded into Phase 3; outer-loop / Picard-option is now moot for the pointwise-friction path) |

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

Map UFEMISM's BC concepts onto DMLabels + `DMAddBoundary`:

1. `BC_prescr_mask_b` / `BC_prescr_u_b` / `BC_prescr_v_b` (prescribed velocity on
   triangles): these arrive on the b-grid; convert to a vertex Dirichlet label
   (a vertex is constrained if all/most incident constrained triangles agree),
   or prescribe on the nearest vertices. Provide the values through the
   Dirichlet callback context.
2. Domain-edge choices already handled by
   `calc_SSA_DIVA_stiffness_matrix_row_BC` (`choice_BC_u/v`: zero, infinite
   slab, ISMIP-HOM periodic, ice-stream periodic). For the FE solver:
   - zero / prescribed: `DM_BC_ESSENTIAL`.
   - periodic: build the DMPlex with periodicity, or add the periodic face pairs
     as a constraint (`find_ti_copy_*` in `mesh_utilities` gives the partner);
     simplest first target is the non-periodic benchmarks.
3. Ice front: derive an ice-front face DMLabel from `geom` masks
   (`mask_cf` / floating vs open ocean) and attach the natural BC there.

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
- [ ] Phase 5: UFEMISM BCs (prescribed, ice front; periodic later).
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
