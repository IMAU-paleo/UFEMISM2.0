# From the SSA equations to the PETSc pointwise callbacks

This note derives, step by step, the functions `f0`, `f1`, `g0`, `g3` that the
`SSA_FEM_PETSc` solver hands to PETSc (in
`src/UFEMISM/ice_dynamics/momentum_balance/SSA_FEM_PETSc/momentum_balance_solver_SSA_FEM_PETSc.f90`).
It is written for glaciologists who know the SSA but not finite elements or
PETSc's weak-form API. Nothing here is PETSc-specific until Section 4.

Contents:

1. The SSA as two coupled PDEs (strong form)
2. Constitutive relations: effective viscosity and the sliding law
3. From the strong form to the weak (variational) form
4. PETSc's pointwise residual convention -> `f0`, `f1`
5. PETSc's pointwise Jacobian convention -> `g0`, `g3`
6. Deriving `g3` (the membrane-stress Jacobian)
7. Deriving `g0` (the basal-drag Jacobian)
8. Index bookkeeping: how the tensors are stored as flat arrays
9. What PETSc does with these functions (the assembly loop)
10. Sanity checks
11. Symbol - code - PETSc dictionary


## 1. The SSA as two coupled PDEs

Let $\mathbf{u}(x,y) = (u_1, u_2) = (u, v)$ be the vertically averaged horizontal
ice velocity, $H$ the ice thickness, $s$ the surface elevation. Write
$\partial_j \equiv \partial/\partial x_j$ and use the summation convention (a
repeated index is summed over $1,2$).

**Strain rate.** The horizontal strain-rate tensor is

$$
\dot\varepsilon_{ij} \;=\; \tfrac12\left(\partial_i u_j + \partial_j u_i\right),
\qquad
\dot{\boldsymbol\varepsilon}
= \begin{pmatrix}
\partial_1 u_1 & \tfrac12(\partial_2 u_1 + \partial_1 u_2)\\[2pt]
\tfrac12(\partial_2 u_1 + \partial_1 u_2) & \partial_2 u_2
\end{pmatrix}.
$$

**Membrane (resistive) stress.** Vertically integrating the deviatoric stress and
eliminating the vertical normal stress with incompressibility
($\dot\varepsilon_{zz} = -\dot\varepsilon_{kk}$) and the plane-stress assumption
gives the SSA membrane-stress tensor

$$
\boxed{\;\mathsf{M} \;=\; 2\,\eta\,H\,\bigl(\dot{\boldsymbol\varepsilon}
        + \operatorname{tr}(\dot{\boldsymbol\varepsilon})\,\mathsf{I}\bigr)
      \;\equiv\; 2\,\eta\,H\,\mathsf{D}\;}
$$

where $\eta$ is the effective viscosity (Section 2) and, written out,

$$
\mathsf{D} =
\begin{pmatrix}
2\,\partial_1 u_1 + \partial_2 u_2 & \tfrac12(\partial_2 u_1 + \partial_1 u_2)\\[4pt]
\tfrac12(\partial_2 u_1 + \partial_1 u_2) & 2\,\partial_2 u_2 + \partial_1 u_1
\end{pmatrix}
=
\begin{pmatrix} 2u_x + v_y & \tfrac12(u_y+v_x)\\[4pt] \tfrac12(u_y+v_x) & 2v_y + u_x \end{pmatrix}.
$$

$\mathsf{D}$ is the code's strain tensor `D` (a length-4 array, row-major:
`D(1)=D_xx`, `D(2)=D_xy`, `D(3)=D_yx`, `D(4)=D_yy`). We also write $N \equiv \eta H$
(the code's `N`), so $\mathsf{M} = 2N\mathsf{D}$.

**Driving stress.**

$$
\boldsymbol\tau_d \;=\; -\,\rho g H \,\nabla s .
$$

(In UFEMISM this is `tau_dx = -ice_density*grav*Hi*dHs_dx`, so the *components*
$\tau_{d,i}$ already carry their sign; on a surface sloping downhill in $+x$,
$\partial_x s < 0$ and $\tau_{d,x} > 0$.)

**Basal drag.**

$$
\boldsymbol\tau_b \;=\; \beta(|\mathbf{u}|)\,\mathbf{u},
$$

with a friction coefficient $\beta \ge 0$ that itself depends on the sliding
speed (Section 2). On floating ice $\beta = 0$.

**Force balance.** The SSA states that membrane-stress divergence, driving stress
and basal drag sum to zero:

$$
\boxed{\;
\partial_j \mathsf{M}_{ij}
\;+\; \tau_{d,i}
\;-\; \beta(|\mathbf{u}|)\,u_i
\;=\; 0
\qquad (i = 1,2).\;}
\tag{SSA}
$$

That is two scalar PDEs, one for each velocity component, coupled through
$\mathsf{M}$ (which mixes $u$ and $v$) and through $\eta$ and $\beta$.


## 2. Constitutive relations

**Effective viscosity (Glen's law).** With

$$
\dot\varepsilon_e^2
= \dot\varepsilon_{xx}^2 + \dot\varepsilon_{yy}^2
+ \dot\varepsilon_{xx}\dot\varepsilon_{yy} + \dot\varepsilon_{xy}^2
= u_x^2 + v_y^2 + u_x v_y + \tfrac14(u_y+v_x)^2 ,
$$

the vertically averaged effective viscosity is

$$
\eta \;=\; \tfrac12\, \bar A^{-1/n}\,
\bigl(\dot\varepsilon_e^2 + \dot\varepsilon_0^2\bigr)^{\frac{1-n}{2n}} ,
$$

with $\bar A$ the vertically averaged flow-rate factor, $n \approx 3$ the Glen
exponent, and $\dot\varepsilon_0^2$ a small regularisation so that $\eta$ stays
finite where the ice is not deforming. This is exactly
`calc_effective_viscosity_Glen_2D` and the helper `SSA_FEM_PETSc_strain`.
$\eta$ depends on $\mathbf{u}$ **only through $\nabla\mathbf{u}$** — this makes the
membrane term non-linear (shear thinning: $\partial\eta/\partial\dot\varepsilon_e^2 < 0$).

**Sliding law (Zoet & Iverson, 2020).** With a regularised speed
$|\mathbf{u}| = \sqrt{\delta_v^2 + u^2 + v^2}$,

$$
\beta(|\mathbf{u}|) \;=\; \tau_c\; |\mathbf{u}|^{\,1/p - 1}\;
\bigl(|\mathbf{u}| + u_t\bigr)^{-1/p} ,
$$

where $\tau_c$ is the till yield stress (here already multiplied by the sub-grid
grounded fraction $f_{gr}^{\,m}$, so friction vanishes under floating ice), $u_t$
the transition velocity and $p$ an exponent. This is `SSA_FEM_PETSc_sliding_beta`
and mirrors `calc_sliding_law_ZoetIverson`. $\beta$ depends on $\mathbf{u}$
**only through $|\mathbf{u}|$**, which makes the basal-drag term non-linear.


## 3. From the strong form to the weak form

Finite elements do not solve (SSA) directly. Instead we require that (SSA), when
multiplied by an arbitrary smooth **test function** $\boldsymbol\phi = (\phi_1,\phi_2)$
and integrated over the domain $\Omega$, gives zero:

$$
\int_\Omega \Bigl[\; \partial_j\mathsf{M}_{ij}
   \;+\;\tau_{d,i}\;-\;\beta u_i \;\Bigr]\,\phi_i \;\mathrm{d}\Omega \;=\; 0
\qquad\text{for every admissible } \boldsymbol\phi .
\tag{W0}
$$

The only troublesome term is $\partial_j\mathsf{M}_{ij}$: it needs a derivative of
the stress, hence a *second* derivative of the velocity, which a piecewise-linear
velocity field does not have. **Integration by parts** (the divergence theorem)
moves that derivative onto the test function:

$$
\int_\Omega \bigl(\partial_j\mathsf{M}_{ij}\bigr)\phi_i \,\mathrm{d}\Omega
=
\underbrace{\oint_{\partial\Omega} \mathsf{M}_{ij}\,n_j\,\phi_i \,\mathrm{d}\Gamma}_{\text{boundary term}}
\;-\;
\int_\Omega \mathsf{M}_{ij}\,\partial_j\phi_i \,\mathrm{d}\Omega .
$$

Two things happen here:

* We now only ever differentiate $\mathbf{u}$ **once** (inside $\mathsf{M}$) and
  $\boldsymbol\phi$ **once**. A continuous piecewise-linear (P1) field has exactly
  one (piecewise-constant) derivative, so it is an admissible solution *and* test
  function. This is the whole point of the weak form.

* The **boundary term** is where boundary conditions live:
  * On a piece of $\partial\Omega$ where the velocity is prescribed (an
    *essential* / Dirichlet condition), the test functions are taken to vanish
    there, $\phi_i = 0$, so the boundary term drops out.
  * On a piece where we do nothing, we are implicitly imposing
    $\mathsf{M}_{ij} n_j = 0$ (zero net membrane traction) - a *natural* boundary
    condition. The current solver does this **everywhere** (an ice front then
    behaves as stress-free; the missing $\tfrac12\rho g H^2$ ocean back-pressure
    is a later addition, and would be an extra boundary integral here).

Dropping the boundary term and substituting back into (W0):

$$
-\int_\Omega \mathsf{M}_{ij}\,\partial_j\phi_i \,\mathrm{d}\Omega
\;+\;\int_\Omega \bigl(\tau_{d,i} - \beta u_i\bigr)\phi_i \,\mathrm{d}\Omega \;=\; 0 .
$$

Multiply by $-1$ (this is a choice of sign; it does not change the solution, but
it makes $\mathsf{f1}$ come out equal to the *physical* stress $\mathsf{M}$):

$$
\boxed{\;
\int_\Omega \mathsf{M}_{ij}\,\partial_j\phi_i \,\mathrm{d}\Omega
\;+\;\int_\Omega \bigl(\beta u_i - \tau_{d,i}\bigr)\phi_i \,\mathrm{d}\Omega
\;=\; 0
\qquad\text{for all } \boldsymbol\phi .\;}
\tag{W}
$$

Equation (W) is the *weak form of the SSA*. Everything below is bookkeeping.


## 4. PETSc's pointwise residual convention -> `f0`, `f1`

PETSc's finite-element assembly (`DMPlexSetSNESLocalFEM`) expects the residual of a
single vector field written in the canonical form

$$
F[\boldsymbol\phi] \;=\;
\int_\Omega \Bigl[\;
   \phi_i\, f0_i(\mathbf{u}, \nabla\mathbf{u}, \mathbf{a}, \dots)
   \;+\;
   \partial_j\phi_i\, f1_{ij}(\mathbf{u}, \nabla\mathbf{u}, \mathbf{a}, \dots)
\;\Bigr]\,\mathrm{d}\Omega
\;=\; 0 ,
$$

and it drives $F = 0$ with Newton's method. **You supply only the algebraic
functions $f0_i$ and $f1_{ij}$**, evaluated at one point, given the local values
of $\mathbf{u}$, $\nabla\mathbf{u}$, the auxiliary data $\mathbf{a}$, and the
uniform constants. PETSc does the integral and the Newton loop.

Matching (W) term by term:

$$
\boxed{\;
f1_{ij} \;=\; \mathsf{M}_{ij} \;=\; 2\,\eta H\,\mathsf{D}_{ij} \;=\; 2N\,\mathsf{D}_{ij}
\;}
\qquad
\boxed{\;
f0_i \;=\; \beta(|\mathbf{u}|)\,u_i \;-\; \tau_{d,i}
\;}
$$

Written out (with $u_x = \partial_1 u_1$, $u_y = \partial_2 u_1$,
$v_x = \partial_1 u_2$, $v_y = \partial_2 u_2$):

| PETSc slot | meaning | value |
| --- | --- | --- |
| `f0[0]` | multiplies $\phi_1$ | $\beta\,u - \tau_{d,x}$ |
| `f0[1]` | multiplies $\phi_2$ | $\beta\,v - \tau_{d,y}$ |
| `f1[0]` = `f1[x,x]` | multiplies $\partial_1\phi_1$ | $2N\,(2u_x + v_y)$ |
| `f1[1]` = `f1[x,y]` | multiplies $\partial_2\phi_1$ | $2N\,\tfrac12(u_y + v_x) = N(u_y+v_x)$ |
| `f1[2]` = `f1[y,x]` | multiplies $\partial_1\phi_2$ | $N(u_y+v_x)$ |
| `f1[3]` = `f1[y,y]` | multiplies $\partial_2\phi_2$ | $2N\,(2v_y + u_x)$ |

In the code:

```fortran
! SSA_FEM_PETSc_f1
call SSA_FEM_PETSc_strain( u_x_values, eps0, n, Abar, D, eps2, eta)   ! D, eps_e^2, eta
N = eta * a_values(i_H)
do m = 1, 4
  f1_values(m) = 2._c_double * N * D(m)
end do

! SSA_FEM_PETSc_f0
call SSA_FEM_PETSc_sliding_beta( u_values(1), u_values(2), a_values(i_tauc), c_values, beta, dbeta_duabs)
f0_values(1) = beta * u_values(1) - a_values(i_taudx)
f0_values(2) = beta * u_values(2) - a_values(i_taudy)
```

$\bar A$, $H$, $\tau_c$, $\tau_{d}$ arrive in the **auxiliary field** `a` (values
of a P1 data field interpolated to the quadrature point); $\dot\varepsilon_0^2$,
$n$ and the sliding parameters arrive in the **constants** array.


## 5. PETSc's pointwise Jacobian convention -> `g0`, `g3`

Newton's method needs the derivative of the residual. PETSc writes the Jacobian
bilinear form (test $\boldsymbol\phi$, trial $\boldsymbol\psi$) as

$$
J[\boldsymbol\phi,\boldsymbol\psi] =
\int_\Omega \Bigl[
  \phi_i\, g0_{ij}\, \psi_j
  + \phi_i\, g1_{ijl}\, \partial_l\psi_j
  + \partial_k\phi_i\, g2_{ijk}\, \psi_j
  + \partial_k\phi_i\, g3_{ijkl}\, \partial_l\psi_j
\Bigr]\mathrm{d}\Omega ,
$$

with

$$
g0_{ij} = \frac{\partial f0_i}{\partial u_j},\quad
g1_{ijl} = \frac{\partial f0_i}{\partial(\partial_l u_j)},\quad
g2_{ijk} = \frac{\partial f1_{ik}}{\partial u_j},\quad
g3_{ijkl} = \frac{\partial f1_{ik}}{\partial(\partial_l u_j)} .
$$

For the SSA:

* $f0_i = \beta(|\mathbf{u}|)u_i - \tau_{d,i}$ depends on $\mathbf{u}$ but **not**
  on $\nabla\mathbf{u}$ $\;\Rightarrow\; g1 = 0$.
* $f1_{ik} = 2\eta(\nabla\mathbf{u})\,H\,\mathsf{D}_{ik}(\nabla\mathbf{u})$ depends
  on $\nabla\mathbf{u}$ but **not** on $\mathbf{u}$ $\;\Rightarrow\; g2 = 0$.

So only `g0` and `g3` are provided (`petsc_ds_set_jacobian(ds, 0, 0, g0, NULL, NULL, g3)`).


## 6. Deriving `g3` (membrane-stress Jacobian)

We need $\displaystyle g3_{ijkl} = \frac{\partial f1_{ik}}{\partial g_{jl}}$, where
$g_{jl} \equiv \partial_l u_j$ is a component of the velocity gradient. Since
$f1_{ik} = 2H\,\eta\,\mathsf{D}_{ik}$ and *both* $\eta$ and $\mathsf{D}$ depend on
$g$, the product rule gives two terms:

$$
\frac{\partial f1_{ik}}{\partial g_{jl}}
= 2H\left[
   \frac{\partial \eta}{\partial g_{jl}}\,\mathsf{D}_{ik}
   \;+\;
   \eta\,\frac{\partial \mathsf{D}_{ik}}{\partial g_{jl}}
\right].
$$

### 6a. The "frozen-viscosity" term $\;\eta\,\partial\mathsf{D}/\partial g$

From $\mathsf{D}_{ik} = \tfrac12(g_{ik} + g_{ki}) + \delta_{ik}\,(g_{11}+g_{22})$,

$$
\frac{\partial \mathsf{D}_{ik}}{\partial g_{jl}}
= \tfrac12\bigl(\delta_{ij}\delta_{kl} + \delta_{il}\delta_{kj}\bigr)
  \;+\; \delta_{ik}\,\delta_{jl} .
$$

This is a **constant** $4\times4$ table (independent of the solution). In the code
it is the parameter `dD_dgradu(m,k)`, with $m\leftrightarrow(i,k)$ and
$k\leftrightarrow(j,l)$ in the layout
$(1,2,3,4) = (xx, xy, yx, yy)$ for $\mathsf{D}$ and
$(1,2,3,4) = (\partial u/\partial x,\ \partial u/\partial y,\ \partial v/\partial x,\ \partial v/\partial y)$
for $g$:

$$
\texttt{dD\_dgradu} =
\begin{pmatrix}
2 & 0 & 0 & 1\\
0 & \tfrac12 & \tfrac12 & 0\\
0 & \tfrac12 & \tfrac12 & 0\\
1 & 0 & 0 & 2
\end{pmatrix}.
$$

Its contribution to `g3` is $2H\eta\;\texttt{dD\_dgradu} = 2N\;\texttt{dD\_dgradu}$.
(This is exactly the Jacobian used before Phase 3, when $\eta$ was frozen.)

### 6b. The shear-thinning term $\;(\partial\eta/\partial g)\,\mathsf{D}$

Let $q \equiv \dfrac{1-n}{2n}$ (so $q = -\tfrac13$ for $n=3$). Then

$$
\frac{\partial \eta}{\partial \dot\varepsilon_e^2}
= \tfrac12 \bar A^{-1/n}\, q\,\bigl(\dot\varepsilon_e^2 + \dot\varepsilon_0^2\bigr)^{q-1}
= \frac{q\,\eta}{\dot\varepsilon_e^2 + \dot\varepsilon_0^2}.
$$

Now the key identity. Differentiating
$\dot\varepsilon_e^2 = g_{11}^2 + g_{22}^2 + g_{11}g_{22} + \tfrac14(g_{12}+g_{21})^2$:

$$
\frac{\partial \dot\varepsilon_e^2}{\partial g_{jl}} \;=\; \mathsf{D}_{jl} .
$$

*(Check: $\partial/\partial g_{11} = 2g_{11}+g_{22} = 2u_x+v_y = \mathsf{D}_{xx}$;
$\partial/\partial g_{12} = \tfrac12(g_{12}+g_{21}) = \mathsf{D}_{xy}$; etc.)* The
derivative of the effective strain rate squared with respect to the velocity
gradient **is the SSA strain tensor** $\mathsf{D}$ itself.

Chain rule:

$$
\frac{\partial \eta}{\partial g_{jl}}
= \frac{\partial \eta}{\partial \dot\varepsilon_e^2}\,
  \frac{\partial \dot\varepsilon_e^2}{\partial g_{jl}}
= \frac{q\,\eta}{\dot\varepsilon_e^2 + \dot\varepsilon_0^2}\,\mathsf{D}_{jl},
$$

so this term's contribution to `g3` is

$$
2H\,\frac{\partial \eta}{\partial g_{jl}}\,\mathsf{D}_{ik}
= \underbrace{\frac{2Hq\,\eta}{\dot\varepsilon_e^2 + \dot\varepsilon_0^2}}_{\displaystyle \texttt{coef}}
  \;\mathsf{D}_{ik}\,\mathsf{D}_{jl}
= \texttt{coef}\;\mathsf{D}(m)\,\mathsf{D}(k) .
$$

It is a **rank-1 (outer-product) update** built from $\mathsf{D}$. Because $q<0$
it is negative - this is the term that can make the full Newton Jacobian
indefinite far from the solution (SNES's line search handles that).

### 6c. Result

$$
\boxed{\;
g3_{ijkl}
= 2N\left[\tfrac12(\delta_{ij}\delta_{kl}+\delta_{il}\delta_{kj}) + \delta_{ik}\delta_{jl}\right]
  \;+\;
  \frac{2Hq\,\eta}{\dot\varepsilon_e^2+\dot\varepsilon_0^2}\;\mathsf{D}_{ik}\,\mathsf{D}_{jl}
\;}
$$

In the code (`SSA_FEM_PETSc_g3`), with `m` the $(i,k)$ (test) index of `D` and
`k_idx` the $(j,l)$ (trial) index:

```fortran
p    = (1.0 - n_glen) / (2.0 * n_glen)          ! q
coef = 2.0 * a_values(i_H) * eta * p / eps2      ! eps2 already includes eps0
...
g3(idx) = 2.0 * N * dD_dgradu(m, k_idx)  +  coef * D(m) * D(k_idx)
```


## 7. Deriving `g0` (basal-drag Jacobian)

We need $g0_{ij} = \partial f0_i/\partial u_j$ with
$f0_i = \beta(|\mathbf{u}|)\,u_i - \tau_{d,i}$ ($\tau_d$ is independent of
$\mathbf{u}$):

$$
g0_{ij}
= \beta\,\frac{\partial u_i}{\partial u_j}
  + u_i\,\frac{\partial \beta}{\partial u_j}
= \beta\,\delta_{ij}
  + u_i\,\frac{\mathrm{d}\beta}{\mathrm{d}|\mathbf{u}|}\,\frac{\partial |\mathbf{u}|}{\partial u_j}.
$$

With $|\mathbf{u}| = \sqrt{\delta_v^2 + u_1^2 + u_2^2}$ we have
$\partial|\mathbf{u}|/\partial u_j = u_j/|\mathbf{u}|$, hence

$$
\boxed{\;
g0_{ij} \;=\; \beta\,\delta_{ij}
\;+\; \frac{1}{|\mathbf{u}|}\frac{\mathrm{d}\beta}{\mathrm{d}|\mathbf{u}|}\; u_i\,u_j
\;}
$$

again $\beta\mathsf{I}$ plus a rank-1 update, this time built from $\mathbf{u}$.

For the Zoet-Iverson $\beta = \tau_c\,|\mathbf{u}|^{a}\,(|\mathbf{u}|+u_t)^{b}$
with $a = 1/p - 1$, $b = -1/p$:

$$
\frac{\mathrm{d}\beta}{\mathrm{d}|\mathbf{u}|}
= \tau_c\Bigl[a\,|\mathbf{u}|^{a-1}(|\mathbf{u}|+u_t)^{b}
           + b\,|\mathbf{u}|^{a}(|\mathbf{u}|+u_t)^{b-1}\Bigr]
= \beta\left(\frac{a}{|\mathbf{u}|} + \frac{b}{|\mathbf{u}|+u_t}\right).
$$

In the code (`SSA_FEM_PETSc_sliding_beta` then `SSA_FEM_PETSc_g0`):

```fortran
beta        = tauc * uabs**aexp * (uabs + ZIut)**bexp        ! aexp = 1/p-1, bexp = -1/p
dbeta_duabs = beta * (aexp/uabs + bexp/(uabs + ZIut))
...
rank1  = dbeta_duabs / uabs
g0(1)  = beta + rank1 * u1*u1          ! g0_11
g0(2)  =        rank1 * u1*u2          ! g0_12
g0(3)  =        rank1 * u2*u1          ! g0_21
g0(4)  = beta + rank1 * u2*u2          ! g0_22
```

(If $\beta$ hits the cap $\beta_{\max}$, $\mathrm{d}\beta/\mathrm{d}|\mathbf{u}|$
is set to zero so the residual and Jacobian stay consistent.)


## 8. Index bookkeeping: flat storage

PETSc passes these tensors as **flat C arrays**. The ordering is the one trap.

* `f0` has length $N_c$ (number of field components $=2$): `f0[i]`.
* `f1` has length $N_c\cdot\mathrm{dim} = 4$, ordered **[component][derivative]**:
  `f1[i*dim + j]` multiplies $\partial_j\phi_i$. So
  `f1[0..3]` = $(\,\mathsf{M}_{x,x},\ \mathsf{M}_{x,y},\ \mathsf{M}_{y,x},\ \mathsf{M}_{y,y}\,)$.
* `g0` has length $N_c^2 = 4$, ordered **[test-comp][trial-comp]**:
  `g0[i*Nc + j]`.
* `g3` has length $N_c^2\cdot\mathrm{dim}^2 = 16$, ordered
  **[test-comp][trial-comp][test-deriv][trial-deriv]**:

$$
\texttt{g3\_flat}\bigl[\,((i\,N_c + j)\,\mathrm{dim} + k)\,\mathrm{dim} + l\,\bigr]
\;=\;
\frac{\partial f1_{ik}}{\partial (\partial_l u_j)} ,
$$

and it multiplies $\partial_k\phi_i \cdot \partial_l\psi_j$ in the Jacobian form.
Note the interleaving: it is **not** "the two indices of `f1`, then the two
indices of $\nabla u$". The test-side pair is $(i,k)$ (split across positions 1
and 3) and the trial-side pair is $(j,l)$ (positions 2 and 4).

The code loop makes this explicit:

```fortran
do ci = 0,1                    ! i  = test component
  do di = 0,1                  ! k  = test derivative direction
    m = 2*ci + di + 1          !    -> row of D / dD_dgradu   (the (i,k) pair)
    do cj = 0,1                ! j  = trial component
      do dj = 0,1              ! l  = trial derivative direction
        k_idx = 2*cj + dj + 1  !    -> row of D / dD_dgradu   (the (j,l) pair)
        idx   = ((ci*2 + cj)*2 + di)*2 + dj + 1     ! = ((i*Nc+j)*dim+k)*dim+l + 1
        g3(idx) = 2*N*dD_dgradu(m,k_idx) + coef*D(m)*D(k_idx)
      end do
    end do
  end do
end do
```

Worked example: the entry $\partial f1_{y,y}/\partial(\partial u/\partial x)$ has
$i=1$ (=y, 0-based), $k=1$ (=y), $j=0$ (=x), $l=0$ (=x), so
$\texttt{idx}_0 = ((1\cdot2 + 0)\cdot2 + 1)\cdot2 + 0 = 10$, i.e. Fortran
`g3(11)`. From the formula: $2N\,\texttt{dD\_dgradu}(4,1) + \texttt{coef}\,\mathsf{D}(4)\mathsf{D}(1)
= 2N\cdot1 + \texttt{coef}\,\mathsf{D}_{yy}\mathsf{D}_{xx}$.


## 9. What PETSc does with these functions

You never write an integral or an element loop. Given the P1 velocity field,
PETSc (through `DMPlexSetSNESLocalFEM`):

1. loops over mesh triangles;
2. at each quadrature point of each triangle, interpolates $\mathbf{u}$,
   $\nabla\mathbf{u}$ and the auxiliary field $\mathbf{a}$ from the P1 nodal
   values, and calls `f0`/`f1` (for the residual) or `g0`/`g3` (for the Jacobian)
   with those pointwise values;
3. multiplies by the quadrature weights and the P1 basis functions / their
   gradients, forming a small element vector (residual) or element matrix
   (Jacobian);
4. scatters the element contributions into the global residual vector / Jacobian
   matrix;
5. SNES (Newton + line search) solves $J\,\delta\mathbf{u} = -F$, updates
   $\mathbf{u}$, and repeats until $\|F\|$ is below tolerance.

The auxiliary field $\mathbf{a}$ (thickness, flow factor, till yield stress,
driving stress) is a separate P1 field on a cloned mesh, filled once per solve
from the UFEMISM vertex arrays; the constants ($\dot\varepsilon_0^2$, $n$, the
sliding parameters) are uniform scalars set with `PetscDSSetConstants`.


## 10. Sanity checks

**Uniform flow / spatially constant coefficients.** If $\nabla\mathbf{u} = 0$
then $\mathsf{D} = 0$, so $f1 = 0$ and the membrane term vanishes. The residual
reduces to $\int (\beta u_i - \tau_{d,i})\phi_i = 0$ for all $\boldsymbol\phi$,
which forces $\beta\mathbf{u} = \boldsymbol\tau_d$ pointwise, i.e.

$$
\mathbf{u} = \boldsymbol\tau_d / \beta .
$$

On grounded ice ($\beta > 0$) this is a finite velocity in the direction of the
driving stress (downhill) - the correct plug-flow limit. (Phase 1 verified this
to machine precision, and this is why the sign in $f0 = \beta u - \tau_d$
matters: the earlier $+\tau_d$ gave $\mathbf{u} = -\boldsymbol\tau_d/\beta$, i.e.
ice flowing uphill.)

**Units.** With SI units and time in years: $[\eta] = \mathrm{Pa\,yr}$,
$[N] = [\eta H] = \mathrm{Pa\,yr\,m}$, $[\mathsf{D}] = \mathrm{yr^{-1}}$, so
$[f1] = \mathrm{Pa\,m}$ (a vertically integrated stress). $[\beta] =
\mathrm{Pa\,yr\,m^{-1}}$, $[u] = \mathrm{m\,yr^{-1}}$, $[\tau_d] = \mathrm{Pa}$,
so $[f0] = \mathrm{Pa}$. PETSc's residual $\int(\phi f0 + \nabla\phi\, f1)$ is
then dimensionally homogeneous ($\mathrm{Pa\,m^2}$ after the area integral, with
$\phi$ dimensionless).


## 11. Symbol - code - PETSc dictionary

| Symbol | Meaning | Code | Source |
| --- | --- | --- | --- |
| $u_1, u_2$ | velocity components | `u_values(1:2)` | PETSc `u` |
| $\partial_l u_k$ | velocity gradient, layout $(u_x, u_y, v_x, v_y)$ | `u_x_values(1:4)` | PETSc `u_x` |
| $\mathsf{D}$ | SSA strain tensor, layout $(xx, xy, yx, yy)$ | `D(1:4)` | `SSA_FEM_PETSc_strain` |
| $\dot\varepsilon_e^2 + \dot\varepsilon_0^2$ | regularised effective strain rate$^2$ | `eps2` | `SSA_FEM_PETSc_strain` |
| $\eta$ | effective viscosity | `eta` | `SSA_FEM_PETSc_strain` |
| $N = \eta H$ | viscosity $\times$ thickness | `N` | computed in `f1`/`g3` |
| $q = (1-n)/2n$ | $\mathrm{d}\ln\eta/\mathrm{d}\ln\dot\varepsilon_e^2$ | `p` | computed in `g3` |
| $\bar A$ | vert. avg. flow factor | `a_values(i_Abar)` | aux field slot 1 |
| $H$ | thickness ($\ge 0.1$ m) | `a_values(i_H)` | aux field slot 2 |
| $\tau_c$ | till yield stress $\times f_{gr}^{\,m}$ | `a_values(i_tauc)` | aux field slot 3 |
| $\tau_{d,x}, \tau_{d,y}$ | driving stress components | `a_values(i_taudx/i_taudy)` | aux field slots 4-5 |
| $\dot\varepsilon_0^2$ | strain-rate regularisation | `c_values(ic_eps0)` | `PetscDSSetConstants` |
| $n$ | Glen exponent | `c_values(ic_nglen)` | `PetscDSSetConstants` |
| $p, u_t, \delta_v, \beta_{\max}$ | Zoet-Iverson params | `c_values(ic_ZIp..ic_betamax)` | `PetscDSSetConstants` |
| $\beta(|\mathbf{u}|)$ | basal friction coefficient | `beta` | `SSA_FEM_PETSc_sliding_beta` |
| $\mathrm{d}\beta/\mathrm{d}|\mathbf{u}|$ | its speed derivative | `dbeta_duabs` | `SSA_FEM_PETSc_sliding_beta` |
| $f0_i$ | residual, undifferentiated part | `f0_values(1:2)` | `SSA_FEM_PETSc_f0` |
| $f1_{ij}$ | residual, gradient part (= $\mathsf{M}$) | `f1_values(1:4)` | `SSA_FEM_PETSc_f1` |
| $g0_{ij}$ | $\partial f0_i/\partial u_j$ | `g0(1:4)` | `SSA_FEM_PETSc_g0` |
| $g3_{ijkl}$ | $\partial f1_{ik}/\partial(\partial_l u_j)$ | `g3(1:16)` | `SSA_FEM_PETSc_g3` |


## References

* MacAyeal, D. R. (1989), Large-scale ice flow over a viscous basal sediment,
  *J. Geophys. Res.*, 94(B4), 4071-4087.
* Zoet, L. K. and Iverson, N. R. (2020), A slip law for glaciers on deformable
  beds, *Science*, 368, 76-78.
* PETSc manual, "PetscFE: Finite Element Discretizations" and the
  `PetscDSSetResidual` / `PetscDSSetJacobian` documentation.
* Reference example: `src/snes/tutorials/ex77.c` in PETSc (nonlinear elasticity),
  which uses the same `f0/f1/g0..g3` structure for a vector field with a
  non-linear material law.
