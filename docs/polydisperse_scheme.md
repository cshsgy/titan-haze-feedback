# Polydisperse microphysics: bimodal 2-moment log-normal scheme

Upgrade of the monodisperse Step-2 model to a **two-mode, two-moment, log-normal**
aerosol microphysics, following Burgalat & Rannou (2017, *J. Aerosol Sci.* 105,
151; "Brownian coagulation of a bi-modal distribution of both spherical and
fractal aerosols") — the same scheme de Trenquelléon (2025) uses. Decisions:
2-moment **fixed-width** log-normal, **bimodal** (spheres + aggregates), with the
**polydisperse optics** re-coupled into Step 3.

## State

Two log-normal modes, each carried by two moments $M_n=\int r^n\,n(r)\,dr$ (over
bulk/mass-equivalent radius $r$): number $M_0$ and volume $M_3$.

| mode | $D_f$ | $E=r_m^{(D_f-3)/D_f}$ | fields |
|---|---|---|---|
| spherical (S) | 3 | 1 | $M_0^S, M_3^S$ |
| fractal (F)   | 2 | $r_m^{-1/2}$ | $M_0^F, M_3^F$ |

So **4 fields** per layer. Bulk radius $r=N^{1/3}r_m$; apparent/mobility radius
$r_a=E\,r^{3/D_f}$ (B&R Eq. 8). Mean radius per mode $r_0=(M_3/M_0)^{1/3}$.

## Log-normal closure (the inter-moment function $\alpha$)

B&R Eq. (3): $M_k/M_0 = r_0^k\,\alpha(k)$, with $\alpha(3)=1$. For a log-normal of
geometric std $\sigma$ (set $s=\ln\sigma$), expanding about $r_0=(M_3/M_0)^{1/3}$:

$$\boxed{\;\alpha(k)=\exp\!\big[\tfrac12\,(k^2-3k)\,s^2\big]\;}\qquad(\alpha(3)=1\ \checkmark,\ \alpha(0)=1).$$

*(Standard log-normal result; confirm convention vs Burgalat et al. 2014 Icarus 231,
which is where B&R import $\alpha$ from. $\sigma$ is a fixed per-mode parameter —
see "Open parameters".)*

Substitution identity (Eq. 2): every $\iint r_i^p r_j^q\,n_X n_Y$ reduces to
$M_p^X M_q^Y$, and any off-order $M_p = M_0\,r_0^p\,\alpha(p)$.

## Sedimentation: moment-averaged settling (B&R Eqs. 29–30)

Single particle (Stokes + first-order slip, $A=1.591$, $\lambda$=mean free path):
$$\omega(r)=\frac{2\rho_p g}{9\eta E}\Big[r^{(3D_f-3)/D_f}+\tfrac{A\lambda}{E}\,r^{(3D_f-6)/D_f}\Big].$$
Moment flux $\Phi_{M_k}=\int r^k n\,\omega\,dr=\frac{2\rho_p g}{9\eta E}\big[M_{(D_f(k+3)-3)/D_f}+\tfrac{A\lambda}{E}M_{(D_f(k+3)-6)/D_f}\big]$, so the
**effective settling velocity differs by moment** $\langle\omega\rangle_k=\Phi_{M_k}/M_k$.
Reducing with $\alpha$ (derived):

$$\langle\omega\rangle_0=\frac{2\rho_p g}{9\eta E}\Big[r_0^{3-3/D_f}\alpha(3-\tfrac3{D_f})+\tfrac{A\lambda}{E}\,r_0^{3-6/D_f}\alpha(3-\tfrac6{D_f})\Big],$$
$$\langle\omega\rangle_3=\frac{2\rho_p g}{9\eta E}\Big[r_0^{3-3/D_f}\alpha(6-\tfrac3{D_f})+\tfrac{A\lambda}{E}\,r_0^{3-6/D_f}\alpha(6-\tfrac6{D_f})\Big].$$

$\langle\omega\rangle_3>\langle\omega\rangle_0$ for $\sigma>1$ (mass settles faster
— gravitational sorting). This is what the monodisperse $\omega(\bar N)$ missed.

## Coagulation closure (B&R Eqs. 6–25)

Per regime $X\in\{$CO (continuum), FM (free-molecular)$\}$, the moment tendencies
are quadratic in the moments (Eqs. 21–24):

$$\dot M_0^S=(p_0^X-2)\gamma_0^{X,SS}(M_0^S)^2-\gamma_0^{X,SF}M_0^S M_0^F,\quad
\dot M_3^S=(p_3^X-1)\gamma_3^{X,SS}(M_3^S)^2-\gamma_3^{X,SF}M_3^S M_3^F,$$
$$\dot M_0^F=(1-p_0^X)\gamma_0^{X,SS}(M_0^S)^2+\gamma_0^{X,FF}(M_0^F)^2,\quad
\dot M_3^F=(1-p_3^X)\gamma_3^{X,SS}(M_3^S)^2+\gamma_3^{X,SF}M_3^S M_3^F.$$

- $\gamma_k^{X,YY}$: scalars assembled by inserting the regime-$X$ kernel power-law
  expansion into the integrals (Eqs. 6–7) and reducing every $\iint$ via Eqs. (2)+(3).
  **CO kernels:** Eqs. 10 (SS), 11 (FF), 12 (SF), $K_{CO}=2k_BT/3\eta$, slip const
  $C=A\lambda$. **FM kernels:** Eqs. 16–20, $K_{FM}=(6k_BT/\rho_p)^{1/2}$, with the
  $\sqrt{r_i^{-3}+r_j^{-3}}\approx b\,(r_i^{-3/2}+r_j^{-3/2})$ approximation.
- $b_k^{(n)}$ ($n=1..5$ for the five G-kernels in Table 2, $\times$ moment order
  $\{0,3\}$ = 10 numbers): the FM bridging factors, each the ratio of the exact FM
  integral to its power-law approximation (Eq. 15) — **computed numerically** over
  the log-normal (functions of $T,P,r_0,\sigma$; tabulate offline).
- **Regime bridge = harmonic mean (Eq. 25)**, per (mode, interaction):
  $\dot M\big|=\dot M\big|^{CO}\dot M\big|^{FM}/(\dot M\big|^{CO}+\dot M\big|^{FM})$,
  then sum interactions.

### Inter-mode transfer (Eq. 26)
$p_k$ = probability an S+S product stays spherical (coalesced radius $\le r_m$):
$p_k=\int_0^{r_m} r^k F_s/\int_0^\infty r^k F_s$, $F_s(r)=\tfrac12\int_0^r\beta_{SS}(u,r{-}u)n(u)n(r{-}u)du$
(volume-conserving pairing $r^3=u^3+(r{-}u)^3$). Precompute as a table in
$(T,P,r_0^S)$. S+S$\to$F feeds the F mode (the $(1-p_k)$ terms); S+F$\to$F always.
This replaces the monodisperse hard $\bar N=1$ switch.

### Production (Eqs. 31–32), into the S mode
$\dot M_3^S=(Q_0/\rho_p)\,G(z)$, $\dot M_0^S=\dot M_3^S/(r_C^3\,\alpha(3))$, with
$r_C=2\times10^{-8}$ m the production (monomer) radius and $G(z)$ the Gaussian.

## Transport BVP (steady 1-D, the architecture carries over)

For each of the 4 moments $M$ with its $\langle\omega\rangle$ and eddy $K$:
$$\partial_z\big[\langle\omega\rangle_M\,M + K\,\partial_z M\big]=\dot M_{\rm coag}+\dot M_{\rm prod},$$
a coupled 4-field two-point BVP (today's 2-field BVP generalized). $K\to0$ gives a
4-ODE master system for the initial guess. BCs: production as a top flux into S;
settling-only deposition at the surface.

## Optics (polydisperse, re-coupled into Step 3)

- **IR absorption (RDG):** $C_{\rm abs}=N\,C_{\rm abs}^{\rm mono}$, additive over
  monomers ⇒ total $=M_1^{\rm monomer\text{-}vol}\,C^{\rm mono}$, i.e. only the
  **volume moment** $M_3$ matters. *Unchanged by polydispersity* (mono already
  exact). The F mode dominates.
- **SW scattering / visible extinction:** needs the size distribution.
  $\int n(r)\,r_a^2 dr = M_{2\cdot 3/D_f}$-type moments per mode → reduce via
  $\alpha$. $\omega_0(\lambda),g(\lambda)$ integrated over both modes' log-normals.
- Re-run the Step-3 coupling with the bimodal haze and check whether the
  bistability persists/shifts.

## Build sequence

1. `microphysics/moments.py` — log-normal $\alpha(k)$, the per-mode/per-moment
   $\langle\omega\rangle_k$, the $\gamma$ coagulation coefficients (CO+FM), the
   $b_k$ FM factors, harmonic-mean bridge, transfer $p_k$. **Unit tests:** recover
   monodisperse as $\sigma\to1$; volume conserved ($\dot M_3^S+\dot M_3^F$ from
   coag $=0$); $\langle\omega\rangle_3>\langle\omega\rangle_0$.
2. `K\to0` 4-ODE master system (initial guess).
3. 4-field eddy-diffusion BVP (generalize `bvp.py`).
4. Polydisperse optics (extend `aggregate_optics.py` / `optics.py`).
5. Cross-validate (extinction scale height, sizes, S/F mass split) + re-run Step 3.

## Results (pieces 4-5; `scripts/polydisperse_compare.py`, `writing/figs/polydisperse_compare.png`)

At `sigma_S=1.5, sigma_F=2.0` (K->0 profiles):

- **Optics.** The bimodal haze is **~4.4x optically thinner** than monodisperse
  (visible column tau **1.9 vs 8.4**) and sits **higher** (centroid 128 vs 107 km).
  Cause: the broad aggregate distribution settles its mass ~5x faster --
  `<w>_3 ~ alpha(4.5) ~ 5x` the single-size speed (gravitational sorting). This is
  a first-order, sigma_F-dependent control on haze opacity the monodisperse model
  cannot represent; tau~2 is below the observed ~8, suggesting Titan's aggregate
  distribution is narrower than sigma_g=2.0 (or production is higher) -- a testable
  inference. **=> the sigma_F sweep is the key follow-up.**
- **Temperature.** The DISORT stratopause is **insensitive** to the 4.4x tau
  difference: ~142 K for both hazes (the lower stratosphere differs, the
  stratopause does not).
- **Bistability RESOLVED.** With the polydisperse haze the radiative-convective
  equilibrium is **monostable** (warm vs cool start: 142 vs 139 K, ~3 K), versus
  the ~31 K split with the monodisperse haze. The absorbing-haze bistability was
  tied to the monodisperse model's overestimate of the haze opacity at altitude;
  the realistic (thinner) polydisperse haze falls below the feedback threshold.
  This is sigma_F-dependent (a thicker haze at smaller sigma_F may restore it).

## Open parameters (need values)

- $\sigma_S,\sigma_F$: fixed geometric widths (Rannou et al. 2004 / Cours et al.
  2011 / Burgalat et al. 2014). Default to a literature value, expose as params.
- $Q_0$ vs the existing $Q_{\rm prod}$: reconcile units (B&R use
  $\mathrm{kg\,m^{-3}\,s^{-1}}$ volumetric; ours is a column rate).
- $r_C$ (production radius) vs the existing seed $r_p$.
