# Microphysics scaling law for Titan's fractal haze

Derivation of a closed law for the steady-state haze profile, in the coordinates
requested:

- `z` — altitude (positive **upward**)
- `N` — number of monomers in an aggregate (size coordinate)
- `C(z,N)` — concentration of aggregates of size `N` at altitude `z`
- `T` — temperature, `d` — monomer radius, `D` — fractal dimension

**Modeling choices** (set with the user):
1. **Monodisperse, self-similar** size closure: `C(z,N) = n(z) · δ(N − N̄(z))`.
   One characteristic size `N̄(z)` and number density `n(z)` per altitude.
2. **Full-transition (Fuchs) kernel** and **slip-corrected settling** — valid
   across the whole column (free-molecular at the top → continuum near surface).
3. **Sedimentation + eddy diffusion** `K(z)` for vertical transport (full BVP,
   §3 — not the advective-only limit).
4. **Two growth phases**: a compact **spherical** phase (`D=3`) for sub-monomer
   particles `N̄<1`, switching to the **fractal** phase (`D` as given, e.g. 2) at
   `N̄=1`.
5. Steady state (microphysics equilibrates within each radiative iteration).
6. Charge inhibition deferred (add `Q̄`, §7, once the base loop closes).

Physical ingredients from Burgalat & Rannou (2017) [BR17] and de Trenquelléon
(2025) [dT25]; see `physics_parameters.md` for symbols/values.

---

## 1. The monodisperse closure

Instead of resolving the full distribution `C(z,N)`, assume that at each altitude
the aggregates are characterized by a single mean monomer-count `N̄(z)`:

```
   C(z,N) = n(z) · δ( N − N̄(z) )
```

so the two unknown fields are the **number density** `n(z)` [m⁻³] and the
**mean size** `N̄(z)` [monomers]. We will write **two conservation laws**
(monomer mass, aggregate number) and close the system. This is exact for the
moments `M₀ = n` and `M₁ = nN̄`; it approximates the coagulation kernel by
evaluating it at the mean size (the standard monodisperse approximation, accurate
to a factor of order unity that can be absorbed into a tunable coefficient — §7).

Here `N̄` is a **continuous monomer-volume count**: it measures particle volume in
units of one monomer's volume, so `N̄<1` is meaningful — a sub-monomer sphere of
radius `r = d N̄^{1/3} < d`. This keeps mass bookkeeping exact through both
coalescence (spheres merging, `N̄<1`) and aggregation (sticking, `N̄≥1`).

Useful derived fields:
- monomer mass `m₁ = ρ_p (4/3)π d³`
- haze mass density `ρ_h(z) = m₁ · n N̄`
- mass radius `r(N̄) = d N̄^{1/3}`; mobility radius `r_a(N̄) = d N̄^{1/D(N̄)}`

**Two-phase fractal dimension** (decided with the user):

```
   D(N̄) = 3   for  N̄ < 1   (compact coalesced spheres; E = d^{0} = 1, r_a = r)
   D(N̄) = D   for  N̄ ≥ 1   (fractal aggregates, e.g. D = 2)
```

Coagulation of sub-monomer spheres conserves volume and produces larger spheres
(`D=3`) up to `N̄=1` (radius `d`, the monomer); beyond that, sticking builds open
aggregates (`D`). Eqs. (1)–(2) below use `D(N̄)`; the switch at `N̄=1` is continuous
in `r_a` (both give `r_a=d`).

---

## 2. Size-dependent transport coefficients (functions of N̄, z, T, d, D)

**Settling velocity** (BR17 Eq. 29, fractal + 1st-order Cunningham slip), with
`r = d N^{1/3}`, `E = d^{(D−3)/D}`:

```
                 2 ρ_p g(z) d²  ⎡              A λ(z)            ⎤
   ω(N,z)  =  ───────────────── ⎢ N^{(D−1)/D} + ────── N^{(D−2)/D}⎥          (1)
                   9 η(z)        ⎣                d               ⎦
```

`η(z)=η(T)`, `λ(z)=λ(T,P)`, `g(z)`, `A=1.591`, and `D=D(N̄)` (the two-phase value
above; `E=d^{(D−3)/D}`, so `E=1` in the spherical phase). Continuum (Kn≪1):
`ω∝N^{(D−1)/D}`; free-molecular (Kn≫1): `ω∝N^{(D−2)/D}`. In the spherical phase
(`D=3`) these are the familiar Stokes-with-slip terms `∝N^{2/3}` and `∝N^{1/3}`.

**Coagulation kernel**, equal-size self-collision `β(N̄,N̄;z)`, full transition via
the Fuchs harmonic bridge of the continuum (BR17 Eq. 10) and free-molecular
(BR17 Eq. 13) forms, with mobility radii `r_a = d N^{1/D}`:

```
   β_CO(N,N) = (2 k_B T / 3η)(r_ai+r_aj)(1/r_ai + 1/r_aj)  · [slip terms]      (2a)
   β_FM(N,N) = (6 k_B T / ρ_p)^{1/2} (r_ai+r_aj)² (r_i^{-3}+r_j^{-3})^{1/2}    (2b)
   β(N,N;z)  = β_CO β_FM / (β_CO + β_FM)                                       (2c)
```

Optionally multiply by the charge-inhibition factor `Q̄ = y/(e^y−1)` (BR17 Eq. 27).
All of `ω`, `β` are explicit functions of `(N̄, z, T, d, D)`.

---

## 3. The governing system (full BVP with eddy diffusion)

Work with two conserved densities: **aggregate number** `n(z)` [m⁻³] and
**monomer-volume density** `M(z) ≡ n N̄` [monomers m⁻³] (mass = `m₁ M`). For any
density `q` with downward settling speed `ω>0` and eddy diffusivity `K(z)`, the
steady 1-D balance (upward-positive flux convention) is

```
   d/dz [ ω(N̄,z) q + K(z) dq/dz ] = L_q(z) − S_q(z)                            (3)
```

where `ω q + K dq/dz` is the **downward flux** of `q`, `S_q` = production, `L_q` =
loss per unit volume. Applied to the two densities (with `N̄ = M/n`, `D=D(N̄)`):

```
   d/dz [ ω(N̄,z) M + K dM/dz ] = − S_M(z)                                      (A)
   d/dz [ ω(N̄,z) n + K dn/dz ] = ½ β(N̄,z) n² − S_n(z)                          (B)
```

- **(A) monomers** are conserved by coagulation (`L=0`), gained only from
  production `S_M`. Coalescence and aggregation both preserve total monomer
  volume, so this holds across both growth phases.
- **(B) aggregate number** is destroyed by self-coagulation at the Smoluchowski
  rate `½ β n²` (each collision removes one particle) and seeded by `S_n`.

This is a coupled **two-point boundary-value problem** (2nd order in each of
`M, n` → 4th order total) for `M(z), n(z)` on `z ∈ [0, z_top]`. The settling and
kernel are evaluated at the local mean size `N̄=M/n` with its two-phase `D(N̄)`.

**Production source** (internal, dT25 Eq. 5): with `G(z)` the normalized Gaussian
`exp[−(z−z₀)²/2Δz²]/(√(2π)Δz)`,

```
   S_M(z) = (𝒬_p / m₁) · G(z) ,   S_n(z) = S_M(z) / N_seed
   𝒬_p = 2.1×10⁻¹³ kg m⁻² s⁻¹ ,  z₀ ≈ 415 km ,  Δz ≈ 20 km ,  N_seed=(r_p/d)³
```

**Boundary conditions:**
- *Top* `z_top` (above the source): vanishing densities, `M, n → 0` (equivalently
  zero incoming flux from above).
- *Surface* `z=0`: deposition by settling only — zero diffusive flux,
  `dM/dz = dn/dz = 0`, so particles leave at the settling flux `ω q`.

Integrating `(A)` once over the source region recovers the intuitive first
integral: well below `z₀`, the **downward monomer flux is constant** and equals
the column production rate,

```
   ω(N̄,z) M + K dM/dz = P ,        P = 𝒬_p / m₁   [monomers m⁻² s⁻¹]           (A′)
```

---

## 4. Advection-only limit → the master scaling ODE (diagnostic / initial guess)

Setting `K→0` in §3 gives a clean closed-form law — useful as the **initial guess**
for the BVP solver and for physical intuition. With `K=0`, `(A′)` reduces to the
classic result that **haze mass density = mass flux / fall speed**:

```
   ω n N̄ = P    ⇒    n(z) = P / (ω N̄)    ⇒    ρ_h(z) = 𝒬_p / ω(N̄,z)            (4)
```

and substituting `ω n = P/N̄` into `(B)` (below the source, `S_n=0`) collapses the
system to a **single first-order ODE for the mean size**:

```
   ┌─────────────────────────────────────────────┐
   │   dN̄/dz  =  −  β(N̄,z) · P / ( 2 ω(N̄,z)² )   │                            (5)
   └─────────────────────────────────────────────┘
        N̄(z₀)=N_seed ;  integrate downward, switching D(N̄): 3→D at N̄=1
```

RHS < 0, so `N̄` **grows as the particle falls**. Then `n=P/(ωN̄)`, `r=dN̄^{1/3}`,
`r_a=dN̄^{1/D}`. Equation (5) is the scaling law in its sharpest form; the §3 BVP
restores eddy mixing around it.

---

## 5. Asymptotic power-law scalings (the "scaling-law" exponents)

Holding the background (`T, η, λ, g, ρ_p, P`) roughly constant over a scale
height, `β` and `ω` become pure powers of `N̄`, and (5) integrates to power laws.
For **equal-size** aggregates:

| Regime | `β(N̄) ∝` | `ω(N̄) ∝` | `dN̄/dz ∝ −N̄^s`, `s =` | `N̄(depth)` for D=2 |
|---|---|---|---|---|
| Continuum (Kn≪1, low alt.) | `N̄⁰` (size-indep., `8k_BT/3η`) | `N̄^{(D−1)/D}` | `−2(D−1)/D` | `s=−1 ⇒ N̄ ∝ depth^{1/2}` |
| Free-molecular (Kn≫1, high alt.) | `N̄^{2/D−1/2}` | `N̄^{(D−2)/D}` | `(6−2D)/D − 1/2` | `s=+1/2 ⇒ N̄ ∝ depth²` |

(`depth ≡ z₀−z`; exponents use the fractal `D`, e.g. 2.) So aggregates grow
**steeply just below the source** (FM, `N̄∝depth²`), then growth slows in the dense
lower atmosphere (continuum, `N̄∝depth^{1/2}`). The continuum size-independence of
`β` (BR17 Eq. 10 leading term `= 8k_BT/3η`) is the well-known Brownian result.
The brief **spherical phase** (`N̄<1`, `D=3`) just below `z₀` follows the same table
with `D=3` (e.g. FM exponent `s = (6−6)/3 − 1/2 = −1/2`). Equation (5) with the
full Fuchs kernel (2c) interpolates smoothly between these limits — no single
exponent, which is why we integrate (5)/§3 numerically rather than quoting one
power law.

---

## 6. Coupling into the radiative transfer (Step 3)

Per layer, the microphysics hands DISORT:
- number density `n(z)` and characteristic radius `r_a(z)=d N̄^{1/D}`;
- extinction `Δτ(λ,z) = n(z) · C_ext(r_a, λ) · Δz`, with `C_ext, ω₀, g` read from
  the aggregate-optics table (dT25 Eqs. 16–18; needs Tomasko 2008a optical
  constants). Aggregate projected area scales `∝ N̄^{2/D} d²`.

The loop: `T(z)` (DISORT) → `η(T), λ(T,P)` → integrate (5)/§3 → `n(z), N̄(z)` →
optics → `Δτ, ω₀, g` → DISORT → new `T(z)`, iterate to a fixed point.

---

## 7. Assumptions, accuracy, and open questions

**Built-in approximations**
- Monodisperse: replaces the true spread by `N̄`. Real distributions coagulate
  somewhat faster; the standard fix is a constant prefactor on `β` in (5)/(B)
  (O(1–2)). Can be calibrated against a full-PBE or two-moment run later.
- Evaluating `β` at the mean size (vs. averaging over the distribution).
- Single monodisperse mode carried through the spherical→fractal switch (rather
  than dT25's two coexisting modes); the switch at `N̄=1` is instantaneous.

**Resolved decisions**
- ✅ **Seed / sub-monomer phase:** carry the spherical phase explicitly — `D=3`
  for `N̄<1`, switching to fractal `D` at `N̄=1` (§1–§2). Seed at `N̄=N_seed=(r_p/d)³`.
- ✅ **Eddy diffusion:** included — the full §3 BVP is solved in
  `src/microphysics/bvp.py` (production imposed as a top flux BC at `z₀`, closed
  by settling-only deposition at the surface); eq. (5) is the `K→0` initial
  guess (`scaling_law.py`).
- ⏸ **Charge inhibition:** deferred. Add the `Q̄` factor (BR17 Eq. 27,
  `n_e≈15 e⁻/µm`) as a multiplier on `β` once the base loop closes.

**Cross-validation (implemented, `scripts/cross_validate.py`).** On the crude
Titan reference column the BVP reproduces the Tomasko (2008) haze extinction
scale height (**62 km** vs. 65 km observed; the `K→0` limit gives 54 km) and a
main-haze characteristic radius of **~0.3–0.5 µm** (cf. dT25 `r_c≈0.46 µm`).
Monomer mass flux is conserved to machine precision; the BVP and master ODE
agree to ~9% in the settling-dominated lower haze. The extinction scale height
is sensitive to `K(z)`, confirming it as the key uncertain input.

**Open inputs (sensible defaults assumed; flag if you disagree)**
- **Background profiles.** `η(T)`, `λ(T,P)` from N₂ kinetic theory; `g(z)` from
  Titan gravity; `K(z)` from a standard Titan eddy-diffusion profile (Vuitton
  2019, as in dT25). Swap in a preferred `K(z)` if you have one.
- **Surface BC.** Zero diffusive flux (settling-only deposition); switch to a
  prescribed deposition velocity if needed.
