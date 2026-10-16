# GeneralizedSasakiNakamura.jl

![license](https://img.shields.io/github/license/ricokaloklo/GeneralizedSasakiNakamura.jl)
[![GitHub release](https://img.shields.io/github/v/release/ricokaloklo/GeneralizedSasakiNakamura.jl.svg)](https://github.com/ricokaloklo/GeneralizedSasakiNakamura.jl/releases)
[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://ricokaloklo.github.io/GeneralizedSasakiNakamura.jl/dev/)

GeneralizedSasakiNakamura.jl solves the frequency-domain radial Teukolsky equation of a Kerr black hole with the Generalized Sasaki-Nakamura (GSN) formalism. It handles *both ingoing and outgoing* radiation of scalar, electromagnetic and gravitational type (spin weight $s = 0, \pm 1, \pm 2$), accepts *complex* frequencies, and uses the $M = 1$ convention throughout. Beyond the source-free (homogeneous) equation, it computes the gravitational ($s = \pm 2$) waveform amplitudes and fluxes radiated by a test particle on a generic stable bound Kerr geodesic.

The angular Teukolsky equation is solved with the accompanying package [SpinWeightedSpheroidalHarmonics.jl](https://github.com/ricokaloklo/SpinWeightedSpheroidalHarmonics.jl) using a spectral decomposition method.

The paper describing both the GSN formalism and the implementation can be found in [2306.16469](https://arxiv.org/abs/2306.16469). A set of Mathematica notebooks deriving all the equations used in the code can be found in [10.5281/zenodo.8080241](https://doi.org/10.5281/zenodo.8080241).

**Contents**
- [What's new in v0.10.0](#whats-new-in-v0100)
- [Installation](#installation)
- [Homogeneous solvers](#homogeneous-solvers)
  - [Solver methods](#solver-methods)
  - [Basic radial-solution example](#basic-radial-solution-example)
  - [Complex frequencies, quasinormal modes, and excitation factors](#complex-frequencies-quasinormal-modes-and-excitation-factors)
- [Inhomogeneous solvers](#inhomogeneous-solvers)
  - [Point-particle sources and single-mode solutions](#point-particle-sources-and-single-mode-solutions)
  - [Total fluxes](#total-fluxes)
- [How to cite](#how-to-cite)
- [License](#license)

More worked examples are on the [Examples](https://ricokaloklo.github.io/GeneralizedSasakiNakamura.jl/dev/examples/) page of the documentation and in the notebooks in [`examples/`](https://github.com/ricokaloklo/GeneralizedSasakiNakamura.jl/tree/main/examples).

## What's new in v0.10.0

Version 0.10.0 adds a GSN-based ISEM (iterative series expansion matching) solver, which solves the generalized Sasaki–Nakamura equation directly in about a millisecond for real and complex frequencies, including exact extremality $a = \pm 1$, and a Kerr quasinormal-mode solver, `qnm`, which returns QNM frequencies for $|a| \le 1$ with their GSN amplitudes and excitation factors (frequency only for the zero-damping modes at $|a| = 1$). Radial solves, QNM calls, and point-particle mode and flux calls are thread-safe, and point-particle calls on the same orbit can share one trajectory through `with_pointparticle_submission` (see [Parallel Calls](https://ricokaloklo.github.io/GeneralizedSasakiNakamura.jl/dev/APIs/#Parallel-Calls)).

Earlier releases: v0.8.0 added waveform amplitudes and fluxes at infinity and at the horizon for a test particle on a _generic (eccentric, inclined) timelike bound orbit_, by solving the inhomogeneous Sasaki–Nakamura (SN) equation using integration by parts; v0.9.0 added the first ISEM solver and the high-level total-flux interface `Teukolsky_pointparticle_flux`.

## Installation
Julia 1.12 or later in the 1.x series is required.

To install the package using the Julia package manager, simply type the following in the Julia REPL:
```julia
using Pkg
Pkg.add("GeneralizedSasakiNakamura")
```

*Note: The dependencies, including [SpinWeightedSpheroidalHarmonics.jl](https://github.com/ricokaloklo/SpinWeightedSpheroidalHarmonics.jl) (v1.4 or later in the 1.x series) and KerrGeodesics.jl (v0.5), are installed automatically from the General registry. To load SpinWeightedSpheroidalHarmonics in your own code, also `Pkg.add("SpinWeightedSpheroidalHarmonics")`.*

## Homogeneous solvers

These solvers construct source-free radial solutions and their asymptotic amplitudes for four boundary conditions: `IN` (purely ingoing at the horizon), `UP` (purely outgoing at infinity), `OUT` (purely outgoing at the horizon) and `DOWN` (purely ingoing at infinity); `OUT` and `DOWN` are linear combinations of `IN` and `UP`. Pass the boundary condition as the sixth argument to get one solution, for example `Teukolsky_radial(-2, 2, 2, 0.7, 0.5, OUT)`; without it, `Teukolsky_radial` and `GSN_radial` return the pair `(IN, UP)` used below. Both accept real and complex frequencies.

### Solver methods

The `method` keyword selects the solver: `"ISEM"`, `"linear"`, `"Riccati"`, or the default `"auto"`, which runs ISEM and, if ISEM throws, logs a warning and falls back to numerical integration (Riccati, then linear, for real frequencies; linear, then Riccati, for complex ones). With an explicit `xm`, `"auto"` raises the ISEM error instead.

#### Numerical integration (`method = "linear"` or `"Riccati"`)

The original GSN solver numerically evolves the radial equation with the `linear` or `Riccati` method and attaches analytical boundary ansatzes near the horizon and infinity:

<p align="center">
  <img width="60%" src="https://github-production-user-asset-6210df.s3.amazonaws.com/55488840/248724944-9707332b-1238-4b3b-b1c0-ac426a1b3dc6.gif">
</p>

At low frequency these methods move their outer boundary outward, and below $|M\omega| = 10^{-8}$ they lose the `UP` solution; use ISEM there.

The on-the-fly benchmark against the Mathematica implementation of the Mano–Suzuki–Takasugi (MST) method is shown below:

<table>
  <tr>
    <th>GeneralizedSasakiNakamura.jl numerical solver</th>
    <th><a href="https://github.com/BlackHolePerturbationToolkit/Teukolsky">Teukolsky</a> Mathematica package using the MST method </th>
  </tr>
  <tr>
    <td><p align="center"><img width="100%" src="https://github-production-user-asset-6210df.s3.amazonaws.com/55488840/248965077-7d216deb-5bae-433f-a699-d40a35f0e35d.gif"></p></td>
    <td><p align="center"><img width="100%" src="https://github-production-user-asset-6210df.s3.amazonaws.com/55488840/248966033-9e7d8027-81ee-4762-98d9-0ad0a1c030ad.gif"></p></td>
  </tr>
</table>

*(No caching: the equation was solved on the fly. The left animation was made with an earlier version of [examples/realtime-demo.ipynb](https://github.com/ricokaloklo/GeneralizedSasakiNakamura.jl/blob/main/examples/realtime-demo.ipynb), which now runs ISEM; the MST side with [examples/realtime-visualization-MST.nb](https://github.com/ricokaloklo/GeneralizedSasakiNakamura.jl/blob/main/examples/realtime-visualization-MST.nb).)*

#### ISEM (`method = "ISEM"`)

ISEM constructs the GSN radial function $X$ from series expansions matched on the generalized Sasaki–Nakamura equation itself: local solutions at the horizon and at infinity (with small-frequency and large-frequency expansions where they apply) are continued by ordinary-point series and matched at an interior point. The Teukolsky function $R$ and the auxiliary function $Y$ used by the point-particle convolutions are transformations of that same solution:

<p align="center">
  <img width="80%" alt="ISEM matching" src="docs/src/isem_matching_original_30fps.gif">
</p>

`method = "ISEM"` never falls back. ISEM chooses its expansions, expansion order and matching point automatically, and the point-particle mode and flux calls use the same solver. Optional controls are the expansion order `N`, the tolerance `tol`, the small- and large-frequency expansion switches `sfe` and `lfe`, and the matching point `xm` in the compact coordinate $x = (r - r_+)/(r - r_-) \in (0, 1)$; see the [API reference](https://ricokaloklo.github.io/GeneralizedSasakiNakamura.jl/dev/APIs/). Leave `xm` unset: an explicit `xm` can introduce a cancellation error in the incidence amplitude that `tol` does not bound.

For $|a| < 1$, at the superradiant threshold $\omega = m a / (2 r_+)$ (exact equality only), the `UP` solution is always built with ISEM in a logarithmic horizon basis, and its incidence and reflection amplitudes are the coefficients of that basis. At $\omega = 0$ only `IN` and `UP` are available: they are the closed-form static solutions (Gauss hypergeometric functions), and their incidence and reflection amplitudes are `missing`.

The GSN-ISEM on-the-fly solve and evaluation speed, recorded with an interactive version of the demo that has a $\log_{10}\omega$ slider:

<p align="center">
  <img width="80%" alt="GSN-ISEM on-the-fly solve and evaluation" src="docs/src/isem_realtime_demo_30s.gif">
</p>

The Pluto notebook [examples/realtime-demo-pluto.jl](https://github.com/ricokaloklo/GeneralizedSasakiNakamura.jl/blob/main/examples/realtime-demo-pluto.jl) shows the same plots and re-solves on every control change (its $\omega$ slider is linear, from 0.05 to 1); [examples/realtime-demo.ipynb](https://github.com/ricokaloklo/GeneralizedSasakiNakamura.jl/blob/main/examples/realtime-demo.ipynb) draws a single frame.

*Change from v0.9.0:* in v0.10.0, `method = "ISEM"` solves the GSN equation directly. The v0.9.0 ISEM matched series solutions of the Teukolsky equation; it has been replaced, together with its `rhom`, `TSinInf`/`TSoutInf`/`TSinHor`/`TSoutHor` controls and its negative `xm` coordinate.

### Basic radial-solution example

The following snippet solves the (source-free) Teukolsky function for the mode $s=-2, \ell=2, m=2, a/M=0.7, M\omega=0.5$ with the purely-ingoing boundary condition at the horizon, $R^{\textrm{in}}$, and the purely-outgoing boundary condition at spatial infinity, $R^{\textrm{up}}$:
```julia
using GeneralizedSasakiNakamura # This is going to take some time to pre-compile, mostly due to DifferentialEquations.jl

# Specify which mode to solve
s=-2; l=2; m=2; a=0.7; omega=0.5;

# NOTE: julia uses 'just-ahead-of-time' compilation. Calling this the first time in each session will take some time
Rin, Rup = Teukolsky_radial(s, l, m, a, omega)
```
Displaying `Rin` in the REPL shows its asymptotic amplitudes:
```julia
julia> Rin
TeukolskyRadialFunction(
    mode = Mode(s = -2, l = 2, m = 2, a = 0.7, omega = 0.5, lambda = 1.6966094016353415),
    boundary_condition = IN,
    transmission_amplitude = 1.0 + 0.0im,
    incidence_amplitude = 6.536587661184236 - 4.941203897065118im,
    reflection_amplitude = -0.1282466191290958 - 0.4404813349644446im,
    normalization_convention = UNIT_TEUKOLSKY_TRANS
)
```

Perturbations go as $e^{-i\omega t + im\phi}$. With `UNIT_TEUKOLSKY_TRANS` the transmission amplitude is 1 and, with $p = \omega - m\Omega_H$, $\Omega_H = a/(2r_+)$, $\Delta = r^2 - 2r + a^2$ and the tortoise coordinate $r_*$,

```math
R^{\rm in} \to
\begin{cases}
\Delta^{-s} e^{-ipr_*}, & r \to r_+,\\
B^{\rm inc} r^{-1} e^{-i\omega r_*} + B^{\rm ref} r^{-2s-1} e^{+i\omega r_*}, & r \to \infty,
\end{cases}
\qquad
R^{\rm up} \to
\begin{cases}
C^{\rm inc} e^{+ipr_*} + C^{\rm ref} \Delta^{-s} e^{-ipr_*}, & r \to r_+,\\
r^{-2s-1} e^{+i\omega r_*}, & r \to \infty.
\end{cases}
```

`incidence_amplitude` and `reflection_amplitude` are $B^{\rm inc}$ and $B^{\rm ref}$ for `IN`, and $C^{\rm inc}$ and $C^{\rm ref}$ for `UP`; `lambda` is the separation constant $\lambda$ of the radial Teukolsky equation. See [2306.16469](https://arxiv.org/abs/2306.16469) for the full conventions.

`Rin(r)` evaluates $R^{\textrm{in}}$ at the Boyer–Lindquist radius $r$, here $r = 10M$:
```julia
julia> Rin(10)
77.57508416830544 - 429.4029095224015im
```

The GSN function $X$ is available in the same way through `GSN_radial`, but it is a function of the tortoise coordinate $r_*$, not of $r$; `rstar_from_r(a, r)` and `r_from_rstar(a, rstar)` convert between the two:
```julia
Xin, Xup = GSN_radial(s, l, m, a, omega)
Xin(rstar_from_r(a, 10))   # X at r = 10M; Xin(10) would be X at r_* = 10M
```
Its amplitudes use `UNIT_GSN_TRANS`: $X^{\rm in} \to e^{-ipr_*}$ at the horizon and $A^{\rm inc} e^{-i\omega r_*} + A^{\rm ref} e^{+i\omega r_*}$ at infinity; $X^{\rm up} \to e^{+i\omega r_*}$ at infinity and $C^{\rm ref} e^{-ipr_*} + C^{\rm inc} e^{+ipr_*}$ at the horizon. These GSN amplitudes are not the Teukolsky amplitudes $B$ and $C$ of `Rin` and `Rup`; the two differ by mode-dependent conversion factors.

`Y_radial(s, l, m, a, omega, bc)` returns the auxiliary function $Y(r)$ that the point-particle convolutions integrate against the source. It needs an explicit boundary condition (there is no `(Yin, Yup)` form), supports only $s = -2$ with `IN` and $s = +2$ with `UP`, and requires $\omega \neq 0$; see the [API reference](https://ricokaloklo.github.io/GeneralizedSasakiNakamura.jl/dev/APIs/).

### Complex frequencies, quasinormal modes, and excitation factors

Finding a quasinormal-mode (QNM) frequency is a root search in complex frequency: `qnm` computes the root and then evaluates the GSN scattering data at it. For the $a=0.68$, $s=-2$, $(\ell,m,n)=(2,2,0)$ mode, the ordinary ($\operatorname{Re}\omega > 0$) and mirror ($\operatorname{Re}\omega < 0$) branches are

```julia
ordinary_mode = qnm(0.68, -2, 2, 2, 0)           # (a, s, l, m, n): the spin comes first here
mirror_mode = qnm(0.68, -2, 2, 2, 0, mirror)
```

<details>
<summary>Output</summary>

```text
QuasiNormalMode(
    parameters = (a = 0.68, s = -2, l = 2, m = 2, n = 0)
    omega = 0.52397510429008398 - 0.081512623631199516im
    lambda = 1.6550030612535795 + 0.36026776076186151im
    incidence_amplitude = -2.9652484956329296e-15 + 1.0066747053005698e-14im
    reflection_amplitude = -0.80073261705109955 + 0.029136276807086189im
    incidence_derivative = -2.143807932934116 + 6.7875557748638862im
    excitation_factor = 0.01983394079017144 + 0.10427056860707475im
    formalism = GSN)

QuasiNormalMode(
    parameters = (a = 0.68, s = -2, l = 2, m = 2, n = 0)
    omega = -0.31116285768546015 - 0.088754663319965579im
    lambda = 5.4219094861401045 + 0.40958044727148446im
    incidence_amplitude = -5.1821813248855796e-15 + 4.7107003677217079e-15im
    reflection_amplitude = -1.2675830722586212 + 0.57434916613852249im
    incidence_derivative = 9.1832313896816231 - 11.990411121525646im
    excitation_factor = 0.13913753384546915 + 0.030226487430463019im
    formalism = GSN)
```

</details>

Check `.status` (`:accepted`, `:estimated` or `:failed`) first, then read the frequency from `.omega`. `qnm` returns a `QNMResult`, displayed as `QuasiNormalMode(...)` above, when every check passes; a `QNMEstimate` (`QuasiNormalModeEstimate(...)`) when the values are finite but one check missed its target, named in `limiting_gate`; a `QNMFailure` when a step fails (for example at the algebraically special frequency $-2i$ at $a = 0$, where the excitation factor is undefined); and, for a zero-damping mode (ZDM) at $a = \pm 1$, a `QNMEndpointResult` that carries the frequency but no amplitudes. At $|a| = 1$ a `QNMResult` displays as `ExtremalQuasiNormalMode(...)` and a `QNMEstimate` as `ExtremalQuasiNormalModeEstimate(...)`, both with `family = DM` (damped mode); a `QNMEndpointResult` displays as `ExtremalQuasiNormalMode(...)` with `family = ZDM`.

Pass the computed frequency itself, not a rounded copy, to the radial interfaces:

```julia
Rin, Rup = Teukolsky_radial(-2, 2, 2, 0.68, ordinary_mode.omega)
```

```julia
julia> Rup
TeukolskyRadialFunction(
    mode = Mode(s = -2, l = 2, m = 2, a = 0.68, omega = 0.523975104290084 - 0.08151262363119952im, lambda = 1.6550030612535795 + 0.3602677607618615im),
    boundary_condition = UP,
    transmission_amplitude = 1.0 + 0.0im,
    incidence_amplitude = -1.6277792171253123e-15 + 2.3406653193234936e-16im,
    reflection_amplitude = 1.101161565720053 + 2.1300599921594445im,
    normalization_convention = UNIT_TEUKOLSKY_TRANS
)
```

At a QNM frequency the incidence amplitude vanishes (here it is about $10^{-15}$): `Rup` is then also purely ingoing at the horizon, so it is proportional to `Rin`.

For a QNM frequency $\omega_q$, let $A^{{\rm GSN},{\rm inc}}$ and $A^{{\rm GSN},{\rm ref}}$ be the amplitudes $A^{\rm inc}$ and $A^{\rm ref}$ of $X^{\rm in}$ above (the result's `incidence_amplitude` and `reflection_amplitude`). The result reports $\alpha$ as `incidence_derivative` and the excitation factor $B$ as `excitation_factor`:

```math
\alpha^{{\rm GSN}}_{\ell m\omega_q}
=\left.\frac{\partial A^{{\rm GSN},{\rm inc}}_{\ell m\omega}}
{\partial\omega}\right|_{\omega=\omega_q},
\qquad
B^{{\rm GSN}}_{\ell m\omega_q}
=\frac{A^{{\rm GSN},{\rm ref}}_{\ell m\omega_q}}
{2\omega_q\alpha^{{\rm GSN}}_{\ell m\omega_q}}.
```

`ordinary` is the default branch. `mirror` gives the mode with negative real frequency: by default its frequency is the reflection $-\bar\omega$ of the ordinary $(s, \ell, -m, n)$ root at the same spin (polished on the negative-real branch from that reflection when the ordinary root is not accepted and did not stop at the negative imaginary axis), and its amplitudes and excitation factor are evaluated at that frequency; pass `mirror_guess` to polish the mirror root from your own starting value instead. Use `detailed=true` only when you need the Teukolsky-normalized values (`result.teukolsky.excitation_factor`, `.incidence_derivative`, `.reflection_amplitude`) or the callable `IN` solutions `X(rstar)`, `Y(r)` and `R(r)` (`Y` uses `UP` for $s = +2$ and is `missing` for $|s| \neq 2$).

When only the frequency is needed, `qnm_frequency` takes the same five numbers and returns the root alone as a `LeaverResult` (no amplitudes): check that `root.status` is `:accepted` and read `root.omega`. `qnm_sequence` follows one mode along a list of spins and returns a vector of `LeaverResult`, stopping after the first spin whose root is not accepted:

```julia
root = qnm_frequency(0.68, -2, 2, 2, 0)                    # (a, s, l, m, n)
root.omega
roots = qnm_sequence([0.0, 0.3, 0.68], -2, 2, 2, 0)        # (spins, s, l, m, n)
[r.omega for r in roots]
```

Overtone labels (`convention`): the default `:overtone` follows the Python [`qnm`](https://github.com/duetosymmetry/qnm) package (Stein 2019) and the Cook–Zalutskiy catalogue: `n` is the Kerr continuation of the `n`-th Schwarzschild overtone. `:complete_spectrum` (Forteza & Mourier, [arXiv:2107.11829](https://arxiv.org/abs/2107.11829)) numbers every root in order of damping. The two differ only for $s = -2$, $\ell = 2$, $|m| = 2$: for co-rotating modes ($a m > 0$ on the ordinary branch, $a m < 0$ on the mirror branch) the Schwarzschild frequency $M\omega = -2i$ (overtone 8) splits at $a \neq 0$ into Cook's $8_0$ (less damped) and $8_1$: both conventions label $8_0$ as $n = 8$, and `:complete_spectrum` labels $8_1$ as $n = 9$, a root that `:overtone` leaves unnumbered. For co- and counter-rotating modes alike, `:complete_spectrum` $n \ge 10$ is `:overtone` $n - 1$. At $a = 0$, $n = 8$ of every $|s| = 2$, $\ell = 2$ mode, and `:complete_spectrum` $n = 9$ of $s = -2$, $\ell = 2$, $|m| = 2$, is exactly $-2i$, where `qnm` returns a `QNMFailure`. For $\ell = 2$ and $n \ge 8$ use $s = -2$. A third value, `:leaver`, makes $n$ the continued-fraction inversion index only: it labels no mode and requires a `guess`.

Segments (`segment`): a few sequences (e.g. $\ell = 2$, $m = 0$, $9 \le n \le 26$) end on the negative imaginary axis at finite spin and reappear at larger spin (Cook's $n_0$ and $n_1$). `segment = 0` (default) is the piece connected to $a = 0$; beyond its end it fails with `stop_reason = :negative_imaginary_axis_endpoint`. `segment = 1` is the reappearing piece ($m = 0$ only), and `segment = :auto` whichever piece exists at $a$; for example `qnm(0.5, -2, 2, 0, 10; segment=1)`. The [QNM documentation](https://ricokaloklo.github.io/GeneralizedSasakiNakamura.jl/dev/QNM/) has the label table and the segment details.

## Inhomogeneous solvers

These solvers compute the radiation from a point particle on a bound Kerr geodesic: single-mode amplitudes and fluxes at infinity and at the horizon, and total fluxes summed over modes.

### Point-particle sources and single-mode solutions

A source is specified by the black-hole spin $a/M$ and the orbit's semi-latus rectum $p$, eccentricity $e$ and inclination parameter $x = \cos\theta_{\rm inc}$. A mode is labeled by $(s, \ell, m, n, k)$, where $n$ and $k$ are the radial and polar harmonic indices; $s = -2$ gives the solution at infinity and $s = +2$ the solution at the horizon.

#### Amplitudes at infinity

For the $s = -2$, $\ell = m = 2$, $n = k = 0$ mode driven by a test particle on a bound geodesic with $a/M = 0.9, p = 6M, e = 0.7, x = \cos(\pi/4)$:
```julia
mode_info = Teukolsky_pointparticle_mode(-2, 2, 2, 0, 0, 0.9, 6, 0.7, cos(π/4))   # (s, l, m, n, k, a, p, e, x)
```
Displaying `mode_info` gives
```julia
julia> mode_info
TeukolskyPointParticleMode(
    mode = Mode(s = -2, l = 2, m = 2, n = 0, k = 0, a = 0.9, omega = 0.06568724726732753, lambda = 3.606789012119982),
    amplitude_inf = 0.00023429507956768906 - 6.558414441160938e-5im,
    energy_flux_inf = 1.0917330113048546e-6,
    angular_momentum_flux_inf = 3.324033375494718e-5,
    Carter_const_flux_inf = 5.890504443058139e-5,
    method = (method = "isem_trapezoidal", radial_method = "ISEM", radial_sfe = false, N = 256, K = 64, truncation_floor = 1.0e-16),
)
```
The `_inf`/`_hor` suffixes appear only in the display: the fields are `mode_info.amplitude`, `.energy_flux`, `.angular_momentum_flux` and `.Carter_const_flux`, together with `.mode`, `.method`, `.trajectory`, `.Y_solution` and `.SWSH`. The default `method = "auto"` runs `"isem_trapezoidal"`: ISEM radial solutions and a trapezoidal convolution on a sampling grid that starts at the recorded `N × K = 256 × 64` and is refined by doubling, up to `2^12 × 2^9`, until the mode's energy flux has converged to the relative tolerance `min(sample_tol, 3√tol)` (`3e-4` with the default keywords `sample_tol = 1e-3` and `tol = 1e-8`) or falls below `truncation_floor`. `radial_sfe` records whether the small-frequency expansion of the ISEM radial solve was switched off (`false`) or left to ISEM (`:auto`).

`mode_info.amplitude` is $Z^{\infty}_{\ell m n k}$, the amplitude of the inhomogeneous radial Teukolsky solution near infinity at that frequency:
```julia
julia> mode_info.amplitude
0.00023429507956768906 - 6.558414441160938e-5im
```

#### Amplitudes at the horizon

For the same parameters, change the sign of $s$ to $2$:
```julia
mode_info = Teukolsky_pointparticle_mode(2, 2, 2, 0, 0, 0.9, 6, 0.7, cos(π/4))
```
gives
```julia
julia> mode_info
TeukolskyPointParticleMode(
    mode = Mode(s = 2, l = 2, m = 2, n = 0, k = 0, a = 0.9, omega = 0.06568724726732753, lambda = -0.39321098788001757),
    amplitude_hor = 0.006089946888790668 - 0.0014130019665083518im,
    energy_flux_hor = -2.843814878429853e-9,
    angular_momentum_flux_hor = -8.658651402627282e-8,
    Carter_const_flux_hor = -1.5343956812851395e-7,
    method = (method = "isem_trapezoidal", radial_method = "ISEM", radial_sfe = false, N = 256, K = 64, truncation_floor = 1.0e-16),
)
```
Here `mode_info.amplitude` is $Z^{\mathrm{H}}_{\ell m n k}$, the amplitude near the horizon:
```julia
julia> mode_info.amplitude
0.006089946888790668 - 0.0014130019665083518im
```

The particle mass $\mu$ is not an input: amplitudes are per unit $\mu$ and fluxes per unit $\mu^2$. The fluxes are time-averaged rates, obtained from the amplitudes as

```math
\dot E^{\infty} = \frac{|Z^{\infty}|^2}{4\pi\omega^2},
\qquad
\dot E^{\rm H} = \frac{\omega\,|Z^{\rm H}|^2}{64\pi\kappa\,(2r_+)^3(\kappa^2 + 4\epsilon^2)},
\qquad
\dot L_z = \frac{m}{\omega}\,\dot E,
```

with $\kappa = \omega - m\Omega_H$ (written $p$ in the radial asymptotics above; in this section $p$ is the semi-latus rectum) and $\epsilon = \sqrt{1-a^2}/(4r_+)$. Horizon fluxes are the rates absorbed by the black hole, so they are negative for superradiant modes ($\omega\kappa < 0$), as in the example above. These formulas fix the normalization of $Z^{\infty}$ and $Z^{\rm H}$.

### Total fluxes

`Teukolsky_pointparticle_flux(a, p, e, x; tol = 1e-8)` sums the mode fluxes to the relative truncation tolerance `tol`. It chooses the summation from the orbit: circular equatorial for `e = 0`, `x = ±1` (no convolution integral), eccentric equatorial for `x = ±1`, inclined (spherical) for `e = 0`, and generic otherwise, which performs two-dimensional convolution integrals and is much slower. It accepts stable bound orbits and marginally stable (critical) bound circular orbits, as classified by KerrGeodesics; for other $(a, p, e, x)$ it warns and returns `nothing`.

Each displayed line is a field of the result (`flux.infinity_energy_flux`, `flux.cost`, …; the `orbital_parameters` line shows `a`, `p`, `e`, `x`), except the `l_reached`, `k_reached_*`, `n_reached` and `convolution_integral` lines, which are stored in `flux.reached` (`flux.reached.l_reached_inf`, `flux.reached.n_reached_hor`, …); `flux.result` holds the raw summation output. `cost` is the wall-clock time of the call in seconds: on the machine that produced these examples, about 0.6–1 ms per mode (`cost / total_modes`) and about 22 minutes for the generic orbit. The first call in a session also compiles the code and takes longer.

<details>
<summary>Circular equatorial orbit</summary>

```julia
julia> flux = Teukolsky_pointparticle_flux(0.9, 6.0, 0.0, 1.0; tol=1e-8)
TeukolskyPointParticleFlux(
    orbital_parameters(a = 0.9, p = 6.0, e = 0.0, x = 1.0),
    orbit_type = circular,
    infinity_energy_flux = 0.0005658659548571375,
    infinity_angular_momentum_flux = 0.008825776472648025,
    infinity_carter_constant_flux = 0.0,
    horizon_energy_flux = -4.177363290666302e-6,
    horizon_angular_momentum_flux = -6.515407815579624e-5,
    horizon_carter_constant_flux = 0.0,
    total_modes = 378,
    l_reached = (infinity = 19, horizon = 12),
    convolution_integral = (strategy = "no convolution integral is needed",),
    tolerance = 1.0e-8,
    truncation_floor = (infinity = 1.0e-16, horizon = 1.0e-16),
    cost = 0.3374300003051758 seconds,
)
```

</details>

<details>
<summary>Spherical (inclined circular) orbit</summary>

```julia
julia> flux = Teukolsky_pointparticle_flux(0.9, 6.0, 0.0, 0.5; tol=1e-8)
TeukolskyPointParticleFlux(
    orbital_parameters(a = 0.9, p = 6.0, e = 0.0, x = 0.5),
    orbit_type = inclined,
    infinity_energy_flux = 0.0006637621566020543,
    infinity_angular_momentum_flux = 0.0060183149433979885,
    infinity_carter_constant_flux = 0.04294455189603736,
    horizon_energy_flux = -3.086574922405825e-6,
    horizon_angular_momentum_flux = -0.00010495330048599916,
    horizon_carter_constant_flux = 8.871154473181861e-5,
    total_modes = 11466,
    k_reached_inf = 15,
    k_reached_hor = 13,
    convolution_integral = (strategy = "ISEM adaptive trapezoidal for all modes",),
    tolerance = 1.0e-8,
    truncation_floor = (infinity = 1.0e-16, horizon = 1.0e-16),
    cost = 10.32890009880066 seconds,
)
```

</details>

<details>
<summary>Eccentric equatorial orbit</summary>

```julia
julia> flux = Teukolsky_pointparticle_flux(0.9, 6.0, 0.7, 1.0; tol=1e-8)
TeukolskyPointParticleFlux(
    orbital_parameters(a = 0.9, p = 6.0, e = 0.7, x = 1.0),
    orbit_type = eccentric,
    infinity_energy_flux = 0.0007457787885708398,
    infinity_angular_momentum_flux = 0.006759517700914085,
    infinity_carter_constant_flux = 0.0,
    horizon_energy_flux = -8.777498432552538e-6,
    horizon_angular_momentum_flux = -7.361786661019827e-5,
    horizon_carter_constant_flux = 0.0,
    total_modes = 45204,
    n_reached = (infinity = 124, horizon = 58),
    convolution_integral = (strategy = "ISEM adaptive trapezoidal for all computed n; tail ISEM adaptive Levin enabled but not triggered", tail_levin_nmin = 50, tail_levin_local_n = 16, tail_levin_max_depth = 8),
    tolerance = 1.0e-8,
    truncation_floor = (infinity = 1.0e-16, horizon = 1.0e-16),
    cost = 26.141580820083618 seconds,
)
```

</details>

<details>
<summary>Generic orbit</summary>

```julia
julia> flux = Teukolsky_pointparticle_flux(0.9, 6.0, 0.7, 0.5; tol=1e-8)
TeukolskyPointParticleFlux(
    orbital_parameters(a = 0.9, p = 6.0, e = 0.7, x = 0.5),
    orbit_type = generic,
    infinity_energy_flux = 0.0011469130759772415,
    infinity_angular_momentum_flux = 0.006380279789657716,
    infinity_carter_constant_flux = 0.04146329448763504,
    horizon_energy_flux = -9.200412363629793e-6,
    horizon_angular_momentum_flux = -0.0003627225703136506,
    horizon_carter_constant_flux = 0.0009653818407530993,
    total_modes = 1319722,
    n_reached = (infinity = 145, horizon = 84),
    convolution_integral = (strategy = "ISEM adaptive trapezoidal for all computed n; tail ISEM adaptive Levin enabled but not triggered", tail_levin_nmin = 50, tail_levin_local_n = 16, tail_levin_max_depth = 8),
    tolerance = 1.0e-8,
    truncation_floor = (infinity = 1.0e-16, horizon = 1.0e-16),
    cost = 1306.2777690887451 seconds,
)
```

</details>

## How to cite
If you have used this code in your research that leads to a publication, please cite the following article:
```
@article{Lo:2023fvv,
    author = "Lo, Rico K. L.",
    title = "{Recipes for computing radiation from a Kerr black hole using a generalized Sasaki-Nakamura formalism: Homogeneous solutions}",
    eprint = "2306.16469",
    archivePrefix = "arXiv",
    primaryClass = "gr-qc",
    doi = "10.1103/PhysRevD.110.124070",
    journal = "Phys. Rev. D",
    volume = "110",
    number = "12",
    pages = "124070",
    year = "2024"
}
```

Additionally, if you have used this code's capability to solve for solutions with complex frequencies, please also cite the following article:
```
@article{Lo:2025njp,
    author = "Lo, Rico K. L. and Sabani, Leart and Cardoso, Vitor",
    title = "{Quasinormal modes and excitation factors of Kerr black holes}",
    eprint = "2504.00084",
    archivePrefix = "arXiv",
    primaryClass = "gr-qc",
    doi = "10.1103/PhysRevD.111.124002",
    journal = "Phys. Rev. D",
    volume = "111",
    number = "12",
    pages = "124002",
    year = "2025"
}
```

If you have used this code's capability to solve for the gravitational waveform amplitudes and fluxes at infinity and at horizon with a test particle orbiting a Kerr black hole in a generic timelike bound and stable orbit (e.g., for extreme mass ratio inspiral waveforms), please cite the following articles:
```
@article{Yin:2025kls,
    author = "Yin, Yucheng and Lo, Rico K. L. and Chen, Xian",
    title = "{Gravitational radiation from Kerr black holes using the Sasaki-Nakamura formalism: waveforms and fluxes at infinity}",
    eprint = "2511.08673",
    archivePrefix = "arXiv",
    primaryClass = "gr-qc",
    doi = "10.1103/9ngz-k1lr",
    journal = "Phys. Rev. D",
    volume = "113",
    number = "12",
    pages = "124007",
    year = "2026"
}

@article{Lo:2025lpo,
    author = "Lo, Rico K. L. and Yin, Yucheng",
    title = "{Near-horizon gravitational perturbations of rotating black holes}",
    eprint = "2512.07937",
    archivePrefix = "arXiv",
    primaryClass = "gr-qc",
    doi = "10.1103/bljh-l413",
    journal = "Phys. Rev. D",
    volume = "113",
    number = "6",
    pages = "L061505",
    year = "2026"
}
```


## License
The package is licensed under the MIT License.
