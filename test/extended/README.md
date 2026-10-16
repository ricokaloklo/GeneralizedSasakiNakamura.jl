# Extended Tests

These tests are opt-in and are not included by `test/runtests.jl`.
They require the package's existing dependencies; the runner does not install anything.

## Running

From the repository root:

```sh
GSN_TEST_GROUPS=list julia test/extended/runtests.jl
GSN_TEST_GROUPS=quick julia --threads=1 test/extended/runtests.jl
GSN_TEST_GROUPS=full julia --threads=1 test/extended/runtests.jl
GSN_TEST_GROUPS=qnm,flux GSN_TEST_LEVEL=extended julia --threads=1 test/extended/runtests.jl
GSN_SOURCE="$PWD" GSN_TEST_GROUPS=quick julia --threads=1 test/extended/runtests.jl
GSN_SOURCE="$PWD" GSN_TEST_GROUPS=threads julia --threads=4 test/extended/runtests.jl
```

Without either selection variable, the runner prints usage and does not load the package.
`GSN_TEST_LEVEL` accepts `quick`, `full`, or `extended`. Setting only the level runs all groups.
`GSN_SOURCE` selects an include-style source checkout; otherwise the installed package is used.
`GSN_SWSH_PATH` optionally selects a local dependency checkout. No machine-specific path is assumed.

## Groups

| Group | Coverage |
|---|---|
| api | Public signatures, return types, display, invalid inputs |
| homogeneous | Independent amplitudes and radial values; identities and transformations |
| isem | Complex contours, endpoints, tolerance handling, explicit junctions |
| qnm | Frequencies, overtone conventions, residues, derivatives, excitation factors |
| pointparticle_mode | Circular, spherical and generic modes; interface consistency |
| flux | Mode sums, Schwarzschild rotation identities, independently summed Kerr grids |
| regression | Historical failures and a snapshot of earlier registered tests |
| threads | Bitwise serial/concurrent radial, QNM and particle-mode calls; shared-result evaluation |

The registered-suite snapshot is not a substitute for current `Pkg.test()`.
`@test_broken` records an unresolved numerical target; it does not count as a successfully solved case.
Reference agreement and equation identities are kept separate from frozen-output regressions.
Quick runs cover the numerical groups on smaller parameter sets.
The thread group skips itself when only one Julia thread is available.

## Reference Conventions

An overtone label is not always a continued-fraction inversion index. The independent Leaver
root checks explicitly select the reference root and inversion index. Separate published-table
checks exercise the default overtone convention.

Detailed Teukolsky excitation factors use
`B_ref / (2 i omega dB_inc/domega)`; the independent qnm(0.68) reference uses `2 omega`.
The test converts that reference by `1/i`, including the corresponding mirror sign.

The reference tables are included in `data/`; running the suite needs no external reference files.
The `src` fields of some tables name the original validation records; they are not runtime dependencies.
The reference-generation scripts are not distributed, except `tools/richartz_extremal_reference.py`.

## Provenance of every reference value

Package-independent ("independent") references are preferred; package outputs are used only where
marked "frozen self-value", with looser tolerances.

**Radial amplitudes (homogeneous, isem, regression)** -- `data/radial_references.jl`
* `STRONG_DAMPING_POINTS` G1-G10: independent mpmath MST and direct complex-path Taylor integration,
  90-160 digits, only agreeing digits kept. E (spheroidal eigenvalue) from the same computation.
* Q_*_f64 (24 points at the exact binary QNM frequencies and +0.001): independent, same solvers; lambda given.
* R_* (27 real points, s = -2, l = m = 2..4, a = 0/0.7/0.9, omega = 1/2/3) and D_* (12 weakly damped):
  independent validation values.
* Tolerance classes: 1e-9 in dominant units (|v - ref| / max(|inc|, |ref|) of the branch); the few wider
  classes are explained in the row's `notes`. The classes were set from measured package errors; endpoint
  changes are checked against the same independent values.
* `GRID_POINTS` (180 parameter sets): the parameter grid of the validation sweep.
* `GRID_REFERENCES` (102 pairs): independent direct Teukolsky Taylor references with independently derived
  endpoint series, 640/896 bits. 896-bit values are used; the 640-bit run disagrees on two deep UP
  incidences (ids 32, 178) and is not used there. `meas_*` columns are measured package errors that select
  tolerance classes only.
* `CLOUD_CELLS` (15 strong-damping cells): independent real-axis R(r) references (horizon Frobenius and exact
  Taylor steps along real r, dps 100/140; lambda by Leaver CF with a spectral check).
* `POSITIVE_AXIS_ROWS`, `POSITIVE_AXIS_POINTS`: independent 384-bit Teukolsky integration near and on the positive
  imaginary axis.
* Hard-coded in `isem` / registered suite: negative-frequency and default-tolerance references
  (`test/runtests.jl`, from the same independent tables), and the real-frequency MST value nu = 2.5 + 1.262i.
* The amplitude Wronskian `C_up (4 i k r_+ + 2 s kappa) = 2 i omega B_inc` (unit transmissions) follows from
  the package's asymptotic forms (`common/helpers.jl`); the suite also checks it on the reference values themselves.

**QNM** -- `data/qnm_references.jl`
* `QNM_CLOUD_ROOTS`: independent Leaver radial + angular CF roots, dps 50/75, 47-57 converged digits.
* `SCHWARZSCHILD_SCAN`: independent Leaver roots, dps 25, inversion index n; mirror roots stored for n >= 9.
* `PUBLISHED_QNM`: Forteza & Mourier arXiv:2107.11829 (10 digits) and Cook, Zenodo 14024959.
* `N20_CF_ROOTS`: six independent radial/angular Leaver-CF roots, at decimal spins
  0.5, 0.68 and 0.9. Increasing precision from 50 to 75 digits, radial depth from
  16000 to 32000 and changing inversion index from 20 to 17/23 preserves at least
  43 digits; spectral continuation checks the angular branch. The corresponding
  Cook rows differ by 8e-9--1.1e-8. The catalogue retains its 2e-8 coarse check;
  the solver must also agree with these independent roots within 1e-10.
* `PUBLISHED_RESIDUES`: arXiv:2609.09531v1 Table I (a = 0.5), converted to the package's residue
  convention. Ordinary rows n >= 9 keep the table's complete-spectrum ordering through explicit
  convention metadata.
* `LEGACY_OVERTONES`: a separate legacy implementation, cross-checked against the published rows at
  n = 8, 9, 20; atol 1e-9. The co-rotating mirror rows at a = -0.5 keep the legacy complete-spectrum labels explicitly.
* `AXISYMMETRIC_N10`: legacy axisymmetric catalogue, confirmed by an independent root.
  It is the re-emerged segment `10₁` (Cook & Zalutskiy), checked as `qnm_frequency(QNMMode(-2, 2, 0, 10), 0.5; segment=1)`
  with atol 1e-9. The NIA endpoints of segment 0 (n = 9, 10, 15) use Cook & Zalutskiy (PRD 94, 104074) Table II at atol 1e-5.
* `Q068`, `Q068_STENCIL`: independent direct solver, dps 40/60. Its excitation factor uses a denominator 2 omega;
  the detailed Teukolsky result uses 2 i omega, so the comparison divides the reference by i.
* Exact-extremal roots: Richartz, PRD 93, 064062 (2016), Tables I, V, VI, cross-checked with
  `tools/richartz_extremal_reference.py`.
* Registered values (overtone n29 at a = 0.5 and n8 branches) as in `test/runtests.jl`.

**Point-particle modes** -- `data/mode_references.jl` and constants in `pointparticle_mode/runtests.jl`
* Schwarzschild circular (2,+-2) at p = 10: independent Zerilli-Moncrief, MST and direct Teukolsky
  integration (40-50 agreeing digits). p = 1000 / 10000: direct Teukolsky integration only, dps 50.
* `MODES_56`: published GSN 0.9.0 values (<= 4e-13 from independent mpmath per-mode values) and frozen
  self-values (used only for the horizon at rtol 1e-9).
* `TABLE_MODES`: rows of the tol = 1e-10 eccentric, generic and spherical summation tables (frozen
  self-values of the summation path; cross-path check at rtol 1e-6/1e-5).
* Schwarzschild orbit constants: closed forms (no data). `ORBITS_DEV`: frozen KerrGeodesics 0.4 values.

**Fluxes** -- `data/flux_references.jl`
* `EQUATORIAL_TOL10`: circular/eccentric totals, normally at tol = 1e-10.
  The Schwarzschild prograde e = 0.2 and 0.5 rows use the lower-sideband stop
  rule at tol = 1e-11.
  They share the radial solver with the tested sums, so they are not independent radial references.
* a = 0 rotation identities (E = E_eq, Lz = x Lz_eq, Q = 2 L (1-x^2) Lz_eq,
  L = p / sqrt(p - 3 - e^2)) applied to `EQUATORIAL_TOL10`.
* `SPHERICAL_FULLGRID`: independent of any stop rule: every mode
  l <= 30, |m| <= 30, k <= 20 at mode tol 1e-10; outer energy layers <= 2e-11 of the totals.
* `GENERIC_FULLGRID`: exhaustive sums over l <= 18, |m| <= 16, |k| <= 18, n <= 16,
  at mode tol 1e-10; outer energy layers <= 7.4e-11 of the totals. These and the
  spherical grids validate stopping rules, not the individual radial/source solver.
* `CANDIDATE_TOL10`: frozen self-values of an alternative stop-rule implementation; used only where no
  full-grid row exists; agrees with the full grid to <= 3e-11 at the infinity totals where both exist.
* `A0_RULEC_TOL10`: informational only.


## Known Limitations

Open targets include poorly conditioned strong-damping
radial values and their error estimates, suppressed reflection coefficients, some explicit
junctions. The tests retain the numerical targets
rather than silently widening them. The frequency-label and excitation-phase comparisons
above are convention conversions, not failed solver targets.

Some high-overtone residue, near-static amplitude and negative-imaginary-axis jump checks lack
an independent reference and remain skipped. Counts of passing, broken and skipped tests must
be reported separately.
