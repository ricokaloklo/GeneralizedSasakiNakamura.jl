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
| isem_direct | Complex contours, endpoints, tolerance handling, explicit junctions |
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
`B_ref / (2 i omega dB_inc/domega)`; the imported qnm(0.68) reference uses `2 omega`.
The test converts that reference by `1/i`, including the corresponding mirror sign.

The tables below are included in `data/`. Source filenames identify the original
validation records, not runtime dependencies. Running the suite needs no external
reference files. The workspace retains the reference-generation scripts separately.

## Provenance of every reference value

Package-independent ("independent") references are preferred; package outputs are used only where
marked "frozen self-value", with looser tolerances.

**Radial amplitudes (homogeneous, isem_direct, regression)** -- `data/radial_references.jl`
* `STRONG_DAMPING_POINTS` G1-G10: `strong_damping_reference/reference_table.tsv`
  (independent: mpmath MST + direct complex-path Taylor integration, 90-160 digits, only agreeing digits;
  README there). E (spheroidal eigenvalue) from the same file.
* Q_*_f64 (24 points at the exact binary QNM frequencies and +0.001): `.../float64_inputs/f64_points.tsv`
  (independent, same solvers; lambda given).
* R_* (27 real points, s = -2, l = m = 2..4, a = 0/0.7/0.9, omega = 1/2/3) and D_* (12 weakly damped):
  `refs/validation.tsv` (independent cloud validation).
* Tolerance classes: 1e-9 in dominant units (|v - ref| / max(|inc|, |ref|) of the branch); exceptions cite
  backlog items 1, 2, 3 (and the n9/n12-off members of the same class). Historical measurements used to
  set the classes: `out/reference73_rd2base.tsv` (2026-10-04 19:39). Later endpoint repairs are checked against the same independent values.
* `GRID_POINTS` (180 parameter sets): copied from `out/grid_root.tsv` (Independent grid driver).
* `GRID_REFERENCES` (102 pairs): Independent direct Teukolsky Taylor references with independently derived
  endpoint series, 640/896 bits (`jq12_grid_reference_{teuk,jq12s1teuk,ycteuk}{640,896}*.tsv`;
  `review.md` "Independent broad-grid reference"). 896-bit values are used; the 640-bit
  run disagrees on two deep UP incidences (ids 32, 178) and is not used there.
  Measured shared-root errors: amplitudes from `out/grid_root.tsv` (2026-10-04 03:49, before the reflection
  batch, so Re omega < 0 cells 28 and 144 show O(1) there), R(r) from `out/rerr_base.tsv` (2026-10-04 20:39).
* `CLOUD_CELLS` (15 strong-damping cells): `cloud_refs/realaxis_R_reference/cells_reference.tsv`
  (independent, dps 100/140; id82 lambda corrected per the README there).
* Item 25: `ITEM25_ROWS` from `backlog25_reference_b384.tsv`, `ITEM25_AXIS` from
  `.../backlog25_wronskian_axis384.tsv` (independent 384-bit Teukolsky integration).
* Hard-coded in `isem_direct` / registered suite: negative-frequency and default-tolerance references
  (`test/runtests.jl`, from the same independent tables), item 18 nu = 2.5 + 1.262i (`review.md` "Round 2").
* The amplitude Wronskian `C_up (4 i k r_+ + 2 s kappa) = 2 i omega B_inc` (unit transmissions) is derived in
  `strong_damping_reference/README.md`; the suite also checks it on the reference values themselves.

**QNM** -- `data/qnm_references.jl`
* `QNM_CLOUD_ROOTS`: `strong_damping_reference/results/qnm_list.json` (independent Leaver radial + angular CF,
  dps 50/75, 47-57 converged digits).
* `SCHWARZSCHILD_SCAN`: `.../results/qnm_schw_scan.json` (independent, dps 25, inversion index n; mirror roots
  stored for n >= 9).
* `PUBLISHED_QNM`: Forteza & Mourier arXiv:2107.11829 (10 digits) and Cook, Zenodo 14024959, as copied in
  the processed published QNM tables (2026-09-24 extraction, not distributed).
* `N20_CF_ROOTS`: six independent radial/angular Leaver-CF roots, at decimal spins
  0.5, 0.68 and 0.9. Increasing precision from 50 to 75 digits, radial depth from
  16000 to 32000 and changing inversion index from 20 to 17/23 preserves at least
  43 digits; spectral continuation checks the angular branch. The corresponding
  Cook rows differ by 8e-9--1.1e-8. The catalogue retains its 2e-8 coarse check;
  the solver must also agree with these independent roots within 1e-10.
* `PUBLISHED_RESIDUES`: arXiv:2609.09531v1 Table I (a = 0.5), via
  the processed published QNM residue table (2026-09-24, not distributed) and the conversion by
  the residue conversion script (2026-09-24, not distributed). Ordinary rows n >= 9
  keep the table's complete-spectrum ordering through explicit convention metadata.
* `LEGACY_OVERTONES`: September-24 legacy bundle (separate implementation, cross-checked against the published
  rows at n = 8, 9, 20): the high-precision overtone reference table (2026-09-24, not distributed); atol 1e-9.
  The co-rotating mirror rows at a = -0.5 keep the legacy complete-spectrum labels explicitly.
* `AXISYMMETRIC_N10`: `.../qnm_axisymmetric_reference_20260924.tsv` (legacy catalogue; the independent comparison confirmed the target).
  It is the re-emerged segment `10₁` (Cook & Zalutskiy), checked as `qnm_frequency(QNMMode(-2, 2, 0, 10), 0.5; segment=1)`
  with atol 1e-9. The NIA endpoints of segment 0 (n = 9, 10, 15) use Cook & Zalutskiy (PRD 94, 104074) Table II at atol 1e-5.
* `Q068`, `Q068_STENCIL`: `cloud_refs/qnm068_bridge_reference/{q068_root.json,stencil_reference.tsv}`
  (independent direct solver, dps 40/60). Its excitation factor uses a denominator 2 omega;
  the detailed Teukolsky result uses 2 i omega, so the comparison divides the reference by i.
* Registered values (overtone n29 at a = 0.5, called n30 before relabeling, and n8 branches) as in `test/runtests.jl`.

**Point-particle modes** -- `data/mode_references.jl` and constants in `pointparticle_mode/runtests.jl`
* Schwarzschild circular (2,+-2) at p = 10: `flux_reference/README.md` (independent:
  Zerilli-Moncrief, MST, direct Teukolsky; 40-50 agreeing digits). p = 1000 / 10000:
  `flux_reference/direct/log_r1000_d50.txt`, `log_r10000_d50.txt` (direct method only, dps 50).
* `MODES_56`: `flux_compare/modes_published.tsv` (published GSN 0.9.0; <= 4e-13 from the cloud's mpmath
  per-mode values, HANDOVER_20261005.md) and `flux_compare/modes_dev.tsv` (frozen self-values, used only for
  the horizon at rtol 1e-9).
* `TABLE_MODES`: rows of `flux_compare/modes_{ecc,gen,sph}_*_tol1e-10.tsv` (frozen self-values of the
  summation path; cross-path check at rtol 1e-6/1e-5).
* Schwarzschild orbit constants: closed forms (no data). `ORBITS_DEV`: frozen KerrGeodesics 0.4 values from
  `flux_compare/orbits_dev.txt`.

**Fluxes** -- `data/flux_references.jl`
* `EQUATORIAL_TOL10`: circular/eccentric totals, normally at tol = 1e-10.
  The Schwarzschild prograde e = 0.2 and 0.5 rows use the corrected sideband
  rule at tol = 1e-11 (`flux_eqref_tol1e-11.tsv`, `flux_eqref05_tol1e-11.tsv`).
  They share the radial solver with the tested sums, so they are not independent radial references.
* a = 0 rotation identities: `flux_compare/check_a0.py` (E = E_eq, Lz = x Lz_eq, Q = 2 L (1-x^2) Lz_eq,
  L = p / sqrt(p - 3 - e^2)) applied to `EQUATORIAL_TOL10`.
* `SPHERICAL_FULLGRID`: `flux_compare/fullgrid_sph_full*.tsv` (independent of any stop rule: every mode
  l <= 30, |m| <= 30, k <= 20 at mode tol 1e-10; outer energy layers <= 2e-11 of the totals).
* `GENERIC_FULLGRID`: exhaustive sums over l <= 18, |m| <= 16, |k| <= 18, n <= 16,
  at mode tol 1e-10; outer energy layers <= 7.4e-11 of the totals. These and the
  spherical grids validate stopping rules, not the individual radial/source solver.
* `CANDIDATE_TOL10`: `flux_compare/fluxP_ref/*_tol1e-10.tsv` (frozen self-values of the private stop-rule
  candidate trees/fluxP; used only where no full-grid row exists; agrees with the full grid to <= 3e-11 at
  the infinity totals where both exist).
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
