# Issue inventory and labeling scheme: cppqed + cppqed-multidiagonal (draft, 2026-10-08)

Temporary IDs: `C-nn` = WignerQuantumOptics/cppqed, `M-nn` = WignerQuantumOptics/cppqed-multidiagonal (both under `vukics` until 2026-10-21), `E-n` = parent/epic issue.
Labels in [brackets] follow the scheme in §2. "Src" points to findings.md (A1…E3), the library ledger (F1, F1b, F2…F11),
CLAUDE.md, or your recollection (U).

---

## 1. Where the numbered findings live today

- **Library CLAUDE.md "Findings ledger".** This was the only list of #1–#11 (now renamed F1, F1b, F2…F11; see below).
- **Library test names.** They carry #1, #1b, #2…#7, and the numbering disagrees with CLAUDE.md:
  - The tests call "colliding pairs" #1b; CLAUDE.md calls it #2.
  - The tests' #2 is "JSON serialization of MultiDiagonal", which is missing from the CLAUDE.md list.
  - #8–#11 appear in no test name.
- **Library commit messages.** No commit references a finding number, although CLAUDE.md requires it.
- **cppqed CLAUDE.md.** It referenced "finding #5" (now F5).
- **findings.md** (this review): A1…E3.
- **Nothing on GitHub, as far as the repos show.** I cannot see GitHub issues; a dump is needed to dedupe.

**Clash:** on GitHub, `#5` in any issue, PR or commit auto-links to issue/PR 5 of that repo. cppqed already has
v2-era issues #4–#32 (commits reference #4, #5, #6, #9, #11, #12, #25–#32), plus SourceForge tickets.
So "finding #5" written in a cppqed issue links to an unrelated 2020 issue.
→ **Done (2026-10-15):** the ledger is renamed **F1, F1b, F2…F11**, with the numbers kept exactly as the tests had them (CLAUDE.md's list had omitted F2). Both CLAUDE.md files are updated. The F-series is closed: new findings get GitHub issue numbers, cited in test names and commits (`#N` in-repo, `WignerQuantumOptics/cppqed-multidiagonal#N` from elsewhere). The test names in the library sources are renamed too; what remains is in M-17.

---

## 2. Labeling scheme (identical families in both repos; only `area:` differs)

| Family | Values | Rule |
|---|---|---|
| ~~`type:`~~ | replaced by the organization's native issue types (2026-10-20), see the mapping below | — |
| `severity:` | `S1-critical` (silently wrong physics or data loss), `S2-high`, `S3-medium`, `S4-low` | bugs/physics only |
| `area:` (cppqed) | `utils/ode`, `utils/trajectory`, `utils/slicing`, `utils/pars`, `utils/random`, `core/quantumdata`, `core/structure`, `core/qjmc`, `core/master`, `core/propagator`, `elements/mode`, `elements/qbit`, `elements/spin`, `elements/composite`, `cmake`, `docs` | one or more |
| `area:` (library) | `multiarray`, `multidiagonal`, `complex`, `tests`, `benchmarks`, `cmake`, `ci`, `docs` | one or more |
| `state:` | `reproduced`, `needs-repro`, `latent` (code currently unreachable or uncompiled) | optional |
| `origin:` | `split` (cppqed↔library), `v2-migration`, `review-2026-10` | optional; useful for tracing |

Other GitHub dimensions:
- **Native issue types** (WignerQuantumOptics, set up 2026-10-20; exactly one per issue). The `type:` values
  used in the tables below map as follows when filing:

  | Issue type | Former `type:` label(s) | Count |
  |---|---|---|
  | Bug | `bug` | 26 |
  | Physics | `physics` | 10 |
  | Feature | `enhancement` | 10 |
  | Task | `build`, `test`, `cleanup`, `performance` (the `area:` labels keep `cmake`/`ci`/`tests` apart) | 23 |
  | Refactor | `refactor` | 9 |
  | Docs | `docs` | 5 |
  | Decision | `discussion` | 11 |

  Task, Bug, Feature and Decision were already defined for Wigner_Time; Physics, Refactor and Docs were added
  for us. The types are organization-wide, so don't rename or delete existing ones.
- **Priority and Effort** are the organization's **issue fields**, not labels or Project-only fields:
  - Priority: Urgent / High / Medium / Low
  - Effort: High / Medium / Low

  Priority changes over time; severity doesn't, so severity stays a label. The organization also has
  Start date and Target date issue fields; we don't plan to use them for now.
- **Milestones** (per repo):
  - cppqed:
    - `M0 split & build`
    - `M1 correct QJMC`
    - `M2 Master restored`
    - `M3 propagators / interaction picture`
    - `M4 ensembles & autostop`
    - `M5 MDoF`
    - `M6 docs`
  - library:
    - `v0.2.1` (fixes)
    - `v0.3.0` (layout change + apply performance)
- **Sub-issues** model the epics below. **Dependencies** ("blocked by") model ordering, e.g. C-01 is blocked
  by M-01 and C-10 is blocked by M-02. Both features exist for personal repos.
- **One organization Project** "C++QED v3" spanning both repos.
  - Created 2026-10-21, number 1; private since 2026-10-22, because the restricted-visibility issue fields
    cannot be added to a public Project. Fields: Status (default), Priority and Effort (the issue fields).
  - Priority and Effort: **decided 2026-10-21: the organization issue fields** (same convention as Wigner_Time),
    set on each issue in Phase 4 via `setIssueFieldValue`; shown in the Project since 2026-10-22.
  - Views: Board by Status for now; by milestone, critical (severity filter) and by repository can be added
    later in the browser.
- **Granularity:** one issue per independently fixable defect. Small same-file items are bundled as a
  checklist issue (marked "bundle" below).

---

## 3. WignerQuantumOptics/cppqed

Filing in progress (2026-10-22): the epics and the E-1 section are filed; their numbers are linked below. Labels use
`area:core/composite` instead of `elements/composite` (Composite.h lives in CPPQEDcore/structure/).

### Epics (parents)
- **E-1 [#14](https://github.com/WignerQuantumOptics/cppqed/issues/14) Switch to cppqed-multidiagonal.** Children: C-01…C-05. [type:enhancement area:cmake origin:split] → M0
- **E-2 [#15](https://github.com/WignerQuantumOptics/cppqed/issues/15) Master equation gives physically wrong / diverging results (U).** [type:physics severity:S1-critical area:core/master] → M2
  - Symptom as you remember it; record the conditions if you can.
  - Note that Master has not compiled since 2bc5d3f8 (2024-04-30), so the observation predates that.
  - Children: C-10, C-11, C-12, C-23, the library fix F5 (arrives via E-1), M-02, and C-13 as acceptance.
- **E-3 [#16](https://github.com/WignerQuantumOptics/cppqed/issues/16) Third-party dependencies via FetchContent (U).** Children: C-07, C-08, and M-11 (cross-repo). [type:build area:cmake] → M0
- **E-4 [#17](https://github.com/WignerQuantumOptics/cppqed/issues/17) Robust unattended cluster runs.** Children: C-40…C-47. [type:bug area:utils/trajectory] → M1
- **E-5 [#18](https://github.com/WignerQuantumOptics/cppqed/issues/18) Interaction-picture propagators (roadmap 1).** Children: C-60, C-14, C-61. [type:enhancement area:core/propagator] → M3

- **E-7 [#20](https://github.com/WignerQuantumOptics/cppqed/issues/20) Replace the own MultiArray(View) machinery with std::mdspan (cross-repo, decision first).** [type:discussion → enhancement, area:core/*, origin:split] → after M1; library side before the v1.0 / paper decision
  - Children: C-80…C-85 here, plus M-14 (the decision) and M-19 (library prototype).
  - Proposal:
    - the view becomes `std::mdspan` used directly through alias templates, with no wrapper class
    - a thin owning type stays (C++ has no standard owning mdarray)
    - the domain-specific retained-axes slicing stays, rebuilt on `submdspan`
  - Toolchain: GCC 13 (local) ships no `<mdspan>`. Use the Kokkos reference implementation (Apache-2.0 WITH LLVM-exception, header-only, FetchContent-pinnable) until the cluster compilers provide `std::mdspan`.

### Split (E-1)
| ID | Title | Labels | Src |
|---|---|---|---|
| C-01 [#21](https://github.com/WignerQuantumOptics/cppqed/issues/21) | Consume cppqed-multidiagonal via FetchContent (SHA-pinned); delete stale MultiArray*/MultiDiagonal copies and MultiDiagonal.cc | type:build, area:cmake, origin:split; blocked by M-01 | A1, A5 |
| C-02 [#22](https://github.com/WignerQuantumOptics/cppqed/issues/22) | Framework Traits.h: keep only the extras on top of `<multidiagonal/Traits.h>`; include it directly where used (consider renaming) | type:build, origin:split | A3 |
| C-03 [#23](https://github.com/WignerQuantumOptics/cppqed/issues/23) | Port Mode/Qbit/Spin operator construction to signed offsets (sign trap: lower → negative) | type:bug, area:elements/*, origin:split | A4 |
| C-04 [#24](https://github.com/WignerQuantumOptics/cppqed/issues/24) | Decide where SliceIterator.h lives (library or cppqed) | type:discussion, area:utils/slicing | A8 |
| C-05 [#25](https://github.com/WignerQuantumOptics/cppqed/issues/25) | Acceptance: jc_qjmc output identical before/after the switch (already shown in a scratch build; make it a test) | type:test, origin:split | A9 |

### Dependencies (E-3)
| ID | Title | Labels | Src |
|---|---|---|---|
| C-07 | Replace submodules with SHA-pinned FetchContent of the **upstream** repos (badaix/popl, imneme/pcg-cpp, Reputeless/Xoshiro-cpp). The forks are not needed: popl and Xoshiro have no own commits, and pcg's only own commit (d33c84d) was superseded upstream by ae2f711 (verified by a round-trip test). Archive the forks; don't delete them, since history .gitmodules reference them. Offline cluster build via `FETCHCONTENT_FULLY_DISCONNECTED` + `FETCHCONTENT_SOURCE_DIR_<NAME>`. Keep Xoshiro (decided 2026-10-19): it is a second independent-stream PRNG; RNG tests cover it alongside pcg64. | type:build, area:cmake | U, E2 |
| C-08 | Drop the unused Boost.Iostreams requirement (breaks minimal cluster Boost installs) | type:build, severity:S4-low | C16 |

Note on C-01/C-07: CMake's built-in `FETCHCONTENT_SOURCE_DIR_CPPQED-MULTIDIAGONAL=<path>` already provides the
local co-development override, so a custom `MULTIDIAGONAL_LOCAL` variable is redundant. Add
`FIND_PACKAGE_ARGS` so an installed library is used when present.

### Master (E-2)
| ID | Title | Labels | Src |
|---|---|---|---|
| C-10 | Port Master.h and Schrödinger.h to the (ha, li, ev[, ex]) API (they reference the deleted QuantumSystemDynamics.h) | type:bug, severity:S2-high, area:core/master, state:reproduced; blocked by M-02 | B4 |
| C-11 | Master: trailing `conj(m.rho)` conjugates ρ every step | type:physics, severity:S1-critical, area:core/master, state:latent | B3 |
| C-12 | Mode superoperators use √rate instead of rate | type:physics, severity:S1-critical, area:elements/mode | B5 |
| C-13 | Driven-qubit Bloch test: Master vs QJMC ensemble vs analytic steady state | type:test, area:core/master | CLAUDE.md |
| C-14 | Master with a non-unitary U needs interaction-picture-aware superoperators (or restrict to the unitary part, as v2 did). Background: v3dev project card "QuantumSystems should report their picture" (Sch / UIP / IP; quoting it: "Master with no specific superoperator", "Integrating QJMC method can be used only for UIP and Sch") | type:physics, area:core/master, core/propagator | B14 |

### QJMC and elements physics (M1)
| ID | Title | Labels | Src |
|---|---|---|---|
| C-20 | ODE abs/rel tolerances swapped in make_controlled | type:bug, severity:S2-high, area:utils/ode | B1 |
| C-21 | H_eff misses ½ΣJ†J identity term → thermal mode from vacuum never jumps | type:physics, severity:S1-critical, area:elements/mode, elements/qbit, state:reproduced | B2 |
| C-22 | Qbit gain channel built with γ_m | type:physics, severity:S2-high, area:elements/qbit | B6 |
| C-23 | SliceIterator: first slice drops parent offset; reads past end of offsets | type:bug, severity:S1-critical, area:utils/slicing, state:reproduced | B7 |
| C-24 | Composite rate lambda rank breaks every 3-subsystem QJMC at compile time | type:bug, severity:S2-high, area:elements/composite | B8 |
| C-25 | Broadcaster offsets assume contiguous strides (cross-rank nesting) | type:bug, severity:S2-high, area:elements/composite | B9 |
| C-26 | Broadcaster captures leaf descriptors by reference (dangling) | type:bug, severity:S3-medium, area:elements/composite | B10 |
| C-27 | QJMC does not normalize initial ψ | type:physics, severity:S3-medium, area:core/qjmc | B11 |
| C-28 | Conventions: qbit::sz sign; Qbit polarization ⟨σ†⟩ vs Mode ⟨a⟩ (bundle) | type:physics, severity:S4-low, area:elements/qbit | B15 |
| C-29 | Mode gain rate counts the truncation edge | type:physics, severity:S4-low, area:elements/mode | B16 |
| C-30 | QJMC bisection discards LogTree → no in-trajectory log for dc>1 | type:bug, severity:S4-low, area:core/qjmc | C12 + CLAUDE.md |
| C-31 | `qbit::make(Pars)` passes δ as the complex frequency: no detuning, no γ_m in H_eff, excited amplitude grows | type:physics, severity:S1-critical, area:elements/qbit | B17 |

### Cluster-run robustness (E-4)
| ID | Title | Labels | Src |
|---|---|---|---|
| C-40 | Continuation silently restarts from t=0 with empty/stale .state | type:bug, severity:S1-critical, area:utils/trajectory, state:reproduced | C1 |
| C-41 | Unknown command-line options must abort | type:bug, severity:S2-high, area:utils/pars, state:reproduced | C2 |
| C-42 | Default --dc overrides --Dt/--NDt | type:bug, severity:S3-medium, area:utils/trajectory, state:reproduced | C3 |
| C-43 | --streamSwitch bit order reversed vs help text | type:bug, severity:S3-medium, area:utils/trajectory | C4 |
| C-44 | NDt mode not idempotent under continuation | type:bug, severity:S3-medium, area:utils/trajectory | C5 |
| C-45 | State-file framing: unchecked gcount, checkpoints never flushed | type:bug, severity:S3-medium, area:utils/trajectory, state:needs-repro | C6 |
| C-46 | ODE edge cases: backward first step → infinite loop; NaN accepted; FSAL guard (bundle) | type:bug, severity:S4-low, area:utils/ode | C11 |
| C-47 | Pars objects must be non-copyable/non-movable | type:bug, severity:S4-low, area:utils/pars, state:needs-repro | C13 |

### Build quality, tests, cleanup
| ID | Title | Labels | Src |
|---|---|---|---|
| C-50 | ODR violations in headers (defaultLogger, observerNoOp, Mode.h/Qbit.h functions); blocks multi-TU tests | type:bug, severity:S3-medium, area:cmake | C7 |
| C-51 | Framework test suite (Catch2): Mode, Qbit, JC, Master invariants, QJMC→Master convergence. Source of analytic cases: the physics tests of the C++QEDv2m9-beta package (MultiLevel tests spelled out in `testsuitePhysics/test.sh`, per a v3dev project card) | type:test; blocked by C-50 | CLAUDE.md |
| C-52 | ~~jc_qjmc.cc~~ (un-ignored and committed 2026-10-20) and benchmarking/ are gitignored | type:build, severity:S2-high | E1 |
| C-53 | Remove or quarantine dead v2 code (Sigma.h, EntanglementMeasures.h, DistributionFunctions, LazyDensityOperatorFFT, Simulated.h, Profiling.h, Mode_.h, CoupledModes, unused frees/interactions). **Before removing `EntanglementMeasures.h`, decide:** `QuantumTrajectory.h` still includes it, and the v2-era card "Resurrect negativity calculation using Eigen" (plus mutual information) is marked Done; keep and wire up, or delete. Also small modernizations: `FormDouble` (`QuantumTrajectory.h`) → `std::format` (a commented draft is in `Trajectory.h`); should `precision` be unsigned? | type:cleanup, origin:v2-migration | C14 |
| C-54 | Structure-layer defects (bundle): SparseMatrix offsets/strides + missing throw; Lindblad without rate; NoOp Hamiltonian; LDO ctor ambiguity; EV concept check; 3-arg H sees τ; static-const ODR; TemporalDataPoint label; superoperatorFromJump | type:bug, severity:S3-medium, area:core/structure | C14, C15, A7 |
| C-55 | Performance bundle: variant-by-value offset cache; checkpoint copies stepper; .state grows unbounded | type:performance, area:utils/trajectory, core/structure | C16 |
| C-56 | ~~CLAUDE.md and README corrections~~ (done 2026-10-20: E2 in CLAUDE.md, README in f1f20c85; library CLAUDE.md layout, D11, also fixed) and the version-banner correction (E3, still open) | type:docs | E2, E3 |

### Roadmap and known problems (from CLAUDE.md)
| ID | Title | Labels | Src |
|---|---|---|---|
| C-60 | UnaryDiagonalPropagator bodies; propagator composition/broadcast; QJMC propagator hook + interaction-picture reset; wire up Qbit and Mode | type:enhancement, area:core/propagator → M3 | CLAUDE.md, B12, C15 |
| C-61 | Reinstate QJMC stepwise (needs full H_eff; anti-Hermitian part exactly once; fix `::time`, calculateRates) | type:bug, area:core/qjmc, state:latent → M3 | CLAUDE.md, B13 |
| C-62 | Reinstate makeEnsemble (dangling average, unnormalized ψ, compile failure, nTraj in stateIO) + RNG stream mapping | type:bug, area:core/qjmc, utils/random, state:latent → M4 | CLAUDE.md, C8, C9 |
| C-63 | Restore autostop | type:bug, area:utils/trajectory, state:latent → M4 | CLAUDE.md, C10 |
| C-64 | structure::concatenate fails for non-constexpr args (derive rank from Lindblad::N_RANK) | type:bug, severity:S4-low, area:core/structure | CLAUDE.md, C15 |
| C-65 | Stubs: compose(), broadcastFor, Placeholder beyond RANK 2, InitializeEnsembleFromArrayOnlyArchive | type:enhancement | CLAUDE.md |
| C-66 | Motional degree of freedom (Particle/ParticleCavity port) | type:enhancement → M5 | CLAUDE.md |
| C-67 | Developer guide and tutorials (the operator-list fix in the guide is done, f1f20c85). Doxygen: MathJax physics extension, no class-hierarchy graphs, CMake switch for internal docs (Milestone 13 project card) | type:docs → M6 | CLAUDE.md |
| C-68 | Discussion: compile-time leading dimensions, OBSERVER concept unification, HDF5 serialization (one issue each) | type:discussion, no milestone | CLAUDE.md |

### Carried over from the v2-era GitHub projects (2026-10-20)
Source: projects #3 "C++QED v3dev" and #4 "C++QED Milestone 13" (dumped 2026-10-20, 115 cards; the rest was done, v2-specific, rejected by the v3 design, or already covered). References are copied from the cards, not verified.

| ID | Title | Labels | Src |
|---|---|---|---|
| C-90 | Benchmarks against QuTiP and QuantumOptics.jl, with and without output; prerequisite for Paper 2 | type:test, origin:v2-migration | #4 (in progress) |
| C-91 | PRNG study: taus2 / mt19937 / pcg64 / xoshiro256, ensemble and time averages vs the exact Master solution (card's case: driven lossy mode, cutoff 100, n_th 5, κ 1, Δ 5, η 3, T 20); home of the Xoshiro decision and the planned CPC paper on PRNGs in MCWF. Blocked by C-62 and E-2 | type:discussion, area:utils/random, origin:v2-migration | #4 |
| C-92 | Time averaging along a single trajectory, as an observer (ties to the OBSERVER discussion in C-68); unmerged 2017 local branch `timeAveraging` as reference | type:enhancement, area:utils/trajectory, origin:v2-migration | #4 |
| C-93 | Steady-state solvers in the QuTiP style (card cites P. D. Nation, arXiv:1504.06768); also a reference for C-13 independent of Master time evolution | type:discussion, origin:v2-migration | #3 |
| C-94 | Spectrum via the "faster-than-the-clock" algorithm (card cites Quantum, q-2022-02-10-649) | type:discussion, origin:v2-migration | #3 |
| C-95 | Master integrators that preserve positivity and trace by construction (card cites Cao & Lu, arXiv:2103.01194; also check QTM's methods), instead of fixing Hermiticity and trace by hand after each step | type:discussion, area:core/master, origin:v2-migration | #3, #4 |
| C-96 | numpy bridge for state files: numpy arrays as initial conditions (`ARRAY_ONLY`), `--sdf` state series into numpy; includes append-vs-replace for state saving and optional cleartext archives for debugging. Child of the serialization issue of C-68 | type:enhancement, area:utils/trajectory, origin:v2-migration | #3, #4 |

### mdspan migration (epic E-7)
| ID | Title | Labels | Src |
|---|---|---|---|
| C-80 | Decision record: mdspan vs own MultiArrayView; consequences for the library paper (joint with M-14) | type:discussion | library CLAUDE.md item 3 |
| C-81 | Add the mdspan dependency: Kokkos reference implementation via FetchContent, behind a `cppqed::mdspan` namespace alias, switching to `std::` when `__cpp_lib_mdspan` is defined | type:build, area:cmake | E-7 |
| C-82 | View aliases: `StateVectorView<R>` etc. become mdspan aliases (layout_left for owners, layout_stride for slices); replace `operator()` indexing with `[]`; Debug bounds checks via a checking accessor policy | type:refactor, area:core/quantumdata | E-7 |
| C-83 | Thin owning array (vector + mapping, move-only, `.view()`, Boost.Serialization + JSON) replacing MultiArray | type:refactor, area:core/quantumdata | E-7 |
| C-84 | Retained-axes slicing (SliceIterator, sliceRange, partialTrace, Broadcaster) rebuilt on submdspan; this structurally removes the C-23/C-25 class of bugs | type:refactor, area:utils/slicing, elements/composite | E-7, B7, B9 |
| C-85 | Follow-up: partly static extents (e.g. `extents<size_t,2,dynamic_extent>` for a qubit⊗mode system) realize "compile-time leading subsystem dimensions" (C-68) | type:enhancement | E-7, C-68 |

### Renames (epic E-6)
- **E-6 [#19](https://github.com/WignerQuantumOptics/cppqed/issues/19) Planned renames.** Children: C-70…C-76, plus M-01 and M-17 in the library. [type:refactor]
- **Timing:** do the structural renames (C-70, C-71, C-73, C-74) in one batch together with the split (M0). The split already touches includes in almost every file, so one round of churn costs less than two.
- **Issues to file before the renames** should cite file:line paths with a commit SHA, so the references stay valid afterwards.
- **Acronym casing (C-72)** is cosmetic and touches the most lines. The library CLAUDE.md says "don't mass-rename". Keep it as a decision issue and do it opportunistically, or in one sweep at a release boundary.

| ID | Title | Labels | Scope today | Src |
|---|---|---|---|---|
| C-70 | Rename `structure/` → `quantumdynamics/` (directory, namespace `structure`, CMake target `CPPQEDcore_structure`, docs) | type:refactor, area:core/structure → M0 | 130 lines in 46 files | CLAUDE.md |
| C-71 | Rename Free → SubSystem: `CPPQEDelements/frees/` → `subsystems/`, target `CPPQEDelements_frees`, docs and comments. The remaining `Free`/`FreeExact` references are in dead v2 code (goes with C-53) | type:refactor, area:elements/* → M0 | directory + 2 CMake files; ~22 textual hits, all in v2 leftovers | CLAUDE.md |
| C-72 | Acronym casing (capitalize first letter only: `Ode`, `Tdp`, `Ldo`, `Qjmc`, `Qsd`, `Ev`, `Json`): decide scope and timing | type:discussion → no milestone | ≈265 identifier occurrences in cppqed (ODE 94, EV 53, TDP 35, JSON 34, LDO 22, QSD 21, QJMC 6); in the library, `toStringJSON` sits in the shared Traits.h, so this is cross-repo | CLAUDE.md (both) |
| C-73 | Rename the framework `Traits.h` (e.g. `FrameworkTraits.h`) so it can never be confused with `<multidiagonal/Traits.h>` | type:refactor, origin:split → M0 | every header that includes "Traits.h" | A3 / C-02 |
| C-74 | Align the reference executable's name and location: `dockerEntry/jc_qjmc.cc` vs `JaynesCummingsQJMC.cc` in CLAUDE.md; consider `dockerEntry/` → `examples/` (do together with C-52, un-ignoring it) | type:refactor, area:cmake → M0 | 1 file, 1 CMake glob | E1, E2 |
| C-75 | Ledger renamed to F1, F1b, F2…F11 in all Markdown (done 2026-10-15); test names renamed as well; the remaining work is in M-17. Listed here only for the record | type:cleanup | cppqed CLAUDE.md | §1 |
| C-76 | Rename `qbit` → `qubit` (namespace, files, option suffixes in examples, docs); cosmetic | type:refactor, area:elements/qbit, origin:v2-migration | namespace + 2 files + docs | v3dev project card |

---

## 4. WignerQuantumOptics/cppqed-multidiagonal

Filed 2026-10-22 as #1–#20 (M-04 is the pre-existing #1, updated in place). Priority, Effort and the milestones of the
unassigned items were proposed during filing; GitHub is now the source of truth for them.

| ID | Title | Labels | Src |
|---|---|---|---|
| M-01 [#2](https://github.com/WignerQuantumOptics/cppqed-multidiagonal/issues/2) | Move headers to include/multidiagonal/ (and install there); blocks cppqed C-01 | type:build, area:cmake, origin:split → v0.3.0 | A2, D3 |
| M-02 [#3](https://github.com/WignerQuantumOptics/cppqed-multidiagonal/issues/3) | Const vectorize/matricize build mutable Eigen::Map from const pointer (regression, ad2836a) | type:bug, severity:S2-high, area:complex, state:reproduced → v0.2.1 | A6, D1 |
| M-03 [#4](https://github.com/WignerQuantumOptics/cppqed-multidiagonal/issues/4) | Tests can't detect a missing conj in MultiDiagonal::hermitianConjugate; plus the test-gap checklist | type:test, area:tests, state:reproduced | D2, D10 |
| M-04 [#1](https://github.com/WignerQuantumOptics/cppqed-multidiagonal/issues/1) | Multi-axis application 2–2.7× slower than Eigen sparse: overhaul apply (U) | type:performance, area:multidiagonal → v0.3.0 | U, CLAUDE.md item 1, D12 |
| M-05 [#5](https://github.com/WignerQuantumOptics/cppqed-multidiagonal/issues/5) | Benchmark honesty: picture1 variants bypass operator(); "eigen_sparse_pernnz" label; inverted plot.py curve; constructJC unvalidated | type:test, area:benchmarks | D7 (agent 10, 16, 17) |
| M-06 [#6](https://github.com/WignerQuantumOptics/cppqed-multidiagonal/issues/6) | README/mainpage: offset-sign formula backwards; (a†)²a² claim; C++20/GCC12; InteractionPictureDiagonal; Doxyfile version; performance claims | type:docs, severity:S2-high (the sign error misleads users) | D7 |
| M-07 [#7](https://github.com/WignerQuantumOptics/cppqed-multidiagonal/issues/7) | Header docs: "applies H/i", "derived type in CPPQEDcore", "homomorphism", stale F5 docstring | type:docs, area:multidiagonal, complex | D8 |
| M-08 [#8](https://github.com/WignerQuantumOptics/cppqed-multidiagonal/issues/8) | twoTimesRealPartOfSelf / hermitianConjugateSelf ignore view offset/strides: assert or document | type:bug, severity:S3-medium, area:complex | D9 |
| M-09 [#9](https://github.com/WignerQuantumOptics/cppqed-multidiagonal/issues/9) | CMake: target_compile_features cxx_std_23; pin Boost/Eigen versions; SameMajorVersion on 0.x; ARCH_INDEPENDENT; no install as subproject; option to skip benchmarks | type:build, area:cmake | D5, D6, D14 |
| M-10 [#10](https://github.com/WignerQuantumOptics/cppqed-multidiagonal/issues/10) | CI: -Werror; fix Clang warnings (MultiArray.h:135, :254); sanitizer job; find_package/FetchContent consumer job | type:build, area:ci, state:reproduced | D4 |
| M-11 [#11](https://github.com/WignerQuantumOptics/cppqed-multidiagonal/issues/11) | Catch2: SHA pin, FIND_PACKAGE_ARGS for offline builds (cross-ref cppqed E-3) | type:build, area:cmake | D6 |
| M-12 [#12](https://github.com/WignerQuantumOptics/cppqed-multidiagonal/issues/12) | Small header items (bundle): greedy variadic ctor; zero-extent diagonal; sqrAbs via std::norm; noInit zero-fills; twoTimesImagPartOf docs; internal-linkage ODR | type:bug, severity:S4-low | D13 |
| M-13 [#13](https://github.com/WignerQuantumOptics/cppqed-multidiagonal/issues/13) | Reconcile the layout docstring (row-major vs first-index-fastest) | type:docs | CLAUDE.md item 2 |
| M-14 [#14](https://github.com/WignerQuantumOptics/cppqed-multidiagonal/issues/14) | Decide std::mdspan migration vs current architecture (joint decision with cppqed E-7/C-80; affects the paper's API) | type:discussion | CLAUDE.md item 3 |
| M-19 [#19](https://github.com/WignerQuantumOptics/cppqed-multidiagonal/issues/19) | Prototype: MultiDiagonal on mdspan (Diagonal storage, operator(), composition), validated by the dense-oracle tests and compared in the benchmarks (also an attempt at M-04) | type:enhancement, area:multidiagonal; blocked by M-14 | E-7 |
| M-15 [#15](https://github.com/WignerQuantumOptics/cppqed-multidiagonal/issues/15) | README usage snippet and blessed construction API (fromStorage/singleDiag vs factory namespace) | type:discussion, docs | CLAUDE.md item 4 |
| M-16 [#16](https://github.com/WignerQuantumOptics/cppqed-multidiagonal/issues/16) | Install/export review; Boost.Serialization as optional dependency | type:build | CLAUDE.md item 5 |
| M-17 [#17](https://github.com/WignerQuantumOptics/cppqed-multidiagonal/issues/17) | Give F8–F11 a named test or `[contract]` tag (F1b, F2 and F7 also lack `[contract]`). The F-renaming of the Markdown and the test names is done (2026-10-15) | type:cleanup, area:tests | §1 |
| M-18 [#18](https://github.com/WignerQuantumOptics/cppqed-multidiagonal/issues/18) | MultiDiagonal JSON is write-only (no value_to) | type:enhancement, area:multidiagonal | D10 |
| M-20 [#20](https://github.com/WignerQuantumOptics/cppqed-multidiagonal/issues/20) | Doxygen site uses the stock theme: no `HTML_EXTRA_STYLESHEET`/`HTML_HEADER` in the Doxyfile. Adopt Doxygen Awesome (the theme the cppqed developer guide §7.8 recommends), vendored at a pinned release or fetched at a pinned tag in `docs.yml`; check its license | type:docs, area:docs | András, 2026-10-21 (verified in Doxyfile) |

Library renames (children of cppqed E-6, cross-repo): M-01 (headers → `include/multidiagonal/`) and M-17 (finding numbers in test names → F-prefix). If acronym casing (C-72) is adopted, `toStringJSON` in the library's Traits.h changes too.

Counts: cppqed 7 epics + 66 leaf issues; library 20. Several are bundles with checklists.

---

## 5. Before filing
1. **Dump existing issues in both repos** (you run it) and dedupe against them. Old v2-era cppqed issues may
   still be open. Decided 2026-10-19: close them and all v2-era projects as obsolete. Issues: closed (2026-10-20). Projects: skimmed 2026-10-20 (C-90…C-96, C-76 and notes in C-14, C-51, C-53, C-67); #3 and #4 to be closed, the empty #5 and #6 deleted.
2. ~~Settle the ledger renaming~~: done. Issue bodies must use F-numbers for old findings.
3. **Order of creation:** labels → milestones → epics → leaf issues (so sub-issue and blocked-by links can
   point at real numbers) → project. Generate the batch script after these decisions.
