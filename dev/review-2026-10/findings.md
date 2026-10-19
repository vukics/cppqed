# C++QED v3 + cppqed-multidiagonal: code review findings (2026-10-08)

Assumption throughout: cppqed consumes cppqed-multidiagonal (library HEAD 4a6d6ac, v0.2.0) instead of its own copies.
cppqed at cd215fbb (v3/dev).

Evidence tags:
- **[R]** reproduced: I built or ran something that shows it.
- **[V]** verified: I read the code and confirmed it.
- **[A]** reported by a review agent after reading the code, not independently checked by me. Where two agents found it independently it says so.

---

## A. The split itself (cppqed ↔ library)

**A1. cppqed does not use the library at all [R].** It has no FetchContent and no `MULTIDIAGONAL_LOCAL`; it compiles its own stale copies:
- `CPPQEDutils/{Traits,MultiArray,MultiArrayComplex}.h`
- `CPPQEDcore/quantumoperator/MultiDiagonal.{h,cc}`

Drift from the library:
- `MultiDiagonal` is still the old upper/lower design; the library uses signed offsets. The diff is about 240 lines.
- The copies still contain library findings that were fixed in the library: F3 (`load()` dangling `dataView`), F4 (`hermitianConjugateSelf` doesn't conjugate the diagonal), F5 (`twoTimesRealPartOfSelf` quadruples the diagonal), F6 (`relativeDeviation(0,0)` = NaN), F7 (JSON complex read), F8 (serialization include), F10 (narrowing), F11 (sign-compare).
- Practical consequence today: the stale `load()` dangles if a state file with different extents is loaded, e.g. via `--initialFileName`.

**A2. The header names collide, and the `<multidiagonal/X.h>` layout both CLAUDE.md files describe does not exist [R].**
- The library exposes and installs a flat `include/{Traits,MultiArray,MultiArrayComplex,MultiDiagonal}.h`.
- Both `Traits.h` files define `toStringJSON`, `passByValue_v`, `ReferenceMF`, `Reference`, `ConstReferenceMF` and `ConstReference`, so including both is a redefinition error.
- Fix in the library: move the headers to `include/multidiagonal/` and install them there. The simulated build below needed exactly this.

**A3. cppqed's `Traits.h` must keep the framework extras, and every user must include it directly [R].**
- The extras are `Empty`, `intro_logger`/`outro_logger`/`logger`, `hana_sequence`, `always_false(_v)`, `overload`, `multilambda`, `compileTimeOrdinals`, `labelled`/`getLabel`, and Boost.Serialization for `LogTree`.
- `Hamiltonian.h`, `Liouvillian.h` and `SliceIterator.h` get `hana_sequence`, `always_false`, `overload` and `multilambda` only indirectly, through `"MultiArray.h"` → `"Traits.h"`. Once that path leads to the library, these break; the simulated build failed exactly here.
- Consider renaming cppqed's file (e.g. `FrameworkTraits.h`) so the two `Traits.h` can never be confused.
- `SliceIterator.h` gets `boost::forward_iterator_helper` only indirectly, via `serialization/vector.hpp` [A].

**A4. Five operator constructions use the old API [R]:**
- `Mode.cc:13` (`aOp`)
- `Qbit.cc:42` (`splus`) and `Qbit.cc:54` (`sz`)
- `Spin.cc:15` (`splus`) and `Spin.cc:35` (`sz`)

Each uses `res.diagonals[k].emplace(Offsets{o}, Diagonal{std::vector})`. The new form is `res.diagonals.emplace(Offsets{±o}, multiarray::fromStorage(v))`.
- **Sign trap:** old index `[1]` means upper (+o) and `[0]` means lower (−o). A naive port of `qbit::splus` and `spin::splus` that keeps `{1}` silently builds σ instead of σ†.

**A5. `quantumoperator/MultiDiagonal.cc` must leave the build [R].** It redefines the library's `inline multidiagonal::identity` and defines a rank-1 `operator|` over the old map.

**A6. The library's const `vectorize`/`matricize` do not compile when instantiated [R].**
- Where: `MultiArrayComplex.h:84-87, 121-125`. They build `Eigen::Map<CVector/CMatrix>` from a `const dcomp*`.
- cppqed's copy had it right (`Map<const …>`); the library regressed in ad2836a, the F10 narrowing fix.
- Breaks `DensityOperator::trace/renorm` (`DensityOperator.h:57`), and therefore Master. Neither cppqed's QJMC path nor any library test instantiates it, which is why nobody noticed.

**A7. `Liouvillian.h:106,115` (`superoperatorFromJump`) calls `copy(view)`, which has no library equivalent [A].** It is also broken in other ways; see C14.

**A8. `SliceIterator.h` exists only in cppqed [V].** The library's CLAUDE.md layout and `MultiArray.h:55` refer to it as part of the library. Decide where it lives.

**A9. Result of the experiment [R].** I copied cppqed to a scratch directory and wired it to the library:
- library headers under a `multidiagonal/` prefix
- a framework-only `Traits.h` that includes `<multidiagonal/Traits.h>`
- the five constructor ports from A4
- `MultiDiagonal.cc` dropped

With that, `jc_qjmc` builds and reproduces the original binary's output to all 10 printed digits, jump time included, on a vacuum-Rabi run. The QJMC path tolerates the drift; Master is where the fixed library matters, and Master does not compile (B4).

---

## B. Physics and numerics in cppqed (present whether or not the split happens)

**B1. ODE tolerances are swapped [V].**
- Where: `ODE.h:213` calls `make_controlled(epsRel, epsAbs, …)`, but Boost's signature is `(abs_error, rel_error, …)` (checked in Boost 1.83).
- Every run so far has used abs = 1e-6 and rel = 1e-12, which is effectively a 1e-6 absolute tolerance. Small quantities (coherences, high Fock amplitudes, the decayed QJMC norm near a jump) get almost no error control.
- `--epsRel` and `--epsAbs` act on each other's tolerance, and the log mislabels both.

**B2. Integrating QJMC: the anti-Hermitian part of H_eff misses its identity-proportional term [R].**
- Mode: ½ΣJ†J = κ(2n_th+1)n + κn_th; the code has only the first term (`Mode.h:39/81`).
- Qbit: ½ΣJ†J = γ_m|1⟩⟨1| + γ_p|0⟩⟨0|; the code has (γ_m−γ_p)|1⟩⟨1|. (Correction 2026-10-19: through `qbit::make(const Pars&)` it actually has no γ term at all, see B17.)
- The integrating algorithm needs d‖ψ‖²/dt = −Σ⟨J†J⟩ exactly.
- Reproduced: a vacuum mode with n_th=1 and κ=1 makes **zero** photon-gain jumps over T=20 for 3 seeds, so ⟨n⟩ stays at 0 instead of relaxing to 1. Defaults (n_th = γ_p = 0) are unaffected.

**B3. `Master.h:72`, `conj(m.rho)`, complex-conjugates ρ after every step [V; found independently by two agents].**
- The sequence U-on-columns → `hermitianConjugateSelf` → U-on-columns already yields UρU†. The trailing `conj` makes it (UρU†)*.
- Without a propagator, ρ → ρ* every step. Coherent dynamics then cancel over pairs of steps.
- This was correct in v2 only because blitz `transposeSelf` re-strided without moving data. Fix: delete the line.

**B4. Master.h and Schrödinger.h have not compiled since 2bc5d3f8 ("restructuring started", 2024-04-30) [R].**
- They include the deleted `QuantumSystemDynamics.h` and use `getHa/getLi/getEx/getEV/getFreqs` and `hamiltonian_ns::broadcast`, none of which exist. Nothing includes them, so the build never noticed.
- So there are no current master-equation results to doubt. The Bloch test cannot run until Master is ported to the `(ha, li, ev[, ex])` API QJMC uses.
- Library finding F5 is also not the only suspect: B3, B5 and B7 each corrupt Master independently.

**B5. Mode superoperators get √rate instead of rate [V].**
- Where: `Mode.cc:81, 97` call `aSuperoperator(sqrt(fact))` and `aDagSuperoperator(sqrt(fact))`, but the superoperator multiplies by its argument once. JρJ† needs 2κ(n_th+1) (and 2κn_th).
- Example: κ=10 recycles 4.47 instead of 20.
- Qbit's superoperators are correct. Master only.

**B6. `Qbit.h:60` builds the gain channel with the loss rate [V].** It calls `gain(gamma_m)` where it should call `gain(gamma_p)`.

**B7. `SliceIterator` first-slice offset bug [R] and read past the end [R].**
- `SliceIterator.h:159` initializes the first slice with `*i` and drops `originalOffset`.
  - Probe: slicing a view with offset 1 gives `0 2 4` for slice 0 instead of `1 3 5`; later slices are right.
  - Affected: Master with any composite Hamiltonian (column slices), nested Broadcasters, `partialTrace` of sliced LDOs.
  - Not affected: top-level QJMC, so `jc_qjmc` is fine.
- `:163` dereferences `++iter` at end on every complete loop. ASan reports a heap-buffer-overflow.

**B8. `Composite.h:35`: the rate lambda takes `StateVectorConstView<RANK-size(ra)>`; it should be `<size(ra)>` [V].** It compiles only when RANK = 2·|ra|, as in JC, so any 3-subsystem QJMC with broadcast Liouvillians fails to build.

**B9. `Composite.h:15`: Broadcaster offsets assume the strides of a contiguous array [A].** Nesting across ranks (a Broadcaster inside a Broadcaster on a sub-view) touches the wrong elements. Together with B7 this contradicts "nesting composites works to any depth" in CLAUDE.md.

**B10. Broadcaster lambdas capture their arguments by reference (`&l`, `&s`) [V pattern, A analysis].**
- Where: `Composite.h:28, 32, 40, 76, 81`.
- Broadcasting a temporary, or returning a composed system from a helper, dangles. jc_qjmc is safe only because main's structured bindings outlive `run()`.

**B11. QJMC never normalizes the initial ψ [A].**
- Where: `QuantumJumpMonteCarlo.h:200/289`, with `bisect_.n1 = 1`.
- A truncated coherent state with ‖ψ₀‖² < 1 gives unnormalized expectation values at t=0, a biased first jump time, and possibly a spurious jump on the first step.

**B12. Integrating QJMC has no exact-propagator hook and no interaction-picture reset [A].**
- Where: line 226 is commented out and there is no `ex` member.
- Harmless today, but needed for roadmap item 1:
  - U must be applied before the norm/bisection check.
  - Backward bisection steps need U(τ<0).
  - `time0` is currently 0 on a fresh run but t_restart after a continuation.

**B13. The stepwise question in CLAUDE.md [A, physics checked].** Stepwise needs the full non-unitary H_eff, not only the unitary part.
- Constraint: the anti-Hermitian part must appear exactly once, in either `ham` or U.
- If U carries it, `−z·n` (and Qbit's `diagonalH`) must leave the ODE Hamiltonian, and drive and coupling terms then need interaction-picture phases.
- Other stepwise compile bugs besides the known `q.ex`/`q.qsd`:
  - `calculateRates(this->li)` has no matching overload.
  - `if (!time)` binds to `::time` from `<ctime>`, so it is always false.

**B14. Master with a non-unitary U (the Mode roadmap) [A, forward-looking].** Time-independent superoperators need interaction-picture factors (e^{−2κ(2n_th+1)τ}). Without them the error is first order, which RKCK54 error control cannot see. v2 refused non-unitary U in Master.

**B15. Conventions [V/A]:**
- `qbit::sz = diag(+½,−½)` gives [sz,s⁺] = −s⁺, the opposite of `spin::sz`. It is unused today.
- Qbit's "polarization" is ⟨σ†⟩ while Mode's "ladder operator" is ⟨a⟩, so JC phases are mutually conjugated.

**B16. Mode gain rate at the truncation edge [A, low].** It counts n = D−1, but the truncated a† jump discards ψ_{D−1}.

**B17. `qbit::make(const Pars&)` passes the real detuning as the complex frequency [A; found 2026-10-19 during the docs fix].**
- Where: `Qbit.h:102` calls `make(p.delta, p.eta, …)`, and the first parameter is `dcomp zSch`, which `diagonalH` uses as `dpsidt(1) -= z*psi(1)`.
- So z = δ (real) instead of z = γ_m − iδ, the analogue of Mode's z = κ(2n_th+1) − iδ (`Mode.h:48`).
- Consequences: the detuning is missing from H; γ_m does not enter H_eff, so the no-jump norm decay of the integrating QJMC is wrong; and with the default δ = −10 the excited-state amplitude grows like e^{10t} between jumps.
- Affects every Qbit run, including `jc_qjmc` (default `--deltaq -10`). Not yet reproduced numerically.
- The fix overlaps with B2 (the γ_p|0⟩⟨0| term cannot be expressed by `diagonalH`, which only touches level 1) and B6.

---

## C. Production driver and robustness (the HPC use case)

**C1. Continuation silently restarts from t=0 [R].**
- Cause: `<o>.state` is created at run start, "continuing" is decided only by `<o>` being non-empty, and an empty `.state` is accepted.
- Reproduced:
  1. run to T=1
  2. truncate `.state` (equivalent to a walltime kill with the default `sdf=0`)
  3. resubmit with T=2
- Result: "Continuing from time 0", the times in the file go 1 → 0.5, and the exit code is 0.
- Variants: a stale `.state` from an earlier run with the same `--o`; with `sdf>0`, the trajectory file runs ahead of the last checkpoint.

**C2. Unknown command-line options only produce a stderr warning; the run proceeds with defaults and exits 0 [R].** A typo in a job-array sweep runs every task with the wrong parameters.

**C3. The default `--dc 10` silently overrides `--Dt` and `--NDt` [R].**

**C4. The `--streamSwitch` bitset is read most-significant-bit first, the reverse of the help text [A].** Passing "0111" to disable the intro disables first-state serialization instead.

**C5. NDt mode is not idempotent under continuation [A].** `timeToReach = t_current + NDt·Dt`, so the final time depends on how often the job was killed.

**C6. State-file reads don't check `gcount`, and checkpoints are never flushed [A, plausible].** A partial last record either kills continuation or silently zero-fills fields.

**C7. ODR violations in headers [A].**
- Non-inline definitions: `ODE.h:117 defaultLogger`, `Trajectory.h:230 observerNoOp`, `QJMC.h:62 defaultLogger`, and the non-inline functions in `Mode.h`/`Qbit.h`.
- These link only because of shared-library interposition. The first multi-TU Catch2 test binary will hit "multiple definition".

**C8. RNG streams [A].**
- QJMC seeds `pcg64{seed, p}`, but `streamOfOrdo` uses `{seed, 1001+p}`.
- Once ensembles return, consecutive ordinals overlap between job-array tasks.
- Not the mt19937 seed+i pattern, and RNG state serialization is exact.

**C9. Ensembles, which are disabled, have further problems beyond the known ones [A]:**
- `AverageTrajectories::_` returns a dangling `const auto&`.
- It averages unnormalized ψ (`averaged(q)` returns raw ψ; `addTo` ignores `weight`).
- Any Ensemble `run()` fails to compile, through `InitializeEnsembleFromArrayOnlyArchive`.
- A qualified `cppqedutils::advance` blocks nested ensembles.
- `stateIO` does not record nTraj.

**C10. Autostop, beyond the known problem [A]:**
- The no-observer `run(traj,p)` overload never compiles (`ranges::fill` on a Hana tuple).
- `replaceZeros` recurses on `t` instead of `v`.
- `valarray<bool>` appears in an `if`.

**C11. ODE edge cases [A]:**
- A backward first step sets dtTry = 0, so `advance()` loops forever.
- A NaN error estimate is accepted.
- Nothing stops someone plugging in an FSAL stepper, which would go stale after jumps, renormalization or the interaction-picture reset.

**C12. QJMC bisection drops the `LogTree` of its recursive `step()` calls (lines 255, 260) [A].** This is the likely cause of the known "no in-trajectory log in DC mode for dc>1".

**C13. `Pars` objects aren't safe to copy or move [A, plausible].** popl stores `&binding`. Copy and move should be deleted.

**C14. Dead or stale code that would fail once instantiated [A]:**
- `Simulated.h`
- `Profiling.h` (`#error` in Boost 1.83)
- `TemporalDataPoint.h:49` `label`
- `EntanglementMeasures.h` (wrong members; `negPT` returns 0; declared functions never defined)
- `SliceIterator.h:257` `transpose` (unqualified `filterIn`)
- `Sigma.h` (v2 types; included by `MLJC.h`)
- `DistributionFunctions.*`, `LazyDensityOperatorFFT.*`
- `Mode_.h`, `CoupledModes.h` (include files that no longer exist)
- `CoupledModes.cc` (fully commented out)
- `superoperatorFromJump` (`auto&` binding to a zip prvalue; ignores offset/strides; assumes a Hermitian input block)

**C15. Structure layer [A]:**
- `SparseMatrix.h:27` ignores offset/strides.
- `SparseMatrix.cc:9` constructs a `range_error` without `throw`.
- A NoOp Hamiltonian hits `static_assert(always_false)`.
- A `Lindblad` without `.rate` throws `bad_function_call`.
- The unconstrained variadic `LDO` constructor makes overloads ambiguous.
- The `ExpectationValues` concept is checked only for state vectors.
- 3-argument H(t,…) silently receives τ = t−t0 ∈ [0,Δt].
- `static const` namespace-scope objects are a formal ODR issue.
- Propagators have no Hana-sequence or Broadcaster branch, so a composite cannot carry several propagators. This blocks roadmap 1 for Qbit+Mode.
- `concatenate` fails because the code relies on P2280 (GCC ≥ 14); derive the rank from `Lindblad::N_RANK` instead.

**C16. Performance [A]:**
- Jump/rate/superoperator variants are passed by value, so the `matrixOffsets` cache is rebuilt on every RHS call.
- Each checkpoint copies the stepper (6-7 state-sized buffers).
- `.state` grows without bound, and continuation deserializes every record.
- Unused `find_package(Boost … iostreams)` can break configure on minimal cluster Boost installs.

---

## D. Library repo (cppqed-multidiagonal)

**D1. Const `vectorize`/`matricize` don't compile; no test instantiates them.** Same as A6 [R].

**D2. Removing `conj(diag)` from `MultiDiagonal::hermitianConjugate` passes all 49 tests [R, mutation test].** Every conjugated operator in the tests is real. CLAUDE.md also credited this to finding F4, which is about `hermitianConjugateSelf` (corrected 2026-10-15).

**D3. Flat `include/` and install; same as A2 [V].** Install rules also run when the library is a subproject [A].

**D4. Warnings and CI.**
- Clang Debug emits `-Wunused-lambda-capture` (`MultiArray.h:135`) and `-Wignored-qualifiers` (`:254`), against the "warning-free on both compilers" rule [R].
- CI has no `-Werror`, so regressions of findings F10 and F11 would pass [A].

**D5. No `target_compile_features(... INTERFACE cxx_std_23)` [A].** Consumers on C++20 get obscure errors.

**D6. Versions not pinned [A].**
- No version on `find_package` for Eigen and Boost; the README says ≥3.4 and ≥1.80.
- `SameMajorVersion` on 0.x lets a request for 0.1 accept the API-breaking 0.2.0.
- `ARCH_INDEPENDENT` is missing.
- Catch2 is pinned by tag and has no `FIND_PACKAGE_ARGS`, so a standalone offline configure fails.

**D7. README and mainpage [A]:**
- The offset formula `|n+i⟩⟨n|` with i>0 is a *lower* element, but the code means upper, so a reader following it builds a† for a.
- "(a†)²a² are further off-diagonals" is wrong: it is n(n−1), on the main diagonal.
- The performance claims overstate what is measured:
  - picture1's `multidiagonal_*` variants never call `operator()`
  - "eigen_sparse_pernnz" is a hand-written CSR loop
  - multi-axis application is 2–2.7× *slower*, which isn't mentioned
- mainpage: says C++20 / GCC ≥ 12; names `InteractionPictureDiagonal`, which doesn't exist.
- `Doxyfile`: says 0.1.0.

**D8. Header docs [V]:**
- `MultiDiagonal.h:49,79` says `operator()` "applies H/i"; it only accumulates diag·ψ.
- `:43-46` mentions a "derived type in CPPQEDcore" with a parallel frequency MultiDiagonal, which contradicts the decided design boundary.
- `:28` uses "homomorphism", against the style rule.
- `MultiArrayComplex.h:182`: the docstring still describes the pre-F5 loop.

**D9. `twoTimesRealPartOfSelf`/`hermitianConjugateSelf` index flat storage as i + N·j and ignore view offset and strides [V].** That is fine for Master's full ρ but wrong on sub-blocks. Assert or document it.

**D10. Test gaps [A]:**
- `directProduct`'s custom-combiner check cannot fail (real input).
- `isEqual` is only tested on equal arrays.
- The serialization round trip doesn't check strides/offset.
- Rank ≥ 3 application and composition are untested, and the "rank-general" test loops over a single element.
- `constructJC` has no validation of its own.
- `[contract]` tags only cover F3–F6.
- The findings numbering differs between CLAUDE.md and the test names (F1b, F2); now reconciled: the ledger is F1, F1b, F2…F11.
- `MultiDiagonal` JSON is one-way only (no `value_to`).
- No complex-valued conjugation test, no sanitizer job, no `find_package`/FetchContent consumer job.

**D11. The library's CLAUDE.md layout section is wrong [V].** It describes an `include/multidiagonal/` subdirectory and a `SliceIterator.h` that don't exist.

**D12. Performance: tested and refuted [R].** I suspected the application loop advances the last axis fastest while axis 0 is contiguous. Flipping the loop order changes `applyJC` by −15%…+33%, about the noise level; Eigen's own timings varied ±20% between runs. Multidiagonal stays 2–2.7× slower than Eigen sparse either way, so the existing diagnosis (per-element address recomputation) stands.

**D13. Smaller header items [A/V]:**
- The unconstrained `explicit MultiDiagonal(auto&&...)` beats the deleted copy constructor for non-const lvalues.
- A zero-extent diagonal (|o| = D) is accepted, then indexes out of bounds on non-leading axes.
- `sqrAbs` = `sqr(abs(v))` is inexact; use `std::norm`.
- `noInit` actually zero-fills.
- `twoTimesImagPartOf` returns A − A† without documenting it.
- Internal-linkage `noInit`, `zeroInit` and `_` are used from external templates (formal ODR).

**D14. Benchmarks always build in a top-level build, with `-march=native` by default and no switch to turn them off [A].**

---

## E. Repo hygiene and CLAUDE.md corrections (cppqed)

**E1. Files outside version control [R].** `dockerEntry/jc_qjmc.cc`, the only end-to-end executable, is gitignored by the `dockerEntry/*.cc` rule, and so is `benchmarking/`. A fresh clone builds no executable.

**E2. cppqed CLAUDE.md is wrong or overstated in places:**
- Dependencies are git submodules, not FetchContent [R].
- The executable is `dockerEntry/jc_qjmc.cc`, not `JaynesCummingsQJMC.cc` [R].
- `Composite.h` is in `CPPQEDcore/structure/` [R].
- Qbit's propagator is not "already wired up": `qbit::make` doesn't use it [A].
- The `std::ranges::fill(dpsidtRaw,0)` item is functionally correct and not a live bug [A, two agents].
- `SliceRangeOwning` "move leaves b_ dangling" is overstated: vector move keeps iterators; the real SliceIterator bugs are B7 [A, two agents].
- "Master results unverified" understates it: Master doesn't compile (B4), and F5 is one of four suspects.

**E3. The version banner still says "2006–2023" [R].**
