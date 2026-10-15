# CLAUDE.md — C++QED v3

C++QED is a framework for simulating open quantum dynamics: quantum-jump Monte Carlo
(QJMC) trajectories, the Lindblad master equation, and ensembles. It targets
**production runs on HPC/HTC clusters**: you compile one executable per physical system
and run it as thousands of job-array tasks.

- It is deliberately *not* an interactive prototyping tool (that is QuTiP's niche).
- v3 rewrites v2's OO/template-metaprogramming design in a C++20 functional style
  (concepts, ranges, Boost.Hana) while keeping v2's performance core: compile-time
  resolution of the tensor-product structure.
- Development resumed in 2026 after a pause in 2023. The first end-to-end executable ran
  in April 2026.

License: BSL-1.0. Developer Guide: `docs/`, published at https://vukics.github.io/cppqed/.

## Repository relationships

- **The array and operator foundations live in a separate repo,
  `cppqed-multidiagonal`.** That covers `Traits.h`, `MultiArray.h`, `MultiArrayComplex.h`,
  `SliceIterator.h`, and `MultiDiagonal.h`.
  - Changes to those files go **there**, not here. If a task needs one, stop and say so.
  - That repo has its own CLAUDE.md.
- This repo's `Traits.h` holds framework-only extras and includes
  `<multidiagonal/Traits.h>`. Never duplicate symbols that the library already defines
  (`passByValue_v`, `ReferenceMF`, …).
- Dependencies come in via FetchContent pinned to **full commit SHAs**: cppqed-multidiagonal
  plus forked popl, xoshiro, and pcg. Each fork should document why it diverges from
  upstream.
- `-DMULTIDIAGONAL_LOCAL=<path>` points the build at a sibling checkout of
  cppqed-multidiagonal during co-development.
- Cluster compute nodes have no outbound network. Builds there use
  `FETCHCONTENT_FULLY_DISCONNECTED=ON` and `FETCHCONTENT_TRY_FIND_PACKAGE_MODE=OPT_IN`
  with system Boost and Eigen. Don't introduce dependencies that break offline builds.

## Layout

```
CPPQEDutils/      physics-agnostic scaffolding: ODE.h (Boost.Odeint, RKCK54), Trajectory.h (run()),
                  StochasticTrajectory.h (Ensemble, AverageTrajectories), Pars.h (popl + JSON), Random.h
CPPQEDcore/
  quantumdata/        StateVector, DensityOperator, LazyDensityOperator (LDO)
  quantumoperator/    framework-side operator helpers (MultiDiagonal itself comes from the library)
  quantumtrajectory/  QuantumTrajectory.h, Master.h, QuantumJumpMonteCarlo.h — header-only
  structure/          Hamiltonian, Liouvillian, ExactPropagator, ExpectationValues (planned rename → quantumdynamics)
CPPQEDelements/
  frees/          SubSystems: Mode, Qbit (planned rename Free → SubSystem)
  interactions/   Composite.h (Broadcaster), JaynesCummings.h
```

The reference executable is `JaynesCummingsQJMC.cc` (qubit on axis 0, mode on axis 1).
It shows the whole user-facing API in about 40 lines.

## Build

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Debug -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
      [-DMULTIDIAGONAL_LOCAL=../cppqed-multidiagonal]
cmake --build build --target <executable>
```

- `compile_commands.json` is needed by clangd. Development happens in VS Code on WSL2
  Ubuntu.
- Compile as C++23, which the multidiagonal headers require.
- Expect long template-heavy builds. Build single targets rather than `all` while
  iterating.

## Design principles: preserve these in every change

- **Compile time first.** Everything known at compile time is resolved at compile time:
  - RANK
  - retained/dummy axes (`retainedAxes<i...>` NTTPs)
  - slice ranks
  - Hamiltonian composition (Hana sequences)
  - ODE/streaming modes (template parameters of `run()`)

  Don't turn any of these into runtime values. Hilbert-space *dimensions* are runtime
  (`Extents<R>`) for now.
- **Functional, not OO.** Use free functions constrained by concepts, and no base-class
  hierarchies.
  - A SubSystem is not an object. `make(pars)` returns a plain tuple
    `(cutoff, hamiltonian, liouvillian, expectationValues, descr)`.
- **Views vs owners.** Views are passed by value; owners are move-only. Watch for
  dangling spans/iterators after moves and reallocations.
- **The H/i convention.**
  - Hamiltonian functionals *add* `H/i · ψ` to `dψ/dt`.
  - H is generally **non-Hermitian**: in QJMC and Master contexts it is
    `H_eff = H_phys − (i/2) Σ J†J`.
  - Three time-dependence levels: `H(psi,dpsidt)`, `H(t−t0,psi,dpsidt)` (interaction
    picture), and `H(t,psi,dpsidt,t0)`.
- **Exact propagators.** These are per-subsystem diagonal propagators U that are **not
  necessarily unitary**. For Mode, `U = diag(1, e^{−zt}, …)` with
  `z = κ(2n_th+1) − iδ`.
  - `applyPropagator` *replaces* ψ by Uψ, whereas Hamiltonians *add*.
- **Continuous interaction-picture reset.** This is numerically essential:
  - After every ODE step, the state is transformed back to the Schrödinger picture and
    `t0 ← t`.
  - Without the reset, matrix elements grow like `e^{NκT}`, which is `e^{10000}` for
    N=100, κ=10, T=10.
  - With it, growth stays bounded by `e^{NκΔt}`. Never remove or defer the reset.
- **The interaction picture belongs to the trajectory, not to operators.** Operators
  (MultiDiagonal) carry no frequencies.
  - v2's `dynamic_cast`-based picture-aware coupling (`expINKX(particle::Ptr, …)`) is
    deliberately not reproduced.
  - Where phase tabulation is needed, use a separate `(alpha, freqs)` MultiDiagonal pair
    in the physics layer.
- **Master equation reuses state-vector code.**
  - It broadcasts `H_eff/i` over the rows of ρ, then applies `twoTimesRealPartOfSelf`
    (X + X†) to obtain the full Lindblad equation minus the recycling terms.
  - The recycling terms `J ρ J†` come from the superoperator loop.
  - ρ is re-symmetrized and renormalized after each step.
- **LDO abstracts over the state type.** Expectation values take a
  `lazy_density_operator`, so the same functional works for a pure state, a mixed state,
  and an ensemble average.
- **Composition.**
  - `Broadcaster<RANK, ra...>` lifts rank-|ra| descriptors to RANK.
  - Hamiltonians and ExpectationValues are recursive Hana sequences, so nesting
    composites works to any depth.
  - Liouvillians are flat `std::vector<Lindblad<RANK>>`; nest them by broadcast plus
    concatenation.
  - ExactPropagators are not composed: each leaf contributes one, broadcast to the full
    rank.
- **RNG.** Use `pcg64` with independent stream ordinals for ensembles. **Never use
  `std::mt19937` with seed+i.**
- **Parameters.**
  - `Pars` structs register through `add(mod, op, tjc, "Title", _("name","desc",default,binding), …)`.
  - `mod` suffixes namespace the options (`--cutoffm`, `--deltaq`).
  - `noJSON` excludes initial conditions from the system descriptor.
- **run() is the production interface.**
  - Options: `--T`/`--NDt`, `--Dt`/`--dc`, `--o`, `--sdf`, `--initialFileName`,
    `--streamSwitch`, `--precision`.
  - Runs continue automatically from `<o>.state`, and every output file documents itself.
  - Keep state files binary-exact; continuation requires full precision.

## Current state and known problems

Treat **master-equation results as unverified**.

- **Suspected cause:** finding F5 in the library (`twoTimesRealPartOfSelf` quadrupled the
  diagonal). That would mean doubled Hamiltonian flow on the diagonal plus artificial
  dephasing from the post-step symmetrization.
- **Confirmation:** a driven-qubit Bloch-equation test comparing three results:
  - Master
  - QJMC ensemble
  - the analytic optical-Bloch steady state

  Run it with patched versus unpatched headers.

Other known problems:
- `UnaryDiagonalPropagator` operator bodies are unimplemented. This is the highest
  priority, because it blocks exact propagation for Qbit (already wired up), Mode (where
  `propagator()` needs reinstating), and the future MDoF subsystem.
- QJMC `stepwise` is broken: it has stale `q.ex`/`q.qsd` references.
  - Open question before reinstating it: stepwise evolves under the physical H and
    accounts for decay through rates. It probably needs only the **unitary part** of the
    propagator, `diag(e^{iδnτ})`. Using Mode's decaying U there would double-count
    dissipation. v2's old code holds the intended answer.
- `makeEnsemble` is commented out, so ensemble runs are disabled.
- Autostop is **broken**: `AutostopHandlerGeneric` assumes a numeric TDP, but TDPs are
  labelled Hana tuples.
  - Callers must pass `trajectory::observerNoOp` explicitly. The default argument is dead
    code because template deduction ignores defaults.
- `structure::concatenate` fails with non-constexpr arguments. Concatenate Liouvillians
  manually.
- `std::ranges::fill(dpsidtRaw, 0)` in QJMC `step()` misuses the returned iterator.
  Separate the fill from the view construction and use `dcomp{}`.
- `SliceRangeOwning` move leaves `b_` dangling into the moved-from vector.
- Stubs:
  - `compose()` in Hamiltonian.h (interactions across composite boundaries)
  - `broadcastFor`
  - `Placeholder` slicing beyond RANK==2
  - `InitializeEnsembleFromArrayOnlyArchive`
- In-trajectory log in DC mode doesn't work for `dc > 1`.

## Roadmap: migration from v2, not new features

1. Interaction-picture propagators through Mode, Qbit, and MDoF.
2. Motional degree of freedom (v2 `Particle`/`ParticleCavity`).
   - Momentum-space representation, `Spatial` grids.
   - The kinetic propagator `exp(−i ω_rec k_j² τ)` becomes a `UnaryDiagonalPropagator`.
   - `expINKX(nK, dim)` becomes a single-diagonal MultiDiagonal at offset `+nK`, and
     `cos`/`sin` follow from it.
   - The 80-line `mfComposition` table is replaced by `operator|`.
3. Restore autostop (`hana::for_each` over the TDP, or flatten it to
   `valarray<double>`).
4. Re-verify the master equation (Bloch test above).
5. Documentation and tutorials.

Later or under discussion, so **don't start these unasked**:
- renames (`structure/` → `quantumdynamics/`, Free → SubSystem, acronym casing)
- compile-time leading subsystem dimensions
- the `OBSERVER` concept unification in Trajectory.h
- HDF5 serialization

## Testing

The framework-level suite is thin. Prefer Catch2 v3, as the library does. Test physics
against analytic results:
- Mode: coherent-state evolution, thermal state
- Qbit: Rabi oscillations
- Jaynes–Cummings: vacuum Rabi splitting
- Master: trace and Hermiticity preservation
- QJMC: ensemble convergence to Master

Operators are truncated, so compare against the truncated dense result, never against
textbook commutators.

## Style and conventions

- camelCase. Keep established physics abbreviations (`dcomp`, `noInit`). Use `FWD(x)`.
  Use `slv<"label">(v)` for labelled output values.
- Doxygen `///` comments. Each file starts with a copyright plus BSL-1.0 line.
- Use physics language in docs and comments ("quanta moved", "phase belongs to the
  diagonal"), not group-theory vocabulary. κ may be set to 1 as the natural unit.
- Be honest about what works. Don't advertise broken features (e.g. autostop) in docs or
  examples.
- Physics facts to keep straight:
  - Generators (a, a†, n, Kerr, two-photon drive, JC coupling, e^{inkx}) are
    multidiagonal.
  - The unitaries they generate (displacement, squeezing) are **not**.
  - The developer guide's operator list (primer §7.11) still has this wrong and needs
    fixing.

## Working rules (Wigner Research Centre institutional policy)

- **Network, credential, and out-of-tree filesystem commands:** show the command, explain
  it, and let András run it. This covers:
  - `git push/fetch/pull`
  - configure steps that download
  - package installs
  - anything outside this repo
  - anything touching keys or tokens

  Building, testing, and editing inside the repo are fine.
- **Destructive operations:** ask for confirmation first. This covers `rm -r`,
  `git reset --hard`, `git clean`, force-push, and rebases of pushed history.
- **Simulation and research data:** process locally only. No telemetry, and no external
  APIs or services unless explicitly requested.
- **Secrets:** never hardcode them.
- **System safety:** never disable antivirus, firewalls, updates, or the VPN, and never
  change network configuration.
