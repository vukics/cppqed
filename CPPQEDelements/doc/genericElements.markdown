Table of supported general-purpose elements {#genericelements}
===========================================

\tableofcontents

\par Abstract
This is a list of general-purpose free and interaction elements (residing in `CPPQEDelements`), that are supported in the framework, meaning that they will always remain an integral part of the framework’s distribution. The supported scripts (residing in `CPPQEDscripts`) rely on them so that they are tested by the testsuite(s). The exact list of these elements does have some historical determination and the list may be extended and supplemented in the future. For how to add custom elements to the framework and write custom scripts using them cf. the directories `CustomElementsExample` and `CustomScriptsExample`.

Frees {#genericelementsfrees}
=====

Reside in `CPPQEDelements/frees`

Mode {#genericelementsfreesmode}
----

A single harmonic-oscillator mode that might be driven and interact with a reservoir of possibly finite temperature. Notation: \f$z\equiv\kappa(2n_\text{Th}+1)-i\delta\f$.

| Class name | Hamiltonian | \f$U\f$ | Liouvillean | Displayed characteristics |
| ---------- | ----------- | ------- | ----------- | ------------------------- |
| Mode | n/a | \f$e^{i\delta\,t\,a^\dagger a}\f$ | n/a | \f$\avr{a^\dagger a},\;\text{var}\lp a^\dagger a\rp,\;\real{a},\;\imag{a},\f$ … |
| ModeSch | \f$-\delta\,a^\dagger a\f$ | n/a | ” | ” |
| PumpedMode | \f$i\lp\eta a^\dagger\,e^{-i\delta\,t}-\hermConj\rp\f$ | = Mode | ” | ” |
| PumpedModeSch | \f$-\delta\,a^\dagger a+i\lp\eta a^\dagger-\hermConj\rp\f$ | n/a | ” | ” |
| LossyMode | n/a | \f$e^{-z\,t\,a^\dagger a}\f$ | \f$2\kappa\lp(n_\text{Th}+1)\,a\rho a^\dagger+n_\text{Th}\,a^\dagger\rho a\rp\f$ | ” |
| LossyModeUIP | \f$-i\kappa(2n_\text{Th}+1)\,a^\dagger a\f$ | = Mode | “ | ” |
| LossyModeSch | \f$-iz\,a^\dagger a\f$ | n/a | ” | ” |
| PumpedLossyMode | \f$i\lp\eta a^\dagger\,e^{z\,t}-\eta^* a\,e^{-z\,t}\rp\f$ | = LossyMode | ” | ” |
| PumpedLossyModeUIP | LossyModeUIP + PumpedMode | = Mode | ” | ” |
| PumpedLossyModeSch | LossyModeUIP + PumpedModeSch | n/a | ” | ” |

Qbit {#genericelementsfreesqbit}
----

Qbit has the same versions as Mode, but finite temperature is not yet implemented in this case. Notation: \f$z\equiv\gamma-i\delta.\f$

|        | Hamiltonian | \f$U\f$ | Liouvillean | Displayed characteristics |
| ------ | ----------- | ------- | ----------- | ------------------------- |
| Qbit …  | \f$-iz\,\sigma^\dagger\sigma+i\lp\eta\sigma^\dagger-\hermConj\rp\f$ | depends on used picture | \f$2\gamma\,\sigma\rho\sigma^\dagger\f$ | \f$\rho_{00},\;\rho_{11},\;\real{\rho_{10}},\;\imag{\rho_{10}}\f$ |

\note The Qbit algebra is identical to a Mode with `cutoff=2`, so that we can reuse much code from Mode when implementing Qbit

Spin {#genericelementsfreesspin}
----

Spin is characterized by a fixed magnitude \f$s,\f$ whereupon its quantum numbers are \f$m=-s\dots s,\f$ so that its dimension is \f$2s+1\f$. In the framework, all state vectors are indexed starting with 0, so that the element \f$\Psi(0)\f$ corresponds to the basis vector \f$\ket{-s}\f$. Notation: \f$z\equiv\gamma-i\delta.\f$

| Class name | Hamiltonian | \f$U\f$ | Liouvillean | Displayed characteristics |
| ---------- | ----------- | ------- | ----------- | ------------------------- |
| Spin | n/a | \f$e^{-z\,t\,S_z}\f$ | n/a | \f$\avr{S_z},\;\avr{S_z^2},\;\real{S_+},\;\imag{S_+},\f$ … |
| LossySpin | n/a | ” | \f$2\gamma\,S_-\rho S_+\f$ | ” |
| SpinSch | \f$-z\,\boldsymbol{\theta}\cdot\mathbf{S}\f$ | n/a | n/a | ” |

\note \f$\boldsymbol{\theta}\cdot\mathbf{S}\f$ is in general not diagonal in the eigenbasis \f$\ket{m}\f$ of \f$S_z,\f$ so that it would not be convenient for defining an interaction picture.

Particle {#genericelementsfreesparticle}
--------

These elements could be more accurately called a 1D motional degree of freedom.

The basic Hamiltonian \f$H=p^2/(2\mu)\f$ is most conveniently implemented in momentum basis. Discrete k-basis amounts to finite quantization length in x-space. Our choice of units is such that the smallest momentum is \f$\Delta k=1\f$, so that the quantisation length in x-space is \f$2\pi\f$. The use of discrete k-basis entails periodic boundary condition in x-basis. Spatial resolution is an integer power of 2 to be able to perform \link fft radix-2 FFT\endlink.

Notation: recoil frequency \f$\omrec\equiv\hbar\,\Delta k^2/(2\mu)=1/(2\mu),\f$ k-operator \f$k\equiv p/(\hbar\,\Delta k)=p\f$. Hence the basic Hamiltonian: \f$H=\omrec k^2.\f$

The Particle elements have conservative dynamics.

| Class name | Hamiltonian | \f$U\f$ | Displayed characteristics |
| ---------- | ----------- | ------- | ------------------------- |
| Particle | n/a | \f$e^{-i\omrec t\,k^2}\f$ | \f$\avr{k},\;\text{var}(k),\;\avr{x},\;\text{dev}(x)\f$ … |
| ParticleSch | \f$\omrec k^2\f$ | n/a | ” |
| PumpedParticle | \f$\vclass\abs{m\Int(x)}^2\f$ | = Particle | ” |
| PumpedParticleSch | \f$\omrec k^2+\vclass\abs{m(x)}^2\f$ | n/a | ” |

Here, \f$m(x)\f$ is the mode function of the pump, which can be \f$\sin(Kx),\;\cos(Kx),\;e^{\pm iKx}\f$ with arbitrary integer \f$K\f$.

\note Simulation of moving particles is inherently hard, since the Schrödinger equation is a partial differential equation, and we inevitably have to deal with both position and momentum representations, which are linked by \link fft Fourier transformation\endlink. In quantum optics, however, the particles are mostly moving in potentials created by electromagnetic fields, mainly standing and running waves. In this case we can stay in momentum space during the whole time evolution. A strange consequence is that in numerical physics the harmonic oscillator seems to be hard, while the cosine potential is easy.

PumpedLossyMultiLevelSch {#genericelementsfreesmultilevel}
------------------------

\see \ref multilevelbundle

Interactions {#genericelementsinteractions}
============

Reside in `CPPQEDelements/interactions`

All the operators are automatically taken in interaction picture, if the underlying free element is in interaction picture.

| Class name | Free elements | Hamiltonian | Displayed characteristics |
| ---------- | ------------- | ----------- | ------------------------- |
| JaynesCummings | (Qbit / Spin) – Mode | \f$i\lp g^*\sigma a^\dagger-\hermConj\rp\f$ | n/a |
| GeneralDicke | Mode – Spin | \f$\displaystyle u\,a^\dagger a\lp S_z+\frac s2\rp+y\lp a+a^\dagger\rp S_x\f$ | n/a |
| NX_CoupledModes | Mode – Mode | \f$u\,a^\dagger a\lp b+b^\dagger\rp\f$ | n/a |
| QbitModeCorrelations | Qbit – Mode | n/a | \f$\real{\avr{\sigma a^\dagger}},\;\imag{\avr{\sigma a^\dagger}},\;\real{\avr{\sigma a}},\;\imag{\avr{\sigma a}},\;\real{\avr{\sigma_z a}},\;\imag{\avr{\sigma_z a}}\f$ |
| ModeCorrelations | Mode – Mode | n/a | covariances of the modes’ quadratures |
| ParticleOrthogonalToCavity | Mode – PumpedParticle | \f$\text{sign}\{U_0\}\sqrt{U_0\vclass}\lp a^\dagger m(x)+\hermConj\rp\f$ | n/a |
| ParticleAlongCavity | Mode – (Pumped)Particle | \f$U_0\abs{m(x)}^2 a^\dagger a+\text{sign}\{U_0\}\sqrt{U_0\vclass}\lp a^\dagger m(x)+\hermConj\rp\f$ | n/a |
| ParticleTwoModes | Mode – Mode – Particle | \f$\sqrt{U_{01}U_{02}}\lp m_1(x)m_2(x)\,a_1^\dagger a_2+\hermConj\rp\f$ | n/a |

\see [This issue](http://sourceforge.net/p/cppqed/bugs/1/)

MLJC {#genericelementsinteractionsmultilevel}
----

<b>m</b>ulti-<b>l</b>evel <b>J</b>aynes-<b>C</b>ummings

\see \ref multilevelbundle
