The MultiLevel bundle {#multilevelbundle}
=====================

\tableofcontents

\par Abstract
The MultiLevel bundle is a framework for representing multi-level quantum systems, e.g. atoms or ions of arbitrary level structure. Coherent driving, loss, and interaction with harmonic-oscillator modes can also be defined. The framework is defined in the header files `CPPQEDelements/frees/MultiLevel.h` and `CPPQEDelements/interactions/MLJC.h`, the latter acronym standing for <b>m</b>ulti-<b>l</b>evel <b>J</b>aynes-<b>C</b>ummings.

Introduction {#multilevelintroduction}
============

An \f$N\f$-level system is characterized by a set of basis vectors
\f[\ket{i},\quad i=0\dots N-1\f]

A generic Hamiltonian with coherent driving – e.g. with laser light – between certain levels and (Jaynes-Cummings) coupling to a (driven) harmonic oscillator with ladder operator \f$a\f$ may read
\f[
H=\sum_i h_i\ket i\bra i
+i\sum_{\underset{\text{driven}}{i,j\in}}\lp\eta_{ij}\,e^{-i\omega_{ij}t}\ket j\bra i-\text{h.c.}\rp
+i\sum_{\underset{\text{coupled}}{i,j\in}}\lp g_{ij}\,a\ket j\bra i-\text{h.c.}\rp
+\omegaM\,\adagger a
+i\lp\eta e^{-i\omega t}\adagger-\text{h.c.}\rp
\f]
with \f$h_i\f$ energy of level \f$i\f$, \f$\omegaM\f$ bare oscillator frequency, \f$\eta_{ij}\f$ driving strength between levels \f$i,j\f$ with frequency \f$\omega_{ij}\f$, \f$g_{i,j}\f$ coupling strength between transition \f$i,j\f$ and the mode, \f$\eta\f$ mode driving strength with frequency \f$\omega\f$.

On passing to the frame rotating with \f$\omega\f$:
\f[
H=\sum_i h_i\ket i\bra i
+i\sum_{\underset{\text{driven}}{i,j\in}}\lp\eta_{ij}\,e^{-i\omega_{ij}t}\ket j\bra i-\text{h.c.}\rp
+i\sum_{\underset{\text{coupled}}{i,j\in}}\lp g_{ij}\,a\,e^{-i\omega t}\ket j\bra i-\text{h.c.}\rp
-\deltaM\,\adagger a
+i\lp\eta \adagger-\text{h.c.}\rp
\f]
it becomes apparent that *from the multilevel system’s* point of view, the driving and coupling terms are algebraically equivalent with \f$\eta_{ij}\,e^{-i\omega_{ij}t}\Leftrightarrow g_{ij}\,a\,e^{-i\omega t}\f$. To eliminate explicit time-dependence, we apply the unitary transformation:
\f[
U=e^{-i\,t\sum_i\Omega_i\ket i\bra i}
\f]
yielding
\f[
H=\sum_i \lp h_i-\Omega_i\rp\ket i\bra i
+i\sum_{\underset{\text{driven}}{i,j\in}}\lp\eta_{ij}\,e^{i\lp\Omega_j-\Omega_i-\omega_{ij}\rp t}\ket j\bra i-\text{h.c.}\rp
+i\sum_{\underset{\text{coupled}}{i,j\in}}\lp g_{ij}\,a\,e^{i\lp\Omega_j-\Omega_i-\omega\rp t}\ket j\bra i-\text{h.c.}\rp
-\deltaM\,\adagger a
+i\lp\eta \adagger-\text{h.c.}\rp
\f]
So the system of equations to be solved for \f$\Omega\f$s in order to eliminate time dependence reads:
\f[\Omega_j-\Omega_i-\omega_{ij}\text{ for }i,j\in\text{set of driven transitions}\f]
\f[\Omega_j-\Omega_i-\omega\text{ for }i,j\in\text{set of coupled transitions}\f]
which can be solved in many cases of interest, especially if, as is often the case in actual experiments, all \f$\omega_{ij}=\omega\f$.

The actual Hamiltonian {#multilevelactualHamiltonian}
======================

Hence we get the simplified Hamiltonian
\f[
H=-\sum_i\delta_i\ket i\bra i
+i\sum_{\underset{\text{driven}}{i,j\in}}\lp\eta_{ij}\ket j\bra i-\text{h.c.}\rp
+i\sum_{\underset{\text{coupled}}{i,j\in}}\lp g_{ij}\,a\ket j\bra i-\text{h.c.}\rp
-\deltaM\,\adagger a
+i\lp\eta \adagger-\text{h.c.}\rp
\f]
with \f$\delta_i=\Omega_i-h_i\f$.

The elements {#multilevelelements}
============

This is the Hamiltonian that can be implemented by the elements in the MultiLevel bundle, in particular,
* PumpedLossyMultiLevelSch (with corresponding maker function multilevel::makePumpedLossySch) implements the first two terms, and adds noise with jump operator \f[J_{j\rightarrow i}=\sqrt{2\gamma_{ij}}\ket i\bra j\f] corresponding to a single radiative transition between two levels with decay rate \f$\gamma_{ij}\f$. \note In case of noise the free part of the Hamiltonian becomes non-Hermitian according to the \ref mcwftrajectory "MCWF method", to wit
\f[-\delta_i\ket i\bra i\Longrightarrow-iz_i\ket i\bra i\text{ with }z_i=-i\delta_i+\sum_k\gamma_{ki}\f]
* MLJC implements the third term.

\note The bundle is designed in such a way that a set of drivings (\f$\eta_{ij}\f$), couplings to harmonic-oscillator modes (\f$g_{ij}\f$), and decays (\f$\gamma_{ij}\f$) can be defined *at compile time* as lists of pairs of \f$i\f$ and \f$j\f$ levels. From these lists, the framework assembles the Hamiltonian of the system *at compile time* to a form corresponding to the last one above, and the Lindblad operators as a set of operators like above. Then, the values of parameters like \f$\eta_{ij}\f$, \f$g_{ij}\f$, and \f$\gamma_{ij}\f$ for all the pairs such defined can be specified at runtime. The rationale of this arrangement is that such pairs of levels as are never expected to be driven or coupled (or lossy), will not pollute the Hamiltonian with terms which would turn out to be zero at runtime. (Or, if the pair is not specified in the list of decays, it will not be considered when trying for quantum jumps.)

\todo Describe phase noise.

Examples {#multilevelexamples}
========

Ladder system {#multilevelexamplesladder}
-------------

In this case, three levels have increasing energy \f$h_0<h_1<h_2\f$, and we have driving (or coupling to a mode or one mode each) between 0,1 and 1,2. Then, the equations can be solved for two different driving strengths. Putting \f$\Omega_0=h_0\f$, we get \f$\Omega_1=\omega_{01}+h_0\f$ and \f$\Omega_2=\omega_{12}+\omega_{01}+h_0\f$. The Hamiltonian reads:
\f[
H=-\left[\omega_{01}-\lp h_1-h_0\rp\right]\ket 1\bra 1-\left[\omega_{01}+\omega_{12}-\lp h_2-h_0\rp\right]\ket 2\bra 2
+i\lp\eta_{01}\ket 1\bra 0-\text{h.c.}\rp+i\lp\eta_{12}\ket 2\bra 1-\text{h.c.}\rp
\f]

Lambda system {#multilevelexampleslambda}
-------------

In this case, three levels are ordered as \f$h_0\approx h_1<h_2\f$, and we have driving (or coupling to a mode or one mode each) between 0,1 and 1,2. Then, the equations can be solved for two different driving strengths. Putting \f$\Omega_0=h_0\f$, we get \f$\Omega_1=\omega_{02}+h_0-\omega_{12}\f$ and \f$\Omega_2=\omega_{02}++h_0\f$. The Hamiltonian reads:
\f[
H=-\left[\omega_{02}-\omega_{12}-\lp h_1-h_0\rp\right]\ket 1\bra 1-\left[\omega_{02}-\lp h_2-h_0\rp\right]\ket 2\bra 2
+i\lp\eta_{02}\ket 2\bra 0-\text{h.c.}\rp+i\lp\eta_{12}\ket 2\bra 1-\text{h.c.}\rp
\f]
In case of \f$\omega_{02}-\omega_{12}=h_1-h_0\f$ we get the Raman scheme, with only the upmost level rotating.

Example scripts {#multilevelexamplesscripts}
---------------

* Cf. `Raman.cc` and `CavityRaman.cc` in `CPPQEDscripts` for actual implementations.

* The script `CPPQEDscripts/Ca40InCavity.cc` demonstrates more complex usage. This script represents a \f$^{40}\text{Ca}^{+}\f$ ion interacting with a cavity mode as described e.g. in Reference \cite Barros09. In this case, the \f$4^{2}\text S_{1/2}\f$, \f$3^{2}\text D_{3/2}\f$, and \f$4^{2}\text P_{1/2}\f$ levels of the ion constitute an eight-level system, with two pumped transitions, six transitions coupled to the cavity mode, and ten decaying transitions.

