Guide on the use of the structure-bundle {#structurebundleguide}
========================================

\tableofcontents

\note The files referred to in this Section are found in directory `examples` in the distribution. The code examples are expected to compile, cf. `examples/Jamfile`.

Implementation of a class representing a harmonic-oscillator mode {#basicoscillator}
=================================================================

We demonstrate how to implement an element representing a pumped lossy mode in a truncated Fock space. In the frame rotating with the pump frequency, it is described by the Hamiltonian:
\f[H=-\delta a^\dagger a+\lp\eta a^\dagger+\hermConj\rp,\f]
where \f$\delta\f$ is the detuning between the oscillator and the pump frequencies. The Liouvillean:
\f[\Liou\rho=\kappa\lp(n_\text{Th}+1)\lp2a\rho a^\dagger-\comm{a^\dagger a}{\rho}_+\rp+n_\text{Th}\lp2a^\dagger\rho a-\comm{a\,a^\dagger}{\rho}_+\rp\rp\f]
The frequency-like parameters are \f$\delta\f$, \f$\kappa\f$, and \f$\eta\f$, representing the mode frequency, loss rate, and pumping Rabi frequency, respectively.
A dimensionless parameter is \f$n\f$ (the average number of thermal photons in the heat bath) and the cutoff of the Fock space.

Using the notation of Sec. \ref mcwftrajectory "Description of the MCWF method" \f[J_0=\sqrt{2\kappa(n_\text{Th}+1)}a,\quad J_1=\sqrt{2\kappa n_\text{Th}}a^\dagger.\f]
The non-Hermitian Hamiltonian for the Monte Carlo wave-function method reads:
\f[\HnH=\lp-\delta-i\kappa(2n_\text{Th}+1)\rp a^\dagger a+\lp\eta a^\dagger+\hermConj\rp\equiv-iz\,a^\dagger a+\lp\eta a^\dagger+\hermConj\rp,\f]
where we have introduced the complex frequency \f[z\equiv\kappa(2n_\text{Th}+1)-i\delta.\f]

The element has to be represented by a class which inherits publicly from the necessary classes in the structure namespace.
In this simple case, it is basically two helper functions returning quantumoperator::Tridiagonal instances, a constructor, and two virtual functions inherited from
ElementAveraged that have to be written.

\see The description of quantumoperator::Tridiagonal

Consider the file `ExampleMode.h`:

\snippet ExampleMode.h basic example mode

Though the NoTime tagging class as a function argument creates some redundancy, it is necessary because of the design of the structure-bundle. Note that ElementLiouvilleanStrategies
(like ElementLiouvillean) assumes that the number of Lindblads is known @ compile time, as is the case here. If this is not the case, Liouvillean has to be used instead.

This will suffice here. Let us look at the implementations in `ExampleMode.cc`:\dontinclude ExampleMode.cc
\until }),
We construct the Free base with the dimension of the system and the name-value-multiplier tuples for the frequency-like parameters of the system, which in this case are all complex
(cf. the explanation @ DynamicsBase). We use C++11 initializer lists and the DynamicsBase::CF typedef.
\until cutoff))),
We construct the time-independent TridiagonalHamiltonian base. This is greatly facilitated by the algebra and helpers of the quantumoperator::Tridiagonal class.
\warning When implementing the Hamiltonian, not \f$H\f$ itself but \f$\frac Hi\f$ has to supplied!

\until ,"photon absorption"}),
We construct the ElementLiouvilleanStrategies base, whose second template argument denotes the number of different quantum jumps, which is 2 in this case.
The constructor takes the strategies for calculating the impact of a jump on a free::StateVectorLow, and for calculating the rate from a free::LazyDensityOperator.
These strategy functions are produced from the free-standing helpers in Lines 10-14 above through argument binding. The strategies are followed by a description of
the lossy element and the decay channels. The number of descriptive labels must agree with the number of strategies.
\until ,"imag(\")"})
We construct the ElementAveraged base, with parameters necessary to produce a simple key of the quantum averages that are communicated towards the user.
Here we calculate only three such averages, the expectation value of the number operator, and the real and imaginary parts of that of the ladder operator.
\until }
With the DynamicsBase::getParsStream function we obtain a stream whereon we can write more information about the object that gets communicated towards the user
in that part of the output which summarizes the parameters of the actual run.

Next, the inherited function Averaged::average_v is implemented:
\skip average_v
\until aJumpRate
the expectation value of the photon number is calculated (where we can reuse our function `aJumpRate`, with unit loss rate).
\until imag(offdiag);
the expectation value of the ladder operator is calculated (real & imaginary parts)
\until }
\until }

The implementation of the helpers is also quite straightforward. It may come to a separate file `ExampleModeImpl.cc`:\dontinclude ExampleModeImpl.cc
\until blitz::tensor::i+1.
\warning In situations like this, it is extremely important to write `blitz::tensor::i+1.`, otherwise the `blitz::sqrt` function will operate within the integers,
which is perhaps a shortcoming of the Blitz++ library in this respect.

\until const Tridiagonal nop
\until }

Exploiting interaction picture {#basicoscillatorip}
------------------------------

In many situations, it pays to transfer to interaction picture defined by the first term of the Hamiltonian \f$\HnH\f$ above. The ladder operator in interaction picture is
\f[a\Int(t)=a e^{-zt},\f] so that the Hamiltonian reads \f[H\Int(t)=\lp\eta a^\dagger e^{zt}+\eta^*ae^{-zt}\rp.\f]

\see [These notes](http://optics.szfki.kfki.hu/~vukics/Pictures.pdf) on how to treat interaction pictures defined by non-unitary transition operators in a consistent way.

In this case, the class representing the element has to be derived from Exact as well, which provides an interface allowing for the transformation between the two pictures.
In addition, instead of TridiagonalHamiltonian`<1,false>`, we need to derive from TridiagonalHamiltonian`<1,true>` because the Hamiltonian is now time-dependent.
\note In general usage, the jump and the averages are calculated in the normal picture also in this case (cf. explanation of classes Hamiltonian, Exact, Liouvillean, and Averaged,
and furthermore quantumtrajectory::MCWF_Trajectory). This allows for reusing the same code in both pictures. (Incidentally, here the Liouvillean remains unchanged anyway.)

\dontinclude ExampleMode.h
\skip FreeExact
\until } // basic
The OneTime tagging class at the same time carries the information about the time instant to which the (diagonal) transformation operator has to be updated.
OneTime can be implicitly converted into a double.

In the implementation, the only difference from the previous case will be the constructor, because the Hamiltonian now also requires furnishing with frequencies
(cf. quantumoperator::furnishWithFreqs), and the implementation of the virtual function FreeExact::updateU.

When “furnished with frequencies”, a quantumoperator::Tridiagonal object will internally take care about the time-dependent phases appearing in \f$a\Int(t)\f$ and \f$H\Int(t)\f$:
\dontinclude ExampleMode.cc
\skip mainDiagonal
\until // PumpedLossyModeIP::average_v exactly the same as PumpedLossyMode::average_v above
FreeExact assumes that the operator transforming between the two pictures is diagonal, and the factors to update are simply its diagonal elements.
If this is not the case, Exact has to be used instead. Here, since there are also real frequency-like parameters, we have to use DynamicsBase::RF as well.

\note Since a lot of the code from the previous case can be reused here, one will usually adopt an inheritence- or class-composition-based solution to implement classes like
`PumpedLossyMode` and `PumpedLossyModeIP` (for an inheritance-based solution, cf. [below](#hierarchicaloscillator); for one based on class-composition, cf. the actual implementation of a harmonic-oscillator 
mode in the framework in `elements/frees/Mode_.h`).

Implementing an X-X interaction {#basicxxinteraction}
===============================

Let us consider the interaction described by the Hamiltonian \f[H_\text{X-X}=g(a+a^\dagger)(b+b^\dagger).\f]

The class implementing this interaction has to be derived from Interaction`<2>` because it is a binary interaction, and TridiagonalHamiltonian`<2,...>`
(note that quantumoperator::Tridiagonal is capable to represent direct products of tridiagonal matrices).

The only thing requiring some care is that once we transform some elements into interaction picture, the whole Hamiltonian is transformed, that is,
\f$a\f$ or \f$b\f$ or both may be in interaction picture. Here, for the sake of simplicity, we assume that both constituents are of the type `PumpedLossyMode`.
Hence, the Hamiltonian is in fact \f[H_{\text{X-X;I}}(t)=g(ae^{-z_at}+a^\dagger e^{z_at})(be^{-z_bt}+b^\dagger e^{z_bt}).\f]

Consider `ExampleInteraction.h`:\dontinclude ExampleInteraction.h
\skip #include
\until } // basic

`ExampleInteraction.cc` then reads \dontinclude ExampleInteraction.cc
\skip #include
\until {}
As we see, the Hamiltonian can be written in a rather straightforward way, and it internally takes care about the time-dependent phases appearing in \f$H_{\text{X-X;I}}(t)\f$,
which result from the use of interaction picture.

The free elements are stored as shared pointers in the interaction element, and it is the task of FreesProxy to turn the constant references supplied to the constructor into (non-owning)
shared pointers. Of course, the free elements have to live in a larger scope than the interaction, otherwise we may run into trouble with dangling pointers.


Using class inheritance {#hierarchicaloscillator}
=======================

For an inheritance-based solution, it pays to define a base class collecting all the common services. Consider the following snippet from `ExampleMode.h`: \dontinclude ExampleMode.h
\skip namespace hierarchical {
\until } // hierarchical
Here, instead of ElementLiouvilleanStrategies, we can rather use ElementLiouvillean, which has as many virtual functions `doActWithJ` and `rate` as there 
are jumps (indicated by the second template argument), distinguished by the tagging classes lindblad::Base::LindbladNo. It results in a compile-time error to instantiate such
a class with an argument not smaller than the number of Lindblads (since the numbering of jumps begins with 0). Via this solution we can get around the awkwardness of specifying
the jump and rate strategies for ElementLiouvilleanStrategies, while retaining a way of controlling the number of Lindblads @ compile time.

Deriving from `ModeBase`, the definition of `PumpedLossyMode` is trivial, while for `PumpedLossyModeIP`, we have to define the virtual functions inherited from FreeExact:
\skip namespace hierarchical {
\until } // hierarchical

\note The correct treatment of frequency-like parameters would require more care in this case, since `ModeBase` does not know about `delta` & `eta`

The implementations come in `ExampleMode.cc`: \dontinclude ExampleMode.cc
\skip hierarchical::
\until // ModeBase::average_v exactly the same as PumpedLossyMode::average_v above
and the derived classes:
\skip hierarchical::PumpedLossyMode
\until // PumpedLossyModeIP::updateU exactly the same as above

We may define a function using runtime-dispatching for a ladder operator either furnished or not furnished with frequencies, depending on the actual type. 
It should be implemented via dynamic casting: \dontinclude ExampleModeImpl.cc
\skip aop(const hierarchical::ModeBase&
\until }

Interaction element in the inheritance-based design {#hierarchicalinteraction}
---------------------------------------------------

Based on the above design of the mode-element and the dispatching ladder-operator function, we may define an interaction element that works correctly with either free element,
if the constructor expects arguments of type `ModeBase`: \dontinclude ExampleInteraction.h
\skip namespace hierarchical {
\until } // hierarchical

The implementation is formally equivalent to the previous: \dontinclude ExampleInteraction.cc
\skip hierarchical::
\until {}

\warning Here, it would cause a hard-to-detect physical error to use TridiagonalHamiltonian`<2,false>` instead of TridiagonalHamiltonian`<2,true>`, because in the former case,
the time update of the binary tridiagonal would not occur even with `PumpedLossyModeIP`.

Other uses of interaction elements {#otherusesofinteraction}
==================================

In the language of the framework, every element that operates on more than one quantum numbers is an *interaction*.

Hence, if we need for instance an element calculating correlations between two free subsystems (or, two quantum numbers in general), it has to be derived from Interaction 
because only an interaction element has a chance to access more than quantum numbers.

Assume we need an element calculating the dynamics of an X-X interaction between two modes as [above](#basicxxinteraction), but it also lets the user monitor the correlations 
\f$\avr{XQ},\;\avr{XP}\;\avr{YQ},\f$ and \f$\avr{YP}\f$ between the modes, where \f$X,\;Y\f$ and \f$Q,\;P\f$ are the quadratures of the two modes, respectively. The element can be
derived from the former interaction element in the [inheritance-based design](#hierarchicalinteraction): \dontinclude ExampleInteraction.h
\skip namespace hierarchical {
\skip } // hierarchical
\skip namespace hierarchical {
\until } // hierarchical
Since now we need to operate on two quantum numbers to calculate the quantum averages, we derived from ElementAveraged`<2>`, which operates on a binary quantumdata::LazyDensityOperator.
The implementation of the averaging function may read \dontinclude ExampleInteraction.cc
\skip InteractionX_X_Correlations
\until averages(3)
\until }
\until }
at this point, the `averages` array contains the real and imaginary parts of \f$\avr{a^\dagger b}\f$ and \f$\avr{a b},\f$ respectively. Note that a quantumdata::LazyDensityOperator of
arity higher than one can be indexed via the auxiliary quantumdata::LazyDensityOperator::Idx type which represents a multi-index of the corresponding arity (and reduces to a single integer
in the unary case).

Now the desired set of quantum averages can be obtained via linear operations:
\until }
