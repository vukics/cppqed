.. _structure:

=============================
The ``structure`` namespace
=============================

.. namespace:: structure

::

  namespace structure

The ``structure`` namespace comprises modules for describing quantum systems. Among them the most important is :class:`~structure::QuantumSystem`, which is an abstract interface every system has to provide to be usable with :ref:`quantum trajectories <quantumtrajectory>` like :class:`~quantumtrajectory::MCWF_Trajectory` or :class:`~quantumtrajectory::Master`. This is why all the :ref:`elementary <generalElements>` and :ref:`composite systems <composites>` are more or less directly derived from :class:`~structure::QuantumSystem`.

Much of the design here depends on the requirements of a step of the Monte-Carlo wave function method, as described in Sec. :ref:`MCWF_Trajectory`, so the reader is asked to have a look at there, too.


Most of the classes in this namespace belong to a single hierarchy, sketched in the following diagram (the broken lines signifying that the inheritance is not direct, due to some classes in between, which can be considered implementation details). [#f1]_

.. image:: figures/structure.png
   :height: 605

We have also indicated how :ref:`elements <generalElements>` like :class:`Mode` and :class:`JaynesCummings`, and :ref:`composite systems <composites>` as :class:`BinarySystem` and :class:`Composite` fit into the hierarchy.

These modules provide a lot of services for implementing new elements in the framework. For examples on how to optimally use these services, cf. :ref:`the structure-bundle tutorial below <structureTutorial>`.

.. _structureTutorial:

---------
Tutorial
---------

.. toctree::

  structureTutorial


-------------------------------
Template argument definitions
-------------------------------

.. _structureTemplates:

.. note::

  Template argument definitions:
  
  ``int RANK``
    Positive integer standing for the number of elementary Hilbert spaces (arity of the state vector)

  ``bool IS_CONST``
    Governs the constness of the corresponding class.

  ``bool IS_TD``
    Governs the time-dependence of the corresponding class.


.. _dynamicsBase:

-----------------------------------
``DynamicsBase``
-----------------------------------

:class:`~structure::DynamicsBase` provides services for dealing with frequency-like parameters, both real and constant, for all elements, frees and interactions alike, which are hence all derived from this class. Such parameters need a special treatment because all such parameters from all the subsystems (either frees or interactions) of the given physical system have to be considered as the largest frequency of the system. This largest frequency is needed for determining the initial time-step of the ODE routine.

On the other hand, these parameters (together with all others) have to be communicated towards the user when the framework summarizes the parameters of a given run. Therefore, for each such parameter, the class stores not only the value, but also the name of the parameter, plus another real number which multiplies the value of the named frequency, to give the actual frequency as appearing in the ODE. An example will make this clear. Consider the Hamiltonian of a free mode:

.. math::

  \omega a^\dagger a

In this case the parameter supplied by the user is :math:`\omega`, but the largest frequency appearing in the ODE is actually this times the dimension of the system (the cutoff of the Fock space). Hence, the tuple stored by the class for this particular frequency-like parameter will be something like::

  "omega",omega,cutoff

In the case of a pumping term:

.. math::

  \eta\lp a^\dagger+a\rp

the multiplier will be different::

   "eta",eta,sqrt(cutoff)

The class also stores an ``std::stringstream`` object, on which the constructor of the given element can write its parameters, and these will in turn be displayed when :func:`~structure::QuantumSystem::displayParameters` is called for the system. Cf. :ref:`tutorial above <structureTutorial>`.


.. py:module:: DynamicsBase.h
   :synopsis: Defines DynamicsBase in namespace structure

.. class:: structure::DynamicsBase

  .. type:: RealFreqs

    This is a standard list storing the name-value-multiplier tuples (cf. `Boost.Tuple <http://www.boost.org/doc/libs/1_45_0/libs/tuple/doc/tuple_users_guide.html>`_) for real frequency-like parameters::

      typedef std::list<boost::tuple<std::string,double,double> > RealFreqs;

  .. type:: ComplexFreqs

    The same for complex::

      typedef std::list<boost::tuple<std::string,dcomp,double> > ComplexFreqs;

  .. function:: explicit DynamicsBase(const RealFreqs& realFreqs=RealFreqs__LP__ __RP__, const ComplexFreqs& complexFreqs=ComplexFreqs__LP__ __RP__ ) 

    The constructor simply expects initializers for the lists. In practice, these are usually created with the ``tuple_list_of`` function from `Boost.Assignment <http://www.boost.org/doc/libs/1_45_0/libs/assign/doc/index.html>`_

  .. function:: double highestFrequency() const

    Calculates the fastest timescale of the given element by simply taking the maximum ``value*multiplier`` value from the lists.

  .. function:: void displayParameters(std::ostream& os) const

    Displaying parameters of the system in two steps:

    #. Converts the stored :member:`paramsStream_` to string and prints it to ``os``

    #. Calls the virtual function :func:`displayMoreParameters`.

  .. function:: std::stringstream& getParsStream()
  
  .. function:: RealFreqs& getRealFreqs()

  .. function:: ComplexFreqs& getComplexFreqs()

    Simple query functions allowing for printing on :member:`paramsStream_` and append to the list of frequencies

  .. function:: void displayMoreParameters(std::ostream& os) const

    Virtual. The default implementation is simply printing ``# name=value`` from the lists.

  .. member:: std::stringstream paramsStream_


-----------------------------------
``QuantumSystem``
-----------------------------------

.. py:module:: QuantumSystem.h
   :synopsis: Defines QuantumSystem in namespace structure

.. class:: structure::QuantumSystem

  ``template <int RANK>`` (cf. :ref:`template parameters <structureTemplates>`); inherits publicly from :class:`DimensionsBookkeeper`

  This class describes an entity that has dimensions in a Hilbert space of arity RANK, it may have some frequencies, and some parameters to communicate towards the user.

  .. type:: Dimensions

    Inherited :type:`DimensionsBookkeeper::Dimensions`

  .. function:: explicit QuantumSystem(const Dimensions& dimensions)

  .. function:: double highestFrequency() const

    Pure virtual.

  .. function:: void displayParameters(std::ostream& os) const

    Pure virtual.

---------------
``Free``
---------------

In the language of the framework, a free system is a system whit arity 1.

.. py:module:: Free.h
   :synopsis: Defines Free in namespace structure and some related types

.. class:: structure::Free

  Inherits publicly from :class:`~structure::QuantumSystem`\ ``<1>`` and :class:`~structure::DynamicsBase`

  .. function:: explicit Free(size_t dimension, const RealFreqs& realFreqs= RealFreqs__LP__ __RP__, const ComplexFreqs& complexFreqs= ComplexFreqs__LP__ __RP__)

  .. function:: double highestFrequency () const

    Simply connects the pure virtual :func:`QuantumSystem::highestFrequency` to the implementation :func:`DynamicsBase::highestFrequency`::

      double highestFrequency () const {return DynamicsBase::highestFrequency();}

  .. function:: void displayParameters(std::ostream& os) const

    The same::

      void displayParameters(std::ostream& os) const {return DynamicsBase::displayParameters(os);}


.. namespace:: structure::free

::

  namespace free

This namespace contains some type definitions for convenience in defining free systems:

.. type:: structure::free::Tridiagonal

  typedef for :class:`quantumoperator::Tridiagonal`\ ``<1>``::

    typedef quantumoperator::Tridiagonal<1> Tridiagonal;

.. type:: structure::free::Frequencies

  typedef for :class:`quantumoperator::Frequencies`\ ``<1>``::

    typedef quantumoperator::Frequencies<1> Frequencies;

.. type:: structure::free::StateVectorLow

  ::

    typedef quantumdata::Types<1>::StateVectorLow StateVectorLow;

.. type:: structure::free::DensityOperatorLow

  ::

    typedef quantumdata::Types<1>::DensityOperatorLow DensityOperatorLow;

.. type:: structure::free::LazyDensityOperator

  ::

    typedef quantumdata::LazyDensityOperator<1> LazyDensityOperator;

.. type:: structure::free::StateVector

  ::

    typedef quantumdata::StateVector<1> StateVector;

.. type:: structure::free::DensityOperator

  ::

    typedef quantumdata::DensityOperator<1> DensityOperator;


----------------
``Interaction``
----------------

This class describes interaction between free systems.

.. py:module:: Interaction.h
   :synopsis: Defines Interaction in namespace structure

.. class:: structure::Interaction

  ``template <int RANK>`` (cf. :ref:`template parameters <structureTemplates>`); inherits publicly from :class:`~structure::DynamicsBase`

  It does not inherit from :class:`~structure::QuantumSystem` because it does not make sense to simulate such an element as describes an interaction alone. However, an interaction can have frequency-like parameters, hence the inheritance from :class:`~structure::DynamicsBase`.

  .. type:: Frees

  ::

    typedef blitz::TinyVector<const Free*,RANK> Frees;

  .. function:: explicit Interaction(const Frees& frees, const RealFreqs& realFreqs=RealFreqs__LP__ __RP__, const ComplexFreqs& complexFreqs=ComplexFreqs__LP__ __RP__)

  .. function:: const Frees& getFrees() const

    The class has knowledge of its :class:`~structure::Free` constituents, and will also communicate them.



---------------------------------------
``Hamiltonian``
---------------------------------------


.. type:: structure::TimeDependence

  An enumeration of the following possibilities for time dependence::

    enum TimeDependence {TWO_TIME, ONE_TIME, NO_TIME};

  ==== ====================== ====================== =========================================================================
  Case ``TimeDependence``     The Hamiltonian
  ==== ====================== ====================== =========================================================================
  1    ``TWO_TIME``           :math:`H(t,t_0)`       Time-dependent problem + :ref:`exact part <Exact>` (:math:`U(t,t_0)`)
  2    ``ONE_TIME``           :math:`H(t)`           Time-dependent problem, no exact part
  3    "                      :math:`H(t-t_0)`       Time-independent problem + :ref:`exact part <Exact>` (:math:`U(t-t_0)`)
  4    ``NO_TIME``            :math:`H(0)`           Time-independent problem, no exact part
  ==== ====================== ====================== =========================================================================

  Where :math:`t_0` is the time instant where the :ref:`two pictures <Exact>` coincide.

.. class:: structure::Hamiltonian

  ``template <int RANK, TimeDependence TD>`` (cf. :ref:`template parameters <structureTemplates>`);

  .. function:: void addContribution(double t, const StateVectorLow& psi, StateVectorLow& dpsidt, double tIntPic0)

    Pure virtual.

    Adds the Hamiltonian contribution :math:`\frac{H(t)}i\ket\Psi` of the given (sub)system to ``dpsidt`` assuming that the time when the Schrödinger picture and interaction picture (if any) coincide is ``tIntPic0``. There are two important points to note:

    1. The contribution has to be *added* to ``dpsidt`` instead of ``dpsidt`` being *replaced*. This is because when the given system is embedded in a larger system, other (sub)systems may also contribute.

      .. todo::

        Somehow signal to the function whether it has to add or replace.

    2. The function has to calculate the effect of :math:`\frac{H(t)}i` and not merely :math:`H`, since it is the former which determines the derivative of the state vector. This is so often missed, that we emphasize it again (although we know that it will still be missed from time to time):

      .. warning::

        When implementing the Hamiltonian, not :math:`H` itself but :math:`\frac Hi` has to supplied!

    For the class, we use partial specializations in the template parameter ``TimeDependent``, and both cases ``ONE_TIME`` and ``NO_TIME`` inherit from ``TWO_TIME``.

    ``ONE_TIME``
      defines ::

        virtual void addContribution(double, const StateVectorLow&, StateVectorLow&) const = 0;

      and implements the :func:`above function <structure::Hamiltonian::addContribution>` as ::

        void addContribution(double t, const StateVectorLow& psi, StateVectorLow& dpsidt, double tIntPic0) const
        {
          addContribution(t-tIntPic0,psi,dpsidt);
        }

    ``NO_TIME``
      defines ::

        virtual void addContribution(const StateVectorLow&, StateVectorLow&) const = 0;

      and implements the :func:`above function <structure::Hamiltonian::addContribution>` as ::

        void addContribution(double, const StateVectorLow& psi, StateVectorLow& dpsidt, double) const
        {
          addContribution(psi,dpsidt);
        }



.. class:: structure::TridiagonalHamiltonian

  ``template <int RANK, bool IS_TD>`` (cf. :ref:`template parameters <structureTemplates>`); inherits publicly from :class:`~structure::Hamiltonian`\ ``<RANK,ONE_TIME>`` when ``IS_TD=true`` and :class:`~structure::Hamiltonian`\ ``<RANK,NO_TIME>`` when ``IS_TD=false``. The present architecture of the classes :class:`~quantumoperator::Tridiagonal` and :class:`~quantumoperator::Frequencies` does not allow to cover the case ``TWO_TIME``.

  It implements the action of a Hamiltonian

  .. math::

    H_\text{tridiagonals}(t)=H_0(t)+H_1(t)+H_2(t)+\dots

  With the :math:`H_i(t)` being all described by :class:`~quantumoperator::Tridiagonal`\ ``<RANK>`` objects, their time-dependence being described by corresponding :class:`~quantumoperator::Frequencies`\ ``<RANK>`` objects.

  Such a class can be constructed with either a list of :class:`~quantumoperator::Tridiagonal`\ ``<RANK>`` (and, in the case when ``IS_TD`` is ``true``, corresponding :class:`~quantumoperator::Frequencies`\ ``<RANK>``) objects, or only one such object when the above sum consists of only one term.


.. _Exact:

---------------------------------------
``Exact``
---------------------------------------

Experience shows that even when a system uses interaction picture (which is automatically the case if any of its subsystems does)---that is, part of its dynamics is solved exactly---it may still want to calculate the jump operators and quantum averages in the normal picture. (Cases 1 & 3 :type:`above <structure::TimeDependence>`.) This is useful e.g. to reuse the code written for the non-interaction-picture case.

In this case, the framework has to be provided with some means to transform between the two pictures. This is fulfilled by the following class, from which classes describing such systems have to inherit.

E.g. if :class:`~quantumtrajectory::MCWF_Trajectory` sees that the simulated system inherits from :class:`~structure::Exact`, then it will make the coherent part of the evolution in interaction picture, and then transform back to normal picture, so that all the rest (jump probabilities, eventual jumps, calculation of quantum averages) can take place in this latter picture. This makes that the two pictures coincide before each timestep. (Cf. also the stages described in Sec. :ref:`MCWF_Trajectory`.)

.. class:: structure::Exact

  ``template <int RANK>`` (cf. :ref:`template parameters <structureTemplates>`).

  .. function:: void actWithU(double dt, StateVectorLow& psi) const

    Pure virtual.

    Describes the operation which transforms from interaction picture to the normal picture.

  .. function:: bool isUnitary() const

    Pure virtual.

    Describes whether the above operation is unitary.

.. class:: structure::FreeExact

  .. function:: void updateU(double dtDid) const


---------------------------------------
``Liouvillean``
---------------------------------------

.. class:: structure::Liouvillean

  ``template <int RANK, bool IS_TD>`` (cf. :ref:`template parameters <structureTemplates>`), we use partial specialization in the second parameter, as with :class:`~structure::Hamiltonian`:

  ``IS_TD=true``
    corresponds to Cases 1 & 2 :type:`above <structure::TimeDependence>`.

  ``IS_TD=false``
    corresponds to Cases 3 & 4 :type:`above <structure::TimeDependence>`. This partial specialization inherits from :class:`~structure::Liouvillean`\ ``<RANK,true>`` and similar function forwarding happens as in :class:`~structure::Hamiltonian`.

  .. type:: Probabilities

    This is just a ``TTD_DArray<1>``.

  .. function:: const Probabilities probabilities(double t, const LazyDensityOperator& matrix) const

    Pure virtual. The first argument is missing when ``IS_TD=false``.

  .. function:: void actWithJ(double t, StateVectorLow& psi, size_t jumpNo) const

    Pure virtual. The first argument is missing when ``IS_TD=false``. The last argument describes which jump is to be performed (a system may define several, assigning an ordinal to each of them).

.. class:: structure::ElementLiouvillean

---------------------------------------
``Averaged``
---------------------------------------

.. class:: structure::Averaged

  ``template <int RANK, bool IS_TD>`` (cf. :ref:`template parameters <structureTemplates>`), the same holds for the time-dependence as with :class:`~structure::Liouvillean`.

  .. type:: Averages

    This is just a ``TTD_DArray<1>``.

  .. function:: const Averages average(double t, const LazyDensityOperator& matrix) const

.. class:: structure::ElementAveraged




-----------------------------------
``MatrixOfHamiltonian``
-----------------------------------

.. note:: 

   This works correctly only in Schrödinger picture!


.. rubric:: Footnotes

.. [#f1] Note that this figure is tentative, and does not fully reflect the actual situation, which cannot be displayed in such a way either, due to the heavy use of templates, partial specializations, and conditional inheritance.
