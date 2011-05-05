.. _structure:

=============================
The ``structure`` namespace
=============================

.. namespace:: structure

::

  namespace structure

The ``structure`` namespace comprises modules for describing quantum systems. Among them the most important is :class:`~structure::QuantumSystem`, which is an abstract interface every system has to provide to be usable with :ref:`quantum trajectories <quantumtrajectory>` like :class:`~quantumtrajectory::MCWF_Trajectory` or :class:`~quantumtrajectory::Master`. This is why all the :ref:`elementary <generalElements>` and :ref:`composite systems <composites>` are more or less directly derived from :class:`~structure::QuantumSystem`.

Most of the classes in this namespace belong to a single hierarchy, sketched in the diagram below (the broken lines signifying that the inheritance is not direct, due to some classes in between, which can be considered implementation details):

.. image:: figures/structure.png
   :height: 605

We have also indicated how :ref:`elements <generalElements>` like :class:`Mode` and :class:`JaynesCummings`, and :ref:`composite systems <composites>` as :class:`BinarySystem` and :class:`Composite` fit into the hierarchy.

This class hierarchy provides a lot of services for implementing new elements in the framework. For examples on how to optimally use these services, cf. :ref:`the structure-bundle tutorial below <structureTutorial>`.

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

  :ref:`template parameters <structureTemplates>`: RANK; inherits publicly from :class:`DimensionsBookkeeper`

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

  Inherits publicly from :class:`~structure::QuantumSystem` with RANK=1 and :class:`~structure::DynamicsBase`

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

  ::

    typedef quantumoperator::Tridiagonal<1> Tridiagonal;

.. type:: structure::free::Frequencies

  ::

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

  :ref:`template parameters <structureTemplates>`: RANK; inherits publicly from :class:`~structure::DynamicsBase`

  It does not inherit from :class:`~structure::QuantumSystem` because it does not make sense to simulate such an element as describes an interaction alone. However, an interaction can have frequency-like parameters, hence the inheritance from :class:`~structure::DynamicsBase`.

  .. type:: Frees

  ::

    typedef blitz::TinyVector<const Free*,RANK> Frees;

  .. function:: explicit Interaction(const Frees& frees, const RealFreqs& realFreqs=RealFreqs__LP__ __RP__, const ComplexFreqs& complexFreqs=ComplexFreqs__LP__ __RP__)

    

  .. function:: const Frees& getFrees() const



---------------------------------------
``Hamiltonian``
---------------------------------------

.. class:: structure::Hamiltonian

.. warning::

  When implementing the Hamiltonian, not :math:`H` itself but :math:`\frac Hi` has to supplied!


.. class:: structure::TridiagonalHamiltonian


---------------------------------------
``Exact``
---------------------------------------

.. class:: structure::Exact

.. class:: structure::FreeExact

  .. function:: void updateU(double dtDid) const


---------------------------------------
``Liouvillean``
---------------------------------------

.. note:: 

  The framework requires that the jump be calculated in the normal picture even in the case when the Hamiltonian is in interaction picture (cf. :func:`quantumtrajectory::MCWF_Trajectory::step`). This allows for reusing the same code in both pictures.

.. class:: structure::Liouvillean

.. class:: structure::ElementLiouvillean

---------------------------------------
``Averaged``
---------------------------------------

.. class:: structure::Averaged

  .. function:: const Averages average(const LazyDensityOperator& matrix) const

.. class:: structure::ElementAveraged

.. note:: 

  The framework requires that the averages be calculated in the normal picture even in the case when the Hamiltonian is in interaction picture (cf. :func:`quantumtrajectory::MCWF_Trajectory::step`). This allows for reusing the same code in both pictures.




-----------------------------------
``MatrixOfHamiltonian``
-----------------------------------

.. note:: 

   This works correctly only in Schrödinger picture!


