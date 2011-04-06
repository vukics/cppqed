Defining metrical transformations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. py:module:: Transformation.h
   :synopsis: Provides the infrastructure for representing metrical transformations in non-orthogonal basesoperations

Any type can represent a metrical transformation, that is, play the role of the template parameter ``TRAFO`` in the non-orthogonal classes, as long as a class :class:`quantumdata::transformation::Traits` can be specialized for the type.

The general requirements on the structure of the traits class is:

.. class:: quantumdata::transformation::Traits

  :ref:`template parameters <quantumdataTemplates>`: TRAFO

  .. c:var:: N_RANK

    ::

      static const int N_RANK

  .. type:: TrafoTypes

    The :class:`~quantumdata::transformation::Composite` transformation stores the underlying elementary transformations as a `Boost.Fusion <http://www.boost.org/doc/libs/1_44_0/libs/fusion/doc/html/fusion/container/vector.html>`_ vector, and this type has to represent such a vector.

  .. type:: StateVectorLow

  .. function:: static void transform(const TRAFO& trafo, const StateVectorLow& in, StateVectorLow& out)

  .. function:: static const TrafoTypes& getTrafos(const TRAFO& trafo)


There is a class which helps with the specialization of :class:`quantumdata::transformation::Traits` for elementary transformations. Traits specializations for elementary transformations can simply inherit from this class:

.. class:: quantumdata::transformation::ElementaryTraits

  :ref:`template parameters <quantumdataTemplates>`: TRAFO

  ::
  
    template<typename TRAFO>
    struct ElementaryTraits
    {
      typedef typename boost::fusion::vector<TRAFO> TrafoTypes;
      // Simply creates a single-element vector out of TRAFO, so that afterwards it can be treated in a unified way with Composites

      static const TrafoTypes getTrafos(const TRAFO& trafo) {return boost::fusion::make_vector(trafo);}
      // The same at runtime

    };


Traits specializations
%%%%%%%%%%%%%%%%%%%%%%%%

For the dummy transformation :class:`~quantumdata::transformation::Identity`:

  ::

    template<int RANK>
    struct Traits<Identity<RANK> > : ElementaryTraits<Identity<RANK> >
    {
      static const int N_RANK=RANK;

      typedef typename Types<N_RANK>::StateVectorLow StateVectorLow;
      static void transform(Identity<RANK> trafo, const StateVectorLow& in, StateVectorLow& out) {out=in;}

    };

  (Most of the times, the dummy ``transform`` function will not be called for such classes because the :class:`~quantumdata::transformation::Identity` components are filtered out at compile time from :class:`~quantumdata::transformation::Composite` transformations.)

A :class:`TTD_CArray`\ ``<2*RANK>`` can represent a transformation of rank ``RANK``:

  ::

    template<int TWO_TIMES_RANK>
    struct Traits<TTD_CArray<TWO_TIMES_RANK> > : ElementaryTraits<TTD_CArray<TWO_TIMES_RANK> >
    {
      static const int N_RANK=TWO_TIMES_RANK/2;

      typedef typename Types<N_RANK>::StateVectorLow StateVectorLow;
      static void transform(const TTD_CArray<TWO_TIMES_RANK>& trafo, const StateVectorLow& in, StateVectorLow& out);

    };

  In this case the ``transform`` function is implemented to represent the application of the matrix.

A pointer to a function with the appropriate signature:

  ::

    template<int RANK>
    struct Traits< void(*)(const TTD_CArray<RANK>&, TTD_CArray<RANK>&) > : ElementaryTraits< void(*)(const TTD_CArray<RANK>&, TTD_CArray<RANK>&) >
    {
      static const int N_RANK=RANK;

      typedef typename Types<N_RANK>::StateVectorLow StateVectorLow;

      typedef void(*TRAFO)(const StateVectorLow&, StateVectorLow&);

      static void transform(const TRAFO trafo, const StateVectorLow& in, StateVectorLow& out) {trafo(in,out);}

    };

