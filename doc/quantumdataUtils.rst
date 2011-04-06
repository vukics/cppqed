``DimensionsBookkeeper``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. py:module:: DimensionsBookkeeper.h
   :synopsis: Defines DimensionsBookkeeper

.. class:: DimensionsBookkeeper

  :ref:`template parameters <quantumdataTemplates>`: RANK, IS_CONST

  This class is designed to store and manipulate dimensions of constructs on composite Hilbert spaces of arity ``RANK``.

  It can be either constant or non-constant depending on the second template parameter.

  It has a straightforward interface which does not need documentation, the only notable point being a technique which is used from time to time throughout the framework. Consider the following constructor:

  .. function:: explicit DimensionsBookkeeper(mpl::bool_<IS_CONST> = mpl::false_ __LP__ __RP__ )

    The aim of the dummy argument with a default value---which creates a nonsensical function signature in the case when ``IS_CONST`` is ``true``---is that this constructor only compiles in the case when ``IS_CONST`` is ``false`` because it is only in the non-constant case that we allow default construction of the class. Since from a template only such parts are compiled as are actually used, a client can use the class in the case when ``IS_CONST`` is ``true`` without problems, getting a compile-time error only when trying to default-construct such an object.

  .. type:: Dimensions

    ::

      typedef TTD_ExtTiny<RANK> Dimensions;



``ArrayBase``
^^^^^^^^^^^^^^^^^^^

.. py:module:: ArrayBase.h
   :synopsis: Defines ArrayBase in namespace quantumdata

.. class:: quantumdata::ArrayBase

  :ref:`template parameters <quantumdataTemplates>`: RANK, IS_CONST


``Types``
^^^^^^^^^^^^^

.. py:module:: Types.h
   :synopsis: Defines the metafunction Types in namespace quantumdata

.. class:: quantumdata::Types

  :ref:`template parameters <quantumdataTemplates>`: RANK, B

  This class is basically only a metafunction.

  .. type:: StateVectorLow

  ::

    typedef TTD_CArray<RANK> StateVectorLow;

  .. type:: DensityOperatorLow

  ::

    typedef TTD_CArray<2*RANK> DensityOperatorLow;

