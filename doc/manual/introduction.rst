
=============
Introduction
=============

This manual describes C++QED: a framework for simulating open quantum dynamics. 

It is mainly for developers who wish to extend the framework for their own further purposes. For a first introduction and the description of some basic ideas underlying the framework cf. the `tutorial <http://cppqed.sourceforge.net/tutorial/tutorial.html>`_.

The organization of the framework may be sketched as follows:

.. image:: figures/organization.png
   :height: 676

The document is accordingly organized as follows:

* In section :ref:`multiarray` we describe the basic organizing concept of the framework, which directly follows from the algebra of composite quantum systems, and entails the need for using template metaprogramming. This section would properly belong to a latter section describing the utility library :ref:`C++Utils <cpputils>`, but due to its importance we put it to the beginning.

* The following sections :ref:`quantumdata`, :ref:`structure`, and :ref:`quantumtrajectory` describe the core of the framework.

* The :ref:`quantumoperator` section describes helper classes which represent quantum operators.

* The :ref:`composites` section describes how composite quantum system are defined and manipulated in the framework.

* There follows a selection of :ref:`generalElements` for demonstrating how elementary quantum systems acting as building blocks for composites can be implemented.

* The next section describes our utility library :ref:`C++Utils <cpputils>`. The scope of this extends far beyond C++QED proper, defining very general (physical) concepts, which may be useful for other projects as well.

* Thorough :ref:`Testing <testing>` of the framework is a highly nontrivial problem, which forms the subject of the next section.

* Finally, the :ref:`appendices` describe some physical considerations underlying the framework, followed by :ref:`rationales` of style, design, and implementation, which have been observed throughout the development, and should be observed by any extension as well.


--------------------
A note on ranges
--------------------

For algorithms we rely throughout on the range concept from the `Boost.Range <http://www.boost.org/doc/libs/1_44_0/libs/range/index.html>`_ library. Since at the moment (Boost version 1.41) the range algorithms are not part of Boost yet (this will be Boost.RangeEx in the future), these are provided in the distribution of the framework under ``C++Utils/include/range_ex/``.


------------------------------
A note on the use of Blitz++
------------------------------

A general problem is that the use of ``int``, ``size_t``, and ``ptrdiff_t`` is not consistent. In the framework we tried to use them consistently, but in Blitz only ``int`` is used in all situations like indexing, extents and even rank template parameters, so in the interaction with Blitz we could not remain consistent.

We have used the main trunk of Blitz throughout, but later it might be desirable to switch to the 64-bit trunk.

Our own extensions to blitz can be found in :ref:`C++Utils <cpputils>` and are defined in :: 

  namespace blitzplusplus

.. todo:: 

   Try to make an array class which acts and feels like a blitz::Array with respect to the functionality needed here. It could be much simpler, leaving every work to the underlying blitz::Array storage, but solving this silly problem of size_t, ptrdiff_t, int conversions. This is not very easy at the moment because of the immense number of constructors a blitz::Array provides, but with the constructor inheritance of C++0x it should be easy.



.. _globalDefs:

-------------------------------
Global typedefs and macros
-------------------------------

The following definitions are in effect throughout the framework and
the present manual:

.. type:: dcomp 

  ::

    typedef std::complex<double> dcomp;

.. c:var:: DCOMP_I

  ::

    const dcomp DCOMP_I(0,1);

.. c:macro:: TTD_DARRAY(r)

  ::

    #define TTD_DARRAY(r) blitz::Array<double,r>

.. c:macro:: TTD_CARRAY(r)

  ::

    #define TTD_CARRAY(r) blitz::Array<dcomp ,r>

.. c:macro:: TTD_EXTTINY(r)

  ::

    #define TTD_EXTTINY(r) blitz::TinyVector<   size_t,r>

.. c:macro:: TTD_IDXTINY(r)

  ::

    #define TTD_IDXTINY(r) blitz::TinyVector<ptrdiff_t,r>

.. type:: linalg::CVector

  ::

    namespace linalg {

    typedef TTD_CARRAY(1) CVector;

    } // linalg

.. type:: linalg::CMatrix

  ::

    namespace linalg {

    typedef TTD_CARRAY(2) CMatrix;

    } // linalg

The prefix ``TTD`` stands for "template typedef" here and throughout the framework.

We will also use the following namespace alias ::

  namespace mpl=boost::mpl;

However, in the framework we are not using this alias globally, as this would lead to all sorts of name clashes.


