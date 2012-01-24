
*************
Introduction
*************

**This manual describes C++QED: a framework for simulating open quantum dynamics.**

It is mainly for developers who wish to extend the framework for their own further purposes. For a first introduction and the description of some basic ideas underlying the framework cf. the `tutorial <http://cppqed.sourceforge.net/tutorial/tutorial.html>`_.

The organization of the framework may be sketched as follows:

.. image:: figures/organization.png
   :height: 676

The document is accordingly organized as follows:

* In section :ref:`multiarray` we describe the basic organizing concept of the framework, which directly follows from the algebra of composite quantum systems, and entails the need for using C++ template metaprogramming (TMP). This section would properly belong to a latter section describing the utility library :ref:`utils <cpputils>`, but due to its importance we put it to the beginning.

* The following sections :ref:`quantumdata`, :ref:`structure`, and :ref:`quantumtrajectory` describe the core of the framework.

* The :ref:`quantumoperator` section describes helper classes which represent quantum operators.

* The :ref:`composites` section describes how composite quantum system are defined and manipulated in the framework.

* There follows a selection of :ref:`generalElements` for demonstrating how elementary quantum systems acting as building blocks for composites can be implemented.

* The next section describes our utility library :ref:`utils <cpputils>`. The scope of this extends far beyond C++QED proper, defining very general (physical) concepts, which may be useful for other projects as well.

* Thorough :ref:`Testing <testing>` of the framework is a highly nontrivial problem, which forms the subject of the next section.

* Finally, the :ref:`appendices` describe some physical considerations underlying the framework, followed by :ref:`rationales` of style, design, and implementation, which have been observed throughout the development, and should be observed by any extension as well.


================
A note on ranges
================

For algorithms we rely throughout on the range concept from the `Boost.Range <http://www.boost.org/doc/libs/1_44_0/libs/range/index.html>`_ library. Since at the moment (Boost version 1.41) the range algorithms are not part of Boost yet (this will be Boost.RangeEx in the future), these are provided in the distribution of the framework under ``C++Utils/include/range_ex/``.


=============================
A note on the use of Blitz++
=============================

A general problem is that the use of ``int``, ``size_t``, and ``ptrdiff_t`` is not consistent. In the framework we tried to use them consistently, but in Blitz only ``int`` is used in all situations like indexing, extents and even rank template parameters, so in the interaction with Blitz we could not remain consistent.

We have used the main trunk of Blitz throughout, but later it might be desirable to switch to the 64-bit trunk.

Our own extensions to blitz can be found in :ref:`C++Utils <cpputils>` and are defined in :: 

  namespace blitzplusplus


.. _globalDefs:

===========================
Global typedefs and macros
===========================

The following definitions are in effect throughout the framework and the present manual:

.. type:: dcomp 

  ::

    typedef std::complex<double> dcomp;

.. c:var:: DCOMP_I

  ::

    const dcomp DCOMP_I(0,1);


-----------------
Template aliases
-----------------

In this manual, we are relying on the C++0x feature of *template aliases* to make documentation easier. However, in the framework, we are not using this feature, as this would limit too much the range of compilers able to cope with the framework.

For instance:

.. class:: TTD_DArray

  ``template <int RANK>``

  is assumed in this document to be a template alias::

    template <int RANK> using TTD_DArray=blitz::Array<double,RANK>;

  In the framework, however, it is simply defined as a macro (with all capitals)::

    #define TTD_DARRAY(r) blitz::Array<double,r>

.. class:: TTD_CArray

  ``template <int RANK>``::

    template <int RANK> using TTD_CArray=blitz::Array<dcomp,RANK>;

.. class:: TTD_ExtTiny

  ``template <int RANK>``::

    template <int RANK> using TTD_ExtTiny=blitz::TinyVector<size_t,RANK>;

.. class:: TTD_IdxTiny::

  ``template <int RANK>``::

    template <int RANK> using TTD_IdxTiny=blitz::TinyVector<ptrdiff_t,RANK>

.. type:: linalg::CVector

  ::

    namespace linalg {

    typedef TTD_CArray<1> CVector;

    } // linalg

.. type:: linalg::CMatrix

  ::

    namespace linalg {

    typedef TTD_CArray<2> CMatrix;

    } // linalg

The prefix ``TTD`` stands for "template typedef" here and throughout the framework.


-------------------
Metafunctions
-------------------

Metafunctions are class templates used to "compute" some type. By convention, the resulting type is stored as a ``type`` typename.

In the framework, metafunctions are named with an ``MF`` suffix. Many of them are simple template aliases, but with some extra functionality, as e.g. compile-time assertions.

------------------
Namespace aliases
------------------

We will also use the following namespace alias ::

  namespace mpl=boost::mpl;

However, in the framework we are not using this alias globally, as this would lead to all sorts of name clashes.

We assume that the following definitions are in effect:

.. c:var:: mpl::constant_true

  ::

    const mpl::true_ mpl::constants_true;

.. c:var:: mpl::constant_false

  ::

    const mpl::false_ mpl::constant_false;
