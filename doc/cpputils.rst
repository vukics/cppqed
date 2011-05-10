.. _cpputils:

=====================================================
Some further general-purpose modules from ``utils``
=====================================================

From ``utils`` what deserves documentation.

.. _cpputils_Evolved:

---------------------------------------------
``Evolved`` and the underlying principle
---------------------------------------------


------------------------------
The ``trajectory`` namespace
------------------------------

.. function:: void trajectory::runDt(trajectory::TrajectoryBase&, double time, double deltaT, bool displayInfo)

.. function:: void trajectory::run  (trajectory::Trajectory<A> &, double time, int          , bool displayInfo)

  template<typename A>


.. _cpputils_Parameters:

----------
``Pars``
----------

.. todo::

   Pars should first read the default values from a file, and update them from the command line afterwards. (Pars should anyway be replaced by the boost thingy)

.. class:: parameters::ParameterTable

.. function:: void parameters::update(parameters::ParameterTable& table, int argc, char** argv, const std::string& mod="--")

-----------------------
``MultiIndexIterator``
-----------------------

.. class:: cpputils::MultiIndexIterator


.. _cpputils_VFMSI:

------------------------------
VectorFromMatrixSliceIterator
------------------------------


------------------
Other
------------------


.. class:: tmptools::OrdinalMF

  ``template <int RANK>``

  .. type:: type


.. class:: tmptools::numerical_contains

  ``template<typename Seq, typename ICW>``


.. class:: tmptools::Vector

  ``template <int V0=TMPTOOLS_VECTOR_DEFAULT_ARG, ...>``


.. class:: linalg::VectorSpace

  ``template<typename T, typename B>``



.. todo::

   Implement a class representing a non-orthogonal vector. It should store the pulled vector, and it should keep track of whether it is up to date. Eg any change in any element makes the pulled vector out of date, and it has to be brought up to date for any operation involving the metric. The same for matrices and indeed tensors of any order. (Also, could make normal tensors of any order, maybe out of CVector using the boost thingy?)

.. todo::

   HermiteCoefficients --- do it better (maybe even with TMP?).
