.. _cpputils:

=====================================================
Some further general-purpose modules from ``utils``
=====================================================

From ``utils`` what deserves documentation.

.. _cpputils_Evolved:

---------------------------------------------
``Evolved`` and the underlying principle
---------------------------------------------

.. class:: evolved::Evolved


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


.. class:: tmptools::Ordinals

  ``template <int RANK>``

  .. type:: type


.. class:: tmptools::numerical_contains

  ``template<typename Seq, typename ICW>``


.. class:: tmptools::Vector

  ``template <int V0=TMPTOOLS_VECTOR_DEFAULT_ARG, ...>``


.. class:: linalg::VectorSpace

  ``template<typename T, typename B>``

