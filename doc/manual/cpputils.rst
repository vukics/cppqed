.. _cpputils:

=====================================================
Some further general-purpose modules from C++Utils
=====================================================

From C++Utils what deserves documentation.

---------------------------------------------
``Evolved`` and the underlying principle
---------------------------------------------


------------------------------
The ``trajectory`` namespace
------------------------------

----------
``Pars``
----------

.. todo::

   Pars should first read the default values from a file, and update them from the command line afterwards. (Pars should anyway be replaced by the boost thingy)

-----------------------
``MultiIndexIterator``
-----------------------

.. class:: cpputils::MultiIndexIterator


--------------------
``vfmsi::Iterator``
--------------------



------------------
Other
------------------

.. todo::

   Implement a class representing a non-orthogonal vector. It should store the pulled vector, and it stould keep track of whether it is up to date. Eg any change in any element makes the pulled vector out of date, and it has to be brought up to date for any operation involving the metric. The same for matrices and indeed tensors of any order. (Also, could make normal tensors of any order, maybe out of CVector using the boost thingy?)

.. todo::

   HermiteCoefficients --- do it better (maybe even with TMP?).
