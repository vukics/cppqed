
********************************************************
C++QED: a framework for simulating open quantum dynamics
********************************************************

**C++QED is a highly flexible framework for simulating open quantum dynamics. It allows users to build arbitrarily complex interacting quantum systems from elementary free subsystems and interactions, and simulate their time evolution with a number of available time-evolution drivers.**


=============
Documentation
=============

To start to understand and use the framework on the highest level, consult the user guide. This will allow for building simulations using the already existing modules and elementary systems as building blocks.

.. toctree::
   :maxdepth: 2

   User Guide <tutorial>

To be able to write new modules representing elementary free subsystems or interactions, consult the following guide:

.. toctree::
   :maxdepth: 2

   elementGuide

For a deeper understanding and to be able to extend the framework with new modules, the upcoming reference manual will have to be consulted. In the meantime feel free to `contact the developers <http://sourceforge.net/project/memberlist.php?group_id=187775>`_.

.. toctree::
   :maxdepth: 2

   Reference Manual <manual>

-------------
ChangeLog
-------------

.. toctree::
   :maxdepth: 2

   If you are using the development version from the Bazaar repository, you will need to consult the ChangeLog from time to time. <changelog>

------------------
Printable versions
------------------

* :download:`User Guide in pdf <_build/latex/C++QED_Tutorial.pdf>`

* :download:`Element Guide in pdf <_build/latex/C++QED_structureTutorial.pdf>`

* `PyCppQED tutorial in pdf <http://github.com/downloads/bastikr/pycppqed/PyCppQED-0.1.1.pdf>`_



========
Download
========

.. highlight:: sh

* `The project summary page <http://sourceforge.net/projects/cppqed/>`_

* `Released packages <http://sourceforge.net/projects/cppqed/files/>`_

* The development version from the Bazaar repository::

    % bzr checkout bzr://cppqed.bzr.sourceforge.net/bzrroot/cppqed C++QED


=========================================================================================
Support
=========================================================================================

* `The project mailing list <https://lists.sourceforge.net/lists/listinfo/cppqed-support>`_ or ...

* `... contact the developers directly <http://sourceforge.net/project/memberlist.php?group_id=187775>`_

Issues with C++QED proper should be addressed to `András Vukics <http://sourceforge.net/users/vukics>`_, those with PyCppQED to `Sebastian Krämer <http://sourceforge.net/users/bastikr>`_.


=====
Links
=====

* `Boost C++ libraries <http://www.boost.org>`_

* `GNU Scientific library (GSL) <http://www.gnu.org/software/gsl/>`_

* `Blitz++ <http://www.oonumerics.org/blitz/>`_, `new homepage <http://blitz.sourceforge.net>`_

* `Flexible Library for Efficient Numerical Solutions (FLENS) <http://flens.sourceforge.net/>`_

* `The first journal article about the framework <http://www.springerlink.com/content/r2237020726t0614/>`_


================
Bazaar access
================

To access the Bazaar repository, configure your Bazaar client as follows::

  bzr://cppqed.bzr.sourceforge.net/bzrroot/cppqed (read-only)

  bzr+ssh://USERNAME@cppqed.bzr.sourceforge.net/bzrroot/cppqed (read/write)


.. highlight:: c++
  :linenothreshold: 10

