
.. |cppqed| replace:: C++QED

..
  image:: figures/logoCppQED.svg
  :alt: C++QED
  :width: 20

********************************************************
C++QED: a framework for simulating open quantum dynamics
********************************************************

**C++QED is a highly flexible framework for simulating open quantum dynamics. It allows users to build arbitrarily complex interacting quantum systems from elementary free subsystems and interactions, and simulate their time evolution with a number of available time-evolution drivers.**

=============
Documentation
=============

To start to understand and use the framework on the highest level, consult the :doc:`User Guide <tutorial>`. This will allow for building simulations using the already existing modules and elementary systems as building blocks.

.. toctree::
   :maxdepth: 2

   User Guide <tutorial>

To be able to write new modules representing elementary free subsystems or interactions, consult the :doc:`Element Guide <elementGuide>`.

.. toctree::
   :maxdepth: 2

   elementGuide

For a deeper understanding and to be able to extend the framework with new modules, consult the constantly developing reference manual. If your problem is not (yet) covered therein, feel free to `contact the developers <http://sourceforge.net/project/memberlist.php?group_id=187775>`_.

.. toctree::
   :maxdepth: 2

   Reference Manual <manual>

-------------
Change Log
-------------

Users of the development version from the Bazaar repository, should consult the :doc:`Change Log <changelog>` from time to time.


-------------
PyCppQED
-------------

`PyCppQED <http://github.com/bastikr/pycppqed>`_ is a Python backend for processing and visualizing data produced by the framework:

* `PyCppQED documentation <http://cppqed.sourceforge.net/pyCppQED_doc/index.html>`_


------------------
Printable versions
------------------

* :download:`User Guide in pdf <_build/latex/C++QED_Tutorial.pdf>`

* :download:`Element Guide in pdf <_build/latex/C++QED_structureTutorial.pdf>`

* `PyCppQED user guide in pdf <http://github.com/downloads/bastikr/pycppqed/PyCppQED-0.1.1.pdf>`_



========
Download
========

.. highlight:: sh

* `The project summary page <http://sourceforge.net/projects/cppqed/>`_

* `Released packages <http://sourceforge.net/projects/cppqed/files/>`_

* The development version from the Bazaar repository::

    bzr checkout bzr://cppqed.bzr.sourceforge.net/bzrroot/cppqed C++QED


=======
Support
=======

* `The project mailing list <https://lists.sourceforge.net/lists/listinfo/cppqed-support>`_ or ...

* `... contact the developers directly <http://sourceforge.net/project/memberlist.php?group_id=187775>`_

Issues with C++QED proper should be addressed to `András Vukics <http://sourceforge.net/users/vukics>`_, those with PyCppQED to `Sebastian Krämer <http://sourceforge.net/users/bastikr>`_.


=====
Links
=====

* `Homepage of A. Vukics at the Researh Institute for Solid State Physics & Optics <http://optics.szfki.kfki.hu/Vukics/Vukics>`_

* Sebastian Krämer's :download:`Master's Thesis <_static/SebastianKraemerThesis.pdf>` at the University of Innsbruck, on simulating harmonic-oscillator modes by an adaptive set of coherent states as basis.

* `Boost C++ libraries <http://www.boost.org>`_

* `GNU Scientific library (GSL) <http://www.gnu.org/software/gsl/>`_

* `Blitz++ <http://www.oonumerics.org/blitz/>`_, `new homepage <http://blitz.sourceforge.net>`_

* `Flexible Library for Efficient Numerical Solutions (FLENS) <http://flens.sourceforge.net/>`_

* `The first journal article about the framework <http://www.springerlink.com/content/r2237020726t0614/>`_


================
Bazaar access
================

To access the Bazaar repository, configure your Bazaar client as follows::

  bzr://cppqed.bzr.sourceforge.net/bzrroot/cppqed

(read-only) ::

  bzr+ssh://vukics@cppqed.bzr.sourceforge.net/bzrroot/cppqed

(read/write—replace ``vukics`` with your username)


.. highlight:: c++
  :linenothreshold: 10

