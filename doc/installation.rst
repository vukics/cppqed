************
Requirements
************

C++QEDv2 depends on a number of open source libraries:

Boost C++ libraries 
  provide indispensable extensions to the C++ standard, and are *de facto* standard by their own right. The framework depends on a number of them, the most notable ones being Fusion, Lambda, MPL, Operators, and Preprocessor. On many systems, (a selection of) the Boost libraries are available. They are packaged for Debian and Mac OS X. Alternatively, they can be downloaded and installed from `the main Boost portal <http://www.boost.org>`_. Version 1.35 or higher is required.

GNU Scientific library (GSL) 
  provides a very wide variety of numerical solutions in a solid object-oriented design (in C!). They are not used directly, but are wrapped into C++ classes and functions, so that they are easily replaced (e.g. if licensing problems arise). Packaged for Debian and Mac OS X, or can be downloaded from `the GSL homepage <http://www.gnu.org/software/gsl/>`_. I haven't thoroughly determined the minimal version, but 1.8 is known to work.

An implementation of BLAS and LAPACK
  For many system, optimized versions of these libraries exist and are even preinstalled. As a fallback, the reference implementation can always be used, which is packaged e.g. for Ubuntu.

These three are best installed on system level. 

The following two libraries are stable, but under more or less steady development.

.. highlight:: sh


Blitz++ 
  provides the fundamental data structure. It hence performs a lot of numerics and lies at the absolute heart of the framework. Blitz++ lives up to its name, as it provides near-Fortran performance in spite of the very high-level abstractions used in the library. This is achieved by TMP, which was discovered in prototype during the development of this very library. More on `the Blitz++ homepage <http://www.oonumerics.org/blitz/>`_.

  .. warning::

    At the time of this writing the released version of Blitz++ is version 0.9. Do not use this release! The CVS version is to be used instead. There is no knowing when a new release will appear, but that future release can be used, of course.

Flexible Library for Efficient Numerical Solutions (FLENS) 
  is a very impressive effort to wrap BLAS-LAPACK functions in a high-level C++ interface. FLENS in turn depends on BLAS-LAPACK (ATLAS). More on `the FLENS homepage <http://flens.sourceforge.net/>`_.

  .. note::

    The use of FLENS is optional because only a very small section of the framework depends on it. Cf. the section :ref:`compilingWithoutFLENS`.

At the corresponding websites, instructions for installing the libraries can be found, while below I have given guidelines concerning installation on an Ubuntu system, cf. also the ``getLibs.sh`` script in directory ``utils``. This will download the CVS versions of both libraries into the subdirectories ``blitz`` and ``FLENS-lite``, respectively, and will also compile them.

I will very much appreciate all feedback regarding also the installation of the framework. 

======================
Installation on Ubuntu
======================

The following is a rudimentary list of prerequisite packages on Ubuntu 9.10:: 

  % sudo apt-get install boost-build libboost1.40-all-dev autoconf automake libtool libgsl0-dev liblapack-dev liblapack-pic liblapack3gf libblas-dev libblas3gf libatlas-headers libatlas3gf-base rcs cvs

This will install the first three libraries listed above on system level.

The following steps to install the remaining two libraries constitute also the script ``utils/getLibs.sh``.

-----------------
Install Blitz++:
-----------------

::

  cvs -d:pserver:anonymous@blitz.cvs.sourceforge.net:/cvsroot/blitz login

Hit enter when prompted for password. ::

  cvs -z3 -d:pserver:anonymous@blitz.cvs.sourceforge.net:/cvsroot/blitz co -P blitz
  cd blitz
  autoreconf -fiv
  ./configure --with-pic
  make lib
  sudo make install

---------------
Install FLENS:
---------------

.. note::

  The use of FLENS is optional because only a very small section of the framework depends on it. Cf. the section :ref:`compilingWithoutFLENS`.

::

  cvs -d:pserver:anonymous@flens.cvs.sourceforge.net:/cvsroot/flens login
  cvs -z3 -d:pserver:anonymous@flens.cvs.sourceforge.net:/cvsroot/flens co -P FLENS-lite
  cd FLENS-lite
  cp config.ubuntu config

Now you have to edit the config file adding to ``CXXFLAGS`` the flag ``-DGSL_CBLAS`` which instructs FLENS to use the CBLAS interface provided by GSL. This is good because hence you don't need a separate package for this. You may also need to remove the flag ``-latlas`` from ``LDFLAGS``. ::

  make
  sudo make install

If the last command issues the error message::

  Makefile.common:19: /config: No such file or directory

you have to edit ``Makefile.common`` replacing the variable ``$(PWD)`` with the path of the current directory.


================
Obtaining C++QED
================

There are two ways, the first being to download the latest package from `<http://sourceforge.net/projects/cppqed/files/>`_. This is only recommended if the package is not too old.

The other is to use the `Bazaar <https://sourceforge.net/scm/?type=bzr&group_id=187775>`_ version::

  bzr checkout bzr://cppqed.bzr.sourceforge.net/bzrroot/cppqed C++QED

Where the last argument can be replaced to the name of the directory for the code to appear in. Alternately, an existing checkout can be updated as::

  bzr pull bzr://cppqed.bzr.sourceforge.net/bzrroot/cppqed

Be aware that C++QED is under development, so changes in the Bazaar version may change the API of certain modules in such a way as breaks your applications. It is advisable to follow the `ChangeLog <http://cppqed.sourceforge.net/changelog.html>`_ of the project. Alternately, the CVS option ``-r date:<date>`` can be used to retrieve the most recent revision no later than ``<date>``. E.g.::

  bzr pull -r date:2010-02-14 bzr://cppqed.bzr.sourceforge.net/bzrroot/cppqed


===========
Compilation
===========

The canonical way to compile the framework is the one using Boost.Build. This is best installed on system level. Typing ::

  bjam 

in the main directory will compile and link the whole framework, creating separate executables from the highest level programs residing in directory ``scripts``. Typing ::

  bjam <script-name-without-extension>

will compile only the given script.


The default compilation mode is ``debug``\ ging mode, meaning that in this case a lot of runtime checks are compiled into the framework, which come from Blitz++, FLENS, and myself. Every time a new script is added it should be compiled and tested in this way because this can detect a *lot* of errors. When we are absolutely sure that everything is all right, for data collection we may compile with ``bjam release``, in which all the checks are omitted and optimisations are used, making the programs *about an order of magnitude faster*.

.. warning::

   Maximum efficiency is achieved only if the framework is compiled with ::

     bjam release 

   or ::

     bjam <script-name-without-extension> release

``bjam`` will put the compiled files into the directories ``bin`` and ``utils/bin``. These directories are the roots of directory structures which mirror the structure of the distribution.

There is a ``Makefile`` which will automatically recognise the executables in directory ``scripts``, compile the framework, and statically link it with necessary libraries. Although with ``make`` it is not easy to provide the same flexibility as with Boost.Build, I am trying to maintain this possibility on an acceptable level. With ``make``, the default compilation mode is optimized mode. Type ::

  make utils
  make <script-name-without-extension>

To switch to debugging mode you need to use ::

  make optimization=no <script-name-without-extension>

Boost.Build, just like ``make``, supports parallel compilation, which can make a significant difference for projects of the magnitude of C++QEDv2. For starting ``n`` threads of compilation use ::

  bjam -j n ... 

C++QEDv2 has been successfully compiled on several Linux platforms and Mac OS X. In all cases the GNU C++ Compiler has been used. It also compiles with the `clang++ <http://clang.llvm.org/>`_ compiler. Portability to other compilers remains to be demonstrated.


.. todo::

   In the make version of build, clarify dynamic linkage.

.. todo::

   Test framework with icc + under Windows

.. _compilingWithoutFLENS:

-------------------------
Compiling without FLENS
-------------------------

There is a compilation feature which can be supplied to Boost.Build::

  bjam with-flens=no <all the rest as before>

In this case, those parts of the framework that rely on FLENS are discreetly disabled. Most notable is the calculation of the negativity of partially transposed density operators, cf. :ref:`assessingEntanglement`.

.. highlight:: c++
  :linenothreshold: 10


=========
``utils``
=========

The content of the directory ``utils`` is a small library of very diverse but quite general tools, that I have abstracted during the development of the framework, and used also in several other projects. This may in time become a project on its own. The reader is encouraged to have a look in there, too: some modules may be useful in themselves. Cf. :ref:`cpputils`.
