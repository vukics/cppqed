************
Requirements
************

C++QEDv2 depends on a number of open source libraries:

Boost C++ libraries 
  provide indispensable extensions to the C++ standard, and are *de facto* standard by their own right. The framework depends on a number of them, the most notable ones being Fusion, Lambda, MPL, Operators, and Preprocessor. On many systems, (a selection of) the Boost libraries are available. They are packaged for Debian and Mac OS X. Alternatively, they can be downloaded and installed from `the main Boost portal <http://www.boost.org>`_. Although version 1.35 is in principle sufficient for the framework, it is advisable to use as new a version as possible, since experience has shown that there can be clashes between old versions of Boost and new versions of compilers.

  All the Boost libraries necessary for the framework are header-only. Therefore, we also provide an alternative package in which the Boost dependencies are placed in a subdirectory under ``utils/include``. When using this package, no separate Boost libraries need to be installed. Cf. sec. :ref:`boostIntegration`.

GNU Scientific library (GSL) 
  provides a very wide variety of numerical solutions in a solid object-oriented design (in C!). They are not used directly, but are wrapped into C++ classes and functions, so that they are easily replaced (e.g. if licensing problems arise). Packaged for Debian and Mac OS X, or can be downloaded from `the GSL homepage <http://www.gnu.org/software/gsl/>`_. I haven't thoroughly determined the minimal version, but 1.8 is known to work.

An implementation of BLAS and LAPACK
  For many systems, optimized versions of these libraries exist and are even preinstalled. As a fallback, the reference implementation can always be used, which is packaged e.g. for Ubuntu.

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

  % sudo apt-get install boost-build libboost1.40-all-dev autoconf automake libtool libgsl0-dev liblapack-dev liblapack-pic liblapack3gf libblas-dev libblas3gf libatlas-headers libatlas3gf-base bzr

This will install the first three libraries listed above on system level.

The following steps to install the remaining two libraries constitute also the script ``utils/getLibs.sh``. The ``sed`` command in the last lines of this script also demonstrates what needs to be changed in ``utils/Jamroot`` if Blitz and FLENS are installed at a location which is not in the search path of the C++ compiler/linker.

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

Under Debian-like operating systems, instead of the last line one can use::

  sudo checkinstall --fstrans=0 -D make install

The ``--fstrans=0`` option has to be used because of a bug in ``checkinstall``.

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

Alternatively, again::

  sudo checkinstall --fstrans=0 -D make install


================
Obtaining C++QED
================

There are two ways, the first being to download the latest package from `<http://sourceforge.net/projects/cppqed/files/>`_. This is only recommended if the package is not too old.

The other is to use the `Bazaar <https://sourceforge.net/scm/?type=bzr&group_id=187775>`_ version::

  bzr checkout bzr://cppqed.bzr.sourceforge.net/bzrroot/cppqed C++QED

Where the last argument can be replaced to the name of the directory for the code to appear in. Alternately, an existing checkout can be updated as::

  bzr pull bzr://cppqed.bzr.sourceforge.net/bzrroot/cppqed

Be aware that C++QED is under development, so changes in the Bazaar version may change the API of certain modules in such a way as breaks your applications. It is advisable to follow the `ChangeLog <http://cppqed.sourceforge.net/changelog.html>`_ of the project. Alternately, the Bazaar option ``-r date:<date>`` can be used to retrieve the most recent revision no later than ``<date>``. E.g.::

  bzr pull -r date:2010-02-14 bzr://cppqed.bzr.sourceforge.net/bzrroot/cppqed


.. _boostIntegration:

=================
Boost integration
=================

To obtain the package with the necessary Boost libraries integrated, download the package file with ``...BoostIntegration...`` in its name. To get the development version, the corresponding Bazaar branch has to be used::

  bzr checkout bzr://cppqed.bzr.sourceforge.net/bzrroot/cppqed/BoostIntegration C++QED

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

A ``Makefile`` is also provided. This will compile the whole framework (together with ``utils``) into a single shared library, and link scripts against this, and the necessary third-party libraries. It automatically recognises the program files in directory ``scripts`` as scripts. The ``Makefile`` also features the option ``with-flens``. All other Makefiles have been removed. Note that in contrast to Boost.Build, ``make`` does not provide the possibility of having several build variants simultaneously. With ``make``, the default compilation mode is optimized mode. Type ::

  make <script-name-without-extension>

To switch to debugging mode you need to use ::

  make optimization=no <script-name-without-extension>

Boost.Build, just like ``make``, supports parallel compilation, which can make a significant difference for projects of the magnitude of C++QEDv2. For starting ``n`` threads of compilation use ::

  bjam -j n ... 

C++QEDv2 has been successfully compiled on several Linux platforms and :ref:`Mac OS X (cf. section below) <installingOnMacOsX>`. In all cases the GNU C++ Compiler has been used. It also compiles with the `clang++ <http://clang.llvm.org/>`_ compiler. Portability to other compilers remains to be demonstrated.


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

and also to ``make``::

  make with-flens=no <all the rest as before>

In this case, those parts of the framework that rely on FLENS are discreetly disabled. Most notable is the calculation of the negativity of partially transposed density operators, cf. :ref:`assessingEntanglement`. The file ``utils/src/DrivenDampedHarmonicOscillator.cc`` is also basically disabled, so that :class:`DrivenDampedHarmonicOscillator` becomes unusable.


=========
``utils``
=========

The content of the directory ``utils`` is a small library of very diverse but quite general tools, that I have abstracted during the development of the framework, and used also in several other projects. This may in time become a project on its own. The reader is encouraged to have a look in there, too: some modules may be useful in themselves. Cf. :ref:`cpputils`.


.. _installingOnMacOsX:

=========
Profiling
=========

With Boost.Build, profiling (e.g. with ``gprof``) will never work in release mode because in this mode it automatically adds the ``--strip-all`` option to ``ld``, which removes the symbols necessary for profiling.

Therefore, for profiling, the ``profile`` variant has to be used. Type::

  bjam profile <script-name-without-extension>

The ``Makefile`` also provides the pertaining option. Type::

  make profiling=yes <script-name-without-extension>

.. note::

  With ``make``, be sure that the whole framework gets recompiled. ``bjam`` will anyway put the binaries into separate directories.

============
Mac OS X
============

Relying on `MacPorts <http://www.macports.org/>`_, installation under Mac OS X is straightforward. There are only a few points to take care about, hence we will describe an example procedure.

In this case, the operating system was Snow Leopard, and initially, there were two instances of GCC installed on the system:

Version ``i686-apple-darwin10-gcc-4.2.1``
  The default version

Version ``gcc-mp-4.4 4.4.6.``
  The version from MacPorts

Since MacPorts uses the second one for its installations, to ensure binary compatibility, we will use this one to compile everything else. Install GSL and Boost.Build::

  port install gsl boost-build

Get Blitz++ from the CVS repository, and perform the steps above with the modification ::

  export CXX=g++-mp-4.4; ./configure --with-pic

For simplicity, use the ``BoostIntegration`` branch::

  bzr checkout bzr://cppqed.bzr.sourceforge.net/bzrroot/cppqed/BoostIntegration C++QED

(Bazaar can be obtained from binary installer for Mac.) Copy ``utils/Jamroot.macosx`` to ``utils/Jamroot``. In your home directory, create a ``user-config.jam`` file with the single line ::

  using darwin : 4.4 : gcc-mp-4.4 ;

(Without this, Boost.Build will pass bad options to ``ld``.) Use ::

  bjam with-flens=no <whatever you want to build>

since FLENS is not present on the system.






.. highlight:: c++
  :linenothreshold: 10
