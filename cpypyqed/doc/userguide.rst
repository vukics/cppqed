.. _userguide:

===================
cpypyqed User Guide
===================

This guide explains how to write C++QED scripts in Python. For a detailed
documentation about all modules, classes and functions look into the
:ref:`reference`.

Concept
=======

Cpypyqed is a thin python wrapper to the C++QED programming framework. The goal is to enable the user
to write python scripts instead of C++ programs, while keeping the structure and syntax very similar.

All the actual work of running teh simulations is done by the C++ libraries under the hood, which means
there is no performance penalty when using cpypyqed over plain C++QED.

C++QED makes heavy use of deeply templated classes and compile-time algorithms. Some of the class templates
are pre-instantiated for cpypyqed, e.g. there is :class:`~.quantumdata.StateVector1`,
:class:`~.quantumdata.StateVector2`, ... :class:`~.quantumdata.StateVector4` as instantiations for
a :core:`quantumdata::StateVector` up to rank 4. For other classes which are templated in a more complex
way, this pre-instantiation is not possible.

Instead, cpypyqed makes use of the cmake build infrastructure of C++QED to build wrappers of these class templates
on demand. The modules resulting out of these on-demand compilation are cached and available for use without
compilation delay from then on.

Installation
============

Easiest installation is to clone the C++QED super repository which comes with everything.
If Python and boost.python is detected, cpypyqed will be built::

  $ git clone --recursive -b python git@ge-c705.uibk.ac.at:cppqed/complete.git C++QED
  $ cd C++QED
  $ make C++QED/build; cd C++QED/build
  $ cmake -DCMAKE_BUILD_TYPE=<Release or Debug> <other cmake options> ..
  $ make cpypyqed
  $ make check_cpypyqed

The last command will run the cpypyqed related tests.

Of course it is also possible to clone the cpypyqed repository as a standalone project.

After building the module, either call :samp:`make install` or point your :samp:`$PYTHONPATH` to the
cpypyqed build directory. As will be explained later, you can also have two cpypyqed build directories installed
or in the :samp:`$PYTHONPATH` simultaneously: one compiled in release mode and one compiled in debug mode.

On-demand compilation
=====================

By default, on-demand compilation is performed in a temporary directory under :file:`~/.cpypyqed` and the
resulting modules are stored in :file:`~/.cpypyqed/cppqedmodules`. The build directory :file:`~/.cpypyqed` can be
changed by setting the environment variable :samp:`$CPYPYQED_BUILDDIR`. Still we will refer to the build directory
as :file:`~/.cpypyqed`.

Configuration file
------------------

By default the configuration file read in by cpypyqed is :file:`~/.cpypyqed/config.txt`. An additional configure
file can be read in by setting the environment variable :samp:`$CPYPYQED_CONFIG`. The configuration file
has the following syntax::

  [Setup]
  # override the location where cpypyqed looks for C++QED libraries (e.g. CPPQEDcore build directory)
  cppqed_dir=<path>
  # override the location where cpypyqed looks for release C++QED libraries (e.g. CPPQEDcore build directory)
  # this also overrides cppqed_dir
  cppqed_dir_debug=<path>
  # override the location where cpypyqed looks for debug C++QED libraries (e.g. CPPQEDcore build directory)
  # this also overrides cppqed_dir
  cppqed_dir_release=<path>

  # Set the C++ compiler path
  compiler=<path>
  # Set additional cmake options
  cmake_opts=-Dsome_option -Dsome_other_option

  # If set to False, the temporary build directory is not deleted after the build
  delete_temp=True

Loading cpypyqed
==================

Because cmake allows different build configurations, there are two versions of cpypyqed:

  * :samp:`import cpypyqed` or :samp:`from cpypyqed import *` for release mode
  * :samp:`import cpypyqed.debug` or :samp:`from cpypyqed.debug import *` for debug mode

If one or both of these import statements succeed depend on which configurations are installed.

Note that importing :samp:`cpypyqed` will fall back to the debug version if the release version is not found,
so this should work in both cases. Use :samp:`cpypyqed.debug` to force the debug version.
If you want to choose between debug and release mode upon calling the script,
here is one possibility with :py:mod:`argparse`::

  import sys
  import argparse

  parser = argparse.ArgumentParser(add_help=False)
  parser.add_argument('--debug', action='store_true')
  (args,remaining)=parser.parse_known_args(sys.argv)

  if vars(args)['debug']:
      from cpypyqed.debug import *
  else:
      from cpypyqed import *

This will look for the command line argument :samp:`--debug` and load :samp:`cpypyqed_d` instead of :samp:`cpypyqed`
if it finds it. All the other command line arguments are stored in :samp:`remaining` for further inspection later.
This is the method all the example scripts of the cpypyqed package use.

Command line arguments
======================

In C++QED, the parameter bundle is responsible to declare parameters, define their defaults
and finally parse and update parameters from the command line. The parameter bundle is exposed in the cpypyqed interface by
:class:`.ParameterTable`, :meth:`.update` and all the relevant :samp:`Pars`-classes like :class:`.mode.Pars`,
:class:`.qbit.Pars`, etc. Usage is pretty much along the lines of the C++QED parameter bundle with two exceptions:
the :meth:`.update` function has a different signature as documented, and
:core:`parameters::ParameterTable::add <parameters::ParameterTable::add(const std::string &, const std::string &, const T &)>`
is missing in cpypyqed for technical reasons. This means additional user-defined commandline parameters cannot be added
in the usual way. However, there is a workaround using the aforementioned :py:mod:`argparse`. This is demonstrated
in the script :file:`1particle1mode.py`::

  import sys
  import argparse

  parser = argparse.ArgumentParser(add_help=False)
  parser.add_argument('--debug', action='store_true')
  parser.add_argument('--1p1mconf', help="System configuration code for 1particle1mode",
                      type=int,default=1)
  (args,remaining)=parser.parse_known_args(sys.argv)

  if vars(args)['debug']:
      from cpypyqed.debug import *
  else:
      from cpypyqed import *

  conf = vars(args)['1p1mconf']

  ...

  parameters.update(p,remaining,'--')

Example scripts
===============

These Python scripts correspond to the examples presented in the `User Guide for the C++ interface <../../cppqed/html/userguide.html>`_

The single-harmonic-oscillator-mode example
-------------------------------------------

Corresponds to `this C++ script <../../cppqed/html/userguide.html#userguideelementaryparameters>`_

.. literalinclude:: /../../../cpypyqed/scripts/tutorialMode.py
  :language: python
  :linenos:

The binary-system example
-------------------------

Corresponds to `this C++ script <../../cppqed/html/userguide.html#userguidebinaryfullfledged>`_

.. literalinclude:: /../../../cpypyqed/scripts/tutorialBinary.py
  :language: python
  :linenos:

Composite examples
------------------

Ring cavity
^^^^^^^^^^^

Corresponds to `this C++ script <../../cppqed/html/userguide.html#userguidemorecomplexringfullfledged>`_

.. literalinclude:: /../../../cpypyqed/scripts/tutorialCompositeRing.py
  :language: python
  :linenos:

Multi-particle example
^^^^^^^^^^^^^^^^^^^^^^

…
