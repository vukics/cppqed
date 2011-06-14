.. _user_guide:

===================
PyCppQED User Guide
===================

This guide explains how to use different features of PyCppQED. For a detailed
documentation about all modules, classes and functions look into the 
:ref:`reference`.

.. contents::
    :depth: 3
    :backlinks: top



Introduction
============

PyCppQED is a python library that helps working with `C++QED`_ - a framework
simulating open quantum dynamics written in C++. Since C++ is not the favorite
programming language of everyone, PyCppQED extends this framework with useful
functionality:

 * Import C++QED output files into python.
 * Export this data as `Matlab`_ ``.mat`` files.
 * Fast and easy visualization/animation of imported data.
 * Generating arbitrary initial conditions for C++QED. 



Requirements
============

**Mandatory:**

    * `Python`_:
        The python interpreter. At least version 2.3 but not 3.x.

    * `NumPy`_:
        A numerical package for python. Provides fast N-dimensional array
        manipulation.

**Optional:**

    * `Matplotlib`_:
        A plotting library. It is needed for any kind of visualization. (For
        3D plots at least 0.99 is needed)

    * `SciPy`_:
        Library providing scientific algorithms. It is needed for exporting
        numpy arrays as ``.mat`` files and importing ``.mat`` files as numpy
        arrays. For exporting multidimensional arrays, at least version 0.7 is
        needed.

    * `PyGTK`_:
        GTK bindings for python. It is needed for animations.



Installation
============

First, the c extensions have to be build::

    $ python setup.py build

This creates a directory like :file:`build/lib.linux-x86_64-2.4/pycppqed/`.
Either this package can be moved somewhere and used directly (you may want to
add it's location to the :envvar:`PYTHONPATH`) or alternatively it can be
installed::

    $ python setup.py install --prefix=PATH



Overview
========

PyCppQED has a strict modular design:

 * Python classes representing objects used in QED:
     * :mod:`pycppqed.statevector` implements state vectors. 

     * :mod:`pycppqed.expvalues` implements classes for working with
       expectation values.

     * :mod:`pycppqed.quantumsystem` implements classes representing
       quantum systems.

 * Functions for generating some useful initial conditions are in
   :mod:`pycppqed.initialconditions`.

 * Everything that has to do with reading and writing C++QED files is in
   the module :mod:`pycppqed.io`.

 * Plotting stuff is in :mod:`pycppqed.visualization` and animation functions
   are implemented in :mod:`pycppqed.animation`.



Usage
=====

PyCppQED can be used either from scripts but also interactively.


Scripts
-------

Scripts are simple text files with valid python code that can be executed by
invoking the python interpreter with the name of the script as first argument::

    $ python myscript.py


Interactive
-----------

To use python interactively just invoke the interpreter without arguments::
    
    $ python
    Python 2.6.2 (r262:71600, Aug 17 2009, 10:52:48)
    [GCC 4.1.2 20080704 (Red Hat 4.1.2-44)] on linux2
    Type "help", "copyright", "credits" or "license" for more information.
    >>>

A good enhancement to the standard python interpreter is `IPython`_::

    $ ipython
    Python 2.4.3 (#1, Jul 27 2009, 17:56:30)
    Type "copyright", "credits" or "license" for more information.

    IPython 0.8.4 -- An enhanced Interactive Python.
    ?         -> Introduction and overview of IPython's features.
    %quickref -> Quick reference.
    help      -> Python's own help system.
    object?   -> Details about 'object'. ?object also works, ?? prints more.

    In [1]:

It provides:

    * **Tab completion** which allows easy inspection of modules and objects::
        
        In [3]: file.re
        file.read       file.readinto   file.readline   file.readlines

        In [3]: file.read

        
    * **Easy Inspection**::

        In [3]: file.read?
        Type:           method_descriptor
        Base Class:     <type 'method_descriptor'>
        String Form:    <method 'read' of 'file' objects>
        Namespace:      Python builtin
        Docstring:
            read([size]) -> read at most size bytes, returned as a string.

            If the size argument is negative or omitted, read until EOF is
            reached. Notice that when in non-blocking mode, less data than what
            was requested may be returned, even if no size parameter was given.

    * **Many more features ...**


Import PyCppQED
---------------

To use PyCppQED you have to import it. This is done using an import statement::

    >>> import pycppqed as qed

From now on all commands will assume that PyCppQED is already imported with
this statement. Now you are ready to do a lot of fancy stuff with PyCppQED!
The next section gives examples how to achieve common tasks.



How to ...
==========

Split up a C++QED output file into standard output and state vectors
--------------------------------------------------------------------

When a C++QED script is invoked using the :option:`svdc` argument, state vectors
are written into the output file between the calculated expectation values.
With PyCppQED it's easy to extract the state vectors into own files and
getting a standard C++QED output file::

    >>> qed.io.split_cppqed("ring.dat", "newring.dat")

This writes the standard output file to :file:`newring.dat` and the state
vectors into separate files named :file:`newring_{time}.dat.sv` where
:token:`time` is substituted by the time when this state vector was reached.


.. _import2python:

Import a C++QED output file
---------------------------

This is done with the function :func:`pycppqed.io.load_cppqed`::

    >>> evs, qs = qed.io.load_cppqed("ring.dat")

This returns two objects which represent the whole information stored
in the C++QED output file:

 * A :class:`pycppqed.expvalues.ExpectationValueCollection` instance which
   holds all expectation values calculated by C++QED.

 * A :class:`pycppqed.quantumsystem.QuantumSystemCompound` instance
   representing the calculated quantum system. This object also stores a 
   :class:`pycppqed.statevector.StateVectorTrajectory` instance which holds
   all calculated state vectors.


Export python data as *.mat* file
---------------------------------

If you want to use `Matlab`_ or `Octave`_ for further processing of the data
you can use PyCppQED to convert a C++QED output file into a *.mat* file.
First, we have to load the file like in :ref:`import2python`. The obtained 
objects (or only parts of it, or any other array ...) can be saved with
the :func:`scipy.io.savemat` function::

    >>> import scipy.io
    >>> scipy.io.savemat("out.mat", {"evs":evs, "svs":qs.statevector})

This file can be used from `Matlab`_ and `Octave`_:

.. code-block:: matlab

    >> load('out.mat')
    >> whos
      Name       Size              Bytes  Class     Attributes

      evs       15x175             21000  double
      svs        4-D              921600  double    complex

    >> size(svs)

    ans =

         9    64    10    10


    >>> size(evs)

    ans =
        
        15   175

.. warning::

    Be aware that old versions of scipy (older than 0.7) can't properly export
    arrays with more than 2 dimensions!


Generate arbitrary initial conditions
-------------------------------------

In the module :mod:`pycppqed.initialconditions` are some convenient functions
that let you easily create common initial conditions. E.g. to create a
gaussian wave packet in the k-space the following command can be used::
    
    >>> sv_p = qed.initialconditions.gaussian(x0=1.1, k0=5, sigma=0.3, fin=7)
    >>> print sv_p
    StateVector(128)

Or a coherent mode::

    >>> sv_m = qed.initialconditions.coherent(alpha=2, N=20)
    >>> print sv_m
    StateVector(20)

To obtain initial conditions for a combined quantum system simply use the
**\^** operator::

    >>> sv = sv_a ^ sv_m
    >>> print sv
    StateVector(128 x 20)

It's easy to create any other initial condition you can think of, by simply
creating a numpy array with the wanted values and then using the array to
build a :class:`pycppqed.statevector.StateVector`::

    >>> import numpy as np
    >>> X = np.linspace(0,10,64) # An array with 64 values between 0 and 10
    >>> Y = np.sin(X)
    >>> sv = qed.statevector.StateVector(Y, norm=True)
    >>> print sv
    StateVector(64)


Export initial conditions as C++QED sv files
--------------------------------------------

Exporting is done with the :func:`pycppqed.io.save_statevector`::

    >>> sv = qed.statevector.StateVector((1,2,3), norm=True)
    >>> qed.io.save_statevector("mystatevector.sv", sv)

The created file then looks like::

    # 0 1
    (0,2)
    [ (0.267261241912,0.0) (0.534522483825,0.0) (0.801783725737,0.0) ]


Calculate standard expectation values
-------------------------------------

By *standard* expectation values we mean values that are also calculated by
C++QED. Automatic calculation for those is implemented in
:mod:`pycppqed.quantumsystem`. All that has to be done is to first create a
proper quantum system object and then call its :meth:`expvalues` method::

    >>> sv = qed.initialconditions.gaussian(x0=0.5, k0=3.2, sigma=0.4, fin=7)
    >>> qs = qed.Particle(sv)
    >>> evs = qs.expvalues()
    >>> print evs
    ExpectationValueCollection('<k>', 'Var(k)', '<x>', 'Std(x)')
    ExpectationValueCollection([ 3.2000+0.j,  1.5625+0.j,  0.5000+0.j, 0.4000+0.j])

It's also possible for combined quantum systems::

    >>> sv_p = qed.initialconditions.gaussian(x0=0.5, k0=3.2, sigma=0.4, fin=7)
    >>> sv_m = qed.initialconditions.coherent(alpha=2, N=20)
    >>> sv = sv_p ^ sv_m
    >>> q = qed.quantumsystem
    >>> qs = q.QuantumSystemCompound(sv, q.Particle, q.Mode)
    >>> print qs
    QuantumSystemCompound(Particle(128), Mode(20))
    >>> print evs
    ExpectationValueCollection('<k>', 'Var(k)', '<x>', 'Std(x)', '<n>', 'Var(n)', 'Re(<a>)', 'Im(<a>)')
    >>> print repr(evs)
    ExpectationValueCollection([  3.19999997e+00+0.j,   1.56250009e+00+0.j,   4.99999995e-01+0.j,  4.00000001e-01+0.j,   3.99999979e+00+0.j,   3.99999747e+00+0.j,   1.99999990e+00+0.j,   6.41996804e-18+0.j])
    

We get a quantum system for free if we load a C++QED output file::

    >>> evs, qs = qed.io.load_cppqed("ring.dat")
    >>> evs_calc = qs.expvalues()


Calculate arbitrary expectation values
--------------------------------------

Expectation values for combined systems are calculated in the following way
(Assuming the operator only acts on first subsystem):

    .. math::

        \langle \Psi | \hat A (k) | \Psi \rangle =
                \sum_{k_1 k_2} \langle k_1 | \hat A (k) | k_2 \rangle
                \sum_m \Psi_{k_1 m}^* \Psi_{k_2 m}

That means the expectation value is determined by specifying the quantity
:math:`A_{k_1 k_2} = \langle k_1 | \hat A (k) | k_2 \rangle`. E.g. let's 
calculate the expectation value of the destruction operator of a combined
system of structure {Particle, Mode}::

    >>> sv_p = qed.initialconditions.gaussian(x0=0.5, k0=3.2, sigma=0.4, fin=7)
    >>> sv_m = qed.initialconditions.coherent(alpha=0.5, N=5)
    >>> sv = sv_p ^ sv_m
    >>> import numpy as np
    >>> a = np.diag(np.sqrt(np.arange(1,5)), 1)
    >>> print a
    [[ 0.          1.          0.          0.          0.        ]
     [ 0.          0.          1.41421356  0.          0.        ]
     [ 0.          0.          0.          1.73205081  0.        ]
     [ 0.          0.          0.          0.          2.        ]
     [ 0.          0.          0.          0.          0.        ]]
    >>> ev_a = sv.expvalue(a, 1)
    >>> print ev_a
    (0.499933315175-7.96953264544e-18j)

The second argument tells the expvalue method that the operator is only working
on the second subsystem. (Python starts counting with 0!)

Let's now consider a slightly more complicated example - a combined system of
the form {Particle, Mode, Mode} and let's try to calculate the expectation
value for the operator :math:`T = a_1^\dagger a_2 + a_2^\dagger a_1`::

    >>> sv_p = qed.initialconditions.gaussian(x0=0.5, k0=3.2, sigma=0.4, fin=7)
    >>> sv_m1 = qed.initialconditions.coherent(alpha=0.5, N=5)
    >>> sv_m2 = qed.initialconditions.coherent(alpha=2, N=20)
    >>> sv = sv_p ^ sv_m1 ^ sv_m2
    >>> import numpy as np
    >>> SV = qed.statevector.StateVector # Define an abbreviation.
    >>> m1_a = SV(np.diag(np.sqrt(np.arange(1,5)), 1))
    >>> m1_at = SV(np.diag(np.sqrt(np.arange(1,5)), -1))
    >>> m2_a = SV(np.diag(np.sqrt(np.arange(1,20)), 1))
    >>> m2_at = SV(np.diag(np.sqrt(np.arange(1,20)), -1))
    >>> T = (m1_at ^ m2_a) + (m1_a ^ m2_at)
    >>> print sv.expvalue(T, (0,1))
    (1.99973315754-1.54328182927e-20j)


Calculate diagonal expectation values
-------------------------------------

If we want to calculate the expectation value of an diagonal operator, we can
use the :meth:`pycppqed.statevector.StateVector.diagexpvalue` method. It takes
only the diagonal elements of the operator's matrix representation and has
the advantage to be faster and to need less memory than the
:meth:`pycppqed.statevector.StateVector.expvalue` method.

A short example::

    >>> sv = qed.initialconditions.gaussian(x0=0.5, k0=3.2, sigma=0.4, fin=7)
    >>> import numpy as np
    >>> K = np.arange(-64, 64)
    >>> print sv.diagexpvalue(K, 0)
    (3.2+0j)


Visualize PyCppQED objects
--------------------------

There are basically 4 different types of PyCppQED objects which might be
interesting to look at:

    * :class:`pycppqed.statevector.StateVector`
    * :class:`pycppqed.statevector.StateVectorTrajectory`
    * :class:`pycppqed.expvalues.ExpectationValueTrajectory`
    * :class:`pycppqed.expvalues.ExpectationValueCollection`

All of them are inheriting from :class:`numpy.ndarray` which means you can
easily plot them using `Matplotlib`_ or `Gnuplot`_. However, PyCppQED
implements some functions to let you take a quick look on these objects. All
but the StateVectorTrajectory class have a :meth:`plot` method::

    >>> sv = qed.initialconditions.coherent(alpha=2.3, N=25)
    >>> sv.plot()

.. image:: media/graph_coherent.png
    :width: 8cm
    :height: 6cm

To change the x-axis we can pass an array of x-coordinated to the plot method::

    >>> import numpy as np
    >>> sv = qed.initialconditions.gaussian(x0=0.5, k0=3.2, sigma=0.05, fin=7)
    >>> K = np.arange(-64,64)
    >>> sv.plot(x=K)

.. image:: media/graph_gaussian.png
    :width: 8cm
    :height: 6cm

Since Matplotlib version 0.99 it's also possible to draw 3D graphs. This
can be used for combined systems::

    >>> sv_p = qed.initialconditions.gaussian(x0=0.5, k0=10, sigma=0.1, fin=6)
    >>> sv_m = qed.initialconditions.coherent(alpha=1.5, N=15)
    >>> sv = sv_p ^ sv_m
    >>> import numpy as np
    >>> K = np.arange(-32,32)
    >>> sv.plot(x=K)

.. image:: media/graph_gaussian&mode.png
    :width: 8cm
    :height: 6cm

The expectation value classes work equivalent. Maybe also useful is the
function :func:`pycppqed.visualization.compare_expvaluecollections`. As its
name says it is used to compare two sets of expectation values::

    >>> evs, qs = qed.io.load_cppqed("ring.dat")
    >>> evs_calc = qs.expvalues()
    >>> qed.visualization.compare_expvalluecollections(evs, evs_calc)

========= ========= =========
|expval1| |expval2| |expval3|
========= ========= =========

.. |expval1| image:: media/graph_expvals1.png
    :width: 5cm
    :height: 8cm

.. |expval2| image:: media/graph_expvals2.png
    :width: 5cm
    :height: 8cm
.. |expval3| image:: media/graph_expvals3.png
    :width: 5cm
    :height: 8cm


The only object that is now left, is the 
:class:`pycppqed.statevector.StateVectorTrajectory` class. It represents the
time evolution of a state vector. For a 1D system there are already three
dimensions to plot. This would be possible, but an alternative is to use
an animation which will also work for 2D systems.
Animations are implemented as interactive window, but it's also possible
to save movies in any format mencoder can write. This functionality
is only very basic and it may need changes on the source code to obtain
professional looking movies. However, here is the code::

    >>> import numpy as np
    >>> time = np.linspace(-np.pi, np.pi, 2**6)
    >>> g = qed.initialconditions.gaussian
    >>> svt = [g(sigma=t) ^ g(sigma=0.15-t) for t in time]
    >>> svt.animate()

And here is the example movie: `animation.avi <_static/animation.avi>`_



.. _C++QED: http://sourceforge.net/projects/cppqed
.. _Python: http://python.org
.. _IPython: http://ipython.scipy.org/moin
.. _Matlab: http://www.mathworks.com
.. _Octave: http://www.gnu.org/software/octave
.. _Gnuplot: http://gnuplot.info
.. _Matplotlib: http://matplotlib.sourceforge.net
.. _SciPy: http://scipy.org
.. _NumPy: http://scipy.org
.. _PyGTK: http://pygtk.org
