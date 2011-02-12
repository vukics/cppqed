.. _quantumoperator:

====================================
The ``quantumoperator`` namespace
====================================

-------------------------------------
Representing "tridiagonal" matrices
-------------------------------------


Synopsis
^^^^^^^^^

Let us consider the following example, which displays the full potential of classes :class:`~quantumoperator::Tridiagonal` and :class:`~quantumoperator::Frequencies`: We have a system consisting of :math:`M` subsystems, each featuring a free (in general, non-Hermitian) Hamiltonian, which is diagonal in the given system's Hilbert space. In addition, we have one more term of a special form in the Hamiltonian coupling all the subsystems:

.. math::

  H=H^\text{free}+H^\text{interaction}=&\sum_{m=0}^{M-1}\bigotimes_{k=0}^{m-1}\mathbf{1}_k\otimes\lp\sum_{n_m=0}^{N_m-1}\omega_{m,n}\ket{n_m}\bra{n_m}\rp\otimes\bigotimes_{k=m+1}^{M-1}\mathbf{1}_k\\&+\bigotimes_{m=0}^{M-1}\lp\sum_{n_m=0}^{N_m-1}\alpha^{(0)}_{m,n}\ket{n_m}\bra{n_m}\right.+\left.\sum_{n_m=0}^{N_m-1-K_m}\lp\alpha^{(+)}_{m,n}\ket{n_m}\bra{n_m+K_m}+\alpha^{(-)}_{m,n}\ket{n_m+K_m}\bra{n_m}\rp\rp.

Here the coefficients :math:`\omega` and :math:`\alpha` are in general complex with the dimension of frequency (:math:`\hbar=1`). The :math:`N` are the dimensions of the subsystems' Hilbert space in which the vectors :math:`\ket{n}` form an orthonormal basis.

The Hamiltonian is indeed in a special form because the interaction term is a direct product of such operators acting on the Hilbert spaces of the individual subsystems, whose matrix contains only three diagonals of nonzero elements. Hence the name of the class :class:`~quantumoperator::Tridiagonal`, although this usually refers to the case when :math:`K=1`.

Now let us transform to the interaction picture defined by :math:`H^\text{free}`. The Hamiltonian in interaction picture reads

.. math::
  :label: HITridiag

  H_I(t)=\bigotimes_{m=0}^{M-1}\lp\sum_{n_m=0}^{N_m-1}\alpha^{(0)}_{m,n}\ket{n_m}\bra{n_m}+\sum_{n_m=0}^{N_m-1-K_m}\lp
  e^{\delta_{m,n}t}\alpha^{(-)}_{m,n}\ket{n_m+K_m}\bra{n_m}+e^{-\delta_{m,n}t}\alpha^{(+)}_{m,n}\ket{n_m}\bra{n_m+K_m}\rp\rp,

where :math:`\delta_{m,n}=i\lp\omega_{m,n+K_m}-\omega_{m,n}\rp`.

Quite generally, classes :class:`~quantumoperator::Tridiagonal` and :class:`~quantumoperator::Frequencies` are designed to store and manipulate Hamiltonians of the form :eq:`HITridiag` for an arbitrary number of subsystems. In particular, the Hamiltonian can be evaluated at a given :math:`t`, applied on a state vector and combined with other :class:`~quantumoperator::Tridiagonal`\ s using algebraic operations.

As their names reflect, a pair of objects of classes :class:`~quantumoperator::Tridiagonal` and :class:`~quantumoperator::Frequencies` correspond to a given Hamiltonian (or a given term of a more complex Hamiltonian) of the form :eq:`HITridiag`, the :class:`~quantumoperator::Tridiagonal` object roughly speaking storing the coefficients :math:`\alpha^{(0,-,+)}_{m,n}`, while the :class:`~quantumoperator::Frequencies` object the frequencies :math:`\delta_{m,n}` for all :math:`m,n`. In itself a :class:`~quantumoperator::Tridiagonal` object therefore reflects the matrix of the operator :eq:`HITridiag` at :math:`t=0`, while it can be updated using the corresponding :class:`~quantumoperator::Frequencies` object to some arbitrary time instant :math:`t`.


:class:`~quantumoperator::Tridiagonal`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _tridiagTemplates:

======== ========== ==================
Template parameter  descriptions
======== ========== ==================
int      RANK       Positive integer standing for the number of elementary Hilbert spaces in :eq:`HITridiag`
======== ========== ==================


.. class:: quantumoperator::Tridiagonal

  :ref:`template parameters <tridiagTemplates>`: RANK

  .. c:var:: LENGTH
  
    The number of :type:`Diagonal`\ s the class has to store::

      static const int LENGTH=tmptools::Power<3,RANK>::value;

  .. type:: Diagonals 
  
    The class is implemented in terms of a :class:`blitzplusplus::TinyOfArrays`, this is the class used to store the :type:`Diagonal`\ s::

      typedef blitzplusplus::TinyOfArrays<dcomp,RANK,LENGTH> Diagonals;

  .. type:: Diagonal

     ::
     
       typedef typename Diagonals::T_numtype Diagonal;

  .. function:: explicit Tridiagonal(const Diagonal& zero =empty, size_t k =0, const Diagonal& minus =empty, const Diagonal& plus =empty, mpl::int_<RANK> = _1_)

    This is the principial way to create an object of this class, which can be used for ``RANK=1`` only, as ensured by the trailing dummy argument. This creates an object corresponding to the elementary operator

    .. math::
       :label: ElemTridiag

       H_I^\text{elem}(0)=\sum_{n=0}^{N-1}\alpha^{(0)}_n\ket{n}\bra{n}+\sum_{n=0}^{N-1-K}\lp\alpha^{(-)}_n\ket{n+K}\bra{n}+\alpha^{(+)}_n\ket{n}\bra{n+K}\rp

    The arguments ``zero``, ``minus``, ``plus``, and ``k`` correspond respectively to :math:`\alpha^{(0)}`, :math:`\alpha^{(-)}`, :math:`\alpha^{(+)}`, and :math:`K`

.. class:: quantumoperator::Frequencies

  :ref:`template parameters <tridiagTemplates>`: RANK



.. warning::

  The following is an outdated version:

Creation
%%%%%%%%%%


\normalsize
where the arrays \code{Zero}, \code{Minus}, and \code{Plus} of type 
\small
\begin{verbatim}
   typedef std::valarray<std::complex<double> > VADC;
\end{verbatim}
\normalsize
have to be previously prepared to contain the coefficients
:math:`\alpha^{(0,-,+)}_n`, respectively. It follows that the size of the first one
is :math:`N`, while that of the other two is :math:`N-K`. If a diagonal is
zero an empty array may be supplied (rather than one of the
appropriate size full of zeros).

Perhaps one of the simplest nontrivial example is the annihilation
operator :math:`a=\sum_{n=0}^{N-1}\sqrt{n+1}\ket n\bra{n+1}` of a harmonic
oscillator. Here clearly :math:`\alpha^{(+)}_n=\sqrt{n+1}` and
:math:`K=1`. Hence a function returning a :class:`~quantumoperator::Tridiagonal` representing
operator :math:`a` for a given cutoff :math:`N` could read
\small
\begin{verbatim}
   Tridiagonal aop(size_t N)
   {
     VADC alphaPlus(N-1);
     for (size_t n=0; n<N-1; n++) alphaPlus[n]=sqrt(n+1.);
     return Tridiagonal(VADC(),1,VADC(),alphaPlus);
   }
\end{verbatim}
\normalsize
The element \code{LossyMode} indeed defines a function in this way.

A copy constructor is defined and a function
\small
\begin{verbatim}
   Tridiagonal Identity(size_t N);
\end{verbatim}
\normalsize
which creates a :class:`~quantumoperator::Tridiagonal` representing identity. In addition, many
elements define :class:`~quantumoperator::Tridiagonal`\ s corresponding to their characteristic
operators. More complex :class:`~quantumoperator::Tridiagonal` s can be created as direct products of
elementary operators.


An elementary :class:`~quantumoperator::Frequencies` object corresponding to the operator
(\ref{eq:ElemTridiag}) can be created using the constructor
\small
\begin{verbatim}
   Frequencies::Frequencies(const VADC& Delta);
\end{verbatim}
\normalsize
where the valarray \code{Delta} has to be previously filled with the
frequencies :math:`\delta_n`.

Operations
%%%%%%%%%%%%

Vector-space operations (addition and multiplication with a complex
number) are defined both in the form \code{+=} and \code{*=} as
member functions and \code{+} and \code{*} as helper functions.
Hermitian transposition is performed by the member function
\small
\begin{verbatim}
   Tridiagonal Tridiagonal::transpose() const;
\end{verbatim}
\normalsize

The tool for creating complex operators is direct product of more
elementary operators, which is also defined in both form \code{*=}
and \code{*}. Together with multiplying :class:`~quantumoperator::Tridiagonal` s using direct
product their corresponding :class:`~quantumoperator::Frequencies` objects have also to be
multiplied in the very same order. (A complex example is given below in
Sec.~\ref{sec:TridiagFreqsExample}.) Operators \code{*=} and
\code{*} are defined accordingly for the :class:`~quantumoperator::Frequencies` class.

Note that normal operator multiplication is not defined because the
product of two operators of the form represented by :class:`~quantumoperator::Tridiagonal` is
in general not of such form. In other words, operators of such form do
not form an algebra. Nonetheless, if really needed, the same effect
can always be achieved by simply \code{Apply}ing (cf below) the
operators succesively --- an example can be found in the
implementation of element \code{ParticleTwoModes}.

The operations described so far are not really meant to be used (as
are not really needed, either) during run-time of the simulation, only
during startup time: In their implementation no special emphasis has
been put on efficiency (eg in many of them :class:`~quantumoperator::Tridiagonal` s are actually
\emph{copied}) but rather on elegance and clarity. There are two
operations, however, which are used abundantly during run-time
(notably by \code{TridiagonalHamiltonian} objects):
The member function
\small
\begin{verbatim}
   void Tridiagonal::Propagate(const Frequencies& F, double Delta_t);
\end{verbatim}
\normalsize uses the frequencies :math:`\delta` stored in the :class:`~quantumoperator::Frequencies` object to update
the :class:`~quantumoperator::Tridiagonal` object representing an operator (\ref{eq:HITridiag}) at time
:math:`t` to time :math:`t+\Delta t`. 

They don't really have to match ()!!!!!! It is now checked for in
\code{Propagate}!!!!!! This solves a lot of problems, eg I don't need
to consider tim dependent and time independent problems separately.


That the structure and dimensions of the
two objects have to match, and this is \emph{not checked for at all}
due to efficiency reasons, therefore violation of this leads to
undefined behaviour.

Finally, the helper function
\small
\begin{verbatim}
   void Apply(const double* Psi, double* dPsidt,
              const Tridiagonal& T, const CPA_Slice&);
\end{verbatim}
\normalsize

Apply, ExpUpdate: runtime 

Example: pumped lossy mode directly and with 'a'.


\code{Apply}: recursive behaviour

Common interface with other operator classes (eg \code{CMatrix})? eg
:class:`~quantumoperator::Tridiagonal` multiplication yielding CMatrix --- we want to avoid
this, because of accidents ---
Does not really make sense, because interaction picture is really
powerful only in case of (very) sparse Hamiltonians. (Not clear how to
define ExpUpdate for CMatrix)

Note that most Hamiltonians in the framework and in quantum
optics in general can be written in such a form.

Exceptions
%%%%%%%%%%%%

If the sizes are
inconsistent, an exception of 

Examples
^^^^^^^^^^

Implementation
^^^^^^^^^^^^^^^^

