.. _quantumoperator:

*********************************
The ``quantumoperator`` namespace
*********************************

===================================
Representing "tridiagonal" matrices
===================================


--------
Synopsis
--------

Let us consider the following situation, which displays the full potential of the :class:`~quantumoperator::Tridiagonal` class: We have a system consisting of :math:`M` subsystems, each featuring a free (in general, non-Hermitian) Hamiltonian, which is diagonal in the given system's Hilbert space. In addition, we have one more term of a special form in the Hamiltonian coupling all the subsystems:

.. math::
  :label: HTridiag

  H=H^\text{free}+H^\text{interaction}=&\sum_{m=0}^{M-1}\bigotimes_{k=0}^{m-1}\mathbf{1}_k\otimes\lp\sum_{n_m=0}^{N_m-1}\omega_{m,n}\ket{n_m}\bra{n_m}\rp\otimes\bigotimes_{k=m+1}^{M-1}\mathbf{1}_k\\&+\bigotimes_{m=0}^{M-1}\lp\sum_{n_m=0}^{N_m-1}\alpha^{(0)}_{m,n}\ket{n_m}\bra{n_m}\right.+\left.\sum_{n_m=0}^{N_m-1-K_m}\lp\alpha^{(+)}_{m,n}\ket{n_m}\bra{n_m+K_m}+\alpha^{(-)}_{m,n}\ket{n_m+K_m}\bra{n_m}\rp\rp.

Here, the coefficients :math:`\omega` and :math:`\alpha` are in general complex with the dimension of frequency (:math:`\hbar=1`). The :math:`N_m`\ s are the dimensions of the subsystems' Hilbert spaces in which the vectors :math:`\ket{n_m}` form an orthonormal basis.

The Hamiltonian is indeed in a special form because the interaction term is a direct product of such operators acting on the Hilbert spaces of the individual subsystems, whose matrix contains only three diagonals of nonzero elements. Hence the name of the class :class:`~quantumoperator::Tridiagonal`, although this usually refers to the case when :math:`K=1`.

Now let us transform to the interaction picture defined by :math:`H^\text{free}`. The Hamiltonian in interaction picture reads

.. math::
  :label: HITridiag

  H_I(t)=\bigotimes_{m=0}^{M-1}\lp\sum_{n_m=0}^{N_m-1}\alpha^{(0)}_{m,n}\ket{n_m}\bra{n_m}+\sum_{n_m=0}^{N_m-1-K_m}\lp
  e^{\delta_{m,n}t}\alpha^{(-)}_{m,n}\ket{n_m+K_m}\bra{n_m}+e^{-\delta_{m,n}t}\alpha^{(+)}_{m,n}\ket{n_m}\bra{n_m+K_m}\rp\rp,

where :math:`\delta_{m,n}=i\lp\omega_{m,n+K_m}-\omega_{m,n}\rp`.

Quite generally, the :class:`~quantumoperator::Tridiagonal` class is designed to store and manipulate Hamiltonians of the form either :eq:`HTridiag` or :eq:`HITridiag` for an arbitrary number of subsystems. In the latter case, it also stores and manipulates the :math:`\delta_{m,n}` frequencies. In particular, the Hamiltonian can be evaluated at a given time :math:`t`, applied on a state vector and combined with other :class:`~quantumoperator::Tridiagonal`\ s using algebraic operations.

Tridiagonal internally bookkeeps to which time instant its state corresponds.

Policy of storing diagonals: only the non-zeros. Policy of storing frequencies: either none or all of them.

--------------
Implementation
--------------

.. _tridiagTemplates:

.. note::

  Template argument definitions:

  ``int RANK``
    Positive integer standing for the number :math:`M` of elementary Hilbert spaces in :eq:`HITridiag`


.. class:: quantumoperator::Tridiagonal

  ``template <int RANK>`` (cf. :ref:`template parameters <tridiagTemplates>`)

  .. c:var:: LENGTH
  
    The number of :type:`Diagonal`\ s the class has to store::

      static const int LENGTH=tmptools::Power<3,RANK>::value;

  .. type:: Diagonals 
  
    The class is implemented in terms of a :class:`blitzplusplus::TinyOfArrays`, this is the class used to store the :type:`Diagonal`\ s::

      typedef blitzplusplus::TinyOfArrays<dcomp,RANK,LENGTH> Diagonals;

  .. type:: Diagonal

     ::
     
       typedef typename Diagonals::T_numtype Diagonal;

  .. function:: explicit Tridiagonal(const Diagonal& zero =empty, size_t k =0, const Diagonal& minus =empty, const Diagonal& plus =empty, mpl::int_<RANK> one=_1_)

    This is the principal way to create an object of this class, which can be used for ``RANK=1`` only, as ensured by the trailing dummy argument. This creates an object corresponding to the elementary operator

    .. math::
      :label: ElemTridiag

      H_I^\text{elem}(0)=\sum_{n=0}^{N-1}\alpha^{(0)}_n\ket{n}\bra{n}+\sum_{n=0}^{N-1-K}\lp\alpha^{(-)}_n\ket{n+K}\bra{n}+\alpha^{(+)}_n\ket{n}\bra{n+K}\rp

    The arguments ``zero``, ``minus``, ``plus``, and ``k`` correspond respectively to :math:`\alpha^{(0)}`, :math:`\alpha^{(-)}`, :math:`\alpha^{(+)}`, and :math:`K`


  .. function:: Tridiagonal(const Tridiagonal& tridiag)

    Copy constructor with deep copy semantics.

  .. function:: Tridiagonal(const Tridiagonal<RANK2>& tridiag1, const Tridiagonal<RANK__MI__RANK2>& tridiag2)

    ``template <int RANK2>``

    Constructing the object as the direct product of ``tridiag1`` and ``tridiag2``.

  .. function:: void apply(const StateVectorLow& psi, StateVectorLow& dpsidt) const

    ``template <int>`` (A dummy template parameter, cf. multiplication with Sigma.)



.. function:: const Tridiagonal<RANK> quantumoperator::furnishWithFreqs(const Tridiagonal<RANK>& tridiag, const Diagonal& mainDiagonal)

.. note::

  A serious limitation of :class:`~quantumoperator::Tridiagonal` is that the composition of two such operators does not in general yield one of the same form. This is one of the reasons why we are planning to deprecate :class:`~quantumoperator::Tridiagonal` in favour of a much more general form 

  .. math::

    H_I(t)=\bigotimes_{m=0}^{M-1}\sum_{i_m\in\mathbb{K}_m}\sum_{n_m=0}^{N_m-1-i_m}e^{\delta^{i_m}_{m,n}t}\alpha^{i_m}_{m,n}\ket{n_m+i_m}\bra{n_m},

  with :math:`\mathbb{K}_m=\left\{K_m^{(0)},K_m^{(1)},\dots\right\}` an arbitrary set, and :math:`\delta_{m,n}^{(i_m)}=i\lp\omega_{m,n+i_m}-\omega_{m,n}\rp`.

======
Sigma
======
