"""
This module provides an implementation for state vectors and density operators.
These two classes are not direct wrappings of their C++QED counterparts to make them
behave more pythonic. In fact, both these classes are derived from `numpy.ndarray`.

The relevant classes are:
    * :class:`StateVector`
    * :class:`DensityOperator`
    * :class:`StateVectorTrajectory`
"""

import numpy
import expvalues
try:
    set()
except NameError:
    from sets import Set as set

class QuantumState(numpy.ndarray):
  r"""This is the base class for :class:`StateVector` and :class:`DensityOperator`.
  It inherits from `numpy.ndarray`.

  :param numpy.ndarray data: Anything a numpy.ndarray can beconstructed from.
  :param double time: (optional) A number defining the point of time when this state vector was
      reached. (Default is 0)

  Any other argument that a numpy array takes. E.g. ``copy=False`` can
  be used so that the QuantumState shares the data storage with the given numpy array.

  Most useful is maybe the tensor product '**' which lets you easily calculate
  state vectors for combined systems.
  """
  def __new__(cls, data, time=None,**kwargs):
    array = numpy.array(data, **kwargs)
    array = numpy.asarray(array,dtype=numpy.complex128).view(cls)
    if time is not None:
        array.time = time
    elif hasattr(data, "time"):
        array.time = data.time
    return array

  def __array_finalize__(self, obj):
    self.time = getattr(obj, "time", 0)

  def __str__(self):
    return numpy.ndarray.__str__(numpy.asarray(self))

  def __repr__(self):
    return numpy.ndarray.__repr__(numpy.asarray(self))

  def __unicode__(self):
    return numpy.ndarray.__unicode__(numpy.asarray(self))

  def _fft_helper(self, axes=None, inverse=False):
    r"""
    Helper function for fft and ifft which performs the actual transformation.
    """
    f = numpy.fft
    norm = numpy.sqrt(numpy.prod(numpy.array(self.shape)[axes]))
    if inverse:
      transform=f.ifftn
    else:
      transform=f.fftn
      norm=1/norm
    array = f.fftshift(transform(f.ifftshift(self, axes=axes), axes=axes), axes=axes)*norm
    return type(self)(array, time=self.time)

class DensityOperator(QuantumState):
  r"""
  A class representing a quantum mechanical density operator. :core2:`quantumdata::DensityOperator` is
  automatically converted to this, but it is not a one to one wrapper.
  """
  def __new__(cls, data, time=None, **kwargs):
    array = QuantumState.__new__(cls, data, time, **kwargs)
    return array

  def __array_finalize__(self, obj):
    QuantumState.__array_finalize__(self,obj)
    ndim=len(self.shape)
    if not ndim % 2 == 0:
      raise ValueError("The number of dimension must be even for a density operator.")
    self.dimensions = self.shape[:ndim/2]
    if not self.shape[:ndim/2]==self.shape[ndim/2:]:
      raise ValueError("The created object is not a valid density operator, dimensions mismatch.")

  def fft(self, subsystems=None):
    r"""
    Return a DensityOperator where the given subsystems are Fourier transformed.
    This is the transformation position space -> momentum space.

    :param subsystems:
            (optional)
            Sequence of ints, subsystems over which the fft is done. (Default is all)
    """
    ndim=len(self.shape)
    if subsystems is None: subsystems = range(ndim/2)
    return self._fft_helper(axes=subsystems,inverse=False)._fft_helper(axes=ndim/2+numpy.array(subsystems),inverse=True)

  def ifft(self, subsystems=None):
    r"""
    Return a DensityOperator where the given subsystems are inversely Fourier transformed.
    This is the transformation momentum space -> position space.

    :param subsystems:
            (optional)
            Sequence of ints, subsystems over which the ifft is done. (Default is all)
    """
    ndim=len(self.shape)
    if subsystems is None: subsystems = range(ndim/2)
    return self._fft_helper(axes=subsystems,inverse=True)._fft_helper(axes=ndim/2+numpy.array(subsystems),inverse=False)


class StateVector(QuantumState):
    r"""
    A class representing a quantum mechanical state. :core2:`quantumdata::StateVector` is
    automatically converted to this, but it is not a one to one wrapper.

    *Usage*
        >>> sv = StateVector((1, 3, 7, 2), time=0.2, norm=True)
        >>> sv = StateVector(numpy.arange(12).reshape(3,4))
        >>> print sv
        StateVector(3 x 4)
        >>> print repr(sv)
        StateVector([[ 0,  1,  2,  3],
               [ 4,  5,  6,  7],
               [ 8,  9, 10, 11]])

    :param data:
            Anything that can be used to create a numpy array, e.g. a nested
            tuple or another numpy array.

    :param double time: (optional)
            A number defining the point of time when this state vector was
            reached. (Default is 0)

    :param bool norm: (optional)
            If set True the StateVector will be automatically normalized.
            (Default is False)

    :param \*\*kwargs: Any other argument that a numpy array takes. E.g. ``copy=False`` can
          be used so that the StateVector shares the data storage with the
          given numpy array.

    Most useful is maybe the tensor product which lets you easily calculate
    state vectors for combined systems::

        >>> sv1 = StateVector((1,2,3))
        >>> sv2 = StateVector((3,4,0), norm=True)
        >>> sv = sv1 ** sv2
        >>> print sv
        StateVector(3 x 3)
        >>> print repr(sv)
        StateVector([[ 0.6,  0.8,  0. ],
               [ 1.2,  1.6,  0. ],
               [ 1.8,  2.4,  0. ]])

    The tensor product is abbreviated by the "**" operator.
    """
    def __new__(cls, data, time=None, norm=False, **kwargs):
        array = QuantumState.__new__(cls, data, time, **kwargs)
        if norm:
            array = normalize(array)
        return array

    def __array_finalize__(self, obj):
        QuantumState.__array_finalize__(self,obj)
        self.dimensions = obj.shape

    def __str__(self):
        clsname = self.__class__.__name__
        return "%s(%s)" % (clsname, " x ".join(map(str, self.dimensions)))

    def norm(self):
        r"""
        Calculate the norm of the StateVector.

        *Usage*
            >>> sv = StateVector((1,2,3,4,5), norm=True)
            >>> print sv.norm()
            1.0
        """
        return norm(self)

    def normalize(self):
        r"""
        Return a normalized StateVector.

        *Usage*
            >>> sv = StateVector((1,2,1,3,1))
            >>> print sv.norm()
            4.0
            >>> nsv = sv.normalize()
            >>> print nsv.norm()
            1.0
        """
        return normalize(self)

    def reduce(self, indices, norm=True):
        r"""
        Return a StateVector where the given indices are reduced.

        *Usage*
            >>> rsv = sv.reduce(1)
            >>> rsv = sv.reduce((1,2))

        :parameter indices:
                An integer or a list of integers specifying over which
                subspaces should be summated.

        :parameter bool norm: (optional)
                If set True the resulting StateVector will be renormalized.

        Reducing means nothing else then summing up over all given indices.
        E.g. a StateVector of rank 4 can be reduced to the first two indices::

            >>> sv1 = StateVector((1,2), norm=True)
            >>> sv2 = StateVector((1,2,3), norm=True)
            >>> sv3 = StateVector((1,2,3,4,5), norm=True)
            >>> sv4 = StateVector((1,2,3,4,5,6), norm=True)
            >>> sv = sv1**sv2**sv3**sv4
            >>> print sv
            StateVector(2 x 3 x 5 x 6)
            >>> print sv.reduce((2,3))
            StateVector(2 x 3)

        This is mathematically equivalent to:

            .. math::

                \Psi_{\alpha \beta} = \frac
                  {\sum_{\gamma \delta} \Psi_{\alpha \beta \gamma \delta}}
                  {\| \sum_{\gamma \delta} \Psi_{\alpha \beta \gamma \delta} \|}

        Reducing is an easy way to find out how subspaces of a high rank
        state vectors behave. Don't use reduced StateVectors for calculating
        expectation values - this will most likely give wrong answers!
        """
        if isinstance(indices, int):
            a = (indices,)
        else:
            a = _sorted_list(indices, True)
        array = self
        if norm:
            for i in a:
                array = array.sum(axis=i).normalize()
        else:
            for i in a:
                array = array.sum(axis=i)
        return array

    def reducesquare(self, indices):
        r"""
        Calculate the reduced Psi-square tensor.

        *Usage*
            >>> sv1 = StateVector((0,1,2,1,0), norm=True)
            >>> sv2 = StateVector((1,0,1), norm=True)
            >>> sv = sv1**sv2
            >>> sqtensor = sv.reducesquare(1)

        :parameter indices:
                An integer or a list of integers specifying over which
                subsystems should be summed up.

        This method calculates the following quantity (simplified for rank 2
        state vectors):

            .. math::

                w_{\alpha_1 \alpha_2} = \sum_\beta \Psi_{\alpha_1 \beta}^*
                                        \Psi_{\alpha_2 \beta}

        Where :math:`\beta` is the reduced index.

        This quantity is useful to calculate expectation values in the
        corresponding subspaces.
        """
        if isinstance(indices, int):
            a = (indices,)
        else:
            a = _sorted_list(indices, True)
        return numpy.tensordot(self, self.conjugate(), (a,a))

    def fft(self, axes=None):
        r"""
        Return a StateVector where the given axes are Fourier transformed.
        This is the transformation position space -> momentum space.

        *Usage*
            >>> sv = StateVector((0,1,1.7,2,1.7,1,0), norm=True)
            >>> print sv.fft()
            StateVector(7)

        :parameter axis: (optional)
                Sequence of ints, axes over which the fft is done. (Default is all)
        """
        return self._fft_helper(axes=axes,inverse=False)

    def ifft(self, axes=None):
        r"""
        Return a StateVector where the given axes are inversely Fourier transformed.
        This is the transformation momentum space -> position space.

        See :func:`StateVector.fft` for details.
        """
        return self._fft_helper(axes=axes,inverse=True)

    def expvalue(self, operator, indices=None, title=None, multi=False):
        r"""
        Calculate the expectation value of the given operator.

        *Usage*
            >>> a = numpy.diag(numpy.ones(3), -1)
            >>> print a
            array([[ 0.,  0.,  0.,  0.],
                   [ 1.,  0.,  0.,  0.],
                   [ 0.,  1.,  0.,  0.],
                   [ 0.,  0.,  1.,  0.]])
            >>> sv = StateVector((1,2,1,2), norm=True)
            >>> print sv.expvalue(a)
            0.6

        :parameter operator:
                A tensor representing an arbitrary operator in the
                basis of the StateVector.

        :parameter indices: (optional)
                Specifies which subsystems should be taken. If None is given
                the whole system is used.

        :parameter bool multi: (optional)
                If multi is True it is assumed that a list of operators is
                given. (Default is False)

        Expectation values for combined systems are calculated in the following
        way (Assuming the operator only acts on first subsystem):

            .. math::

                \langle \Psi | \hat A (k) | \Psi \rangle =
                    \sum_{k_1 k_2} \langle k_1 | \hat A (k) | k_2 \rangle
                    \sum_m \Psi_{k_1 m}^* \Psi_{k_2 m}

        The second sum is exactly what :meth:`reducesquare` does while the
        first expression is the matrix representation of the given operator
        in the same basis as the StateVector.
        """
        if indices is not None:
            A = self.reducesquare(_conjugate_indices(indices, self.ndim))
        else:
            A = self**self.conjugate()
        length = A.ndim
        index = range(0, length, 2) + range(1, length, 2)
        if multi:
            evs = [(A*op.transpose(index)).sum() for op in operator]
            return expvalues.ExpectationValueCollection(evs, self.time, title)
        else:
            return (A*operator.transpose(index)).sum()

    def diagexpvalue(self, operator, indices=None, title=None, multi=False):
        r"""
        Calculate the expectation value for the given diagonal operator.

        *Usage*
            >>> a = numpy.arange(4)
            >>> print a
            array([ 0.,  1.,  2.,  3.])
            >>> sv = StateVector((1,2,1,4), norm=True)
            >>> print sv.diagexpvalue(a)
            2.45454545455

        :parameter operator:
                The diagonal elements of a tensor representing an arbitrary
                diagonal operator in the basis of the StateVector.

        :parameter indices: (optional)
                Specifies which subsystems should be taken. If None is given
                the whole system is used.

        :parameter bool multi: (optional)
                If multi is True it is assumed that a list of operators is
                given. (Default is False)

        Expectation values for combined systems are calculated in the following
        way (Assuming the operator only acts on first subsystem):

            .. math::

                \langle \Psi | \hat A (k) | \Psi \rangle =
                    \sum_k \langle k | \hat A (k) | k \rangle
                    \sum_m \Psi_{k m}^* \Psi_{k m}

        Other than the general :meth:`expvalue` method :meth:`diagexpvalue`
        only works for diagonal operators and only needs the diagonal elements
        of the matrix representation.
        """
        if isinstance(indices, int):
            indices = (indices,)
        A = self*self.conjugate()
        if indices is not None:
            indices = _sorted_list(_conjugate_indices(indices, self.ndim),
                                   True)
            for index in indices:
                A = A.sum(index)
        if multi:
            evs = [(A*op).sum() for op in operator]
            return expvalues.ExpectationValueCollection(evs, self.time, title)
        else:
            return (A*operator).sum()

    def dyad(self):
      r"""
      Calculate the dyadic product with itself.

      :returns: The DensityOperator :math:`| \Psi \rangle \langle \Psi |`.
      """
      return DensityOperator(self**numpy.conjugate(self), time=self.time)


    def outer(self, array):
        r"""
        Return the outer product between this and the given StateVector.

        *Usage*
            >>> sv = StateVector((0,1,2), norm=True)
            >>> print repr(sv.outer(StateVector((3,4), norm=True)))
            StateVector([[ 0.        ,  0.        ],
                   [ 0.26832816,  0.35777088],
                   [ 0.53665631,  0.71554175]])
            >>> print sv.outer((3,4)) # Not normalized!
            StateVector([[ 0.        ,  0.        ],
                   [ 1.34164079,  1.78885438],
                   [ 2.68328157,  3.57770876]])

        :param array:
                Some kind of array (E.g. StateVector, numpy.array, list, ...).

        As abbreviation ``sv1**sv2`` can be written instead of
        ``sv1.outer(sv2)``.
        """
        return StateVector(numpy.multiply.outer(self, array))

    def __pow__(self, other):
      if type(other)==type(self):
        return self.outer(other)
      elif type(other)==int:
        return numpy.ndarray.__pow__(self,other)


class StateVectorTrajectory(numpy.ndarray):
    """
    A class holding StateVectors for different points of time.

    :param data:
            Some nested structure which holds state vector like arrays for
            different points of time.

    :param double time: (optional)
            An array which specifies the point of time for every state
            vector. This array must have as many entries as there are state
            vectors.

    :param \*\*kwargs:
            Any other argument that a numpy array takes. E.g. ``copy=False`` can
            be used so that the StateVectorTrajectory shares the data storage
            with the given numpy array.

    Most methods are simple mapped to all single StateVectors. For more
    documentation regarding these methods look into the docstrings of the
    corresponding :class:`StateVector` methods.
    """
    def __new__(cls, data, time=None, **kwargs):
        array = numpy.array(data, **kwargs)
        array = numpy.asarray(array).view(cls)
        if time is None:
            array.time = numpy.array([sv.time for sv in data])
        else:
            array.time = time
        svs = [None]*array.shape[0]
        for i, entry in enumerate(array):
            svs[i] = StateVector(entry, time=array.time[i], copy=True)
        array.statevectors = svs
        return array

    def __array_finalize__(self, obj):
        self.dimensions = obj.shape[1:]

    def map(self, func, svt=True):
        """
        Apply the given function to every single StateVector.

        *Usage*
            >>> norm = svt.map(lambda sv:sv.norm())

        :paramter func:
                Function that takes a StateVector as argument.

        :parameter bool svt: (optional)
                If svt is True, the return value will be an instance of
                StateVectorTrajectory.
        """
        svs = [None]*self.shape[0]
        for i, sv in enumerate(self.statevectors):
            svs[i] = func(sv)
        if svt:
            return StateVectorTrajectory(svs)
        else:
            return svs

    def norm(self):
        """
        Return a list of norms for every single StateVector.

        See also: :meth:`StateVector.norm`
        """
        return self.map(lambda sv:sv.norm(), False)

    def normalize(self):
        """
        Return a StateVectorTrajectory where all StateVectors are normalized.

        See also: :meth:`StateVector.normalize`
        """
        return self.map(lambda sv:sv.normalize())

    def reduce(self, indices, norm=True):
        """
        Return a StateVectorTrajectory where all StateVectors are reduced.

        See also: :meth:`StateVector.reduce`
        """
        return self.map(lambda sv:sv.reduce(indices, norm=norm))

    def fft(self, axis=0):
        """
        Return a StateVectorTrajectory whith Fourier transformed StateVectors.

        See also: :meth:`StateVector.fft`
        """
        return self.map(lambda sv:sv.fft(axis))

    def expvalue(self, operator, indices=None, multi=False, titles=None):
        """
        Calculate the expectation value of the operator for all StateVectors.

        :returns:
                An :class:`.expvalues.ExpectationValuesTrajectory`
                instance.

        See also: :meth:`StateVector.expvalue`
        """
        evs = self.map(lambda sv:sv.expvalue(operator, indices, multi=multi),
                       False)
        if not multi:
            return expvalues.ExpectationValueTrajectory(evs, self.time, titles)
        evs = numpy.array(evs).swapaxes(0,1)
        return expvalues.ExpectationValueCollection(
                            evs, self.time, titles, copy=False)

    def diagexpvalue(self, operator, indices=None, multi=False, titles=None):
        """
        Calculate the expectation value of the diagonal operator for all SVs.

        :returns:
                An :class:`.expvalues.ExpectationValuesTrajectory`
                instance.

        See also: :meth:`StateVector.diagexpvalue`
        """
        evs = self.map(lambda sv:sv.diagexpvalue(operator, indices,
                            multi=multi), False)
        if not multi:
            return expvalues.ExpectationValueTrajectory(evs, self.time, titles)
        evs = numpy.array(evs).swapaxes(0,1)
        return expvalues.ExpectationValueCollection(
                            evs, self.time, titles, copy=False)


def norm(array):
    """
    Return the norm of the array.
    """
    return numpy.sqrt((array*array.conj()).sum())

def normalize(array):
    """
    Return a normalized array.
    """
    return array/norm(array)

def adjust(array, length):
    """
    Adjust the dimensionality of a 1D array.
    """
    import scipy.interpolate
    X_old = numpy.linspace(0,1,len(array))
    f = scipy.interpolate.interp1d(X_old, array)
    X_new = numpy.linspace(0,1,length)
    return StateVector(f(X_new))

def _dim2str(dimensions):
    """
    Return the corresponding dimension string for the given nested tuple.
    """
    dims = []
    for d in dimensions:
        dims.append("(%s,%s)" % d)
    return " x ".join(dims)

def _conjugate_indices(indices, ndim):
    """
    Return all numbers from 0 to ndim which are not in indices.
    """
    if isinstance(indices, int):
        indices = (indices,)
    return set(range(ndim)).difference(indices)

def _sorted_list(iterable, reverse=False):
    """
    Transform an iterable to a sorted list.
    """
    a = list(iterable)
    a.sort()
    if reverse:
        a.reverse()
    return a
