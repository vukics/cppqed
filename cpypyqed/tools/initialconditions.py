"""
This module provides convenient methods for creating initial conditions.

At the moment the following initial conditions are implemented:
    * :func:`gaussian`
    * :func:`coherent`
"""

import numpy
import numpy.fft as fft
import quantumstate


def gaussian(x0=0, k0=0, sigma=0.5, fin=6, isItInK=False, cppqed=False):
    r"""
    Generate a StateVector with a normal distribution.

    *Usage*
        >>> sv = gaussian(x0=0.3, k0=4, sigma=0.6, fin=7)
        >>> print sv
        StateVector(128)

    :param double x0: (optional)
            Center in the real space.

    :param double k0: (optional)
            Center in the k-space.

    :param double sigma: (optional)
            Width in the real space :math:`(\sigma = \sqrt{Var(x)})`.

    :param int fin: (optional)
            :math:`2^{fin}` determines the amount of sample points.

    :param bool isItInK: (optional)
            if true, sigma is interpreted as width of the wavepacket in k-space

    :param bool cppqed: (optional)
            C++QED compatibility flag, if set to true then x0 (but not sigma) is expected
            in units of Pi. Defaults to False.

    :returns sv:
            A :class:`.quantumstate.StateVector` representing this
            gaussian wave packet in the k-space.

    The generated StateVector is normalized and given in the k-space. It
    is the discretized and Fourier transformed of the following expression:

        .. math::

            \Psi(x) = \frac {1} {\sqrt[4]{2 \pi}} * \frac {1} {\sqrt{\Delta x}}
                            e^{-\frac {x^2} {4*{\Delta x}^2}}*e^{i k_0x}

    The discretization is given by:

        .. math::

            \Psi_j = \sqrt{\frac{L}{N}} \Psi(x_j),\qquad x_j=-\pi,-\pi+dx,...,\pi-dx
    """

    N = 2**fin
    if cppqed: x0=x0*numpy.pi
    if isItInK:
        L = N
        offset1=k0
        offset2=-x0
    else:
        L = 2*numpy.pi
        offset1=x0
        offset2=k0
    if 6.*sigma > L:
        print "Warning: Sigma is maybe too big."
    dx = L/float(N)
    if sigma < dx:
        print "Warning: Sigma is maybe too small."
    array = numpy.linspace(-L/2., L/2., N, endpoint=False)
    Norm = numpy.sqrt(L/N)/(2*numpy.pi)**(1./4)/numpy.sqrt(sigma)
    Psi = quantumstate.StateVector(Norm*numpy.exp(-(array-offset1)**2/(4*sigma**2))*numpy.exp(1j*array*offset2))
    if not isItInK: Psi=Psi.fft()
    return Psi.normalize()

def coherent(alpha=2, N=20):
    r"""
    Generate a coherent StateVector in the Fock space.

    *Usage*
        >>> sv = coherent(alpha=2, N=20)
        >>> print sv
        StateVector(20)

    :param alpha: (optional)
            A complex number specifying the coherent state. (Default is 2)

    :param N: (optional)
            A number determining the dimension of the Fock space. (Default is
            20)

    :returns sv:
            A :class:`.quantumstate.StateVector`.

    The coherent state is given by the formula:

        .. math::

            |\alpha\rangle = e^{-\frac {|\alpha|^2} {2}} \sum_{n=0}^{N}
                                \frac {\alpha^n} {\sqrt{n!}} |n\rangle

    Calculation is done using the recursive formula:

        .. math::

            a_0 = e^{- \frac {|\alpha|^2} {2}}

        .. math::

            a_n = a_{n-1} * \frac {\alpha} {\sqrt n}
    """
    x = numpy.empty(N)
    x[0] = 1
    for n in range(1,N):
        x[n] = x[n-1]*alpha/numpy.sqrt(n)
    x = numpy.exp(-numpy.abs(alpha)**2/2.)*x
    return quantumstate.StateVector(x)

def fock(dim,i):
    r"""
    Generate a Fock space basis vector.

    *Usage*
        >>> sv = fock(dim=8, i=4)
        >>> print sv
        StateVector(8)

    :param int dim:
            Dimension of the Fock space.

    :param int i:
            Genereate the i-th basis vector

    :returns sv:
            A :class:`.quantumstate.StateVector`.
    """
    psi=numpy.zeros(dim)
    psi[i]=1.
    return quantumstate.StateVector(psi)

def momentum_eigen(k,fin=6):
    r"""
    Generate a momentum eigenstate with momentum p in a space with finesse fin.

    *Usage*
        >>> sv = momentum_eigen(0, fin=4)
        >>> print sv
        StateVector(16)

    :param int k:
            Wavenumber of the resulting state.

    :param int fin: (optional)
            Finesse of the resulting state (default 4). 

    :returns sv:
            A :class:`.quantumstate.StateVector`.
    """

    num=2**fin
    return fock(num,k+num/2)