"""
This module provides convenient methods for creating operators.

At the moment the following initial conditions are implemented:
    * :func:`gaussian`
    * :func:`coherent`
    * :func:`fock`
    * :func:`momentum_eigen`
"""

import numpy as np


def create(dim):
  r"""
  Generate the creation operator for a mode with `N` states

  :param int dim:
    Dimension of the space.
  :returns sv:
    A :class:`numpy.array`.
  """

  assert dim>0
  return np.diag([np.sqrt(n) for n in range(1,dim)],-1)

def destroy(dim):
  r"""
  Generate the destruction operator for a mode with `N` states

  :param int dim:
    Dimension of the space.
  :returns sv:
    A :class:`numpy.array`.
  """

  assert dim>0
  return np.diag([np.sqrt(n) for n in range(1,dim)],1)