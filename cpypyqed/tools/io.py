"""
This module provides convenient methods for accessing C++QED data (expectation values and state files)
"""

import numpy as np
import expvalues
import quantumstate
from ..config import build_type
if build_type=="debug":
  from ..io_d import read
else:
  from ..io import read

def load_evs(filename, maxevs=None):
  r"""
  Load a trajectory file into an :class:`.expvalues.ExpectationValueCollection`.

  :param filename: The name of the trajectory file.
  :type filename: str
  :param maxevs: Maximal nunmber of columns to read in.
  :type maxevs: int
  :rtype: expvalues.ExpectationValueCollection
  """

  evs = np.loadtxt(filename).swapaxes(0,1)
  time = evs[0,:]
  titles = []
  evstraj = expvalues.ExpectationValueCollection(evs, time=time, copy=False)
  return evstraj


def load_statevector(filename):
  """
  Load a C++QED state vector file from the given location.

  *Usage*
      >>> sv = load_statevector("ring.sv")

  :param filename:
          Path to the C++QED state vector file that should be loaded.
          :class:`IOError` is raised if the required `cpypyqed.io` module is not
          available.
  :rtype: quantumstate.StateVectorTrajectory
  """
  states, times = read(filename)[1:]
  return quantumstate.StateVectorTrajectory(states,times)

