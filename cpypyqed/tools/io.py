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

def _open_file(filename):
  if filename.endswith("bz2"):
    try:
      import bz2
    except ImportError:
      raise IOError("{filename} is compressed, but the bz2 python module could not be loaded.".format(filename=filename))
    f = bz2.BZ2File(filename)
  else:
    f = open(filename)
  return f

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
  Returns a tuple with a list of state vectors and a list of times.

  *Usage*
      >>> svs,times = load_statevector("ring.sv")

  :param filename:
          Path to the C++QED state vector file that should be loaded.
  :rtype: tuple
  """
  svs,times = read(filename)[1:]
  return [quantumstate.StateVector(sv[0],time=sv[1]) for sv in zip(svs,times)],times

def load_densityoperator(filename):
  """
  Load a C++QED state vector file from the given location.
  Returns a tuple with a list of densityoperators and a list of times.

  *Usage*
      >>> svs,times = load_statevector("ring.sv")

  :param filename:
          Path to the C++QED state vector file that should be loaded.
  :rtype: tuple
  """
  svs,times = read(filename)[1:]
  return [quantumstate.DensityOperator(sv[0],time=sv[1]) for sv in zip(svs,times)],times

def evs_header(filename):
  """
  Load the header of a C++QED trajectory file from the given location. The file
  can be bzip2-compressed if the bz2 module is available.

  *Usage*
      >>> header = evs_header("ring.sv")

  :param filename:
          Path to the C++QED state vector file that should be loaded.
          :class:`IOError` is raised if the required `cpypyqed.io` module is not
          available.
  :rtype: string
  """

  f=_open_file(filename)

  result=""
  notDone=True
  while notDone:
    try:
      line = f.next()
    except StopIteration:
      break
    if line.strip(" \n") and not line.startswith("#"):
      notDone=False
    else:
      result+=line
  return result

def load_jumps(filename):
  """
  Load the jump times from a  C++QED trajectory file. The file
  can be bzip2-compressed if the bz2 module is available. Returns a
  two-dimensional :class:`np.ndarray` where the first column is the time
  and the second column is the dissipation channel.

  *Usage*
      >>> jumps = load_jumps("ring.out.1001.bz2")
      >>> # every jump of decay channel 1
      >>> print(jumps[jumps[:,1]==1,0])
      array([ 3.24715])

  :param filename:
          Path to the C++QED trajectory file from which the jumps should be loaded.
  :rtype: np.ndarray
  """

  f=_open_file(filename)

  for num, line in enumerate(f, 1):
    if '# MCWF Trajectory:' in line:
      f.close()
      return np.loadtxt(filename,skiprows=num,comments='',usecols=(1,2))
  return None