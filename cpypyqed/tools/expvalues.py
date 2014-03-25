"""
This module provides classes for working with expectation values.

The relevant classes are:
    * :class:`ExpectationValueTrajectory`
    * :class:`ExpectationValueCollection`

The :class:`ExpectationValueTrajectory` stores **one** expectation value at
different points of time while the :class:`ExpectationValueCollection` stores
**several** expectation values at different points of time.
"""

import numpy
import utils

class ExpectationValueTrajectory(numpy.ndarray):
    r"""
    A class representing one expectation value at different points of time.

    *Usage*
        >>> T = numpy.linspace(0,10)
        >>> X = numpy.sin(T)
        >>> ev = ExpectationValueTrajectory(X, T, "<x>")

    *Arguments*
        * *data*
            A 1d array or similar data structure representing an expectation
            value at different points of time.

        * *time* (optional)
            A 1d array or list specifying the points of time.

        * *title* (optional)
            The name of the expectation value. (Can be any unicode string)

        * Any other argument that a numpy array can use for creation. E.g.
          ``copy = False`` can be used so that the ExpectationValueTrajecory
          shares the data storage with the given numpy array.
    """
    def __new__(cls, data, time=None, title=None, **kwargs):
        array = numpy.array(data, **kwargs).view(cls)
        if time is not None:
            array.time = time
        else:
            array.time = getattr(data, "time", None)
        if title is not None:
            array.title = title
        else:
            array.title = getattr(data, "title", None)
        return array

    def __array_finalize__(self, obj):
        self.time = getattr(obj, "time", None)
        self.title = getattr(obj, "title", None)

    def __str__(self):
        clsname = self.__class__.__name__
        if self.title is None:
            title="?"
        else:
            title = self.title
        return "%s('%s')" % (clsname, title)


class ExpectationValueCollection(numpy.ndarray):
    r"""
    A class representing several expectation values at different points of time.

    *Usage*
        >>> import numpy as np
        >>> T = np.linspace(0,10)
        >>> X = np.sin(T)
        >>> Y = np.cos(T)
        >>> ev = ExpectationValueTrajectoryCollection((X,Y), T, ("<x>", "<y>"))

    *Arguments*
        * *data*
            A list of 1D arrays or some similar data structure.

        * *time* (optional)
            A 1d array or list specifying the points of time.

        * *titles* (optional)
            A list of names for the expectation values. (The names can be any
            unicode string)

        * *subsystems* (optional)
            Dictionary specifying subsystems. E.g. ``{"Mode1" : (1,3)}``.

        * Any other argument that a numpy array can use for creation. E.g.
          ``copy = False`` can be used so that the
          ExpectationValueCollection shares the data storage with the given
          numpy array.
    """
    def __new__(cls, data, time=None, titles=None, subsystems=None, **kwargs):
        if isinstance(data, ExpectationValueTrajectory):
            array = numpy.asarray(data).reshape((1,-1)).view(cls)
            array.evtrajectories = (data,)
        else:
            array = numpy.array(data, **kwargs).view(cls)
            if titles is None:
                titles = []
            else:
                titles = list(titles)
            titles = titles + [None]*(len(data) - len(titles))
            traj = [None]*len(data)
            for i, col in enumerate(data):
                traj[i] = ExpectationValueTrajectory(col, time, titles[i])
            array.evtrajectories = tuple(traj)
        if time is not None:
            array.time = time
        else:
            array.time = getattr(data, "time", None)
        array.subsystems = utils.OrderedDict()
        if subsystems is not None:
            for key, value in subsystems.iteritems():
                array.subsystems[key] = cls(
                        array[value[0]:value[1]:], array.time,
                        array.titles[value[0]:value[1]:], copy=True)
        return array

    def __array_finalize__(self, obj):
        self.time = getattr(obj, "time", None)

    def titles(self):
        """
        Return a list of titles of the holded ExpectationValueTrajectories.

        Unknown titles are represented by "?".
        """
        titles = []
        for evt in self.evtrajectories:
            if evt.title is None:
                titles.append("?")
            else:
                titles.append(evt.title)
        return tuple(titles)

    titles = property(titles)

    def __str__(self):
        clsname = self.__class__.__name__
        return "%s('%s')" % (clsname, "', '".join(self.titles))

