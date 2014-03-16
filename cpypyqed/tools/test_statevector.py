import unittest
import numpy
import quantumstate
import initialconditions
import expvalues

eps = 1e-15

class StateVectorTestCase(unittest.TestCase):
    def test_creation(self):
        sv = quantumstate.StateVector((1,2,3,4), 1)
        self.assertEqual(sv.dimensions, (4,))
        self.assertEqual(sv.time, 1)
        sv = quantumstate.StateVector(numpy.empty((12,10)))
        self.assertEqual(sv.dimensions, (12, 10))
        self.assertEqual(sv.time, 0)

    def test_normalize(self):
        sv = quantumstate.StateVector((5,1,6), norm=True)
        N = numpy.sqrt((sv*sv.conjugate()).sum())
        self.assert_(1-eps<N<1+eps)
        self.assert_(1-eps<sv.norm()<1+eps)

    def test_reduce(self):
        sv1 = quantumstate.StateVector((2,1,1), norm=True)
        sv2 = quantumstate.StateVector((3,2,6,7), norm=True)
        sv3 = quantumstate.StateVector((1,0,3,4), norm=True)
        sv = sv1 ** sv2 ** sv3
        self.assert_(((sv.reduce((1,2))-sv1)<eps).all())
        self.assert_(((sv.reduce((0,2))-sv2)<eps).all())
        self.assert_(((sv.reduce((0,1))-sv3)<eps).all())
        self.assert_(((sv.reduce(2)-(sv1**sv2))<eps).all())

    def test_reducesquare(self):
        sv1 = quantumstate.StateVector((2,1,1), norm=True)
        sv2 = quantumstate.StateVector((3,2,6,7), norm=True)
        sv3 = quantumstate.StateVector((1,0,3,4), norm=True)
        sv = sv1 ** sv2 ** sv3
        rsv = sv.reducesquare((2,1))
        self.assertEqual(rsv.shape, (3,3))
        rsv = sv.reducesquare(2)
        self.assertEqual(rsv.shape, (3,4,3,4))

    def test_fft(self):
        sv_k = initialconditions.gaussian(0.5,2,0.3)
        sv_x = sv_k.fft()
        X = numpy.linspace(-numpy.pi, numpy.pi, 64, endpoint=False)
        ev_x = sv_x.diagexpvalue(X)*2*numpy.pi/64
        self.assert_(numpy.abs(ev_x-0.5)<eps)

    def test_expvalue(self):
        sv1 = quantumstate.StateVector((1,1,1), norm=True)
        sv2 = quantumstate.StateVector((1,2,3,4), norm=True)
        sv = sv1**sv2
        X = numpy.diag((1,2,3))
        ev1_X = sv1.expvalue(X)
        self.assert_(not isinstance(ev1_X, numpy.ndarray))
        self.assert_((ev1_X-2)<eps)
        ev_X = sv.expvalue(X, indices=(0))
        self.assert_(not isinstance(ev_X, numpy.ndarray))
        self.assert_((ev_X-2)<eps)
        Y = numpy.diag((0,1,2))
        ev = sv1.expvalue((X,Y), multi=True)
        self.assert_(isinstance(ev, expvalues.ExpectationValueCollection))
        self.assertEqual(ev.shape, (2,))

    def test_diagexpvalue(self):
        sv1 = quantumstate.StateVector((1,1,1), norm=True)
        sv2 = quantumstate.StateVector((1,2,3,4), norm=True)
        sv = sv1**sv2
        X = (1,2,3)
        ev1_X = sv1.diagexpvalue(X)
        self.assert_(not isinstance(ev1_X, numpy.ndarray))
        self.assert_((ev1_X-2)<eps)
        ev_X = sv.diagexpvalue(X, indices=(0))
        self.assert_(not isinstance(ev_X, numpy.ndarray))
        self.assert_((ev_X-2)<eps)
        Y = (0,1,2)
        ev = sv1.diagexpvalue((X,Y), multi=True)
        self.assert_(isinstance(ev, expvalues.ExpectationValueCollection))
        self.assertEqual(ev.shape, (2,))

    def test_outer(self):
        sv1 = quantumstate.StateVector((1,2), norm=False)
        sv2 = quantumstate.StateVector((3,4), norm=False)
        sv = quantumstate.StateVector(((3,4),(6,8)))
        self.assert_((sv1**sv2==sv).all())

    def test_adjust(self):
        sv = quantumstate.StateVector(numpy.sin(numpy.linspace(0,10)))
        sv_new = quantumstate.adjust(sv, 10)
        self.assertEqual(len(sv_new), 10)


class StateVectorTrajectoryTestCase(unittest.TestCase):
    def sv(self, t, points=16):
        X = numpy.linspace(-numpy.pi, numpy.pi, points)
        return quantumstate.StateVector(numpy.sin(X)*numpy.exp(-t), time=t)

    def test_creation(self):
        t = numpy.linspace(0,5,20)
        sv = quantumstate.StateVectorTrajectory(map(self.sv, t))
        self.assertEqual(sv.shape, (20,16))

    def test_diagexpvalue(self):
        t = numpy.linspace(0,5,20)
        sv = quantumstate.StateVectorTrajectory(map(self.sv, t))
        X = numpy.arange(16)
        ev = sv.diagexpvalue(X)
        self.assert_(isinstance(ev, expvalues.ExpectationValueTrajectory))
        self.assertEqual(ev.shape, (20,))

        ev = sv.diagexpvalue((X,X+2), multi=True)
        self.assert_(isinstance(ev, expvalues.ExpectationValueCollection))
        self.assertEqual(ev.shape, (2,20))


def suite():
    load = unittest.defaultTestLoader.loadTestsFromTestCase
    suite = unittest.TestSuite([
            load(StateVectorTestCase),
            load(StateVectorTrajectoryTestCase)
            ])
    return suite


if __name__ == "__main__":
    unittest.main()
