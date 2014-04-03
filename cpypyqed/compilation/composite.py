# Copyright Raimar Sandner 2012-2014. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
# -*- coding: utf-8 -*-

r"""This module is the wrapper for the C++QED :core2:`composite` bundle.

Typically the enduser only has to call the convenience maker function :func:`makeComposite`.
"""

from ondemand import OnDemand
import os

def makeComposite(acts):
    r"""makeComposite(acts) -> Composite<VA> :
    Wrapper of the C++QED function :core2:`composite::make` with a modified signature better suited for Python.

    This function takes a dictionary where the keys describe an :class:`Act` object and
    values are :class:`InteractionN <.structure.Interaction2>`-objects (with N the rank of
    the interaction). The key has to be a N-tuple of integers corresponding to the involved subsystems.
    For example, the key value pair :samp:`{(0,2):poc}` will be a :samp:`Act<0,2>(poc)`, where :samp:`poc`
    could be an instance of :class:`.ParticleOrthogonalToCavity` and the interaction acts on the subsystems 0 and 2.

    Example usage (where :samp:`poc` is a :class:`.ParticleOrthogonalToCavity`-object, subsytem 0 is a mode and subsystems
    1 and 2 are particles)::

        >>> c = makeComposite({(0,2):poc,(0,1):poc})

    :param dict acts: A dictionary holding the subsystem tuples as keys and the interaction objects as values.
    :returns: A wrapped C++QED :core2:`Composite`-object, compiled on demand if necessary.
    :rtype: :samp:`Composite<VA>`
    """
    return Composite(acts).maker()

class Composite(OnDemand):
    r"""A :class:`.OnDemand`-class wrapping the C++QED class :core2:`Composite`.

    The constructor takes a dictionary where the keys describe an :class:`Act` object and
    values are :class:`InteractionN <.structure.Interaction2>`-objects (with N the rank of
    the interaction). The key has to be a N-tuple of integers corresponding to the involved subsystems.
    For example, the key value pair `{(0,2):poc}` will be a `Act<0,2>(poc)`, where `poc` could be an instance of
    :class:`.ParticleOrthogonalToCavity` and the interaction acts on the subsystems 0 and 2.

    Note that the constructor will only yield an :class:`.OnDemand`-object, not a real C++QED composite object.
    However, calling the :meth:`.OnDemand.maker`-method on the resulting object will return a true
    composite object, while instantiating and compiling the class template as needed. The convinience maker function
    :func:`makeComposite` has the same signature as the constructor here and automatically calls the
    underlying :meth:`.OnDemand.maker` function, thus returning a true composite object which can be
    passed to :func:`.evolve`.

    :param dict acts: A dictionary holding the subsystem tuples as keys and the interaction objects as values.
    """

    def __init__(self, acts):
        self.acts = acts
        self.actslist = sorted(self.acts.keys())
        OnDemand.__init__(self, "Composite", self._make_id(), makerfunction="compositeMake")

    def _make_from_actlist(self,prefix,postfix,enum_large,enum_small):
        return enum_large.join([prefix+enum_small.join(map(str,i))+postfix for i in self.actslist])

    def _make_id(self):
        return self._make_from_actlist('', '', '__', '_')

    def _make_actlist(self):
        return self._make_from_actlist('Act<', '>', ',', ',')

    def _make_const_act_arguments(self):
        return self._make_from_actlist('const Act<', '>&', ',', ',')

    def _make_custodian_and_wards(self, depth=1):
        if depth>len(self.actslist): return ""
        if depth == 1:
          return "with_custodian_and_ward_postcall<0,1{0}>()".format(self._make_custodian_and_wards(depth+1))
        else:
          return ",with_custodian_and_ward_postcall<0,{0}{1}>".format(depth,self._make_custodian_and_wards(depth+1))

    def generate_source(self,builddir):
        r"""Overriding :meth:`.OnDemand.generate_source`

        This takes care of creating the proper source file to compile a boost.Python module, which holds
        the correct C++QED :core2:`Composite` class.
        """
        with open(os.path.join(self.cpypyqeddir,"CompositeTemplate.cc")) as f:
            composite_template = f.read().format(actlist=self._make_actlist(),
                                                 const_act_arguments=self._make_const_act_arguments(),
                                                 classname=self.classname,
                                                 custodian_and_wards=self._make_custodian_and_wards(),
                                                 modulename=self.modulename)
        with open(os.path.join(builddir,self.sourcefile),"w") as f:
            f.write(composite_template)

    def maker(self):
        r"""Overriding :meth:`.OnDemand.maker`

        This function compiles and creates the needed C++QED :core2:`Act`-objects on demand and passes them to
        :meth:`.OnDemand.maker`, where they are passed into :core2:`composite::make` to create the right
        :core2:`Composite` object.
        """
        args=[]
        for i in self.actslist:
            args.append(Act(*i).maker(self.acts[i]))
        return OnDemand.maker(self,*args)

class Act(OnDemand):
    r"""A :class:`.OnDemand`-class wrapping the C++QED class :core2:`Act`.

    The template parameters denoting the involved subsystems are given as constructor arguments to this class.
    The interaction which this object should hold has to be passed in as an argument to :meth:`.OnDemand.maker`.
    """
    def __init__(self,*args):
        self.subsys=map(str,args)
        self.act="Act<"+','.join(self.subsys)+">"
        OnDemand.__init__(self, "Act", '_'.join(self.subsys))

    def generate_source(self,builddir):
        with open(os.path.join(self.cpypyqeddir,"ActTemplate.cc")) as f:
            act_template = f.read().format(act=self.act,modulename=self.modulename,classname=self.classname)
        with open(os.path.join(builddir,self.sourcefile),"w") as f:
            f.write(act_template)