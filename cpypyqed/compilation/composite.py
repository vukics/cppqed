# -*- coding: utf-8 -*-

from ondemand import OnDemand
import os

def makeComposite(acts):
    return Composite(acts).maker()

class Composite(OnDemand):

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

    def _make_custodian_and_wards(self):
        return ','.join(['with_custodian_and_ward_postcall<0,'+str(i)+'>()' for i in range(1,len(self.actslist)+1)])

    def generate_source(self,builddir):
        with open(os.path.join(self.cpypyqeddir,"CompositeTemplate.cc")) as f:
            composite_template = f.read().format(actlist=self._make_actlist(),
                                                 const_act_arguments=self._make_const_act_arguments(),
                                                 classname=self.classname,
                                                 custodian_and_wards=self._make_custodian_and_wards(),
                                                 modulename=self.modulename)
        with open(os.path.join(builddir,self.sourcefile),"w") as f:
            f.write(composite_template)

    def maker(self):
        args=[]
        for i in self.actslist:
            args.append(Act(*i).maker(self.acts[i]))
        return OnDemand.maker(self,*args)

class Act(OnDemand):

    def __init__(self,*args):
        self.subsys=map(str,args)
        self.act="Act<"+','.join(self.subsys)+">"
        OnDemand.__init__(self, "Act", '_'.join(self.subsys))

    def generate_source(self,builddir):
        with open(os.path.join(self.cpypyqeddir,"ActTemplate.cc")) as f:
            act_template = f.read().format(act=self.act,modulename=self.modulename,classname=self.classname)
        with open(os.path.join(builddir,self.sourcefile),"w") as f:
            f.write(act_template)