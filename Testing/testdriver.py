import logging
from optparse import OptionParser
import ConfigParser
import sys
import os
import errno
import subprocess
import numpy as np

logging.basicConfig(level=logging.DEBUG, format="%(asctime)s %(levelname)s %(message)s")

# helpers
def mkdir_p(path):
  """http://stackoverflow.com/questions/600268/mkdir-p-functionality-in-python
  """
  try:
    os.makedirs(path)
  except OSError, exc:
    if exc.errno == errno.EEXIST:
      pass
    else: raise

def rm_f(filename):
  """http://stackoverflow.com/a/10840586
  """
  try:
    os.remove(filename)
  except OSError as e: # this would be "except OSError, e:" before Python 2.6
    if e.errno != errno.ENOENT: # errno.ENOENT = no such file or directory
      raise # re-raise exception if a diff

class OptionsManager(object):
  def __init__(self, options, cp):
    self.options,self.cp = options,cp
    self.test = options.test
    if not self.test: sys.exit('--test missing')
    if self.cp.has_option(self.test,'import'):
      for item in self.cp.items(self.cp.get(self.test,'import')):
        if not self.cp.has_option(self.test,item[0]):
          self.cp.set(self.test, *item)


class OutputManager(OptionsManager):
  def __init__(self, *args, **kwargs):
    OptionsManager.__init__(self, *args, **kwargs)
    self.outputdir   = self.cp.get('Setup','outputdir')
    self.expecteddir = self.cp.get('Setup','expecteddir')
    mkdir_p(self.outputdir)
    self.script = self.options.script
    if not self.script: sys.exit('--script missing')
    self.runmodes = self.cp.get(self.test, 'runmodes').split(',')

  def output(self, runmode):
    return os.path.join(self.outputdir, self.test+'_'+runmode)

  def state_output(self, runmode):
    return self.output(runmode)+'.state'

  def clean(self, runmode):
    for file in (self.output(runmode),self.state_output(runmode)):
      rm_f(file)


# The test classes

class Runner(OutputManager):
  def run(self, clean=True):
    for runmode in self.runmodes:
      if clean: self.clean(runmode)
      command = self._build_commandline(runmode)
      logging.debug(subprocess.list2cmdline(command))
      ret = subprocess.call(command)
      if not ret==0: sys.exit(ret)

  def _extend_opts(self, options, section, option_prefix):
    for option in [ item[0] for item in self.cp.items(section) if item[0].startswith(option_prefix)]:
      options.extend(self.cp.get(section,option).split())

  def _build_commandline(self, runmode):
    result = [self.options.script]
    self._extend_opts(result, 'Scripts','opts')
    self._extend_opts(result, self.test,'opts')
    self._extend_opts(result, self.test,runmode)

    result.extend(('--evol',runmode))
    result.extend(('--o',self.output(runmode)))
    return result

class Verifier(OutputManager):
  def run(self):
    for runmode in self.runmodes:
      svfile = self.output(runmode)
      statefile = self.state_output(runmode)
      self._verify_ev(svfile,os.path.join(self.expecteddir,os.path.basename(svfile)))
      self._verify_state(statefile,os.path.join(self.expecteddir,os.path.basename(statefile)))
  def _load_sv(self,fname):
    return np.genfromtxt(fname)
  def _load_state(self,fname):
    return cpypyqed.io.read(fname)
  def _differ(self,res,exp):
    sys.exit("Error: {} and {} differ.".format(res,exp))
  def _equiv(self,res,exp):
    logging.debug("{} and {} are equivalent.".format(res,exp))
  def _verify_ev(self,res,exp):
    if not np.allclose(self._load_sv(res),self._load_state(exp)):
      self._differ(res,exp)
    else:
      self._equiv(res,exp)
  def _verify_state(self,res,exp):
    r_state,r_time = self._load_state(res)
    e_state,e_time = self._load_state(exp)
    if not (np.allclose(r_state,e_state) and np.allclose(r_time,e_time)):
      self._differ(res,exp)
    else:
      self._equiv(res,exp)

class Continuer(Runner):
  def run(self):
    self.cp.set(self.test,'opts_thisrun',self.cp.get(self.test,'firstrun'))
    Runner.run(self)
    self.cp.set(self.test,'opts_thisrun',self.cp.get(self.test,'secondrun'))
    Runner.run(self,clean=False)

def main():
  op = OptionParser()
  cp = ConfigParser.SafeConfigParser()

  op.add_option("--test")
  op.add_option("--testclass")
  op.add_option("--script")

  (options,args) = op.parse_args()

  if len(args) != 1:
    op.error("incorrect number of arguments")

  cp.read(args[0])
  sys.path.append(cp.get('Setup','modulepath'))
  # we can only load pycppqed after we know where to look for the cpypyqed module
  global cpypyqed
  import cpypyqed.io
  logging.info("Taking cpypyqed from {}".format(cpypyqed.io.__file__))

  if options.testclass:
    constructor = globals()[options.testclass]
    myTestclass = constructor(options,cp)
    myTestclass.run()

if __name__ == '__main__':
  main()