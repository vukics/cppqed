import logging
from optparse import OptionParser
import ConfigParser
import sys
import os
import errno
import subprocess
import numpy as np
import shutil

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

    # if the current section has the key 'import=othersection', import all keys from 'othersection'
    # if they are not present already
    import_section=self.get_option('import')
    if import_section:
      for item in self.cp.items(import_section):
        if not self.cp.has_option(self.test,item[0]):
          self.cp.set(self.test, *item)

  def get_option(self, name, default=None, required=False):
    if self.cp.has_option(self.test,name):
      return self.cp.get(self.test,name)
    else:
      if not required: return default
      else: sys.exit("Error: required option {} not found in section {}.".format(name,self.test))

class OutputManager(OptionsManager):
  def __init__(self, *args, **kwargs):
    OptionsManager.__init__(self, *args, **kwargs)
    self.outputdir   = self.cp.get('Setup','outputdir')
    self.expecteddir = self.cp.get('Setup','expecteddir')
    mkdir_p(self.outputdir)

  def runmodes(self,section=None):
    if section is None: section=self.test
    return self.cp.get(section, 'runmodes').split(',')

  def output(self, runmode, section=None, statefile=False):
    if section is None: section=self.test
    output = os.path.join(self.outputdir, section+'_'+runmode)
    if statefile: output+=".state"
    return output

  def clean(self, runmode):
    rm_f(self.output(runmode))
    rm_f(self.output(runmode,statefile=True))

# The test classes

class Runner(OutputManager):
  def run(self, clean=True, extra_opts=None, *args, **kwargs):
    for runmode in self.runmodes():
      if clean: self.clean(runmode)
      command = self._build_commandline(runmode,extra_opts)
      logging.debug(subprocess.list2cmdline(command))
      ret = subprocess.call(command, *args, **kwargs)
      if not ret==0: sys.exit(ret)

  def _extend_opts(self, options, section, option_prefix):
    for option in sorted([ item[0] for item in self.cp.items(section) if item[0].startswith(option_prefix)]):
      options.extend(self.cp.get(section,option).split())

  def _build_commandline(self, runmode, extra_opts=None):
    result = [self.options.script]
    if extra_opts: result+=extra_opts
    self._extend_opts(result, 'Scripts','opts')
    self._extend_opts(result, self.test,'opts')
    self._extend_opts(result, self.test,runmode)

    result.extend(('--evol',runmode))
    result.extend(('--o',self.output(runmode)))
    return result

class PythonRunner(Runner):
  def run(self, clean=True):
    cpypyqed_builddir = self.options.cpypyqed_builddir
    cpypyqed_config   = self.options.cpypyqed_config
    env = os.environ.copy()
    if cpypyqed_builddir:
      env['CPYPYQED_BUILDDIR']=cpypyqed_builddir
      if clean: shutil.rmtree(os.path.join(cpypyqed_builddir,'cppqedmodules'),ignore_errors=True)
    if cpypyqed_config: env['CPYPYQED_CONFIG']=cpypyqed_config
    env['PYTHONPATH']=self.cp.get('Setup','modulepath')
    if extra_opts is None: extra_opts = []
    if self.options.configuration.lower()=="debug": extra_opts += ['--debug']
    Runner.run(self,clean=clean,extra_opts=extra_opts,env=env,*args,**kwargs)

class Verifier(OutputManager):
  def __init__(self,*args,**kwargs):
    OutputManager.__init__(self,*args,**kwargs)
    self.thisSection  = self.get_option('this',default=self.test)
    self.otherSection = self.get_option('other')

  def run(self):
    mode=self.get_option('verify')
    if mode is None or mode=='full':
      self._verify_full()
    elif mode=='outcome':
      self._verify_outcome()
  def _verify_full(self):
    for runmode in self.runmodes(section=self.thisSection):
      self._verify_ev(self._thisOutput(runmode),self._otherOutput(runmode))
      self._verify_state(self._thisOutput(runmode,statefile=True),self._otherOutput(runmode,statefile=True))
  def _thisOutput(self,runmode,statefile=False):
    return self.output(runmode,section=self.thisSection,statefile=statefile)
  def _otherOutput(self,runmode,statefile=False):
    if self.otherSection is None:
      return os.path.join(self.expecteddir,os.path.basename(self._thisOutput(runmode,statefile)))
    else:
      return self.output(runmode,section=self.otherSection,statefile=statefile)
  def _load_sv(self,fname):
    return np.genfromtxt(fname)
  def _load_state(self,fname):
    return io.read(fname)
  def _differ(self,this,other):
    sys.exit("Error: {} and {} differ.".format(this,other))
  def _equiv(self,this,other):
    logging.debug("{} and {} are equivalent.".format(this,other))
  def _verify_ev(self,this,other):
    if not np.allclose(self._load_sv(this),self._load_sv(other)): self._differ(this,other)
    else: self._equiv(this,other)
  def _verify_state(self,this,other):
    _,r_state,r_time = self._load_state(this)
    _,e_state,e_time = self._load_state(other)
    if not (np.allclose(r_state,e_state) and np.allclose(r_time,e_time)): self._differ(this,other)
    else: self._equiv(this,other)
  def _verify_outcome(self,this,other):
    _,r_state,r_time=self._load_state(this)
    _,e_state,e_time=self._load_state(other)
    if not (np.allclose(r_state[-1],e_state[-1]) and np.allclose(r_time[-1],e_time[-1])):
      self._differ(this,other)
    else: self._equiv(this,other)

class VerifiedRunner(Runner,Verifier):
  def run(self):
    Runner.run(self)
    Verifier.run(self)

def ContinueFactory(base,*args,**kwargs):
  class GenericContinuer(base):
    def run(self):
      self.cp.set(self.test,'opts_thisrun',self.get_option('firstrun',default=''))
      base.run(self)
      self.cp.set(self.test,'opts_thisrun',self.get_option('secondrun',default=''))
      base.run(self,clean=False)
  return GenericContinuer(*args,**kwargs)
class Continuer(Runner, GenericContinuer):
  """!
  @ingroup Testclasses
  GenericContinuer version of Runner.

  \ref GEnericContinuer_keys "Configuration file keys" this class understands.
  """

  ## @addtogroup TestclassKeys
  #
  # @anchor Continuer
  # ## Continuer configuration file keys
  # See \ref GenericContinuer_keys "GenericContinuer keys".

  def run(self, *args, **kwargs):
    """!
    Delegates to GenericContinuer::continued_run().
    """
    GenericContinuer.continued_run(self, Runner.run, *args, **kwargs)

class PythonContinuer(PythonRunner, GenericContinuer):
  """!
  @ingroup Testclasses
  GenericContinuer version of PythonRunner.

  \ref GEnericContinuer_keys "Configuration file keys" this class understands.
  """

  ## @addtogroup TestclassKeys
  #
  # @anchor PythonContinuer
  # ## PythonContinuer configuration file keys
  # See \ref GenericContinuer_keys "GenericContinuer keys".

  def run(self, *args, **kwargs):
    """!
    Delegates to GenericContiuer::continued_run().
    """
    GenericContinuer.continued_run(self, PythonRunner.run, *args, **kwargs)

class CompileFail(OptionsManager):
  def run(self):
    error=self.options.error
    cmake=self.cp.get('Setup','cmake')
    builddir=self.cp.get('Setup','builddir')
    command=[cmake,'--build',builddir,'--target',self.options.script]
    logging.debug(subprocess.list2cmdline(command))
    p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (std,err) = p.communicate()
    returncode = p.returncode
    if returncode == 0:
      sys.exit("Compilation was successfull, but failure was expected.")
    if not error in std:
      sys.exit("Compilation failed as expected, but {} was not found in the error message.".format(error))

def main():
  op = OptionParser()
  cp = ConfigParser.SafeConfigParser()

  op.add_option("--test", help="the name of the test, and the name of the section in the config file")
  op.add_option("--testclass", help="the name of the testclass to use, must implement run()")
  op.add_option("--script", help="the script to run or the target to compile")
  op.add_option("--configuration", help="debug or release")
  op.add_option("--cpypyqed_builddir", help="directory for on-demand module compilation")
  op.add_option("--cpypyqed_config", help="configure file for on-demand module compilation")
  op.add_option("--error", metavar='STRING', help="string to expect in the compilation failure for CompileFail class")

  (options,args) = op.parse_args()

  if len(args)==0: op.error("Need configuration file(s) as argument(s).")
  cp.read(args)
  sys.path.insert(0,cp.get('Setup','modulepath'))
  # we can only load the io module after we know where to look for the cpypyqed package
  global io
  if options.configuration.lower()=="release":
    import cpypyqed.io as io
  elif options.configuration.lower()=="debug":
    import cpypyqed_d.io_d as io
  logging.info("Taking cpypyqed from {}".format(io.__file__))

  if options.testclass:
    constructor = globals()[options.testclass]
    myTestclass = constructor(options,cp)
    myTestclass.run()

if __name__ == '__main__':
  main()