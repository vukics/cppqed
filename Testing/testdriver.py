## @package testdriver
# This is the Python testdriver for the \ref testsuite.
#
# It is intended to be used with the CMake CTest utility.
# When called with the parameter `--testclass=<TESTCLASS>`, it calls the `run`
# method of the specified runner class. Success of a test is indicated by the
# return value 0.

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

## @name Helper functions
# @{

## Create a directory with parent directories.
# @param path The path to create.
#
# From this [stackoverflow question](http://stackoverflow.com/questions/600268/mkdir-p-functionality-in-python)
def mkdir_p(path):
  try:
    os.makedirs(path)
  except OSError, exc:
    if exc.errno == errno.EEXIST:
      pass
    else: raise

## Remove a file without error if it doesn't exist.
# @param filename The file to delete.
#
# From this [stackoverflow question](http://stackoverflow.com/a/10840586)
def rm_f(filename):
  try:
    os.remove(filename)
  except OSError as e: # this would be "except OSError, e:" before Python 2.6
    if e.errno != errno.ENOENT: # errno.ENOENT = no such file or directory
      raise # re-raise exception if a diff

## Loads a trajectory file.
# \param fname File name to load from.
# \return array Numpy array.
def load_sv(fname):
  return np.genfromtxt(fname)

## @}

## @defgroup TestclassHelpers Helpers
# @ingroup Testclasses
# \brief Helper base classes to test classes.
# These classes cannot be used as a test class directly, but serve as base to other test classes
# and define some \ref TestclassKeys "configuration file keys" and \ref TestclassOptions "command line options".

class OptionsManager(object):
  """!
  @ingroup TestclassHelpers
  \brief Stores command line options and configuration file keys.

  Each OptionsManager instance has its own section in the configuration file, named after
  the current test name (OptionsManager::test). If the current section has the key
  `import=othersection`, import all keys from `othersection` if they are not present already
  (doesn't work recursively).

  \ref OptionsManager_options "Command line" options this class understands.
  """

  ## @addtogroup TestclassOptions
  #
  # @anchor OptionsManager_options
  # ## OptionsManager command line options
  # * `--test=<testname>`: The name of the test. This defines the section in the configuration file
  #     and also ends up in output files etc.

  def __init__(self, options, cp):
    """!
    @param options optparse.Values: object holding all the command line options.
    @param cp ConfigParser: ConfigParser instance holding all configuration file keys.
    """

    ## optparse.Values: command line options
    self.options = options
    ## ConfigParser: configuration file keys
    self.cp = cp
    ## The name of the current test
    self.test = options.test
    if not self.test: sys.exit('--test missing')

    import_section=self.get_option('import')
    if import_section:
      for item in self.cp.items(import_section):
        if not self.cp.has_option(self.test,item[0]):
          self.cp.set(self.test, *item)

  def get_option(self, name, default=None, required=False):
    """!
    Get configuration file keys in a safe way.
    \param name Name of the key.
    \param default Default value to return if key does not exist.
    \param required Fail if True and key does not exist.
    \return The value to the key.

    This methods looks up the key `name` in the section name OptionsManager::test.
    """
    if self.cp.has_option(self.test,name):
      return self.cp.get(self.test,name)
    else:
      if not required: return default
      else: sys.exit("Error: required option {} not found in section {}.".format(name,self.test))

class OutputManager(OptionsManager):
  """!
  @ingroup TestclassHelpers
  \brief Manages output files for different run modes.

  \ref OutputManager_keys "Configuration file keys" this class understands.
  """

  ## @addtogroup SetupKeys
  #
  # * `outuptdir`: All output files end up here.
  # * `expecteddir`: Where to look for pre-run simulations to compare test runs to.


  ## @addtogroup TestclassKeys
  #
  # @anchor OutputManager_keys
  # ## OutputManager configuration file keys
  # * `runmodes`: comma separated list of runmodes (single, master ensemble)


  def __init__(self, *args, **kwargs):
    """!
    Arguments are passed through to OptionsManager.
    """
    OptionsManager.__init__(self, *args, **kwargs)
    ## All output files end up here.
    self.outputdir   = self.cp.get('Setup','outputdir')
    ## Where to look for pre-run simulations to compare test runs to.
    self.expecteddir = self.cp.get('Setup','expecteddir')
    mkdir_p(self.outputdir)

  def runmodes(self,section=None):
    """!
    Return runmodes.
    \param section (optional) String: Where to look up the runmodes, take current test section if not specified.
    \return A list of runmodes in this section.
    """
    if section is None: section=self.test
    return self.cp.get(section, 'runmodes').split(',')

  def output(self, runmode, section=None, statefile=False):
    """!
    The name of the output file for a given runmode.
    \param runmode String: The runmode for which the filename should be generated.
    \param section (optional) String: Output file name for which section, current test section if left empty.
    \param statefile (optional) Boolean: By default generate the file name for a trajectory file. If set to true
      generate the file name for a state file.
    \return Full path including OutputManager::outputdir.
    """
    if section is None: section=self.test
    if runmode == "generic":
      output = os.path.join(self.outputdir, section)
    else:
      output = os.path.join(self.outputdir, section+'_'+runmode)
    if statefile: output+=".state"
    return output

  def clean(self, runmode):
    """!
    Delete the trajectory file and state file for a given runmode.
    \param runmode String: The runmode for which output files should be deleted.
    """
    rm_f(self.output(runmode))
    rm_f(self.output(runmode,statefile=True))


# The test classes

class Runner(OutputManager):
  """!
  @ingroup Testclasses
  Runs a script repeatedly for all declared runmodes and succeeds if the scripts do.

  \ref Runner_keys "Configuration file keys" this class understands.
  """
  def run(self, clean=True, extra_opts=None, *args, **kwargs):
    """!
    The method to run the test.
    \param clean (optional) `Boolean`: Whether to remove old output before running the test.
    \param extra_opts (optional) `List`: Additional command line options appended to the script call.
    \param args passed through to `subprocess.call`
    \param kwargs passed through to `subprocess.call`

    This method terminates the test driver with a return value equal to that of the script call
    if one of the scripts fail.
    """
    for runmode in self.runmodes():
      if clean: self.clean(runmode)
      command = self._build_commandline(runmode,extra_opts)
      logging.debug(subprocess.list2cmdline(command))
      ret = subprocess.call(command, *args, **kwargs)
      if not ret==0: sys.exit(ret)

  ## @addtogroup TestclassKeys
  #
  # @anchor Runner_keys
  # ## Runner configuration file keys
  # * `opts*`: The command line options used for running the script, multiple keys matching `opts*` can be given
  # * `single*`, `master*`, `ensemble*`: Additional options for the specific runmodes. Multiple keys
  #     matching `<runmode>*` can be given.
  #
  # Example usage:
  #
  #     # The options used for running the scripts, multiple keys can be given if they match opts*
  #     opts=--etat 8 --sdf 3
  #     opts1=--dc 0 --Dt 0.1 --NDt 10
  #
  #     # runmode specific options
  #     single=...
  #     single1=...
  #     ensemble=...
  #     master=...

  def _extend_opts(self, options, section, option_prefix):
    for option in sorted([ item[0] for item in self.cp.items(section) if item[0].startswith(option_prefix)]):
      options.extend(self.cp.get(section,option).split())

  def _build_commandline(self, runmode, extra_opts=None):
    result = [self.options.script]
    if extra_opts: result+=extra_opts

    ## @addtogroup SetupKeys
    #
    # * `opts`: Script command line options added to all scripts

    self._extend_opts(result, 'Setup','opts')
    self._extend_opts(result, self.test,'opts')
    self._extend_opts(result, self.test,runmode)

    if not runmode=="generic": result.extend(('--evol',runmode))
    result.extend(('--o',self.output(runmode)))
    return result

class PythonRunner(Runner):
  """!
  @ingroup Testclasses
  Runs a cpypyqed script repeatedly for all declared runmodes and succeeds if the scripts do.

  \ref PythonRunner_options "Configuration file keys" this class understands.
  """

  ## @addtogroup TestclassOptions
  #
  # @anchor PythonRunner_options
  # ## PythonRunner command line options
  # * `--cpypyqed_builddir=<dir>`: Directory for on-demand compilation
  # * `--cpypyqed_config=<config-file>`: Configuration file for on-demand compilation

  def run(self, clean=True, extra_opts=None, *args, **kwargs):
    """!
    The method to run the test.
    \param clean (optional) `Boolean`: Whether to remove old output before running the test.
    \param extra_opts (optional) `List`: Additional command line options appended to the script call.
    \param args passed through to Runner.run()
    \param kwargs passed through to Runner.run()

    This method terminates the test driver with a return value equal to that of the script call
    if one of the scripts fail.
    """
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
  """!
  @ingroup Testclasses
  Verifies the output of a script 'this' to an expected output or the output of some other test run 'other'

  \ref Verifier_keys "Configuration file keys" this class understands.
  """

  ## @addtogroup TestclassKeys
  #
  # @anchor Verifier_keys
  # ## Verifier configuration file keys
  # The Verifier compares some test 'this' to another test 'other'.
  # * `this`: Test name of 'this', by default the current test if missing
  # * `other`: Testname of 'other', by default the results from the directory of expected results
  #     (OutputManager::expecteddir)
  # * `verify`: Verify that both trajectories are exactly equal (default if this key is missing or
  #     `verify=full`), or verify that the last outcome of the simulation is equal, e.g. timesteps may differ
  #     (`verify=outcome`)
  #
  # If `this=some_test` is specified, it is probably also a good idea to `import=some_test` to keep
  # the runmodes in sync. Currently the directory of expected results is `Testing/expected`, it is kept
  # under version control so that changes in the output of the scripts are noticed.


  def __init__(self,*args,**kwargs):
    """!
    \param args passed through to OutputManager
    \param kwargs passed through to OutputManager
    """
    OutputManager.__init__(self,*args,**kwargs)
    self.thisSection  = self.get_option('this',default=self.test)
    self.otherSection = self.get_option('other')

  def run(self):
    """!
    Run the test.
    """
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
  def _differ(self,this,other):
    sys.exit("Error: {} and {} differ.".format(this,other))
  def _equiv(self,this,other):
    logging.debug("{} and {} are equivalent.".format(this,other))
  def _verify_ev(self,this,other):
    if not np.allclose(load_sv(this),load_sv(other)): self._differ(this,other)
    else: self._equiv(this,other)
  def _verify_state(self,this,other):
    _,r_state,r_time = io.read(this)
    _,e_state,e_time = io.read(other)
    if not (np.allclose(r_state,e_state) and np.allclose(r_time,e_time)): self._differ(this,other)
    else: self._equiv(this,other)
  def _verify_outcome(self,this,other):
    _,r_state,r_time=io.read(this)
    _,e_state,e_time=io.read(other)
    if not (np.allclose(r_state[-1],e_state[-1]) and np.allclose(r_time[-1],e_time[-1])):
      self._differ(this,other)
    else: self._equiv(this,other)

class VerifiedRunner(Runner,Verifier):
  """!
  @ingroup Testclasses
  Combines the functionality of Runner and Verifier to a single test.
  """

  def run(self):
    """!
    Run the test.
    """
    Runner.run(self)
    Verifier.run(self)

class GenericContinuer(OptionsManager):
  """!
  @ingroup TestclassHelpers
  This class hosts continued_run(), which will run and then continue a script.
  """

  ## @addtogroup TestclassKeys
  #
  # @anchor GenericContinuer_keys
  # ## GenericContinuer configuration file keys
  # * `firstrun`: script options for the first run
  # * `secondrun`: script options for the second run

  def continued_run(self, runfn, *args, **kwargs):
    """!
    Run, then continue a script.
    \param runfn Function: The run function to call.
    \param args passed through to `runfn`
    \param kwargs passed through to `runfn`
    """
    runfn(self, extra_opts=self.get_option('firstrun',default='').split(), *args, **kwargs)
    runfn(self, clean=False, extra_opts=self.get_option('secondrun',default='').split(), *args, **kwargs)

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

class CompileTarget(OptionsManager):
  """!
  @ingroup Testclasses
  \brief This test tries to compile a %CMake target.

  If the `--error` option is not given,
  the test succeeds if the target can be compiled, otherwise the test succeeds if the
  compilation fails and the string specified together with `--error` is found.

  \ref CompileTarget_options "Command line options" this class understands.
  """

  ## @addtogroup SetupKeys
  #
  # * `cmake`: Path of the cmake executable
  # * `builddir`: Top-level build directory
  # * `

  ## @addtogroup TestclassOptions
  #
  # @anchor CompileTarget_options
  # ## CompileTarget command line options
  # * `--script`: The name of the target to compile.

  ## @addtogroup TestclassKeys
  #
  # @anchor CompileTarget_keys
  # ## CompileTarget configuration file keys
  # * `error`: Turn on "failure mode". The error message which is expected in the output.
  # * `dependencies`: Space separated list of dependencies to compile first. These are
  #     always required to succeed, independent of the presence of `error`.

  def run(self):
    """!
    Runs the test.
    """
    error=self.get_option('error')
    cmake=self.cp.get('Setup','cmake')
    builddir=self.cp.get('Setup','builddir')
    command=[cmake,'--build',builddir,'--target']
    dependencies=self.get_option('dependencies',default="").split()
    for dep in dependencies:
      logging.debug(subprocess.list2cmdline(command+[dep]))
      p = subprocess.Popen(command+[dep], stdout=subprocess.PIPE,stderr=subprocess.PIPE)
      (std,err) = p.communicate()
      if not p.returncode==0:
        sys.exit("Compilation of dependency {} for {} failed.".format(dep,self.options.script))
    logging.debug(subprocess.list2cmdline(command+[self.options.script]))
    p = subprocess.Popen(command+[self.options.script], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (std,err) = p.communicate()
    returncode = p.returncode
    if error is None:
      if returncode != 0:
        sys.exit("Compilation of {} failed.".format(self.options.script))
    else:
      if returncode == 0:
        sys.exit("Compilation was successful, but failure was expected.")
      if not error in std:
        sys.exit("Compilation failed as expected, but {} was not found in the error message.".format(error))

def main():
  """!
  \brief Main function of the Python test driver.

  Command line options are defined here. It is responsible of loading the right `cpypyqed` module
  (release or debug) as well as instantiating and running the test class.
  """
  op = OptionParser()
  cp = ConfigParser.SafeConfigParser()

  op.add_option("--test", help="the name of the test, and the name of the section in the config file")
  op.add_option("--testclass", help="the name of the testclass to use, must implement run()")
  op.add_option("--script", help="the script to run or the target to compile")
  op.add_option("--configuration", help="debug or release")
  op.add_option("--cpypyqed_builddir", help="directory for on-demand module compilation")
  op.add_option("--cpypyqed_config", help="configure file for on-demand module compilation")

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