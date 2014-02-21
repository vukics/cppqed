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
import ast
import scipy.interpolate
from scipy.integrate import quadrature
from scipy import exp

try:
  import matplotlib
  matplotlib.use('Agg')
  import matplotlib.pyplot as plt
  from matplotlib.backends.backend_pdf import PdfPages
  from matplotlib.font_manager import FontProperties

  plot=True
except ImportError:
  plot=False

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
def load_sv(fname, format=None):
  if format is None: return np.genfromtxt(fname)

  floatingReString=r'([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)'
  complexReString =r'\(\s*'+floatingReString+'\s*,\s*'+floatingReString+'\s*\)'

  return np.fromregex(fname,format.replace(r'+',r'\s*').replace('f',floatingReString).replace('c',complexReString),np.float)

## @}

def PTLA_postprocess(input):
  result=np.zeros((input.shape[0],6))
  result[:,[0,1]]=input[:,[0,1]]
  result[:,2]=(1+input[:,2])/2
  result[:,3]=(1-input[:,2])/2
  result[:,4]=input[:,3]/2
  result[:,5]=input[:,4]/2
  return result

def PLM_Evolved_postprocess(input):
  result=np.zeros((input.shape[0],5))
  result[:,[0,1]]=input[:,[0,1]]
  result[:,2]=input[:,2]**2+input[:,3]**2
  result[:,3]=input[:,2]
  result[:,4]=input[:,3]
  return result

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
  (works recursively). Values which end in `_local` are never imported.

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

  def _import_section(self,section=None):
    if section is None: section = self.test
    if self.cp.has_option(section,'import'):
      import_section=self.cp.get(section,'import')
      self._import_section(section=import_section) # import recursively
      for item in self.cp.items(import_section):
        if not self.cp.has_option(section,item[0]) and not item[0].endswith('_local'):
          self.cp.set(section, *item)
      self.cp.remove_option(section,'import')

  def get_option(self, name, default=None, required=False, section=None):
    """!
    Get configuration file keys in a safe way.
    \param name Name of the key.
    \param default Default value to return if key does not exist.
    \param required Fail if True and key does not exist.
    \param section The section name to look in, defaults to OptionsManager::test if None.
    \return The value to the key.

    This methods looks up the key `name` in the section name OptionsManager::test.
    """
    if section is None: section=self.test
    self._import_section(section=section)
    if self.cp.has_option(section,name):
      return self.cp.get(section,name)
    else:
      if not required: return default
      else: sys.exit("Error: required option \"{}\" not found in section {}.".format(name,section))

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
    return self.get_option('runmodes', section=section, default='generic').split(',')

  def _filter_runmodes(self, section):
    filter_runmodes=self.get_option('runmodes_'+self.test+'_local',section=section)
    if not filter_runmodes is None: filter_runmodes=filter_runmodes.split(',')
    for mode in self.runmodes(section=section):
      if not filter_runmodes is None and not mode in filter_runmodes: continue
      yield(mode)

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

class Plotter(OutputManager):
  """!
  \brief This is a helper class which helps with plotting functions to a pdf file.

  If the global variable `plot` is False, all functions are a no-op.
  """

  def _plot(self):
    return plot and not self.get_option('pdf') is None
  def start_pdf(self):
    """!
    \brief Initialize a new pdf file.

    The file is read from the configuration key `pdf`.
    """
    if not self._plot(): return
    self.pdf = PdfPages(os.path.join(self.outputdir,self.get_option('pdf')))
  def close_pdf(self):
    """!
    \brief Saves the pdf file to disc after all plots are finished.
    """
    if not self._plot(): return
    for n in plt.get_fignums():
      plt.figure(num=n)
      self._place_legend()
      self.finish_plot()
    self.pdf.close()
  def finish_plot(self):
    """!
    \brief Adds the current plot to the pdf file.
    """
    if not self._plot(): return
    self.pdf.savefig()
    plt.close()
  def figureLegendRight(self,ylabel,title,n):
    """!
    \brief Creates a new plot with figure legend right of the plot.
    \param ylabel The label of the y axis.
    \param title The title of the plot
    \param n The value number.
    """
    if not self._plot(): return
    if n in plt.get_fignums():
      plt.figure(num=n)
      return
    f = plt.figure(num=n,figsize=(11.6,8.2))
    f.add_axes([0.09, 0.1, 0.6, 0.75])
    plt.title(title)
    plt.ylabel(ylabel)
    plt.xlabel('t')
  def _place_legend(self):
    if not self._plot(): return
    fontP = FontProperties()
    fontP.set_size('small')
    leg=plt.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.,prop = fontP)
    llines=leg.get_lines()
    plt.setp(llines, linewidth=1.5)
  def plot(self,time,data,**kwargs):
    """!
    \brief Wraps matplotlibs plot function.
    \param time An array of time values.
    \param data An array of data values.
    \param **kwargs These are passed to `matplotlib.plot`.
    """
    if not self._plot(): return
    plt.plot(time,data,**kwargs)


def final_temperature(nTh):
  def fn(states):
    state=states[-1]
    n=np.arange(state.shape[0],dtype=float)
    expected_rho=np.diag(nTh**n/(1.+4)**(n+1))
    return np.sqrt(np.sum(np.abs(state-expected_rho)**2))
  return fn

class StateComparer(OutputManager):
  """!
  @ingroup Testclasses
  Tests final states of several trajectories by applying a given function.

  \ref StateComparer_keys "Configuration file keys" this class understands.
  """

  ## @addtogroup TestclassKeys
  #
  # @anchor StateComparer_keys
  # ## StateComparer configuration file keys
  # * `trajectories`: List of comma-separated trajectories which should be tested.
  # * `function`: A meta-function which should return the actual test function. The actual test function
  #     should accept the state array and return some epsilon value (the measure of the test).
  # * `parameters`: Tuple of function parameters passed to the meta function.
  #
  # The following configuration keys are read from the 'target'-sections.
  # * `runmodes_<test>`: For the compare test <test>, only use these runmodes.
  # * `epsilon_<runmode>_<test>`: Acceptable deviation for the given runmode and comparison test.

  def run(self):
    trajectories=self.get_option('trajectories',required=True).split(',')
    function=globals()[self.get_option('function',required=True)]
    parameters=ast.literal_eval(self.get_option('parameters'))
    if parameters is None: parameters=[]
    failure=False
    for traj in trajectories:
      for runmode in self._filter_runmodes(section=traj):
        statefile=self.output(runmode=runmode,section=traj,statefile=True)
        _,states,_=io.read(statefile)
        logging.debug("Evaluating {}.".format(os.path.basename(statefile)))
        eps=float(self.get_option('epsilon_'+runmode+'_'+self.test,section=traj,required=True))
        value=function(*parameters)(states)
        logging.debug("Value: {}, epsilon: {}".format(value,eps))
        if not value<eps:
          failure=True
          logging.debug("====== FAILED ======")
    if failure: sys.exit(-1)


class TrajectoryComparer(Plotter):
  """!
  @ingroup Testclasses
  Compares several trajectories to a reference trajectory by using function interpolation.

  \ref TrajectoryComparer_keys "Configuration file keys" this class understands.
  """

  ## @addtogroup TestclassKeys
  #
  # @anchor TrajectoryComparer_keys
  # ## TrajectoryComparer configuration file keys
  # * `pdf`: Save plots to this pdf file.
  # * `reference`: Section of reference trajectory
  # * `trajectories`: List of comma-separated trajectories which should be compared to the reference.
  #
  # The following configuration keys are read from the 'target'-sections.
  # * `runmodes_<test>`: For the compare test <test>, only use these runmodes.
  # * `columns_<test>`: Use these columns of the output files for the comparison.
  # * `epsilon_<runmode>_<test>`: List of acceptable deviations for the given runmode and comparison test.
  # * `postprocess_local`: Name of a global function which expects the data array as input and postprocesses the data.
  # * `format_local`: specifies which columns are floats (`f`) and which are complex numbers (`c`). Example:
  #     "f+f+c+c" will result in 6 columns, the two complex number columns are split into real and imaginary parts.
  # * `start_<test>`: The first row of the data lines to consider for the comparison test `<test>`.
  # * `length_<test>`: How many lines of data to consider for the comparison test `<test>`.


  def run(self):
    """!
    Runs the test.
    """
    trajectories=self.get_option('trajectories',required=True).split(',')
    failure=False
    self.start_pdf()
    reference_plotted=dict()
    for traj in trajectories:
      for runmode in self._filter_runmodes(section=traj):
        for n in range(len(self._get_columns(traj,runmode))):
          self.figureLegendRight(ylabel='value '+str(n+1), title=self.test, n=n)

          data,timeArray,data_label=self._get_data(section=traj,runmode=runmode,n=n)
          reference,reference_label=self._get_reference(section=traj,runmode=runmode,n=n)
          if not reference_plotted.has_key((reference_label,n)):
            self.plot(timeArray,reference(timeArray),label=reference_label)
            reference_plotted[(reference_label,n)]=True
          self.plot(timeArray,data(timeArray),label=data_label)
          logging.debug("Evaluating {}, value number {}.".format(data_label,n+1))
          eps=self._get_eps(runmode, traj, n)
          if not self._regression(reference,data,timeArray,eps):
            logging.debug("====== FAILED ======")
            failure=True
    self.close_pdf()
    if failure:
      sys.exit(-1)

  def _get_eps(self, runmode, section, n):
    return float(self.get_option('epsilon_'+runmode+'_'+self.test,section=section,required=True).split(',')[n])

  def _get_columns(self,section,runmode):
    return map(int,self.get_option('columns_'+self.test,section=section,required=True).split(','))

  def _get_reference(self,section,runmode,n):
    reference=self.get_option('reference',required=True)
    reference_runmode=self.runmodes(section=reference)[0]
    result=self._get_data(section=reference,runmode=reference_runmode,n=n)
    return result[0],result[2]

  def _get_data(self,section,runmode,n):
    fname=self.get_option('postprocess_local',section=section)
    format=self.get_option('format_local',section=section)
    length=self.get_option('length_'+self.test,section=section)
    start=self.get_option('start_'+self.test,section=section)
    postprocess=globals()[fname] if not fname is None else lambda x: x
    result=postprocess(load_sv(self.output(runmode=runmode,section=section),format=format))
    if not start is None:  result=result[int(start):]
    if not length is None: result=result[:int(length)]
    timeArray = result[:,0]
    data = result[:,self._get_columns(section,runmode)[n]]
    return self._interpolate(timeArray,data),timeArray,os.path.basename(self.output(runmode,section))

  def _interpolate(self,timeArray,array):
    return scipy.interpolate.interp1d(timeArray,array)

  def _regression(self, f1, f2, timeArray, eps) :
    t0=timeArray[ 0]
    t1=timeArray[-1]
    res=quadrature(lambda t : (f1(t)-f2(t))**2,t0,t1,maxiter=100)[0]
    logging.debug("Quadrature: {}, epsilon: {}".format(res,eps))
    return res<eps

def exponential(a,l):
  def fn(t):
    return a*exp(-l*t)
  return fn,"{}*exp(-{}*t)".format(a,l)

def FreeParticleX(x0,p0):
  def fn(t):
    return x0+2*p0*t
  return fn, "{}+2*{}*t".format(x0,p0)

def FreeParticleVarX(dx0,dp0):
  def fn(t):
    return (dx0+4.*dp0*t**2)**.5
  return fn, "({}+(4*{}*t)^2)^0.5"

class FunctionComparer(TrajectoryComparer):
  """!
  @ingroup Testclasses
  Compares several trajectories to a reference function by using function interpolation.

  \ref FunctionComparer_keys "Configuration file keys" this class understands.
  """

  ## @addtogroup TestclassKeys
  #
  # @anchor FunctionComparer_keys
  # ## FunctionComparer configuration file keys
  # * `reference_function`: Name of a global function, which should return a tuple of a unary function and a label used in plots.
  #
  # The following configuration keys are read from the 'target'-sections.
  # * `paramters_<test>`: List of tuples or single tuple which are passed to the reference function.
  #     Example: `[(1,5,3),(2,2,1)]` or `(1,5,3)`. If this is a list, each entry corresponds to a column of the data file,
  #     otherwise the same parameters are used for all columns.
  def _get_reference(self, section, runmode, n):
    reference = globals()[self.get_option('reference_function', required=True)]
    parameters=self.get_option('parameters_'+self.test, section=section)
    parameters=() if parameters is None else ast.literal_eval(parameters)
    if type(parameters)==list:parameters=parameters[n]
    return reference(*parameters)

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