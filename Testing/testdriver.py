import logging
from optparse import OptionParser
import ConfigParser
import sys
import os
import errno
import subprocess

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
    self.section = self.test
    if self.cp.has_option(self.section,'import'):
      for item in self.cp.items(self.cp.get(self.section,'import')):
        if not self.cp.has_option(self.section,item[0]):
          self.cp.set(self.section, *item)


class OutputManager(OptionsManager):
  def __init__(self, *args, **kwargs):
    OptionsManager.__init__(self, *args, **kwargs)
    self.outputdir = self.cp.get('Setup','outputdir')
    mkdir_p(self.outputdir)
    self.script = self.options.script
    if not self.script: sys.exit('--script missing')
    self.runmodes = self.cp.get(self.section, 'runmodes').split(',')

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
    self._extend_opts(result, self.section,'opts')
    self._extend_opts(result, self.section,runmode)

    result.extend(('--evol',runmode))
    result.extend(('--o',self.output(runmode)))
    return result

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
  import cpypyqed.io as cppio

  if options.testclass:
    constructor = globals()[options.testclass]
    myTestclass = constructor(options,cp)
    myTestclass.run()

if __name__ == '__main__':
  main()