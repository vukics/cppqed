import logging
from optparse import OptionParser
import ConfigParser
import sys

class Runner:
  def __init__(self, op, cp):
    self.op,self.cp = op,cp

  def run(self):
    pass

def main():
  op = OptionParser()
  cp = ConfigParser.SafeConfigParser()
  cp.read('testdriver.conf')
  sys.path.append(cp.get('Setup','modulepath'))
  import cpypyqed.io as cppio

  op.add_option("--testname")
  op.add_option("--testclass")
  op.add_option("--script")

  (options,args) = op.parse_args()

  if options.testclass:
    constructor = globals()[options.testclass]
    myTestclass = constructor(op,cp)
    myTestclass.run()

if __name__ == '__main__':
  main()