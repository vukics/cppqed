import logging
from optparse import OptionParser
import ConfigParser
import sys

op = OptionParser()
cp = ConfigParser.SafeConfigParser()


def main():
  cp.read('testdriver.conf')
  sys.path.append(cp.get('Setup','modulepath'))
  import cpypyqed.io as cppio



if __name__ == '__main__':
  main()