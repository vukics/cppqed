# -*- coding: utf-8 -*-

import ConfigParser
import os
import sys
import errno
import tempfile
import subprocess
import shutil

from cpypyqed_config import cppqed_build_type,cppqed_module_suffix
if cppqed_build_type=="release":
  from ..core import core_git
else:
  from ..core_d import core_git

# The compilation directory can be customized by changing ondemand.cpypyqed_builddir
cpypyqed_builddir = "~/.cpypyqed"

def mkdir_p(path):
    """http://stackoverflow.com/questions/600268/mkdir-p-functionality-in-python
    """
    try:
        os.makedirs(path)
    except OSError, exc:
        if exc.errno == errno.EEXIST:
            pass
        else: raise

class OnDemand(object):

    def __init__(self, basename, classid, makerfunction=None):
        self.config = ConfigParser.SafeConfigParser(
                        dict(delete_temp='True',cmake_opts='',compiler=''))
        self.dir = os.path.expanduser(cpypyqed_builddir)
        self.configfile = os.path.join(self.dir,"config.txt")
        self.config.read(self.configfile)
        self.logfile = os.path.join(self.dir,"build.log")
        self.cpypyqeddir = os.path.dirname(__file__)
        self.cmaketemplate = os.path.join(self.cpypyqeddir,"CMakeListsTemplate.txt")
        self.packagename = "cppqedmodules"
        self.modulepath = os.path.join(self.dir,self.packagename)
        self.basename = basename
        self.classid = classid
        self.classname = self.basename+self.classid
        self.modulename = self.classname+cppqed_module_suffix
        self.sourcefile = self.modulename + ".cc"
        self.library = self.modulename+".so"
        self.fullclass = self.packagename+"."+self.modulename+"."+self.classname
        if makerfunction:
          self.makerfunction = self.packagename+"."+self.modulename+"."+makerfunction
        else:
          self.makerfunction = None
        if not sys.path.count(self.dir):
            sys.path.insert(0,self.dir)
        mkdir_p(self.modulepath)
        f = open(os.path.join(self.modulepath,"__init__.py"),"w")
        f.close()

    def import_class(self,s):
        r"""Import a class specified by the string s.

        :param s: The name of the class to import, e.g. 'mypackage.mymodule.myclass'
        :returns: The class.
        """
        components = s.split('.')
        modulename = '.'.join(components[:-1])
        classname = components[-1]
        module = __import__(modulename, globals(), locals(), [classname])
        return getattr(module,classname)

    def maker(self,*args,**kwargs):
        try:
          thisClass = self.import_class(self.fullclass)
        except ImportError:
          self.build()
          thisClass = self.import_class(self.fullclass)
        if not thisClass.core_git == core_git:
          os.remove(os.path.join(self.modulepath,self.library))
          sys.exit("Error: {} was compiled with different library version, removing.\nPlease restart the script.\n".format(self.fullclass))
        if self.makerfunction:
            thisClass = self.import_class(self.makerfunction)
        return thisClass(*args,**kwargs) 

    def _check_return_value(self,ret, errormsg=None):
        if not ret==0:
          if errormsg: sys.stderr.write(errormsg)
          raise(RuntimeError("Could not compile the on-demand python module. Refer to "+self.logfile+" for error messages."))

    def build(self):
        builddir = tempfile.mkdtemp(dir=self.dir)
        try:
          self.generate_source(builddir)
          with open(self.cmaketemplate) as f:
              cmake = f.read()
          cmake = cmake.format(modulename=self.modulename, sourcefile=self.sourcefile)
          with open(os.path.join(builddir,"CMakeLists.txt"),"w") as f:
              f.write(cmake)
          with open(self.logfile, "w") as log:
              print("Configuring the build project for {}, please stand by...".format(self.classname))
              opts=self.config.get('Setup','cmake_opts').split()
              compiler=self.config.get('Setup','compiler')
              if compiler: opts.append('-DCMAKE_CXX_COMPILER={}'.format(compiler))
              if cppqed_build_type=="debug":
                opts.append("-DCMAKE_BUILD_TYPE=Debug")
              elif cppqed_build_type=="release":
                opts.append("-DCMAKE_BUILD_TYPE=Release")
              ret = subprocess.call(args=['cmake','.']+opts,cwd=builddir,stdout=log,stderr=subprocess.STDOUT)
              self._check_return_value(ret, errormsg="Error: cmake failed")
              print("Building {}, please stand by...".format(self.classname))
              ret = subprocess.call(args=('make',),cwd=builddir,stdout=log,stderr=subprocess.STDOUT)
              self._check_return_value(ret, errormsg="Error: make failed")
          shutil.copy(os.path.join(builddir,self.library),self.modulepath)
        finally:
          if self.config.getboolean('Setup','delete_temp'):
            shutil.rmtree(builddir, ignore_errors=True)

    def generate_source(self, builddir):
        pass

