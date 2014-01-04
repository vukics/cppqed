# -*- coding: utf-8 -*-

r"""This module provides the infrastructure to instantiate and build C++QED wrapper classes on demand.
The main class to implementing this functionality is :class:`OnDemand`.
"""

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
    """Equivalent of `mkdir -p`.

    http://stackoverflow.com/questions/600268/mkdir-p-functionality-in-python

    :param str path: The directory to create.
    """
    try:
        os.makedirs(path)
    except OSError, exc:
        if exc.errno == errno.EEXIST:
            pass
        else: raise

class OnDemand(object):
    r"""This class serves as a base to classes which need to be compiled on the fly in cpypyqed.

    This is needed for class templates which cannot be pre-instantiated at cpypyqed compile time
    because there are too many possibilities of template parameters. A typical example is the class
    :core:`Composite`.

    Classes deriving from :class:`OnDemand` need to implement :meth:`generate_source`.

    :param str basename: The descriptive name of the underlying class
    :param str classid: An id which encodes the template parameters needed to instantiate the C++ class
    :param makerfunction: Use this makerfunction instead of the class constructor to create an instances of the
      underlying class (default None)
    :type makerfunction: function or None
    """

    cpypyqeddir=None
    r"""This will be set to the directory where the C++ source templates are found,
    typically used in the implementation of :meth:`generate_source`."""

    def __init__(self, basename, classid, makerfunction=None):
        self.config = ConfigParser.SafeConfigParser(
                        dict(delete_temp='True',cmake_opts='',compiler='',
                             cppqed_dir='',cppqed_dir_release='',cppqed_dir_debug=''))
        if not self.config.has_section("Setup"):
          self.config.add_section("Setup")
        self.dir = os.path.expanduser(cpypyqed_builddir)
        if os.environ.has_key('CPYPYQED_BUILDDIR'):
          self.dir = os.path.expanduser(os.environ['CPYPYQED_BUILDDIR'])
        self.configfiles = [os.path.join(self.dir,"config.txt")]
        if os.environ.has_key('CPYPYQED_CONFIG'):
          self.configfiles.append(os.path.expanduser(os.environ['CPYPYQED_CONFIG']))
        self.config.read(self.configfiles)
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

    def _import_class(self,s):
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
        r"""Returns an instance of the underlying class represented by this :class:`OnDemand` object.
        The arguments of this function are passed through to the constructor of the underlying class.
        If the on-demand module holding the class has not been built yet or the C++QED library version
        does not fit, the module is compiled on the fly.
        """
        try:
          thisClass = self._import_class(self.fullclass)
        except ImportError:
          self._build()
          thisClass = self._import_class(self.fullclass)
        if not thisClass.core_git == core_git:
          os.remove(os.path.join(self.modulepath,self.library))
          sys.exit("Error: {} was compiled with different library version, removing.\nPlease restart the script.\n".format(self.fullclass))
        if self.makerfunction:
            thisClass = self._import_class(self.makerfunction)
        return thisClass(*args,**kwargs) 

    def _check_return_value(self,ret, errormsg=None):
        if not ret==0:
          if errormsg: sys.stderr.write(errormsg)
          raise(RuntimeError("Could not compile the on-demand python module. Refer to "+self.logfile+" for error messages."))

    def _build(self):
        r"""Build the module holding the underlying class.
        """
        builddir = os.path.join(self.dir,self.modulename)
        shutil.rmtree(builddir,ignore_errors=True)
        mkdir_p(builddir)
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
              cppqed_dir=os.path.expanduser(self.config.get('Setup','cppqed_dir'))
              if compiler: opts.append('-DCMAKE_CXX_COMPILER={}'.format(compiler))
              if cppqed_build_type=="debug":
                opts.append("-DCMAKE_BUILD_TYPE=Debug")
                cppqed_dir_debug=os.path.expanduser(self.config.get('Setup','cppqed_dir_debug'))
                if cppqed_dir_debug: cppqed_dir=cppqed_dir_debug
              elif cppqed_build_type=="release":
                opts.append("-DCMAKE_BUILD_TYPE=Release")
                cppqed_dir_release=os.path.expanduser(self.config.get('Setup','cppqed_dir_release'))
                if cppqed_dir_release: cppqed_dir=cppqed_dir_release
              if cppqed_dir: opts.append('-DCPPQED_DIR={}'.format(cppqed_dir))
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
        r"""This creates the C++ source file of the python module. It is a stub and has to be implmented by
        deriving classes.

        Implementations have to save the result to `builddir`. They can access template files in the directory
        :attr:`cpypyqeddir`.

        :param str builddir: The directory where the module is going to be built. The result has to be saved here.
        """
        pass

