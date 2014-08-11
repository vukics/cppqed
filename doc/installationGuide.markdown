Installation guide {#installationguide}
==================

# Quick start: Installation from packages (Ubuntu, Arch Linux) {#installationguidequickstart}

Binary packages have been prepared for Arch Linux and [Ubuntu](https://launchpad.net/~raimar-sandner/+archive/cppqed-development),
which might also work on Debian (not tested).  Installation from binary packages is recommended for users
who want to try C++QED or use it exclusively on its \ref userguide "highest level", i.e. writing
scripts which use \ref genericelements "elements and interactions" already implemented in C++QED. If
you are not using Ubuntu or want to develop your own elements and interactions
you have to install C++QED from \ref installationguidefromsource "source".

### Ubuntu

All packages are compiled for Ubuntu Precise (12.04 LTS) and Saucy (13.10). You can choose between the latest release and daily builds. For the release:

    $ sudo add-apt-repository ppa:raimar-sandner/cppqed-development
    $ sudo apt-get update
    $ sudo apt-get install libcppqedelements-2.10-dev cppqedscripts-2.10 python-cpypyqed-2.10

For the daily builds:

    $ sudo add-apt-repository ppa:raimar-sandner/cppqed-daily
    $ sudo apt-get update
    $ sudo apt-get install libcppqedelements-daily-dev cppqedscripts python-cpypyqed

The current documentation can be installed (Ubuntu Saucy only) with

    $ sudo apt-get install cppqed-doc-2.10
    # or
    $ sudo apt-get install cppqed-doc-daily

Here is an overview of the available packages and what they install (the exact version numbers may change over time).

| Package | Description |
|---------|-------------|
| `cppqedscripts-2.10` | Installs the compiled example scripts to `/usr/bin` and their sourcecode for inspection to `/usr/share/doc/cppqedscripts-2.10/examples/`.|
| `cppqedscripts-debug-2.10` | Installs the example scripts compiled in debug mode to `/usr/bin`. This is a separate package because they need a lot of disc space.|
| `libcppqedcore-2.10-2` | Installs the core library to `/usr/lib/x86_64-linux-gnu`. |
| `libcppqedcore-2.10-dev` | Installs the core header files to `/usr/include/CPPQED-2.10` and cmake auxilliary files to `/usr/lib/x86_64-linux-gnu/cmake`. It also Installs example projects for custom elements and scripts to `/usr/share/doc/libcppqedcore-2.10-dev/examples` which can be copied and used as starting point for own projects|
| `libcppqedelements-2.10-2` | Installs the elements library to `/usr/lib/x86_64-linux-gnu`.|
| `libcppqedelements-2.10-dev` | Installs the elements header files to `/usr/include/CPPQED-2.10`.|
| `python-cpypyqed-2.10` | Installs the Python wrapper for C++QED to `/usr/lib/python2.7/dist-packages/cpypyqed`. Example scripts which demonstrate how to use the wrapper are installed in `/usr/share/doc/python-cpypyqed-2.10/examples`.|
| `cppqed-doc-2.10` | (Only Ubuntu >= 13.10) Installs the documentation, also online [here](http://cppqed.sourceforge.net).|

You can uninstall our packages by calling

    $ sudo apt-get purge libcppqedelements-2.10-dev cppqedscripts-2.10 python-cpypyqed-2.10
    $ sudo apt-get autoremove
    $ sudo apt-add-repository --remove ppa:raimar-sandner/cppqed-development
    $ sudo apt-get update

### Arch

The C++QED package is in the [Arch User Repository (AUR)](https://aur.archlinux.org/packages/cppqed-git), you can install it for example with the `yaourt` package manager (`pacman` alternative which can install from AUR):

    $ yaourt -S cppqed-git

This will install all dependencies and compile C++QED.

# Installation from source {#installationguidefromsource}

## Requirements {#installationguiderequirements}

C++QEDv2 depends on a number of open-source libraries:

#### C++11 compatible compiler

C++QED uses some \ref cppelevenfeatures "C++11 features", therefore only modern C++ compilers are supported. These compilers are known to work:

* GCC, g++ >= 4.7
* LLVM/CLang, clang++ >= 3.1
* Intel compiler 14.0
* Apple LLVM version 5.0 (comes with the developer tools for Mac Os X Mavericks)

#### %CMake

%CMake is the build system of C++QED, version 2.8.9 or newer is needed.

#### Boost C++ libraries
The boost libraries provide indispensable extensions to the C++ standard, and are *de facto* standard by their own right. The framework depends on a number of them, the most notable ones being Fusion, Lambda, MPL, Operators, and Preprocessor. On many systems, (a selection of) the Boost libraries are available. They are packaged for Debian and Mac OS X. Alternatively, they can be downloaded and installed from [the main Boost portal](http://www.boost.org). Although version 1.46 is in principle sufficient for the framework, it is advisable to use as new a version as possible, since experience has shown that there can be clashes between old versions of Boost and new versions of compilers, especially when template metaprogramming is concerned.

#### GNU Scientific library (GSL)
GSL provides a very wide variety of numerical solutions in a solid object-oriented design (in C!). They are not used directly, but are wrapped into C++ classes and functions, so that they are easily replaced (e.g. if licensing problems arise). Packaged for Debian and Mac OS X, or can be downloaded from [the GSL homepage](http://www.gnu.org/software/gsl/). I haven't thoroughly determined the minimal version, but 1.8 is known to work.

The following two libraries are stable, but under more or less steady development.

#### Blitz++

\note Installation of Blitz++ can be omitted if the build system is later configure with `-DBUNDLED_BLITZ=On`. This will download and compile Blitz++ automatically when the framework is built.

Blitz++ provides the fundamental data structure. It hence performs a lot of numerics and lies at the absolute heart of the framework. Blitz++ lives up to its name, as it provides near-Fortran performance in spite of the very high-level abstractions used in the library. This is achieved by template metaprogramming, which was discovered in prototype during the development of this very library. More on [the Blitz++ homepage] (http://blitz.sourceforge.net).

C++QEDv2 needs some additional patches, therefore we host [our own version](http://sourceforge.net/p/cppqed/blitz/ci/default/tree/) which is backward compatible.

#### Optional: Flexible Library for Efficient Numerical Solutions (FLENS)

\note If FLENS should be used, the build system can be configured with `-DBUNDLED_FLENS=On`. This will download FLENS automatically when the framework is built.

[FLENS](http://flens.sourceforge.net) is a very impressive effort to wrap BLAS-LAPACK functions in a high-level C++ interface.
The use of FLENS is optional because only a very small tract of the framework depends on it.

#### Optional: Python

The \ref testsuite "testsuite" and the [CpypyQED](\cpypyqedMainPage) Python wrapper depends on Python >= 2.6 including NumPy and SciPy modules. If your system Python version is too old, we recommend the [Enthought Python Distribution](https://www.enthought.com/products/epd/free/).

## Installation of requirements {#installationguideinstallationrequirements}

On modern systems, the compiler, %CMake, Boost and GSL should be installed on system level through the package manager.
For example, on Ubuntu 13.10 it is sufficient to call:

    % sudo apt-get install cmake libboost-all-dev autoconf automake libtool libgsl0-dev \
      bzr git python-scipy

For Blitz++ and FLENS it is recommended to use the bundled versions (configure the framework with `-DBUNDLED_BLITZ=On` and `-DBUNDLED_FLENS=On`).

\note The rest of this section is only needed if any of the components is too old on your system, otherwise you can skip to the \ref cppqed_cmake "next section". By following these instructions, C++QED could successfully be installed on a quiet ancient CentOS 5.10.

Each requirement can be installed manually. The following sections give a hint how to do so. We will assume the manually installed components should be built inside `$BUILDDIR` (for example `$HOME/build`) and installed to the sub directory `$PREFIX` (for example `$HOME/local`). In order to find the components later in the actual build, add the following to the file `~/.bash_profile` and log in again:

    export BUILDDIR=$HOME/build
    export PREFIX=$HOME/local
    export PATH=$PREFIX/bin:$PATH
    export LD_LIBRARY_PATH=$PREFIX/lib64:$LD_LIBRARY_PATH
    export CMAKE_PREFIX_PATH=$PREFIX
    export NPROC=4

The environment variable `$NPROC` refers to the number of cores on your system (adjust appropriately), we will use that to speed up compilation.

Because some packages install libraries to `$PREFIX/lib`, some to `$PREFIX/lib64`, we make some preparations in `$PREFIX`:

    mkdir -p $PREFIX/lib64
    cd $PREFIX
    ln -s lib lib64

#### Compiler

If your compiler is too old, the following commands will install gcc-4.8.2 to `$PREFIX`:

    cd $BUILDDIR
    wget ftp://ftp.gnu.org/gnu/gcc/gcc-4.8.2/gcc-4.8.2.tar.bz2
    tar -xvjf gcc-4.8.2.tar.bz2
    cd gcc-4.8.2
    ./contrib/download_prerequisites
    cd ..
    mkdir gcc-4.8.2-build
    cd gcc-4.8.2-build
    $PWD/../gcc-4.8.2/configure --prefix=$PREFIX --enable-languages=c,c++ --disable-multilib
    make -j$NPROC
    make install

The compilation might take a long time. To use the new compiler, add the following to your `~/.bash_profile`:

    export CXX=$PREFIX/bin/g++
    export CC=$PREFIX/bin/gcc

#### %CMake

If your %CMake version is too old, the following commands will install cmake-2.8.12.2 to `$PREFIX`:

    cd $BUILDDIR
    wget http://www.cmake.org/files/v2.8/cmake-2.8.12.2.tar.gz
    tar -xvf cmake-2.8.12.2.tar.gz
    cd cmake-2.8.12.2
    ./bootstrap --prefix=$PREFIX
    make -j$NPROC
    make install

#### Python (optional, but needed for test suite)

If your Python version is too old, the following commands will install the Enthought Python Distribution "Canopy" to `$PREFIX`:

    wget https://store.enthought.com/downloads/canopy/rh5/64/free/ -O canopy.sh
    chmod a+x canopy.sh
    ./canopy.sh

After agreeing to the license, you can choose the installation directory (I will refer to this choice as `<CANOPY>`, replace it by the real directory in later commands). We will need the Python root installation directory later, you can display it by calling

    ls <CANOPY>/appdata/canopy*

It will display something like `<CANOPY>/appdata/canopy-1.3.0.1715.rh5-x86_64` maybe with a different version number (I will refer to this directory by `<CANOPY_ROOT>`, replace it by the real directory in later commands).

Call `canopy` for the first time, choose an environment directory and agree to set this python version as the default. Log out and log in again so that the changes can take effect.


#### Boost

If your Boost version is too old or if you installed Canopy Python above, you need to compile Boost. The following commands will install boost-1.55.0 to `$PREFIX`:

    cd $BUILDDIR
    wget http://sourceforge.net/projects/boost/files/boost/1.55.0/boost_1_55_0.tar.bz2/download -O boost_1_55_0.tar.bz2
    tar -xvf boost_1_55_0.tar.bz2
    cd boost_1_55_0

If you have installed Canopy Python earlier, call

    ./bootstrap.sh --prefix=$PREFIX --with-libraries=serialization,python,test --with-python-root=<CANOPY_ROOT>

and remember to replace `<CANOPY_ROOT>` by the real directory from above. Otherwise, if you are using the system version of Python, call

    ./bootstrap.sh --prefix=$PREFIX --with-libraries=serialization,python,test

In both cases, continue with

    ./b2 -j$NPROC
    ./b2 install


#### GSL

If your GSL version is too old, compile GSL manually. The following commands will install gsl-1.16 to `$PREFIX`:

    cd $BUILDDIR
    wget ftp://ftp.gnu.org/gnu/gsl/gsl-1.16.tar.gz
    tar -xvf gsl-1.16.tar.gz
    cd gsl-1.16
    ./configure --prefix=$PREFIX
    make -j$NPROC
    make install

## Building C++QED with CMake {#cppqed_cmake}

This page describes the C++QED build system. [CMake](http://www.cmake.org/) is a
cross-platform build system which checks the system for the location of required
and optional dependencies and enables/disables features based on what is found.
It is similar to the usual GNU autotools (aka `./configure; make; make install`)
and furthermore offers easy integration into popular IDE‘s, e.g.
[Eclipse](http://www.eclipse.org/),
[kdevelop](http://www.kdevelop.org/) or
[Xcode](http://developer.apple.com/xcode/).


### Short building guide {#cmake_short_guide}

This section is for the impatient, it describes the minimal steps to build the
framework. If %CMake cannot find installed requirements or if you wish to fine-tune
the build system, please also read the following sections.

Under the top level C++QED directory, create a build directory (e.g. `build`)

    mkdir build; cd build

Next, configure the build system by invoking %CMake

    cmake -DBUNDLED_BLITZ=On -DBUNDLED_FLENS=On -DCMAKE_BUILD_TYPE=<type> ..

where `<type>` is either `Debug` or `Release`. The default is `Release` if this
flag is missing. In debugging mode, a lot of runtime
checks are compiled into the framework, which come from Blitz++, FLENS, and
C++QED. Every time a new script is added it should be compiled and tested
in this way because this can detect a lot of errors. When we are absolutely
sure that everything is all right, for data collection we may compile in
release mode, where all the checks are omitted and optimisations are used,
making the programs about an order of magnitude faster. To switch the build type
create a new build directory or call `cmake` again with different `<type>`
parameter.

\warning
Maximum efficiency is achieved only if the build system is configured with

    cmake -DCMAKE_BUILD_TYPE=Release ..

To build and install the framework, type

    make -jN

where `N` is the number of parallel build processes. In order to install everything
on the system (which is not strictly necessary, you can run scripts out of the build
directory), call

    sudo make install

Every script in the scripts directory has its own target, so to compile only a single
script, it is possible to call

    make <script-name-without-extension>

### Fine tuning the build process {#cmake_fine_tuning}

If %CMake fails to locate a dependency, it will inform the user about the
disabled feature or fail if the library was required. If some libraries are
installed in non-standard locations, the user has to specify the directory with
`-DCMAKE_PREFIX_PATH="/some/path;/some/other/path"`. This will look for libraries
in `/some/path/lib`, `/some/other/path/lib` and for include files in
`/some/path/include`, `/some/other/path/include`. It is also possible to specify
additional library paths and include paths separately with `-DCMAKE_LIBRARY_PATH=`
and `-DCMAKE_INCLUDE_PATH=`.

\note
%CMake caches all detected library locations. If you want to be sure that the library
checks are performed, you have to delete the file `CMakeCache.txt` before calling `cmake` again.

By default, calling `make install` will install everything to the system directory
`/usr/local`. You can change this prefix by adding the option `-DCMAKE_INSTALL_PREFIX=<install path>`.

#### Variables which influence the build process ### {#cmake_build_defines}

There are a couple of %CMake options defined by C++QED itself, by which the build
process can be fine tuned:

* `-DCOMPILE_SCRIPTS=Off`: don't compile the example scripts
* `-DCOMPILE_CPYPYQED=Off`: don't compile Python modules
* `-DPYIO=Off`: don't compile the python module which allows reading and writing state files with Python.
* `-DSERIALIZATION=Off`: don't compile state vector serialization, even if boost serialization is found.
Note that without this feature, state files cannot be written to or read from disk.
* `-DFLENS=Off`: disable `FLENS` support, even if the library is detected.
* `-DREGISTRY=Off`: don't write information about the build trees to the
* `-DBUNDLED_FLENS=On/Off`: download FLENS automatically
* `-DBUNDLED_BLITZ=On/Off`: download and build Blitz++ automatically
[cmake package registry][cmake-registry].

See also the (possibly more complete) list of %CMake options in the CMake namespace ("Project options"), which is generated from the project
build system documentation.

#### How Cmake finds the C++QED components ### {#cmake_find_components}

The C++QED framework is organized in several sub-projects (core,elements,scripts,cpypyqed). For every component,
%Cmake has to figure out where the needed libraries are located. For example, the elements component needs
the libraries from core, the scripts and cpypyqed need the libraries from elements and core.

In the monolithic project layout (clone from the C++QED repository), all the component libraries are found
automatically in the common build directory, without user intervention. This section describes the case
when the different components are built in standalone projects.

##### Libraries are installed #####

If `make install` was called for C++QEDcore and the installation prefix is a standard system path
(e.g. `/usr`, `/usr/local`), then other sub-projects will find the core libraries automatically. If the libraries
are installed to a non-system prefix (e.g. `~/local`), this prefix can be added to the search path in
other sub-projects by calling %CMake with `-DCMAKE_PREFIX_PATH=~/local`. The same is true for C++QEDelements.

##### Libraries are not installed #####

It is possible to use C++QED libraries directly from the build directory where they were compiled. Libraries from
build directories are even found automatically most of the time, because these build directories are registered
in the [Cmake registry][cmake-registry]. You can override those locations by giving hints to %CMake where to
find the build directories:

* `-DCPPQED_DIR=<path>`: override path to the build directory of the core component
* `-DCPPQEDelements_DIR=<path>`: override path to the build directory of the elements component
* `-DCPPQEDcustomelements_DIR=<path>`: override path to a custom elements project named `customelements`

This is useful if you are using several build directories for different build configurations (Release/Debug).
In such a case, or if the libraries are intended to be installed with `make install`, you might consider to
call %CMake with `-DREGISTRY=Off` in order to not write information about build directories to the registry,
which might interfere with the installed libraries.

\note
To clear the %CMake registry, just delete `~/.cmake/packages` or individual directories there. Before running
%CMake again, also remove the `CMakeCache.txt` in the build directory to erase the cached locations.

### Building the documentation {#cmake_documentation}

The main target for building the documentation is `doc`. To build or rebuild only parts of the documentation,
the targets `<component>_doc` exist, where `<componenet>` is one of

- `cppqed`: Overall C++QED documentation, build system, test suite.
- `core`: C++QED core API
- `elements`: C++QED elements API
- `cpypyqed`: Python wrapper documentation (in Sphinx)

The targets mentioned above are only available if `doxygen` and `dot` is installed on the system. The
`cpypyqed` documentation additionally depends on [Sphinx](http://sphinx-doc.org/) and the
[doxylink](https://pypi.python.org/pypi/sphinxcontrib-doxylink) extension.

### Integrated Development Environments (IDE) {#cmake_ide}

Using an IDE can make the development of scripts and the framework much more convenient. %CMake has generators for project files of several popular IDE's, check the %CMake man page to see which are supported on your system. For example, to generate an Eclipse CDT project, call ::

    cmake -G "Eclipse CDT4 - Unix Makefiles" [other options]

in the C++QED root directory. This will generate the necessary Eclipse project files, afterwards you can import the project into Eclipse by choosing "File->Import->Existing Project into workspace". See the [CMake wiki] (http://www.cmake.org/Wiki/Eclipse_CDT4_Generator) for details.

For kdevelop4 you can just open the top level ``CMakeLists.txt`` file as a project.

Under Mac OS X the ``-G Xcode`` generator will create a project which can be opened by Xcode.

# Mac OS X {#installationguidemacos}

Compilation of C++QED under Mac OS X 10.9 is straightforward once you have installed XCode and the prerequisite libraries. If you do not intend to use XCode for development, it is sufficient to install only the development command line tools.

XCode can be obtained from the [Apple Developer site] (https://developer.apple.com/downloads) (you need to register an Apple ID).

## Prerequisite Libraries

The easiest way to satisfy the dependencies is to use [Homebrew] (http://mxcl.github.com/homebrew), a package manager for Mac Os X.
Unfortunately Boost cannot be compiled with Homebrew and Apple's Clang compiler in Mac Os X 10.8. It might still be possible to install a different compiler and succeed, but this is not tested.


1. Install [Homebrew](http://mxcl.github.com/homebrew) and don't forget to run `brew update` afterwards.

2. Install prerequisites:

       brew install libtool autoconf automake gsl cmake mercurial

3. Install Boost (this works in Mac Os X 10.9 but not 10.8):

       brew install boost --c++11

4. If you want to use the Python parts of C++QED (testsuite and Python wrapper), you also need NumPy and Scipy for Mac Os X:

       brew tap Homebrew/python
       brew install scipy

If you installed Scipy, follow the instructions to add `/usr/local/lib/python2.7/site-packages` to your `PYTHONPATH`. After installing the required libraries, continue \ref cppqed_cmake "as usual" building the project with %CMake.

[cmake-registry]: http://www.cmake.org/Wiki/CMake/Tutorials/Package_Registry
