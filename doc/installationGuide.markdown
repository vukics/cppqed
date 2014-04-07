Installation guide {#installationguide}
==================

# Quick start: Installation from packages (Ubuntu, Arch Linux) {#installationguidequickstart}

Binary packages have been prepared for Arch Linux and [Ubuntu](http://www.launchpad.net/~raimar-sandner/cppqed),
which might also work on Debian (not tested).  Installation from binary packages is recommended for users
who want to try C++QED or use it exclusively on its highest level, i.e. writing
scripts which use elements and interactions already implemented in C++QED. If
you are not using Ubuntu or want to develop your own elements and interactions
you have to install C++QED from \ref installationguidefromsource "source".

### Ubuntu

All packages are compiled for Ubuntu Precise (12.04 LTS) and Saucy (13.10). You can choose between the latest release and daily builds. For the release:

    $ sudo add-apt-repository ppa:raimar-sandner/cppqed-development
    $ sudo apt-get update
    $ sudo apt-get install libcppqedelements-2.10-dev cppqedscripts python-cpypyqed

For the daily builds:

    $ sudo add-apt-repository ppa:raimar-sandner/cppqed-daily
    $ sudo apt-get update
    $ sudo apt-get install libcppqedelements-daily-dev cppqedscripts python-cpypyqed

The current documentation can be installed with

    $ sudo apt-get install cppqed-doc

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

## Installation of requirements

On modern systems, the compiler, %CMake, Boost and GSL should be installed on system level through the package manager.
For example, on Ubuntu 13.10 it is sufficient to call:

    % sudo apt-get install cmake libboost-all-dev autoconf automake libtool libgsl0-dev \
      bzr git python-scipy

For Blitz++ and FLENS it is recommended to use the bundled versions (configure the framework with `-DBUNDLED_BLITZ=On` and `-DBUNDLED_FLENS=On`). You can then skip the rest of this section.

If any of the components is too old on your system, you can install it manually. The following sections give a hint how to do so. We will assume the manually installed components should be built inside `$BUILDDIR` (for example `$HOME/build`) and installed to the sub directory `$PREFIX` (for example `$HOME/local`). In order to find the components, add the following to the file `~/.bash_profile`:

    export BUILDDIR=$HOME/build
    export PREFIX=$HOME/local
    export PATH=$PREFIX/bin:$PATH
    export LD_LIBRARY_PATH=$PREFIX/lib:$PREFIX/lib64:$LD_LIBRARY_PATH
    export CMAKE_PREFIX_PATH=$PREFIX
    export NPROC=4

The environment variable `$NPROC` refers to the number of cores on your system (adjust appropriately), we will use that to speed up compilation.

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
    $PWD/../gcc-4.8.2/configure --prefix=$PREFIX --enable-languages=c,c++ --disable-multilib --libdir=$PREFIX/lib64
    make -j$NPROC
    make install

This might take a long time.

#### %CMake

If your %CMake version is too old, the following commands will install cmake-2.8.12.2 to `$PREFIX`:

    cd $BUILDDIR
    wget http://www.cmake.org/files/v2.8/cmake-2.8.12.2.tar.gz
    tar -xvf cmake-2.8.12.2.tar.gz
    cd cmake-2.8.12.2
    ./bootstrap --prefix=$PREFIX
    make -j$NPROC
    make install

Note that you probably also have to compile Boost and GSL with the new compiler to have binary compatible libraries.

#### Python (optional, but needed for test suite)

If your Python version is too old, the following commands will install the Enthought Python Distribution to `$PREFIX`:

    wget https://www.enthought.com/downloads/canopy/rh5-64/free/ -O canopy.sh
    chmod a+x canopy.sh
    ./canopy.sh

After agreeing to the license, you can choose the installation directory. Call `canopy` and choose the environment directory (I will refer to this choice as `<CANOPY>`). Add the following to your `$CMAKE_PREFIX_PATH` variable in `.bash_profile`:

    export CMAKE_PREFIX_PATH=<CANOPY>/System:$PREFIX

Log out and log in again so that the changes can take effect.


#### Boost

If your Boost version is too old or you installed a new compiler, the following commands will install boost-1.55.0 to `$PREFIX`:

    cd $BUILDDIR
    wget http://sourceforge.net/projects/boost/files/boost/1.55.0/boost_1_55_0.tar.bz2/download -O boost_1_55_0.tar.bz2
    tar -xvf boost_1_55_0.tar.bz2
    cd boost_1_55_0
    ./bootstrap.sh --prefix=$PREFIX --with-libraries=serialization,python,test
    ./b2 -j$NPROC
    ./b2 install


#### GSL

If your GSL version is too old or you installed a new compiler, the following commands will install gsl-1.16 to `$PREFIX`:

    cd $BUILDDIR
    wget ftp://ftp.gnu.org/gnu/gsl/gsl-1.16.tar.gz
    tar -xvf gsl-1.16.tar.gz
    cd gsl-1.16
    ./configure --prefix=$PREFIX
    make -j$NPROC
    make install


