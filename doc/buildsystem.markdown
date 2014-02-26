# The cmake build system # {#cppqed_cmake}

\tableofcontents

This page describes the C++QED build system. [CMake](http://www.cmake.org/) is a
cross-platform build system which checks the system for the location of required
and optional dependencies and enables/disables features based on what is found.
It is similar to the usual GNU autotools (aka `./configure; make; make install`)
and furthermore offers easy integration into popular IDE‘s, e.g.
[Eclipse](http://www.eclipse.org/),
[kdevelop](http://www.kdevelop.org/) or
[Xcode](http://developer.apple.com/xcode/).

C++QED requires CMake version 2.8.9 or later.

# Building the framework {#cmake_building}

## Short user guide {#cmake_short_guide}

This section is for the impatient, it describes the minimal steps to build the
framework. If CMake cannot find installed requirements or if you wish to fine-tune
the build system, please also read the following sections.

Under the top level C++QED directory (the monolithic clone of all project components),
create a build directory (e.g. `build`)

    mkdir build; cd build

Next, configure the build system by invoking CMake

    cmake -DCMAKE_BUILD_TYPE=<type> ..

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

## Fine tuning the build process {#cmake_fine_tuning}

If CMake fails to locate a dependency, it will inform the user about the
disabled feature or fail if the library was required. If some libraries are
installed in non-standard locations, the user has to specify the directory with
`-DCMAKE_PREFIX_PATH="/some/path;/some/other/path"`. This will look for libraries
in `/some/path/lib`, `/some/other/path/lib` and for include files in
`/some/path/include`, `/some/other/path/include`. It is also possible to specify
additional library paths and include paths separately with `-DCMAKE_LIBRARY_PATH=`
and `-DCMAKE_INCLUDE_PATH=`.

\note
CMake caches all detected library locations. If you want to be sure that the library
checks are performed, you have to delete the file `CMakeCache.txt` before calling `cmake` again.

By default, calling `make install` will install everything to the system directory
`/usr/local`. You can change this prefix by adding the option `-DCMAKE_INSTALL_PREFIX=<install path>`.

### Variables which influence the build process ### {#cmake_build_defines}

There are a couple of CMake options defined by C++QED itself, by which the build
process can be fine tuned:

* `-DCOMPILE_SCRIPTS=Off`: don't compile the example scripts
* `-DCOMPILE_CPYPYQED=Off`: don't compile Python modules
* `-DPYIO=Off`: don't compile the python module which allows reading and writing state files with Python.
* `-DSERIALIZATION=Off`: don't compile state vector serialization, even if boost serialization is found.
Note that without this feature, state files cannot be written to or read from disk.
* `-DFLENS=Off`: disable `FLENS` support, even if the library is detected.
* `-DREGISTRY=Off`: don't write information about the build trees to the
[cmake package registry][cmake-registry].
* `-DEXAMPLES=Off`: con't compile some examples in the C++QEDcore repository.

See also the (possibly more complete) list of %CMake options in the CMake namespace ("Project options"), which is generated from the project
build system documentation.

### How Cmake finds the C++QED components ### {#cmake_find_components}

The C++QED framework is organized in several sub-projects (core,elements,scripts,cpypyqed). For every component,
Cmake has to figure out where the needed libraries are located. For example, the elements component needs
the libraries from core, the scripts and cpypyqed need the libraries from elements and core.

In the monolithic project layout (clone from the C++QED repository), all the component libraries are found
automatically in the common build directory, without user intervention. This section describes the case
when the different components are built in standalone projects (e.g. clones
from the repositories C++QEDcore, C++QEDelements etc.).

#### Libraries are installed ####

If `make install` was called for C++QEDcore and the installation prefix is a standard system path
(e.g. `/usr`, `/usr/local`), then other sub-projects will find the core libraries automatically. If the libraries
are installed to a non-system prefix (e.g. `~/local`), this prefix can be added to the search path in
other sub-projects by calling Cmake with `-DCMAKE_PREFIX_PATH=~/local`. The same is true for C++QEDelements.

#### Libraries are not installed ####

It is possible to use C++QED libraries directly from the build directory where they were compiled. Libraries from
build directories are even found automatically most of the time, because these build directories are registered
in the [Cmake registry][cmake-registry]. You can override those locations by giving hints to Cmake where to
find the build directories:

* `-DCPPQED_DIR=<path>`: override path to the build directory of the core component
* `-DCPPQEDelements_DIR=<path>`: override path to the build directory of the elements component
* `-DCPPQEDcustomelements_DIR=<path>`: override path to a custom elements project named `customelements`

This is useful if you are using several build directories for different build configurations (Release/Debug).
In such a case, or if the libraries are intended to be installed with `make install`, you might consider to
call CMake with `-DREGISTRY=Off` in order to not write information about build directories to the registry,
which might interfere with the installed libraries.

\note
To clear the CMake registry, just delete `~/.cmake/packages` or individual directories there. Before running
Cmake again, also remove the `CMakeCache.txt` in the build directory to erase the cached locations.

# Building the documentation {#cmake_documentation}

The main target for building the documentation is `doc`. To build or rebuild only parts of the documentation,
the targets `<component>_doc` exist, where `<componenet>` is one of

- `cppqed`: Overall C++QED documentation, build system, test suite.
- `core`: C++QED core API
- `elements`: C++QED elements API
- `cpypyqed`: Python wrapper documentation (in Sphinx)

The targets mentioned above are only available if `doxygen` and `dot` is installed on the system. The
`cpypyqed` documentation additionally depends on [Sphinx](http://sphinx-doc.org/) and the
[doxylink](https://pypi.python.org/pypi/sphinxcontrib-doxylink) extension.

The documentation can only be installed if the project is built in monolithic mode. In this case the
destination directory is set to `/usr/share/doc/cppqed-doc-${CPPQED_ID}`. In standalone mode the documentation
can still be built (it resides in `${CMAKE_BINARY_DIRECTORY}/doc`) and also inter-project links should work,
however the documentation cannot be moved to another place and thus cannot be installed system wide.

[cmake-registry]: http://www.cmake.org/Wiki/CMake/Tutorials/Package_Registry