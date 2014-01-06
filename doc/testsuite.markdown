# Test suite {#testsuite}

\tableofcontents

# Short user guide  {#shortguide}

After [cloning](#cloning) the C++QED repository, there are several ways to run the test suite
fully or partially.

### `make check*` targets

Make targets have the advantage that they build everything the tests depend on automatically.
They will run all tests with verbose output and offer the possibility to run fine grained groups
of tests. Example usage:

    # Run the full test suite
    $ make check

    # Perform the tests of the "continue" directory
    $ make check_continue

    # Perform only a single test by name
    $ make fail_CompositeError2

Success is indicated by the make target completing without errors. Parallelization is possible to
a limited extent by running `make -jN <target>`. Only the different make targets, not the individual tests
will be run in parallel. This method is suited best for usage from within kdevelop, for example while
debugging some feature one can add the relevant test to the build sequence and just hit F8, which will
recompile the dependencies and run the test.

### `make test` target

This target is added by cmake automatically. It has no dependency checking, which means the framework has to be
built in full before calling it. It features a nicer and compact output (no verbose messages from the tests) and
a summary of all completed and failed tests. No parallelization is supported.

### `ctest`

`ctest` is the tool that manages the tests, it is internally called by all the make targets. It can also be called
directly. Example usage:

    # Run the full test suite on 2 cores
    $ ctest -j2

    # Run all tests matching run_ on 2 cores
    $ ctest -j2 -R run_

Again there is no dependency checking, the framework has to be built in full before calling `ctest`. However,
tests depending on other tests will be run in the right order. Running `ctest` directly has the best parallelization
support, because individual tests are run in parallel, making this the fastest method. Verbose output from
tests can be enabled by the `-V` switch.

# Test suite layout {#testsuite_layout}

## The `cmake` part

Tests are grouped into sub-directories of `Testing`. To add a new directory, the `CMakeLists.txt` has to include
the statement

    testdir(<dirname> [dependency1 dependency2 ...])

Among other things, this will create the `check_<dirname>` target. The optional dependencies are `cmake` targets
which are required to be built before all of the tests in this directory can be run. For example, if one of the tests
requires `1particle1mode` one has to add `1particle1mode` to the dependencies.

Each test can use one of two test drivers: the [Python](#testsuite_python) or [boost](#testsuite_boost) test driver.
The method how to add a new test depends on the test driver. In any case, before a new test can be used, the
project has to be re-configured with cmake.

## The Python test driver {#testsuite_python}

To add a new test, one has to call **delcaretest**, **add_test** and optionally **set_tests_properties**.
Here is an example how to add a new test in the `CMakeLists.txt`:

    declaretest(<unique test name>)
    add_test(NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT} <parameters>)
    # the following line is optional, adds a dependency:
    set_tests_properties(${TESTNAME} PROPERTIES DEPENDS <some other test>)

The variable `${TESTNAME}` is set to the current test name by `declaretest`. This will call the
Python script **testdriver.py** as test driver. The variable `${TESTSCRIPT}` points to `testdriver.py`
and already adds some needed parameters to pass in the build configuration (debug or release) and the
configuration file. The `<parameters>` are:

* `--test=TEST` (required): the name of the test, usually passed in as `--test=${TESTNAME}`. This name also defines
the section in the configuration file.
* `--testclass=TESTCLASS` (required): the name of the python class which defines the mode of operation of the
test driver
* `--script=SCRIPT`: the script to run by the test driver or the target to compile. Typically passed in for example as
`--script=$<TARGET_FILE:1particle1mode>` which will be evaluated by cmake to the full path of either `1particle1mode`
(release configuration) or `1particle1mode_d` (debug configuration).
* `--error=ERROR`: only needed for the `CompileFail` testclass. The test will be considered successful if `ERROR`
is found in the output.

### Configuration files {#testsuite_py_configuration}

The Python test driver reads in a global configuration file (**Testing/testdriver.conf.in**, pre-processed by cmake) and
a per-directory configuration file (**testdriver.conf**). In the `[Scripts]` section of the global configuration file,
one can define options which will be in effect for all scripts:

    [Scripts]
    opts=--precision 6 --epsAbs 1e-5 --logLevel 2

In the per-directory configuration files, there is a section for every test name. The available configuration keys
of the sections depend on the used test class and are documented [there](#testsuite_py_classes).
One can import the options from another section by an `import` statement:

    [run_PTLA_CPPQED_dc]
    ...

    [run_PTLA_CPPQED_dt]
    import=run_PTLA_CPPQED_dc
    ...

Note that configuration keys which are also defined in the current section will not be overwritten by an import.

### Test classes {#testsuite_py_classes}

This section documents all the available test classes, which can be used with the
`--testclass=TESTCLASS` parameter. Some classes inherit from others, in this case the configuration
options from the base class can be used in the derived one.

#### OutputManager

This cannot be used as a test class directly, but it serves as base to other test classes.
Configuration options:

    # List of run modes to use for this test
    runmodes=single,ensemble,master

#### Runner (derives from OutputManager)

Runs a script and succeeds if the script does.
Configuration options:

    # The options used for running the scripts, multiple keys can be given if they match opts*
    opts=--etat 8 --sdf 3
    opts1=--dc 0 --Dt 0.1 --NDt 10

    # The options which are additionally added for the specified run mode
    # multiple keys can be given if they match <runmode>*
    single=...
    ensemble=...
    master=...


#### Verifier (derives from OutputManager)

Verifies the output of a script 'this' to an expected output or the output of some other test run 'other'.
By default, 'this' is the current test. Configuration options:

    # Verify that both trajectories are exactly equal (default if this key is missing)
    verify=full
    # OR verify that the last outcome of the simulation is equal (e.g. timesteps may differ)
    verify=outcome

    # Testname of 'this', by default the current test if missing
    this=some_test
    # Testname of 'other', by default the results from the directory of expected results
    other=some_other_test

If `this=some_test` is specified, it is probably also a good idea to `import=some_test` to keep
the runmodes in sync. The directory of expected results is `Testing/expected`, it is kept
under version control so that changes in the output of the scripts are noticed.

#### VerifiedRunner (derives from Runner, Verifier)

Combines the functionality of `Runner` and `Verifier` to a single test.

#### Continuer (derives from Runner)

Runs a script, then continues the trajectory. Succeeds if both runs do. Use `Verifier` to compare the outcome
to another test or an expected output. Configuration options:

    # Options for the first run and the second run, respectively
    firstrun=...
    secondrun=...

#### CompileFail

This test succeeds if the compilation of a target fails in an expected way. For this test class, the
`--script=SCRIPT` parameter is not the path to a binary but a cmake target. The expected failure message
is passed to the test driver by the `--error=ERROR` command line parameter. The cmake targets which should
fail to compile are typically defined in `Testing/compilefail`.

## The boost test driver {#testsuite_boost}

Boost comes with its own test driver, which can be easily used with this test suite. All boost tests are
declared in the directory `Testing/boost`.

To add a new test, a C++ source file, e.g. `new_test.cc`, has to be added to `Testing/boost` which implements
the test:

    #include <boost/test/unit_test.hpp>
    ...
    BOOST_AUTO_TEST_CASE( UNIQUE_TEST_NAME )
    {
    ...
    // for example: BOOST_CHECK(...)
    }

In `Testing/boost/CMakeLists.txt`, one has to add `new_test.cc` to `BOOST_TEST_SOURCES` and `UNIQUE_TEST_NAME`
to `BOOST_TESTS`:

    ...
    set(BOOST_TEST_SOURCES ... new_test.cc)
    ...
    set(BOOST_TESTS ... UNIQUE_TEST_NAME)
    ...

Afterwards, the new test will be integrated into the test suite.

# Cloning the C++QED repository {#cloning}

The test suite is available in the Development branch of the C++QED repository, which
contains the whole framework. It can be cloned for read-write access with

    $ git clone --recursive -b Development ssh://<username>@git.code.sf.net/p/cppqed/cppqed C++QED

where <username> is the SourceForge user name. Note that the git submodules are cloned from read-only
URLs, in each submodule you have to set the URL to a read-write one, for example:

    $ cd C++QED/CPPQEDcore
    $ git remote set-url origin ssh://<username>@git.code.sf.net/p/cppqed/cppqed_core

and likewise for the other submodules.