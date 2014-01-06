# Test suite {#testsuite}

\tableofcontents

# Short user guide  {#shortguide}

After [cloning](#cloning) the C++QED repository, there are several ways to run the test suite
in full or partial.

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

To add a new test, one has to call **delcaretest**, **add_test** and optionally **set_tests_properties**.
Here is an example how to add a new test in the `CMakeLists.txt`:

    declaretest(<unique test name>)
    add_test(NAME ${TESTNAME} COMMAND <command with parameters>)
    # the following line is optional, adds a dependency:
    set_tests_properties(${TESTNAME} PROPERTIES DEPENDS <some other test>)

The command can be anything which indicates success by a return value of 0, typically this will be the
[Python](#testsuite_python) or boost test driver. The variable `${TESTNAME}` is set to the current test
name by `declaretest`.

## The Python part {#testsuite_python}

The Python script **testdriver.py** can serve as a command to `add_test`, the syntax is:

    add_test(NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT} <parameters>)

# Cloning the C++QED repository {#cloning}

The test suite is available in the Development branch of the C++QED repository, which
contains the whole framework. It can be cloned for read-write access with

    $ git clone --recursive -b Development ssh://<username>@git.code.sf.net/p/cppqed/cppqed C++QED

where <username> is the SourceForge user name. Note that the git submodules are cloned from read-only
URLs, in each submodule you have to set the URL to a read-write one, for example:

    $ cd C++QED/CPPQEDcore
    $ git remote set-url origin ssh://<username>@git.code.sf.net/p/cppqed/cppqed_core

and likewise for the other submodules.