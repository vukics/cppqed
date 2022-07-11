* [An introductory talk](https://youtu.be/yunR73LNJ1M) given at GPU Day 2019
* [Video tutorial part 1](https://youtu.be/Ozj4sDNdlmM) — How to download, configure and build C++QED, and run the testsuite on a fresh Ubuntu20.04 installation. The commands used in the tutorial can be found in `helperscripts/cloudScript.sh`. They are as follows:
```
# Starting from fresh Ubuntu 20.04 LTS installation
# Install prerequisites:
sudo apt install libgsl-dev libboost-all-dev libeigen3-dev clang cmake cmake-curses-gui git python3-scipy # ipython3
# Clone repo + submodules:
git clone --recurse-submodules --remote-submodules https://github.com/vukics/cppqed.git; cd cppqed
# Checkout Ubuntu branch: 
git checkout Ubuntu20.04LTS
mkdir build; cd build
# Configure build:
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=/usr/bin/clang++
# Build everything needed for tests:
make -j4 boostTester fewer_scripts compileTests cpypyqed 
# Tests can be run with
# ctest -j4 # from the main build folder
```

C++QED is a C++/Python framework for simulating open quantum dynamics. It allows users to build arbitrarily complex interacting quantum systems from elementary free subsystems and interactions, and simulate their time evolution with a number of available time-evolution drivers.

Earlier, C++QED was located at sourceforge: http://cppqed.sf.net. The (at the moment somewhat outdated) documentation can still be found there.
* [User guide](http://cppqed.sourceforge.net/cppqed/html/userguide.html)

![Example trajectories generated for a recent usecase, cf. https://quantum-journal.org/papers/q-2019-06-03-150/](/doc/exampleTrajectories.png)
