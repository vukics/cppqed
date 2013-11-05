* boost::tuples => std::tuples --- deprecate DynamicsBase::RealFreqs,ComplexFreqs FREQS, etc.

* Make headers self-contained & minimal

* MultiLevel probably severely outdated
  
* C++11 features currently used:

  * variadic templates – variadic parameter lists
  * template aliases
  * new function declaration syntax
  * auto keyword => new for syntax
  * rvalue references
  * new initialization syntax + initialization lists
  
  => g++ >= 4.7 ; clang++ >= 3.1

  * To be adopted
  
    * delegating constructors 4.7,3.0
    * inherited constructors 4.8,3.3
   (* rvalue references for StateVectorLow&)
    * initialization of static constants within the class
    * forward-declared enums 4.6,3.1
    * Non-static data member initializers
    * lambda 4.5,3.1
    * decltype 4.8,2.9


# Code organization + packaging

  * DONE hierarchical: core <- elements <- scripts ==> git sub-repositories (hierarchy within elements: frees <- interactions)
  
  * ubuntu packages should also reflect this three-fold structure (separate dev & bin, dev depends on lower-level dev & bin, bin depends only on lower-level bin)
  
  * Not needed, incorporated into core. (4th sub-repository + package for build infrastructure)
  
  * DONE cmake should observe structural dependencies (e.g. quantumdata is not allowed to use anything from quantumtrajectory, only the other way round) – can be achieved by specific cmake files in these directories


# Build system

* DONE make a stronger connection between scripts and core in cmake, scripts should basically inherit the same compile-configuration as was used for core

* DONE quantumdata may include only from utils; quantumtrajectory may include quantumdata; structure may include both quantumtrajectory & quantumdata – TridiagonalHamiltonian should rather go to quantumoperator?

* Git Branch & Commit details to be displayed by scripts on call with --version (note, this information must be captured @ compile time, and not simply @ configuration time, that is, the query must be performed by make and not by cmake)

* Deprecate Boost.Build

  
# New release

  * Documentation
  
    * Finish basic API documentation (doxygen branch)
    * Update structure-bundle documentation (highlight c++11 features – give a name, link, etc.)
    * User guide
      * trajectory-state serialization
      * trajectory logging
     (* heavier use of smart pointers)
    * Installation ???
    
    * Python interface
    * pycppqed ???
    
  * Form of distribution (less reliance on packages, more integration?)
  
    * Switch to new FLENS (header only) & distribute it together with the project
    * Try to turn Blitz++ also into a header-only dependence
    * Boost – distribute together with the project those parts that are needed, if it gets installed, then enable serialization
      * or, try to use Boost’s cmake to integrate its compilation
    
  * How much Python to include??? (perhaps only proofs of principle)


# Python i/o

  * state-saving saves 0 for dtTry, and then when the initial-condition file is read in, the actual dtTry should be calculated through the highest frequency
  
  * in trajectory::Trajectory, introduce a new virtual (setInitalDtTry) – this can be implemented by quantumtrajectories by calling the highestFrequency function of the system (put it in a base class QuantumTrajectory common to Master & MCWF_Trajectory)
  
  * make SerializationMetaData a data member of trajectory::Trajectory, then, the virtual interface can be as it was before (no writeMeta & no extra argument for readState is needed)
  
  * Rationale: each archive is self-contained, so that it contains its own metadata. Then, the state files can eventually even be split into several archives, each containing a single self-contained state.
  
