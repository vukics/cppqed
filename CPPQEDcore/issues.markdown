Issues
======

* \b DONE boost::tuples => std::tuples --- deprecate DynamicsBase::RealFreqs,ComplexFreqs FREQS, etc.
* (re)consider use of Boost.Math instead of GSL throughout

* (re)consider use of Boost.MultiArray instead of Blitz

* define Ptr types outside classes as template aliases eg QuantumSystem::Ptr => QuantumSystemPtr


* Make headers self-contained & minimal

* \b DONE MultiLevel probably severely outdated
  

# Code organization + packaging

  * \b DONE hierarchical: core <- elements <- scripts ==> git sub-repositories (hierarchy within elements: frees <- interactions)
  
  * ubuntu packages should also reflect this three-fold structure (separate dev & bin, dev depends on lower-level dev & bin, bin depends only on lower-level bin)
  
  * <b>Not needed, incorporated into core.</b> (4th sub-repository + package for build infrastructure)
  
  * \b DONE cmake should observe structural dependencies (e.g. quantumdata is not allowed to use anything from quantumtrajectory, only the other way round) – can be achieved by specific cmake files in these directories


# Build system

* \b DONE make a stronger connection between scripts and core in cmake, scripts should basically inherit the same compile-configuration as was used for core

* \b DONE quantumdata may include only from utils; quantumtrajectory may include quantumdata; structure may include both quantumtrajectory & quantumdata – TridiagonalHamiltonian should rather go to quantumoperator?

* \b DONE Git Branch & Commit details to be displayed by scripts on call with --version (note, this information must be captured @ compile time, and not simply @ configuration time, that is, the query must be performed by make and not by cmake)

* Deprecate Boost.Build

  
# New release

  * Documentation
  
    * Finish basic API documentation (doxygen branch)
    * \b DONE Update structure-bundle documentation (highlight c++11 features – give a name, link, etc.)
    * User guide
      * trajectory-state serialization
      * trajectory logging (especially EnsembleMCWF)
      * time averaging
      * version information
      * (heavier use of smart pointers)
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
