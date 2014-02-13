Issues
======

* (re)consider use of Boost.Math instead of GSL throughout

* (re)consider use of Boost.MultiArray instead of Blitz

* define Ptr types outside classes as template aliases eg QuantumSystem::Ptr => QuantumSystemPtr

* Make headers self-contained & minimal

# New release

  * Documentation
  
    * Finish basic API documentation (doxygen branch)
    * User guide
      * trajectory-state serialization
      * trajectory logging (especially EnsembleMCWF)
      * time averaging
      * version information
      * (heavier use of smart pointers)

    * pycppqed ???
    
  * Form of distribution (less reliance on packages, more integration?)

    * Switch to new FLENS (header only) & distribute it together with the project
    * Try to turn Blitz++ also into a header-only dependence
    * Boost – distribute together with the project those parts that are needed, if it gets installed, then enable serialization
      * or, try to use Boost’s cmake to integrate its compilation
    