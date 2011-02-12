-----------------
Design rationales
-----------------

* All information available at compile time (in particular, the complete layout of the simulated system) should be processed as such, with the help of template metaprogramming.

  * Note that whether the system is Hamiltonian, Liouvillean, etc is generally NOT known at compile time since it depends on parameters. (e.g. if eta is nonzero, Lossy mode is Hamiltonian, but otherwise not)

  * The fact whether the system is NonOrthogonal, however, is known at compile time, so this can be treated eg. in Trajectories with the help of a boolean template parameter.


* Interactions should not bother with shifting frequencies as in version 1, because these frequencies appear at so many places that it is hardly possible to track them all. Instead, the shift of frequencies, when required at all, should be performed on the highest level, that is, in scripts like ``1particle1mode``.

* The basic data structures should be built on Blitz++ --- memory block of doubles, state vectors, complex matrices, etc. 
  
* Every non-trivial numerics should be done by GSL with bindings developed in the C++Utils. Alternative implementations are possible. Basic linear algebra should be tackled using the Blitz++ operators, while more involved ones using LAPACK.

* Rely heavily on STL, TR1, Boost both in design and implementation, however, use Boost instead of TR1: TR1 support for gcc is scarce, and what is implemented doesn't interoperate well with other parts of Boost. (E.g. boost::lambda and tr1::function)

* Elements should declare their parameter structs in their own namespaces, and split them into separate header & implementation files. Eg. ParsParticle -> particle::Pars.


* Elementary systems can select from which class to derive (Hamiltonian etc.), because here we literally DON'T KNOW how to implement Hamiltonian if they are not indeed Hamiltonian. Composite derive from all (we know how to implement: iterate over subsystems AND dummies). This means we have to put up with a lot of ``dynamic_cast``s as in ``MCWF_Trajectory``, and the constant queries whether the system at hand features this and that interface. Alternative: see below.

* In object hierarchies make the objects noncopyable, but allow them to be Cloneable. Then they can be conveniently manipulated by ptr_containers.


Future
^^^^^^
 
* Somehow try to separate the DEFINITION of the layout of systems from the functionalities like ``Hamiltonian``. Then ``Composite`` could be only defining the layout, but not derived from ``Hamiltonian``.

* Implement a more elaborate set of trajectories like `ExactTrajectory`, `HamiltonianTrajectory` AND `MCWF_Trajectory`
