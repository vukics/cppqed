Rationales {#rationales}
==========

\tableofcontents

\par Abstract
Failure to supply contemporaneous rationale for design decisions is a major defect in many software projects. Lack of accurate rationale causes issues to be revisited endlessly,
causes maintenance bugs when a maintainer changes something without realizing it was done a certain way for some purpose, and shortens the useful lifetime of software.
**Rationale is fairly easy to provide at the time decisions are made, but very hard to accurately recover even a short time later.**


Design rationales {#designrationales}
=================

* For polymorphic classes, use the [non-virtual-interface idiom](http://www.gotw.ca/publications/mill18.htm).
Down the hierarchy, virtual functions should almost always be `private` — `protected` only in the case when a derived class uses the class's implementation.

* All information available at compile time (in particular, the complete layout of the simulated system) should be processed as such, with the help of template metaprogramming.

  * Note that whether the system is \link structure::Hamiltonian Hamiltonian\endlink, \link structure::Liouvillean Liouvillean\endlink, etc.
is generally *not* known at compile time since it depends on parameters. (E.g. if there is pumping, Mode is \link structure::Hamiltonian Hamiltonian\endlink, but otherwise not.)

  * The fact whether the system is NonOrthogonal, however, is known at compile time, so this can be treated e.g. in trajectories with the help of a boolean template parameter.

* Interactions should not bother with shifting frequencies as in C++QEDv1, because these frequencies appear at so many places that it is hardly possible to track them all.
Instead, the shift of frequencies, when required at all, should be performed on the highest level, that is, in scripts like `1particle1mode`.

* The basic data structures should be built on [Blitz++](http://blitz.sourceforge.net/) — memory block of doubles, state vectors, complex matrices, etc.
\todo Consider switching to \refBoost{Boost.MultiArray,multi_array}. This would also provide const views which are not present in Blitz++.

* All non-trivial numerics should be performed by [GSL](http://www.gnu.org/software/gsl/) with bindings developed in `utils`.
Alternative implementations are possible, consider in particular [Boost.Math](http://www.boost.org/doc/libs/?view=category_Math).

  Basic linear algebra should be tackled using the Blitz++ operators, while more involved ones using [FLENS/LAPACK](http://flens.sourceforge.net/‎). \todo Consider switching to [Eigen](http://eigen.tuxfamily.org/).

* Rely heavily on STL and Boost both in design and implementation, however, use Boost instead of TR1:
TR1 support for gcc appears scarce, and what is implemented doesn't interoperate well with other parts of Boost. (E.g. `boost::lambda` and `tr1::function`)
\todo Revisit this issue (might not exist anymore in C++11)

* Elements should declare their parameter `struct`s in their own namespaces, and split them into separate header & implementation files. Eg. ParsParticle -> particle::Pars.

* In object hierarchies make the objects `noncopyable`, but allow them to be \refBoost{cloneable,ptr_container/doc/tutorial.html#cloneability}.
Then, they can be conveniently manipulated by \refBoost{Boost.PointerContainer,ptr_container}s.


Future {#designrationalesfuture}
------
 
* Implement a more elaborate set of trajectories like `ExactTrajectory`, `HamiltonianTrajectory` *and* \link quantumtrajectory::MCWF_Trajectory MCWF_Trajectory\endlink.
Then, a maker function can decide on the basis of the dynamical characteristics of the quantum system which one to use. \todo Implement this



Implementation rationales {#implementationrationales}
=========================

General robust style + Effective C++ + Exceptional C++

* Be aware that \refStdCppConstruct{assert,cassert/assert} works only in debug mode (probably a macro switches it off in `NDEBUG`). It will not even evaluate its argument in `NDEBUG` mode!

* Note that function template partial specializations are not allowed, but the same effect can almost always be achieved by function template *overloads*.

  * Be aware of the possibility of “partial template-argument-inference” of *trailing* parameters, as exempified e.g. by

    * quantumoperator::Tridiagonal::doApply, where only the `REMAINING` argument is inferred

    * evolve, where only the `RANK` paramter is inferred

    * quantumdata::partialTrace, where the arguments `V` and `T` are not inferred, etc.

* When using \refBoostConstruct{BOOST_PP_ITERATE,preprocessor/doc/ref/iterate.html} to generate code, use *the same file* for reentering with the include guard 
`BOOST_PP_IS_ITERATING` *outside* the normal include guards (if a header file in is question, but a `.cc` file can be reentered like this as well).
Cf. the example at \refBoost{Boost.Preprocessor,preprocessor/doc/ref/iterate.html}. \see Composite.h, Sigma.cc, Tridiagonal.tcc, BlitzArraySliceIterator.tcc, etc.

  The advantage of this technique is that the generating pattern gets placed in the same file where it belongs.

* Note that declaring RANK template arguments as unsigned *does not work* because of template argument deduction problems. Unfortunately, Blitz++ declares them as int…

* Try to separate template declaration from definition as much as possible (without `export`ing, of course). One technique is explicit instatiation,
but this can be used only when the class is not very deeply templated (e.g. when there is only a single `RANK` argument).
For this, it is very important to factor out parameter-independent code from templates as much as possible. Explicit instatiation will be essential for a Python interface.

* When using containers of pointers to abstract base together with bind expessions always use `boost::bind` instead of `boost::lambda::bind`! See `c++essays/lambda_smartPtr` for a demonstration. Be *extremely* careful that Boost.Bind introduces the placeholders into the global namespace, so all kinds of clashes with `boost::lambda::bind` are easily created. In particular

      bind(foo,_1) 

  is correct, but

      bind(foo,bll::_1) 

  will never compile because bind is taken from boost as this is searched before `boost::lambda`.

      bll::bind(foo,bll::_1) 

  *must* be used instead!!!

  Rule is: use lambda when *really* needed, but bind in all other cases.

  From the above it also follows that `using namespace boost::lambda` should never really occur, use the namespace alias `bll` instead.

  Note that Boost.Bind takes bound arguments by value, and unbound arguments (`_n`) by non-const reference (same for Boost.Lambda). \todo Merge \refBoost{Boost.Bind,bind} and \refBoost{Boost.Lambda,lambda} into \refBoost{Boost.Phoenix,phoenix}

* Templates taking `int RANK` arguments may implement specializations, when needed, in a brute-force way, when the number of needed specializations depend linearly on RANK, just as in Blitz++. Examples:

  * BlitzArraySliceIterator.h -> blitzplusplus::basi::Transposer, blitzplusplus::basi::Indexer

  * ComplexArrayExtensions.h -> blitzplusplus::HermitianConjugateHelper

  * Tridiagonal.h -> quantumoperator::Tridiagonal::apply

  A notable exception is blitzplusplus::basi::Iterator itself because here the number of possibilities goes as `RANK` choose the size of `V`. Here obviously one has to resort to more sophisticated TMP techniques.

  Such brute force implementations should still rely on \refBoost{Boost.Preprocessor,preprocessor}.

* All offset (indexing) arguments should be of type `std::size_t` (`std::ptrdiff_t`) and all the extent arguments of type `std::size_t`. In template parameters use only `int`.

* Avoid hand-crafted loops whenever practical with the help of STL.Algorithm & \refBoost{Boost.Functional,functional} & \refBoost{Boost.Phoenix,phoenix}.

  * Use \refBoost{Boost.Range,range} whenever possible, together with range adaptors to avoid nested bind statements \todo Switch to range adaptors in existing codebase

* Consider design alternatives to dynamic_cast

* Plain pointers should practically never be used, but replaced with a smart pointer. Newed objects directly to `shared_ptr` in class factories. Use the `dynamic_pointer_cast` etc. templates.

* Avoid newed arrays, use `std::vector` instead together with `&v[0]`

* Effective C++ item 37!!!

* Do not inherit from STL containers.

* Prefer preincrements over postincrements, using the latter only when the original value is needed.

* If `xns` is a namespace, prefer the syntax `xns::` instead of re-opening the namespace when defining a function declared in that namespace. Cf. Stroustrup 8.2.9.3

* In return statements, be aware that the semantics of value return is the same as those of initialization.

* Use unnamed namespaces instead of `details` namespaces in `.cc` files.

* Classes not meant to be used from the outside should have protected constructors. Then, it's more convenient to use virtual bases which can have trivial protected constructors.
These (if not really meant to be used), can emit a runtime error, when called. Cf. `c++_essays/basics/virtualBaseSolution.cc`

* Identify uses of `boost::assign::repeat_fun()`. 

* Note the possibility of manual overload resolution by *statically* casting to the necessary function signature. This may enable the use of `boost::function` and `boost::bind` even in situations when they could not resolve the overload by themselves.

* In the factory idiom, encapsulate derived classes in makers (even maker functions). \see evolved::MakerGSL::_


C++11 features currently used {#cppelevenfeatures}
-----------------------------

* variadic templates – variadic parameter lists
* template aliases
* new function declaration syntax
* `auto` keyword => new `for` syntax
* rvalue references
* new initialization syntax + initialization lists
* delegating constructors

\note This sets the compiler requirement to g++ >= 4.7 ; clang++ >= 3.1

* To be adopted

  * inherited constructors 4.8, 3.3
  * initialization of static constants within the class
  * forward-declared enums 4.6, 3.1
  * Non-static data member initializers
  * lambda 4.5, 3.1
  * decltype 4.8, 2.9



Coding-style conventions {#codingstyleconventions}
========================

-# Reserve `struct` for TMP (traits classes, metafunctions, etc.) and stateless classes, otherwise always use `class`

-# Names representing types must be in mixed case starting with upper case

-# Variable names must be in mixed case starting with lower case

-# Private class variables must have underscore suffix (this makes them easier to spot than the underscore prefix.)

-# Named constants (including enumeration values, and template parameters of integral const type) must be all uppercase using underscore to separate words — Enumeration constants must be prefixed by a common type name

-# Names representing methods or functions should be verbs and written in mixed case starting with lower case — They are distinguished from variable names anyway by their special syntax

-# Names representing template metafunctions must be named as if they were functions, but in mixed case starting with upper case as if they were types

-# Names representing namespaces should be all lowercase

-# Names representing template types should be a single uppercase letter

-# Generic variables should have the same name as their type — Non-generic variables have a role: these variables can often be named by combining role and type

-# Variables with a large scope should have long names, variables with a small scope can have short names

-# The name of the object is implicit, and should be avoided in a method name

-# The terms `get`/`set` must be used where an attribute is accessed directly — This, together with items 2, 3, 4, 6, and 10, allows for very convenient naming, eg::

       class Foo
       {
       public:
         Foo(Bar bar) : bar_(bar) {}  

         bar  getBar(       ) const {return bar_;}
         void setBar(Bar bar)       {bar_=bar   ;}
       private:
         Bar bar_;
       };

-# Plural form should be used on names representing a collection of objects

-# The prefix `n` should be used for variables representing a number of objects

-# The suffix `No` should be used for variables representing an entity number

-# Iterator variables should be called `i`, `j`, `k`, etc.

-# The prefix `is` should be used for boolean variables and methods

-# Naming pointers specifically should be avoided

-# Negated boolean variable names must be avoided

-# Exception classes should be suffixed with `Exception`

-# Functions (methods returning something) should be named after what they return, procedures (void methods) after what they do

-# C++ header files should have the extension `.h`, source files can have the extension `.cpp` or `.cc`

-# Header files must contain an include guard. The guard macro name should be the all-uppercase filename including the path relative to the project root with `.` and `/` replaced by `_` and with the suffix  `_INCLUDED`. It must not start with an underscore.

-# Include statements should be sorted and grouped.  Sorted by their hierarchical position in the system with high level files included first. Leave an empty line between groups of include statements. Cf. below in code organization

-# Include statements must be located at the top of a file only

-# The parts of a class must be sorted `public`, `protected`, and `private`.  All sections must be identified explicitly. Not applicable sections should be left out.

-# Abbreviations and acronyms must not be uppercase when used as name

-# Complement names must be used for complement operations

-# C++ pointers and references should have their reference symbol next to the type rather than to the name

-# The nominal case should be put in the if-part and the exception in the else-part of an if statement

-# Use alignment wherever it enhances readability

-# Use // for all comments, including multi-line comments — then any larger passage can be commented out with `/* … */`

-# The function return type can be put in the left column immediately above the function name

-# Consider using a beautifier (bcpp?)



Code organization {#codeorganizationconventions}
=================

* Organize headers around a single concept: either a class (e.g. Evolved.h) or a coherent set of services (e.g. BlitzArrayExtensions.h).

* For headers that declare classes provide also a forward-declaring header with suffix `Fwd`. (Never manually forward declare anything, because this hardwires the forward declaration at several places, which is then hard to change if we want to migrate e.g. from a class to a template with default arguments.)

* Keep headers minimal but idempotent (apply include guards, etc.).

* Source files first include the corresponding header, next any headers developed for the project, then any other third party headers, and finally standard and system headers. This ensures that headers are genuinely self-contained, not accidentally relying on features they do not themselves include. Within each of the groupings use alphabetical ordering.

* Rely on forward declarations whenever possible. 

  Some cases when it is not enough: base classes, data members stored by value, `throw` expressions, etc. 

  `typedef`s and `enum`s may not be forward declared!

* A header should *always* include the corresponding `Fwd` header if it exists, because this helps to keep the two in synchron.

* Avoid the use of inline functions out of laziness.

* Template definitions must of course also be in header files with `.tcc` extension.

  The following is the preferred method to include implementation headers, because it reduces compilation-dependencies and enables IDEs to parse the headers correctly:

  * The `.tcc` header includes the corresponding header (also in the case of non-classes to keep the two in synchron) and further headers needed for the implementation.

  * The `.tcc` header must then be included at the point where the templates are actually instantiated, typically in `.cc` files.
In this case it is indicated for `.tcc` headers to include further `.tcc` headers for classes that are subsequently also instantiated.
This is for convenience of use. *However*: script developers should not be required to include `<some-implementation>.tcc` files.
Instead, high level include files should be provided which present a specific concept towards the user (e.g. Evolution.h, MCWF.h etc.).
Each element `Element` with template classes and/or a corresponding `ParsElement` class should be defined in an `Element_.h` and 
provide a high level `Element.h` which bundles `Element.tcc` (or `Element_.h` for non-template classes) and `ParsElement.h`.



Namespaces {#namespacescodeorganizationconventions}
----------

* The content of `utils` is treated as a library, so that here everything should be wrapped in namespaces

  * well-defined bundles in their own namespaces (this results in names like trajectory::Trajectory, but so be it)

  * all the rest in namespace `cpputils`

* Higher-level components of the framework can come in the global namespace, while their helpers come in namespaces named after their module (like Composite, but composite::Base).


Miscellaneous rationales, ideas, notes {#miscellaneousrationales}
======================================

* Instead of the current utils/Blitz2FLENS.h Blitz and FLENS could be better integrated by making a `GeneralStorageView` of a `blitz::Array` and passing it to `GEMatrix` as template parameter.

* Transposition of `blitz::Array` only permutes strides, does not touch storage, as we very well know. What is maybe not so clear is that *this remains the case after deep copying as well*, since also in this case the strides are just copied so that the result does not have obvious ordering.
To get an array with some conventional storage order, the way to go is to first create an array with the desired ordering and then assign to it from the original array.

  A more general solution would be to implement an Array class, which is *noncopyable* and inherits publicly from `blitz::Array`.
Then, in “copy-like” constructors and transposition, one must in each and every case specify whether deep or shallow semantics is meant.
In this case, since no real copy constructor is available, one should resort to another way to return these arrays from functions:
A straightforward solution is to return the underlying `blitz::Array`, and reconstruct the Array from this.
It should also provide something like a “cloning” member, which allows for dressing up a chunk of memory with the very same Array interface, with the correct storage layout, etc.
Clearly, this new Array class obsoletes blitzplusplus::TinyOfArrays. 

  The class could leave all work to the underlying `blitz::Array` storage, but solve this silly problem of `size_t`, `ptrdiff_t`, `int` conversions.
The difficulty is the construction, because `blitz::Array` provides an immense number of constructors, but with the constructor inheritance of C++11 it should be possible.

* Passing dcomp by reference is slightly faster if no optimizaton is used, but if the function is `inline`d, there is no difference.

* In structure::Interaction, it would be an interesting possibility to supply the type of the base class plugs as a compile-time vector (e.g. JaynesCummings could be derived from structure::Interaction`<mpl::vector<QbitBase,ModeBase> >`).
Then, Composite could do much more static sanity checks. The problem with this is that in this case these would appear in composite::_ as well and creep into the user interface, but this could be perhaps avoided somehow?
