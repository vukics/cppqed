
-------------------------
Implementation rationales
-------------------------

General robust style + Effective C++ + Exceptional C++

* Be aware that ``assert`` works only in debug mode (probably a macro switches it off in ``NDEBUG``). It will not even evaluate its argument in ``NDEBUG`` mode!

* Note that function template partial specializations are not allowed, but the same effect can almost always be achieved by function template *overloads* which are of course allowed!

* Note that declaring RANK template arguments as unsigned *does not work* because of template argument deduction problems. Blitz declares them as int... No interoperation...

* Try to separate template declaration from definition as much as possible (without ``export``\ ing, of course). One technique is explicit instatiation, but this can be used only when the class is not very deeply templated (e.g. when there is only a single RANK parameter). For this, it is very important to factor out parameter-independent code from templates as much as possible. Explicit instatiation will be essential for a Python interface.

* When using containers of pointers to abstract base together with bind expessions always use ``boost::bind`` instead of ``boost::lambda::bind``! See ``c++essays/lambda_smartPtr`` for a demonstration. Be *extremely* careful that ``boost::bind`` introduces the placeholders into the global namespace, so all kinds of clashes with ``boost::lambda::bind`` are easily created. In particular :: 

	bind(foo,_1) 

  is correct, but ::

	bind(foo,bll::_1) 

  will never compile because bind is taken from boost as this is searched before boost::lambda. ::

	bll::bind(foo,bll::_1) 

  *must* be used instead!!!

  Rule is: use lambda when *really* needed, but bind in all other cases.

  From the above it also follows that ``using namespace boost::lambda`` should never really be used, use the namespace alias ``bll`` instead.

  Note that ``boost::bind`` takes bound arguments by value, and unbound arguments (_n) by non-const reference (same for boost::lambda).

  .. note::

    Such problems should disappear when switching to Boost.Phoenix, which replaces both Boost.Bind and Boost.Lambda.

* With templates taking ``int RANK`` arguments may implement specializations, when needed, in a brute-force way, when the number of needed specializations depend linearly on RANK, just as in Blitz++.

  Examples: 

  ``BlitzArraySliceIterator.h`` -> :class:`~blitzplusplus::basi::Transposer`, :class:`~blitzplusplus::basi::Indexer`

  ``ComplexArrayExtensions.h`` -> :func:`blitzplusplus::HermitianConjugateHelper`

  ``Tridiagonal.h`` -> :func:`~quantumoperator::Tridiagonal::apply`

  A notable exception is :class:`blitzplusplus::basi::Iterator`\ ``<RANK,V>`` itself because here the number of possibilities goes as ``RANK`` choose the size of ``V``. Here obviously one has to resort to more sophisticated TMP techniques.

  Such brute force implementations should still rely on Boost.Preprocessor.

* All the offset (indexing) args of type ``std::size_t`` (``std::ptrdiff_t``) and all the extent args of type ``std::size_t``. In template parameters use only int!

* Avoid loops whenever practical with the help of STL.Algorithm & Boost.Functional & Boost.Bind & Boost.Lambda.

* consider alternatives to dynamic_cast---static_cast, but also design alternatives

* newed objects directly to ``shared_ptr`` (factory functions!)

* avoid newed arrays, use ``std::vector`` instead together with ``&v[0]``

* Effective C++ item 37!!!

* No inheritance from STL containers.

* Prefer preincrements over postincrements, using the latter only when the original value is needed.

* If ``xns`` is a namespace, prefer the syntax ``xns::`` instead of re-opening the namespace when defining a function declared in that namespace. Cf BS. 8.2.9.3

* In return statements, be aware that the semantics of value return is the same as those of initialization.

* Use unnamed namespaces instead of detail namespaces in ``.cc`` files.

* Classes not meant to be used from the outside should have protected ctors. Then, it's more convenient to use protected bases which can have trivial protected ctors. These (if not really meant to be used), can emit a runtime error, when called. Cf. ``c++_essays/basics/virtualBaseSolution.cc``

* Identify uses of ``boost::assign::repeat_fun()``. 

* Note the possibility of manual overload resolution by *statically* casting to the necessary function signature. This may enable the use of ``boost::function`` and ``boost::bind`` even in situations when they could not resolve the overload by themselves.

C++0x
^^^^^^^^^

* Template typedefs are at the moment "defined" in the form of this ``...MF`` craziness, in the case of stateless classes by inheritance, and often by macros. Problem will be obsoleted in C++0x.

* Functions cannot be defined in function scope, though this would be often very convenient for defining helpers. However, the effect can be emulated because classes *can* be defined in function scope. These must not have friend functions, but they can have static functions. Problem will be obsoleted in C++0x.

Further useful features:

* Variadic templates

* Rvalue references --- identify uses

* ``auto`` + ``decltype`` keyword --- several uses

* Defaulted and deleted functions
