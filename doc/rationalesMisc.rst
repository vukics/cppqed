
-----------------------------------------
Miscellaneous rationales, ideas, notes
-----------------------------------------

* Instead of the current ``utils/Blitz2FLENS.h`` Blitz and FLENS could be better integrated by making a ``GeneralStorageView`` of a ``blitz::Array`` and passing it to ``GEMatrix`` as template parameter.

* Transposition of ``blitz::Array`` only permutes strides, does not touch storage, as we very well know. What is maybe not so clear is that *this remains the case after deep copying as well*, since also in this case the strides are just copied so that the result does not have obvious ordering. The way to go is to first create an array with the desired ordering and then assign to it from the original array. 

  A more general solution would be to implement an Array class, which is *noncopyable* and inherits publicly from blitz::Array. Then in "copy-like" constructors and transposition, one must in each and every case specify whether deep or shallow semantics is meant. In this case, since no real copy constructor is available, one should resort to another way to return these arrays from functions: A straightforward solution is to return the underlying ``blitz::Array``, and reconstruct the Array from this. It should also provide something like a "cloning" member, which allows for dressing up a chunk of memory with the very same Array interface, with the correct storage layout, etc. Clearly, This new Array class obsolesces :class:`blitzplusplus::TinyOfArrays`. The difficulty is the construction, of course, as ``blitz::Array`` has an immense number of constructors.

* Passing ``dcomp`` by reference is slightly faster if no optimizaton is used, but if the function is ``inline``\ d, there is no difference.

* In :class:`~structure::Interaction`, it would be an interesting possibility to supply the type of the base class plugs as a compile-time vector (e.g. :class:`JaynesCummings` could be derived from :class:`~structure::Interaction`\ ``<mpl::vector<QbitBase,ModeBase> >``). Then, :class:`Composite` could do much more static sanity checks. The problem with this is that in this case these would appear in :class:`Act` as well and creep into the user interface, but this could be perhaps avoided somehow?

* The following code-organization scheme: header files are in pairs (``XFwd.h`` and ``X.h``, for templates also ``impl/X.tcc``), but the header files *do not* ``#include`` *anything*. Rather, files to be compiled ``#include`` the necessary headers *in the right order*. Fwd-declaration headers may sometimes be necessary to break loops.

  .. warning::

    *Do not do this, it is crazy!!!* Just imagine if a header starts to rely on a new ``#include``, how many times this has to be written into different files…
