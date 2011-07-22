.. _codeOrganization:

*****************
Code organization
*****************


* Organize headers around a single concept: either a class (e.g. ``Evolved.h``) or a coherent set of services (eg ``BlitzArrayExtensions.h``)

* For headers that declare classes provide also a forward-declaring header with suffix ``Fwd``. (Never manually forward declare anything, because this hardwires the forward declaration at several places, which is then hard to change if we want to migrate eg from a class to a template with default arguments.)

* Keep headers minimal but idempotent (apply include guards, etc.)

* Source files first include the corresponding header, next any headers developed for the project, then any other third party headers, and finally standard and system headers. This ensures that headers are genuinely self-contained, not accidentally relying on features they do not themselves include. Within each of the groupings use alphabetical ordering.

* Rely on forward declarations whenever possible. 

  Some cases when it is not enough: base classes, data members that are stored by value, throw expressions, etc. 

  Typedefs and enums may not be forward declared!

* A header should *always* include the corresponding ``Fwd`` header if it exists, because this helps to keep the two in synchron.

* Avoid the use of inline functions out of laziness.

* Template definitions must of course also be in header files with suffix ``Impl`` or a subdirectroy ``impl``. They should have ``.tcc`` extension.

  There are two alternatives at the moment:

  1. The ``Impl`` header includes the corresponding header (also in the case of non-classes to keep the two in synch) and further headers needed for the implementation. The ``Impl`` header must then be included at the point where the templates are actually instantiated, typically in ``.cpp`` files. In this case it is probably indicated for ``Impl`` headers to include further ``Impl`` headers for classes that are subsequently also instantiated. This is for convenience of use.

  2. The header includes the ``Impl`` header at the end of the file. This is even more convenient for users, but not really a *minimal* approach because definitions get included everywhere where declarations would have sufficed.


