C++QED: a framework for simulating open quantum dynamics {#mainpage}
========================================================

\par Abstract

C++QED is a highly flexible framework for simulating open quantum dynamics. It allows users to build arbitrarily complex interacting quantum systems from elementary free subsystems and interactions, and simulate their time evolution with a number of available time-evolution drivers.

The framework is very sensitive to performance, which can be increased by physics and computational ways. Among the physics ways to increase performance, the most notable is the maximal use of interaction picture, which may help to get rid of very separate timescales in the problem. Among the computational ones we can mention

- Maximal exploitation of special operator structures, i.e., sparse and tridiagonal matrices.

- The use of adaptive-stepsize methods for evolving ordinary differential equations.

- Judicious use of memory. Guideline: if it is necessary to copy something, ask first whether it is *really* necessary.


[<strong>The SourceForge.net project summary page</strong>](http://sourceforge.net/projects/cppqed/).

\tableofcontents

\note C++QED is a free open-source project, the only thing the developers ask is to please cite **both** our \ref mainpagepublications "papers" in publications where C++QED has any role, as this is the only way for us to obtain credit for this work.


Documentation {#mainpagedocumentation}
=============

[Click here for the Milestone 9 documentation](http://cppqed.sourceforge.net/oldSphinx/)

\par Installation of the framework
\ref installationguide

\par Background and fundamental design principles
\ref userguideintroduction

\par A quick start on writing and running actual physical simulations in C++
\ref userguide

\par The Python frontend to the framework
[CpypyQED](\cpypyqedMainPage)

\par The out-of-the-box physical elements
\ref genericelements

\par How to implement new physical systems
\ref structurebundleguide

\par Walkthrough for the API documentation
…

\par Some technical terms used throughout the documentation
\ref glossary


Download {#mainpagedownload}
========

[Click here for the Milestone 9 downloads](http://cppqed.sourceforge.net/oldSphinx/#download)

- Milestone 10 released packages

  - SourceForge file releases (coming soon!)
  - Debian packages (coming soon!)
  - AUR package (coming soon!)

- The Git version:

      git clone git://git.code.sf.net/p/cppqed/monolithic <directory name>

\note This online documentation is meant to reflect the development version, but delays compared to the actual state are possible. The released packages include the corresponding documentation.


Support {#mainpagesupport}
=======

- [the project mailing list](http://sourceforge.net/p/cppqed/mailman/cppqed-support/)

- [tracker system](http://sourceforge.net/p/cppqed/_list/tickets/)

We offer full support for the framework both in terms of writing and testing new elements and scripts on demand from the quantum optics community, and of advising with the use of the existing software in the framework or with the development of new software (elements, scripts, time-evolution drivers). In the first case we may require to become co-author in the publications stemming from the work.


Developers {#mainpagedevelopers}
==========

C++QED has been originally conceived and created by [András Vukics](http://optics.szfki.kfki.hu/Vukics/Vukics), who developed v1 between 2006–2008, and has developed and maintained v2 since 2008.

[Raimar Sandner](http://www.uibk.ac.at/th-physik/people/staffdb/660275.xml) joined the project in early 2012, and has since made substantial contributions. The [CpypyQED](\cpypyqedMainPage) Python frontend to C++QED is his work. He contributed the elaborate \ref cppqed_cmake "build system" and \ref testsuite "testsuite" of the framework. In addition, he made many important contributions in the design and implementation of the C++ codebase, such as the \ref userguideiotrajectorystate "trajectory-state input/output" feature. He is the maintainer of the Debian and AUR packages.


Publications {#mainpagepublications}
============

- András Vukics. *C++QEDv2: The multi-array concept and compile-time algorithms in the definition of composite quantum systems.* **Comp. Phys. Comm.**, 183(6):1381–1396, 2012. [(link)](http://www.sciencedirect.com/science/article/pii/S0010465512000562)

- A. Vukics and H. Ritsch. *C++QED: an object-oriented framework for wave-function simulations of cavity QED systems.* **Eur. Phys. J. D**, 44:585–599, 2007. [(link)](http://link.springer.com/article/10.1140%2Fepjd%2Fe2007-00210-x)

Links {#mainpagelinks}
=====

- [Homepage of A. Vukics at the Wigner Research Centre for Physics](http://optics.szfki.kfki.hu/Vukics/Vukics)

- [Boost C++ libraries](http://www.boost.org)

- [GNU Scientific library (GSL)](http://www.gnu.org/software/gsl)

* [Blitz++](http://blitz.sourceforge.net)

* [Flexible Library for Efficient Numerical Solutions (FLENS)](http://flens.sourceforge.net)


Acknowledgement {#mainpageacknowledgement}
===============

First of all, we acknowledge the developers of Boost for making C++ an even more powerful language than it originally was. Without the Boost libraries, the framework could not have taken form.

We would like to thank the developers of GSL, LAPACK, Blitz++, and FLENS for their effort, without which scientific computing in C++ in general, and the present framework in particular would not look as nice as it looks today.
