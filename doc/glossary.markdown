Glossary {#glossary}
========

\tableofcontents

element {#glossaryelement}
=======

<em>Elements</em> are building blocks of composite systems. \see \ref genericelements "Supported general-purpose elements", [Main page of the <tt>elements</tt> module](\elementsMainPage)

They come in two flavours:

free {#glossaryelementfree}
----

Elementary free subsystem, where “free” means that the system has only one quantum number \see \ref genericelementsfrees "Supported general-purpose free elements", structure::Free

This kind of element can be used for simulating in itself.

interaction {#glossaryelementinteraction}
-----------

An interaction of arbitrary arity between such free systems. \see \ref genericelementsinteractions "Supported general-purpose interaction elements", structure::Interaction

composite {#glossarycomposite}
=========

System composed of several \ref glossaryelementfree elements linked by \ref glossaryelementinteraction elements

\see BinarySystem, Composite, \ref userguidemorecomplex

script {#glossaryscript}
======

High-level C++ program assembling the modules necessary for a given problem domain. E.g. specification of a particular physical system and of what to do with it – how to simulate it. It uses the framework as a library.

\see \ref userguide, directory `CPPQEDscripts`