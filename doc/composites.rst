.. _composites:

==============
Composites
==============

-------------------
``BinarySystem``
-------------------

.. class:: BinarySystem

-----------------
``Composite``
-----------------

.. class:: Composite

   Assume that we have a system composed of modes and particles layed out in the following way: (Free No. 0) mode (1) particle (2) mode (3) particle (4) mode (5) particle (6) mode (7) particle. Assume that we need the interaction :class:`ParticleTwoModes2D` acting between (4) (0) (7) (1). Furthermore assume that ::

     const int RANK=8;
     structure::Types<8>::StateVectorLow psi, dpsidt;
     
   represent state vectors of the *complete* system. Then the action of the Hamiltonian of this interaction in this complete Hilbert space can be calculated as ::
   
     for_each(fullRange(psi,Vector<4,0,7,1>()),basi::begin(dpsidt,Vector<4,0,7,1>()),PTM2D_Hamiltonian);
     
   ...
   
   It is the task of this class to keep track of the elements of the system together with their corresponding slices, and implement the operations of the full system in terms of such loops...