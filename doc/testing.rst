.. _testing:


====================================
Principles and examples of testing
====================================


A simplified testsuite::

  compile -j3 PTLA_Evolved PTLA_C++QED PumpedLossyQbit PTLA_EvolvedHL PumpedLossyMode_Evolved PumpedLossyMode_C++QED QbitMode_C++QED QbitMode_Evolved QbitMode_Matrix release

requiring only the elements::

  interactions/JaynesCummings.cc frees/PumpedTwoLevelAtom.cc frees/TimeIndependentMatrixHamiltonian.cc frees/QM_Picture.cc [ glob frees/*Mode*.cc frees/*Qbit*.cc composites/*.cc ] 
