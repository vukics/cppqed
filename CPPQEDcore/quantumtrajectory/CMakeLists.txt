# Copyright Raimar Sandner 2012–2020.
# Copyright András Vukics 2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

add_library(CPPQEDcore_quantumtrajectory SHARED Evolution.cc MCWF_TrajectoryLogger.cc ParsMCWF_Trajectory.cc)
add_library(CPPQED::CPPQEDcore_quantumtrajectory ALIAS CPPQEDcore_quantumtrajectory)

target_link_libraries(CPPQEDcore_quantumtrajectory PUBLIC CPPQED::CPPQEDcore_structure)
