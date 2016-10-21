// Copyright András Vukics 2006–2016. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-
#ifndef   CPPQEDCORE_UTILS_TRAJECTORYFWD_H_INCLUDED
#define   CPPQEDCORE_UTILS_TRAJECTORYFWD_H_INCLUDED


namespace trajectory {

struct ParsRun;

struct ParsEvolved;

class StoppingCriterionReachedException;
class TrajectoryFileOpeningException;
class StateFileOpeningException;

struct SerializationMetadata;

class Trajectory;

template<typename>
class AdaptiveIO;

template<typename>
class Adaptive;

} // trajectory

#endif // CPPQEDCORE_UTILS_TRAJECTORYFWD_H_INCLUDED
