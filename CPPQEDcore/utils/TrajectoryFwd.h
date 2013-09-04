// -*- C++ -*-
#ifndef   UTILS_INCLUDE_TRAJECTORYFWD_H_INCLUDED
#define   UTILS_INCLUDE_TRAJECTORYFWD_H_INCLUDED


namespace trajectory {

struct ParsRun;

struct ParsEvolved;

class StoppingCriterionReachedException;
class TrajectoryFileOpeningException;
class StateFileOpeningException;

struct SerializationMetadata;

class Trajectory;

template<typename>
class Adaptive;

} // trajectory

#endif // UTILS_INCLUDE_TRAJECTORYFWD_H_INCLUDED
