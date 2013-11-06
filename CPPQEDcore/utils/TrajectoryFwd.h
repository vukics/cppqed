// -*- C++ -*-
#ifndef   UTILS_TRAJECTORYFWD_H_INCLUDED
#define   UTILS_TRAJECTORYFWD_H_INCLUDED


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

#endif // UTILS_TRAJECTORYFWD_H_INCLUDED
