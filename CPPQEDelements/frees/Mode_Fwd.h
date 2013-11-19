// -*- C++ -*-

// Implements the following Hamiltonian: 
// (-delta-i*kappa)*adagger*a+i*(eta*adagger-h.c.)
// with the usual Liouvillean

#ifndef   FREES_MODE_FWD_H_INCLUDED
#define   FREES_MODE_FWD_H_INCLUDED


namespace mode {

class Averaged;
// A basic, extensible Averaging class

class AveragedQuadratures;

template<typename Base=Averaged> class AveragedMonitorCutoff;

struct DoNotAverage {};

} // mode


class ModeBase;

template<typename AveragingType=mode::Averaged> class Mode;
template<typename AveragingType=mode::Averaged> class ModeSch;

template<typename AveragingType=mode::Averaged> class PumpedMode;
template<typename AveragingType=mode::Averaged> class PumpedModeSch;

// When not lossy, IP and UIP coincides

template<bool IS_FINITE_TEMP=false, typename AveragingType=mode::Averaged> class LossyMode;
template<bool IS_FINITE_TEMP=false, typename AveragingType=mode::Averaged> class LossyModeSch;
template<bool IS_FINITE_TEMP=false, typename AveragingType=mode::Averaged> class LossyModeUIP;

template<bool IS_FINITE_TEMP=false, typename AveragingType=mode::Averaged> class PumpedLossyMode;
template<bool IS_FINITE_TEMP=false, typename AveragingType=mode::Averaged> class PumpedLossyModeSch;
template<bool IS_FINITE_TEMP=false, typename AveragingType=mode::Averaged> class PumpedLossyModeUIP;



namespace mode {

struct Pars;
struct ParsPumped;
struct ParsLossy;
struct ParsPumpedLossy;

struct FockStatePreparationError_CheckYourCutoffAgainstDesiredFockState;

} // mode

#endif // FREES_MODE_FWD_H_INCLUDED
