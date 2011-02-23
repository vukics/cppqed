// -*- C++ -*-

// Implements the following Hamiltonian: 
// (-delta-i*kappa)*adagger*a+i*(eta*adagger-h.c.)
// with the usual Liouvillean

#ifndef   _MODE_FWD_H
#define   _MODE_FWD_H


namespace mode {

class Averaged;
// A basic, extensible Averaging class

class AveragedQuadratures;

template<typename Base=Averaged> class AveragedMonitorCutoff;

struct DoNotAverage {};

} // mode


class ModeBase;

template<typename A=mode::Averaged> class Mode;
template<typename A=mode::Averaged> class ModeSch;

template<typename A=mode::Averaged> class PumpedMode;
template<typename A=mode::Averaged> class PumpedModeSch;

// When not lossy, IP and UIP coincides

template<bool IS_FINITE_TEMP=false, typename A=mode::Averaged> class LossyMode;
template<bool IS_FINITE_TEMP=false, typename A=mode::Averaged> class LossyModeSch;
template<bool IS_FINITE_TEMP=false, typename A=mode::Averaged> class LossyModeUIP;

template<bool IS_FINITE_TEMP=false, typename A=mode::Averaged> class PumpedLossyMode;
template<bool IS_FINITE_TEMP=false, typename A=mode::Averaged> class PumpedLossyModeSch;
template<bool IS_FINITE_TEMP=false, typename A=mode::Averaged> class PumpedLossyModeUIP;



namespace mode {

struct Pars;
struct ParsPumped;
struct ParsLossy;
struct ParsPumpedLossy;

struct PrepError;

} // mode

#endif // _MODE_FWD_H
