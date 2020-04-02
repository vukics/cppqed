// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-

// Implements the following Hamiltonian: 
// (-delta-i*kappa)*adagger*a+i*(eta*adagger-h.c.)
// with the usual Liouvillean

#ifndef   CPPQEDELEMENTS_FREES_MODE_FWD_H_INCLUDED
#define   CPPQEDELEMENTS_FREES_MODE_FWD_H_INCLUDED


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

template<bool TEMPERATURE=false, typename AveragingType=mode::Averaged> class LossyMode;
template<bool TEMPERATURE=false, typename AveragingType=mode::Averaged> class LossyModeSch;
template<bool TEMPERATURE=false, typename AveragingType=mode::Averaged> class LossyModeUIP;

template<bool TEMPERATURE=false, typename AveragingType=mode::Averaged> class PumpedLossyMode;
template<bool TEMPERATURE=false, typename AveragingType=mode::Averaged> class PumpedLossyModeSch;
template<bool TEMPERATURE=false, typename AveragingType=mode::Averaged> class PumpedLossyModeUIP;



namespace mode {

struct Pars;
struct ParsPumped;
struct ParsLossy;
struct ParsPumpedLossy;

struct FockStatePreparationError_CheckYourCutoffAgainstDesiredFockState;

} // mode

#endif // CPPQEDELEMENTS_FREES_MODE_FWD_H_INCLUDED
