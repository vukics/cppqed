// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef   CPPQEDELEMENTS_FREES_QBIT_FWD_H_INCLUDED
#define   CPPQEDELEMENTS_FREES_QBIT_FWD_H_INCLUDED

class QbitBase;

class Qbit;
class QbitSch;

class PumpedQbit;
class PumpedQbitSch;

class LossyQbit;
class LossyQbitSch;
class LossyQbitUIP;

class PumpedLossyQbit;
class PumpedLossyQbitSch;
class PumpedLossyQbitUIP;


namespace qbit {

struct Pars;

struct ParsPumped;

template <typename BASE=Pars>
struct ParsLossy;

typedef ParsLossy<ParsPumped> ParsPumpedLossy;

template <typename BASE=ParsLossy<>>
struct ParsLossyPhaseNoise;

typedef ParsLossyPhaseNoise<ParsPumpedLossy> ParsPumpedLossyPhaseNoise;

// struct PrepError;

} // qbit


#endif // CPPQEDELEMENTS_FREES_QBIT_FWD_H_INCLUDED
