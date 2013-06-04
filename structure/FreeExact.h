// -*- C++ -*-
#ifndef STRUCTURE_FREEEXACT_H_INCLUDED
#define STRUCTURE_FREEEXACT_H_INCLUDED

#include "FreeExactFwd.h"

#include "Exact.h"


namespace structure {

class FreeExact : public Exact<1>
{
public:
  typedef CArray<1> Factors;

protected:
  explicit FreeExact(size_t dim) : dtDid_(0), factors_(int(dim)) {}

  Factors& getFactors() const {return factors_;}

private:
  void actWithU_v(double, StateVectorLow&) const;

  virtual void updateU(double) const = 0;

  // auxiliaries for U
  mutable double  dtDid_;
  mutable Factors factors_;

};

} // structure


#endif // STRUCTURE_FREEEXACT_H_INCLUDED
