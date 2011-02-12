// -*- C++ -*-
#ifndef _FREE_EXACT_H
#define _FREE_EXACT_H

#include "FreeExactFwd.h"

#include "Exact.h"


namespace structure {

class FreeExact : public Exact<1>
{
public:
  typedef TTD_CARRAY(1) Factors;

protected:
  explicit FreeExact(size_t dim) : dtDid_(0), factors_(int(dim)) {}

  Factors& getFactors() const {return factors_;}

private:
  void actWithU(double, StateVectorLow&) const;

  virtual void updateU(double) const = 0;

  // auxiliaries for U
  mutable double  dtDid_;
  mutable Factors factors_;

};

} // structure


#endif // _FREE_EXACT_H
