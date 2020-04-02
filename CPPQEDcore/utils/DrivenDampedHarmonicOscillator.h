// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-

/*
  
  Solution of the ODE

  diff(y(t),t,t)+2*gamma*diff(y(t),t)+y(t)-exp(I*omega*t)=0

  gamma==1 is critical damping

*/

#ifndef   CPPQEDCORE_UTILS_DRIVENDAMPEDHARMONICOSCILLATOR_H_INCLUDED
#define   CPPQEDCORE_UTILS_DRIVENDAMPEDHARMONICOSCILLATOR_H_INCLUDED

#include "DrivenDampedHarmonicOscillatorFwd.h"

#include "ComplexExtensions.h"

#include <flens/flens.h>

#include <boost/shared_ptr.hpp>


typedef boost::shared_ptr<DrivenDampedHarmonicOscillator> DDHO_Ptr;

const DDHO_Ptr makeDDHO(double gamma, double omega, dcomp ampTI, dcomp ampDerivTI, double tInit=0);


class DrivenDampedHarmonicOscillator
{
public:
  virtual const dcomp amp     (double t) const=0;
  virtual const dcomp ampDeriv(double t) const=0;

  virtual ~DrivenDampedHarmonicOscillator() {}

protected:
  typedef flens::GeMatrix<flens::FullStorage<dcomp,flens::ColMajor> > Matrix;
  typedef flens::DenseVector<flens::Array<dcomp> >                    Vector;

  static const Vector calculateA(Matrix, dcomp ampTI, dcomp ampDerivTI, dcomp cTInit, double omega);
  
  DrivenDampedHarmonicOscillator(double gamma, double omega, const Matrix&, dcomp ampTI, dcomp ampDerivTI, double tInit=0);
  // The initial condition is supplied at time tInit.

  const dcomp c(double t) const;

  const double gamma_, omega_;

  const Vector a_;

  static const Matrix makeMatrix(const dcomp&, const dcomp&, const dcomp&, const dcomp&);

};


#endif // CPPQEDCORE_UTILS_DRIVENDAMPEDHARMONICOSCILLATOR_H_INCLUDED
