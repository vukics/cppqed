// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

/*
  
  Solution of the ODE

  diff(y(t),t,t)+2*gamma*diff(y(t),t)+y(t)-exp(I*omega*t)=0

  gamma==1 is critical damping

*/

#ifndef   CPPQEDCORE_UTILS_DRIVENDAMPEDHARMONICOSCILLATOR_H_INCLUDED
#define   CPPQEDCORE_UTILS_DRIVENDAMPEDHARMONICOSCILLATOR_H_INCLUDED

#ifdef    EIGEN3_FOUND

#include "ComplexExtensions.h"

#include <Eigen/Dense>

#include <memory>


namespace cppqedutils::ddho {

typedef Eigen::Matrix2cd Matrix;
typedef Eigen::Vector2cd Vector;

class _
{
public:
  virtual dcomp amp     (double t) const=0;
  virtual dcomp ampDeriv(double t) const=0;

  virtual ~_() {}

protected:
  _(double gamma, double omega, Matrix, dcomp ampTI, dcomp ampDerivTI, double tInit=0); // The initial condition is supplied at time tInit.

  dcomp c(double t) const; // The particular solution of the inhomogeneous equation

  const double gamma_, omega_;

  const Vector a_;

};


typedef std::shared_ptr<_> Ptr;

Ptr make(double gamma, double omega, dcomp ampTI, dcomp ampDerivTI, double tInit=0);

} // cppqedutils::ddho

#endif // EIGEN3_FOUND

#endif // CPPQEDCORE_UTILS_DRIVENDAMPEDHARMONICOSCILLATOR_H_INCLUDED
