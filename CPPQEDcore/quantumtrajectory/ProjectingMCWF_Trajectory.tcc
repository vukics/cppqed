// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-
#ifndef   CPPQEDCORE_QUANTUMTRAJECTORY_PROJECTINGMCWF_TRAJECTORY_TCC_INCLUDED
#define   CPPQEDCORE_QUANTUMTRAJECTORY_PROJECTINGMCWF_TRAJECTORY_TCC_INCLUDED


#include "ProjectingMCWF_Trajectory.h"

#include "MCWF_Trajectory.tcc"

#include "Blitz2FLENS.tcc"


namespace quantumtrajectory {


template<int RANK>
const linalg::CMatrix
ProjectingMCWF_Trajectory<RANK>::help() const
{
  int dim=basis_.size();

  linalg::CMatrix res(blitz::shape(dim,dim));

  for (int i=0; i<dim; i++) {
    res(i,i)=basis_[i].norm();
    for (int j=0; j<dim; j++)
      res(j,i)=conj(res(i,j)=braket(basis_[i],basis_[j]));
  }

  // now we have the dd components of the metric, and we have to invert it to get the uu components:
  if (dim) {
    using namespace blitz2flens;

    typedef GeMatrixOf<dcomp,RowMajor> GeMatrix;

    GeMatrix a(matrix<RowMajor>(res));
    flens::DenseVector<flens::Array<int> > pivots(dim);
    
    flens::lapack::trf(a,pivots);
    flens::lapack::tri(a,pivots);
  }

  return res;

}


template<int RANK>
std::ostream&
ProjectingMCWF_Trajectory<RANK>::display_v(std::ostream& os, int precision) const
{
  using namespace formdouble;

  const StateVector& psi=this->getPsi();
  const FormDouble fd(precision);
  
  Base::display_v(os,precision);
  
  if (int dim=basis_.size()) {
    os<<"\t";
    double sumR=0;
    for (int i=0; i<dim; i++) {
      double temp=0;
      for (int j=0; j<dim; j++) temp+=real(braket(psi,basis_[j])*braket(basis_[i],psi)*metricTensor_uu_(j,i));
      os<<fd(mathutils::sqrAbs(braket(basis_[i],psi)));
      sumR+=temp;
    }
    os<<'\t'<<fd(sumR);
  }
  return os;
}


template<int RANK>
std::ostream&
ProjectingMCWF_Trajectory<RANK>::displayKey_v(std::ostream& os, size_t& i) const
{
  return Base::displayKey_v(os,i)<<"# ProjectingMCWF_Trajectory\n# "
                             <<i<<'-'<<i+basis_.size()-1<<". Overlap of Monte-Carlo state vector with basis vectors\n# "
                             <<i+basis_.size()<<". Sum of overlaps\n";
}


} // quantumtrajectory


#endif // CPPQEDCORE_QUANTUMTRAJECTORY_PROJECTINGMCWF_TRAJECTORY_TCC_INCLUDED
