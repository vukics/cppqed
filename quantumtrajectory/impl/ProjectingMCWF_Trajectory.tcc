// -*- C++ -*-
#ifndef   PROJECTING_MCWF_TRAJECTORY_IMPL_INCLUDED
#define   PROJECTING_MCWF_TRAJECTORY_IMPL_INCLUDED


#include "ProjectingMCWF_Trajectory.h"

#include "Blitz2FLENS.h"


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

  // now we have the dd components of the metric, and we have to
  // invert it to get the uu components:
  if (dim) {
    using namespace blitz2flens;
    typedef GeMatrixMF<dcomp,RowMajor>::type GeMatrix;

    GeMatrix a(matrix(res,RowMajorTag()));
    DenseVector<Array<int> > pivots(dim);
    
    trf(a,pivots);
    tri(a,pivots);
  }

  return res;

}


template<int RANK>
ProjectingMCWF_Trajectory<RANK>::ProjectingMCWF_Trajectory(
							   StateVector& psi,
							   const Basis& basis,
							   const QuantumSystem& sys,
							   const ParsMCWF_Trajectory& p,
							   const StateVectorLow& scaleAbs
							   )
  : TrajectoryBase(p), Base(psi,sys,p,scaleAbs),
    basis_(basis), metricTensor_uu_(help())
{
}


template<int RANK>
void
ProjectingMCWF_Trajectory<RANK>::displayEvenMore(int precision) const
{
  using namespace formdouble;

  const StateVector& psi=getPsi();
  const FormDouble fd(precision);
 
  if (int dim=basis_.size()) {
    getOstream()<<"\t";
    double sumR=0;
    for (int i=0; i<dim; i++) {
      double temp=0;
      for (int j=0; j<dim; j++) temp+=real(braket(psi,basis_[j])*braket(basis_[i],psi)*metricTensor_uu_(j,i));
      getOstream()<<fd(mathutils::sqrAbs(braket(basis_[i],psi)));
      sumR+=temp;
    }
    getOstream()<<'\t'<<fd(sumR);
  }
}


template<int RANK>
size_t
ProjectingMCWF_Trajectory<RANK>::displayMoreKey() const
{
  size_t res=Base::displayMoreKey();
  getOstream()<<"# ProjectingMCWF_Trajectory "<<res<<'-'<<res+basis_.size()-1<<". Overlap of Monte-Carlo state vector with basis vectors "<<res+basis_.size()<<". Sum of overlaps"<<std::endl;
  return res;
}


} // quantumtrajectory


#endif // PROJECTING_MCWF_TRAJECTORY_IMPL_INCLUDED
