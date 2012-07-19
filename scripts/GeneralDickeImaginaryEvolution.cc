#include "EvolutionHigh.h"
#include "ProjectingMCWF_Trajectory.h"

#include "GeneralDicke.h"

#include "BinarySystem.h"

#include "MatrixOfHamiltonian.h"

#include "MathExtensions.h"

#include<flens/flens.h>

#include<boost/lambda/bind.hpp>


using namespace std;
using namespace linalg;
using namespace blitz2flens;


typedef quantumtrajectory::ProjectingMCWF_Trajectory<2> MCWF;

typedef structure::Averaged<2> Averaged;

typedef quantumdata::StateVector<1> StateVector1;
typedef quantumdata::StateVector<2> StateVector2;

typedef TTD_DARRAY(1) DARRAY;

int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem
  try {

  ParameterTable p;

  ParsEvolution pe(p); // Driver Parameters

  mode::ParsLossy plm(p);

  spin::Pars ps(p);

  // generaldicke::Pars pgd(p);

  dcomp& u=p.add("u","General Dicke interaction u parameter",dcomp(1.,0));
  dcomp& y=p.add("y","General Dicke interaction y parameter",dcomp(1.,0));

  unsigned& numberOfProjectingStates=p.add("nops","Number of states to project the MCWF state vector on",6u);

  bool& noEigenStates=p.add("noEigenStates","Don't calculate eigenstates",false);

  pe.precision=6;

  // Parameter finalization
  update(p,argc,argv,"--");
  
  // ****** ****** ****** ****** ****** ******
  // ****** ****** ****** ****** ****** ******
  // ****** ****** ****** ****** ****** ******

  u/=ps.twoS+1; y/=sqrt(ps.twoS+1);

  plm.kappa=ps.gamma=0;

  const size_t dim=plm.cutoff*(ps.twoS+1);

  numberOfProjectingStates=noEigenStates ? 0 : min(size_t(numberOfProjectingStates),dim);
    
  MCWF::Basis eigenStates;
  DARRAY      eigenValues(dim);

  if (!noEigenStates) {
    LossyModeSch<> mode(plm);
    SpinSch        spin(ps );
    // Important that everything is in Sch picture here.

    GeneralDicke<> gd(mode,spin,u,y);
    BinarySystem sys(gd);

    static_cast<structure::QuantumSystem<2>&>(sys).displayParameters(cout);

    CMatrix 
      hamiltonian (calculateMatrix(sys)), 
      eigenVectors(hamiltonian.copy()  );

    assert(!mathutils::fcmp(1-max(abs(hamiltonian.copy()-blitzplusplus::hermitianConjugate(hamiltonian))),1,1e-12));
	
    {
      typedef HeMatrixMF<RowMajor>::type HeMatrix;
      HeMatrix a(hermitianMatrix(eigenVectors,RowMajorTag()));

      DenseVectorMF<double>::type v(blitz2flens::vector(eigenValues));
    
      ev(true,a,v);
    }

    for (int i=0; i<dim; i++) {
      eigenStates.push_back(new StateVector2(sys.getDimensions()));
      eigenStates[i].vectorView()=eigenVectors(i,blitz::Range::all());
    }

  }

  {
    // ****** ****** ****** ****** ****** ******
    // Mapping for imaginary time propagation:

    u*=-DCOMP_I;
    y*=-DCOMP_I;
    swap(plm.delta,plm.kappa); plm.kappa*=-1;
    swap(ps .omega,ps .gamma); ps .omega*=-1;

    pe.noise=false;

    mode::SmartPtr mode(mode::make(plm,QMP_IP));

    Spin spin(ps);

    StateVector1 psiMode(mode::init(plm)), psiSpin(spin.getDimension());

    psiSpin()(0)=1; psiSpin()(1)=1;

    StateVector2 psi(psiMode*psiSpin);

    psi.renorm();

    GeneralDicke<> gd(mode,spin,u,y);
    BinarySystem sys(gd);

    
    MCWF traj(psi,eigenStates,sys,pe);

    if (!noEigenStates) {
      using namespace boost;
      namespace bll=lambda;

      ostream& os=traj.getOstream();

      for (pair<MCWF::Basis::const_iterator,DARRAY::const_iterator> i(eigenStates.begin(),eigenValues.begin()); 
	   i.first!=eigenStates.end(); 
	   ++i.first, ++i.second) {
	os<<formdouble::high()(*i.second); Averaged::display(0,*i.first,os,pe.precision,&sys); os<<endl;
      }

    }

    evolve(traj,pe);

    // evolve(psi,BinarySystem(GeneralDicke(mode.get(),&spin,u,y)),pe,tmptools::Vector<0>());

  }

  } catch (const ParsNamedException& pne) {cerr<<"Pars named error: "<<pne.getName()<<endl;}



}

