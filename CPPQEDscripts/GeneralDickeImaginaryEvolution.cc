// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "EvolutionBinary.h"

#include "Mode.h"

#include "GeneralDicke.h"

#include "MatrixOfHamiltonian.h"

#include "MathExtensions.h"

#include "ProjectingMCWF_Trajectory.tcc"

#include<flens/flens.h>


using namespace std;
using namespace linalg;
using namespace blitz2flens;


typedef quantumtrajectory::ProjectingMCWF_Trajectory<2> MCWF;

typedef structure::Averaged<2> Averaged;

typedef quantumdata::StateVector<1> StateVector1;
typedef quantumdata::StateVector<2> StateVector2;

typedef DArray<1> DARRAY;

int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem

  ParameterTable p;

  evolution::Pars<> pe(p); // Driver Parameters

  mode::ParsDissipative plm(p);

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
    DissipativeModeSch<> mode(plm);
    SpinSch        spin(ps );
    // Important that everything is in Sch picture here.

    GeneralDicke<> gd(mode,spin,u,y);
    
    binary::Ptr sys(binary::make(gd));

    sys->streamParameters(cout);

    CMatrix 
      hamiltonian (calculateMatrix(*sys)), 
      eigenVectors(hamiltonian.copy()  );

    assert(!mathutils::fcmp(1-max(abs(hamiltonian.copy()-blitzplusplus::hermitianConjugate(hamiltonian))),1,1e-12));

    {
      HeMatrixOf<RowMajor> a(hermitianMatrix<RowMajor>(eigenVectors));

      DenseVectorOf<double> v(blitz2flens::vector(eigenValues));
    
      // Hermitian eigenvalues are apparently not implemented in FLENS-LAPACK at the moment – should be replaced by generic complex eigenvalues
      // flens::lapack::ev(true,a,v);
    }

    for (int i=0; i<dim; i++) {
      eigenStates.push_back(new StateVector2(sys->getDimensions()));
      eigenStates[i].vectorView()=eigenVectors(i,blitz::Range::all());
    }

  }

  {
    // ****** ****** ****** ****** ****** ******
    // Mapping for imaginary time propagation:

    u*=-1i;
    y*=-1i;
    swap(plm.delta,plm.kappa); plm.kappa*=-1;
    swap(ps .omega,ps .gamma); ps .omega*=-1;

    pe.noise=false;

    mode::Ptr mode(mode::make(plm,QMP_IP));

    Spin spin(ps);

    StateVector1 psiMode(mode::init(plm)), psiSpin(spin.getDimension());

    psiSpin(0)=1; psiSpin(1)=1;

    StateVector2 psi(psiMode*psiSpin);

    psi.renorm();

    GeneralDicke<> gd(mode,spin,u,y);
    
    binary::Ptr sys(binary::make(gd));
    
    MCWF traj(psi,eigenStates,sys,pe);

    if (!noEigenStates) {
      using namespace boost;

      ostream& os=cout;

      for (pair<MCWF::Basis::const_iterator,DARRAY::const_iterator> i(eigenStates.begin(),eigenValues.begin()); 
           i.first!=eigenStates.end();
           ++i.first, ++i.second) {
        os<<formdouble::high()(*i.second); structure::stream(structure::qsa<2>(sys),0,*i.first,os,pe.precision); os<<endl;
      }

    }

    run(traj,pe);

    // evolve<0>(psi,BinarySystem(GeneralDicke(mode.get(),&spin,u,y)),pe);

  }

}

