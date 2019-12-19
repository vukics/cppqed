//#define BOOST_MPL_CFG_NO_PREPROCESSED_HEADERS // cf. http://www.boost.org/doc/libs/1_51_0/libs/mpl/doc/refmanual/limit-vector-size.html
//#define BOOST_MPL_LIMIT_VECTOR_SIZE 30 // must be rounded to the nearest multiple of 10
//#define FUSION_MAX_VECTOR_SIZE 21 // for 10qbits
#define BLITZ_ARRAY_LARGEST_RANK 14
#define QUANTUMOPERATOR_TRIDIAGONAL_MAX_RANK 5 // 14 as maximal rank would be an overkill here

#include "EvolutionComposite.h"

#include "AveragingUtils.tcc"
#include "JaynesCummings.h"
#include "QbitModeCorrelations.h"

using namespace std;

typedef quantumdata::StateVector<7> StateVector;


int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem

  ParameterTable p;

  evolution::Pars<> pe(p); // Driver Parameters
  qbit::ParsPumpedLossy pplqb(p);
  mode::ParsPumpedLossy pplm (p); 
  jaynescummings::Pars  pjc  (p);

  // Parameter finalization
  auto qmp=updateWithPicture(p,argc,argv,"--");
  
  // ****** ****** ****** ****** ****** ******

  const qbit::Ptr qbit(qbit::make(pplqb,qmp));

  const mode::Ptr mode(mode::make<mode::AveragedQuadratures>(pplm ,qmp));

  const jaynescummings::Ptr jc(jaynescummings::make(qbit,mode,pjc)),
                            jcRDO(jaynescummings::make<ReducedDensityOperatorNegativity<2,tmptools::Vector<0> > >(qbit,mode,pjc,"JaynesCummings1-6",jc->getDimensions())),
                            jcCorr(jaynescummings::make<QbitModeCorrelations>(qbit,mode,pjc));
                            
  
  StateVector psi(qbit::init(pplqb)*qbit::init(pplqb)*qbit::init(pplqb)*qbit::init(pplqb)*qbit::init(pplqb)*qbit::init(pplqb)*mode::init(pplm));
  psi.renorm();

  evolve(psi,composite::make(Act<0,6>(jcCorr),
                             Act<1,6>(jcRDO),
                             Act<2,6>(jc),
                             Act<3,6>(jc),
                             Act<4,6>(jc),
                             Act<5,6>(jc)),
         pe);
 



}
