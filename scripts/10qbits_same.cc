#define BOOST_MPL_CFG_NO_PREPROCESSED_HEADERS // cf. http://www.boost.org/doc/libs/1_51_0/libs/mpl/doc/refmanual/limit-vector-size.html
#define BOOST_MPL_LIMIT_VECTOR_SIZE 30 // must be rounded to the nearest multiple of 10
#define FUSION_MAX_VECTOR_SIZE 23
#define BLITZ_ARRAY_LARGEST_RANK 23
#define QUANTUMOPERATOR_TRIDIAGONAL_MAX_RANK 2

#include "EvolutionComposite.h"
#include "JaynesCummings.h"
#include "QbitModeCorrelations.h"

using namespace std;

typedef quantumdata::StateVector<11> StateVector;


int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem
  try {

  ParameterTable p;

  ParsEvolution pe(p); // Driver Parameters
  qbit::ParsPumpedLossy pplqb(p);
  mode::ParsPumpedLossy pplm (p); 
  jaynescummings::Pars  pjc  (p);

  QM_Picture& qmp=p.add("picture","Quantum mechanical picture",QMP_IP);

  // Parameter finalization
  update(p,argc,argv,"--");
  
  // ****** ****** ****** ****** ****** ******

  if ((pe.evol==EM_MASTER /* || pe.evol==EM_convergence */) && qmp==QMP_IP) qmp=QMP_UIP;

  const qbit::Ptr qbit(qbit::make(pplqb,qmp));

  const mode::Ptr mode(mode::make(pplm ,qmp));

  const jaynescummings::Ptr jcCorr(jaynescummings::make<QbitModeCorrelations>(qbit,mode,pjc)), jc(jaynescummings::make(qbit,mode,pjc));
  
  StateVector psi(qbit::init(pplqb)*qbit::init(pplqb)*qbit::init(pplqb)*qbit::init(pplqb)*qbit::init(pplqb)*qbit::init(pplqb)*qbit::init(pplqb)*qbit::init(pplqb)*qbit::init(pplqb)*qbit::init(pplqb)*mode::init(pplm));
  psi.renorm();

  evolve<tmptools::Vector<0> >(psi,composite::make(Act<0,10>(jcCorr),
						   Act<1,10>(jc),
						   Act<2,10>(jc),
						   Act<3,10>(jc),
						   Act<4,10>(jc),
						   Act<5,10>(jc),
						   Act<6,10>(jc),
						   Act<7,10>(jc),
						   Act<8,10>(jc),
						   Act<9,10>(jc)),pe);
 
  } catch (const cpputils::TaggedException& te) {cerr<<"Caught exception with tag: "<<te.getTag()<<endl; exit(1);}


}
