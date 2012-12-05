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
  qbit::ParsPumpedLossy pplqb2(p,"2");
  qbit::ParsPumpedLossy pplqb3(p,"3"); 
  qbit::ParsPumpedLossy pplqb4(p,"4");
  qbit::ParsPumpedLossy pplqb5(p,"5");
  qbit::ParsPumpedLossy pplqb6(p,"6");
  qbit::ParsPumpedLossy pplqb7(p,"7");
  qbit::ParsPumpedLossy pplqb8(p,"8"); 
  qbit::ParsPumpedLossy pplqb9(p,"9");
  qbit::ParsPumpedLossy pplqb10(p,"10"); 

  mode::ParsPumpedLossy pplm (p); 
  jaynescummings::Pars  pjc  (p);
  jaynescummings::Pars  pjc2  (p,"2");
  jaynescummings::Pars  pjc3  (p,"3");
  jaynescummings::Pars  pjc4  (p,"4"); 
  jaynescummings::Pars  pjc5  (p,"5");
  jaynescummings::Pars  pjc6  (p,"6");
  jaynescummings::Pars  pjc7  (p,"7");
  jaynescummings::Pars  pjc8  (p,"8");
  jaynescummings::Pars  pjc9  (p,"9");
  jaynescummings::Pars  pjc10  (p,"10");

  QM_Picture& qmp=p.add("picture","Quantum mechanical picture",QMP_IP);

  // Parameter finalization
  update(p,argc,argv,"--");
  
  // ****** ****** ****** ****** ****** ******

  if ((pe.evol==EM_MASTER /* || pe.evol==EM_convergence */) && qmp==QMP_IP) qmp=QMP_UIP;

  qbit::Ptr qbit(qbit::make(pplqb,qmp));
  qbit::Ptr qbit2(qbit::make(pplqb2,qmp));
  qbit::Ptr qbit3(qbit::make(pplqb3,qmp));
  qbit::Ptr qbit4(qbit::make(pplqb4,qmp));
  qbit::Ptr qbit5(qbit::make(pplqb5,qmp));
  qbit::Ptr qbit6(qbit::make(pplqb6,qmp));
  qbit::Ptr qbit7(qbit::make(pplqb7,qmp));
  qbit::Ptr qbit8(qbit::make(pplqb8,qmp));
  qbit::Ptr qbit9(qbit::make(pplqb9,qmp));
  qbit::Ptr qbit10(qbit::make(pplqb10,qmp));

  mode::Ptr mode(mode::make(pplm ,qmp));

  StateVector psi(qbit::init(pplqb)*qbit::init(pplqb2)*qbit::init(pplqb3)*qbit::init(pplqb4)*qbit::init(pplqb5)*qbit::init(pplqb6)*qbit::init(pplqb7)*qbit::init(pplqb8)*qbit::init(pplqb9)*qbit::init(pplqb10)*mode::init(pplm));
  psi.renorm();

  evolve<tmptools::Vector<0> >(psi,composite::make(Act<0,10>(jaynescummings::make<QbitModeCorrelations>(qbit,mode,pjc)),
						   Act<1,10>(jaynescummings::make(qbit2,mode,pjc2)),
						   Act<2,10>(jaynescummings::make(qbit3,mode,pjc3)),
						   Act<3,10>(jaynescummings::make(qbit4,mode,pjc4)),
						   Act<4,10>(jaynescummings::make(qbit5,mode,pjc5)),
						   Act<5,10>(jaynescummings::make(qbit6,mode,pjc6)),
						   Act<6,10>(jaynescummings::make(qbit7,mode,pjc7)),
						   Act<7,10>(jaynescummings::make(qbit8,mode,pjc8)),
						   Act<8,10>(jaynescummings::make(qbit9,mode,pjc9)),
						   Act<9,10>(jaynescummings::make(qbit10,mode,pjc10))),pe);
 
  } catch (const cpputils::TaggedException& te) {cerr<<"Caught exception with tag: "<<te.getTag()<<endl; exit(1);}


}
