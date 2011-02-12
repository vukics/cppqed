#include "EvolutionHigh.h"
#include "Composite.h"

#include "Qbit.h"
#include "Mode.h"

#include "Interaction.h"
#include "TimeIndependentHamiltonian.h"



class MeiserToyModel : public structure::Interaction<2>, public structure::TimeIndependentHamiltonian<2>
{
public:
  MeiserToyModel(mode::SmartPtr cavity, mode::SmartPtr oscillator, double g);

  void addContribution(const StateVectorLow&, StateVectorLow&) const;

private:
  const linalg::CMatrix matrix_;

};


typedef structure::Interaction<3> I3;

int main(int argc, char* argv[])
{
  ParameterTable p;

  ParsEvolution pe(p);
  mode::ParsPumpedLossy pmC(p,"C");
  mode::ParsPumpedLossy pmO(p,"O");

  qbit::ParsPumpedLossy pq(p);

  double& g=p.add("g","Meiser's coupling",1.);

  update(p,argc,argv,"--");


  mode::SmartPtr cavity    (maker(pmC,QMP_IP));
  mode::SmartPtr oscillator(maker(pmO,QMP_IP));

  qbit::SmartPtr qbit(maker(pq,QMP_IP));

  MeiserToyModel mtm(cavity,oscillator,g);

  I3 dummy(I3::Frees(cavity.get(),oscillator.get(),qbit.get()));

  quantumdata::StateVector<3> psi(init(pmC)*init(pmO)*init(pq));

  evolve(psi,
	 makeComposite(	         
		       Act<0,1>  (mtm),
		       Act<0,1,2>(dummy)
				 ),
	 pe);

}


#include<boost/assign/list_of.hpp>
#include<boost/assign/list_inserter.hpp>

#include<boost/tuple/tuple_io.hpp>


using namespace boost::assign;

using namespace blitz;


MeiserToyModel::MeiserToyModel(mode::SmartPtr cavity, mode::SmartPtr oscillator, double g)
  : structure::Interaction<2>(Frees(cavity.get(),oscillator.get()),
			      tuple_list_of("g",g,cavity->getDimension())), 
    matrix_(shape(oscillator->getDimension(),oscillator->getDimension()))
{
  getParsStream()<<"# Meiser's toy model\n";
}


void MeiserToyModel::addContribution(const StateVectorLow& psi, StateVectorLow& dpsidt) const
{
  using tensor::i; using tensor::j;
  typedef mode::StateVectorLow SV1;

  for (int n=0; n<psi.extent(0); n++) {
    const SV1 psiSlice(psi(n,Range::all())); 
    SV1 dpsidtSlice(dpsidt(n,Range::all()));
    dpsidtSlice+=double(n)*sum(matrix_(i,j)*psiSlice(j),j);
  }

}
