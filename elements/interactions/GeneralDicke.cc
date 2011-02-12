#include "GeneralDicke.h"

#include<boost/assign/list_of.hpp>


using namespace boost::assign;

using namespace structure;
using namespace spin;
using namespace mode;


namespace {

size_t spinN(const SpinBase* spin)
{
  return spin->getDimension()-1;
}

}


generaldicke::Base::Base(const ModeBase* mode, const SpinBase* spin, dcomp u, dcomp y)
  : Interaction<2>(Frees(mode,spin),
		   RealFreqs(),
		   tuple_list_of("u",u,mode->getDimension()*spin->getDimension())("y",y,
										  sqrt(mode->getDimension()*spin->getDimension()))
		   ),
    TridiagonalHamiltonian<2,true>((u*nop(mode)*(sz(spinN(spin))+spinN(spin)/2.*quantumoperator::identity(spin->getDimension()))
				   +
				    y*(aop(mode).dagger()+aop(mode))*sx(spinN(spin)))/DCOMP_I,
				   freqs(mode)*freqs(spin))
{
  getParsStream()<<"# General Dicke interaction\n";
  // getParsStream()<<sx(spin->getDimension()-1)<<std::endl;
}
