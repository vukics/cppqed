#include "GeneralDicke.h"

#include "impl/Tridiagonal.tcc"

#include<boost/assign/list_of.hpp>


using namespace boost::assign;

using namespace structure;
using namespace spin;
using namespace mode;



generaldicke::Base::Base(mode::Ptr mode, spin::Ptr spin, dcomp u, dcomp y)
  : Interaction<2>(Frees(mode,spin),
		   RealFreqs(),
		   tuple_list_of("u",u,mode->getDimension()*spin->getDimension())("y",y,
										  sqrt(mode->getDimension()*spin->getDimension()))
		   ),
    TridiagonalHamiltonian<2,true>((u*nop(mode)*(sz(spin)+spin->getTwoS()/2.*quantumoperator::identity(spin->getDimension()))
				   +
				    y*(aop(mode).dagger()+aop(mode))*sx(spin))/DCOMP_I)
{
  getParsStream()<<"# General Dicke interaction\n";
  // getParsStream()<<sx(spin->getDimension()-1)<<std::endl;
}
