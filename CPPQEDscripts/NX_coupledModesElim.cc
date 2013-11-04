#include "Evolution.h"

#include "ParsMode.h"
#include "Particle.h"

#include "ElementLiouvillean.h"

#include "FFT.tcc"


using namespace std;
using namespace particle;
using namespace mathutils;

class NX_coupledModesElim 
  : private boost::base_from_member<Spatial>, 
    public structure::Free, 
    public Averaged, 
    public structure::Hamiltonian<1,structure::NO_TIME>/*, public structure::ElementLiouvillean<1,1>*/
{
public:
  typedef boost::base_from_member<Spatial> SpatialStorage;

  typedef particle::StateVectorLow StateVectorLow;

  typedef Spatial::Array DArray;
  typedef StateVectorLow CArray;

  NX_coupledModesElim(particle::Pars, const mode::ParsPumpedLossy&, double u);

  const Spatial& getSpace() const {return space_;}

private:
  void addContribution_v(const StateVectorLow&, StateVectorLow&) const;

  const Spatial& space_;

  const DArray& x_, y_;

  const CArray twoOmrecY_SqrOverI_, hOverI_X_;


};

/*

  array0_  2*omrec*y^2/I

  array1_  (2*omrec*x^2+eta^2/kappa*atan

 */


int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem

  ParameterTable p;

  ParsEvolution pe   (p); // Driver Parameters
  ParsPumped    ppart(p); 

  mode::ParsPumpedLossy pmode(p);

  double& u=p.add("u","",1.);

  // Parameter finalization
  update(p,argc,argv,"--");
  
  // ****** ****** ****** ****** ****** ******

  // if (pe.evol==EM_MASTER && qmp==QMP_IP) qmp=QMP_UIP;

  NX_coupledModesElim system(ppart,pmode,u);

  /*
  Ptr part(make(ppart,qmp));

  if (!ppart.init.getSig() && !ppart.vClass) {cerr<<"Incorrect initial condition"<<endl; abort();}
  */

  StateVector psi(wavePacket(ppart.init,system.getSpace()));

  evolve(psi,system,pe);




}




NX_coupledModesElim::NX_coupledModesElim(particle::Pars ppart, const mode::ParsPumpedLossy& pmode, double u) 
  : SpatialStorage(ppart.fin,sqrt(2*PI/(1<<ppart.fin))),
    Free(member.getDimension()),
    Averaged(member),
    space_(member),
    x_(member.getX()), y_(member.getK()),
    twoOmrecY_SqrOverI_(2*ppart.omrec*blitz::sqr(y_)/DCOMP_I),
    hOverI_X_((2*ppart.omrec*blitz::sqr(x_)+sqr(pmode.eta)/pmode.kappa*atan((u*x_-pmode.delta)/pmode.kappa))/DCOMP_I)
{
  getParsStream()<<"# NX_coupledModesElim\n";
  member.header(getParsStream());
}



void NX_coupledModesElim::addContribution_v(const StateVectorLow& psi, StateVectorLow& dpsidt) const
{
  dpsidt+=hOverI_X_*psi;

  StateVectorLow psiY(psi.copy());

  transform(psiY,fft::DIR_XK); psiY*=twoOmrecY_SqrOverI_; transform(psiY,fft::DIR_KX); dpsidt+=psiY;
}