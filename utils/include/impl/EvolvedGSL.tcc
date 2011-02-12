// -*- C++ -*-
#ifndef   EVOLVED_GSL_IMPL_INCLUDED
#define   EVOLVED_GSL_IMPL_INCLUDED


namespace evolved {


template<typename A>
const typename Maker<A>::SmartPtr MakerGSL<A>::operator()(
							  A& a,
							  Derivs derivs,
							  double dtInit,
							  double epsRel,
							  double epsAbs,
							  const A& scaleAbs
							  ) const
{
  return SmartPtr(new details::EvolvedGSL<A>(a,derivs,dtInit,epsRel,epsAbs,scaleAbs));
}



namespace details {

template<typename A>
int auxFunction(double t, const double* y, double* dydt, void* aux)
{
  EvolvedGSL<A>* e=static_cast<EvolvedGSL<A>*>(aux);

  const A    yInterfaceA(cpputils::ArrayMemoryTraits<A>::create(y   ,e->getA()));
        A dydtInterfaceA(cpputils::ArrayMemoryTraits<A>::create(dydt,e->getA()));

  e->getDerivs()(t,yInterfaceA,dydtInterfaceA);

  return onSuccess;
}



template<typename A>
EvolvedGSL<A>::EvolvedGSL(
			  A& a,
			  Derivs derivs,
			  double dtInit,
			  double epsRel,
			  double epsAbs,
			  const A& scaleAbs
			  )
  : Base(a,derivs,dtInit,epsRel,epsAbs),
    pImpl_(createImpl(this,Traits::size(getA()),auxFunction<A>,epsRel,epsAbs,Traits::data(scaleAbs)))
{
  if (!Traits::isStorageContiguous(a)) throw (NonContiguousStorageException());
}



template<typename A> void EvolvedGSL<A>::step(double deltaT)
{
  double time=getTime(), dtTry=getDtTry();
  apply(pImpl_,&time,time+deltaT,&dtTry,Traits::data(getA()));
  Base::update(time,dtTry);
}


} // details


} // evolved

#endif // EVOLVED_GSL_IMPL_INCLUDED
