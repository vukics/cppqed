// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef   CPPQEDCORE_COMPOSITES_COMPOSITE_TCC_INCLUDED
#define   CPPQEDCORE_COMPOSITES_COMPOSITE_TCC_INCLUDED


// template<int>
struct CompositeConsistencyException : std::logic_error
{
  CompositeConsistencyException(int idxx, int ii) : std::logic_error(std::to_string(idxx)+" "+std::to_string(ii)), idx(idxx), i(ii) {}

  const int idx, i;
};






//////////////
//
// Liouvillean
//
//////////////




////////
//
// Maker
//
////////


#define DISPATCHER(EX,HA,LI) (all(systemCharacteristics==SystemCharacteristics{EX,HA,LI})) return std::make_shared<Composite<Acts,EX,HA,LI> >(frees,acts)

template<typename Acts>
const typename composite::Base<Acts>::Ptr composite::doMake(const Acts& acts)
{
  using namespace structure;
  
  const typename Base<Acts>::Frees frees(fillFrees(acts));

  SystemCharacteristics systemCharacteristics{false,false,false};
  
  mpl::for_each<typename Base<Acts>::Ordinals>([&](auto t) -> void {
    const SubSystemFree& free=frees[std::decay_t<decltype(t)>::value];
    systemCharacteristics|=SystemCharacteristics{free.getEx()!=0,free.getHa()!=0,free.getLi()!=0};
  });
  
  systemCharacteristics=boost::fusion::fold(acts,systemCharacteristics,[&](SystemCharacteristics sc, const auto& act) -> SystemCharacteristics {
    return sc || SystemCharacteristics{act.getEx()!=0,act.getHa()!=0,act.getLi()!=0};
  });
    
  if      DISPATCHER(true ,true ,true ) ;
  else if DISPATCHER(true ,true ,false) ;
  else if DISPATCHER(true ,false,true ) ;
  else if DISPATCHER(true ,false,false) ;
  else if DISPATCHER(false,true ,true ) ;
  else if DISPATCHER(false,true ,false) ;
  else if DISPATCHER(false,false,true ) ;
  else return std::make_shared<Composite<Acts,false,false,false> >(frees,acts);
}


#undef DISPATCHER

#undef CALL_composite_worker
#endif // CPPQEDCORE_COMPOSITES_COMPOSITE_TCC_INCLUDED
