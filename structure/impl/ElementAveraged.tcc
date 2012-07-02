// -*- C++ -*-
#ifndef STRUCTURE_ELEMENT_AVERAGED_IMPL_INCLUDED
#define STRUCTURE_ELEMENT_AVERAGED_IMPL_INCLUDED


#include "Algorithm.h"

#include <boost/iterator/transform_iterator.hpp>

#define TRANSFORMED_iterator(beginend) boost::make_transform_iterator(collection.beginend(),boost::bind(getLabels,_1))

template<int RANK>
structure::averaged::Collecting<RANK>::Collecting(const Collection& collection)
  : Base(collection.begin()->getTitle(),
	 cpputils::concatenateGrow(make_iterator_range(TRANSFORMED_iterator(begin),TRANSFORMED_iterator(end)),KeyLabels())),
    collection_(collection.clone())
{}


template<int RANK>
structure::averaged::Collecting<RANK>::Collecting(const Collecting& collecting)
  : Base(collecting),
    collection_(collecting.collection_.clone())
{}

#undef TRANSFORMED_iterator



#define TRANSFORMED_iterator(beginend) boost::make_transform_iterator(collection_.beginend(),boost::bind(&Element::average,_1,boost::cref(ldo)))

template<int RANK>
const typename structure::averaged::Collecting<RANK>::Averages
structure::averaged::Collecting<RANK>::average(const LazyDensityOperator& ldo) const
{
  Averages res(nAvr()); res=0;
  return cpputils::concatenate(make_iterator_range(TRANSFORMED_iterator(begin),TRANSFORMED_iterator(end)),res);

}

#undef TRANSFORMED_iterator



template<int RANK>
void
structure::averaged::Collecting<RANK>::process(Averages& avr) const
{
  struct Helper
  {
    static void doIt(const Element& eav, Averages& avr, ptrdiff_t& l, ptrdiff_t& u)
    {
      using blitz::Range;
      if ((u=l+eav.nAvr())>l) {
	Averages temp(avr(Range(l+1,u)));
	eav.process(temp);
      }
      std::swap(l,u);
    }
  };
 
  ptrdiff_t l=-1, u=0;
  for_each(collection_,boost::bind(Helper::doIt,_1,avr,l,u));
}





#endif // STRUCTURE_ELEMENT_AVERAGED_IMPL_INCLUDED
