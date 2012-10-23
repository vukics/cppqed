// -*- C++ -*-
#ifndef STRUCTURE_IMPL_ELEMENTAVERAGED_TCC_INCLUDED
#define STRUCTURE_IMPL_ELEMENTAVERAGED_TCC_INCLUDED

#include "ElementAveraged.h"

#include "Algorithm.h"

#include <boost/iterator/transform_iterator.hpp>

#define TRANSFORMED_iterator(beginend) boost::make_transform_iterator(collection.beginend(),boost::bind(&Element::getLabels,_1))

template<int RANK, bool IS_TD>
structure::averaged::Collecting<RANK,IS_TD>::Collecting(const Collection& collection)
  : Base(collection.begin()->getTitle(),
	 cpputils::concatenateGrow(make_iterator_range(TRANSFORMED_iterator(begin),TRANSFORMED_iterator(end)),KeyLabels())),
    collection_(collection.clone())
{}


template<int RANK, bool IS_TD>
structure::averaged::Collecting<RANK,IS_TD>::Collecting(const Collecting& collecting)
  : Base(collecting),
    collection_(collecting.collection_.clone())
{}

#undef TRANSFORMED_iterator



#define TRANSFORMED_iterator(beginend) boost::make_transform_iterator(collection_.beginend(),boost::bind(&Element::average,_1,t,boost::cref(matrix)))

template<int RANK, bool IS_TD>
const typename structure::averaged::Collecting<RANK,IS_TD>::Averages
structure::averaged::Collecting<RANK,IS_TD>::average_v(double t, const LazyDensityOperator& matrix) const
{
  Averages res(nAvr()); res=0;
  return cpputils::concatenate(make_iterator_range(TRANSFORMED_iterator(begin),TRANSFORMED_iterator(end)),res);

}

#undef TRANSFORMED_iterator



template<int RANK, bool IS_TD>
void
structure::averaged::Collecting<RANK,IS_TD>::process_v(Averages& avr) const
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





#endif // STRUCTURE_IMPL_ELEMENTAVERAGED_TCC_INCLUDED
