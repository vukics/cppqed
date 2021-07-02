// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

#ifndef   CPPQEDELEMENTS_UTILS_AVERAGINGUTILS_TCC_INCLUDED
#define   CPPQEDELEMENTS_UTILS_AVERAGINGUTILS_TCC_INCLUDED

#include "AveragingUtils.h"

#include "DensityOperator.h"
#include "LazyDensityOperator.h"
#include "NegPT.tcc"

#include "Algorithm.h"
#include "MathExtensions.h"
#include "MultiIndexIterator.h"


#include <boost/bind.hpp>

#include <boost/range/adaptor/transformed.hpp>

#include <algorithm>


template<int RANK>
auto
ReducedDensityOperator<RANK>::helper(const Dimensions& dim, bool offDiagonals, const KeyLabels& subsequent) -> const KeyLabels
{
  typedef cppqedutils::MultiIndexIterator<RANK> Iterator;
  using std::ostringstream;

  KeyLabels res;
  Iterator i(dim-1,cppqedutils::mii::begin);

  // Diagonals
  for (; i!=i.getEnd(); ++i) {
    ostringstream ss;
    ss<<"rho_"<<*i<<';'<<*i;
    res.push_back(ss.str());
  }

  if (!offDiagonals) return res;

  // Offdiagonals
  for (i.setToBegin(); i!=i.getEnd(); ++i) 
    for (Iterator j(std::next(i)); j!=j.getEnd(); ++j) {
      {
        ostringstream ss;
        ss<<"real["<<"rho_"<<*i<<';'<<*j<<']';
        res.push_back(ss.str());
      }
      {
        ostringstream ss;
        ss<<"imag["<<"rho_"<<*i<<';'<<*j<<']';
        res.push_back(ss.str());
      }
  }

  res.insert(res.end(),subsequent.begin(),subsequent.end());

  return res;

}


template<int RANK>
ReducedDensityOperator<RANK>::ReducedDensityOperator(const std::string& label, const Dimensions& dim, bool offDiagonals, const KeyLabels& subsequent) :
  DimensionsBookkeeper<RANK>(dim),
  Base(label,helper(getDimensions(),offDiagonals,subsequent)),
  offDiagonals_(offDiagonals)
{
}


template<int RANK>
const structure::Averages 
ReducedDensityOperator<RANK>::average_v(structure::NoTime, const LazyDensityOperator& matrix) const
{
  return deflate(matrix,offDiagonals_);
}


template<int RANK, typename V>
const structure::Averages 
ReducedDensityOperatorNegativity<RANK,V>::average_v(structure::NoTime t, const LazyDensityOperator& matrix) const
{
  auto res(Base::average_v(t,matrix));
  res.resize(res.size()+1);
  return res;
}


template<int RANK, typename V>
void ReducedDensityOperatorNegativity<RANK,V>::process_v(structure::Averages& averages) const
{
  quantumdata::DensityOperator<RANK> rho(getDimensions());
  linalg::CMatrix matrix(rho.matrixView()); // references the same data
  const size_t dim=getTotalDimension();
  
  {
    int idx=0;  
    for (int i=0; i<dim; ++i, ++idx) matrix(i,i)=averages(idx);
    for (int i=0; i<dim; ++i) for (int j=i+1; j<dim; ++j, idx+=2) matrix(j,i)=conj(matrix(i,j)=dcomp(averages(idx),averages(idx+1)));
  }
  
  averages(averages.size()-1)=
#ifndef   DO_NOT_USE_FLENS
    quantumdata::negPT(rho,V())
#else  // DO_NOT_USE_FLENS
    0
#endif // DO_NOT_USE_FLENS
    ;
  
}


template<int RANK, bool IS_TIME_DEPENDENT>
averagingUtils::Collecting<RANK,IS_TIME_DEPENDENT>::Collecting(const Collection& collection)
  : Base(collection.front().getTitle(),
         cppqedutils::concatenateGrow<KeyLabels>( collection | boost::adaptors::transformed(bind(&Element::getLabels,_1)) )),
    collection_(collection.clone())
{}


template<int RANK, bool IS_TIME_DEPENDENT>
averagingUtils::Collecting<RANK,IS_TIME_DEPENDENT>::Collecting(const Collecting& collecting)
  : Base(collecting),
    collection_(collecting.collection_.clone())
{}


namespace averagingUtils { namespace details {

double convert(structure::OneTime t) {return t ;}
double convert(structure:: NoTime  ) {return 0.;}

} } // averagingUtils::details


template<int RANK, bool IS_TIME_DEPENDENT>
auto
averagingUtils::Collecting<RANK,IS_TIME_DEPENDENT>::average_v(Time t, const LazyDensityOperator& matrix) const -> const Averages
{
  Averages res(this->nAvr()); res=0;
  return cppqedutils::concatenate( collection_ | boost::adaptors::transformed(bind(&Element::average,_1,details::convert(t),boost::cref(matrix))) , res);

}


template<int RANK, bool IS_TIME_DEPENDENT>
void
averagingUtils::Collecting<RANK,IS_TIME_DEPENDENT>::process_v(Averages& avr) const
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


template<int RANKFROM, int RANKTO, bool IS_TIME_DEPENDENT>
auto
averagingUtils::Transferring<RANKFROM,RANKTO,IS_TIME_DEPENDENT>::average_v(Time t, const LazyDensityOperator&) const -> const Averages
{
  return averaged_->average(details::convert(t),ldo_);
}


/* The earlier solution fails with the C++11 standard due to a "bug" in Boost.Assign, cf. https://svn.boost.org/trac/boost/ticket/7364

template<int RANK>
struct ReducedDensityOperator<RANK>::Helper
{
  typedef cppqedutils::MultiIndexIterator<RANK> Iterator;
  
  Helper(const Dimensions& dim) : i_(Dimensions(size_t(0)),dim-1,cppqedutils::mii::begin), j_(i_), real_(true), offDiagonals_(false) {++i_; ++j_;}
  
  Helper() : i_(Dimensions(size_t(0)),Dimensions(size_t(0)),cppqedutils::mii::begin), j_(i_), real_(true), offDiagonals_(false) {}

  const std::string operator()()
  {
    using namespace std;
    stringstream ss(stringstream::out);

    // Diagonals
    if (!offDiagonals_ && i_!=i_.getEnd()) {
      ss<<"rho_"<<*i_<<';'<<*i_;
      ++i_;
    }
    else if (!offDiagonals_ && i_==i_.getEnd()) {
      offDiagonals_=true;
      i_.setToBegin();
    }
    
    // Offdiagonals
    if (offDiagonals_) {
      ss<<(real_ ? "real[" : "imag[")<<"rho_"<<*i_<<','<<*j_<<']';
      if (real_)
        real_=false;
      else {
        real_=true;
        ++j_;
        if (j_==j_.getEnd())
          ++(j_=++i_);
      }
    }
    
    return ss.str();
  }

private:
  Iterator i_, j_;
  bool real_, offDiagonals_;
  
};



template<int RANK>
ReducedDensityOperator<RANK>::ReducedDensityOperator(const std::string& label, const Dimensions& dim, bool offDiagonals, const KeyLabels& subsequent) :
  DimensionsBookkeeper<RANK,true>(dim),
  Base(label,boost::assign::repeat_fun((offDiagonals ? mathutils::sqr(getTotalDimension()) : getTotalDimension())-1,
                                       Helper(getDimensions())).range(subsequent)),
  offDiagonals_(offDiagonals)
{
}

*/

#endif // CPPQEDELEMENTS_UTILS_AVERAGINGUTILS_TCC_INCLUDED
