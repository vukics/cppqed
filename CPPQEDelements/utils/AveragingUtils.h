// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef   CPPQEDELEMENTS_UTILS_AVERAGINGUTILS_H_INCLUDED
#define   CPPQEDELEMENTS_UTILS_AVERAGINGUTILS_H_INCLUDED

#include "ElementAveraged.h"

#include "DensityOperator.h"
#include "LazyDensityOperator.h"
#include "NegPT.tcc"

#include "DimensionsBookkeeper.h"

#include "Algorithm.h"
#include "MathExtensions.h"
#include "MultiIndexIterator.h"

#include <boost/assign/list_of.hpp>

#include <boost/ptr_container/ptr_list.hpp>

#include <boost/bind.hpp>

#include <boost/range/adaptor/transformed.hpp>

#include <algorithm>


template<int RANK>
class ReducedDensityOperator : private DimensionsBookkeeper<RANK>, public structure::ClonableElementAveraged<RANK>
{
private:
  typedef structure::ClonableElementAveraged<RANK> Base;

public:
  typedef typename Base::KeyLabels KeyLabels;
  
  typedef typename DimensionsBookkeeper<RANK>::Dimensions Dimensions;
  
  using DimensionsBookkeeper<RANK>::getDimensions; using DimensionsBookkeeper<RANK>::getTotalDimension;

  ReducedDensityOperator(const std::string& label, const Dimensions& dim, bool offDiagonals = false, const KeyLabels& subsequent = KeyLabels{}) :
    DimensionsBookkeeper<RANK>(dim),
    Base{label,[&] () {
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
      for (i.setToBegin(); i!=i.getEnd(); ++i) {
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
      }

      res.insert(res.end(),subsequent.begin(),subsequent.end());
        
      return res;
        
    } () },
    offDiagonals_(offDiagonals) {}

private:
  Base*const do_clone() const override {return new ReducedDensityOperator(*this);}

protected:
  const structure::Averages average_v(structure::NoTime, const quantumdata::LazyDensityOperator<RANK>& matrix) const override
  {
    return deflate(matrix,offDiagonals_);
  }

private:
  const bool offDiagonals_;

};


template<int RANK, typename V>
class ReducedDensityOperatorNegativity : public ReducedDensityOperator<RANK>
{
public:
  typedef ReducedDensityOperator<RANK> Base;
  typedef typename Base::Dimensions Dimensions;
  
  using Base::getDimensions; using Base::getTotalDimension;
  
  ReducedDensityOperatorNegativity(const std::string& label, const Dimensions& dim)
    : Base(label,dim,true,boost::assign::list_of("negativity")) {}

private:
  const structure::Averages average_v(structure::NoTime t, const quantumdata::LazyDensityOperator<RANK>& matrix) const override
  {
    auto res(Base::average_v(t,matrix));
    res.resize(res.size()+1);
    return res;
  }

  void process_v(structure::Averages& averages) const override
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

  
};


namespace averagingUtils {


using structure::Averages;

  
template<int RANK, bool IS_TIME_DEPENDENT=false>
class Collecting : public structure::ClonableElementAveraged<RANK,IS_TIME_DEPENDENT>
{
private:
  typedef structure::ClonableElementAveraged<RANK,IS_TIME_DEPENDENT> Base;

public:
  typedef Base Element;
  typedef boost::ptr_list<Element> Collection;

  typedef typename Base::KeyLabels KeyLabels;

  typedef typename Base::Time Time;
  
  Collecting(const Collection& collection)
    : Base(collection.front().getTitle(),
           cppqedutils::concatenateGrow<KeyLabels>( collection | boost::adaptors::transformed(bind(&Element::getLabels,_1)) )),
      collection_(collection.clone()) {}

  Collecting(const Collecting& collecting)
    : Base(collecting), collection_(collecting.collection_.clone()) {}

private:
  Base*const do_clone() const {return new Collecting(*this);}

  const Averages average_v(Time t, const quantumdata::LazyDensityOperator<RANK>& matrix) const override
  {
    Averages res(this->nAvr()); res=0;
    return cppqedutils::concatenate( collection_ | boost::adaptors::transformed(bind(&Element::average,_1,structure::time::convert(t),boost::cref(matrix))) , res);
  }

  void process_v(Averages& avr) const override
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

  const Collection collection_;

};



template<int RANKFROM, int RANKTO, bool IS_TIME_DEPENDENT=false>
class Transferring : public structure::AveragedTimeDependenceDispatched<RANKFROM,IS_TIME_DEPENDENT>
// Transfers the calculation of averages to another Averaged class, possibly with different RANK.
// The LazyDensityOperator for the other class should reference the same data.
// For a usage example cf. scripts/QbitMode_Matrix.cc
{
private:
  typedef structure::AveragedTimeDependenceDispatched<RANKFROM,IS_TIME_DEPENDENT> Base;
  
public:
  typedef typename Base::Time Time;
  
  using AveragedToPtr = structure::AveragedPtr<RANKTO> ;
  
  typedef quantumdata::LazyDensityOperator<RANKTO> LazyDensityOperatorTo;

  Transferring(AveragedToPtr averaged, const LazyDensityOperatorTo& ldo)
    : averaged_(averaged), ldo_(ldo) {}

private:
  std::ostream& streamKey_v(std::ostream& os, size_t& n) const override {return averaged_->streamKey(os,n);}

  void process_v(Averages& averages) const override {averaged_->process(averages);}

  std::ostream& stream_v(const Averages& averages, std::ostream& os, int n) const override {return averaged_->stream(averages,os,n);}

  size_t nAvr_v() const override {return averaged_->nAvr();}

  const Averages average_v(Time t, const quantumdata::LazyDensityOperator<RANKFROM>&) const override
  {
    return averaged_->average(structure::time::convert(t),ldo_);
  }

  const AveragedToPtr averaged_;
  const LazyDensityOperatorTo& ldo_;

};


} // averagingUtils


#endif // CPPQEDELEMENTS_UTILS_AVERAGINGUTILS_H_INCLUDED
