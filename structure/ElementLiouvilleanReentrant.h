// -*- C++ -*-

#define COND_ARG(IS_TD) BOOST_PP_EXPR_IIF(ISTD,double) BOOST_PP_COMMA_IF(ISTD)
#define COND_ARG_T(IS_TD) BOOST_PP_EXPR_IIF(ISTD,double t) BOOST_PP_COMMA_IF(ISTD)


template<int RANK, int NOJ>
class ElementLiouvillean<RANK,NOJ,BOOST_PP_IIF(ISTD,true,false)> : public Liouvillean<RANK,BOOST_PP_IIF(ISTD,true,false)>
{
public:
  typedef Liouvillean<RANK,BOOST_PP_IIF(ISTD,true,false)> Base;

  typedef typename Base::     StateVectorLow      StateVectorLow;

  typedef typename Base::LazyDensityOperator LazyDensityOperator;

  typedef typename LiouvilleanCommon::Probabilities Probabilities;

  typedef boost::function<void  (COND_ARG(ISTD)       StateVectorLow&     )> JumpStrategy;
  typedef boost::function<double(COND_ARG(ISTD) const LazyDensityOperator&)> JumpProbabilityStrategy;

  typedef blitz::TinyVector<JumpStrategy           ,NOJ> JumpStrategies;
  typedef blitz::TinyVector<JumpProbabilityStrategy,NOJ> JumpProbabilityStrategies;

protected:
  ElementLiouvillean(const JumpStrategies& jumps, const JumpProbabilityStrategies& jumpProbas) 
    : jumps_(jumps), jumpProbas_(jumpProbas) {}

private:
  size_t nJumps() const {return NOJ;}

  const Probabilities probabilities(COND_ARG(ISTD) const LazyDensityOperator&) const;

  void actWithJ(COND_ARG_T(ISTD) StateVectorLow& psi, size_t jumpNo) const
  {jumps_(jumpNo)(BOOST_PP_EXPR_IIF(ISTD,t) BOOST_PP_COMMA_IF(ISTD) psi);}

  const JumpStrategies            jumps_     ;
  const JumpProbabilityStrategies jumpProbas_;

};


template<int RANK>
class ElementLiouvillean<RANK,1,BOOST_PP_IIF(ISTD,true,false)> : public Liouvillean<RANK,BOOST_PP_IIF(ISTD,true,false)>
// This specialization can use the virtual-function technique of old
{
public:
  typedef Liouvillean<RANK,BOOST_PP_IIF(ISTD,true,false)> Base;

  typedef typename Base::    StateVectorLow     StateVectorLow;

  typedef typename Base::LazyDensityOperator LazyDensityOperator;

  typedef typename LiouvilleanCommon::Probabilities Probabilities;

private:
  size_t nJumps() const {return 1;}

  const Probabilities probabilities(COND_ARG_T(ISTD) const LazyDensityOperator& matrix) const 
  {Probabilities probas(1); probas(0)=probability(BOOST_PP_EXPR_IIF(ISTD,t) BOOST_PP_COMMA_IF(ISTD) matrix); return probas;}

#ifndef   NDEBUG
  struct ElementLiouvilleanException : cpputils::Exception {};
#endif // NDEBUG

  void actWithJ(COND_ARG_T(ISTD) StateVectorLow& psi, size_t 
#ifndef   NDEBUG
		  jumpNo
#endif // NDEBUG
		  ) const {
#ifndef   NDEBUG
    if (jumpNo) throw ElementLiouvilleanException(); 
#endif // NDEBUG
    doActWithJ(BOOST_PP_EXPR_IIF(ISTD,t) BOOST_PP_COMMA_IF(ISTD) psi);
  }


  virtual void   doActWithJ (COND_ARG(ISTD)       StateVectorLow     &) const = 0;
  virtual double probability(COND_ARG(ISTD) const LazyDensityOperator&) const = 0;

};




template<int RANK, int NOJ>
const LiouvilleanCommon::Probabilities ElementLiouvillean<RANK,NOJ,BOOST_PP_IIF(ISTD,true,false)>::probabilities(COND_ARG_T(ISTD) const LazyDensityOperator& matrix) const
{
  Probabilities probas(NOJ);
  // Note that this cannot be anything like static because of the by-reference semantics of blitz::Array

  boost::transform(jumpProbas_,probas.begin(),
		   bind(&JumpProbabilityStrategy::operator(),_1,BOOST_PP_EXPR_IIF(ISTD,t) BOOST_PP_COMMA_IF(ISTD) boost::cref(matrix)));
  return probas;
}


#undef COND_ARG
#undef COND_ARG_T

