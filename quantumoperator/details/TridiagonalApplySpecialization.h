// -*- C++ -*-

#define rank BOOST_PP_ITERATION()

#define AT_helper(V,i) ranges(i)(at_c<V,i>::type::value)

#define DPSIDT_print(z,m,data) AT_helper(V_DPSIDT,m)
#define      A_print(z,m,data) AT_helper(V_A     ,m)
#define    PSI_print(z,m,data) AT_helper(V_PSI   ,m)


#define INDEXER_print2(z,m,data) boost::fusion::at_c<m>(Base::cache_)

template<> template<int START, typename V_DPSIDT, typename V_A, typename V_PSI>
void quantumoperator::Tridiagonal<rank>::doApply(mpl::int_<0>,
						 const Ranges& ranges, const StateVectorLow& psi, StateVectorLow& dpsidt) const
{
  // std::cerr<<START<<std::endl;

  using mpl::at_c;
  if (diagonals_(START).size()) 
    dpsidt             (BOOST_PP_ENUM(rank,DPSIDT_print,~))
      +=
      diagonals_(START)(BOOST_PP_ENUM(rank,     A_print,~))
      *
      psi              (BOOST_PP_ENUM(rank,   PSI_print,~))
      ;
}

#undef    PSI_print
#undef      A_print
#undef DPSIDT_print

#undef AT_helper

#undef rank
