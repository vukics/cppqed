#include "impl/Composite.tcc"

const int v=composite::result_of::Make<Act<0,1>,Act<0,2>,Act<0,3>,Act<3,2,5,1> >::type::RANK ; 
