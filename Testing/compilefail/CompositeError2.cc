#include "Composite.tcc"

// Does not compile: (NEGATIVE_ELEMENT_in_NONNEGATIVE_VECTOR)

template struct composite::result_of::Make<Act<1,0>,Act<2,0>,Act<3,4>,Act<1,2,0>,Act<2,-1,0>,Act<3,2,0,1> >;
