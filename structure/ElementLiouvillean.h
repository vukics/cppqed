// -*- C++ -*-
#ifndef STRUCTURE_ELEMENTLIOUVILLEAN_H_INCLUDED
#define STRUCTURE_ELEMENTLIOUVILLEAN_H_INCLUDED

#include "ElementLiouvilleanFwd.h"

#include "Liouvillean.h"

#include "KeyPrinter.h"
#include "Range.h"

#ifndef   NDEBUG
#include "Exception.h"
#endif // NDEBUG

#include <boost/function.hpp>
#include <boost/bind.hpp>

#include <boost/preprocessor/control/expr_iif.hpp>
#include <boost/preprocessor/control/iif.hpp>
#include <boost/preprocessor/punctuation/comma_if.hpp>

namespace structure {


#define ISTD 0
#include "details/ElementLiouvilleanReentrant.h"
#undef  ISTD
#define ISTD 1
#include "details/ElementLiouvilleanReentrant.h"
#undef  ISTD


} // structure



#endif // STRUCTURE_ELEMENTLIOUVILLEAN_H_INCLUDED
