// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef CPPQEDELEMENTS_UTILS_QM_PICTURE_H_INCLUDED
#define CPPQEDELEMENTS_UTILS_QM_PICTURE_H_INCLUDED

#include "Pars.h"

#include <iosfwd>


enum QM_Picture {QMP_IP, QMP_UIP, QMP_SCH};


std::ostream& operator<<(std::ostream&, QM_Picture);
std::istream& operator>>(std::istream&, QM_Picture&);


namespace picture {
  
/// Convenience version of parameters::update meant to tackle the problem described in Sec. \ref masterequationlimitations
QM_Picture& updateWithPicture(parameters::Table& p, int argc, char* argv[], const std::string& prefix="--");

}

#endif // CPPQEDELEMENTS_UTILS_QM_PICTURE_H_INCLUDED
