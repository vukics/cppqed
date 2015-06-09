// Copyright András Vukics 2006–2015. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-
#ifndef CPPQEDELEMENTS_UTILS_QM_PICTURE_H_INCLUDED
#define CPPQEDELEMENTS_UTILS_QM_PICTURE_H_INCLUDED

#include "QM_PictureFwd.h"

#include "Pars.h"

#include<iosfwd>

std::ostream& operator<<(std::ostream&, QM_Picture);
std::istream& operator>>(std::istream&, QM_Picture&);


namespace picture {
/// Convenience version of parameters::update meant to tackle the problem described in Sec. \ref masterequationlimitations
QM_Picture& updateWithPicture(parameters::ParameterTable& p, int argc, char* argv[], const std::string& prefix="--");
}

#endif // CPPQEDELEMENTS_UTILS_QM_PICTURE_H_INCLUDED
