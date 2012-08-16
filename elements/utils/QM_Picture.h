// -*- C++ -*-
#ifndef ELEMENTS_UTILS_QM_PICTURE_H_INCLUDED
#define ELEMENTS_UTILS_QM_PICTURE_H_INCLUDED

#include "QM_PictureFwd.h"

#include<iosfwd>

std::ostream& operator<<(std::ostream&, QM_Picture);
std::istream& operator>>(std::istream&, QM_Picture&);

#endif // ELEMENTS_UTILS_QM_PICTURE_H_INCLUDED
