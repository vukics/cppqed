// -*- C++ -*-
#ifndef UTILS_QM_PICTURE_H_INCLUDED
#define UTILS_QM_PICTURE_H_INCLUDED

#include "QM_PictureFwd.h"

#include<iosfwd>

std::ostream& operator<<(std::ostream&, QM_Picture);
std::istream& operator>>(std::istream&, QM_Picture&);

#endif // UTILS_QM_PICTURE_H_INCLUDED
