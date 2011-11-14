// -*- C++ -*-
#ifndef _QM_PICTURE_H
#define _QM_PICTURE_H

#include "QM_PictureFwd.h"

#include<iosfwd>

std::ostream& operator<<(std::ostream&, QM_Picture);
std::istream& operator>>(std::istream&, QM_Picture&);

#endif // _QM_PICTURE_H
