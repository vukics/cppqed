// -*- C++ -*-
#ifndef   STRUCTURE_LIOUVILLEANFWD_H_INCLUDED
#define   STRUCTURE_LIOUVILLEANFWD_H_INCLUDED


namespace structure {


class LiouvilleanCommon;

#ifndef   NDEBUG
struct LiouvilleanNumberMismatchException;
#endif // NDEBUG

template<int, bool IS_TD=true>
class Liouvillean;


} // structure


#endif // STRUCTURE_LIOUVILLEANFWD_H_INCLUDED
