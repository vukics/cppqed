// Copyright Raimar Sandner 2012–2014. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-

#ifndef CPYPYQED_INCLUDE_PARSPROPERTYMACROS_H_INCLUDED
#define CPYPYQED_INCLUDE_PARSPROPERTYMACROS_H_INCLUDED

#define PARS_GETTER(type, Class, val) type get##val(const Class *p){return p->val;}

#define PARS_GETTER_SETTER(type, Class, val) \
  PARS_GETTER(type, Class, val)\
  void set##val(Class *p, type v){p->val=v;}
  
#define PARS_PROPERTY(val) add_property(#val, &get##val, &set##val)
#define PARS_RO_PROPERTY(val) add_property(#val, &get##val)

#endif // CPYPYQED_INCLUDE_PARSPROPERTYMACROS_H_INCLUDED