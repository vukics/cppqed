// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{The parameter-bundle}
#ifndef CPPQEDCORE_UTILS_PARS_H_INCLUDED
#define CPPQEDCORE_UTILS_PARS_H_INCLUDED

#include "ParsFwd.h"

#include "BooleanNegatedProxy.h"

#include <boost/ptr_container/ptr_list.hpp>

#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <typeinfo>


/// The parameter-bundle
namespace parameters {


struct Empty
{
  Empty(Table&, const std::string& ="") {}
};


inline const size_t maxTypeLabelLength=24;


//////////////////
// 
// Pars Exceptions
//
//////////////////


/// Thrown if a certain Typed name is not found in Table when \link Table::operator[] subscripting\endlink.
struct UnrecognisedParameterException : std::runtime_error {using std::runtime_error::runtime_error;};

/// Thrown if a given Typed name is attempted to be used more than once when \link Table::add pushing\endlink to Table
struct AttemptedRecreationOfParameterException : std::runtime_error {using std::runtime_error::runtime_error;};


/// Common, template-parameter independent base class of Typed, to be stored in a pointer-container in Table
/**
 * \note This class is not really part of the API, it is included here simply because \link Typed parameters\endlink are stored in Table as a pointer container (`boost::ptr_list<Base>`).
 * The inclusion here could be avoided by resorting to a pointer-to-implementation technique in Table, however, … \see … Typed
 */
class Base
{
public:
  Base(const std::string& s, ///< parameter name (as appears in command line)
                const std::string& d  ///< parameter description (to be listed in help)
               ) : s_(s), d_(d) {}
  virtual ~Base() {}

  const std::string& getS() const {return s_;}
  const std::string& getD() const {return d_;}

  void print(size_t smw, size_t tmw, size_t dmw) const {print_v(smw,tmw,dmw);}
  void read(std::istream& is) {read_v(is);}

protected:
  virtual void  read_v(std::istream&)                = 0;

private:
  virtual void print_v(size_t, size_t, size_t) const = 0;

  const std::string s_; // as appearing in command line
  const std::string d_; // short description

};


/// Template containing value for the given parameter
/** \note This class is not really part of the API, it is included here simply in order that specializations of it might be declared elsewhere in the framework */
template<typename T>
class Typed : public Base 
{
public:
  Typed(const std::string& s, ///< parameter name (as appears in command line)
            const std::string& d, ///< parameter description (to be listed in help)
            const T& v ///< parameter default value
           ) : Base(s,d), v_(v) {}

  const T& get() const {return v_;}
        T& get()       {return const_cast<T&>(static_cast<const Typed*>(this)->get());}

protected:
  void read_v(std::istream& is) {is>>v_;}

private:
  void print_v(size_t smw, size_t tmw, size_t dmw) const
  {
    using namespace std;
    const string name(typeid(T).name());
    cout<<setw(smw+3)<<left<<getS()
        <<setw(tmw+3)<<left<<(name.size()>maxTypeLabelLength-3 ? name.substr(0,maxTypeLabelLength-3)+"..." : name)
        <<setw(dmw+3)<<left<<getD()<<v_<<endl;
  }
  
  T v_;
  
};


/// The parameter table according to which the command line will be parsed by update()
/** \see \ref userguideelementaryparameters in the \ref userguide */
class Table
{
public:
  typedef std::shared_ptr<std::string> ParsedCommandLine; ///< Type for passing the \link getParsedCommandLine parsed command line\endlink throughout the framework

  Table();

  /// \name Subscription
  //@{
  const Base& operator[](const std::string&  ) const; ///< Subscription by parameter name
        Base& operator[](const std::string& s)        {return const_cast<Base&>(static_cast<const Table&>(*this)[s]);} ///< ” non-const version
  //@}

  /// \name Adders
  //@{
  /// generic version
  /**
   * In function scope, if `p` is a Table, then, e.g. a call like
   * 
   *       double & omega=p.add("omega","harmonic oscillator mode frequency",1.);
   * 
   * will add a parameter of type `double` to the Table, so that if on the highest level of the application the command line gets parsed with the update() function like
   * 
   *       int main(int argc, char* argv[])
   *       {
   *         // Table p declared and filled with omega with p.add(...) as above
   * 
   *         // parsing the command line according to p
   *         update(p,argc,argv,"--");
   *         
   *         // ... operations involving omega
   *       }
   * 
   * then, when starting the program from the command line, the value of `omega` can be specified with e.g. `--omega 10.`, while if unspecified, the value will be the default `1`.
   * When starting the program with `--help`, `omega`, with all other parameters will be listed in a nice tabulated format with the specified description and default value.
   * 
   * In addition, the reference `omega` serves as a convenient handle for further operations involving omega in the program.
   * 
   * \see e.g. the implementation of trajectory::ParsEvolved::ParsEvolved()
   * 
   */
  template<typename T> 
  T& add(const std::string& s, ///< parameter name (as appears in command line)
         const std::string& d, ///< parameter description (to be listed in help)
         const T& v ///< parameter default value
        )
  {
    using namespace std;
    auto i=findSubscript(s);
    if (i==table_.end()) {
      Typed<T>* pptr=new Typed<T>(s,d,v);
      table_.push_back(pptr);
      smwidth_=max(smwidth_,s.length());
      tmwidth_=max(tmwidth_,min(strlen(typeid(T).name()),maxTypeLabelLength));
      dmwidth_=max(dmwidth_,d.length());
      return pptr->get();
    }
    else throw AttemptedRecreationOfParameterException(s);
  }
  
  /// Overload of add(const std::string &, const std::string &, const T &) for boolean
  /**
   * This will create a switch, that does not require any value. It will also add to the Table an additional parameter of type BooleanNegatedProxy which is the opposite switch with a `no_` prefix.
   * 
   * \see trajectory::ParsRun::streamInfo – if `--streamInfo` is present in the command line, the switch will have `true` value, but if `--no_streamInfo` is also present subsequently, the value will be `false`.
   */
  bool& add(const std::string& s, const std::string& d, bool v);

  /// Adds the parameter with a modifier suffix – useful when many parameters are distinguished only by suffixes (e.g. numerical suffix)
  /** \see  \ref userguidemorecomplexring "This" and \ref userguidemorecomplexmultiparticle "this" section of the \ref userguide */
  template<typename T> T& add(const std::string& s, const std::string& mod, const std::string& d, const T& v)
  {
    return add(s+mod,d,v);
  }

  /// This adds a dummy “parameter” whose only effect is to cause a newline and a title to be printed for grouping parameters belonging to different modules in parameter \link printList listing\endlink
  Table& addTitle(const std::string& s, const std::string& mod="");
  //@}
  
  /// Print a full list of Parameters in help with name, typeID, description, and (default) value
  /** Invoked by update() if the switch `--help` is found in the command line */
  void printList() const;
  
  /// Getter for the full command line set by the update() function
  const ParsedCommandLine getParsedCommandLine(                    ) const {return  parsedCommandLine_;}

  /// Setter for the full command line used by the update() function
  void                    setParsedCommandLine(const std::string& s)       {        *parsedCommandLine_=s;}

private:
  typedef boost::ptr_list<Base> Impl;

  Impl::const_iterator findSubscript(const std::string& s) const;

  Impl table_;

  size_t smwidth_; // maximal width of s_ entries
  size_t tmwidth_; // maximal width of typeIDs
  size_t dmwidth_; // maximal width of d_ entries

  ParsedCommandLine parsedCommandLine_;

};


/// The function parsing the command line and putting the parameters into the Table \related Table
void update(Table& p,
            int argc, ///< number of words in the command line
            char* argv[], ///< content of the command line word by word
            const std::string& prefix="--" ///< prefix of parameters in the command line
           );


} // parameters


#endif // CPPQEDCORE_UTILS_PARS_H_INCLUDED
