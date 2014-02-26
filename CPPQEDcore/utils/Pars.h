/// \briefFile{The parameter-bundle}
// -*- C++ -*-
#ifndef CPPQEDCORE_UTILS_PARS_H_INCLUDED
#define CPPQEDCORE_UTILS_PARS_H_INCLUDED

#include "ParsFwd.h"

#include "BooleanNegatedProxy.h"
#include "Exception.h"

#include <boost/ptr_container/ptr_list.hpp>

#include <string>
#include <sstream>


/// The parameter-bundle
namespace parameters {

const size_t maxTypeLabelLength=24;


/// The function parsing the command line and putting the parameters into the ParameterTable \related ParameterTable
/** \todo Somehow echo the whole command line into the output file – difficulty: the parameter-bundle has no knowledge about the output file */
void update(ParameterTable& p,
            int argc, ///< number of words in the command line
            char* argv[], ///< content of the command line word by word
            const std::string& prefix="--" ///< prefix of parameters in the command line
           );


//////////////////
// 
// Pars Exceptions
//
//////////////////


class Exception : public cpputils::Exception
{
public:
  virtual ~Exception() {}
};


class NamedException : public Exception 
{
public:
  NamedException(const std::string& name);

  const std::string& getName() const {return name_;}

private:
  const std::string name_;
};


/// Thrown if a certain Parameter name is not found in ParameterTable when \link ParameterTable::operator[] subscripting\endlink.
class UnrecognisedParameterException : public NamedException 
{
public: 
  UnrecognisedParameterException(const std::string& name) : NamedException(name) {}
};


/// Thrown if a given Parameter name is attempted to be used more than once when \link ParameterTable::add pushing\endlink to ParameterTable
class AttemptedRecreationOfParameterException : public NamedException 
{
public: 
  AttemptedRecreationOfParameterException(const std::string& name);
};



/// Common, template-parameter independent base class of Parameter, to be stored in a pointer-container in ParameterTable
/**
 * \note This class is not really part of the API, it is included here simply because \link Parameter parameters\endlink are stored in ParameterTable as a pointer container (`boost::ptr_list<ParameterBase>`).
 * The inclusion here could be avoided by resorting to a pointer-to-implementation technique in ParameterTable, however, … \see … Parameter
 */
class ParameterBase
{
public:
  ParameterBase(const std::string& s, ///< parameter name (as appears in command line)
                const std::string& d  ///< parameter description (to be listed in help)
               ) : s_(s), d_(d) {}
  virtual ~ParameterBase() {}

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
class Parameter : public ParameterBase 
{
public:
  Parameter(const std::string& s, ///< parameter name (as appears in command line)
            const std::string& d, ///< parameter description (to be listed in help)
            const T& v ///< parameter default value
           ) : ParameterBase(s,d), v_(v) {}
  const T& get() const {return v_;}
        T& get()       {return const_cast<T&>(static_cast<const Parameter*>(this)->get());}
  
  void print(size_t smw, size_t tmw, size_t dmw) const;
  void read(std::istream& is) {is>>v_;}


protected:
  void read_v(std::istream& is) {is>>v_;}

private:
  void print_v(size_t smw, size_t tmw, size_t dmw) const;
  
  T v_;
  
};


/// The parameter table according to which the command line will be parsed by update()
/** \see \ref userguidescriptselementaryparameters in the \ref userguide */
class ParameterTable
{
public:
  ParameterTable() : table_(), smwidth_(0), tmwidth_(6), dmwidth_(0), stream_() {} // tmwidth_ cf bool!

  /// \name Subscription
  //@{
  const ParameterBase& operator[](const std::string&  ) const; ///< Subscription by parameter name
        ParameterBase& operator[](const std::string& s)        {return const_cast<ParameterBase&>(static_cast<const ParameterTable&>(*this)[s]);} ///< ” non-const version
  //@}

  /// \name Adders
  //@{
  /// generic version
  /**
   * In function scope, if `p` is a ParameterTable, then, e.g. a call like
   * 
   *       double & omega=p.add("omega","harmonic oscillator mode frequency",1.);
   * 
   * will add a parameter of type `double` to the ParameterTable, so that if on the highest level of the application the command line gets parsed with the update() function like
   * 
   *       int main(int argc, char* argv[])
   *       {
   *         // ParameterTable p declared and filled with omega with p.add(...) as above
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
  template<typename T> T& add(const std::string& s, ///< parameter name (as appears in command line)
                              const std::string& d, ///< parameter description (to be listed in help)
                              const T& v ///< parameter default value
                             );
  
  
  /// Overload of add(const std::string &, const std::string &, const T &) for boolean
  /**
   * This will create a switch, that does not require any value. It will also add to the ParameterTable an additional parameter of type BooleanNegatedProxy which is the opposite switch with a `no_` prefix.
   * 
   * \see trajectory::ParsRun::displayInfo – if `--displayInfo` is present in the command line, the switch will have `true` value, but if `--no_displayInfo` is also present subsequently, the value will be `false`.
   */
  bool& add(const std::string& s, const std::string& d, bool v);

  /// Adds the parameter with a modifier suffix – useful when many parameters are distinguished only by suffixes (e.g. numerical suffix)
  /** \see  \ref userguidescriptsmorecomplexring "This" and \ref userguidescriptsmorecomplexmultiparticle "this" section of the \ref userguide */
  template<typename T> T& addMod(const std::string& s, const std::string& mod, const std::string& d, const T& v)
  {
    return add(s+mod,d,v);
  }

  /// This adds a dummy “parameter” whose only effect is to cause a newline and a title to be printed for grouping parameters belonging to different modules
  ParameterTable& addTitle(const std::string& s, const std::string& mod="");
  //@}
  
  /// Print a full list of Parameters in help with name, typeID, description, and (default) value
  /** Invoked by update() if the switch `--help` is found in the command line */
  void printList() const;
  
  /// The stream whereon parameter i/o occurs (an intermediary \refStdCppConstruct{stringstream,sstream/stringstream/}) is exposed in order that i/o manipulators can be applied
  std::iostream& getStream() {return stream_;}

private:
  typedef boost::ptr_list<ParameterBase> Impl;

  Impl table_;

  size_t smwidth_; // maximal width of s_ entries
  size_t tmwidth_; // maximal width of typeIDs
  size_t dmwidth_; // maximal width of d_ entries

  std::stringstream stream_;
  
};


} // parameters


#endif // CPPQEDCORE_UTILS_PARS_H_INCLUDED
