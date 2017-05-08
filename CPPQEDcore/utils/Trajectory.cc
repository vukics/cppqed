// Copyright András Vukics 2006–2017. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "ArrayTraits.h"
#include "BlitzArray.h"
#include "Trajectory.tcc"
#include "SmartPtr.h"

#include "ParsTrajectory.h"
#include "FormDouble.tcc"

#include "core_config.h"

#ifndef DO_NOT_USE_BOOST_COMPRESSION
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/device/file.hpp>
#endif // DO_NOT_USE_BOOST_COMPRESSION
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <iostream>
#include <list>


using namespace std;


namespace {

bool isbz2(const std::string filename)
{
  using namespace std;
  ifstream file(filename, ios_base::in | ios_base::binary);
  if (file.peek() == ifstream::traits_type::eof())
    return false;
  string header; header.resize(3);
  file.read(&header[0],3);
  file.close();
  return header=="BZh";
}

}

void trajectory::run(Trajectory & traj, double time, double deltaT, unsigned sdf, const std::string& ofn, const std::string& initialFileName, int precision,
                     bool displayInfo, bool firstStateDisplay,
                     double autoStopEpsilon, unsigned autoStopRepetition, const std::string& parsedCommandLine)
{details::run(traj,time,deltaT,sdf,ofn,initialFileName,precision,displayInfo,firstStateDisplay,autoStopEpsilon,autoStopRepetition,parsedCommandLine);}

void trajectory::run(Trajectory & traj, long   nDt , double deltaT, unsigned sdf, const std::string& ofn, const std::string& initialFileName, int precision,
                     bool displayInfo, bool firstStateDisplay,
                     double autoStopEpsilon, unsigned autoStopRepetition, const std::string& parsedCommandLine)
{details::run(traj,nDt ,deltaT,sdf,ofn,initialFileName,precision,displayInfo,firstStateDisplay,autoStopEpsilon,autoStopRepetition,parsedCommandLine);}

void trajectory::run(Trajectory& traj, const ParsRun& p)
{ 
  if (!p.Dt) {cerr<<"Nonzero Dt required!"<<endl; return;}
  if (p.NDt) run(traj,p.NDt,p.Dt,p.sdf,p.ofn,p.initialFileName,p.precision,p.displayInfo,p.firstStateDisplay,p.autoStopEpsilon,p.autoStopRepetition,p.getParsedCommandLine());
  else       run(traj,p.T  ,p.Dt,p.sdf,p.ofn,p.initialFileName,p.precision,p.displayInfo,p.firstStateDisplay,p.autoStopEpsilon,p.autoStopRepetition,p.getParsedCommandLine());
}


ostream& trajectory::Trajectory::display(ostream& os, int precision) const
{
  const FormDouble fd(formdouble::positive(precision));
  return display_v( os<<fd(getTime())<<fd(getDtDid()) , precision)<<endl; // Note: endl flushes the buffer
}


ostream& trajectory::Trajectory::displayParameters(ostream& os) const
{
  size_t i=3;
  return displayKey_v(displayParameters_v(os)<<endl<<"# Key to data:\n# Trajectory\n#  1. time\n#  2. dtDid\n" , i);
}

boost::shared_ptr<istream> trajectory::openStateFileReading(const std::string &filename)
{
#ifdef DO_NOT_USE_BOOST_COMPRESSION
  boost::shared_ptr<ifstream> ifs = boost::make_shared<ifstream>(filename, ios_base::in | ios_base::binary);
  if (!ifs->is_open()) throw StateFileOpeningException(filename);
  return ifs;
#else
  using namespace boost::iostreams;
  boost::shared_ptr<filtering_istream> in = boost::make_shared<filtering_istream>();
  file_source file(filename, std::ios_base::in | std::ios_base::binary);
  if (!file.is_open()) throw StateFileOpeningException(filename);
  if (isbz2(filename))
    in->push(bzip2_decompressor());
  in->push(file);
  return in;
#endif // DO_NOT_USE_BOOST_COMPRESSION
}

boost::shared_ptr<ostream> trajectory::openStateFileWriting(const std::string &filename, const ios_base::openmode mode)
{
#ifdef DO_NOT_USE_BOOST_COMPRESSION
  boost::shared_ptr<ofstream> ofs = boost::make_shared<ofstream>(filename, mode);
  if (!ofs->is_open()) throw StateFileOpeningException(filename);
  return ofs;
#else
  using namespace boost::iostreams;
  boost::shared_ptr<filtering_ostream> out = boost::make_shared<filtering_ostream>();
  file_sink file(filename, mode);
  if (!file.is_open()) throw StateFileOpeningException(filename);
  if (isbz2(filename)) {
    std::cerr << "Appending to compressed state files is not supported because of a boost bug." << std::endl;
    throw StateFileOpeningException(filename);
    // Appending to compressed state files does not work because of a bug in boost when handling multi-stream bz2 files.
    // The files can be written just fine, but cannot be read in afterwards. Hopefully this gets solved.
    // The patch in https://svn.boost.org/trac/boost/ticket/9749 actually does fix the problem.
    // See also: http://stackoverflow.com/q/32870991/1132850
    out->push(bzip2_compressor());
  }
  out->push(file);
  return out;
#endif // DO_NOT_USE_BOOST_COMPRESSION
}

bool trajectory::details::restoreState(Trajectory& traj, const string& trajectoryFileName, const string& stateFileName, const string& initialFileName)
{

  if (trajectoryFileName!="") {

    ifstream trajectoryFile(trajectoryFileName.c_str());

    if (trajectoryFile.is_open() && (trajectoryFile.peek(), !trajectoryFile.eof()) ) {
      boost::shared_ptr<istream> stateFile = openStateFileReading(stateFileName);
      while ( (stateFile->peek(), !stateFile->eof()) ) readViaSStream(traj,stateFile.get());
      return true;
    }

  }

  if (initialFileName!="") {
    boost::shared_ptr<istream> initialFile = openStateFileReading(initialFileName);
    while ( (initialFile->peek(), !initialFile->eof()) ) readViaSStream(traj,initialFile.get());
  }

  return false;
}


void trajectory::details::writeNextArchive(ostream* ofs, const ostringstream &oss)
{
  const string& buffer=oss.str();
  *ofs<<buffer.size(); ofs->write(&buffer[0],buffer.size());
}

void trajectory::details::readNextArchive(istream* ifs, istringstream &iss)
{
  string buffer;
  streamsize n; *ifs>>n; buffer.resize(n);
  ifs->read(&buffer[0],n);
  iss.str(buffer);
}


auto trajectory::readMeta(istream* ifs) -> SerializationMetadata
{
  istringstream iss(ios_base::binary);
  details::readNextArchive(ifs,iss);
  cpputils::iarchive archive(iss);
  SerializationMetadata meta;
  archive >> meta;
  return meta;
}

const std::string trajectory::SerializationMetadata::UNSPECIFIED = "Unspecified";
const std::string trajectory::SerializationMetadata::ARRAY_ONLY  = "ArrayOnly";


std::ostream& trajectory::details::DisplayAndAutostopHandler::display(std::ostream& os, int precision) const {return traj_.display(os,precision);}

auto trajectory::details::makeDisplayAndAutostopHandler(const Trajectory& traj, double autoStopEpsilon, unsigned autoStopRepetition) -> const DisplayAndAutostopHandler::Ptr
{
  class _ : public DisplayAndAutostopHandler
  {
  public:
    _(const Trajectory& traj, double autoStopEpsilon, unsigned autoStopRepetition)
      : DisplayAndAutostopHandler(traj), autoStopEpsilon_(autoStopEpsilon), autoStopRepetition_(autoStopRepetition) {}

  private:
    typedef DArray<1> Averages;

    std::ostream& display(std::ostream& os, int precision) const
    {
      istringstream iss;

      {
        ostringstream oss;
        DisplayAndAutostopHandler::display(oss,precision);

        os<<oss.str();

        iss.str(oss.str());
      }

      { double t, dtDid; iss>>t>>dtDid; } // eat time and dtDid

      if (!averages_.size()) { // This means that no display has yet occured: the number of averages must be determined
        size_t size(0);
        for (double dummy; (iss.peek(), !iss.eof()); (iss>>dummy, ++size) ) ;
        averages_.resize(size-1); averages_=0.;
      }
      else {
        Averages newItem(averages_.size());

        for (auto i=newItem.begin(); i!=newItem.end(); iss>>*i++); // Fill the new item

        if (queue_.size()==autoStopRepetition_) {
          if (max(abs(averages_-newItem)/(abs(averages_)+abs(newItem)))<autoStopEpsilon_) throw StoppingCriterionReachedException();
          averages_=averages_+(newItem-queue_.front())/autoStopRepetition_; // update the averages set for next step
          queue_.pop_front(); // remove obsolate first element of the queue
        }
        else averages_=(queue_.size()*averages_+newItem)/(queue_.size()+1); // build the initial averages set to compare against

        queue_.push_back(newItem); // place the new item to the back of the queue

      }

      return os;
    }

    const double autoStopEpsilon_;
    const unsigned autoStopRepetition_;

    mutable std::list<Averages> queue_;
    mutable Averages averages_;

  };

  if (autoStopRepetition) return std::make_shared<_>(traj,autoStopEpsilon,autoStopRepetition);
  else                    return std::make_shared<DisplayAndAutostopHandler>(traj);

}

