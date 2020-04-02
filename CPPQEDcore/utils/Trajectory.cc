// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "ArrayTraits.h"
#include "BlitzArray.h"
#include "Trajectory.tcc"

#include "ParsTrajectory.h"
#include "FormDouble.tcc"

#include <iostream>
#include <deque>


using namespace std;


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


bool trajectory::details::restoreState(Trajectory& traj, const string& trajectoryFileName, const string& stateFileName, const string& initialFileName)
{
  
  if (trajectoryFileName!="") {

    ifstream trajectoryFile(trajectoryFileName.c_str());

    if (trajectoryFile.is_open() && (trajectoryFile.peek(), !trajectoryFile.eof()) ) {
      
      ifstream stateFile(stateFileName.c_str()/*, ios_base::binary*/);// stateFile.exceptions(ifstream::eofbit);
      
      if (!stateFile.is_open()) throw StateFileOpeningException(stateFileName);

      // readViaSStream(traj,stateFile,true);
      while ( (stateFile.peek(), !stateFile.eof()) ) readViaSStream(traj,stateFile);
      
      return true;
    }
    
  }
  
  if (initialFileName!="") {
    ifstream initialFile(initialFileName.c_str(), ios_base::binary);
    if (!initialFile.is_open()) throw StateFileOpeningException(initialFileName);
    while ( (initialFile.peek(), !initialFile.eof()) ) readViaSStream(traj,initialFile);
  }
  
  return false;
}


void trajectory::details::writeNextArchive(ofstream* ofs, const ostringstream &oss)
{
  const string& buffer=oss.str();
  *ofs<<buffer.size(); ofs->write(&buffer[0],buffer.size());
}

void trajectory::details::readNextArchive(ifstream& ifs, istringstream &iss)
{
  string buffer;
  streamsize n; ifs>>n; buffer.resize(n);
  ifs.read(&buffer[0],n);
  iss.str(buffer);
}


auto trajectory::readMeta(ifstream& ifs) -> SerializationMetadata
{
  istringstream iss(ios_base::binary);
  streampos pos = ifs.tellg();
  details::readNextArchive(ifs,iss);
  ifs.seekg(pos);
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

      { double t; iss>>t; } // eat time

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

    mutable std::deque<Averages> queue_;
    mutable Averages averages_;

  };

  if (autoStopRepetition) return std::make_shared<_>(traj,autoStopEpsilon,autoStopRepetition);
  else                    return std::make_shared<DisplayAndAutostopHandler>(traj);

}

