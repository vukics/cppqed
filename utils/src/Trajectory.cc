#include "impl/Trajectory.tcc"

#include "ParsTrajectory.h"
#include "impl/FormDouble.tcc"

#include <iostream>


using namespace std;


namespace trajectory {


void run(Trajectory& traj, const ParsRun& p)
{ 
  if (!p.Dt) {cerr<<"Nonzero Dt required!"<<endl; return;}
  if (p.NDt) run(traj,p.NDt,p.Dt,p.sdf,p.ofn,p.precision,p.displayInfo,p.firstStateDisplay);
  else       run(traj,p.T  ,p.Dt,p.sdf,p.ofn,p.precision,p.displayInfo,p.firstStateDisplay);
}


ostream& Trajectory::display(ostream& os, int precision) const
{
  const FormDouble fd(formdouble::positive(precision));
  return display_v( os<<fd(getTime())<<fd(getDtDid()) , precision)<<endl; // Note: endl flushes the buffer
}


ostream& Trajectory::displayKey(ostream& os) const
{
  size_t i=3;
  return displayKey_v( os<<"# Trajectory\n#  1. time\n#  2. dtDid\n" , i);
}

bool details::restoreState(Trajectory& traj, const string& trajectoryFileName, const string& stateFileName)
{
  if (trajectoryFileName!="") {

    ifstream trajectoryFile(trajectoryFileName.c_str());

    if (trajectoryFile.is_open() && (trajectoryFile.peek(), !trajectoryFile.eof()) ) {
      
      ifstream stateFile(stateFileName.c_str(), ios_base::binary); stateFile.exceptions(ifstream::eofbit);
      
      if (!stateFile.is_open()) throw StateFileOpeningException(stateFileName);
      
      cpputils::iarchive stateArchive(stateFile);
      
      try {
        while (true) /* for (int i=0; i<2; ++i) */ traj.readState(stateArchive);
      }
      catch (boost::archive::archive_exception /*ifstream::failure*/) {}
      
      return true;
    }
  }
  return false;
}

} // trajectory
