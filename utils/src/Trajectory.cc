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
      
      ifstream stateFile(stateFileName.c_str()/*, ios_base::binary*/);// stateFile.exceptions(ifstream::eofbit);
      
      if (!stateFile.is_open()) throw StateFileOpeningException(stateFileName);

      readViaSStream(traj,stateFile,true);
      
      return true;
    }
    
  }
  
  return false;
}


void writeViaSStream(const Trajectory& traj, ofstream* ofs)
{
  if (ofs && ofs->is_open()) {
    ostringstream oss(ios_base::binary);
    cpputils::oarchive stateArchive(oss);
    traj.writeState(stateArchive);
    const string& buffer=oss.str();
    *ofs<<buffer.size(); ofs->write(&buffer[0],buffer.size());
  }
}


void readViaSStream(Trajectory& traj, ifstream& ifs, bool fromEnd)
{ // scope of buffer
  string buffer;
  streamsize n; ifs>>n; buffer.resize(n);
  if (fromEnd) ifs.seekg(-n,ios_base::end); 
  ifs.read(&buffer[0],n);
  istringstream iss(buffer,ios_base::binary);
  cpputils::iarchive stateArchive(iss);
  traj.readState(stateArchive);
}


} // trajectory
