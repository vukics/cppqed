#include "Trajectory.tcc"

#include "ParsTrajectory.h"
#include "FormDouble.tcc"

#include <iostream>


using namespace std;


namespace trajectory {


void run(Trajectory& traj, const ParsRun& p)
{ 
  if (!p.Dt) {cerr<<"Nonzero Dt required!"<<endl; return;}
  if (p.NDt) run(traj,p.NDt,p.Dt,p.sdf,p.ofn,p.initialFileName,p.precision,p.displayInfo,p.firstStateDisplay);
  else       run(traj,p.T  ,p.Dt,p.sdf,p.ofn,p.initialFileName,p.precision,p.displayInfo,p.firstStateDisplay);
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


bool details::restoreState(Trajectory& traj, const string& trajectoryFileName, const string& stateFileName, const string& initialFileName)
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


namespace details {
  
void writeNextArchive(ofstream* ofs, const ostringstream &oss)
{
  const string& buffer=oss.str();
  *ofs<<buffer.size(); ofs->write(&buffer[0],buffer.size());
}

void readNextArchive(ifstream& ifs, istringstream &iss)
{
  string buffer;
  streamsize n; ifs>>n; buffer.resize(n);
  ifs.read(&buffer[0],n);
  iss.str(buffer);
}

} //details

SerializationMetadata readMeta(ifstream& ifs)
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

void writeViaSStream(const Trajectory& traj, ofstream* ofs)
{
  if (ofs && ofs->is_open()) {
    ostringstream oss(ios_base::binary);
    cpputils::oarchive stateArchive(oss);
    traj.writeState(stateArchive);
    details::writeNextArchive(ofs,oss);
  }
}


void readViaSStream(Trajectory& traj, ifstream& ifs)
{
  istringstream iss(ios_base::binary);
  details::readNextArchive(ifs,iss);
  cpputils::iarchive stateArchive(iss);
  traj.readState(stateArchive);
}

const std::string SerializationMetadata::UNSPECIFIED = "Unspecified";
const std::string SerializationMetadata::ARRAY_ONLY  = "ArrayOnly";

} // trajectory
