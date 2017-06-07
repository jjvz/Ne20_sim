#include <iostream>
#include "myDetData.h"

// Member functions definitions including constructor
myDetData::myDetData(void) {
//   std::cout << "\nObject myDetData is being created." << std::endl;
}

myDetData::~myDetData(void) {
//   std::cout << "Object myDetData is being deleted.\n" << std::endl;
}

void myDetData::ClearEvent()
{
  DetEnergy.clear();
  DetTheta.clear();
  DetNr.clear();
  DetChan.clear();
  HitPos.clear();
}

void   myDetData::SetEnergy(double ener) { DetEnergy.push_back(ener); }
void   myDetData::SetTheta(double thet)  { DetTheta.push_back(thet); }
void   myDetData::SetHitDet(int det) 	 { DetNr.push_back(det); }
void   myDetData::SetHitChan(int ch) 	 { DetChan.push_back(ch); }
void   myDetData::SetChanPos(double x0)  { ChanPos.push_back(x0); }	// only x0 positions of strips - 4 detectors x 16 strips each = 64 strip positions

void   myDetData::SetHitPos(double x, double y)  
{ 
  HitPos.push_back(x); 
  HitPos.push_back(y); 
}

double myDetData::GetEnergy(int i)       { return DetEnergy.at(i); }
double myDetData::GetTheta(int i)        { return DetTheta.at(i); }
int    myDetData::GetHitDet(int det) 	 { return DetNr.at(det); }
int    myDetData::GetHitChan(int ch) 	 { return DetChan.at(ch); }
double myDetData::GetChanPos(int i)  	 { return ChanPos.at(i); }	// only x0 positions of strips - 4 detectors x 16 strips each = 64 strip positions

std::vector<double> myDetData::GetHitPos(int pos) 
{ 
	std::vector<double> xy;
	xy.clear();
		
	xy.push_back( HitPos.at(pos) );
	xy.push_back( HitPos.at(pos+1) );	
    return xy;
}

unsigned int myDetData::SizeOfEvent(char c)
{
  if(c=='p') { return HitPos.size(); }
  else if(c=='e') { return DetEnergy.size(); }
  else if(c=='s') { return DetChan.size(); }
  else if(c=='c') { return ChanPos.size(); }
  else if(c=='d') { return DetNr.size(); }
  else { return DetEnergy.size(); }
}

