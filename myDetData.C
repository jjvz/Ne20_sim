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
  if(c=='p') {
    	unsigned int result = HitPos.size();
  	return result;
  }
  else if(c=='e') { 
  	unsigned int result = DetEnergy.size();
  	return result;
  }
  else if(c=='s') { 
  	unsigned int result = DetChan.size();
  	return result;
  }
  else if(c=='c') { 
  	unsigned int result = ChanPos.size();
  	return result;
  }
  else if(c=='d') { 
  	unsigned int result = DetNr.size();
  	return result;
  }
  else { 
  	unsigned int result = DetEnergy.size();
  	return result;
  }
}

