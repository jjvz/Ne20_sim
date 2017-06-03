#ifndef __myDetData__
#define __myDetData__ 1

#include <vector>
#include <stdio.h>
#include <iostream>
#include <TSystem.h> 
#include <TObject.h>


class myDetData : public TObject {
//class myDetData {
  public :
    myDetData();	// this is the class contructor
    virtual ~myDetData();	// this is the class destructor
  
  private :
    std::vector<double> DetEnergy;
    std::vector<double> DetTheta;
    std::vector<int> DetNr;
    std::vector<int> DetChan;
    std::vector<double> HitPos;
    std::vector<double> ChanPos;
  
  public :
    std::vector<double> Detx0;
    std::vector<double> Dety0;
    std::vector<double> Detx1;
    std::vector<double> Dety1;

    void SetEnergy(double ener);
    void SetTheta(double thet);
    void SetHitDet(int det);
    void SetHitChan(int ch);
    void SetChanPos(double x0);
	void SetHitPos(double x, double y);  

    double GetEnergy(int i);
    double GetTheta(int i);
    int GetHitDet(int det);
    int GetHitChan(int ch);
	std::vector<double> GetHitPos(int pos); 
	double GetChanPos(int i);	// only x0 positions of strips - 4 detectors x 16 strips each = 64 strip positions
    void ClearEvent();
    unsigned int SizeOfEvent(char c);

  ClassDef(myDetData,1)
    
};

ClassImp(myDetData);

#endif
