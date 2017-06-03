//****************************************************************************************************
//****************************************************************************************************
//Monte Carlo type program for incident gammas of specific energy and relative
//intencity.
//The program receives a list of gammas and rel. intensities and generates
//(a)  a Monte Carlo simulation,
//(b)  a Gaussian "convolution" with a user defined resolution
//
//Usage: in ROOT:>> .x MCRoot.C)
//
//Written by JJ van Zyl, 30-04-2004, addapted to ROOT: Sept 2013, last updated: July 2015
//******************************************************************************
//------------- MCROOT() -----------------------------------------
#include "myMCFuncts.h"
#include "myDetData.h"

void MCRoot(Float_t Egam0, Int_t pindx, myDetData *gam, std::string loglable);
void INTERACT();
float COMPTON();
float MU(float iEmu, int type, float *MuE);
void ZeroTTreeVariablesMC();
float gooi_dice_kn(); 
float convolv(float imu);
float WHERE(float iEpos, int type);
void BOOST(float bthet, float bphi);
float *ROTAT(float xx, float yy, float angl);	//usage: newcoord = rotat(old_x, old_y, angle), then newx=newcoord[0] and newy = newcoord[1]


