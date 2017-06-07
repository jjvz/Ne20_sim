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
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

/* root includes */
#include <TH1F.h>
#include <TH2F.h>     
#include <TTree.h>
#include <TFile.h>
#include <TRandom3.h>
#include <TROOT.h>
#include <TStyle.h>

#include "MCGAMMA.h"
#include "myMCFuncts.h"
#include "myDetData.h"

#define mulable   "mu.dat"
//#define LOG_FILE		//uncomment to write logfile for GAMMA TRANSPORT - this may slow down code !!!! 

// FOR MCROOT - GAMMA DECAY:
const float Res   = 60.;       	// resolution % (float) of Hagar
const float D     = 30.0;     	// distance from source (cm) (float)
const float H0    = 11.9;     	// radius of detector face (cm)(float)
const float maxdx = 35.6;       // depth of detector (cm) (float)
const int maxkan  = 18000;      // number of channels (integer)
const int NG      = 500;		// # gamma evrnts to print in logfile

float MuE[200], Mum[200][4];
float plength, E,theta_f,phi,rho_f,dx;
float *coord;
const int muentries=40;	// number of att. coeff. entries in mu.dat file

Int_t   numesc, numcomp, numpaar, numfoto;	
Float_t THETA, PHI, Edx; 	

// General variables for TTree
Int_t b_Ngam;
Int_t b_isgam=0;
Int_t b_esc;
Int_t b_comp1=0;
Int_t b_annihal=0;
Int_t b_fot1=0;
Int_t b_fot2=0;
Double_t b_compt;       
Double_t b_foto; 
Double_t b_paar; 
Double_t b_gamraw=0.; 
Double_t b_dist, b_gthet, b_gthet0;

TH1F *h_GamThet   = new TH1F("h_GamThet","Scatter angle",1000,0,180);
TH1F *h_GamDist   = new TH1F("h_GamDist","Free Paths",1000,0,maxdx);
TH1F *h_GamMca    = new TH1F("h_GamMca","Gamma rays; Energy [keV]",maxkan,0,maxkan);

#ifdef LOG_FILE
    std::ofstream logfile;
#endif

//------------- BEGIN MCROOT() -----------------------------------------
void MCRoot(Float_t Egam0, Int_t pindx, myDetData *gam, std::string loglable)
{
    Int_t i=0,  j, l, bini, NI;
    Float_t RHO, DX0, DY0, THETA_MAX;
   	TRandom3 *rgaus = new TRandom3(0);
    Char_t  buff[1000];

// **************** Read Mu file (att. coefficients) **************************************
// Data file with gamma ray mass att. coeff's (in cm2/g): 0=Compt.; 1=Photo; 2=Pair; 3=Total 
// as function of gamma ray energy (MuE) in keV:    
#ifndef MUFILE_H	// open and read only once!
#define MUFILE_H
	FILE *mufile;
    mufile = fopen(mulable, "r");
    if(mufile){
    	fgets(buff, 80, mufile);
    	fgets(buff, 80, mufile);
    	for(j=0; j<muentries; j++) {
        	fscanf(mufile,"%e\t%e\t%e\t%e\t%e",&MuE[j],&Mum[j][0],&Mum[j][1],&Mum[j][2],&Mum[j][3]);
    	}
    	fclose (mufile);
    }
    else {
    	std::cout<< "\nERROR: "<< mulable <<" file does not exist!!!!"<<std::endl;
    	return;
    }
#endif
// ****************************************************************************************

#ifdef LOG_FILE
	if(!logfile.is_open()) {
		logfile.open (loglable.c_str(), std::ios::out | std::ios::trunc);
	}
#endif

//	b_Ngam=0;
    b_isgam=1;
    THETA_MAX = atan(H0/D);	        //max angle (deg) for detector face area.
                                                            //anders vat die program eeue!!
    numesc=0;
    numcomp=0;
    numpaar=0;
    numfoto=0;
        
// Begin event:>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    b_gamraw=Egam0;	//stores initial gamma energy
    Edx=0.;
  
	THETA = THETA_MAX*rgaus->Rndm();    //gamma incident angle, theta (rad)
	PHI = PI/180.*360.*rgaus->Rndm();                  //gamma incident angle, phi (rad)
	RHO = D*tan(THETA);                              //gamma incident radial position, rho
	plength = WHERE(Egam0,3);      //path length inside detector
	dx = plength*cos(THETA);
	b_dist = plength;
	h_GamDist->Fill(plength);
	b_gthet0 = THETA*180./PI;	// initial gamma incident angle on HAGAR
	phi = PHI;
	rho_f = RHO+(dx*tan(THETA));                              //new radial distance after plength

	if(dx>maxdx || rho_f>H0) {                      //if outside detector
		++numesc;
		b_esc=1;
	}
	if(dx<=maxdx || rho_f<=H0) INSIDE=true;
	else INSIDE=false;

#ifdef LOG_FILE
	logfile<<"MCRoot called:---------------------------------------\n";
	logfile<<"energy = "<<Egam0<<" keV ; b_Ngam = "<<b_Ngam<<"\n";
#endif
		
// ****************************** INTERACTION************************************
	for(E=Egam0,theta_f=THETA; ((INSIDE && (E!=0.)) || (b_annihal==1) || (b_annihal==2)); i++) { //loop until out of detector or fully absorbed (E=0.)

#ifdef LOG_FILE
		logfile<<"\tNow in MCRoot(): "<<i<<" ---------------------------------------\n";
		logfile<<"\tenergy = "<<E<<" keV ; b_Ngam = "<<b_Ngam<<"\n";
		logfile<<"\tINSIDE? "<<INSIDE<<" \n";
#endif

		INTERACT();         //At point B
		
#ifdef LOG_FILE
  		if(b_Ngam<NG) 
			logfile<<" foto abs.: E = "<<E<<" ("<<b_annihal<<") ; numfoto = "<<numfoto<<"\n";
#endif

		h_GamThet->Fill(theta_f*180./PI);
	}
	if(Edx>0.) {
		gam->SetEnergy(Edx);	
   
		h_GamMca->Fill(Edx);
	}

#ifdef LOG_FILE
    if(b_Ngam<NG) {
//    	fprintf(logfile, "Event: %d , gam evt: %d E* = %.3f MeV, gamma %.0f detected with Energy = %.0f keV \n", b_N0, b_Ngam, b_energSA, b_gamraw,Edx);
		logfile<<" Gamma evt: "<<b_Ngam<<", gamma "<<b_gamraw<<" keV detected with Energy = "<<Edx<<" keV\n";
		logfile<<" _________________________________________________________________________________ \n\n";
    }
#endif

	Edx=0.;
//	ZeroTTreeVariablesMC();
    b_Ngam++;

//#ifdef LOG_FILE	// Can't close file before whole program is done! Else file will be cleared.
//	logfile.close();
//#endif
}
//-------------END MCROOT()-------------------------------------------
//---------------------------INTERACT-----------------------------
void INTERACT()
{
    float R, p_abs, p_foto, p_compt, p_paar, dEpaar, paarrho=0, paarx=0;       
   	TRandom3 *rgaus = new TRandom3(0);

    p_compt=1-exp(-MU(E,0,MuE)*plength);
    p_foto=1-exp(-MU(E,1,MuE)*plength);
    if(E>=1022.) p_paar=1-exp(-MU(E,2,MuE)*plength);
    else p_paar=0.;

    p_abs=p_foto+p_compt+p_paar;		     //absorbsie

    R=rgaus->Rndm();
//------------------- FOTO ABSORBPTION -------------------------------
    if(R <= (p_foto/p_abs)) {   // Foto electriese absorbsie
        ++numfoto;
        if(b_comp1==0) b_fot1=1;
        else b_fot2=1;
#ifdef LOG_FILE
  		if(b_Ngam<NG) 
			logfile<<" foto abs.: E = "<<E<<" ("<<b_annihal<<") ; numfoto = "<<numfoto<<"\n";
#endif
        Edx=Edx+convolv(E);   // sum energy within plength - convoluted, may also already contain... 
        b_foto=Edx;           // ... previous compton Edx !!!!                     
        if(b_annihal==1) {	// if gamma 1 just interacted, switch to gamma 2
        	E=511.;
        	dx=paarx;		// back where the annihilation took place
        	rho_f=paarrho;	// back where the annihilation took place
            BOOST(theta_f+PI, phi+PI);	// second 511 gamma goes in opposite diretion
        	if(INSIDE) b_annihal=2; 	// gamma 2 must still interact
        	else {
        		b_annihal=0;	// gamma 2 outside, done with
        		E=0.;
        	}
        }
        else {
        	b_annihal=0;	// gamma 2 done with
        	E=0.;			// gamma 2 just interacted, done with paar, Event stops here - gamma fully absorbed!
        }
    }
//------------------- COMPTON SCATTERING -------------------------------
    else if(R <= p_compt/(p_compt+p_paar)) {
        ++numcomp;
        if(b_comp1==0) b_comp1=1;
        E = COMPTON();            // secondary gamma; new dx, theta_f, rho_f and phi
        BOOST(theta_f, phi);
    }
//------------------- PAIR PRODUCTION -------------------------------
    else {  // Pair Production: we assume (for now) both annihilations gammas escape, so initial 
            // gamma dissapears, using 1022 keV for e+ and e- creation, energy balance (E-1022) 
            // shared by e+ and e- deposited in detector as e+ and e- is slowed down. 
            // e+ then annihilates, creating 2x 511 keV gammas, which both mostly escapes small detectors 
            // - thus peak at E-1022. If only one escapes, extra 511 deposited in detector - peak at 
            // E-1022+511. If both absorbed - full energy deposited  - full energy peak.
        ++numpaar;	// new energy of (e-,e+) pair to be absorbed... or not. The e- takes ~half energy and is usually immediately absorbed and constitutes Edx. 
#ifdef LOG_FILE
		logfile<<" pair: E = "<<E<<" keV; numpair = "<<numpaar<<"\n";
#endif

        Edx=Edx+convolv((E-1022.)/2.);	 	// electron absorbed taking half the remainder (Egam0-1022):
        Edx=Edx+convolv((E-1022.)/2.);		// positron kinetic energy absorbed as it is slowed down before being annihilated:
       	E = 511.;		// gamma 1: the e+ annihilates, creating two opposite 511 keV gammas.
        b_paar=Edx;
        paarx=dx;
        paarrho=rho_f;
        BOOST(PI/180.*360.*rgaus->Rndm(), PI/180.*360.*rgaus->Rndm());	// 511 goes anywhere.
        if(INSIDE) b_annihal=1;		// gamma 1 must still interact
        else {				// switch to gamma 2:
        	dx=paarx;		// back where the annihilation took place
        	rho_f=paarrho;	// back where the annihilation took place
            BOOST(theta_f+PI, phi+PI);	// second 511 gamma goes in opposite diretion		
	        if(INSIDE) b_annihal=2;		// gamma 2 must still interact
			else b_annihal=0; 			// gamma 2 lost, exit	
		}
    }
}
//---------------------------COMPTON-------------------------------
//This function calculates compton scattering secondary gamma energy
float COMPTON()
{
    float Egam, Ee;
   	TRandom3 *rgaus = new TRandom3(0);

    theta_f=gooi_dice_kn();	                //Klein-Nishina angular cross-section
    b_gthet=theta_f*180./PI;              //!! this overrides previous compt. scat. theta!!!
    phi=PI/180.*360.*rgaus->Rndm();                //isotropic
    Egam=E/(1.+E*(1.-cos(theta_f))/511.);       //scattered/secondary gamma energy
    Ee=E-Egam;                                  //energy deposited to electrons

    Edx=Edx+convolv(Ee);
    b_compt=Edx;    // CAREFULL!! - I'm overwriting previous comptons during multiple compt. scat.
#ifdef LOG_FILE
    if(b_Ngam<NG) {
		logfile<<" compt. scatt. of "<<E<<" ("<<b_annihal<<") at thet = "<<theta_f*180./PI<<" deg, phi = "<<phi*180./PI<<" deg; numcomp = "<<numcomp<<"\n";
		logfile<<" \tEgam out = "<<Egam<<" keV ("<<b_annihal<<"), "<<Ee<<" keV deposited\n";
    }
#endif
    return Egam;   // remember Egam is UNconvoluted (100% resolution)!!
}
//---------------------------COMPTON-------------------------------
//--------------------------- ZeroTTreeVariablesMC() -----------------------------
void ZeroTTreeVariablesMC()
{
   b_dist=0.;
   b_gthet=0.;
   b_gthet0=0.;     
   b_foto=-1.0,
   b_compt=-1.0;
   b_esc=0;
   b_comp1=0;
   b_fot1=0;    
   b_fot2=0;    
   b_paar=-1.0;
   b_isgam=1;
   b_gamraw=0.;
   b_annihal=0;
}

//---------------------------KLEIN-NISHINA ANGLE-------------------
float gooi_dice_kn()        //hoek in radians
{
    float thet_i, Rkn, integral, kn, prob, facto, Ykn;
    int ikn, succ2;
   	TRandom3 *rgaus = new TRandom3(0);

    for(integral=0., ikn=0; ikn<=1800; ikn++) {
	thet_i=ikn*PI/1800.;
    facto=1.+E*(1.-cos(thet_i))/511.;
	kn=pow((1./facto),2)*(facto+(1./facto)-pow(sin(thet_i),2));
	integral=integral+kn;
    }
    for(succ2=0; succ2!=1; ) {
        thet_i = PI*rgaus->Rndm();              //random thet between 0 and pi radians     
        facto=1.+E*(1.-cos(thet_i))/511.;
        Ykn=pow((1./facto),2)*(facto+(1./facto)-pow(sin(thet_i),2))/integral;
        Rkn=rgaus->Rndm();
        if(Rkn <= Ykn) succ2=1;
        else {};
    }
    return thet_i;
}
//---------------------------KLEIN-NISHINA ANGLE-------------------
//---------------------------CONVOLVE------------------------------
//This function spreads the energy using a gaussian random distribution 
// according to the input selected resolution.
//It returns an energy value, distributed accordingly. 

float convolv(float imu)
{
    int succ;
    float Xg, Yg, rg;     
   	TRandom3 *rgaus = new TRandom3(0);
    
    for(succ=0; succ!=1; ) {
        Xg=(imu-3.322*Res)+0.006644*Res*(1000.*rgaus->Rndm());                   //full-width-at-1%-max
//        Yg=(1./(1.0645*Res))*exp(-pow((Xg-imu)/(0.6006*Res),2.0));       //gauss function normalised to 1
        Yg=1.*exp(-pow((Xg-imu)/(0.6006*Res),2.0));       //gauss function normalised to 1
        rg=rgaus->Rndm();
        if(rg <= Yg) succ=1;
        else {};
    }
	return Xg;
}
//---------------------------CONVOLVE------------------------------
//--------------------------- MU(E, type, MuE); -----------------------------
// Calculates the linear attenuation coefficient at given gamma ray energy
// from mass att. coeff. [in cm^2/g] data file ("mu.dat"):
float MU(float iEmu, int type, float *MuE)
{
  float MM,CC,DENS=3.67;               //density of NaI=3.67 g/cm3
  int imu;
  for(imu=0;imu<muentries;imu++){
     if(MuE[imu]<iEmu) continue;
     else if(MuE[imu]==iEmu) return Mum[imu][type]*DENS;        //interpolated linear attenuation coeff.
     else if(MuE[imu]>iEmu){
        MM=(Mum[imu][type]-Mum[imu-1][type])/(MuE[imu]-MuE[imu-1]);
        CC=Mum[imu][type]-MM*MuE[imu];
        return (MM*iEmu+CC)*DENS;        //interpolated linear attenuation coeff.
     }
     else { 
        std::cout<<"problem in MU().... no coefficient found!!"<<std::endl;
        return 1000000.;
     }
  }
}


//--------------------------- WHERE(E,mu) -----------------------------
// Determines the dx (path length) into the medium where the interaction took place:
float WHERE(float iEpos, int type)
{  
   	TRandom3 *rgaus = new TRandom3(0);
	float YY,XX=0.;
	for(YY=0.;YY==0.;){
		YY=rgaus->Rndm();
		XX=-log(1-YY)/MU(iEpos,type,MuE);    //XX is path length in cm, if density is in g/cm3
	}
	return XX;
}
//--------------------------- WHERE(E,mu) -----------------------------

//--------------------------- BOOST(path,theta,phi) -----------------------------
// Transports gamma to new distance and angle:
void BOOST(float bthet, float bphi)
{  
	plength = WHERE(E,3);          	// secondary gam. now travels to point C
   	theta_f = bthet;              	//gamma incident angle, isotropic (rad)
   	phi = bphi;                   	//gamma incident angle, isotropic (rad)
    coord = ROTAT(plength*cos(theta_f),tan(theta_f),theta_f);
    rho_f = rho_f + coord[1];
    dx = dx + coord[0];				// dx is projection of plength in z-direction
  	if(dx<=maxdx || rho_f<=H0) INSIDE=true;
	else INSIDE=false;

}
//--------------------------- BOOST(path,theta,phi) -----------------------------

//--------------------------- ROTATE -----------------------------
// rotates the (x,y) coordinates in 2D through angle 'angl':

float *ROTAT(float xx, float yy, float angl)
{
    static float newcoord[2];
    newcoord[0] = xx*cos(angl)- yy*sin(angl);
    newcoord[1] = xx*sin(angl)+ yy*cos(angl);
    return newcoord;
}
//--------------------------- ROTATE -----------------------------
