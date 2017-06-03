//=========Display Scattering Geometry for 20Ne gas cell target setup
//=========  JJvZ, Jul 7, 2014

// Usage:
// To Compile with g++ in Ubuntu:
// 1) In order to add the new TClass "myDetData" to Root's Dictionary of classes, use "rootint" as follows:
// ~/codes$ rootcint -f ne20_Dict.C -c myDetData.h my_si_LinkDef.h
// The file "my_det_LinkDef.h" is shown below and must be in the same directory. The files "my_det_Dict.C" and "my_det_Dict.h" are created by rootint
// 2) Compile the files into a library:
// ~/codes$ g++ Ne20_gascell.C BREAKUP.C E_OUT.C myMCFuncts.C myDetData.C ne20_Dict.C `root-config --libs --cflags` -o Ne20_gascell
// 3) Run the executable:
// ~/codes$ ./Ne20_gascell
//
// This following bash script will do all the above:
// source ./run_ne20.sh 

// Running:
// ./Ne20_gascell

// To compile in ROOT with CINT:
// root [0] .L BREAKUP.C
// root [1] .L E_OUT.C
// root [2] .L myMCFuncts.C
// root [3] .L myDetData.C++
// root []  .L MCGAMMA.C
// root [4] .L Ne20_gascell.C
// root [5] ne20_gascell()

// .cd E:\Navorsing\iThemba labs\20Ne\Sims\Ne02_sim-master

// To view events/branches/leafs, run:
// root [3] DATA->Scan("evntno:GamData->DetEnergy[0]");	//e.g. for leaf evntno and GamData->DetEnergy[0]
// root [4] DATA->Scan("evntno:GamData->DetEnergy[2]:SiData->DetEnergy[2]", "evntno==456");	// or just a specific entry:	
// root [4] DATA->Scan("evntno:GamData->DetEnergy[2]:SiData->DetEnergy[2]", "evntno==456", "colsize=13 precision=3 col=13::15.10");

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include <TH1F.h>
#include <TH2F.h>     
#include <TTree.h>
#include <TFile.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h> 
#include <TObject.h>
#include <TCanvas.h>
#include <TWbox.h>
#include <TEllipse.h>
#include <TLine.h>
#include <TBox.h>
#include <TApplication.h>

#include "E_OUT.h"
#include "BREAKUP.h"

// *************** Find relative energy between breakup particles************************
float E_REL(float ec, float ed, float tc, float td);

#define gamlable1 "ne20-gammas.dat"
#define gamlable2 "o16-gammas.dat"
#define gamlable3 "c12-gammas.dat"

FILE *gamfile1; 
FILE *gamfile2; 
FILE *gamfile3; 

const int NN=50;

// INPUT:
const double E0=200.0 ;        // [in MeV]
const double THscat=0.0 ;          // outgoing angle of ejectile [in deg]
const double mp=4.00260 ;      // mass of projectile [in amu]  // 4.00260 for 4He; 1.00782 for H; 15.99491 for...
const double p_Ne20=0.65;	// probability for 20Ne -> a + 16O
const double p_Ne20Be=0.20;	// probability for 20Ne -> 8Be + 12C
const double p_C12=0.15;	// probability for 12C -> a + 8Be

const float XFWHM = 0.08;		// fwhm (in MeV) of focal plane x-axis

const int numlines1=69;		// # states for Ne20 datfile
const int numlines2=14;	// # states for O16 datfile
const int numlines3=4;	// # states for C12 datfile   

/*------------ Constants----------------------------------------------*/
const double PI=3.14159265359;
const float amu = 931.5016;  // 1 antomic mass unit [in MeV]
const int   q_charge=2;
const float spec_rad=8.14;
const float regit=276.74;
const float m0=2809.4;
const float Epeak=200.0;         // [in MeV]


/*-- GLOBAL DECLARATIONS --------------------------------------------*/

Float_t I[30],E_rel=-1.;
Int_t   counterj;		
Char_t  flagt, IN[28], lable[28];
Float_t m_b, m_B, E_tot, p_A, p_c, p_d, p_b, p_B, E_b=0., E_B, rawEb=0, rawEB=0, rawEc=0, rawEd=0, E_rel=-1., E_c, E_d;

Int_t cnt=0, ti, flag=0;
float *breakdat;


// ************************************************************************************

int breakup_test()
{
    Int_t miss, i, j, p, Nj, numscat=0, nummis=0, ind=0, indx;
    Float_t rg, x1, y1, x2b, y2b, y_si1, y_si2, z1, z2, ran, ran2, ran3, ran4;
    Float_t rgE, rgGam, rgE2, x2B, y2B, Estar_A=0, Estar_B, S_b, m_A, m_c=0, m_d=0, S_2a;
    Float_t energ_b, energ_B, E_A, thet_b, thet_c=0, thet_d=0, thet_B, phi, PHI_MAX, Xposition, thet_max, THET_MAX;
	Float_t thet_in=0;
    Char_t  buff[1000];
	Int_t b_Ne20, b_Ne20Be, b_C12;
	Int_t evt;

	TRandom3 *rgaus = new TRandom3(0);
	TRandom3 r(0);

// **************************************************************************************
// ********** Start event ***************************************************************
// **************************************************************************************
 
    for(Nj=0; Nj<NN; Nj++) {

		b_Ne20=0;
		b_Ne20Be=0;
		b_C12=0;
    	thet_in=0;

	    std::cout<<"\n ***************** EVENT: "<<Nj+1 <<" ***********************"<<std::endl;
	    std::cout<<" ------------------------------------------"<<std::endl;

    	rgE = rgaus->Rndm();;    
// ********************************* Ne-20 state - a + 16O ******************************
    	if(rgE<=p_Ne20) {
        	Estar_A = 13.58;
        	S_b = -4.730;               // separation energy for cluster from tgt
        	m_A = 19.99244;             // breakup nucleus: 19.99244 for 20Ne
        	m_b = 4.00260;              // light breakup ion (alpha)
        	m_B = 15.99491;             // heavy breakup ion (16O)

        	E_A = E0-Estar_A-E_OUT(E0, mp, m_A, mp, m_A, 0, Estar_A, THscat);
        	Int_t b_Ne20 = 1;
			Estar_B = 6.049;           
		    std::cout<<" ***************** Ne-20 -> a + 16O :  ***********************"<<std::endl;
    	}
    // ********************************* Ne-20 state - 8Be + 12C *****************************************************
    	else if(rgE>p_Ne20 && rgE<=(p_Ne20+p_Ne20Be)){

		    Estar_A = 18.8;
		    S_b = -11.984;               // separation energy for cluster from tgt
		    m_A = 19.99244;             // breakup nucleus: 19.99244 for 20Ne
		    m_b = 8.00531;              // light breakup ion (8Be)
		    m_B = 12.00000;             // heavy breakup ion (12C)
		    m_c = 4.00260;              // breakup alphas (8Be->a+a)
		    m_d = 4.00260;              // breakup alphas (8Be->a+a)
			S_2a = 0.092;       // separation energy for cluster from tgt

		    E_A = E0-Estar_A-E_OUT(E0, mp, m_A, mp, m_A, 0, Estar_A, THscat);
		    Int_t b_Ne20Be = 1;
		    Estar_B = 4.43891;        
		    std::cout<<" ***************** Ne-20 -> 8Be + 12C :  ***********************"<<std::endl;
		} 

// ********************************* C-12 state - a + 8Be ***********************************************************
		else {
			Estar_A=16.45;
		    S_b = -7.367;               // separation energy for cluster from tgt
		    m_A = 12.00000;             // breakup nucleus: 19.99244 for 20Ne
		    m_b = 4.00260;              // light breakup ion (alpha)
		    m_B = 8.00531;             // heavy breakup ion (8Be)
		    m_c = 4.00260;              // breakup alphas (8Be->a+a)
		    m_d = 4.00260;              // breakup alphas (8Be->a+a)
			S_2a = 0.092;       // separation energy for cluster from tgt

		    E_A = E0-Estar_A-E_OUT(E0, mp, m_A, mp, m_A, 0, Estar_A, THscat);
		    Int_t b_C12 = 1;
	  
		    Estar_B = 0.;
		    std::cout<<" ***************** C-12 -> a + 8Be :  ***********************"<<std::endl;
		}        
// ****************************************************************************************************************
		
		m_A = m_A*amu;
		m_b = m_b*amu;
		m_B = m_B*amu;
		m_c = m_c*amu;
		m_d = m_d*amu;
// ****************** Scattering (breakup) angle of light ion: **************************
		p_A = sqrt(2.*m_A*E_A);

	    std::cout<<"E_A = "<<E_A<<" ,m_A = "<<m_A<<", p_A = "<<p_A<<std::endl;

		breakdat = breakup(p_A, m_A, m_b, m_B, E_A, (Estar_A-Estar_B), S_b);  // Call this action ONCE to calc breakup kinematics!
		if(breakdat[7]==1) continue;

		p_b = breakdat[0];
		p_B = breakdat[1];
		E_b = breakdat[2];
		E_B = breakdat[3];
		thet_b = breakdat[4]; // Breakup angle of heavy (16O) ion (CLOCKwise from +x-axis) [radians]
		thet_B = breakdat[5]; // Breakup angle of heavy (16O) ion (CLOCKwise from +x-axis) [radians]
		thet_max = breakdat[6]; // Breakup angle of heavy (16O) ion (CLOCKwise from +x-axis) [radians]

	    std::cout<<"p_b = "<<p_b<<", p_B = "<<p_B<<std::endl;
	    std::cout<<"thet_b = "<<thet_b*180./PI<<", thet_B = "<<thet_B*180./PI<<std::endl;
	    std::cout<<"thet_max = "<<thet_max*180./PI<<std::endl;

//******************************* breakup of 8Be -> a + a *******************************************
		if(b_Ne20Be==1) {
			breakdat = breakup(p_b, m_b, m_c, m_d, E_b, 0, S_2a);  // Call this action ONCE to calc breakup kinematics!
			if(breakdat[7]==1) continue;

			p_c = breakdat[0];
			p_d = breakdat[1];
			E_c = breakdat[2];
			E_d = breakdat[3];
			thet_c = breakdat[4]; // Breakup angle of heavy ion (CLOCKwise from +x-axis) [radians]
			thet_d = breakdat[5]; // Breakup angle of heavy ion (CLOCKwise from +x-axis) [radians]
			thet_max = breakdat[6]; // Breakup angle of heavy (16O) ion (CLOCKwise from +x-axis) [radians]

			E_rel = E_REL(E_c, E_d, thet_c, thet_d); 

			std::cout<<" ***************** 8Be -> a + a :  ***********************"<<std::endl;
		    std::cout<<"E_0 = "<<E_b<<" ,m_0 = "<<m_b<<", p_b = "<<p_b<<std::endl;
			std::cout<<"p_c = "<<p_c<<", p_d = "<<p_d<<std::endl;
			std::cout<<"thet_c = "<<thet_c*180./PI<<", thet_d = "<<thet_d*180./PI<<std::endl;
			std::cout<<"E_rel = "<<E_rel<<std::endl;
			std::cout<<"thet_max = "<<thet_max*180./PI<<std::endl;
		}

//******************************* breakup of 8Be -> a + a *******************************************
		else if(b_C12==1) {
			breakdat = breakup(p_B, m_B, m_c, m_d, E_B, 0, S_2a);  // Call this action ONCE to calc breakup kinematics!
			if(breakdat[7]==1) continue;

			p_c = breakdat[0];
			p_d = breakdat[1];
			E_c = breakdat[2];
			E_d = breakdat[3];
			thet_c = breakdat[4]; // Breakup angle of heavy ion (CLOCKwise from +x-axis) [radians]
			thet_d = breakdat[5]; // Breakup angle of heavy ion (CLOCKwise from +x-axis) [radians]
			thet_max = breakdat[6]; // Breakup angle of heavy (16O) ion (CLOCKwise from +x-axis) [radians]

			E_rel = E_REL(E_c, E_d, thet_c, thet_d); 

			std::cout<<" ********************** from 12Cjj: 8Be -> a + a :  *****************"<<std::endl;
		    std::cout<<"E_0 = "<<E_B<<" ,m_0 = "<<m_B<<", p_b = "<<p_B<<std::endl;
			std::cout<<"p_c = "<<p_c<<", p_d = "<<p_d<<std::endl;
			std::cout<<"thet_c = "<<thet_c*180./PI<<", thet_d = "<<thet_d*180./PI<<std::endl;
			std::cout<<"E_rel = "<<E_rel<<std::endl;
			std::cout<<"thet_max = "<<thet_max*180./PI<<std::endl;
		}

	}  // end of event
    std::cout<<std::endl;

	return 0;
}

/*
 ***************** C-12 -> a + 8Be :  ***********************
E_A = 0.128643 ,m_A = 11178, p_A = 53.6279
p_b = 207.698, p_B = 226.061
thet_b = 103.206, thet_B = 63.4404
thet_max = 360
 ********************** from 12C: 8Be -> a + a :  *****************
E_0 = 3.42658 ,m_0 = 7456.96, p_b = 226.061
p_c = 96.1411, p_d = 130.363
thet_c = 4.17171, thet_d = 3.07535
E_rel = 0.0919533
thet_max = 9.42833

*/

// *************** Find relative energy between breakup particles************************
float E_REL(float ec, float ed, float tc, float td)
{
	return (0.5*(ec+ed) - sqrt(ec*ed)*cos(tc + td));
}


