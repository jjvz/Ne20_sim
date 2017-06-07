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
#include <TArrow.h>
#include <TLatex.h>
#include <TApplication.h>

#include "myMCFuncts.h"
#include "myDetData.h"
#include "E_OUT.h"
#include "BREAKUP.h"
#include "MCGAMMA.h"

#define gamlable1 "ne20-gammas.dat"
#define gamlable2 "o16-gammas.dat"
#define gamlable3 "c12-gammas.dat"
#define diaglable "ne20-diag.dat"
#define loglable  "ne20-gammas-log.dat"
#define plotlable "gascell_draw.dat"

#define PRNT_DIAG		//uncomment to print scatter diagram 
//#define DIAG_FILE		//uncomment to write diagnosis file - this may slow down code !!!! 

FILE *gamfile1; 
FILE *gamfile2; 
FILE *gamfile3; 

const int nSi=4;
const int n_si_det=4;
const int n_si_ch=16;
const int n_gam_det=1;
const int n_gam_ch=7;
const int NN   = 200000;	// number of event iterations 
const int PRNT = 2000;	// number of scatter events to plot in diagram
// INPUT:
const double E0=200.0 ;        // [in MeV]
const double THscat=0.0 ;          // outgoing angle of ejectile [in deg]
const double mp=4.00260 ;      // mass of projectile [in amu]  // 4.00260 for 4He; 1.00782 for H; 15.99491 for...
const double p_Ne20=0.75;	// probability for 20Ne -> a + 16O
const double p_Ne20Be=0.05;	// probability for 20Ne -> 8Be + 12C
const double p_Ne20p=0.10;	// probability for 20Ne -> p + 19F
const double p_O16=0.05;	// probability for 16O -> a + 12C
const double p_C12=0.05;	// probability for 12C -> a + 8Be
const double p_Ne20s0=0.334;// probability for 20Ne -> a0 + 16Og.s.
const double p_Ne20s1=0.333;	// probability for 20Ne -> a1 + 16O*
const double p_Ne20s2=0.333;	// probability for 20Ne -> a2 + 16O*
const double p_Ne20Bes0=0.667;	// probability for 20Ne -> a0 + 16Og.s.
const float  XFWHM = 0.08;		// fwhm (in MeV) of focal plane x-axis
const float  SiFWHM = 0.25;		// fwhm (in MeV) of Si-detectorts

const int numlines1=69;		// # states for Ne20 datfile
const int numlines2=14;	// # states for O16 datfile
const int numlines3=4;	// # states for C12 datfile

// Gas cell geometry (in mm):
const float L = 155.;
const float H = 84.;
const float orgx = 0;
const float orgy = 0;        // position of source (origin, mm)

// Si detector dimentions:
// Dimentions of cell and Si detectors (in mm)
const float SiX0 = 13.;
const float SiL  = 48.;
const float SiH  = 48.;
const float SiT  = 1.;        // active area (mm)
const float SiSepX  = 12.;    // separation of Si active areas (specify in mm)
const float SiSepY1 = 15.;    // y-distance from centre line to si-detectors (specify in mm)
const float SiSepY2 = 19.;    // y-distance from centre line to si-detectors (specify in mm)
const float SiSW    = 3.;     // strip width
const float SiSSep  = 0.133;

// Canvas dimentions:
// Length conversion (factor to convert mm to fraction of Canvas):
const float Lconv = 0.80/L;        // convert to Canvas dimentions: x: 200mm * x = 0.65 -> x = 0.65/200
const float Hconv = 0.80/H;        // convert to Canvas dimentions: y: 50mm * y = 0.35 -> x = 0.35/50

// Gas cell geometry (in fractions of Canvas size):
const float cOrgx = 0.15;
const float cOrgy = 0.5;        // position of source (origin, fraction of Canvas)
float cellL, cellH, cSiL, cSiH, cSiT, cSiSepX, cSiSepY1, cSiSepY2, cSiSW, cSiSSep;       

/*------------ Constants----------------------------------------------*/
const double PI=3.14159265359;
const float amu = 931.5016;  // 1 antomic mass unit [in MeV]
const int   q_charge=2;
const float spec_rad=8.14;
const float regit=276.74;
const float m0=2809.4;
const float Epeak=200.0;         // [in MeV]

/*-- Module declaration --------------------------------------------*/
void ZeroTTreeVariables();
void print_evt(myDetData *si);
std::vector<float> get_xy(float x, float y ,float a, myDetData *si); 
float E_REL(float ec, float ed, float tc, float td);
void set_si_det(myDetData *si, float x, float y, float a, float en);

/*-- GLOBAL DECLARATIONS --------------------------------------------*/	
Float_t I[30];
Int_t   counterj;		
Char_t  flagt, IN[28], lable[28];
Float_t m_b, m_B, p_A, p_b, p_B;
std::vector<float> pos1xy;
std::vector<float> pos2x;
std::vector<float> pos2y;
std::vector<float> thet2;
std::vector<float> rawE;
std::vector<float> mass;
std::vector<float> Erel;
std::vector<float> En;
std::vector<float> Qsep;
Int_t cnt=0, ti, flag=0;
float *breakdat;
bool INSIDE = true;
bool hit = false;

// ROOT TREE VARIABLE
Int_t b_N0=-10, b_Nbreak=-10;
Float_t b_energA=-10., b_energSA=-10., b_Ex=-10., b_energSB=-10.;
Float_t b_Xpos=-1;
Int_t b_Ne20=0, b_O16=0, b_C12=0, b_Ne20p=0, b_Ne20Be=0, b_spn=0;

/*
// General variables for TTree
Int_t b_isgam=0;
Int_t b_esc;
Int_t b_comp1=0;
Int_t b_annihal=0;
Int_t b_fot1=0;
Int_t b_fot2=0;
Double_t b_compt;       
Double_t b_foto; 
Double_t b_paar; 
//Double_t b_sumenerg=-100.0;  // to prevent ch 0 from increasing
Double_t b_gamraw=0.; 
Double_t b_dist, b_gthet, b_gthet0;
*/

DetBox Si1box;          // Declare Si1box of type Box
DetBox Si2box;          // Declare Si2box of type Box
DetBox Si3box;          // Declare Si3box of type Box
DetBox Si4box;          // Declare Si4box of type Box
    
// Note, channel -1 contains all other events, so Draw() with added condition >=0 to
// only see relevant events!! OR just add "esc==0" 
//e.g. Draw("sumenerg>>hMca","sumenerg>=0 && paar<0",""); will show sum events without pair events.
//e.g. Draw("sumenerg>>hMca","sumenerg>=0 && paar<0 && foto<0",""); will show only comptpon events. 
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

// ************************************************************************************
// **************** Begin MAIN() ********************************************************************
// ************************************************************************************

int ne20_gascell()
{
	gROOT->ProcessLine("#include <vector>");	// don't know why I have to do this!!!!
//	gROOT->ProcessLine("#include \"myDetData.h\"");
//	gROOT->ProcessLine("#include \"my_si_data.h\"");
//	gROOT->ProcessLine(".L myDetData.C++");
//	gROOT->ProcessLine(".L my_si_data.C++");

    Int_t miss, i, j, p, Nj, numscat=0, nummis=0, ind=0, indx;
	Int_t b_evntno=0;
    Float_t rg, x1, y1, z2, ran, ran2, ran3, ran4;
    Float_t cx1, cy1;
    Float_t rgE, rgGam, rgE2, Estar_A=0, Estar_B, S_b, m_A, m_c=0, m_d=0, S_2a;
    Float_t E_A, thet_b, thet_c=0, thet_d=0, thet_B, phi, PHI_MAX, Xposition, thet_max, THET_MAX, Erel_tmp=-1;
	Float_t E_b=-1., E_B=-1.0, E_c=-1.0, E_d=-1.0;
	Float_t thet_in=0;
	std::vector<float> pos_tmp;
	bool lyn0=true;
    Char_t  buff[1000];

	
#ifdef DIAG_FILE
    std::ofstream diagfile;
	diagfile.open (diaglable, std::ios::out);
#endif

#ifdef PRNT_DIAG
    std::ofstream posfile;
	if(!posfile.is_open()) {
		posfile.open (plotlable, std::ios::out | std::ios::trunc);
	}

	TCanvas *c1=new TCanvas("c1","Scattering Geometry",900,480,660,350);
#endif

	myDetData *si = new myDetData;
	myDetData *gam = new myDetData;

	si->ClearEvent();
	gam->ClearEvent();

	TRandom3 *rgaus = new TRandom3(0);
	TRandom3 r(0);
//	rgaus->SetSeed(0); 	

// Create a ROOT Tree
// !!! NB !!! MUST DECLARE THE FILE AND TREE FIRST BEFORE DECLAIRING THE HISTOS AND BRANCHES !!!!!!
	TFile hfile("ne20.root","RECREATE","Demo ROOT file with histograms & trees");
	TTree *t1 = new TTree("DATA","An example of ROOT tree with a few branches");

// Create some histograms and a profile histogram
//	TH1F *h_SiE     = new TH1F("h_SiE","Si Energy",1000,0,20);
//	TH1F *h_GamE    = new TH1F("h_GamE","Gam Energy",1000,0,20);
//	TH1F *h_GamThet = new TH1F("h_GamThet","Gam Angle",1000,0,360);
/*
    t1->Branch("isgam",&b_isgam,"b_isgam/I");
    t1->Branch("esc",&b_esc,"b_esc/I");
    t1->Branch("comp1",&b_comp1,"b_comp1/I");
    t1->Branch("fot1",&b_fot1,"b_fot1/I");
    t1->Branch("fot2",&b_fot2,"b_fot2/I");
    t1->Branch("dist",&b_dist,"b_dist/D");
    t1->Branch("gthet",&b_gthet,"b_gthet/D");
    t1->Branch("gthet0",&b_gthet0,"b_gthet0/D");
    t1->Branch("foto",&b_foto,"foto/D");
    t1->Branch("compt",&b_compt,"compt/D");
    t1->Branch("paar",&b_paar,"paar/D");
    t1->Branch("annihal",&b_annihal,"annihal/I");
    t1->Branch("sumenerg",&b_sumenerg,"b_sumenerg/D");
    t1->Branch("gamraw",&b_gamraw,"b_gamraw/D");
*/
//	t1->Branch(branchname, className, &p_object, bufsize, splitlevel);
	t1->Branch("evntno",&b_evntno,"b_evntno/I");
	t1->Branch("SiData","myDetData",&si);
	t1->Branch("GamData","myDetData",&gam);
    t1->Branch("N0",&b_N0,"b_N0/I");
    t1->Branch("Nbreak",&b_Nbreak,"b_Nbreak/I");
    t1->Branch("energA",&b_energA,"b_energA/F");
    t1->Branch("energSA",&b_energSA,"b_energSA/F");
    t1->Branch("Qsep",&Qsep);
    t1->Branch("Erel",&Erel);
    t1->Branch("Ex",&b_Ex,"b_Ex/F");
    t1->Branch("energSB",&b_energSB,"b_energSB/F");
    t1->Branch("rawE",&rawE);
    t1->Branch("mass",&mass);
    t1->Branch("pos1xy",&pos1xy);
    t1->Branch("pos2x",&pos2x);
    t1->Branch("pos2y",&pos2y);
    t1->Branch("Xpos",&b_Xpos,"b_Xpos/F");
    t1->Branch("Ne20",&b_Ne20,"b_Ne20/I");
    t1->Branch("Ne20p",&b_Ne20p,"b_Ne20p/I");
    t1->Branch("Ne20Be",&b_Ne20Be,"b_Ne20Be/I");
    t1->Branch("O16",&b_O16,"b_O16/I");
    t1->Branch("C12",&b_C12,"b_C12/I");
    t1->Branch("spn",&b_spn,"b_spn/I");
    
    TH1F *hX1pos = new TH1F("hX1pos","Simulated FP position spectrum; position, x (mm); counts",2000,0,800);
    TH1F *hx0 = new TH1F("hx0","x0 of original breakup",500,0,1.1*L);
    TH1F *hE = new TH1F("hE","Energy of breakup particles; Energy [MeV]; counts",200,0,20);
    TH1F *hEx = new TH1F("hEx","Excitation energy spectrum; E_{x} [MeV]",2000,0,30);
    TH1F *hQsep = new TH1F("hQsep","Separation Q-value of breakup; Energy [MeV]; counts",500,-30.,30.);
    TH1F *hErel = new TH1F("hErel","Relative kinetic energy of breakup particles; Energy [MeV]; counts",1000,0.,30.);
    TH2F *hE_bB = new TH2F("hE_bB","Energy_#alpha vs. Energy_B; E_B [MeV]; E_#alpha [MeV]",2000,0,20,2000,0,20);
    TH2F *hEvsThet = new TH2F("hEvsThet","Energy vs. Scattering angle; #theta [deg]; E [MeV]",720,0,360,199,0,20);
    TH2F *hESivsEx = new TH2F("hESivsEx","Energy in Si's vs. Excitation Energy; E_{x} [MeV]; E_{Si} [keV]",2000,0,30,1000,0,20000);
    TH2F *hESivsXpos = new TH2F("hESivsXpos","Energy in Si's vs. FP position; Position [mm]; E_{Si} [keV]",2000,0,800,1000,0,20000);
    TH2F *hEvsThet180 = new TH2F("hEvsThet180","Energy vs. Scattering angle; #theta [deg]; E [MeV]",360,0,180,199,0,20);
    TH2F *hThetbVsThetB = new TH2F("hThetbVsThetB","#theta_#alpha vs. #theta_B; #theta_B [deg]; #theta #alpha [deg]",360,0,360,360,0,360);
    TH2F *hEHagarvsEx = new TH2F("hEHagarvsEx","Energy in Hagar vs. Excitation Energy; E_{x} [MeV]; E_{Hagar} [keV]",2000,0,30,1800,0,18000);

    gStyle->SetOptStat(kFALSE);
    hEvsThet180->GetYaxis()->SetTitleOffset(1.5);
    
// ******************** Detector setup ********************************************************
//*********************************************************************************************    
//  one push_back for each si det. in order 1, 2, 3, 4:
// 4 si det. start x positions:
 	si->Detx0.push_back(SiX0);
 	si->Detx0.push_back(SiX0+SiL+SiSepX);
 	si->Detx0.push_back(SiX0);
 	si->Detx0.push_back(SiX0+SiL+SiSepX);

//  one push_back for each si det. in order 1, 2, 3, 4:
// 4 si det. end x positions:
 	si->Detx1.push_back(SiX0+SiL);
 	si->Detx1.push_back(SiX0+SiL+SiSepX+SiL);
 	si->Detx1.push_back(SiX0+SiL);
 	si->Detx1.push_back(SiX0+SiL+SiSepX+SiL);

//  one push_back for each si det. in order 1, 2, 3, 4:
// 4 si det. start y positions:
 	si->Dety0.push_back(SiSepY1);
 	si->Dety0.push_back(SiSepY2);
 	si->Dety0.push_back(-SiSepY1);
 	si->Dety0.push_back(-SiSepY2);

//  one push_back for each si det. in order 1, 2, 3, 4:
// 4 si det. end y positions:
 	si->Dety1.push_back(SiSepY1+SiT);
 	si->Dety1.push_back(SiSepY2+SiT);
 	si->Dety1.push_back(-SiSepY1-SiT);
 	si->Dety1.push_back(-SiSepY2-SiT);
 
    PHI_MAX=atan(SiH/2./SiSepY1);			// max angle (radians) for detector width - I don't use this yet.

    std::cout<<"\n ***************** SI DETECTORS SETUP ****************"<<std::endl;
    std::cout<<"Si-1/3x: "<<(si->Detx0.at(0)-orgx)<<" to "<<(si->Detx1.at(0)-orgx)<<" mm from gas cell face"<<std::endl;
    std::cout<<"Si-2/4x: "<<(si->Detx0.at(1)-orgx)<<" to "<<(si->Detx1.at(1)-orgx)<<" mm from gas cell face"<<std::endl;
    std::cout<<"Si-1/2y: "<<(si->Dety0.at(0)-orgy)<<" to "<<(si->Dety1.at(0)-orgy)<<" mm and "<<(si->Dety0.at(1)-orgy)<<" to "<<(si->Dety1.at(1)-orgy)<<" mm "<<std::endl;
    std::cout<<"Si-3/4y: "<<(si->Dety0.at(2)-orgy)<<" to "<<(si->Dety1.at(2)-orgy)<<" mm and "<<(si->Dety0.at(3)-orgy)<<" to "<<(si->Dety1.at(3)-orgy)<<" mm "<<std::endl;

// *******************Drawing gas cell canvas: ******************************************
#ifdef PRNT_DIAG
	cellL = L*Lconv;       // length of gas cell (rel. to canvas)
	cellH = H/2.*Hconv;    // half height of gas cell (rel. to canvas)
	cSiL = SiL*Lconv;
	cSiH = SiH*Lconv;
	cSiT = SiT*Hconv;        // active area (in canvas dim.)
	cSiSepX = SiSepX*Lconv;    // separation of Si active areas (specify in mm)
	cSiSepY1 = SiSepY1*Hconv;   // y-distance from centre line to si-detectors (specify in mm)
	cSiSepY2 = SiSepY2*Hconv;   // y-distance from centre line to si-detectors (specify in mm)
	cSiSW = SiSW*Lconv;             // strip width
	cSiSSep = SiSSep*Lconv;

    c1->cd();
// Draw Source(origin):
// TEllipse(x1,y1,r1,r2,phimin,phimax,theta-rotate)
    TEllipse *ellipse = new TEllipse(cOrgx,cOrgy,0.005,0.005,0,360,0);
    ellipse->Draw();
    c1->Update();

// Draw Gas Cell: (bottom_left: cOrgx,y0; top_right: x1,y1)
    TWbox *box = new TWbox(cOrgx,cOrgy-cellH,cOrgx+cellL,cOrgy+cellH);
    box->SetLineColor(1);
    box->SetLineStyle(1);
    box->SetLineWidth(1);
    box->SetBorderSize(1);
//    box->SetFillStyle(0);  
    box->SetFillColor(19);  
    box->Draw();
    
	TLatex *tex = new TLatex(0.3155704,0.67,"Si 1");
	tex->SetLineWidth(2);
	tex->Draw();
	  tex = new TLatex(0.63,0.71,"Si 2");
	tex->SetLineWidth(2);
	tex->Draw();
	  tex = new TLatex(0.32,0.29,"Si 3");
	tex->SetLineWidth(2);
	tex->Draw();
	  tex = new TLatex(0.63,0.26,"Si 4");
	tex->SetLineWidth(2);
	tex->Draw();
	  tex = new TLatex(0.48,0.91,"Gas cell");
	tex->SetLineWidth(2);
	tex->Draw();
	  tex = new TLatex(0.04,0.52,"beam");
	tex->SetLineWidth(2);
	tex->Draw();
	TArrow *arrow = new TArrow(0.026,0.5,0.132,0.5,0.02,">");
	arrow->SetFillColor(1);
	arrow->SetFillStyle(1001);
	arrow->SetLineWidth(3);
	arrow->SetAngle(52);
	arrow->Draw();

    c1->Update();
  
// Draw Centre axis:
    TLine *line = new TLine(0.01,0.5,0.99,0.5);
    line->SetLineStyle(3);
    line->Draw();
    c1->Update();
    
	// Si1box specification
	Si1box.sx0 = cOrgx + SiX0*Lconv; 
	Si1box.sy0 = cOrgy+cSiSepY1+cSiT; 
	Si1box.sx1 = Si1box.sx0+cSiL;
	Si1box.sy1 = cOrgy+cSiSepY1;

	// Si2box specification
	Si2box.sx0 = Si1box.sx1+cSiSepX; 
	Si2box.sy0 = cOrgy+cSiSepY2+cSiT; 
	Si2box.sx1 = Si2box.sx0+cSiL;
	Si2box.sy1 = cOrgy+cSiSepY2;

	// Si3box specification
	Si3box.sx0 = cOrgx + SiX0*Lconv; 
	Si3box.sy0 = cOrgy-cSiSepY1-cSiT; 
	Si3box.sx1 = Si3box.sx0+cSiL;
	Si3box.sy1 = cOrgy-cSiSepY1;

	// Si4box specification
	Si4box.sx0 = Si3box.sx1+cSiSepX; 
	Si4box.sy0 = cOrgy-cSiSepY2-cSiT; 
	Si4box.sx1 = Si4box.sx0+cSiL;
	Si4box.sy1 = cOrgy-cSiSepY2;

    TBox *si_det1 = new TBox(Si1box.sx0,Si1box.sy0,Si1box.sx1,Si1box.sy1);
    TBox *si_det2 = new TBox(Si2box.sx0,Si2box.sy0,Si2box.sx1,Si2box.sy1);
    TBox *si_det3 = new TBox(Si3box.sx0,Si3box.sy0,Si3box.sx1,Si3box.sy1);
    TBox *si_det4 = new TBox(Si4box.sx0,Si4box.sy0,Si4box.sx1,Si4box.sy1);
	//    si_det1->SetLineColor(1);
	//    si_det1->SetLineStyle(1);
	//    si_det1->SetLineWidth(1);
    si_det1->SetFillStyle(0);  
    si_det2->SetFillStyle(0);   
    si_det3->SetFillStyle(0);  
    si_det4->SetFillStyle(0);  
    si_det1->Draw();
    si_det3->Draw();
    si_det2->Draw();
    si_det4->Draw();
    c1->Update();
#endif

// *********************************************************
// Subdivide Si's into 16 strips each, 3mm wide:
// Set channel x0 positions:

	for(int k=0;k<n_si_det;k++) {
		if(k==0) std::cout<<"Si1/3 strips: \n";
		else if(k==1) std::cout<<"\n\nSi2/4 strips: \n";
// ----------------------------------------------------------------------------------------
		for(int i=0;i<n_si_ch;i++) { 
			si->SetChanPos(si->Detx0.at(k) + i*SiSW + i*SiSSep); 	// 4*16 push_back's
// ----------------------------------------------------------------------------------------
			if(k==0 || k==1) {
				if(i<8) printf("|%4.1f+3|", (si->Detx0.at(k) + i*SiSW + i*SiSSep));         
				else if(i==8) printf("\n|%4.1f+3|", (si->Detx0.at(k) + i*SiSW + i*SiSSep));          
				else printf("|%4.1f+3|", (si->Detx0.at(k) + i*SiSW + i*SiSSep));          
			}
		}			// si det. 1, strip 1-16 (x0)
	}

	std::cout<<std::endl;
	
// ********************* 20Ne states: ***********************************   
    Float_t cent1[numlines1], gam1[numlines1][3], gamI1[numlines1][3];
//    Float_t sig1[18]={0.657,0.362,0.51,0.14,0.143,0.12,0.123,0.692,0.581,0.726,0.695,0.976,0.252,0.304,
//                        0.017,0.465,0.66,0.327};
    Float_t width1[numlines1]={0.,0.,0.,0.000028,0.019,0.0082,0.0034,0.0151,0.002,0.8,0.0021,0.019,0.0032,0.029,0.0003,0.08,
	0.024,0.045,0.013,0.58,0.,0.175,0.0003,0.04,0.0011,0.046,0.03,0.001,0.001,0.39,0.0373,0.0244,0.03,0.038,0.162,0.048,
	0.04,0.08,0.053,0.024,0.195,0.024,0.061,0.076,0.012,0.009,0.017,0.08,0.136,0.175,0.00014,0.074,0.0035,0.079,0.07,0.14,
	0.042,0.033,0.068,0.116,0.025,0.1,0.066,0.035,0.08,0.002,0.086,0.00950,0.1};

    Int_t spn1[numlines1];
    Int_t numgam1[numlines1];
    memset(&numgam1, 0, sizeof(numgam1));

// *************** Get Ne20 gamma data, energies, intensities, etc. ***********************  
    gamfile1 = fopen(gamlable1, "r");
    if(gamfile1){
    	fgets(buff, 80, gamfile1);
    	fgets(buff, 80, gamfile1);
    	for(j=0; j<numlines1; j++) {
        	fscanf(gamfile1,"%f\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f",&cent1[j],&spn1[j],&numgam1[j],&gam1[j][0],&gamI1[j][0],&gam1[j][1],&gamI1[j][1],&gam1[j][2],&gamI1[j][2]);
    	}
    	fclose (gamfile1);
    }
    else {
    	std::cout<< "\nERROR: "<< gamlable1 <<" file does not exist!!!!"<<std::endl;
    	return 0;
    }
    
//	valarray<int> myvalarray(numgam1,numlines1);
// ****************************************************************************************
// ******************  16O states: **********************************    
    Float_t cent2[numlines2], gam2[numlines2][3], gamI2[numlines2][3];
//    Float_t sig2[numlines2]={0.94,0.081,0.014,0.287,0.48,0.322,0.261,0.198,0.069,0.062,0.085};
    Float_t width2[numlines2]={0.,0.,0.,0.,0.,0.420,0.0006,0.071,0.0015,0.091,0.150,0.130,0.185,0.166};
    Int_t spn2[numlines2];
    Int_t numgam2[numlines2];
    memset(&numgam2, 0, sizeof(numgam2));
// *************** Get O16 gamma data, energies, intensities, etc. ***********************  
    gamfile2 = fopen(gamlable2, "r");
    if(gamfile2){
    	fgets(buff, 80, gamfile2);
    	fgets(buff, 80, gamfile2);
    	for(j=0; j<numlines2; j++) {
        	fscanf(gamfile2,"%f\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f",&cent2[j],&spn2[j],&numgam2[j],&gam2[j][0],&gamI2[j][0],&gam2[j][1],&gamI2[j][1],&gam2[j][2],&gamI2[j][2]);
    	}
    	fclose (gamfile2);
    }
    else {
    	std::cout<< "\nERROR: "<< gamlable2 <<" file does not exist!!!!"<<std::endl;
    	return 0;
    }

// ********************* 12C states: ********************************   
    Float_t cent3[numlines3], gam3[numlines3][3], gamI3[numlines3][3];
//    Float_t sig3[numlines3]={0.869,0.342,0.987,0.570};
    Float_t width3[numlines3]={0.,0.00001,0.0085,0.034};
    Int_t spn3[numlines3];
    Int_t numgam3[numlines3];
    memset(&numgam3, 0, sizeof(numgam3));
// *************** Get C12 gamma data, energies, intensities, etc. ***********************  
    gamfile3 = fopen(gamlable3, "r");
    if(gamfile3){
    	fgets(buff, 80, gamfile3);
    	fgets(buff, 80, gamfile3);
    	for(j=0; j<numlines3; j++) {
        	fscanf(gamfile3,"%f\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f",&cent3[j],&spn3[j],&numgam3[j],&gam3[j][0],&gamI3[j][0],&gam3[j][1],&gamI3[j][1],&gam3[j][2],&gamI3[j][2]);
    	}
    	fclose (gamfile3);
    }
    else {
    	std::cout<< "\nERROR: "<< gamlable3 <<" file does not exist!!!!"<<std::endl;
    	return 0;
    }
	    
	b_evntno=0;

// **************************************************************************************
// ********** Start event ***************************************************************
// **************************************************************************************
 
    for(Nj=0; Nj<NN; Nj++) {
  		printf("\r...progress: %.2f%%",100.*Nj/NN);
    	b_N0=Nj;
    	thet_in=0;

		si->ClearEvent();
		gam->ClearEvent();
		
		b_evntno=Nj;		
    	ZeroTTreeVariables();
    		  	
// **************************************************************************************
// ********** Determine where in the cell the reaction occured: (x1,y1) *****************
// cOrgx is entrance of cell.  (x1,y1) is interaction point. 
// **************************************************************************************
		rg = rgaus->Rndm();     	// Find x position of breakup
		
		x1 = orgx + rg*L;       	// x-pos of breakup (rel. to left side of Canvas)
		y1 = orgy + noise(0.,0.5);	// noise in [mm] - exagerated on diagram
		hx0->Fill(x1);
        pos1xy.push_back(x1-orgx);       // x-position [in mm] from gas cell face
        pos1xy.push_back(y1-orgy);       // y-position [in mm] from gas cell face

		cx1 = cOrgx + x1*Lconv;     	// x-pos of breakup (rel. to left side of Canvas)
		cy1 = cOrgy + y1*Hconv;
    
// **************************************************************************************
// ********************* choose E* of tgt ***********************************************
// **************************************************************************************
    	rgE = rgaus->Rndm(); 

// ********************************* Ne-20 state - a + 16O ******************************
    	if(rgE<=p_Ne20) {
        	indx=choose_state(sizeof(cent1)/4);
        	Estar_A = cent1[indx];
        	S_b = -4.730;               // separation energy for cluster from tgt
        	m_A = 19.99244;             // breakup nucleus: 19.99244 for 20Ne
        	m_b = 4.00260;              // light breakup ion (alpha)
        	m_B = 15.99491;             // heavy breakup ion (16O)
        	// E_A = Erecoil1[indx];         // Ne-20 state
        	// E_OUT(200, 4.00260, 8.00531, 4.00260, 8.00531, -4.23, 0., 0.)
        	// Erec = E0-Eout-Estar;
        	E_A = E0-Estar_A-E_OUT(E0, mp, m_A, mp, m_A, 0, Estar_A, THscat);
        	b_Ne20 = 1;
        	b_spn=spn1[indx];
  
		   	for(p=0;p<numgam1[indx];p++) {	// do for every gamma of Estar
		   		if(p>2) {printf("\nerror....................... p>3"); break;}
		   		rgGam=rgaus->Rndm();    		// probability of gamma to be emitted
	   			if(rgGam<=gamI1[indx][p]) MCRoot(gam1[indx][p], p, gam, loglable);	// do MCRoot() for gamma energy gam1[][]  
		   	}

		    rgE2=rgaus->Rndm();   // choose E* of heavy ion: 0.0, 6.05 or 6.92 MeV
		    if(rgE2<=p_Ne20s0) { Estar_B = 0.; }
		    else if(rgE2>p_Ne20s0 && rgE2<=(p_Ne20s0+p_Ne20s1)) {
		    	Estar_B = 6.049;        
		    	indx=1;	// state index 1 in data file of heavy ion
			   	for(p=0;p<numgam2[indx];p++) {	// do for every gamma of Estar
			   		rgGam=rgaus->Rndm();    		// probability of gamma to be emitted
		   			if(rgGam<=gamI2[indx][p]) MCRoot(gam2[indx][p], p, gam, loglable);	// do MCRoot() for gamma energy gam1[][]  
			   	}
			}
		    else {
		    	Estar_B = 6.917;
		    	indx=3;	// state index 3 in data file of heavy ion
			   	for(p=0;p<numgam2[indx];p++) {	// do for every gamma of Estar
			   		rgGam=rgaus->Rndm();    		// probability of gamma to be emitted
		   			if(rgGam<=gamI2[indx][p]) { MCRoot(gam2[indx][p], p, gam, loglable); }	// do MCRoot() for gamma energy gam1[][]  
			   	}
		    }        
    	}
    // ********************************* Ne-20 state - 8Be + 12C *****************************************************
    	else if(rgE>p_Ne20 && rgE<=(p_Ne20+p_Ne20Be)){
			for(Estar_A=0.;Estar_A<12.0;) {	// repeat search untill E* > E_sep = 11.984 MeV
				indx=choose_state(sizeof(cent1)/4);	
				Estar_A = cent1[indx];
				//std::cout<<"******* test Estar_A = "<<Estar_A<<std::endl; 
			}
		    Estar_A = cent1[indx];
		    S_b = -11.984;               // separation energy for cluster from tgt
		    m_A = 19.99244;             // breakup nucleus: 19.99244 for 20Ne
		    m_b = 8.00531;              // light breakup ion (8Be)
		    m_B = 12.00000;             // heavy breakup ion (12C)
		    m_c = 4.00260;              // breakup alphas (8Be->a+a)
		    m_d = 4.00260;              // breakup alphas (8Be->a+a)
			S_2a = 0.092;       // separation energy for cluster from tgt
	//        E_A = Erecoil1[indx];         // Ne-20 state
		    E_A = E0-Estar_A-E_OUT(E0, mp, m_A, mp, m_A, 0, Estar_A, THscat);
		    b_Ne20Be = 1;
		    b_spn=spn1[indx];
	  
		   	for(p=0;p<numgam1[indx];p++) {	// do for every gamma of Estar
		   		if(p>2) {printf("\nerror....................... p>3"); break;}
		   		rgGam=rgaus->Rndm();    		// probability of gamma to be emitted
		   		
	   			if(rgGam<=gamI1[indx][p]) MCRoot(gam1[indx][p], p, gam, loglable);	// do MCRoot() for gamma energy gam1[][]  
		   	}

		    rgE2=rgaus->Rndm();   // choose E* of heavy ion: 0.0 or 4.44 MeV
		    if(rgE2<=p_Ne20Bes0) Estar_B = 0.;
		    else {
		    	Estar_B = 4.43891;        
		    	indx=1;	// state index 1 in data file of heavy ion
			   	for(p=0;p<numgam3[indx];p++) {	// do for every gamma of Estar
			   		rgGam=rgaus->Rndm();    		// probability of gamma to be emitted
		   			if(rgGam<=gamI3[indx][p]) MCRoot(gam3[indx][p], p, gam, loglable);	// do MCRoot() for gamma energy gam1[][]  
			   	}
			} 
    	}
// ********************************* Ne-20 state - p + F ***********************************************************
		else if(rgE>(p_Ne20+p_Ne20Be) && rgE<=(p_Ne20+p_Ne20Be+p_Ne20p)){	// Ne-20 state - breakup into p + F-19
			for(Estar_A=0.;Estar_A<12.9;) {	// repeat search untill E* > E_sep = 12.844 MeV
				indx=choose_state(sizeof(cent1)/4);	
				Estar_A = cent1[indx];
				//std::cout<<"******* test Estar_A = "<<Estar_A<<std::endl; 
			}
			//std::cout<<"******* Estar_A = "<<Estar_A<<std::endl; 
			S_b = -12.844;               // separation energy for cluster from tgt
			m_A = 19.99244;              // breakup nucleus: 19.99244 for 20Ne
			m_b = 1.007825;              // light breakup ion (p)
			m_B = 18.998403;             // heavy breakup ion (19F)
	//		E_A = Erecoil1[indx];        // Ne-20 state - breakup into p + F-19
		    E_A = E0-Estar_A-E_OUT(E0, mp, m_A, mp, m_A, 0, Estar_A, THscat);
			b_Ne20p = 1;
		    Estar_B = 0.;
		    b_spn=spn1[indx];
	  
			for(p=0;p<numgam1[indx];p++) {	// do for every gamma of Estar
				rgGam=rgaus->Rndm();    		// probability of gamma to be emitted
		   		if(rgGam<=gamI1[indx][p]) MCRoot(gam1[indx][p], p, gam, loglable);	// do MCRoot() for gamma energy gam1[][]  
			}
		}        
// ********************************* O-16 state - a + 12C ***********************************************************
		else if(rgE>(p_Ne20+p_Ne20Be+p_Ne20p) && rgE<=(p_Ne20+p_Ne20Be+p_Ne20p+p_O16)){
	 		for(Estar_A=0.;Estar_A<7.2;) {	// repeat search untill E* > E_sep = 7.162 MeV
				indx=choose_state(sizeof(cent2)/4);	
				Estar_A = cent2[indx];
			}
		    S_b = -7.162;               // separation energy for cluster from tgt
		    m_A = 15.99491;             // breakup nucleus: 19.99244 for 20Ne
		    m_b = 4.00260;              // light breakup ion (alpha)
		    m_B = 12.00000;             // heavy breakup ion (12C)
	//        E_A = Erecoil2[indx];       // O-16 state
		    E_A = E0-Estar_A-E_OUT(E0, mp, m_A, mp, m_A, 0, Estar_A, THscat);
		    b_O16 = 1;
		    b_spn=spn2[indx];
	  
		   	for(p=0;p<numgam2[indx];p++) {	// do for every gamma of Estar
		   		rgGam=rgaus->Rndm();   		// probability of gamma to be emitted
	   			if(rgGam<=gamI2[indx][p]) MCRoot(gam2[indx][p], p, gam, loglable);	// do MCRoot() for gamma energy gam1[][]  
		   	}
		    rgE2=rgaus->Rndm();   // choose E* of heavy ion: 0 or 4.44 MeV
		    if(rgE2<=p_Ne20s0) Estar_B = 0.;
		    else if(rgE2>p_Ne20s0 && rgE2<=(p_Ne20s0+p_Ne20s1)){
		    	Estar_B = 4.43891;
		    	indx=1;	// state index 1 in data file of heavy ion
	   	       	for(p=0;p<numgam3[indx];p++) {	// do for every gamma of Estar
			   		rgGam=rgaus->Rndm();    		// probability of gamma to be emitted
		   			if(rgGam<=gamI3[indx][p]) MCRoot(gam3[indx][p], p, gam, loglable);	// do MCRoot() for gamma energy gam1[][]  
			   	}
		    }        
		    else {
		    	Estar_B = 7.6542;
		    	indx=2;	// state index 2 in data file of heavy ion
	   	       	for(p=0;p<numgam3[indx];p++) {	// do for every gamma of Estar
			   		rgGam=rgaus->Rndm();    		// probability of gamma to be emitted
		   			if(rgGam<=gamI3[indx][p]) MCRoot(gam3[indx][p], p, gam, loglable);	// do MCRoot() for gamma energy gam1[][]  
			   	}
		    }        
		}        
// ********************************* C-12 state - a + 8Be ***********************************************************
		else {
			for(Estar_A=0.;Estar_A<7.4;) {	// repeat search untill E* > E_sep = 7.367 MeV
				indx=choose_state(sizeof(cent3)/4);	
				Estar_A = cent3[indx]; 
			}
		    S_b = -7.367;               // separation energy for cluster from tgt
		    m_A = 12.00000;             // breakup nucleus: 19.99244 for 20Ne
		    m_b = 4.00260;              // light breakup ion (alpha)
		    m_B = 8.00531;             // heavy breakup ion (8Be)
		    m_c = 4.00260;              // breakup alphas (8Be->a+a)
		    m_d = 4.00260;              // breakup alphas (8Be->a+a)
			S_2a = 0.092;       // separation energy for cluster from tgt
	//        E_A = Erecoil3[indx];         // C-12 state
		    E_A = E0-Estar_A-E_OUT(E0, mp, m_A, mp, m_A, 0, Estar_A, THscat);
		    b_C12 = 1;
		    b_spn=spn3[indx];
	  
		   	for(p=0;p<numgam3[indx];p++) {	// do for every gamma of Estar
		   		rgGam=rgaus->Rndm();    		// probability of gamma to be emitted
	   			if(rgGam<=gamI3[indx][p]) MCRoot(gam3[indx][p], p, gam, loglable);	// do MCRoot() for gamma energy gam1[][]  
		   	}
		    Estar_B = 0.;
		}        
// ****************************************************************************************************************

		b_energSA = Estar_A;	//once the gamma MC for all targets are done, move earliers assignment here.
		b_energSB = Estar_B;

		m_A = m_A*amu;
		m_b = m_b*amu;
		m_B = m_B*amu;
		m_c = m_c*amu;
		m_d = m_d*amu;

// Now fill excitation energy spectrum: I just invert Xpos calibration, thus retaining shape of peaks
// **************************************************************************************
		if(b_Ne20==1 || b_Ne20p==1 || b_Ne20Be==1) b_Ex = noise(Estar_A,(1-width1[indx])*XFWHM/2.35);	//b_Ex is now broadened statistically
		else if(b_O16==1) b_Ex = noise(Estar_A,(1-width2[indx])*XFWHM/2.35);	//b_Ex is now broadened statistically
		else if(b_C12==1) b_Ex = noise(Estar_A,(1-width3[indx])*XFWHM/2.35);	//b_Ex is now broadened statistically
		hEx->Fill(b_Ex);

// DATA->Draw("GamData->DetEnergy:Ex>>hEHagarvsEx","","col");
		for ( int gi=0; gi<gam->SizeOfEvent('e'); gi++ ) {
			hEHagarvsEx->Fill(b_Ex,gam->GetEnergy(gi));  // here Ex are correct, but may miss showing a few gammas!     
		}
		//	ZeroTTreeVariablesMC();
		// Now fill FP position spectrum (wih calibration)
		//    Xposition = -0.2012*Estar_A*Estar_A - 26.31*Estar_A + 911.35;
		Xposition = -30.656*b_Ex + 987.52;
		//    Xposition = noise(Xposition,(1-width1[j])*XFWHM/2.35);
		hX1pos->Fill(Xposition);

// **************************************************************************************
// **************** Now consider breakup of excited state: ******************************
// **************************************************************************************

// ****************** Scattering (breakup) angle of light ion: **************************
		p_A = sqrt(2.*m_A*E_A);

		breakdat = breakup(p_A, m_A, m_b, m_B, E_A, (Estar_A-Estar_B), S_b);  // Call this action ONCE to calc breakup kinematics!
		if(breakdat[7]==1) continue;
		p_b = breakdat[0];
		p_B = breakdat[1];
		E_b = breakdat[2];
		E_B = breakdat[3];
		thet_b = breakdat[4]; // Breakup angle of one ion (anti-CLOCKwise from +x-axis) [radians]
		thet_B = breakdat[5]; // Breakup angle of second ion (CLOCKwise from +x-axis) [radians]
		thet_max = breakdat[6]; 

		hThetbVsThetB->Fill( thet_B*180./PI, thet_b*180./PI );  // NB! angle_B is stil CLOCKwise from +x-axis

		thet_B = 2.*PI - thet_B;     // NB! Now convert All second particle angles to ANTI-clockwise from +x-axis
		//	thet_B = noise(thet_B,5.*PI/180.); // ~1mm Si strip width -> ~1 deg

/*
		#ifdef DIAG_FILE
			diagfile<<"\n*********** Event: "<< Nj <<" ****************************\n";
			diagfile<<" *********** Ne20 -> a + B ****************************\n";
			diagfile<<" p0 = "<<p_A<<", E0 = "<<E_A<<" MeV \n";
			diagfile<<" breakdat[0] = p1 = "<<breakdat[0]<<" \n";
			diagfile<<" p1x = "<<breakdat[0]*cos(thet_b)<<" \n";
			diagfile<<" p2x = "<<breakdat[1]*cos(thet_B)<<" \n";
			diagfile<<" breakdat[1] = p2 = "<<breakdat[1]<<" \n";
			diagfile<<" p1y = "<<breakdat[0]*sin(thet_b)<<" \n";
			diagfile<<" p2y = "<<breakdat[1]*sin(thet_B)<<" \n";
			diagfile<<" breakdat[2] = E1 = "<<breakdat[2]<<" \n";
			diagfile<<" breakdat[3] = E2 = "<<breakdat[3]<<" \n";
		 	diagfile<<" Max. breakup cone angle (of one particle) is "<< thet_max*180./PI<<" deg!!!\n";
			diagfile<<" breakdat[5] = thet1 = "<<breakdat[4]*180./PI<<" deg\n";
			diagfile<<" breakdat[5] = thet2 = "<<breakdat[5]*180./PI<<" deg\n";
		#endif 
*/
//******************************* breakup of 8Be -> a + a *******************************************
		if(b_Ne20Be==1) {
					
			breakdat = breakup(p_b, m_b, m_c, m_d, E_b, 0, S_2a);  // Call this action ONCE to calc breakup kinematics!
			if(breakdat[7]==1) continue;
			
			E_c = breakdat[2];
			E_d = breakdat[3];
			thet_c = breakdat[4]; // Breakup angle of heavy ion (CLOCKwise from +x-axis) [radians]
			thet_d = breakdat[5]; // Breakup angle of heavy ion (CLOCKwise from +x-axis) [radians]
			thet_max = breakdat[6]; // Breakup angle of heavy (16O) ion (CLOCKwise from +x-axis) [radians]

			mass.push_back(m_c);
			mass.push_back(m_d);
			mass.push_back(m_B);

		    rawE.push_back(E_c);   // E without noise!
		    rawE.push_back(E_d);   // E without noise!
		    rawE.push_back(E_B);   // E without noise!

		    En.push_back(noise(E_c, SiFWHM/2.35));   // E WITH noise!
		    En.push_back(noise(E_d, SiFWHM/2.35));   // E WITH noise!
		    En.push_back(noise(E_B, SiFWHM/2.35));   // E WITH noise!

			Erel_tmp = E_REL(En.at(0), En.at(1), thet_c, thet_d); 
			Erel.push_back(Erel_tmp); 
	        hErel->Fill(Erel_tmp);    

/*
			#ifdef DIAG_FILE
				diagfile<<"\n*********** Event: "<< Nj <<" ****************************\n";
				diagfile<<" *********** Ne20 -> Be8 + C12 ****************************\n";
				diagfile<<" p0 = "<<p_b<<", E0 = "<<E_b<<" MeV \n";
				diagfile<<" breakdat[0] = p1 = "<<breakdat[0]<<" \n";
				diagfile<<" breakdat[1] = p2 = "<<breakdat[1]<<" \n";
				diagfile<<" breakdat[2] = E1 = "<<breakdat[2]<<" \n";
				diagfile<<" breakdat[3] = E2 = "<<breakdat[3]<<" \n";
			 	diagfile<<" Max. breakup cone angle (of one particle) is "<< thet_max*180./PI<<" deg!!!\n";
				diagfile<<" breakdat[4] = thet1 = "<<breakdat[4]*180./PI<<" deg\n";
				diagfile<<" thet1_C = "<<atan((breakdat[0]/m_c*sin(breakdat[4]))/(breakdat[0]/m_c*cos(breakdat[4])-p_b/m_b))*180./PI<<" deg\n"; // [v1 cos(thet1) - v0]/v1
				diagfile<<" breakdat[5] = thet2 = "<<breakdat[5]*180./PI<<" deg\n";
				diagfile<<" thet2_C = "<<atan((breakdat[1]/m_d*sin(breakdat[5]))/(p_b/m_b-breakdat[1]/m_d*cos(breakdat[5])))*180./PI<<" deg\n";	// [V0 - v2 cos(thet1)]/v1
			#endif 
*/

			thet_c = thet_b + thet_c;
			thet_d = thet_b - thet_d;

			thet2.push_back(thet_c);
			thet2.push_back(thet_d);
			thet2.push_back(thet_B);

/*
			#ifdef DIAG_FILE
				diagfile<<" thet_c = "<<thet_c*180./PI<<" deg\n";
				diagfile<<" thet_d = "<<thet_d*180./PI<<" deg\n";
			#endif 
*/
			for ( int i=0;i<thet2.size();i++ ) {
				pos_tmp = get_xy(x1,y1,thet2.at(i), si);	// <pos> is cleared inside function 
				set_si_det(si, pos_tmp.at(0), pos_tmp.at(1), thet2.at(i), En.at(i));
				pos2x.push_back(pos_tmp.at(0));
				pos2y.push_back(pos_tmp.at(1)); 
			}
//			std::cout<<"\n******************** size of chan evt = "<<si->SizeOfEvent('s')<<"\n";
		}

//*******************************breakup of 12C -> a + 8Be -> a + a ****************************************************
		else if(b_C12==1) {
			breakdat = breakup(p_B, m_B, m_c, m_d, E_B, 0, S_2a);  // Call this action ONCE to calc breakup kinematics!
			if(breakdat[7]==1) continue;
			
			E_c = breakdat[2];
			E_d = breakdat[3];
			thet_c = breakdat[4]; // Breakup angle of heavy ion (CLOCKwise from +x-axis) [radians]
			thet_d = breakdat[5]; // Breakup angle of heavy ion (CLOCKwise from +x-axis) [radians]

			mass.push_back(m_b);
			mass.push_back(m_c);
			mass.push_back(m_d);

		    rawE.push_back(E_b);   // E without noise!
		    rawE.push_back(E_c);   // E without noise!
		    rawE.push_back(E_d);   // E without noise!

		    En.push_back(noise(E_b, SiFWHM/2.35));   // E WITH noise!
		    En.push_back(noise(E_c, SiFWHM/2.35));   // E WITH noise!
		    En.push_back(noise(E_d, SiFWHM/2.35));   // E WITH noise!

			Erel_tmp = E_REL(En.at(1), En.at(2), thet_c, thet_d); 
			Erel.push_back(Erel_tmp); 
	        hErel->Fill(Erel_tmp);    

			thet_c = thet_B + thet_c;
			thet_d = thet_B - thet_d;

			thet2.push_back(thet_b);
			thet2.push_back(thet_c);
			thet2.push_back(thet_d);

			for ( int i=0;i<thet2.size();i++ ) {
				pos_tmp = get_xy(x1,y1,thet2.at(i), si); 
				set_si_det(si, pos_tmp.at(0), pos_tmp.at(1), thet2.at(i), En.at(i));	// push_back element for b
				pos2x.push_back(pos_tmp.at(0));
				pos2y.push_back(pos_tmp.at(1));
			}
//			std::cout<<"\n************* size of chan evt = "<<si->SizeOfEvent('s')<<"\n";
//			std::cout<<"\n************* thet2.size(): "<<thet2.size()<<", # chan hits: "<<si->SizeOfEvent('s')<<"\n";
		}
		else {
			mass.push_back(m_b);
			mass.push_back(m_B);

		    rawE.push_back(E_b);   // E without noise!
		    rawE.push_back(E_B);   // E without noise!

		    En.push_back(noise(E_b, SiFWHM/2.35));   // E WITH noise!
		    En.push_back(noise(E_B, SiFWHM/2.35));   // E WITH noise!

			Erel_tmp = E_REL(En.at(0), En.at(1), thet_b, thet_B); 
			Erel.push_back(Erel_tmp); 
	        hErel->Fill(Erel_tmp);    

			thet2.push_back(thet_b);
			thet2.push_back(thet_B);

			for ( int i=0;i<thet2.size();i++ ) {
				pos_tmp = get_xy(x1,y1,thet2.at(i), si); 	// determines the hit pos. on si detectors
				set_si_det(si, pos_tmp.at(0), pos_tmp.at(1), thet2.at(i), En.at(i));	// fills si det. structure, etc. - only valid si det. hits
				pos2x.push_back(pos_tmp.at(0));
				pos2y.push_back(pos_tmp.at(1));
			}
			hE_bB->Fill(En.at(1),En.at(0));
		}

// **************************************************************************************
// ********* Now Si detecor array hit pattern and geometry : ****************************
// **************************************************************************************
// Not really that vital tqo include random phi angle...? All it does is bring 
// more randomness, i.e. a lesser probability of detected breakup events, since
// half of the events will miss the top and bottom Si detectors.    
// If tan(phi)<0, either z2 < 0 or y_si < 0.    
		phi=89.*PI/180.;
		while((phi>PHI_MAX && phi<PI-PHI_MAX) || (phi>PI+PHI_MAX && phi<2*PI-PHI_MAX)) {
		    phi = 2.*PI*(rgaus->Rndm());     // Random, isotropic phi angle [in radians]
		}
		z2 = SiSepY1*tan(phi);

        numscat++; 
        b_Xpos = Xposition;
        b_Nbreak = numscat;     
        b_energA = E_A;
// ************ note I convert to keV for Si energies to compare with experiment ***************
		for(Int_t ii=0; ii < si->SizeOfEvent('e'); ii++){
		    hESivsEx->		Fill(b_Ex,si->GetEnergy(ii)*1000.);
		    hESivsXpos->	Fill(Xposition,si->GetEnergy(ii)*1000.);
	        hEvsThet->		Fill(si->GetTheta(ii),si->GetEnergy(ii));       // NB! angle_b is 0 - 360 deg.
	        hEvsThet180->	Fill(si->GetTheta(ii),si->GetEnergy(ii));   // NOTE angle_b is now only from 0 - 180 deg!
	        hE->			Fill(si->GetEnergy(ii));

			Qsep.push_back(si->GetEnergy(ii)*1.25024223 - b_Ex);	// energies in MeV
			hQsep->			Fill(Qsep.at(ii));  
		}
        
        t1->Fill();         // I fill t1 here to record only valid events.

        if(numscat==1234 && flag==0) print_evt(si);

		#ifdef PRNT_DIAG
//		std::cout<<"\n******************** size of thet evt = "<<thet2.size()<<"\n";
		if(numscat<=PRNT) {       // plots only first 100 detected events
			Int_t jcount = 2*si->SizeOfEvent('e');
			
			i=0;
			
			TEllipse *ellipse = new TEllipse(cx1,cy1,0.001,0.003,0,360,0);
			ellipse->Draw();

 			posfile << jcount << std::endl;
 			posfile << x1 <<"\t"<< y1 << std::endl;

			for(Int_t j=0; j < jcount; j=j+2) {
				// pos = si->GetHitPos(0); x = pos.at(0); y = pos.at(1);
				posfile << si->GetHitPos(j).at(0)<<"\t"<< si->GetHitPos(j).at(1) <<"\t"<< si->GetTheta(i) << std::endl;

				line = new TLine(cx1,cy1, cOrgx+(si->GetHitPos(j).at(0))*Lconv, cOrgy+(si->GetHitPos(j).at(1))*Hconv);
				line->SetLineStyle(1);
				if(i==0) line->SetLineColor(2);
				else if(i==1) line->SetLineColor(kBlue); 
				else if(i==2) line->SetLineColor(kGreen); 
				line->Draw();
				i++;
			}
		}    
		#endif

		#ifdef DIAG_FILE
			diagfile<<"\n\n*********** Event: "<< Nj <<" ****************************\n";
			diagfile<<" ------------- Breakup of ";
			if(b_Ne20==1) diagfile<<"Ne-20:"<<std::endl; 
			else if(b_O16==1) diagfile<<"O-16:"<<std::endl; 
			else if(b_C12==1) diagfile<<"C-12:"<<std::endl; 
			else if(b_Ne20p==1)diagfile<<"Ne-20 into p + F-19:"<<std::endl; 
			else if(b_Ne20Be==1) diagfile<<"Ne-20 into 8Be + 12C:"<<std::endl; 
			diagfile<<"pos1 = ("<<pos1xy.at(0)<<", "<<pos1xy.at(1)<<") mm from gas cell face"<<std::endl;
			diagfile<<"energA = "<<b_energA<<" MeV"<<std::endl;
			diagfile<<"energSA = "<<b_energSA<<" MeV"<<std::endl;
			diagfile<<"energSB = "<<b_energSB<<" MeV"<<std::endl;
			diagfile<<" --------------- Si detector(s) hit: ------------------\n";
			for(int i=0; i<pos2x.size();i++) {
				diagfile<<"pos2 = ("<<pos2x.at(i)<<", "<<pos2y.at(i)<<") mm, breakup angle = "<<thet2.at(i)*180./PI<<" deg\n";
			}
			for(int i=0; i < si->SizeOfEvent('d');i++){
				diagfile<<"Si "<<si->GetHitDet(i)<<", strip "<<si->GetHitChan(i)<<", E = "<<si->GetEnergy(i)<<" MeV (E_raw = "<<rawE.at(i)<<" MeV)" <<std::endl;
			}
		#endif 

	}  // end of event
//**********************************************************************************************************************	
//**********************************************************************************************************************	

	std::cout<< "\nNumber scattered: "<< numscat << "/"<<NN<<" ("<<100*numscat/NN<<"%)\n"<<std::endl;

	#ifdef DIAG_FILE
	  diagfile.close();
	#endif

	#ifdef PRNT_DIAG
	  c1->Write();
	  posfile.close();
	#endif
	 	
//	Choose event:	
//	evt=(int)(NN*r.Rndm());		//	j=(int)(0.01*(rand()%1001));
//	Char_t txt[128];
//	sprintf (txt,  "evntno>=%d && evntno<%d", evt, evt+10); 		
//	t1->Scan("GamData->DetEnergy[2]:SiData->DetEnergy[2]", txt, "colsize=16 precision=3");	//	t1->Show(evt);

	hfile.Write();
	hfile.Close();

	return 0;
}

//----------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------
 
# ifndef __CINT__
int main()
{
  ne20_gascell();
  return 0;
}
# endif
//----------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------

//------------------ZERO VARIABLES ------------------------------
void ZeroTTreeVariables()
{
	pos1xy.clear();
	pos2x.clear();
	pos2y.clear();
	thet2.clear();
	rawE.clear();
	En.clear();
 	Erel.clear();
 	Qsep.clear();
	mass.clear();
    b_Ne20=0;
    b_O16=0;
    b_C12=0;
    b_Ne20p=0;
    b_Ne20Be=0;
    hit = false;
}

//--------------- PRINT EVENT -----------------------------------
void print_evt(myDetData *si)
{
    std::cout<<"\n ***************** PRINT EVENT : "<< b_N0 <<" ***********************"<<std::endl;
//    std::cout<<"N0 = "<<b_N0<<std::endl;
    std::cout<<"Nbreak = "<<b_Nbreak<<std::endl;
    std::cout<<" ------------- Breakup of ";
    if(b_Ne20==1) std::cout<<"Ne-20:"<<std::endl; 
    else if(b_O16==1) std::cout<<"O-16:"<<std::endl; 
    else if(b_C12==1) std::cout<<"C-12:"<<std::endl; 
    else if(b_Ne20p==1) std::cout<<"Ne-20 into p + F-19:"<<std::endl; 
    else if(b_Ne20Be==1) std::cout<<"Ne-20 into 8Be + 12C:"<<std::endl; 
    std::cout<<"pos1 = ("<<pos1xy.at(0)<<", "<<pos1xy.at(1)<<") mm from gas cell face"<<std::endl;
    std::cout<<"energA = "<<b_energA<<" MeV"<<std::endl;
    std::cout<<"energSA = "<<b_energSA<<" MeV"<<std::endl;
    std::cout<<"energSB = "<<b_energSB<<" MeV"<<std::endl;

//	std::cout<<"\n******************** size of chan evt = "<<si->SizeOfEvent('s')<<"\n";
//	std::cout<<"\n******************** size of pos evt = "<<si->SizeOfEvent('p')<<"\n";
//	std::cout<<"\n******************** size of energy evt = "<<si->SizeOfEvent('e')<<"\n";
//	std::cout<<"\n******************** size of det evt = "<<si->SizeOfEvent('d')<<"\n";

    std::cout<<" --------------- Si detector(s) hit: ------------------\n";
    for(int i=0; i<pos2x.size();i++) {
    	std::cout<<"pos2 = ("<<pos2x.at(i)<<", "<<pos2y.at(i)<<") mm, breakup angle = "<<thet2.at(i)*180./PI<<" deg\n";
	}
//    std::cout<<std::endl;
    for(int i=0; i < si->SizeOfEvent('d');i++){
    	std::cout<<"Si "<<si->GetHitDet(i)<<", strip "<<si->GetHitChan(i)<<", E = "<<si->GetEnergy(i)<<" MeV (E_raw = "<<rawE.at(i)<<" MeV)" <<std::endl;
    }

//    std::cout<<"Hagar energy (sumenerg) = "<<b_sumenerg<<"\t of gamma "<<b_gamraw<<" at E* = "<<b_energSA<<std::endl;
    std::cout<<" *****************************************************"<<std::endl;
    flag=1;     // this way print only 1 event!
}

// *************** Find relative energy between breakup particles************************
float E_REL(float ec, float ed, float tc, float td)
{
	return (0.5*(ec+ed) - sqrt(ec*ed)*cos(tc + td));
}

// *************** Find (x,y) hit points at Si detectors ************************************
// Determines pos. on si det. - NOT neccessarily a valid si hit!!
std::vector<float> get_xy(float x, float y ,float a, myDetData *si) 
{ 
	float xx=0,yy=0;
	std::vector<float> pos;
	pos.clear();
		
	if(a<=PI) { 
		xx = x + SiSepY1/tan(a);
		if(xx >= si->Detx0.at(1)) yy = SiSepY2-y;
		else yy = SiSepY1-y;
		pos.push_back( x + yy/tan(a) );
		pos.push_back( y + yy );	
	}
    else { 	// angle a > PI:
		xx = x - SiSepY1/tan(a);
		if(xx >= si->Detx0.at(1)) yy = SiSepY2+y;
		else yy = SiSepY1+y;
		pos.push_back( x - yy/tan(a) );
		pos.push_back( y - yy );	
	}    
    return pos;
}

// *************** Identify the SI detector triggered ***********************************
// ***************** SI DETECTORS SETUP ****************
//    std::cout<<"Si-1/3: "<<(si->Detx0.at(0)-orgx)<<" mm to "<<(si->Detx1.at(0)-orgx)<<" mm from gas cell face"<<std::endl;
//    std::cout<<"Si-2/4: "<<(si->Detx0.at(1)-orgx)<<" mm to "<<(si->Detx1.at(1)-orgx)<<" mm from gas cell face"<<std::endl;
// Si-1/3: 13 mm to 61 mm from gas cell face
// Si-2/4: 73 mm to 121 mm from gas cell face
void set_si_det(myDetData *si, float x, float y, float a, float en)
{
//	std::cout<<"\n************************** pos_tmp = ("<<pos.at(0)<<","<<pos.at(1)<<")\n";
	if(x>=si->Detx0.at(0) && x<si->Detx1.at(0)) {	// front 
		if(y>=si->Dety0.at(0)) { 
//------ now find channel ---------------
			for(Int_t ch=0; ch<n_si_ch; ch++) {
				if(x >= si->GetChanPos(0*n_si_ch+ch) && x < si->GetChanPos(0*n_si_ch+ch+1)) {
			        si->SetHitChan(ch+1);
					si->SetHitDet(1); 	// top front
					si->SetHitPos(x,y);	// only if si dets. has been hit
					si->SetTheta(a*180./PI);
			    	si->SetEnergy(en);
			        hit=true;
			    }
			    else hit=false;
			}
		}
		else { 
//------ now find channel ---------------
			for(Int_t ch=0; ch<n_si_ch; ch++) {
				if(x >= si->GetChanPos(2*n_si_ch+ch) && x < si->GetChanPos(2*n_si_ch+ch+1)) {
			        si->SetHitChan(ch+1);
					si->SetHitDet(3); 	// bottom front
					si->SetHitPos(x,y);	// only if si dets. has been hit
					si->SetTheta(a*180./PI);
			    	si->SetEnergy(en);
			        hit=true;
			    }
			    else hit=false;
			}
		} 
	}
	else if(x>=si->Detx0.at(1) && x<si->Detx1.at(1)) {	// back
		if(y>=si->Dety0.at(1)) { 
//------ now find channel ---------------
			for(Int_t ch=0; ch<n_si_ch; ch++) {
				if(x >= si->GetChanPos(1*n_si_ch+ch) && x < si->GetChanPos(1*n_si_ch+ch+1)) {
			        si->SetHitChan(ch+1);
					si->SetHitDet(2); 	// top back
					si->SetHitPos(x,y);	// only if si dets. has been hit
					si->SetTheta(a*180./PI);
			    	si->SetEnergy(en);
			        hit=true;
			    }
			    else hit=false;
			}
		}
		else { 
//		std::cout<<"\n********************** y in bottom: "<<y<<"\n"; 
//------ now find channel ---------------
			for(Int_t ch=0; ch<n_si_ch; ch++) {
				if(x >= si->GetChanPos(1*n_si_ch+ch) && x < si->GetChanPos(1*n_si_ch+ch+1)) {
			        si->SetHitChan(ch+1);
					si->SetHitDet(4); 	// bottom back
					si->SetHitPos(x,y);	// only if si dets. has been hit
					si->SetTheta(a*180./PI);
			    	si->SetEnergy(en);
			        hit=true;
			    }
			    else hit=false;
			}
		}
	}
	else { hit=false; /* std::cout<<"************************** no si hit !! \n"; */}
}
 
