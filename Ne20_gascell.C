//=========Display Scattering Geometry for 20Ne gas cell target setup
//=========  JJvZ, Jul 7, 2014
// Usage: .L Ne20_gascell.C
//        main(NN)      // default: NN=20000 
// (where NN is number of repeats, e.g. main(200)) 

#include <stdio.h>
#include <iostream>
//#include <valarray>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

/* root includes */
#include <TH1F.h>
#include <TH2F.h>     
#include <TTree.h>
#include <TFile.h>

#define PI 3.141592654
#define muentries 40	// number of att. coeff. entries in mu.dat file
#define nSi 4 
#define mulable "/home/jjvz/codes/montecarlo/mu.dat"
//#define mulable "mu.dat"
#define gamlable1 "ne20-gammas.dat"
#define gamlable2 "o16-gammas.dat"
#define gamlable3 "c12-gammas.dat"
#define loglable "ne20-gammas-log.dat"

#define XFWHM   0.08		// fwhm (in MeV) of focal plane x-axis
#define Res     60        	// resolution % (float) of Hagar
#define D       30.0      	// distance from source (cm) (float)
#define H0      11.9      	// radius of detector face (cm)(float)
#define maxdx   35.6       	// depth of detector (cm) (float)
#define maxkan  18000       // number of channels (integer)
#define logflag 1
#define NG		500			// # gamma evrnts to print in logfile
#define numlines1	71		// # states for Ne20 datfile
#define numlines2	14		// # states for O16 datfile
#define numlines3	4		// # states for C12 datfile

FILE *gamfile1; 
FILE *gamfile2; 
FILE *gamfile3; 
FILE *mufile;
FILE *logfile;

/*------------ Constants----------------------------------------------*/
const int   q_charge=2;
const float spec_rad=8.14;
const float regit=276.74;
const float m0=2809.4;

/*-- Module declaration --------------------------------------------*/
void MCRoot(Float_t E0, Int_t pindx);
void INTERACT();
int round(float x);
int gooi_dice(float energ, float rho_i, float theta_i);
int lpos(char *str, char ch);
float gooi_dice_kn();
float COMPTON();
float convolv(float imu);
float MU(float iEmu, int type, float *MuE);
float WHERE(float iEpos, int type);
float *ROTAT(float xx, float yy, float angl);	//usage: newcoord = rotat(old_x, old_y, angle), then newx=newcoord[0] and newy = newcoord[1]
void ZeroTTreeVariables();
void ZeroTTreeVariablesMC();
void print_evt();
void BOOST(float bthet, float bphi);
Float_t *breakup(Float_t angl_a);
Float_t noise(Float_t val, Float_t dval);
Int_t choose_state(Int_t cent);

/*-- GLOBAL DECLARATIONS --------------------------------------------*/
Float_t rho_f,theta_f,phi,THETA,PHI,Edx,plength,E,dx, scalx, orgx=0., orgy=0., x2, y2; 	
Float_t MuE[200],Mum[200][4], I[30];
Float_t *coord;
Int_t   numesc,numcomp,numpaar,numfoto;	
Int_t   counterj;		
Char_t  flagt, IN[28], lable[28];
Float_t m_a, m_B, E_tot, p_A, p_a, p_B, E_alph=0., E_B, angl_B, rawEa=0, rawEB=0;
Float_t posx=0;
Int_t NN, cnt=0, ti, flag=0;
Float_t amu    = 931.5016;  // 1 antomic mass unit [in MeV]
Float_t *breakdat;
TRandom3 *rgaus = new TRandom3();
bool INSIDE = true;


// ROOT TREE VARIABLE
Int_t b_N0=-10, b_Nbreak=-10, b_Si1=0, b_Si2=0, b_Si3=0, b_Si4=0, b_Si5=0, b_Si6=0;
Int_t Sis1[16], Sis2[16], Sis3[16], Sis4[16], Sis5[16], Sis6[16];
Float_t b_energA=-10., b_energSA=-10., b_Excit=-10., b_energSB=-10., b_energB=-10., b_energa=-10., b_hoeka=-10., b_hoekB=-10.;
Float_t b_energBraw=-10., b_energaraw=-10., b_posxa=0, b_posxB=0, b_posya=0, b_posyB=0;
Float_t energy[2], hoek[2], StripPos[50][2], b_Xpos=-1, b_Qsep=0.;
Int_t b_Ne20=0, b_O16=0, b_C12=0, b_Ne20p=0;
// general variables for TTree
Int_t b_Ngam=0;
Int_t b_isgam=0;
Int_t b_esc;
Int_t b_comp1=0;
Int_t b_annihal=0;
Int_t b_fot1=0;
Int_t b_fot2=0;
Double_t b_compt;       
Double_t b_foto; 
Double_t b_paar; 
Double_t b_sumenerg=-100.0;  // to prevent ch 0 from increasing
Double_t b_gamraw=0.; 
Double_t b_dist, b_gthet, b_gthet0;
Double_t b_gamEnerg[3]={-1.,-1.,-1.};

// CLASSES
class Box
{
    public:
        double sx0;   // Length of a box
        double sy0;   // Breadth of a box
        double sx1;   // Height of a box
        double sy1;   // Height of a box
};

    TFile *f1=new TFile("Ne20.root","RECREATE"); 
    TTree *t1=new TTree("DATA","Tree data for good events");
    TTree *t2=new TTree("MCDATA","Tree data for good gamma events");

    TH1F *h_GamThet   = new TH1F("h_GamThet","Scatter angle",1000,0,180);
    TH1F *h_GamDist   = new TH1F("h_GamDist","Free Paths",1000,0,maxdx);
    TH1F *h_GamMca    = new TH1F("h_GamMca","Gamma rays; Energy [keV]",maxkan,0,maxkan);
// Note, channel -1 contains all other events, so Draw() with added condition >=0 to
// only see relevant events!! OR just add "esc==0" 
//e.g. Draw("sumenerg>>hMca","sumenerg>=0 && paar<0",""); will show sum events without pair events.
//e.g. Draw("sumenerg>>hMca","sumenerg>=0 && paar<0 && foto<0",""); will show only comptpon events. 
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

// ************************************************************************************
// **************** Begin MAIN() ********************************************************************
// ************************************************************************************
void ne20_gascell(Int_t NN=20000) 
{
    Int_t miss, i, j, p, Nj, numscat=0, nummis=0, ind=0, indx;
    Float_t rg, x0, x1, y1, x2, y2, x_start, x_end, y_si1, y_si2, z1, z2;
    Float_t rgE, rgGam, rgE2, x2B, y2B, Estar_A, Estar_B, S_a, m_A;
    Float_t Si1x0, Si1x1, Si2x0, Si2x1, Si3x0, Si3x1;
    Float_t Sis1x0[16], Sis1x1[16], Sis2x0[16], Sis2x1[16], Sis3x0[16], Sis3x1[16];
    Float_t energ_a, energ_B, E_A, thet, thet_B, phi, PHI_MAX, Xposition;

	bool lyn0=true;
    Char_t  buff[1000];
    
    Box Si1box;          // Declare Si1box of type Box
    Box Si2box;          // Declare Si2box of type Box
    Box Si3box;          // Declare Si3box of type Box
    Box Si4box;          // Declare Si4box of type Box
    
    memset(&StripPos, 0, sizeof(StripPos));

    t2->Branch("isgam",&b_isgam,"b_isgam/I");
    t2->Branch("esc",&b_esc,"b_esc/I");
    t2->Branch("comp1",&b_comp1,"b_comp1/I");
    t2->Branch("fot1",&b_fot1,"b_fot1/I");
    t2->Branch("fot2",&b_fot2,"b_fot2/I");
    t2->Branch("dist",&b_dist,"b_dist/D");
    t2->Branch("gthet",&b_gthet,"b_gthet/D");
    t2->Branch("gthet0",&b_gthet0,"b_gthet0/D");
    t2->Branch("foto",&b_foto,"foto/D");
    t2->Branch("compt",&b_compt,"compt/D");
    t2->Branch("paar",&b_paar,"paar/D");
    t2->Branch("annihal",&b_annihal,"annihal/I");
    t2->Branch("sumenerg",&b_sumenerg,"b_sumenerg/D");
    t2->Branch("gamraw",&b_gamraw,"b_gamraw/D");
    t2->Branch("gamEnerg",b_gamEnerg,"b_gamEnerg[3]/D");

// **************** Read Mu file (att. coefficients) **************************************
// Data file with gamma ray mass att. coeff's (in cm2/g): 0=Compt.; 1=Photo; 2=Pair; 3=Total 
// as function of gamma ray energy (MuE) in keV:    
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
    	cout<< "ERROR: "<< mulable <<" file does not exist!!!!"<<endl;
    	return;
    }
// ****************************************************************************************

// ********************* 20Ne states: ***********************************   
    Float_t cent1[numlines1], gam1[numlines1][3], gamI1[numlines1][3];
//    Float_t sig1[18]={0.657,0.362,0.51,0.14,0.143,0.12,0.123,0.692,0.581,0.726,0.695,0.976,0.252,0.304,
//                        0.017,0.465,0.66,0.327};
    Float_t Erecoil1[numlines1]={0.000,0.001,0.009,0.011,0.012,0.015,0.015,0.016,0.018,0.021,0.022,0.021,0.023,0.024,0.024,0.026,
	0.030,0.031,0.032,0.034,0.033,0.034,0.034,0.035,0.036,0.036,0.036,0.038,0.040,0.041,0.043,0.043,0.043,0.044,0.044,0.047,
	0.047,0.049,0.049,0.050,0.049,0.050,0.051,0.052,0.052,0.053,0.053,0.052,0.053,0.053,0.053,0.054,0.055,0.056,0.055,0.055,
	0.056,0.056,0.057,0.058,0.059,0.060,0.061,0.062,0.064,0.065,0.077,0.079,0.080,0.086,0.098};
    Float_t width1[numlines1]={0.,0.,0.,0.000028,0.019,0.0082,0.0034,0.0151,0.002,0.8,0.0021,0.019,0.8,0.0032,0.,0.029,0.0003,0.08,
	0.024,0.045,0.013,0.,0.58,0.,0.175,0.0003,0.04,0.0011,0.046,0.03,0.001,0.001,0.39,0.0373,0.0244,0.03,0.038,0.162,0.048,
	0.04,0.08,0.053,0.024,0.195,0.024,0.061,0.076,0.012,0.009,0.017,0.08,0.136,0.175,0.00014,0.074,0.0035,0.079,0.07,0.14,
	0.042,0.033,0.068,0.116,0.025,0.1,0.066,0.035,0.08,0.002,0.086,0.0095};

    Int_t numgam1[numlines1];
    memset(&numgam1, 0, sizeof(numgam1));

// *************** Get Ne20 gamma data, energies, intensities, etc. ***********************  
    gamfile1 = fopen(gamlable1, "r");
    if(gamfile1){
    	fgets(buff, 80, gamfile1);
    	fgets(buff, 80, gamfile1);
    	for(j=0; j<numlines1; j++) {
        	fscanf(gamfile1,"%f\t%d\t%f\t%f\t%f\t%f\t%f\t%f",&cent1[j],&numgam1[j],&gam1[j][0],&gamI1[j][0],&gam1[j][1],&gamI1[j][1],&gam1[j][2],&gamI1[j][2]);
    	}
    	fclose (gamfile1);
    }
    else {
    	cout<< "ERROR: "<< gamlable1 <<" file does not exist!!!!"<<endl;
    	return;
    }
    
//	valarray<int> myvalarray(numgam1,numlines1);
// ****************************************************************************************
    if(logflag==1) logfile=fopen(loglable, "w");

// ******************  16O states: **********************************    
    Float_t cent2[numlines2], gam2[numlines2][3], gamI2[numlines2][3];
//    Float_t sig2[numlines2]={0.94,0.081,0.014,0.287,0.48,0.322,0.261,0.198,0.069,0.062,0.085};
    Float_t Erecoil2[numlines2]={0.,0.013,0.013,0.017,0.017,0.032,0.034,0.046,0.051,0.054,0.060,0.060,0.069,0.081};
    Float_t width2[numlines2]={0.,0.,0.,0.,0.,0.420,0.0006,0.071,0.0015,0.091,0.150,0.130,0.185,0.166};
    Int_t numgam2[numlines2];
    memset(&numgam2, 0, sizeof(numgam2));
// *************** Get O16 gamma data, energies, intensities, etc. ***********************  
    gamfile2 = fopen(gamlable2, "r");
    if(gamfile2){
    	fgets(buff, 80, gamfile2);
    	fgets(buff, 80, gamfile2);
    	for(j=0; j<numlines2; j++) {
        	fscanf(gamfile2,"%f\t%d\t%f\t%f\t%f\t%f\t%f\t%f",&cent2[j],&numgam2[j],&gam2[j][0],&gamI2[j][0],&gam2[j][1],&gamI2[j][1],&gam2[j][2],&gamI2[j][2]);
    	}
    	fclose (gamfile2);
    }
    else {
    	cout<< "ERROR: "<< gamlable2 <<" file does not exist!!!!"<<endl;
    	return;
    }

// ********************* 12C states: ********************************   
    Float_t cent3[numlines3], gam3[numlines3][3], gamI3[numlines3][3];
//    Float_t sig3[numlines3]={0.869,0.342,0.987,0.570};
    Float_t Erecoil3[numlines3]={0.,0.009,0.027,0.043};
    Float_t width3[numlines3]={0.,0.00001,0.0085,0.034};
    Int_t numgam3[numlines3];
    memset(&numgam3, 0, sizeof(numgam3));
// *************** Get C12 gamma data, energies, intensities, etc. ***********************  
    gamfile3 = fopen(gamlable3, "r");
    if(gamfile3){
    	fgets(buff, 80, gamfile3);
    	fgets(buff, 80, gamfile3);
    	for(j=0; j<numlines3; j++) {
        	fscanf(gamfile3,"%f\t%d\t%f\t%f\t%f\t%f\t%f\t%f",&cent3[j],&numgam3[j],&gam3[j][0],&gamI3[j][0],&gam3[j][1],&gamI3[j][1],&gam3[j][2],&gamI3[j][2]);
    	}
    	fclose (gamfile3);
    }
    else {
    	cout<< "ERROR: "<< gamlable3 <<" file does not exist!!!!"<<endl;
    	return;
    }
    
// Dimentions of cell and Si detectors (in mm)
// Length conversion (factor to convert mm to fraction of Canvas):
    Float_t L = 155.;
    Float_t H = 84.;
    Float_t Lconv = 0.80/L;        // convert to Canvas dimentions: x: 200mm * x = 0.65 -> x = 0.65/200
    Float_t Hconv = 0.80/H;        // convert to Canvas dimentions: y: 50mm * y = 0.35 -> x = 0.35/50
    Float_t cellL = L*Lconv;       // length of gas cell (specify in mm)
    Float_t cellH = H/2.*Hconv;    // half height of gas cell (specify in mm)
    Float_t SiL = 48.*Lconv;
    Float_t SiH = 48.*Lconv;
    Float_t SiT = 1.*Hconv;        // active area (in canvas dim.)
    Float_t SiSepX = 12.*Lconv;    // separation of Si active areas (specify in mm)
    Float_t SiSepY1 = 15.*Hconv;   // y-distance from centre line to si-detectors (specify in mm)
    Float_t SiSepY2 = 19.*Hconv;   // y-distance from centre line to si-detectors (specify in mm)
    Float_t SiSW = 3.*Lconv;             // strip width
    Float_t SiSSep = 0.133*Lconv;

// Gas cell geometry (in fractions of Canvas size):
    Float_t orgx = 0.15, orgy = 0.5;        // position of source (origin, fraction of Canvas)
    Float_t Dcell = 0;      // distance from source to cell face (entrance window)
    
    x0=orgx+Dcell;                          // canvas x pos. of gas cell face
    x_start=x0+13.*Lconv;                // canvas x pos. of first si-detector

    PHI_MAX=atan(SiH/2./SiSepY1);	// max angle (radians) for detector width.
    
    t1->Branch("N0",&b_N0,"b_N0/I");
    t1->Branch("Nbreak",&b_Nbreak,"b_Nbreak/I");
    t1->Branch("Si1",&b_Si1,"b_Si1/I");
    t1->Branch("Si2",&b_Si2,"b_Si2/I");
    t1->Branch("Si3",&b_Si3,"b_Si3/I");
    t1->Branch("Si4",&b_Si4,"b_Si4/I");
    t1->Branch("Sis1",Sis1,"Sis1[16]/I");
    t1->Branch("Sis2",Sis2,"Sis2[16]/I");
    t1->Branch("Sis3",Sis3,"Sis3[16]/I");
    t1->Branch("Sis4",Sis4,"Sis4[16]/I");
    t1->Branch("energA",&b_energA,"b_energA/F");
    t1->Branch("energSA",&b_energSA,"b_energSA/F");
    t1->Branch("Qsep",&b_Qsep,"b_Qsep/F");
    t1->Branch("Excit",&b_Excit,"b_Excit/F");
    t1->Branch("energSB",&b_energSB,"b_energSB/F");
    t1->Branch("energy",energy,"energy[2]/F");
    t1->Branch("hoek",hoek,"hoek[2]/F");
    t1->Branch("energa",&b_energa,"b_energa/F");
    t1->Branch("energB",&b_energB,"b_energB/F");
    t1->Branch("energaraw",&b_energaraw,"b_energaraw/F");
    t1->Branch("energBraw",&b_energBraw,"b_energBraw/F");
    t1->Branch("hoeka",&b_hoeka,"b_hoeka/F");
    t1->Branch("hoekB",&b_hoekB,"b_hoekB/F");
    t1->Branch("posxa",&b_posxa,"b_posxa/F");
    t1->Branch("posxB",&b_posxB,"b_posxB/F");
    t1->Branch("posya",&b_posya,"b_posya/F");
    t1->Branch("posyB",&b_posyB,"b_posyB/F");
    t1->Branch("Xpos",&b_Xpos,"b_Xpos/F");
    t1->Branch("Ne20",&b_Ne20,"b_Ne20/I");
    t1->Branch("Ne20p",&b_Ne20p,"b_Ne20p/I");
    t1->Branch("O16",&b_O16,"b_O16/I");
    t1->Branch("C12",&b_C12,"b_C12/I");

    TCanvas *c1=new TCanvas("c1","Scattering Geometry",900,480,660,350);
    TCanvas *c2=new TCanvas("c2","Data distribution",900,80,660,350);
    c2->Divide(2,1);
    
    TH1F *hX1pos = new TH1F("hX1pos","Simulated FP position spectrum; position, x (mm); counts",2000,0,800);
    TH1F *hE = new TH1F("hE","Energy of breakup particles; Energy [MeV]; counts",200,0,20);
    TH1F *hEx = new TH1F("hEx","Excitation energy spectrum; E_{x} [MeV]",2000,0,30);
    TH1F *hQsep = new TH1F("hQsep","Separation Q-value of breakup; Energy [MeV]; counts",500,-30.,30.);
    TH2F *hE_aB = new TH2F("hE_aB","Energy_#alpha vs. Energy_B; E_B [MeV]; E_#alpha [MeV]",200,0,20,200,0,20);
    TH2F *hX = new TH2F("hX","x_^{16}O vs. x_#alpha (rel. to gas cell window)",400,0,cellL/Lconv,400,0,cellL/Lconv);
    TH2F *hEvsThet = new TH2F("hEvsThet","Energy vs. Scattering angle",720,0,360,199,0,20);
    TH2F *hESivsEx = new TH2F("hESivsEx","Energy in Si's vs. Excitation Energy; E_{x} [MeV]; E_{Si} [keV]",2000,0,30,1000,0,20000);
    TH2F *hESivsXpos = new TH2F("hESivsXpos","Energy in Si's vs. FP position; Position [mm]; E_{Si} [keV]",2000,0,800,1000,0,20000);
    TH2F *hEvsThet180 = new TH2F("hEvsThet180","Energy vs. Scattering angle",360,0,180,199,0,20);
    TH2F *hThetaVsThetB = new TH2F("hThetaVsThetB","#theta_#alpha vs. #theta_B",360,0,180,360,0,180);
    TH2F *hEHagarvsEx = new TH2F("hEHagarvsEx","Energy in Hagar vs. Excitation Energy; E_{x} [MeV]; E_{Hagar} [keV]",2000,0,30,1800,0,18000);

    gStyle->SetOptStat(kFALSE);
    hX->GetXaxis()->SetTitle("x_#alpha");
    hX->GetYaxis()->SetTitle("x_B");
    hX->GetYaxis()->SetTitleOffset(1.5);
    hEvsThet->GetXaxis()->SetTitle("#theta [deg]");
    hEvsThet->GetYaxis()->SetTitle("E [MeV]");
    hEvsThet180->GetXaxis()->SetTitle("#theta [deg]");
    hEvsThet180->GetYaxis()->SetTitle("E [MeV]");
    hEvsThet180->GetYaxis()->SetTitleOffset(1.5);
    hThetaVsThetB->GetXaxis()->SetTitle("#theta_B [deg]");
    hThetaVsThetB->GetYaxis()->SetTitle("#theta #alpha [deg]");
    c1->cd();
    
// *******************Drawing gas cell: ******************************************
// Draw Source(origin):
// TEllipse(x1,y1,r1,r2,phimin,phimax,theta-rotate)
    TEllipse *ellipse = new TEllipse(orgx,orgy,0.005,0.005,0,360,0);
    ellipse->Draw();
    c1->Update();

// Draw Gas Cell: (bottom_left: x0,y0; top_right: x1,y1)
    TWbox *box = new TWbox(x0,orgy-cellH,x0+cellL,orgy+cellH);
    box->SetLineColor(1);
    box->SetLineStyle(1);
    box->SetLineWidth(1);
    box->SetBorderSize(1);
//    box->SetFillStyle(0);  
    box->SetFillColor(19);  
    box->Draw();
    c1->Update();
  
// Draw Centre axis:
    TLine *line = new TLine(0.01,0.5,0.99,0.5);
    line->SetLineStyle(3);
    line->Draw();
    c1->Update();
    
// Si1box specification
   Si1box.sx0 = x_start; 
   Si1box.sy0 = orgy+SiSepY1+SiT; 
   Si1box.sx1 = Si1box.sx0+SiL;
   Si1box.sy1 = orgy+SiSepY1;

// Si2box specification
   Si2box.sx0 = Si1box.sx1+SiSepX; 
   Si2box.sy0 = orgy+SiSepY2+SiT; 
   Si2box.sx1 = Si2box.sx0+SiL;
   Si2box.sy1 = orgy+SiSepY2;
    
// Si3box specification
   Si3box.sx0 = x_start; 
   Si3box.sy0 = orgy-SiSepY1-SiT; 
   Si3box.sx1 = x_start+SiL;
   Si3box.sy1 = orgy-SiSepY1;

// Si4box specification
   Si4box.sx0 = Si3box.sx1+SiSepX; 
   Si4box.sy0 = orgy-SiSepY2-SiT; 
   Si4box.sx1 = Si4box.sx0+SiL;
   Si4box.sy1 = orgy-SiSepY2;

    TBox *si_det1 = new TBox(Si1box.sx0,Si1box.sy0,Si1box.sx1,Si1box.sy1);
    TBox *si_det2 = new TBox(Si2box.sx0,Si2box.sy0,Si2box.sx1,Si2box.sy1);
    TBox *si_det3 = new TBox(Si3box.sx0,Si3box.sy0,Si3box.sx1,Si3box.sy1);
    TBox *si_det4 = new TBox(Si4box.sx0,Si4box.sy0,Si4box.sx1,Si4box.sy1);
    cout<<"\n ***************** SI DETECTORS SETUP ****************"<<endl;
    cout<<"Si-1/3: "<<(Si1box.sx0-x0)/Lconv<<" mm to "<<(Si1box.sx1-x0)/Lconv<<" mm from gas cell face"<<endl;
    cout<<"Si-2/4: "<<(Si2box.sx0-x0)/Lconv<<" mm to "<<(Si2box.sx1-x0)/Lconv<<" mm from gas cell face"<<endl;
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

// *********************************************************
// Subdivide Si's into 16 strips each, 3mm wide:
  Sis1x0[0]=Si1box.sx0;
  Sis1x1[0]=Sis1x0[0]+SiSW;
  for(Int_t strip=1;strip<16;strip++) {
    Sis1x0[strip]=Sis1x1[strip-1]+SiSSep;
    Sis1x1[strip]=Sis1x0[strip]+SiSW;
  }
  printf("\nSi1/3 strips:\n");         
  for(Int_t strip=0;strip<16;strip++) {
    StripPos[strip][0]=Sis1x0[strip];
    StripPos[strip][1]=Sis1x1[strip];
    printf("%3.0f ", (StripPos[strip][0]-x0)/Lconv);         
  }      
  Sis2x0[0]=Si2box.sx0;
  Sis2x1[0]=Sis2x0[0]+SiSW;
  for(Int_t strip=1;strip<16;strip++) {
    Sis2x0[strip]=Sis2x1[strip-1]+SiSSep;
    Sis2x1[strip]=Sis2x0[strip]+SiSW;
  }
  StripPos[16][0]=Si1box.sx1;        // Gap between Si 1 and 2
  StripPos[16][1]=Si2box.sx0;
  printf("\nSi2/4 strips:\n");         
  for(Int_t strip=0;strip<16;strip++) {
    StripPos[17+strip][0]=Sis2x0[strip];
    StripPos[17+strip][1]=Sis2x1[strip];
    printf("%3.0f ", (StripPos[17+strip][0]-x0)/Lconv);         
  }      
    memset(&Sis1, 0, sizeof(Sis1));
    memset(&Sis2, 0, sizeof(Sis2));
    memset(&Sis3, 0, sizeof(Sis3));
    memset(&Sis4, 0, sizeof(Sis4));

    cout <<"\n0%---------25%---------50%----------75%----------100%"<<endl;
    srand(time(0));

// **************************************************** START EVENTS *************************************************
  for(Nj=0; Nj<NN; Nj++) { 
    b_N0=Nj;

    ZeroTTreeVariables();
    ZeroTTreeVariablesMC();
    
    rgE=0.001*(rand()%1001);    // choose E* of tgt
// ********************************* Ne-20 state - a + 16O ***********************************************************
    if(rgE<=0.8) {
        indx=choose_state(sizeof(cent1)/4);
        E_A = Erecoil1[indx];         // Ne-20 state
        Estar_A = cent1[indx];
        S_a = -4.730;               // separation energy for cluster from tgt
        m_A = 19.99244;             // breakup nucleus: 19.99244 for 20Ne
        m_a = 4.00260;              // light breakup ion (alpha)
        m_B = 15.99491;             // heavy breakup ion (16O)
        b_Ne20 = 1;
  
       	for(p=0;p<numgam1[indx];p++) {	// do for every gamma of Estar
       		rgGam=0.001*(rand()%1001);    		// probability of gamma to be emitted
   			if(rgGam<=gamI1[indx][p]) MCRoot(gam1[indx][p],p);	// do MCRoot() for gamma energy gam1[][]  
       	}

        rgE2=0.001*(rand()%1001);   // choose E* of heavy ion: 0.0, 6.05 or 6.92 MeV
        if(rgE2<=0.34) Estar_B = 0.;
        else if(rgE2>0.34 && rgE2<=0.67) {
        	Estar_B = 6.049;        
        	indx=1;	// state index 1 in data file of heavy ion
	       	for(p=0;p<numgam2[indx];p++) {	// do for every gamma of Estar
	       		rgGam=0.001*(rand()%1001);    		// probability of gamma to be emitted
	   			if(rgGam<=gamI2[indx][p]) MCRoot(gam2[indx][p],p);	// do MCRoot() for gamma energy gam1[][]  
	       	}
	    }
        else {
        	Estar_B = 6.917;
        	indx=3;	// state index 3 in data file of heavy ion
	       	for(p=0;p<numgam2[indx];p++) {	// do for every gamma of Estar
	       		rgGam=0.001*(rand()%1001);    		// probability of gamma to be emitted
	   			if(rgGam<=gamI2[indx][p]) MCRoot(gam2[indx][p],p);	// do MCRoot() for gamma energy gam1[][]  
	       	}
        }        
    }
// ********************************* Ne-20 state - p + F ***********************************************************
    else if(rgE>0.8 && rgE<=0.9){	// Ne-20 state - breakup into p + F-19
		for(Estar_A=0.;Estar_A<12.9;) {	// repeat search untill E* > E_sep = 12.844 MeV
			indx=choose_state(sizeof(cent1)/4);	
		    Estar_A = cent1[indx];
			//cout<<"******* test Estar_A = "<<Estar_A<<endl; 
		}
		//cout<<"******* Estar_A = "<<Estar_A<<endl; 
		E_A = Erecoil1[indx];        // Ne-20 state - breakup into p + F-19
		S_a = -12.844;               // separation energy for cluster from tgt
		m_A = 19.99244;              // breakup nucleus: 19.99244 for 20Ne
		m_a = 1.007825;              // light breakup ion (p)
		m_B = 18.998403;             // heavy breakup ion (19F)
		b_Ne20p = 1;
        Estar_B = 0.;
  
	    for(p=0;p<numgam1[indx];p++) {	// do for every gamma of Estar
	    	rgGam=0.001*(rand()%1001);    		// probability of gamma to be emitted
	   		if(rgGam<=gamI1[indx][p]) MCRoot(gam1[indx][p],p);	// do MCRoot() for gamma energy gam1[][]  
	    }
    }        
// ********************************* O-16 state - a + 12C ***********************************************************
    else if(rgE>0.9 && rgE<=0.95){
 		for(Estar_A=0.;Estar_A<7.2;) {	// repeat search untill E* > E_sep = 7.162 MeV
			indx=choose_state(sizeof(cent2)/4);	
		    Estar_A = cent2[indx];
		}
        E_A = Erecoil2[indx];       // O-16 state
        S_a = -7.162;               // separation energy for cluster from tgt
        m_A = 15.99491;             // breakup nucleus: 19.99244 for 20Ne
        m_a = 4.00260;              // light breakup ion (alpha)
        m_B = 12.00000;             // heavy breakup ion (12C)
        b_O16 = 1;
  
       	for(p=0;p<numgam2[indx];p++) {	// do for every gamma of Estar
       		rgGam=0.001*(rand()%1001);    		// probability of gamma to be emitted
   			if(rgGam<=gamI2[indx][p]) MCRoot(gam2[indx][p],p);	// do MCRoot() for gamma energy gam1[][]  
       	}
        rgE2=0.001*(rand()%1001);   // choose E* of heavy ion: 0 or 4.44 MeV
        if(rgE2<=0.34) Estar_B = 0.;
        else if(rgE2>0.34 && rgE2<=0.67){
        	Estar_B = 4.43891;
        	indx=1;	// state index 1 in data file of heavy ion
   	       	for(p=0;p<numgam3[indx];p++) {	// do for every gamma of Estar
	       		rgGam=0.001*(rand()%1001);    		// probability of gamma to be emitted
	   			if(rgGam<=gamI3[indx][p]) MCRoot(gam3[indx][p],p);	// do MCRoot() for gamma energy gam1[][]  
	       	}
        }        
        else {
        	Estar_B = 7.6542;
        	indx=2;	// state index 2 in data file of heavy ion
   	       	for(p=0;p<numgam3[indx];p++) {	// do for every gamma of Estar
	       		rgGam=0.001*(rand()%1001);    		// probability of gamma to be emitted
	   			if(rgGam<=gamI3[indx][p]) MCRoot(gam3[indx][p],p);	// do MCRoot() for gamma energy gam1[][]  
	       	}
        }        
    }        
// ********************************* C-12 state - a + 8Be ***********************************************************
    else {
		for(Estar_A=0.;Estar_A<7.4;) {	// repeat search untill E* > E_sep = 7.367 MeV
			indx=choose_state(sizeof(cent3)/4);	
		    Estar_A = cent3[indx]; 
		}
        E_A = Erecoil3[indx];         // C-12 state
        S_a = -7.367;               // separation energy for cluster from tgt
        m_A = 12.00000;             // breakup nucleus: 19.99244 for 20Ne
        m_a = 4.00260;              // light breakup ion (alpha)
        m_B = 8.00531;             // heavy breakup ion (8Be)
        b_C12 = 1;
  
       	for(p=0;p<numgam3[indx];p++) {	// do for every gamma of Estar
       		rgGam=0.001*(rand()%1001);    		// probability of gamma to be emitted
   			if(rgGam<=gamI3[indx][p]) MCRoot(gam3[indx][p],p);	// do MCRoot() for gamma energy gam1[][]  
       	}
        Estar_B = 0.;
    }        
// ****************************************************************************************************************

    b_energSA = Estar_A;	//once the gama MC for all targets are done, move earliers assignment here.
    b_energSB = Estar_B;
    
    E_tot = Estar_A + S_a + E_A - Estar_B;    // E_tot = E_alph + E_B = 15.891;
    if((Estar_A + E_A)<=S_a || E_tot<=0) continue;

    m_A = m_A*amu;
    m_a = m_a*amu;
    m_B = m_B*amu;
    
    p_A = sqrt(2*m_A*E_A);
    
// **************************************************************************************
// Excitation energy has been selected:
// **************************************************************************************
// Now fill excitation energy spectrum: I just invert Xpos calibration, thus retaining shape of peaks
    if(b_Ne20==1 || b_Ne20p==1) b_Excit = noise(Estar_A,(1-width1[indx])*XFWHM/2.35);	//b_Excit is now broadened statistically
    else if(b_O16==1) b_Excit = noise(Estar_A,(1-width2[indx])*XFWHM/2.35);	//b_Excit is now broadened statistically
    else if(b_C12==1) b_Excit = noise(Estar_A,(1-width3[indx])*XFWHM/2.35);	//b_Excit is now broadened statistically
    hEx->Fill(b_Excit);
	hEHagarvsEx->Fill(b_Excit,b_gamEnerg[0]);  // here the Excit energies are correct, but may miss showing a few gammas!     
	hEHagarvsEx->Fill(b_Excit,b_gamEnerg[1]);  // here the Excit energies are correct, but may miss showing a few gammas!     
	hEHagarvsEx->Fill(b_Excit,b_gamEnerg[2]);  // here the Excit energies are correct, but may miss showing a few gammas!     
//	hEHagarvsEx->Fill(b_Excit,b_sumenerg);  // here the Excit energies are correct, but may miss showing a few gammas!     
//	ZeroTTreeVariablesMC();
// Now fill FP position spectrum (wih calibration)
//    Xposition = -0.2012*Estar_A*Estar_A - 26.31*Estar_A + 911.35;
    Xposition = -30.656*b_Excit + 987.52;
//    Xposition = noise(Xposition,(1-width1[j])*XFWHM/2.35);
    hX1pos->Fill(Xposition);

// ********** Determine where in the cell the reaction occured: (x1,y1) ********
    rg=0.001*(rand()%1001);     // Find x position of breakup
    x1=x0+rg*cellL;                 // x-pos of breakup (rel. to left side of Canvas)
    y1=orgy;
    y1=y1+noise(0.,0.5)*Hconv;  // noise in [mm] - exagerated on diagram

// **************************************************************************************
// Now consider breakup of excited state:
// **************************************************************************************

// ************* Scattering (breakup) angle of light ion: **********************
    thet=PI/180.*0.01*(rand()%36001);     // Random, isotropic [in radians]
    breakdat = breakup(thet);  // Call this action ONCE to calc breakup kinematics!
    thet_B = breakdat[4]; // Breakup angle of heavy (16O) ion (CLOCKwise from +x-axis) [radians]
    thet_B = noise(thet_B,5.*PI/180.); // ~1mm Si strip width -> ~1 deg

//------------------------------------------------------------------------------
// Not really that vital to include random phi angle...? All it does is bring 
// more randomness, i.e. a lesser probability of detected breakup events, since
// half of the events will miss the top and bottom Si detectors.    
// If tan(phi)<0, either z2 < 0 or y_si < 0.    
    phi=89.*PI/180.;
    while((phi>PHI_MAX && phi<PI-PHI_MAX) || (phi>PI+PHI_MAX && phi<2*PI-PHI_MAX)) 
        phi = PI/180.*0.01*(rand()%36001);     // Random, isotropic phi angle [in radians]
    z2 = SiSepY1*tan(phi);
//------------------------------------------------------------------------------
            
    x2 = x1+SiSepY1/Hconv*Lconv/tan(thet); // x-pos of Si detection (rel. to left side of Canvas)
    y2 = orgy+SiSepY1;
    if(x2>=Si2box.sx0) y2 = orgy+SiSepY2;    
    x2B = x1+SiSepY2/Hconv*Lconv/tan(thet_B);
    y2B = orgy-SiSepY1;
    if(x2B>=Si2box.sx0) y2B = orgy-SiSepY2;    
    if(thet>PI) {
        x2 = x1-SiSepY1/Hconv*Lconv/tan(thet);
        y2 = orgy-SiSepY1;
        if(x2>=Si2box.sx0) y2 = orgy-SiSepY2;    
        x2B = x1-SiSepY2/Hconv*Lconv/tan(thet_B);
        y2B = orgy+SiSepY1;
        if(x2B>=Si2box.sx0) y2B = orgy+SiSepY2;    
    }    
// Adding some statistical noise to Si detector position and energy:
// **************************************************************
//    x2 = x2 + noise(0., 1.)*Lconv;     // simulate discrete Si strip resolution, noise in [mm]
//    x2B = x2B + noise(0., 1.)*Lconv;
// I should actually here recalculate thet and thet_B with included uncertainty/noise 
// in x2 and x2B..., but I prefer to add noise to thet_B so I don't have to recalc. the x2's! 
    if((x2>=Si1box.sx0 && x2<=Si1box.sx1) || (x2>=Si2box.sx0 && x2<=Si2box.sx1)) {
        if((x2B>=Si1box.sx0 && x2B<=Si1box.sx1) || (x2B>=Si2box.sx0 && x2B<=Si2box.sx1)) {
// *********************** Find Si strips triggered: *************************          
            miss=0;            
//            cout<<"\n begin evt: "<<b_N0<<", xa = "<<(x2-x0)/Lconv<<" x2B = "<<(x2B-x0)/Lconv<<endl;
            // Si 1/3
            for(Int_t strip=0;strip<16;strip++) {
                if(x2<=StripPos[strip][1] && x2>=StripPos[strip][0]) {
                    if(y2>orgy) {
                        Sis1[strip]=1;
                        b_Si1=1;
                        miss++;
//                        cout<<"loop si 1-3-a1"<<endl;
                        break;
                    }
                    else if(y2<orgy) {
                        Sis3[strip]=1;
                        b_Si3=1;
//                        cout<<"loop si 1-3-a2"<<endl;
                        miss++;
                        break;
                    }
                }
            }
            for(Int_t strip=0;strip<16;strip++) {
                if(x2B<=StripPos[strip][1] && x2B>=StripPos[strip][0]) {
                    if(y2B>orgy) {
                        Sis1[strip]=1;
                        b_Si1=1;
//                        cout<<"loop si 1-3-B1"<<endl;
                        miss++;
                        break;
                    }
                    else if(y2B<orgy) {
                        Sis3[strip]=1;
                        b_Si3=1;
//                        cout<<"loop si 1-3-B2"<<endl;
                        miss++;
                        break;
                    }
                }
            }
            // Si 2/4:
            for(Int_t strip=0;strip<16;strip++) {
                if(x2<=StripPos[17+strip][1] && x2>=StripPos[17+strip][0]) {
                    if(y2>orgy) {
                        Sis2[strip]=1;
                        b_Si2=1;
//                        cout<<"loop si 2-4-a1"<<endl;
                        miss++;
                        break;
                    }
                    else if(y2<orgy) {
                        Sis4[strip]=1;
                        b_Si4=1;
//                        cout<<"loop si 2-4-a2"<<endl;
                        miss++;
                        break;
                    }
                }
            }
            for(Int_t strip=0;strip<16;strip++) {
                if(x2B<=StripPos[17+strip][1] && x2B>=StripPos[17+strip][0]) {
                    if(y2B>orgy) {
                        Sis2[strip]=1;
                        b_Si2=1;
//                        cout<<"loop si 2-4-B1"<<endl;
                        miss++;
                        break;
                    }
                    else if(y2B<orgy) {
                        Sis4[strip]=1;
                        b_Si4=1;
//                        cout<<"loop si 2-4-B2"<<endl;
                        miss++;
                        break;
                    }
                }
            }
//            cout<<"detected: "<<miss<<endl;
            if(miss!=2) continue;            
//******************************************************************************88888888
            numscat++;
            posx = (x1-x0)/Lconv;       // x-position [in mm] from gas cell face
            b_posya = (y2-orgy)/Hconv;       // y-position [in mm] from gas cell face
            b_posyB = (y2B-orgy)/Hconv;       // y-position [in mm] from gas cell face
            b_Xpos = Xposition;
            E_alph = breakdat[2];
            rawEa = E_alph;   // E without noise!
            E_alph = noise(E_alph, 0.15);
            E_B = breakdat[3];
            rawEB = E_B;   // E without noise!
            E_B = noise(E_B, 0.15);
            hEvsThet->Fill(thet*180./PI,E_alph);       // NB! angle_a is 0 - 360 deg.
// ************ note I convert to keV for Si energies to compare with experiment ***************
            hESivsEx->Fill(b_Excit,E_alph*1000.);
            hESivsEx->Fill(b_Excit,E_B*1000.);
            hESivsXpos->Fill(Xposition,E_alph*1000.);
            hESivsXpos->Fill(Xposition,E_B*1000.);

// ******************** angle_B still CLOCKwise from +x axis, angle_a < 180 deg: ***************
            if(thet>PI) b_hoeka = 360 - thet*180./PI;
            else b_hoeka = thet*180./PI;
            hEvsThet180->Fill(b_hoeka,E_alph);   // NOTE angle_a is now only from 0 - 180 deg!
            hoek[0] = b_hoeka;

            if(thet_B<0) b_hoekB = -thet_B*180./PI;      // only positive angles_B and a
            else b_hoekB = thet_B*180./PI;      // only positive angles_B and a
            hEvsThet180->Fill(b_hoekB,E_B);         // NB! angle_B is CLOCKwise from +x-axis
            hoek[1] = b_hoekB;

            hThetaVsThetB->Fill(b_hoekB,b_hoeka);  // NB! angle_B is still CLOCKwise from +x-axis
            b_posxa=(x2-x0)/Lconv;
            b_posxB=(x2B-x0)/Lconv;
            hX->Fill((x2-x0)/Lconv,(x2B-x0)/Lconv);     //  offset to cell face/window
            hE->Fill(E_alph);
            hE->Fill(E_B);
            hE_aB->Fill(E_B,E_alph);

// ******************** angle_B now ANTI-clockwise from +x axis: ******
            if(thet_B<=0) thet_B = -1.*thet_B;     // NB! All angles ANTI-clockwise from +x-axis
            else if(thet_B>0) thet_B = 2*PI - thet_B;     // NB! All angles ANTI-clockwise from +x-axis
            hEvsThet->Fill(thet_B*180./PI,E_B);

            b_energaraw = rawEa;
            b_energBraw = rawEB;
            b_energA = E_A;
            b_energa = E_alph;
            energy[0]= E_alph;
            energy[1]= E_B;
            b_energB = E_B;
            b_Nbreak = numscat;        
            b_Qsep = E_alph*1.25024223 - b_Excit;	// energies in MeV
            hQsep->Fill(b_Qsep);  
// from E_tot = E* + E_recoil + Q_separation = (E_a + E_B) = E_a + E_a*(M_a/M_B) = E_a*(1.25024223)
//         e.g. 13.5 + 0.051 - 4.730 = 8.821 --->  8.821 = E_a*(1.25024223)
// we want to plot a constant, i.e. Q_sep = (E_a + E_B) - E* - E_recoil = E_a*1.25024223 - E* - E_recoil = -4.730
// to use as filter gate for only a_0 events in Si vs. E* plot.
// we neglect E_recoil as small, and assume 16O in g.s., i.e. alpha_0
// 
// to plot: DataTree->Draw("Qsep>>hQsep","","")
// and use as cut: DataTree->Draw("energa*1000.:Excit>>hESivsEx","Qsep>-6 && Qsep<-3","col")
/*
//---------------------------------------------------------------------------------------
    cout<<" ***************** PRINT EVENT ***********************"<<endl;
    cout<<"Nbreak = "<<b_Nbreak<<endl;
    cout<<"pos x0 = "<<posx<<" mm from gas cell face"<<endl;
    cout<<"pos2 a = ("<<b_posxa<<", "<<b_posya<<") mm from gas cell face, middle"<<endl;
    cout<<"pos2 B = ("<<b_posxB<<", "<<b_posyB<<") mm from gas cell face, middle"<<endl;
    cout<<"hoeka = "<<b_hoeka<<endl;
    cout<<"hoekB = "<<b_hoekB<<endl;
    cout<<"Si1 = "<<b_Si1<<"\tSi2 = "<<b_Si2<<"\tSi3 = "<<b_Si3<<endl;
    cout<<"Si4 = "<<b_Si4<<"\tSi5 = "<<b_Si5<<"\tSi6 = "<<b_Si6<<endl;
    cout<<"Sis1 strips: ";
    for(Int_t strip=0;strip<16;strip++) cout<<Sis1[strip]<<" ";
    cout<<endl;
    cout<<"Sis2 strips: ";
    for(strip=0;strip<16;strip++) cout<<Sis2[strip]<<" ";
    cout<<endl;
    cout<<"Sis3 strips: ";
    for(strip=0;strip<16;strip++) cout<<Sis3[strip]<<" ";
    cout<<endl;
    cout<<"Sis4 strips: ";
    for(strip=0;strip<16;strip++) cout<<Sis4[strip]<<" ";
    cout<<endl;
    cout<<"Sis5 strips: ";
    for(strip=0;strip<16;strip++) cout<<Sis5[strip]<<" ";
    cout<<endl;
    cout<<"Sis6 strips: ";
    for(strip=0;strip<16;strip++) cout<<Sis6[strip]<<" ";
    cout<<endl;
//---------------------------------------------------------------------------------------
*/
            t1->Fill();         // I fill t1 here to record only valid events.
            t2->Fill();         // I fill t1 here to record only valid events.
            if(numscat==1234 && flag==0) print_evt();
	ZeroTTreeVariablesMC();

        }
// *** if no Si detector was hit **************        
        else {
            x2=x1;
            y2=y1;
            x2B=x1;
            y2B=y1;
        }
    }
    else {
        x2=x1;
        y2=y1;
        x2B=x1;
        y2B=y1;
    }
    if(numscat<=100) {       // plots only first 100 detected events
        TEllipse *ellipse = new TEllipse(x1,y1,0.001,0.003,0,360,0);
        ellipse->Draw();

        line = new TLine(x1,y1,x2,y2);
        line->SetLineStyle(1);
        line->SetLineColor(2);
        line->Draw();
        TLine *line_B = new TLine(x1,y1,x2B,y2B);	//    TLine *line = new TLine(0.01,0.5,0.99,0.5);
        line_B->SetLineStyle(1);
//        line_B->SetLineWidth(2);
        line_B->SetLineColor(kBlue);
        line_B->Draw();
        c1->Update();       // uncomment to see individual events
    }    

  }  // end of event
//**********************************************************************************************************************	

  cout<< "\nNumber scattered: "<< numscat << "/"<<NN<<" ("<<100*numscat/NN<<"%)\n"<<endl;

  c2->cd(1);
  hX->Draw("col");  
//  hE_aB->Draw("col");  
  c2->cd(2);
  hEvsThet180->Draw("col"); 
//  c2->Update();
  if(logflag==1) fclose(logfile);

  c1->Write();
  f1->Write();

//  f1->Close();  
} // end of main()


//************************** BREAKUP ***************************
// Determines the breakup kinematics (energies(2,3), momenta(0,1), angles(4)), 
// given the (i) breakup angle of light ion, 
//          (ii) total energy of two ions, and
//         (iii) momentum of initial tgt nucleus A.
// NOTE: angl_a is ANTI-clockwize from +x-axis!
//       angl_B is CLOCKwize from +x-axis!
//
Float_t *breakup(Float_t angl_a)
{
    Float_t det, alpha, beta;
    static Float_t bu[5];
    
    memset(&bu, 0, sizeof(bu));

    alpha = 2*p_A*cos(angl_a);
    det = 4*p_A*p_A*cos(angl_a)*cos(angl_a) - 4*(1+m_B/m_a)*(p_A*p_A-2*m_B*E_tot);
    if(det<0) {
        cout<<"ERROR - negative determinant - check kinematics!..."<<endl;
        return 0;
    }
    else beta = sqrt(det);
// Found beta always > alpha, i.e. always use p_a = alpha + beta     
    if(beta>alpha) p_a = alpha + beta;
    else p_a = alpha - beta;
    p_a = p_a/(2*(1+m_B/m_a));  

    E_alph = p_a*p_a/2/m_a;
    E_B = E_tot - E_alph;          // def: E_tot = E_alph + E_B
    if(E_B<0) {
        cout<<"ERROR - negative determinant - check kinematics!..."<<endl;
        return 0;
    }
    else p_B = sqrt(2*m_B*E_B);
    angl_B = asin(p_a/p_B*sin(angl_a));     // breakup angle of heavy ion, CLOCKWISE from +x-axis

    if(angl_a<=PI/2. && (p_A*p_A+p_B*p_B)<p_a*p_a) 
        angl_B = PI - angl_B; // This is to change angle_B from acute to obtuse
    else if(angl_a>3.*PI/2. && (p_A*p_A+p_B*p_B)<p_a*p_a) 
        angl_B = -PI - angl_B;

    bu[0] = p_a;
    bu[1] = p_B;
    bu[2] = E_alph;
    bu[3] = E_B;
    bu[4] = angl_B;
  
    return bu;
}
//---------------------------------------------------------------
//----------------- ADD NOISE -----------------------------------
Float_t noise(Float_t val, Float_t dval)        //rgaus is [-1,1]
{
    val = val + dval*rgaus->Gaus(0,1);         // +- dval [mm]
    return val;
}
//---------------------------------------------------------------
//------------------ZERO VARIABLES ------------------------------
void ZeroTTreeVariables()
{
    b_Si1=0;
    b_Si2=0;
    b_Si3=0;
    b_Si4=0;
    posx=0;
    b_posxa=0;
    b_posxB=0;
    b_posya=0;
    b_posyB=0;
    memset(&Sis1, 0, sizeof(Sis1));
    memset(&Sis2, 0, sizeof(Sis2));
    memset(&Sis3, 0, sizeof(Sis3));
    memset(&Sis4, 0, sizeof(Sis4));
    b_Ne20=0;
    b_O16=0;
    b_C12=0;
    b_Ne20p=0;
    b_Qsep=0;

}
//--------------- PRINT EVENT -----------------------------------
void print_evt()
{
    cout<<" ***************** PRINT EVENT ***********************"<<endl;
    cout<<"N0 = "<<b_N0<<endl;
    cout<<"Nbreak = "<<b_Nbreak<<endl;
    cout<<"pos x0 = "<<posx<<" mm from gas cell face"<<endl;
    cout<<"pos2 a = ("<<b_posxa<<", "<<b_posya<<") mm from gas cell face, middle"<<endl;
    cout<<"pos2 B = ("<<b_posxB<<", "<<b_posyB<<") mm from gas cell face, middle"<<endl;
    cout<<"hoeka = "<<b_hoeka<<endl;
    cout<<"hoekB = "<<b_hoekB<<endl;
    cout<<"Si1 = "<<b_Si1<<"\tSi2 = "<<b_Si2<<endl;
    cout<<"Si3 = "<<b_Si3<<"\tSi4 = "<<b_Si4<<endl;
    cout<<"Sis1 strips: ";
    for(Int_t strip=0;strip<16;strip++) cout<<Sis1[strip]<<" ";
    cout<<endl;
    cout<<"Sis2 strips: ";
    for(Int_t strip=0;strip<16;strip++) cout<<Sis2[strip]<<" ";
    cout<<endl;
    cout<<"Sis3 strips: ";
    for(Int_t strip=0;strip<16;strip++) cout<<Sis3[strip]<<" ";
    cout<<endl;
    cout<<"Sis4 strips: ";
    for(Int_t strip=0;strip<16;strip++) cout<<Sis4[strip]<<" ";
    cout<<endl;
    cout<<"energA = "<<b_energA<<endl;
    cout<<"energSA = "<<b_energSA;
    if(b_Ne20==1) cout<<"   --- breakup of Ne-20"<<endl; 
    else if(b_O16==1) cout<<"   --- breakup of O-16"<<endl; 
    else if(b_C12==1) cout<<"   --- breakup of C-12"<<endl; 
    else if(b_Ne20p==1) cout<<"   --- breakup of Ne-20 into p + F-19"<<endl; 
    cout<<"energSB = "<<b_energSB<<endl;
    cout<<"energa = "<<b_energa<<"\t energaraw = "<<b_energaraw<<endl;
    cout<<"energB = "<<b_energB<<"\t energBraw = "<<b_energBraw<<endl;
    cout<<" *****************************************************"<<endl;
    cout<<"Hagar energy (sumenerg) = "<<b_sumenerg<<"\t of gamma "<<b_gamraw<<" at E* = "<<b_energSA<<endl;
    cout<<" *****************************************************"<<endl;
    flag=1;     // this way print only 1 event!
}
//--------------- CHOOSE STATE -----------------------------------
// selects a random item in array with "cent" elements
Int_t choose_state(Int_t cent)
{
    Float_t ran = 0.001*(rand()%1001);
    return int(cent*ran);
}
//---------------------------------------------------------------
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
//------------- BEGIN MCROOT() -----------------------------------------
void MCRoot(Float_t E0, Int_t pindx)
{
    Int_t   j, l, bini, NI;
    Float_t RHO, DX0, DY0, THETA_MAX;
    b_Ngam++;
    b_isgam=1;
    
    THETA_MAX=180.*atan(H0/D)/PI;	        //max angle (deg) for detector face area.
                                                            //anders vat die program eeue!!
    numesc=0;
    numcomp=0;
    numpaar=0;
    numfoto=0;
        
// Begin event:>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    b_gamraw=E0;	//stores initial gamma energy
    Edx=0.;
  
	THETA=PI*0.01*(rand()%int(100.*THETA_MAX+1))/180.;    //gamma incident angle, theta (rad)
	PHI=PI*0.01*(rand()%36001)/180.;                  //gamma incident angle, phi (rad)
	RHO=D*tan(THETA);                              //gamma incident radial position, rho

	plength=WHERE(E0,3);      //path length inside detector
	dx=plength*cos(THETA);
	b_dist=plength;
	h_GamDist->Fill(plength);
	b_gthet0=THETA*180./PI;	// initial gamma incident angle on HAGAR
	phi=PHI;
	rho_f=RHO+(dx*tan(THETA));                              //new radial distance after plength

	if(dx>maxdx || rho_f>H0) {                      //if outside detector
		++numesc;
		b_esc=1;
	}
	if(dx<=maxdx || rho_f<=H0) INSIDE=true;
	else INSIDE=false;
		
	for(E=E0,theta_f=THETA;((INSIDE && (E!=0.)) || (b_annihal==1) || (b_annihal==2));) { //loop until out of detector or fully absorbed (E=0.)
// ****************************** INTERACTION************************************
		INTERACT();         //At point B
		h_GamThet->Fill(theta_f*180./PI);
	}
	if(Edx>0.) {
		b_sumenerg=Edx;
		h_GamMca->Fill(Edx);
	}
	Edx=0.;

    if(logflag==1 && b_Ngam<NG) {
    	fprintf(logfile, "Event: %d , gam evt: %d E* = %.3f MeV, gamma %.0f detected with Energy = %.0f keV \n", b_N0, b_Ngam, b_energSA, b_gamraw, b_sumenerg);
    	fprintf(logfile, "_________________________________________________________________________________ \n");
    }
    b_gamEnerg[pindx]=b_sumenerg;
//	t1->Fill();   
//	ZeroTTreeVariablesMC();
}
//-------------END MCROOT()-------------------------------------------
//---------------------------INTERACT-----------------------------
void INTERACT()
{
    float R, p_abs, p_foto, p_compt, p_paar, dEpaar, paarrho, paarx;       

    p_compt=1-exp(-MU(E,0,MuE)*plength);
    p_foto=1-exp(-MU(E,1,MuE)*plength);
    if(E>=1022.) p_paar=1-exp(-MU(E,2,MuE)*plength);
    else p_paar=0.;

    p_abs=p_foto+p_compt+p_paar;		     //absorbsie

    R=0.001*(rand()%1001);
//------------------- FOTO ABSORBPTION -------------------------------
    if(R <= (p_foto/p_abs)) {   // Foto electriese absorbsie
        ++numfoto;
        if(b_comp1==0) b_fot1=1;
        else b_fot2=1;
  		if(logflag==1 && b_Ngam<NG) 
    		fprintf(logfile, "foto abs.: E = %.0f (%d); numfoto = %d\n", E, b_annihal, numfoto);
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
        if(logflag==1 && b_Ngam<NG) fprintf(logfile, "paar: E = %.0f; numpaar = %d\n", E, numpaar);

        Edx=Edx+convolv((E-1022.)/2.);	 	// electron absorbed taking half the remainder (E0-1022):
        Edx=Edx+convolv((E-1022.)/2.);		// positron kinetic energy absorbed as it is slowed down before being annihilated:
       	E = 511.;		// gamma 1: the e+ annihilates, creating two opposite 511 keV gammas.
        b_paar=Edx;
        paarx=dx;
        paarrho=rho_f;
        BOOST(PI*0.01*(rand()%36001)/180., PI*0.01*(rand()%36001)/180.);	// 511 goes anywhere.
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

    theta_f=gooi_dice_kn();	                //Klein-Nishina angular cross-section
    b_gthet=theta_f*180./PI;              //!! this overrides previous compt. scat. theta!!!
    phi=PI*0.01*(rand()%36001)/180.;                //isotropic
    Egam=E/(1.+E*(1.-cos(theta_f))/511.);       //scattered/secondary gamma energy
    Ee=E-Egam;                                  //energy deposited to electrons

    Edx=Edx+convolv(Ee);
    b_compt=Edx;    // CAREFULL!! - I'm overwriting previous comptons during multiple compt. scat.
    if(logflag==1 && b_Ngam<NG) {
    	fprintf(logfile, "compt. scatt. of %.0f (%d) at thet=%.2f deg, phi=%.2f deg; numcomp = %d\n", E, b_annihal, theta_f*180./PI, phi*180./PI, numcomp);
    	fprintf(logfile, "\tEgam out = %.0f keV (%d), %.0f keV deposited\n", Egam, b_annihal, Ee);         
    }
    return Egam;   // remember Egam is UNconvoluted (100% resolution)!!
}
//---------------------------COMPTON-------------------------------

//---------------------------KLEIN-NISHINA ANGLE-------------------
float gooi_dice_kn()        //hoek in radians
{
    float thet_i, Rkn, integral, kn, prob, facto, Ykn;
    int ikn, succ2;

    for(integral=0., ikn=0; ikn<=1800; ikn++) {
	thet_i=ikn*PI/1800.;
    facto=1.+E*(1.-cos(thet_i))/511.;
	kn=pow((1./facto),2)*(facto+(1./facto)-pow(sin(thet_i),2));
	integral=integral+kn;
    }
    for(succ2=0; succ2!=1; ) {
        thet_i=PI*0.01*(rand()%18001)/180.;              //random thet between 0 and pi radians     
        facto=1.+E*(1.-cos(thet_i))/511.;
        Ykn=pow((1./facto),2)*(facto+(1./facto)-pow(sin(thet_i),2))/integral;
        Rkn=0.001*(rand()%1001);
        if(Rkn <= Ykn) succ2=1;
        else ;
    }
    return thet_i;
}
//---------------------------KLEIN-NISHINA ANGLE-------------------

//---------------------------ROUND---------------------------------
int round(float x)
{
    int X;
    X=int(x);
    if((x-X*1.0)<0.5)
	return X;
    else
	return (X+1);
}
//---------------------------ROUND---------------------------------

//---------------------------CONVOLVE------------------------------
//This function spreads the energy using a gaussian random distribution 
// according to the input selected resolution.
//It returns an energy value, distributed accordingly. 

float convolv(float imu)
{
    int succ;
    float Xg, Yg, rg;     
    
    for(succ=0; succ!=1; ) {
        Xg=(imu-3.322*Res)+0.006644*Res*(rand()%1000);                   //full-width-at-1%-max
//        Yg=(1./(1.0645*Res))*exp(-pow((Xg-imu)/(0.6006*Res),2.0));       //gauss function normalised to 1
        Yg=1.*exp(-pow((Xg-imu)/(0.6006*Res),2.0));       //gauss function normalised to 1
        rg=0.001*(rand()%1001);
        if(rg <= Yg) succ=1;
        else ;
    }
    return Xg;
}
//---------------------------CONVOLVE------------------------------

//---------------------------GET INDEX-----------------------------
// get left string index to a character

int lpos(char *str, char ch)
{
 int loc;
    for (loc = 0; *str != 0; loc++)
        if (*str++ == ch)
           return loc;
}
//---------------------------GET INDEX-----------------------------

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
        cout<<"problem in MU().... no coefficient found!!"<<endl;
        return 1000000.;
     }
  }
}

//--------------------------- WHERE(E,mu) -----------------------------
// Determines the dx (path length) into the medium where the interaction took place:
float WHERE(float iEpos, int type)
{  
  float YY,XX=0.;
  for(YY=0.;YY==0.;){
      YY=0.001*(rand()%1001);
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

//--------------------------- ZeroTTreeVariablesMC() -----------------------------
void ZeroTTreeVariablesMC()
{
   b_dist=0.;
   b_gthet=0.;
   b_gthet0=0.;
   b_sumenerg=-1.0;      
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
   b_gamEnerg[0]=-1.;
   b_gamEnerg[1]=-1.;
   b_gamEnerg[2]=-1.;
}

//--------------------------- ROTATE -----------------------------
// rotates the (x,y) coordinates in 2D through angle theta:

float *ROTAT(float xx, float yy, float angl)
{
    static float newcoord[2];
    newcoord[0] = xx*cos(angl)- yy*sin(angl);
    newcoord[1] = xx*sin(angl)+ yy*cos(angl);
    return newcoord;
}
//--------------------------- ROTATE -----------------------------
//****************************************************************************************************


    
    
