 // Two-body Kinematics code to calculate: 
 // outgoing energy 
 // based on input: mass of projectile, target, ejectile and residual nucleus,  
 // incident energy, Q-value, excitation energy and scattering angle [in radians].
 // 
 //                                              mb
 //                                            /
 // reaction is: a + A -> B + b    --ma->--- mA
 //                                            \\
 //                                             \\     
 //                                              mB
 // JJvZ - Jan 2013
 //

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
 
#include "E_OUT.h" 
 
//float E_OUT(200, 4.00260, 8.00531, 4.00260, 8.00531, 0, 0., 0.)
float E_OUT(Float_t Ein, Float_t mp, Float_t mtgt, Float_t me, Float_t mB, Float_t Qgs, Float_t Estar, Float_t THscat)
{
// VARIABLES:     
 Float_t Eout=0, b3LI, a, b, c, d2, betaC, ynew, e3cm, ecmf, QR, ecmi; 
     
// CONSTANTS:
 Float_t amu    = 931.5016;  // 1 antomic mass unit [in MeV]
 Float_t amu2kg = 1.66e-27;
 Float_t mev2j  = 1.60e-13;
 Float_t hbar   = 197.33;         //hbar [in MeV.fm/c] 
 Float_t PI 	= 3.14159265359;
 
// INPUT:
// Float_t Ein    = 200.;          //[in MeV]
// Float_t THscat = 0.0;           // outgoing angle of ejectile [in deg]
// Float_t mp    = 4.00260;       // mass of projectile [in amu]  // 4.00260 for 4He; 1.00782 for H; 15.99491 for...
// Float_t me    = 4.00260;        // mass of ejectile [in amu]  // 3.01603 for 3He
// Float_t Estar  = 0.;          // excitation energy of tgt [in MeV]  // 20.5; 2.28; 6.012
// mtgt   = 8.00531;       //[in amu] 
// mB     = 8.00531;       //[in amu]
// Qgs    = -4.23;                  // for inelastic scatter

 mp=mp*amu;     // E0=m0c^2 [in MeV]
 mtgt=mtgt*amu;
 me=me*amu;
 mB=mB*amu;
 QR=Qgs-Estar;
 
 betaC = sqrt(Ein*(Ein+2*mp))/(Ein+mp+mtgt);
 ecmi = sqrt((mp+mtgt)*(mp+mtgt)+2*Ein*mtgt);
 ecmf = ecmi+Qgs-mp-mtgt+me+mB;
 e3cm = (ecmf*ecmf+(me+mB+Estar)*(me-mB-Estar))/(2*ecmf);
 ynew = (e3cm*e3cm/me/me)*(1-betaC*betaC);
 c = 1-ynew;
 b = cos(THscat)*(-betaC);
 a = ynew+b*b;
 d2 = b*b-a*c; 
 b3LI = (-b+sqrt(d2))/a;
    
 Eout = me*(1/sqrt(1-b3LI*b3LI)-1);
// Erec = Ein-Eout-Estar;
 return Eout;
}

