//************************** BREAKUP ***************************
// Determines the breakup kinematics (energies(2,3), momenta(0,1), angles(4)), 
// given the (i) breakup angle of light ion, 
//          (ii) total energy of two ions, and
//         (iii) momentum of initial tgt nucleus A.
// NOTE: angl1 is ANTI-clockwize from +x-axis!
//       angl2 is CLOCKwize from +x-axis!
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

#include "BREAKUP.h"

const double PI=3.14159265359;
int countp, countpn;

float *breakup(Float_t p_0, Float_t m_0, Float_t m_1, Float_t m_2, Float_t E_0, Float_t E_st, Float_t Q)
{
    Float_t det, alpha=0, beta=0, angl_1=0, angl_2=0, p_1=0, p_2=0, E_1=0, E_2=0, Etot=0, th_m=0;
    Int_t err1=0, err2=0, err0=0;
    static Float_t bu[8];
	TRandom3 *rgaus = new TRandom3(0);
    
    memset(&bu, 0, sizeof(bu));

//	E_tot = E_A + Estar_A - Estar_B + S_b;    // E_tot = E_b + E_B = 15.891;
	Etot = E_0 + E_st + Q;    // E_tot = E_b + E_B = 15.891;
	if(Etot>0 && p_0 <= sqrt(2.*(m_1+m_2)*Etot)) {	// breakup possible
//	    tm = (1+m2/m1)*2.*m2*Et/p0/p0 - m2/m1;
//		only possible if(tm>0 && tm<=1)

	    th_m = MAX_CONE(m_1, m_2, Etot, p_0);	// neg. if E_tot > 4*E_b
		angl_1 = th_m*rgaus->Rndm();     // Random, isotropic [in radians]
//		angl_1=205.*PI/180.;
//			std::cout<<" deg, and angl_1 = "<<angl_1*180./PI<<" deg\n";

		alpha = p_0*cos(angl_1);
		det = p_0*p_0*cos(angl_1)*cos(angl_1) - (1.+m_2/m_1)*(p_0*p_0-2.*m_2*Etot);
		if(det>=0) beta = sqrt(det);
		else {
			err1++;
			//std::cout<<"\nERROR 1 - negative determinant - check kinematics!...\n";
			det=0;
		}

// p_1 = alpha + beta (~96%, always > 0)  
//	OR   
// p_1 = alpha - beta (4%, only if beta < alpha ) - here only when Theta_b < 90 deg, and p_b > p_A.
// beta >= 0
// alpha+beta > 0
// Generally beta > alpha, i.e. generally p1 = alpha + beta     

		if(beta>=alpha) { 
			p_1 = alpha + beta;  
		}
		else { 
			p_1 = alpha - beta;	 	// however, both + and - solutons are valid, though I found the + and - solutions equal	
		}	
//		if(angl_1>PI/2.) std::cout<<" thet_b = "<<angl_1*180./PI<<" deg\n"; 

		p_1 = p_1/(1+m_2/m_1);  
		E_1 = p_1*p_1/2./m_1;
		E_2 = Etot - E_1;          // def: E_tot = E_b + E_B
		if(E_2>0) {
			p_2 = sqrt(2.*m_2*E_2);

			angl_2 = asin(p_1/p_2*sin(angl_1));     // breakup angle of heavy ion, CLOCKWISE from +x-axis
//			if(p_1<0) std::cout<<"\n ** p0 = "<<p_0<<"\talpha = "<<alpha<<"\t thet_B = "<<angl_2*180./PI<<" deg\n"; 
			if(p_1*sin(angl_1)<0) angl_2 = 2.*PI + angl_2;	// this is to correct for the incorrect acute angle limit of asin() function

			bu[0] = p_1;
			bu[1] = p_2;
			bu[2] = E_1;
			bu[3] = E_2;
			bu[4] = angl_1;
			bu[5] = angl_2;
			bu[6] = th_m;
		}
		else {	//std::cout<<"\n..... breakup NOT possible - p_0 > sqrt(2.*(m_1+m_2)*Etot) or Etot < 0 ...\n"; 
			err0=1; 
		}
	}
	else {	//std::cout<<"\n..... breakup NOT possible - p_0 > sqrt(2.*(m_1+m_2)*Etot) or Etot < 0 ...\n"; 
		err0=1; 
	}
	
	bu[7] = err0;
	return bu;
}

// *************** Find max breakup angle cone of one particle ************************
float MAX_CONE(float m1, float m2, float Et, float p0)
{
    float tm=0;
    tm = (1+m2/m1)*2.*m2*Et/p0/p0 - m2/m1;	// for m1=m2=m0/2 -> tm = 4.*m2*Et/p0/p0 - 1 = 2.*m2/m0*(E_0+E*+Q)/E_0 - 1 = (E_0+E*+Q)/E_0 - 1 = (E*+Q)/E_0
	if(tm>0 && tm<=1) {
		tm = asin(sqrt(tm));
	}
	else if(tm>1) tm = 2.*PI;
	else tm=0;	// no breakup possible!!

	return tm;
}

