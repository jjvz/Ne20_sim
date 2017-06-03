#ifndef CONSTANTS_H
#define CONSTANTS_H

extern float MuE[200], Mum[200][4];
extern float plength, E,theta_f,phi,rho_f,dx;
extern float *coord;
extern bool INSIDE;
extern const int muentries;	// number of att. coeff. entries in mu.dat file
extern const double PI;
extern const int NN;
extern const int n_si_det;
extern const int n_si_ch;
extern const int n_gam_det;
extern const int n_gam_ch;

extern const double p_Ne20;	// probability for 20Ne -> a + 16O
extern const double p_Ne20Be;	// probability for 20Ne -> 8Be + 12C
extern const double p_Ne20p;	// probability for 20Ne -> p + 19F
extern const double p_O16;	// probability for 16O -> a + 12C
extern const double p_C12;	// probability for 12C -> a + 8Be
extern const double p_Ne20s0;// probability for 20Ne -> a0 + 16Og.s.
extern const double p_Ne20s1;	// probability for 20Ne -> a1 + 16O*
extern const double p_Ne20s2;	// probability for 20Ne -> a2 + 16O*
extern const double p_Ne20Bes0;	// probability for 20Ne -> a0 + 16Og.s.extern const double E0;        // [in MeV]
extern const double THscat;          // outgoing angle of ejectile [in deg]
extern const double mp;      // mass of projectile [in amu]  // 4.00260 for 4He; 1.00782 for H; 15.99491 for...

// FOR MCROOT - GAMMA DECAY:
extern const float XFWHM;		// fwhm (in MeV) of focal plane x-axis
extern const float Res;       	// resolution % (float) of Hagar
extern const float D;     	// distance from source (cm) (float)
extern const float H0;     	// radius of detector face (cm)(float)
extern const float maxdx;       	// depth of detector (cm) (float)
extern const int maxkan;      // number of channels (integer)
extern const int NG;	// # gamma evrnts to print in logfile
extern const int numlines1;		// # states for Ne20 datfile
extern const int numlines2;	// # states for O16 datfile
extern const int numlines3;	// # states for C12 datfile

/*
// Si detector dimentions:
// Dimentions of cell and Si detectors (in mm)
// Length conversion (factor to convert mm to fraction of Canvas):
extern const float L;
extern const float H;
extern const float Lconv;        // convert to Canvas dimentions: x: 200mm * x = 0.65 -> x = 0.65/200
extern const float Hconv;        // convert to Canvas dimentions: y: 50mm * y = 0.35 -> x = 0.35/50
extern const float cellL;       // length of gas cell (specify in mm)
extern const float cellH;    // half height of gas cell (specify in mm)
extern const float SiL;
extern const float SiH;
extern const float SiT;        // active area (in canvas dim.)
extern const float SiSepX;    // separation of Si active areas (specify in mm)
extern const float SiSepY1;   // y-distance from centre line to si-detectors (specify in mm)
extern const float SiSepY2;   // y-distance from centre line to si-detectors (specify in mm)
extern const float SiSW;             // strip width
extern const float SiSSep;

// Gas cell geometry (in fractions of Canvas size):
extern const float orgx;             // strip width
extern const float orgy;             // strip width
*/
//TRandom3 *rgaus = new TRandom3(0);

int round(float x);
Int_t choose_state(Int_t cent);
int lpos(char *str, char ch);
Float_t noise(Float_t val, Float_t dval);

// CLASSES
class DetBox
{
    public:
        double sx0;   // Length of a box
        double sy0;   // Breadth of a box
        double sx1;   // Height of a box
        double sy1;   // Height of a box
};
 
#endif

