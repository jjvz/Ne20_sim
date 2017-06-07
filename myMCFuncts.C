#include <iostream>
#include <TRandom3.h>

#include "myMCFuncts.h"

//---------------------------ROUND---------------------------------
int round(float x)
{
    int X;
    X=int(x);
    if((x-X*1.0)<0.5) return X;
    else return (X+1);
}
//---------------------------------------------------------------
//---------------------------GET INDEX-----------------------------
// get left string index to a character
int lpos(char *str, char ch)
{
 int loc;
    for (loc = 0; *str != 0; loc++)
        if (*str++ == ch)
           return loc;
}
//---------------------------------------------------------------
//--------------- CHOOSE STATE -----------------------------------
// selects a random item in array with "cent" elements
Int_t choose_state(Int_t cent)
{
   	TRandom3 *rgaus = new TRandom3(0);
  	float ran = rgaus->Rndm();
    return int(cent*ran);
}
//---------------------------------------------------------------
//----------------- ADD NOISE -----------------------------------
// Every call to the rgaus->Gaus(0,1) function gives a random number between approix. -10 to +10 (centered around 0)
// with an occurance/frequency proportional to the gaus function with sigma = 1, around the mean 0, i.e. values
// close to 0 will occur more than values towards + or -10. -> 68% of occurance will be values between -1 and +1. 
Float_t noise(Float_t val, Float_t dval)        //rgaus is [-1,1]
{
   	TRandom3 *rgaus = new TRandom3(0);
	val = val + dval*rgaus->Gaus(0,1);         // +- dval [mm]
    return val;
}
//---------------------------------------------------------------

