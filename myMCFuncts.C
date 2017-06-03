#include <iostream>
#include <TRandom3.h>

#include "myMCFuncts.h"



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
Float_t noise(Float_t val, Float_t dval)        //rgaus is [-1,1]
{
   	TRandom3 *rgaus = new TRandom3(0);
	val = val + dval*rgaus->Gaus(0,1);         // +- dval [mm]
    return val;
}
//---------------------------------------------------------------

