// To compile in ROOT with CINT:
void ne20_cint()
{
  gROOT->ProcessLine(".L BREAKUP.C++");
  gROOT->ProcessLine(".L E_OUT.C++");
  gROOT->ProcessLine(".L MCGAMMA.C++");
  gROOT->ProcessLine(".L myMCFuncts.C++");
  gROOT->ProcessLine(".L myDetData.C++"); 
  gROOT->ProcessLine(".L Ne20_gascell.C");
//  gROOT->ProcessLine("ne20_gascell()");
}
