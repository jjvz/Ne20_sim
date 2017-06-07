========= Display Scattering Geometry for 20Ne gas cell target setup ========= 
========= JJvZ, Jun 2016 =========

Simulates K600 spectrometer DAQ for 20Ne(a,a') scatter including alpha breakup reactions plus Si-array and NaI anciliry detectors. Code is based on simple monte carlo calculations and creates a ROOT tree file with event branches similar to the K600 DAQ run output.

// Usage:
// To Compile with g++ in Ubuntu14.04 LTS:
// 1) In order to add the new TClass "myDetData" to Root's Dictionary of classes, use "rootint" as follows:
// ~/codes$ rootcint -f ne20_Dict.C -c myDetData.h my_si_LinkDef.h
// The file "my_det_LinkDef.h" is shown below and must be in the same directory. The files "my_det_Dict.C" and "my_det_Dict.h" are created by rootint
// 2) Compile the files into a library:
// ~/codes$ g++ Ne20_gascell.C BREAKUP.C E_OUT.C myMCFuncts.C myDetData.C ne20_Dict.C `root-config --libs --cflags` -o Ne20_gascell
// 3) Run the executable:
// $ ./Ne20_gascell
//
// This following bash script will do all the above:
// source ./run_ne20.sh 

// Running executable in Ubuntu:
// ./Ne20_gascell

// To compile in ROOT with CINT:
// root [0] .L BREAKUP.C
// root [1] .L E_OUT.C
// root [2] .L myMCFuncts.C
// root [3] .L myDetData.C++
// root []  .L MCGAMMA.C
// root [4] .L Ne20_gascell.C
// root [5] ne20_gascell()

To view the drawing of the Si detector box and breakup scattering, run the foillowing in Root:
(make sure the original Ne20_Gascell program was run with the "#define PRNT_DIAG" uncommented, anf the number of events PRNT set to something like 2000)
root [1] .L gascell_draw.C
root [2] gascell_draw(2000)


Typical Draw() calls:
 root [12] DATA->Draw("SiData->HitPos[1]:SiData->HitPos[0]>>h04(500,0,160,500,-35,35)","","col")	// all hits at Si faces, actually hitting the Si detectors
 root [10] DATA->Draw("pos1xy[1]:pos1xy[0]>>h03(500,0,160,500,-35,35)","","colsame")
