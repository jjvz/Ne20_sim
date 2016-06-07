# Ne20_sim
Simulates K600 spectrometer DAQ for 20Ne(a,a') scatter including alpha breakup reactions plus Si-array and NaI anciliry detectors. Code is based on simple monte carlo calculations and creates a ROOT tree file with event branches similar to the K600 DAQ run output.

Things to do:
1) The file "lyne-ne20.txt" contains the info for simulating the gamma decays into the large volume NaI detector HAGAR. User can change items herein.

Usage:
 root[] .L Ne20_gascell.C
 root[] ne20_gascell(NN)      // default: NN=20000 
        (where NN is number of repeats / incident projectiles, e.g. ne20_gascell(20000)) 
        


========= Display Scattering Geometry for 20Ne gas cell target setup =========
========= JJvZ, Jun 2016 =========
