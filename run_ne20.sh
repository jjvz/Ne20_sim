#!/bin/bash
# usage: source ./run_ne20.sh

rootcint -f ne20_Dict.C -c myDetData.h my_si_LinkDef.h
echo "Root Dictionary added"
g++ -g Ne20_gascell.C BREAKUP.C E_OUT.C MCGAMMA.C myMCFuncts.C myDetData.C ne20_Dict.C `root-config --libs --cflags --glibs` -o Ne20_gascell
echo "Code compiled with g++ - Now running executable:"
./Ne20_gascell

