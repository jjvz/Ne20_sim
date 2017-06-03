//=========Display Scattering Geometry for 20Ne gas cell target setup
//=========  JJvZ, May 2017

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include <TH1F.h>
#include <TH2F.h>     
#include <TTree.h>
#include <TFile.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h> 
#include <TObject.h>
#include <TCanvas.h>
#include <TWbox.h>
#include <TEllipse.h>
#include <TLine.h>
#include <TBox.h>

#include "myMCFuncts.h"

#define diaglable "gascell_draw.dat"
#define PRNT_DIAG		//uncomment to print scatter diagram 'online'

FILE *posfile; 

const int NN   = 5500;	// number of event iterations 
//const int PRNT = 10;	// number of scatter events to plot in diagram

// Gas cell geometry (in mm):
const float L = 155.;
const float H = 84.;
const float orgx = 0;
const float orgy = 0;        // position of source (origin, mm)

// Si detector dimentions:
// Dimentions of cell and Si detectors (in mm)
const float SiX0 = 13.;
const float SiL  = 48.;
const float SiH  = 48.;
const float SiT  = 1.;        // active area (mm)
const float SiSepX  = 12.;    // separation of Si active areas (specify in mm)
const float SiSepY1 = 15.;    // y-distance from centre line to si-detectors (specify in mm)
const float SiSepY2 = 19.;    // y-distance from centre line to si-detectors (specify in mm)
const float SiSW    = 3.;     // strip width
const float SiSSep  = 0.133;

// Canvas dimentions:
// Length conversion (factor to convert mm to fraction of Canvas):
const float Lconv = 0.80/L;        // convert to Canvas dimentions: x: 200mm * x = 0.65 -> x = 0.65/200
const float Hconv = 0.80/H;        // convert to Canvas dimentions: y: 50mm * y = 0.35 -> x = 0.35/50

// Gas cell geometry (in fractions of Canvas size):
const float cOrgx = 0.15;
const float cOrgy = 0.5;        // position of source (origin, fraction of Canvas)
float cellL, cellH, cSiL, cSiH, cSiT, cSiSepX, cSiSepY1, cSiSepY2, cSiSW, cSiSSep;       

/*------------ Constants----------------------------------------------*/
const double PI=3.14159265359;

Int_t b_N0=-10, b_Nbreak=-10, b_Si1=0, b_Si2=0, b_Si3=0, b_Si4=0, b_Si5=0, b_Si6=0;
Float_t b_posxb=0, b_posxB=0, b_posyb=0, b_posyB=0;

DetBox Si1box;          // Declare Si1box of type Box
DetBox Si2box;          // Declare Si2box of type Box
DetBox Si3box;          // Declare Si3box of type Box
DetBox Si4box;          // Declare Si4box of type Box
   
// ************************************************************************************

int gascell_draw(Int_t PRNT=100)
{
    Int_t i=1, j=0, p=0,q, Nj, jc;
    Float_t x1=0, y1=0, x2b=0, y2b=0, x2B=0, y2B=0, hoekb=0, hoekB=0;
    std::vector<float> cx1, cy1, cx2b, cy2b, choekb;
    std::vector<int> jcount;
    Char_t  buff[1000], answ='c';
		
    posfile = fopen(diaglable, "r");
    if(posfile){
//    	fgets(buff, 80, posfile);
//    	fgets(buff, 80, posfile);
		loop:
		while(i > NULL){
			i = fscanf(posfile,"%d",&jc);
			i = fscanf(posfile,"%f\t%f",&x1,&y1);
			jcount.push_back(jc);
			cx1.push_back(cOrgx + x1*Lconv);
			cy1.push_back(cOrgy + y1*Hconv);
//			std::cout<<"\n"<<x1<<" (of cell lenght 155mm)\t"<<cx1.at(p)<<" (of canvas l=0,8)";

			for(j=0;j<jcount.at(p);j=j+2) { 	// jcount=4 -> j=0: x0,y1,h -> j=2: x2,y3,h -> j=4: stop
				i = fscanf(posfile,"%f\t%f\t%f",&x2b,&y2b,&hoekb);
				cx2b.push_back(cOrgx + x2b*Lconv);
				cy2b.push_back(cOrgy + y2b*Hconv);
				choekb.push_back(hoekb);
//				std::cout<<"\n"<<x2b<<"\t"<<y2b<<"\t"<<hoekb;
			}
			p++;
			goto loop;
		}
	    std::cout<<"\n "<<cx1.size()<<" events read."<<std::endl;
    	fclose (posfile);
    }
    else {
    	std::cout<< "\nERROR: "<< diaglable <<" file does not exist!!!!"<<std::endl;
    	return 0;
    }

#ifdef PRNT_DIAG
	TCanvas *c1=new TCanvas("c1","Scattering Geometry",900,480,660,350);
#endif

// *******************Drawing gas cell canvas: ******************************************
#ifdef PRNT_DIAG
	cellL = L*Lconv;       // length of gas cell (rel. to canvas)
	cellH = H/2.*Hconv;    // half height of gas cell (rel. to canvas)
	cSiL = SiL*Lconv;
	cSiH = SiH*Lconv;
	cSiT = SiT*Hconv;        // active area (in canvas dim.)
	cSiSepX = SiSepX*Lconv;    // separation of Si active areas (specify in mm)
	cSiSepY1 = SiSepY1*Hconv;   // y-distance from centre line to si-detectors (specify in mm)
	cSiSepY2 = SiSepY2*Hconv;   // y-distance from centre line to si-detectors (specify in mm)
	cSiSW = SiSW*Lconv;             // strip width
	cSiSSep = SiSSep*Lconv;

    c1->cd();
// Draw Source(origin):
    TEllipse *ellipse = new TEllipse(cOrgx,cOrgy,0.005,0.005,0,360,0);
    ellipse->Draw();
    c1->Update();

// Draw Gas Cell: (bottom_left: cOrgx,y0; top_right: x1,y1)
    TWbox *box = new TWbox(cOrgx,cOrgy-cellH,cOrgx+cellL,cOrgy+cellH);
    box->SetLineColor(1);
    box->SetLineStyle(1);
    box->SetLineWidth(2);
    box->SetBorderSize(1);
    box->SetFillColor(19);  
    box->Draw("l");
    
       TLatex *   tex = new TLatex(0.3155704,0.67,"Si 1");
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.63,0.71,"Si 2");
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.32,0.29,"Si 3");
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.63,0.26,"Si 4");
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.48,0.91,"Gas cell");
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.04,0.52,"beam");
   tex->SetLineWidth(2);
   tex->Draw();
   TArrow *arrow = new TArrow(0.026,0.5,0.132,0.5,0.02,">");
   arrow->SetFillColor(1);
   arrow->SetFillStyle(1001);
   arrow->SetLineWidth(3);
   arrow->SetAngle(52);
   arrow->Draw();

    c1->Update();
  
// Draw Centre axis:
    TLine *line = new TLine(0.01,0.5,0.99,0.5);
    line->SetLineStyle(3);
    line->Draw();
    c1->Update();
    
// Si1box specification
	Si1box.sx0 = cOrgx + SiX0*Lconv; 
	Si1box.sy0 = cOrgy+cSiSepY1+cSiT; 
	Si1box.sx1 = Si1box.sx0+cSiL;
	Si1box.sy1 = cOrgy+cSiSepY1;

// Si2box specification
	Si2box.sx0 = Si1box.sx1+cSiSepX; 
	Si2box.sy0 = cOrgy+cSiSepY2+cSiT; 
	Si2box.sx1 = Si2box.sx0+cSiL;
	Si2box.sy1 = cOrgy+cSiSepY2;

// Si3box specification
	Si3box.sx0 = cOrgx + SiX0*Lconv; 
	Si3box.sy0 = cOrgy-cSiSepY1-cSiT; 
	Si3box.sx1 = Si3box.sx0+cSiL;
	Si3box.sy1 = cOrgy-cSiSepY1;

// Si4box specification
	Si4box.sx0 = Si3box.sx1+cSiSepX; 
	Si4box.sy0 = cOrgy-cSiSepY2-cSiT; 
	Si4box.sx1 = Si4box.sx0+cSiL;
	Si4box.sy1 = cOrgy-cSiSepY2;

    TBox *si_det1 = new TBox(Si1box.sx0,Si1box.sy0,Si1box.sx1,Si1box.sy1);
    TBox *si_det2 = new TBox(Si2box.sx0,Si2box.sy0,Si2box.sx1,Si2box.sy1);
    TBox *si_det3 = new TBox(Si3box.sx0,Si3box.sy0,Si3box.sx1,Si3box.sy1);
    TBox *si_det4 = new TBox(Si4box.sx0,Si4box.sy0,Si4box.sx1,Si4box.sy1);
 
    si_det1->SetFillStyle(0);  
    si_det2->SetFillStyle(0);   
    si_det3->SetFillStyle(0);  
    si_det4->SetFillStyle(0);  
    si_det1->Draw();
    si_det3->Draw();
    si_det2->Draw();
    si_det4->Draw();
    c1->Update();
#endif

//**********************************************************************************************************************	
	q=0;
    for(Nj=0; Nj<std::min(cx1.size(),PRNT); Nj++) {	// for every read line:
  		printf("\r...progress: %.2f%%",100.*Nj/std::min(cx1.size(),PRNT));

#ifdef PRNT_DIAG
		TEllipse *ellipse = new TEllipse(cx1.at(Nj),cy1.at(Nj),0.001,0.003,0,360,0);
		ellipse->Draw();

		c1->Modified();
		c1->Update();       // uncomment to see individual events
		if(answ!='R') {  
			std::cout<<"\n.... press any key to STEP, 'R' to RUN all, or 'Q' to QUIT...\t";
			answ = std::cin.get();
			if(answ=='Q') return 0;
		}

		for(j=0; j < jcount.at(Nj); j=j+2) {	// for every (x,y)-set in line Nj:
			line = new TLine(cx1.at(Nj),cy1.at(Nj),cx2b.at(q),cy2b.at(q));
			line->SetLineStyle(1);
			if(j==0) line->SetLineColor(kRed);
			else if(j==2) line->SetLineColor(kBlue); 
			else if(j==4) line->SetLineColor(kGreen);  
			if(jcount.at(Nj)>0) { 	// select which events to draw:
				line->Draw();
				c1->Modified();
				c1->Update();       // uncomment to see individual events
			} 
			q++;
		}
	//		std::cout<<"\nb: "<<choekb.at(Nj)<<" deg\t B: "<<choekB.at(Nj)<<" deg \n";
//		std::cout<<" \t x1: "<<cx2b.at(Nj)<<"\t y1: "<<cy2b.at(Nj)<<"\t x2: "<<cx2B.at(Nj)<<"\t y2: "<<cy2B.at(Nj)<<"\n";
#endif
	}  // end of event
//**********************************************************************************************************************	
	std::cout<<std::endl; 	
	return 0;
}

