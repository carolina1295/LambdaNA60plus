#include <Riostream.h>
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TMath.h"
#include "TClonesArray.h"
#include "Riostream.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH1F.h"
#include "TAxis.h"
#include "TString.h"
#include "TCanvas.h"
#include "TH3D.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include <string>
#include <vector>
#include <TNtuple.h>
#include "TH3F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TMultiGraph.h"

void efficiency(TH1F* hRec, TH1F* hMC, TH1F* hGen, char* xaxis, TFile *f);
void efficiencyGen(TH3F *hMCGenRec[3],char* title, TFile* f);

void EfficiencyLambda(){
    TFile *fLambda = new TFile("LambdastandconBKG-Signal-histos.root");
	TFile *fLambdaBar = new TFile("LambdaBarstandconBKG-Signal-histos.root");
    TFile *fLambdaMoreLayer = new TFile("LambdaMoreLayerconBKG-Signal-histos.root");
	TFile *fLambdaBarMoreLayer = new TFile("LambdaBarMoreLayerconBKG-Signal-histos.root");
	TFile *fLambdaMoreLayer5NoBKG = new TFile("LambdaMoreLayer5NoBKG-Signal-histos.root");
	TFile *fLambdaBarMoreLayer5NoBKG = new TFile("LambdaBarMoreLayer5NoBKG-Signal-histos.root");
	TFile *fLambdaMoreLayer7NoBKG = new TFile("LambdaMoreLayer7NoBKG-Signal-histos.root");
	TFile *fLambdaBarMoreLayer7NoBKG = new TFile("LambdaBarMoreLayer7NoBKG-Signal-histos.root");
	TFile *fLambdaMoreLayer10NoBKG = new TFile("LambdaMoreLayer10NoBKG-Signal-histos.root");
	TFile *fLambdaBarMoreLayer10NoBKG = new TFile("LambdaBarMoreLayer10NoBKG-Signal-histos.root");
	TFile *fPion = new TFile("PION-Signal-histos.root");
	TFile *fPionpos = new TFile("PIONpos-Signal-histos.root");
	TFile *fProt = new TFile("PROTON-Signal-histos.root");
	TFile *fAntiProt = new TFile("AntiPROTON-Signal-histos.root");
	TFile *fD0 = new TFile("D0senzaBKG-Signal-histos.root");
	TFile *fD0MoreLayer = new TFile("D0MoreLayer-Signal-histos.root");

	TFile *foutLambda=new TFile("Efficiency.root","RECREATE");
    //Lambda
    //Distribuzione 3D per efficienza
	TH3F* hYPtLzMC[2], *hYPtLzGen[2], *hYPtLzRec[2] ;
	 TH3F* hYPtLzGen0fake[2], *hYPtLzRec0fake[2] ;
	
	hYPtLzMC[0]=(TH3F*)fLambda->Get("hYPtLzMC");
	hYPtLzGen[0]=(TH3F*)fLambda->Get("hYPtLzGen");
	hYPtLzRec[0]=(TH3F*)fLambda->Get("hYPtLzRec");
	hYPtLzGen0fake[0]=(TH3F*)fLambda->Get("hYPtLzGen0fake");
	hYPtLzRec0fake[0]=(TH3F*)fLambda->Get("hYPtLzRec0fake");

	hYPtLzMC[1]=(TH3F*)fLambdaMoreLayer->Get("hYPtLzMC");
	hYPtLzGen[1]=(TH3F*)fLambdaMoreLayer->Get("hYPtLzGen");
	hYPtLzRec[1]=(TH3F*)fLambdaMoreLayer->Get("hYPtLzRec");
	hYPtLzGen0fake[1]=(TH3F*)fLambdaMoreLayer->Get("hYPtLzGen0fake");
	hYPtLzRec0fake[1]=(TH3F*)fLambdaMoreLayer->Get("hYPtLzRec0fake");

	TH1F* hYMC[2],*hYGen[2], *hYRec[2], *hPtMC[2], *hPtGen[2], *hPtRec[2], *hLzMC[2], *hLzGen[2], *hLzRec[2];
	TH1F* hYGen0fake[2], *hYRec0fake[2], *hPtGen0fake[2], *hPtRec0fake[2], *hLzGen0fake[2], *hLzRec0fake[2];

	for(int i=0;i<2;i++){
		//efficiency for y
		hYMC[i]=(TH1F*)hYPtLzMC[i]->ProjectionX();
		hYGen[i]=(TH1F*)hYPtLzGen[i]->ProjectionX();
		hYRec[i]=(TH1F*)hYPtLzRec[i]->ProjectionX();
		hYGen0fake[i]=(TH1F*)hYPtLzGen0fake[i]->ProjectionX();
		hYRec0fake[i]=(TH1F*)hYPtLzRec0fake[i]->ProjectionX();

		//efficiency for Pt 	
		hPtMC[i]=(TH1F*)hYPtLzMC[i]->ProjectionY();
		hPtGen[i]=(TH1F*)hYPtLzGen[i]->ProjectionY();
		hPtRec[i]=(TH1F*)hYPtLzRec[i]->ProjectionY();
		hPtGen0fake[i]=(TH1F*)hYPtLzGen0fake[i]->ProjectionY();
		hPtRec0fake[i]=(TH1F*)hYPtLzRec0fake[i]->ProjectionY();
		
		//efficiency for Lz
		hLzMC[i]=(TH1F*)hYPtLzMC[i]->ProjectionZ();
		hLzGen[i]=(TH1F*)hYPtLzGen[i]->ProjectionZ();
		hLzRec[i]=(TH1F*)hYPtLzRec[i]->ProjectionZ();
		hLzGen0fake[i]=(TH1F*)hYPtLzGen0fake[i]->ProjectionZ();
		hLzRec0fake[i]=(TH1F*)hYPtLzRec0fake[i]->ProjectionZ();
		
		if(i==0){
			efficiency(hYRec[i],hYMC[i],hYGen[i],"y #Lambda",foutLambda);
			efficiency(hPtRec[i],hPtMC[i],hPtGen[i],"P_{t} #Lambda",foutLambda);
			efficiency(hLzRec[i],hLzMC[i],hLzGen[i],"Lz #Lambda",foutLambda);
				
			efficiency(hYRec0fake[i],hYMC[i],hYGen0fake[i],"y #Lambda 0fake",foutLambda);
			efficiency(hPtRec0fake[i],hPtMC[i],hPtGen0fake[i],"P_{t} #Lambda 0fake",foutLambda);
			efficiency(hLzRec0fake[i],hLzMC[i],hLzGen0fake[i],"Lz #Lambda 0fake",foutLambda);
		}

		if(i==1){
			efficiency(hYRec[i],hYMC[i],hYGen[i],"y #Lambda 10 layers",foutLambda);
			efficiency(hPtRec[i],hPtMC[i],hPtGen[i],"P_{t} #Lambda 10 layers",foutLambda);
			efficiency(hLzRec[i],hLzMC[i],hLzGen[i],"Lz #Lambda 10 layers",foutLambda);
			
			efficiency(hYRec0fake[i],hYMC[i],hYGen0fake[i],"y #Lambda 10 layers 0fake",foutLambda);
			efficiency(hPtRec0fake[i],hPtMC[i],hPtGen0fake[i],"P_{t} #Lambda 10 layers 0fake",foutLambda);
			efficiency(hLzRec0fake[i],hLzMC[i],hLzGen0fake[i],"Lz #Lambda 10 layers 0fake",foutLambda);
		}
	}

	TH3F* hYPtLz5NoBkg[3], *hYPtLz7NoBkg[3], *hYPtLz10NoBkg[3];
	hYPtLz5NoBkg[0]=(TH3F*)fLambdaMoreLayer5NoBKG->Get("hYPtLzMC");
	hYPtLz5NoBkg[1]=(TH3F*)fLambdaMoreLayer5NoBKG->Get("hYPtLzGen");
	hYPtLz5NoBkg[2]=(TH3F*)fLambdaMoreLayer5NoBKG->Get("hYPtLzRec");
	efficiencyGen(hYPtLz5NoBkg,"for minimum 5 hits", foutLambda);

	hYPtLz7NoBkg[0]=(TH3F*)fLambdaMoreLayer7NoBKG->Get("hYPtLzMC");
	hYPtLz7NoBkg[1]=(TH3F*)fLambdaMoreLayer7NoBKG->Get("hYPtLzGen");
	hYPtLz7NoBkg[2]=(TH3F*)fLambdaMoreLayer7NoBKG->Get("hYPtLzRec");
	efficiencyGen(hYPtLz7NoBkg,"for minimum 7 hits", foutLambda);

	hYPtLz10NoBkg[0]=(TH3F*)fLambdaMoreLayer10NoBKG->Get("hYPtLzMC");
	hYPtLz10NoBkg[1]=(TH3F*)fLambdaMoreLayer10NoBKG->Get("hYPtLzGen");
	hYPtLz10NoBkg[2]=(TH3F*)fLambdaMoreLayer10NoBKG->Get("hYPtLzRec");
	efficiencyGen(hYPtLz10NoBkg,"for minimum 10 hits", foutLambda);

	TH3F* hYPtLz5NoBkgBar[3], *hYPtLz7NoBkgBar[3], *hYPtLz10NoBkgBar[3];
	hYPtLz5NoBkgBar[0]=(TH3F*)fLambdaBarMoreLayer5NoBKG->Get("hYPtLzMC");
	hYPtLz5NoBkgBar[1]=(TH3F*)fLambdaBarMoreLayer5NoBKG->Get("hYPtLzGen");
	hYPtLz5NoBkgBar[2]=(TH3F*)fLambdaBarMoreLayer5NoBKG->Get("hYPtLzRec");
	efficiencyGen(hYPtLz5NoBkgBar,"for minimum 5 hits Bar", foutLambda);

	hYPtLz7NoBkgBar[0]=(TH3F*)fLambdaBarMoreLayer7NoBKG->Get("hYPtLzMC");
	hYPtLz7NoBkgBar[1]=(TH3F*)fLambdaBarMoreLayer7NoBKG->Get("hYPtLzGen");
	hYPtLz7NoBkgBar[2]=(TH3F*)fLambdaBarMoreLayer7NoBKG->Get("hYPtLzRec");
	efficiencyGen(hYPtLz7NoBkgBar,"for minimum 7 hits Bar", foutLambda);

	hYPtLz10NoBkgBar[0]=(TH3F*)fLambdaBarMoreLayer10NoBKG->Get("hYPtLzMC");
	hYPtLz10NoBkgBar[1]=(TH3F*)fLambdaBarMoreLayer10NoBKG->Get("hYPtLzGen");
	hYPtLz10NoBkgBar[2]=(TH3F*)fLambdaBarMoreLayer10NoBKG->Get("hYPtLzRec");
	efficiencyGen(hYPtLz10NoBkgBar,"for minimum 10 hits Bar", foutLambda);

	TH1F* h5eff[3], *h7eff[3], *h10eff[3], *h5effBar[3], *h7effBar[3], *h10effBar[3];
	for(int j=0;j<3;j++){
		if(j==0){
		 	h5eff[j]=(TH1F*)foutLambda->Get("Gen over GenMC for y for minimum 5 hits");
			h7eff[j]=(TH1F*)foutLambda->Get("Gen over GenMC for y for minimum 7 hits");
			h10eff[j]=(TH1F*)foutLambda->Get("Gen over GenMC for y for minimum 10 hits");
			h5effBar[j]=(TH1F*)foutLambda->Get("Gen over GenMC for y for minimum 5 hits Bar");
			h7effBar[j]=(TH1F*)foutLambda->Get("Gen over GenMC for y for minimum 7 hits Bar");
			h10effBar[j]=(TH1F*)foutLambda->Get("Gen over GenMC for y for minimum 10 hits Bar");
		}
		if(j==1){
            h5eff[j]=(TH1F*)foutLambda->Get("Gen over GenMC for Pt for minimum 5 hits");
			h7eff[j]=(TH1F*)foutLambda->Get("Gen over GenMC for Pt for minimum 7 hits");
			h10eff[j]=(TH1F*)foutLambda->Get("Gen over GenMC for Pt for minimum 10 hits");
			h5effBar[j]=(TH1F*)foutLambda->Get("Gen over GenMC for Pt for minimum 5 hits Bar");
			h7effBar[j]=(TH1F*)foutLambda->Get("Gen over GenMC for Pt for minimum 7 hits Bar");
			h10effBar[j]=(TH1F*)foutLambda->Get("Gen over GenMC for Pt for minimum 10 hits Bar");
		}
		if(j==2){
            h5eff[j]=(TH1F*)foutLambda->Get("Gen over GenMC for Lz for minimum 5 hits");
			h7eff[j]=(TH1F*)foutLambda->Get("Gen over GenMC for Lz for minimum 7 hits");
			h10eff[j]=(TH1F*)foutLambda->Get("Gen over GenMC for Lz for minimum 10 hits");
			h5effBar[j]=(TH1F*)foutLambda->Get("Gen over GenMC for Lz for minimum 5 hits Bar");
			h7effBar[j]=(TH1F*)foutLambda->Get("Gen over GenMC for Lz for minimum 7 hits Bar");
			h10effBar[j]=(TH1F*)foutLambda->Get("Gen over GenMC for Lz for minimum 10 hits Bar");
      	}
	}


	TCanvas *cPoint = new TCanvas("cPoint");
	TCanvas *cPointBar = new TCanvas("cPointBar");
	cPoint->Divide(3,1);
	cPointBar->Divide(3,1);
	for(int i=0;i<3;i++){
		cPoint->cd(i+1);
		gPad->SetGrid();
		h5eff[i]->SetDirectory(0);
		h7eff[i]->SetDirectory(0);
		h10eff[i]->SetDirectory(0);
		h5eff[i]->SetLineColor(2);
		h7eff[i]->SetLineColor(4);
		h10eff[i]->SetLineColor(93);
		h5eff[i]->SetMarkerStyle(20);
		h5eff[i]->SetMarkerColor(2);
		h7eff[i]->SetMarkerStyle(22);
		h7eff[i]->SetMarkerColor(4);
		h10eff[i]->SetMarkerStyle(21);
		h10eff[i]->SetMarkerColor(93);
		h5eff[i]->Draw();
		h7eff[i]->Draw("SAME");
		h10eff[i]->Draw("SAME");
		cPointBar->cd(i+1);
		gPad->SetGrid();
		h5effBar[i]->SetLineColor(2);
		h7effBar[i]->SetLineColor(4);
		h10effBar[i]->SetLineColor(93);
		h5effBar[i]->SetDirectory(0);
		h7effBar[i]->SetDirectory(0);
		h10effBar[i]->SetDirectory(0);
		h5effBar[i]->SetMarkerStyle(20);
		h5effBar[i]->SetMarkerColor(2);
		h7effBar[i]->SetMarkerStyle(22);
		h7effBar[i]->SetMarkerColor(4);
		h10effBar[i]->SetMarkerStyle(21);
		h10effBar[i]->SetMarkerColor(93);
		h5effBar[i]->Draw();
		h7effBar[i]->Draw("SAME");
		h10effBar[i]->Draw("SAME");
	}
	cPoint->cd(1);
	auto legendP = new TLegend(0.1,0.7,0.48,0.9);
   	//legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
   	legendP->AddEntry(h5eff[0],"#Lambda efficiency for 5 min hits","lep");
  	legendP->AddEntry(h7eff[0],"#Lambda efficiency for 7 min hits","lep");
   	legendP->AddEntry(h10eff[0],"#Lambda efficiency for 10 min hits","lep");	   
   	legendP->Draw();

	cPointBar->cd(1);
	auto legendPBar = new TLegend(0.1,0.7,0.48,0.9);
	legendPBar->AddEntry(h5effBar[0],"#LambdaBar efficiency for 5 min hits","lep");
  	legendPBar->AddEntry(h7effBar[0],"#LambdaBar efficiency for 7 min hits","lep");
   	legendPBar->AddEntry(h10effBar[0],"#LambdaBar efficiency for 10 min hits","lep");	   
   	legendPBar->Draw();
	h5eff[0]->SetTitle("Efficiency on y");
	h5eff[1]->SetTitle("Efficiency on Pt");
	h5eff[2]->SetTitle("Efficiency on Lz");
	h5effBar[0]->SetTitle("Efficiency on y");
	h5effBar[1]->SetTitle("Efficiency on Pt");
	h5effBar[2]->SetTitle("Efficiency on Lz");

	


	//Distribuzione 3D per efficienza lambda bar
	TH3F* hYPtLzMCBar[2], *hYPtLzGenBar[2], *hYPtLzRecBar[2];
	TH3F* hYPtLzGenBar0fake[2], *hYPtLzRecBar0fake[2];
	hYPtLzMCBar[0]=(TH3F*)fLambdaBar->Get("hYPtLzMC");
	hYPtLzGenBar[0]=(TH3F*)fLambdaBar->Get("hYPtLzGen");
	hYPtLzRecBar[0]=(TH3F*)fLambdaBar->Get("hYPtLzRec");
	hYPtLzGenBar0fake[0]=(TH3F*)fLambda->Get("hYPtLzGen0fake");
	hYPtLzRecBar0fake[0]=(TH3F*)fLambda->Get("hYPtLzRec0fake");

	hYPtLzMCBar[1]=(TH3F*)fLambdaBarMoreLayer->Get("hYPtLzMC");
	hYPtLzGenBar[1]=(TH3F*)fLambdaBarMoreLayer->Get("hYPtLzGen");
	hYPtLzRecBar[1]=(TH3F*)fLambdaBarMoreLayer->Get("hYPtLzRec");
	hYPtLzGenBar0fake[1]=(TH3F*)fLambdaMoreLayer->Get("hYPtLzGen0fake");
	hYPtLzRecBar0fake[1]=(TH3F*)fLambdaMoreLayer->Get("hYPtLzRec0fake");


	TH1F* hYMCBar[2],*hYGenBar[2], *hYRecBar[2], *hPtMCBar[2], *hPtGenBar[2], *hPtRecBar[2], *hLzMCBar[2], *hLzGenBar[2], *hLzRecBar[2];
	TH1F* hYGenBar0fake[2], *hYRecBar0fake[2], *hPtGenBar0fake[2], *hPtRecBar0fake[2], *hLzGenBar0fake[2], *hLzRecBar0fake[2];
	//efficiency for y
	for(int i=0;i<2;i++){
		hYMCBar[i]=(TH1F*)hYPtLzMCBar[i]->ProjectionX();
		hYGenBar[i]=(TH1F*)hYPtLzGenBar[i]->ProjectionX();
		hYRecBar[i]=(TH1F*)hYPtLzRecBar[i]->ProjectionX();

		hYGenBar0fake[i]=(TH1F*)hYPtLzGenBar0fake[i]->ProjectionX();
		hYRecBar0fake[i]=(TH1F*)hYPtLzRecBar0fake[i]->ProjectionX();
		
		//efficiency for Pt 	
		hPtMCBar[i]=(TH1F*)hYPtLzMCBar[i]->ProjectionY();
		hPtGenBar[i]=(TH1F*)hYPtLzGenBar[i]->ProjectionY();
		hPtRecBar[i]=(TH1F*)hYPtLzRecBar[i]->ProjectionY();

		hPtGenBar0fake[i]=(TH1F*)hYPtLzGenBar0fake[i]->ProjectionY();
		hPtRecBar0fake[i]=(TH1F*)hYPtLzRecBar0fake[i]->ProjectionY();
		
		//efficiency for Lz
		hLzMCBar[i]=(TH1F*)hYPtLzMCBar[i]->ProjectionZ();
		hLzGenBar[i]=(TH1F*)hYPtLzGenBar[i]->ProjectionZ();
		hLzRecBar[i]=(TH1F*)hYPtLzRecBar[i]->ProjectionZ();

		hLzGenBar0fake[i]=(TH1F*)hYPtLzGenBar0fake[i]->ProjectionZ();
		hLzRecBar0fake[i]=(TH1F*)hYPtLzRecBar0fake[i]->ProjectionZ();

		
		if(i==0){
			efficiency(hYRecBar[i],hYMCBar[i],hYGenBar[i],"y #Lambda Bar",foutLambda);
			efficiency(hPtRecBar[i],hPtMCBar[i],hPtGenBar[i],"P_{t} #Lambda Bar",foutLambda);
			efficiency(hLzRecBar[i],hLzMCBar[i],hLzGenBar[i],"Lz #Lambda Bar",foutLambda);

			efficiency(hYRecBar0fake[i],hYMC[i],hYGenBar0fake[i],"y #Lambda Bar 0fake",foutLambda);
			efficiency(hPtRecBar0fake[i],hPtMC[i],hPtGenBar0fake[i],"P_{t} #Lambda Bar 0fake",foutLambda);
			efficiency(hLzRecBar0fake[i],hLzMC[i],hLzGenBar0fake[i],"Lz #Lambda Bar 0fake",foutLambda);
			
		}

		if(i==1){
			efficiency(hYRecBar[i],hYMCBar[i],hYGenBar[i],"y #Lambda Bar 10 layers",foutLambda);
			efficiency(hPtRecBar[i],hPtMCBar[i],hPtGenBar[i],"P_{t} #Lambda Bar 10 layers",foutLambda);
			efficiency(hLzRecBar[i],hLzMCBar[i],hLzGenBar[i],"Lz #Lambda Bar 10 layers",foutLambda);

			efficiency(hYRecBar0fake[i],hYMC[i],hYGenBar0fake[i],"y #Lambda Bar 10 layers 0fake",foutLambda);
			efficiency(hPtRecBar0fake[i],hPtMC[i],hPtGenBar0fake[i],"P_{t} #Lambda Bar 10 layers 0fake",foutLambda);
			efficiency(hLzRecBar0fake[i],hLzMC[i],hLzGenBar0fake[i],"Lz #Lambda Bar 10 layers 0fake",foutLambda);
		}
	}

	TH1F* hPtGenD0=(TH1F*)fD0->Get("hPtGen");
   	TH1F* hPtRecoAllD0=(TH1F*)fD0->Get("hPtRecoAll");
  	TH1F* hPtGenRecoAllD0=(TH1F*)fD0->Get("hPtGenRecoAll");
	TH1F* hPtGenMLD0=(TH1F*)fD0MoreLayer->Get("hPtGen");
   	TH1F* hPtRecoAllMLD0=(TH1F*)fD0MoreLayer->Get("hPtRecoAll");
  	TH1F* hPtGenRecoAllMLD0=(TH1F*)fD0MoreLayer->Get("hPtGenRecoAll");
	TH1F* hYGenD0=(TH1F*)fD0->Get("hYGen");
   	TH1F* hYRecoAllD0=(TH1F*)fD0->Get("hYRecoAll");
  	TH1F* hYGenRecoAllD0=(TH1F*)fD0->Get("hYGenRecoAll");
	TH1F* hYGenMLD0=(TH1F*)fD0MoreLayer->Get("hYGen");
   	TH1F* hYRecoAllMLD0=(TH1F*)fD0MoreLayer->Get("hYRecoAll");
  	TH1F* hYGenRecoAllMLD0=(TH1F*)fD0MoreLayer->Get("hYGenRecoAll");
	efficiency(hPtRecoAllD0,hPtGenD0,hPtGenRecoAllD0, "p_{t} D_{0}",foutLambda);
	efficiency(hYRecoAllD0,hYGenD0,hYGenRecoAllD0, "y D_{0}",foutLambda);
	efficiency(hPtRecoAllMLD0,hPtGenMLD0,hPtGenRecoAllMLD0, "p_{t} D_{0} 10 layers",foutLambda);
	efficiency(hYRecoAllMLD0,hYGenMLD0,hYGenRecoAllMLD0,"y D_{0} 10 layers",foutLambda);    

	TCanvas *effD0_c = new TCanvas("effD0_c");

	TH1F *effD0[2], *effD010lay[2] ;

	effD0[0]=(TH1F*)foutLambda->Get("Gen over GenMC for y D_{0}");
	effD0[1]=(TH1F*)foutLambda->Get("Gen over GenMC for p_{t} D_{0}");
	effD010lay[0]=(TH1F*)foutLambda->Get("Gen over GenMC for y D_{0} 10 layers");
	effD010lay[1]=(TH1F*)foutLambda->Get("Gen over GenMC for p_{t} D_{0} 10 layers");

	effD0_c->Divide(2,1);
	for(int i=0;i<2;i++){
		
		effD0[i]->SetDirectory(0);
		effD0[i]->SetMarkerStyle(20);
		effD0[i]->SetMarkerSize(0.8);
		effD0[i]->SetMarkerColor(1);
		effD010lay[i]->SetDirectory(0);
		effD010lay[i]->SetMarkerStyle(20);
		effD010lay[i]->SetMarkerSize(0.8);
		effD010lay[i]->SetMarkerColor(2);
		effD0_c->cd(i+1);
		gPad->SetGrid();
		effD0[i]->Draw();
		effD010lay[i]->Draw("SAME");
	}	

	TCanvas *effLambda = new TCanvas("effLambda");
	TCanvas *effLambdaBar = new TCanvas("effLambdaBar");
	TCanvas *effLambda0fake = new TCanvas("effLambda0fake");
	TCanvas *effLambdaBar0fake = new TCanvas("effLambdaBar0fake");
	
	
	TH1F* eff[3], *effBar[3], *eff10lay[3], *effBar10lay[3], *eff10lay0fake[3], *effBar10lay0fake[3];
	eff[0]=(TH1F*)foutLambda->Get("Gen over GenMC for y #Lambda");
	eff[1]=(TH1F*)foutLambda->Get("Gen over GenMC for P_{t} #Lambda");	
	eff[2]=(TH1F*)foutLambda->Get("Gen over GenMC for Lz #Lambda");
	effBar[0]=(TH1F*)foutLambda->Get("Gen over GenMC for y #Lambda Bar");
	effBar[1]=(TH1F*)foutLambda->Get("Gen over GenMC for P_{t} #Lambda Bar");	
	effBar[2]=(TH1F*)foutLambda->Get("Gen over GenMC for Lz #Lambda Bar");

	eff10lay[0]=(TH1F*)foutLambda->Get("Gen over GenMC for y #Lambda 10 layers");
	eff10lay[1]=(TH1F*)foutLambda->Get("Gen over GenMC for P_{t} #Lambda 10 layers");	
	eff10lay[2]=(TH1F*)foutLambda->Get("Gen over GenMC for Lz #Lambda 10 layers");
	effBar10lay[0]=(TH1F*)foutLambda->Get("Gen over GenMC for y #Lambda Bar 10 layers");
	effBar10lay[1]=(TH1F*)foutLambda->Get("Gen over GenMC for P_{t} #Lambda Bar 10 layers");	
	effBar10lay[2]=(TH1F*)foutLambda->Get("Gen over GenMC for Lz #Lambda Bar 10 layers");

	eff10lay0fake[0]=(TH1F*)foutLambda->Get("Gen over GenMC for y #Lambda 10 layers 0fake");
	eff10lay0fake[1]=(TH1F*)foutLambda->Get("Gen over GenMC for P_{t} #Lambda 10 layers 0fake");	
	eff10lay0fake[2]=(TH1F*)foutLambda->Get("Gen over GenMC for Lz #Lambda 10 layers 0fake");
	effBar10lay0fake[0]=(TH1F*)foutLambda->Get("Gen over GenMC for y #Lambda Bar 10 layers 0fake");
	effBar10lay0fake[1]=(TH1F*)foutLambda->Get("Gen over GenMC for P_{t} #Lambda Bar 10 layers 0fake");	
	effBar10lay0fake[2]=(TH1F*)foutLambda->Get("Gen over GenMC for Lz #Lambda Bar 10 layers 0fake");



	effLambda->Divide(3,1);
	effLambdaBar->Divide(3,1);
	for(int i=0;i<3;i++){
		eff[i]->SetDirectory(0);
	eff10lay[i]->SetDirectory(0);
		eff[i]->SetMarkerStyle(20);
		eff10lay[i]->SetMarkerStyle(21);
		eff[i]->SetMarkerSize(0.8);
eff10lay[i]->SetMarkerSize(0.8);
		eff[i]->SetMarkerColor(2);
		eff10lay[i]->SetMarkerColor(4);
		eff[i]->SetLineColor(2);
		eff10lay[i]->SetLineColor(4);
		effLambda->cd(i+1);
		gPad->SetGrid();
		eff[i]->Draw();
		eff10lay[i]->Draw("SAME");
		

		
				effBar[i]->SetDirectory(0);
		effBar10lay[i]->SetDirectory(0);
		effBar[i]->SetMarkerStyle(20);
		effBar[i]->SetMarkerSize(0.8);
		effBar10lay[i]->SetMarkerStyle(20);
		effBar10lay[i]->SetMarkerSize(0.8);
		effBar[i]->SetMarkerColor(2);
		effBar10lay[i]->SetMarkerColor(4);
		effBar[i]->SetLineColor(2);
		effBar10lay[i]->SetLineColor(4);

		effLambdaBar->cd(i+1);
		gPad->SetGrid();
		effBar[i]->Draw();
		effBar10lay[i]->Draw("SAME");
	}
	effLambda->cd(1);
	auto legendL = new TLegend(0.1,0.7,0.48,0.9);
   	//legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
   	legendL->AddEntry(eff[0],"#Lambda efficiency for 5 layers","lep");
  	legendL->AddEntry(eff10lay[0],"#Lambda efficiency for 10 layers","lep");
   	legendL->Draw();
	eff[0]->SetTitle("Efficiency on y");
	eff[1]->SetTitle("Efficiency on Pt");
	eff[2]->SetTitle("Efficiency on Lz");
	effBar[0]->SetTitle("Efficiency on y");
	effBar[1]->SetTitle("Efficiency on Pt");
	effBar[2]->SetTitle("Efficiency on Lz");

	effLambdaBar->cd(1);
	auto legendLB = new TLegend(0.1,0.7,0.48,0.9);
   	//legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
   	legendLB->AddEntry(effBar[0],"#bar{#Lambda} efficiency for 5 layers","lep");
  	legendLB->AddEntry(effBar10lay[0],"#bar{#Lambda} efficiency for 10 layers","lep");
   	legendLB->Draw();

	//Effetto delle fake

	effLambda0fake->Divide(3,1);
	effLambdaBar0fake->Divide(3,1);
	for(int i=0;i<3;i++){

		eff10lay[i]->SetDirectory(0);
		eff10lay0fake[i]->SetDirectory(0);
		eff10lay[i]->SetMarkerStyle(20);
		eff10lay0fake[i]->SetMarkerStyle(20);
		eff10lay[i]->SetMarkerSize(0.8);
		eff10lay0fake[i]->SetMarkerSize(0.8);
		eff10lay[i]->SetMarkerColor(1);
		eff10lay0fake[i]->SetMarkerColor(2);

		effLambda0fake->cd(i+1);
		gPad->SetGrid();
		eff10lay[i]->Draw();
		eff10lay0fake[i]->Draw("SAME");

		effBar10lay[i]->SetDirectory(0);
		effBar10lay0fake[i]->SetDirectory(0);
		effBar10lay[i]->SetMarkerStyle(20);
		effBar10lay0fake[i]->SetMarkerStyle(20);
		effBar10lay[i]->SetMarkerSize(0.8);
		effBar10lay[i]->SetMarkerColor(1);
		effBar10lay0fake[i]->SetMarkerSize(0.8);
		effBar10lay0fake[i]->SetMarkerColor(2);

		effLambdaBar0fake->cd(i+1);
		gPad->SetGrid();
		effBar10lay[i]->Draw();
		effBar10lay0fake[i]->Draw("SAME");
	}
	effLambda0fake->cd(1);
	auto legendL0fake = new TLegend(0.1,0.7,0.48,0.9);
   	//legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
   	legendL0fake->AddEntry(eff10lay[0],"#Lambda efficiency for 10 layers","lep");
  	legendL0fake->AddEntry(eff10lay0fake[0],"#Lambda efficiency for 10 layers 0fake","lep");
   	legendL0fake->Draw();
	eff10lay[0]->SetTitle("Efficiency on y");
	eff10lay[1]->SetTitle("Efficiency on Pt");
	eff10lay[2]->SetTitle("Efficiency on Lz");
	effBar10lay[0]->SetTitle("Efficiency on y");
	effBar10lay[1]->SetTitle("Efficiency on Pt");
	effBar10lay[2]->SetTitle("Efficiency on Lz");

	effLambdaBar0fake->cd(1);
	auto legendLB0fake = new TLegend(0.1,0.7,0.48,0.9);
   	//legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
   	legendLB0fake->AddEntry(effBar10lay[0],"#bar{#Lambda} efficiency for 10 layers","lep");
  	legendLB0fake->AddEntry(effBar10lay0fake[0],"#bar{#Lambda} efficiency for 10 layers 0fake","lep");
   	legendLB0fake->Draw();

	

/*	TH1D* hPion, *hPionpos, *hProt, *hAntiProt;
	hPion = (TH1D*)fPion->Get("hEfficiency");
	hPionpos = (TH1D*)fPionpos->Get("hEfficiency");
	hProt = (TH1D*)fProt->Get("hEfficiency");
	hAntiProt = (TH1D*)fAntiProt->Get("hEfficiency");
	TCanvas *EffPurePart = new TCanvas("EffPurePart");
	gPad->SetGrid();
	hPion->SetMarkerColor(1);
	hPion->SetTitle("Efficiency #frac{n_{rec}}{n_{gen}}");

	hPionpos->SetMarkerColor(2);
	hProt->SetMarkerColor(4);
	hAntiProt->SetMarkerColor(6);
	hPion->SetLineWidth(0);
	hPionpos->SetLineWidth(0);
	hProt->SetLineWidth(0);
	hAntiProt->SetLineWidth(0);

	hPion->Draw();
	hPionpos->Draw("SAME");
	hProt->Draw("Same");
	hAntiProt->Draw("SAME");
	auto legend = new TLegend(0.1,0.7,0.48,0.9);
   	//legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
   	legend->AddEntry(hPion,"#pi- efficiency","lep");
  	legend->AddEntry(hPionpos,"#pi+ efficiency","lep");
   	legend->AddEntry(hProt,"p efficiency","lep");
   	legend->AddEntry(hAntiProt,"#bar{p} efficiency","lep");	   
   	legend->Draw();
	*/

	foutLambda->Close();
}

void efficiency(TH1F* hRec, TH1F* hMC, TH1F* hGen, char* xaxis, TFile* f){
	TH1F *Rec_Gen = (TH1F*)hRec->Clone("Rec_Gen");
	Rec_Gen->SetDirectory(0);
	Rec_Gen->Divide(hRec,hMC,1,1,"B");
	//TCanvas *eff_y=new TCanvas(Form("eff_%s",xaxis));
	Rec_Gen->SetTitle(Form("Rec over GenMC for %s",xaxis));
	Rec_Gen->GetXaxis()->SetTitle(Form("%s",xaxis));
	Rec_Gen->GetYaxis()->SetTitle("Efficiency #epsilon");
	//Rec_Gen->Draw("E1");
	TH1F *Gen_GenMC = (TH1F*)hGen->Clone("Gen_GenMC");
	Gen_GenMC->SetDirectory(0);
	Gen_GenMC->Divide(hGen,hMC,1,1,"B");
	//TCanvas *effGen=new TCanvas(Form("effGen%s",xaxis));
	Gen_GenMC->SetTitle(Form("Gen over GenMC for %s",xaxis));
	Gen_GenMC->GetXaxis()->SetTitle(Form("%s",xaxis));
	Gen_GenMC->GetYaxis()->SetTitle("Efficiency #epsilon");
	//Gen_GenMC->Draw("E1");
	f->cd();
	Gen_GenMC->SetName(Form("Gen over GenMC for %s",xaxis));
	Gen_GenMC->Write();
	Rec_Gen->SetName(Form("Rec over GenMC for %s",xaxis));
	Rec_Gen->Write();

}
void efficiencyGen(TH3F *hMCGenRec[3],char* title, TFile* f){
	TH1F* hY[3],* hPt[3], *hLz[3];//0=MC,1=Gen,2=Rec
	for(int i=0;i<3;i++){
		hY[i]=(TH1F*)hMCGenRec[i]->ProjectionX();
		hPt[i]=(TH1F*)hMCGenRec[i]->ProjectionY();
		hLz[i]=(TH1F*)hMCGenRec[i]->ProjectionZ();
	}
	TH1F *Rec_GenY = (TH1F*)hY[2]->Clone("Rec_GenY");
	Rec_GenY->SetDirectory(0);
	Rec_GenY->Divide(hY[2],hY[0],1,1,"B");
	//TCanvas *eff_y=new TCanvas(Form("eff_%s",xaxis));
	Rec_GenY->SetTitle(Form("Rec over GenMC for y %s",title));
	Rec_GenY->GetXaxis()->SetTitle("y");
	Rec_GenY->GetYaxis()->SetTitle("Efficiency #epsilon");
	//Rec_Gen->Draw("E1");
	TH1F *Gen_GenMCY = (TH1F*)hY[1]->Clone("Gen_GenMC");
	Gen_GenMCY->SetDirectory(0);
	Gen_GenMCY->Divide(hY[1],hY[0],1,1,"B");
	//TCanvas *effGen=new TCanvas(Form("effGen%s",xaxis));
	Gen_GenMCY->SetTitle(Form("Gen over GenMC for y %s",title));
	Gen_GenMCY->GetXaxis()->SetTitle("y");
	Gen_GenMCY->GetYaxis()->SetTitle("Efficiency #epsilon");
	//Gen_GenMC->Draw("E1");
	


		TH1F *Rec_GenPt = (TH1F*)hPt[2]->Clone("Rec_GenPt");
	Rec_GenPt->SetDirectory(0);
	Rec_GenPt->Divide(hPt[2],hPt[0],1,1,"B");
	//TCanvas *eff_Pt=new TCanvas(Form("eff_%s",xaxis));
	Rec_GenPt->SetTitle(Form("Rec over GenMC for Pt %s",title));
	Rec_GenPt->GetXaxis()->SetTitle("Pt");
	Rec_GenPt->GetYaxis()->SetTitle("Efficiency #epsilon");
	//Rec_Gen->Draw("E1");
	TH1F *Gen_GenMCPt = (TH1F*)hPt[1]->Clone("Gen_GenMC");
	Gen_GenMCPt->SetDirectory(0);
	Gen_GenMCPt->Divide(hPt[1],hPt[0],1,1,"B");
	//TCanvas *effGen=new TCanvas(Form("effGen%s",xaxis));
	Gen_GenMCPt->SetTitle(Form("Gen over GenMC for Pt %s",title));
	Gen_GenMCPt->GetXaxis()->SetTitle("Pt");
	Gen_GenMCPt->GetYaxis()->SetTitle("Efficiency #epsilon");
	//Gen_GenMC->Draw("E1");
	

	TH1F *Rec_GenLz = (TH1F*)hLz[2]->Clone("Rec_GenLz");
	Rec_GenLz->SetDirectory(0);
	Rec_GenLz->Divide(hLz[2],hLz[0],1,1,"B");
	//TCanvas *eff_Lz=new TCanvas(Form("eff_%s",xaxis));
	Rec_GenLz->SetTitle(Form("Rec over GenMC for Lz %s",title));
	Rec_GenLz->GetXaxis()->SetTitle("Lz");
	Rec_GenLz->GetYaxis()->SetTitle("Efficiency #epsilon");
	//Rec_Gen->Draw("E1");
	TH1F *Gen_GenMCLz = (TH1F*)hLz[1]->Clone("Gen_GenMC");
	Gen_GenMCLz->SetDirectory(0);
	Gen_GenMCLz->Divide(hLz[1],hLz[0],1,1,"B");
	//TCanvas *effGen=new TCanvas(Form("effGen%s",xaxis));
	Gen_GenMCLz->SetTitle(Form("Gen over GenMC for Lz %s",title));
	Gen_GenMCLz->GetXaxis()->SetTitle("Lz");
	Gen_GenMCLz->GetYaxis()->SetTitle("Efficiency #epsilon");
	//Gen_GenMC->Draw("E1");
	f->cd();
	Gen_GenMCY->SetName(Form("Gen over GenMC for y %s",title));
	Gen_GenMCY->Write();
	Rec_GenY->SetName(Form("Rec over GenMC for y %s",title));
	Rec_GenY->Write();
	Gen_GenMCPt->SetName(Form("Gen over GenMC for Pt %s",title));
	Gen_GenMCPt->Write();
	Rec_GenPt->SetName(Form("Rec over GenMC for Pt %s",title));
	Rec_GenPt->Write();
	Gen_GenMCLz->SetName(Form("Gen over GenMC for Lz %s",title));
	Gen_GenMCLz->Write();
	Rec_GenLz->SetName(Form("Rec over GenMC for Lz %s",title));
	Rec_GenLz->Write();

}
