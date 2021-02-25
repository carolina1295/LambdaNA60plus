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

void EfficiencyLambda(){
    TFile *fLambda = new TFile("LambdastandconBKG-Signal-histos.root");
	TFile *fLambdaBar = new TFile("LambdaBarstandconBKG-Signal-histos.root");
    TFile *fLambdaMoreLayer = new TFile("LambdaMoreLayerconBKG-Signal-histos.root");
	TFile *fLambdaBarMoreLayer = new TFile("LambdaBarMoreLayerconBKG-Signal-histos.root");
	//TFile *fLambdaMoreLayerBKG = new TFile("LambdaMoreLayerconBKG-Signal-histos.root");
	//TFile *fLambdaBarMoreLayerBKG = new TFile("LambdaBarMoreLayerconBKG-Signal-histos.root");
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
		effBar[i]->SetDirectory(0);
		eff[i]->SetMarkerStyle(20);
		effBar[i]->SetMarkerStyle(20);
		eff[i]->SetMarkerSize(0.8);
		effBar[i]->SetMarkerSize(0.8);
		eff[i]->SetMarkerColor(1);
		effBar[i]->SetMarkerColor(1);

		eff10lay[i]->SetDirectory(0);
		effBar10lay[i]->SetDirectory(0);
		eff10lay[i]->SetMarkerStyle(20);
		effBar10lay[i]->SetMarkerStyle(20);
		eff10lay[i]->SetMarkerSize(0.8);
		effBar10lay[i]->SetMarkerSize(0.8);
		eff10lay[i]->SetMarkerColor(2);
		effBar10lay[i]->SetMarkerColor(2);
		effLambda->cd(i+1);
		gPad->SetGrid();
		eff[i]->Draw();
		eff10lay[i]->Draw("SAME");
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

	effLambdaBar0fake->cd(1);
	auto legendLB0fake = new TLegend(0.1,0.7,0.48,0.9);
   	//legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
   	legendLB0fake->AddEntry(effBar10lay[0],"#bar{#Lambda} efficiency for 10 layers","lep");
  	legendLB0fake->AddEntry(effBar10lay0fake[0],"#bar{#Lambda} efficiency for 10 layers 0fake","lep");
   	legendLB0fake->Draw();

	

	TH1D* hPion, *hPionpos, *hProt, *hAntiProt;
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
