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
    TFile *hLambda=new TFile("Lambdastand-Signal-histos.root");
	TFile *hLambdaBar=new TFile("LambdaBarstand-Signal-histos.root");
    TFile *hLambdatrasl=new TFile("Lambdatrasl-Signal-histos.root");
	TFile *hLambdaBartrasl=new TFile("LambdaBartrasl-Signal-histos.root");
    TFile *hLambdaMoreLayer=new TFile("LambdaMoreLayer-Signal-histos.root");
	TFile *hLambdaBarMoreLayer=new TFile("LambdaBarMoreLayer-Signal-histos.root");

	TFile *foutLambda=new TFile("Efficiency.root","RECREATE");
    //Lambda
    //Distribuzione 3D per efficienza
	TH3F* hYPtLzMC[3], *hYPtLzGen[3], *hYPtLzRec[3];
	hYPtLzMC[0]=(TH3F*)hLambda->Get("hYPtLzMC");
	hYPtLzGen[0]=(TH3F*)hLambda->Get("hYPtLzGen");
	hYPtLzRec[0]=(TH3F*)hLambda->Get("hYPtLzRec");

	hYPtLzMC[1]=(TH3F*)hLambdatrasl->Get("hYPtLzMC");
	hYPtLzGen[1]=(TH3F*)hLambdatrasl->Get("hYPtLzGen");
	hYPtLzRec[1]=(TH3F*)hLambdatrasl->Get("hYPtLzRec");

	hYPtLzMC[2]=(TH3F*)hLambdaMoreLayer->Get("hYPtLzMC");
	hYPtLzGen[2]=(TH3F*)hLambdaMoreLayer->Get("hYPtLzGen");
	hYPtLzRec[2]=(TH3F*)hLambdaMoreLayer->Get("hYPtLzRec");

	TH1F* hYMC[3],*hYGen[3], *hYRec[3], *hPtMC[3], *hPtGen[3], *hPtRec[3], *hLzMC[3], *hLzGen[3], *hLzRec[3];
	//efficiency for y
	for(int i=0;i<3;i++){
		hYMC[i]=(TH1F*)hYPtLzMC[i]->ProjectionX();
		hYGen[i]=(TH1F*)hYPtLzGen[i]->ProjectionX();
		hYRec[i]=(TH1F*)hYPtLzRec[i]->ProjectionX();
		
		//efficiency for Pt 	
		hPtMC[i]=(TH1F*)hYPtLzMC[i]->ProjectionY();
		hPtGen[i]=(TH1F*)hYPtLzGen[i]->ProjectionY();
		hPtRec[i]=(TH1F*)hYPtLzRec[i]->ProjectionY();
		
		//efficiency for Lz
		hLzMC[i]=(TH1F*)hYPtLzMC[i]->ProjectionZ();
		hLzGen[i]=(TH1F*)hYPtLzGen[i]->ProjectionZ();
		hLzRec[i]=(TH1F*)hYPtLzRec[i]->ProjectionZ();
		
		if(i==0){
			efficiency(hYRec[i],hYMC[i],hYGen[i],"y #Lambda",foutLambda);
			efficiency(hPtRec[i],hPtMC[i],hPtGen[i],"P_{t} #Lambda",foutLambda);
			efficiency(hLzRec[i],hLzMC[i],hLzGen[i],"Lz #Lambda",foutLambda);
		}

		if(i==1){
			efficiency(hYRec[i],hYMC[i],hYGen[i],"y #Lambda 5 layers traslated",foutLambda);
			efficiency(hPtRec[i],hPtMC[i],hPtGen[i],"P_{t} #Lambda 5 layers traslated",foutLambda);
			efficiency(hLzRec[i],hLzMC[i],hLzGen[i],"Lz #Lambda 5 layers traslated",foutLambda);
		}

		if(i==2){
			efficiency(hYRec[i],hYMC[i],hYGen[i],"y #Lambda 10 layers",foutLambda);
			efficiency(hPtRec[i],hPtMC[i],hPtGen[i],"P_{t} #Lambda 10 layers",foutLambda);
			efficiency(hLzRec[i],hLzMC[i],hLzGen[i],"Lz #Lambda 10 layers",foutLambda);
		}
	}

	    //Distribuzione 3D per efficienza
	TH3F* hYPtLzMCBar[3], *hYPtLzGenBar[3], *hYPtLzRecBar[3];
	hYPtLzMCBar[0]=(TH3F*)hLambdaBar->Get("hYPtLzMC");
	hYPtLzGenBar[0]=(TH3F*)hLambdaBar->Get("hYPtLzGen");
	hYPtLzRecBar[0]=(TH3F*)hLambdaBar->Get("hYPtLzRec");

	hYPtLzMCBar[1]=(TH3F*)hLambdaBartrasl->Get("hYPtLzMC");
	hYPtLzGenBar[1]=(TH3F*)hLambdaBartrasl->Get("hYPtLzGen");
	hYPtLzRecBar[1]=(TH3F*)hLambdaBartrasl->Get("hYPtLzRec");

	hYPtLzMCBar[2]=(TH3F*)hLambdaBarMoreLayer->Get("hYPtLzMC");
	hYPtLzGenBar[2]=(TH3F*)hLambdaBarMoreLayer->Get("hYPtLzGen");
	hYPtLzRecBar[2]=(TH3F*)hLambdaBarMoreLayer->Get("hYPtLzRec");

	TH1F* hYMCBar[3],*hYGenBar[3], *hYRecBar[3], *hPtMCBar[3], *hPtGenBar[3], *hPtRecBar[3], *hLzMCBar[3], *hLzGenBar[3], *hLzRecBar[3];
	//efficiency for y
	for(int i=0;i<3;i++){
		hYMCBar[i]=(TH1F*)hYPtLzMCBar[i]->ProjectionX();
		hYGenBar[i]=(TH1F*)hYPtLzGenBar[i]->ProjectionX();
		hYRecBar[i]=(TH1F*)hYPtLzRecBar[i]->ProjectionX();
		
		//efficiency for Pt 	
		hPtMCBar[i]=(TH1F*)hYPtLzMCBar[i]->ProjectionY();
		hPtGenBar[i]=(TH1F*)hYPtLzGenBar[i]->ProjectionY();
		hPtRecBar[i]=(TH1F*)hYPtLzRecBar[i]->ProjectionY();
		
		//efficiency for Lz
		hLzMCBar[i]=(TH1F*)hYPtLzMCBar[i]->ProjectionZ();
		hLzGenBar[i]=(TH1F*)hYPtLzGenBar[i]->ProjectionZ();
		hLzRecBar[i]=(TH1F*)hYPtLzRecBar[i]->ProjectionZ();
		
		if(i==0){
			efficiency(hYRecBar[i],hYMCBar[i],hYGenBar[i],"y #Lambda Bar",foutLambda);
			efficiency(hPtRecBar[i],hPtMCBar[i],hPtGenBar[i],"P_{t} #Lambda Bar",foutLambda);
			efficiency(hLzRecBar[i],hLzMCBar[i],hLzGenBar[i],"Lz #Lambda Bar",foutLambda);
		}

		if(i==1){
			efficiency(hYRecBar[i],hYMCBar[i],hYGenBar[i],"y #Lambda Bar 5 layers traslated",foutLambda);
			efficiency(hPtRecBar[i],hPtMCBar[i],hPtGenBar[i],"P_{t} #Lambda Bar 5 layers traslated",foutLambda);
			efficiency(hLzRecBar[i],hLzMCBar[i],hLzGenBar[i],"Lz #Lambda Bar 5 layers traslated",foutLambda);
		}

		if(i==2){
			efficiency(hYRecBar[i],hYMCBar[i],hYGenBar[i],"y #Lambda Bar 10 layers",foutLambda);
			efficiency(hPtRecBar[i],hPtMCBar[i],hPtGenBar[i],"P_{t} #Lambda Bar 10 layers",foutLambda);
			efficiency(hLzRecBar[i],hLzMCBar[i],hLzGenBar[i],"Lz #Lambda Bar 10 layers",foutLambda);
		}
	}
	TCanvas *effLambda=new TCanvas("effLambda");
	TCanvas *effLambdaBar=new TCanvas("effLambdaBar");
	

	TCanvas *effLnotL=new TCanvas("effLnotL");
	TCanvas *effLnotLtrasl=new TCanvas("effLnotLtrasl");
	TCanvas *effLnotL10lay=new TCanvas("effLnotL10lay");
	TH1F* eff[3], *effBar[3], *efftrasl[3], *effBartrasl[3], *eff10lay[3], *effBar10lay[3];
	eff[0]=(TH1F*)foutLambda->Get("Gen over GenMC for y #Lambda");
	eff[1]=(TH1F*)foutLambda->Get("Gen over GenMC for P_{t} #Lambda");	
	eff[2]=(TH1F*)foutLambda->Get("Gen over GenMC for Lz #Lambda");
	effBar[0]=(TH1F*)foutLambda->Get("Gen over GenMC for y #Lambda Bar");
	effBar[1]=(TH1F*)foutLambda->Get("Gen over GenMC for P_{t} #Lambda Bar");	
	effBar[2]=(TH1F*)foutLambda->Get("Gen over GenMC for Lz #Lambda Bar");
	efftrasl[0]=(TH1F*)foutLambda->Get("Gen over GenMC for y #Lambda 5 layers traslated");
	efftrasl[1]=(TH1F*)foutLambda->Get("Gen over GenMC for P_{t} #Lambda 5 layers traslated");	
	efftrasl[2]=(TH1F*)foutLambda->Get("Gen over GenMC for Lz #Lambda 5 layers traslated");
	effBartrasl[0]=(TH1F*)foutLambda->Get("Gen over GenMC for y #Lambda Bar 5 layers traslated");
	effBartrasl[1]=(TH1F*)foutLambda->Get("Gen over GenMC for P_{t} #Lambda Bar 5 layers traslated");	
	effBartrasl[2]=(TH1F*)foutLambda->Get("Gen over GenMC for Lz #Lambda Bar 5 layers traslated");
	eff10lay[0]=(TH1F*)foutLambda->Get("Gen over GenMC for y #Lambda 10 layers");
	eff10lay[1]=(TH1F*)foutLambda->Get("Gen over GenMC for P_{t} #Lambda 10 layers");	
	eff10lay[2]=(TH1F*)foutLambda->Get("Gen over GenMC for Lz #Lambda 10 layers");
	effBar10lay[0]=(TH1F*)foutLambda->Get("Gen over GenMC for y #Lambda Bar 10 layers");
	effBar10lay[1]=(TH1F*)foutLambda->Get("Gen over GenMC for P_{t} #Lambda Bar 10 layers");	
	effBar10lay[2]=(TH1F*)foutLambda->Get("Gen over GenMC for Lz #Lambda Bar 10 layers");
	effLnotL->Divide(3,1);
	effLnotLtrasl->Divide(3,1);
	effLnotL10lay->Divide(3,1);
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
		effLnotL->cd(i+1);
		gPad->SetGrid();
		eff[i]->Draw();
		effBar[i]->Draw("SAME");
		efftrasl[i]->SetDirectory(0);
		effBartrasl[i]->SetDirectory(0);
		efftrasl[i]->SetMarkerStyle(20);
		effBartrasl[i]->SetMarkerStyle(20);
		efftrasl[i]->SetMarkerSize(0.8);
		effBartrasl[i]->SetMarkerSize(0.8);
		efftrasl[i]->SetMarkerColor(2);
		effBartrasl[i]->SetMarkerColor(2);
		effLnotLtrasl->cd(i+1);
		gPad->SetGrid();
		efftrasl[i]->Draw();
		effBartrasl[i]->Draw("SAME");
		eff10lay[i]->SetDirectory(0);
		effBar10lay[i]->SetDirectory(0);
		eff10lay[i]->SetMarkerStyle(20);
		effBar10lay[i]->SetMarkerStyle(20);
		eff10lay[i]->SetMarkerSize(0.8);
		effBar10lay[i]->SetMarkerSize(0.8);
		eff10lay[i]->SetMarkerColor(4);
		effBar10lay[i]->SetMarkerColor(4);
		effLnotL10lay->cd(i+1);
		gPad->SetGrid();
		eff10lay[i]->Draw();
		effBar10lay[i]->Draw("SAME");
		effLambda->cd(i+1);
		gPad->SetGrid();
		eff[i]->Draw();
		efftrasl[i]->Draw("SAME");
		eff10lay[i]->Draw("SAME");
		effLambdaBar->cd(i+1);
		gPad->SetGrid();
		effBar[i]->Draw();
		effBartrasl[i]->Draw("SAME");
		effBar10lay[i]->Draw("SAME");
	}

    //Distribuzione 3D per efficienza Lambda Bar
/*	TH3F* hYPtLzMCBar=(TH3F*)hLambdaBar->Get("hYPtLzMC");
	TH3F* hYPtLzGenBar=(TH3F*)hLambdaBar->Get("hYPtLzGen");
	TH3F* hYPtLzRecBar=(TH3F*)hLambdaBar->Get("hYPtLzRec");
	//efficiency for y
	TH1F* hYMCBar=(TH1F*)hYPtLzMCBar->ProjectionX();
	TH1F* hYGenBar=(TH1F*)hYPtLzGenBar->ProjectionX();
	TH1F* hYRecBar=(TH1F*)hYPtLzRecBar->ProjectionX();
	efficiency(hYRecBar,hYMCBar,hYGenBar,"y #Lambda Bar",foutLambda);
	//efficiency for Pt 	
	TH1F* hPtMCBar=(TH1F*)hYPtLzMCBar->ProjectionY();
	TH1F* hPtGenBar=(TH1F*)hYPtLzGenBar->ProjectionY();
	TH1F* hPtRecBar=(TH1F*)hYPtLzRecBar->ProjectionY();
	efficiency(hPtRecBar,hPtMCBar,hPtGenBar,"P_{t} #Lambda Bar",foutLambda);
	//efficiency for Lz
	TH1F* hLzMCBar=(TH1F*)hYPtLzMCBar->ProjectionZ();
	TH1F* hLzGenBar=(TH1F*)hYPtLzGenBar->ProjectionZ();
	TH1F* hLzRecBar=(TH1F*)hYPtLzRecBar->ProjectionZ();
	efficiency(hLzRecBar,hLzMCBar,hLzGenBar,"Lz #Lambda Bar",foutLambda);*/
	
	/*TCanvas *effLnotL=new TCanvas("effLnotL");
	TH1F* eff[9], *effBar[9];
	eff[0]=(TH1F*)foutLambda.Get("Gen over GenMC for y #Lambda");
	eff[1]=(TH1F*)foutLambda.Get("Gen over GenMC for P_{t} #Lambda");	
	eff[2]=(TH1F*)foutLambda.Get("Gen over GenMC for Lz #Lambda");
	effBar[0]=(TH1F*)foutLambda.Get("Gen over GenMC for y #Lambda Bar");
	effBar[1]=(TH1F*)foutLambda.Get("Gen over GenMC for P_{t} #Lambda Bar");	
	effBar[2]=(TH1F*)foutLambda.Get("Gen over GenMC for Lz #Lambda Bar");

	effLnotL->Divide(3,1);
	for(int i=0;i<3;i++){
		eff[i]->SetMarkerStyle(20);
		effBar[i]->SetMarkerStyle(20);
		eff[i]->SetMarkerSize(0.9);
		effBar[i]->SetMarkerSize(0.9);
		eff[i]->SetMarkerColor(1);
		effBar[i]->SetMarkerColor(2);
		effLnotL->cd(i+1);
		eff[i]->Draw();
		effBar[i]->Draw("SAME");
	}
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
