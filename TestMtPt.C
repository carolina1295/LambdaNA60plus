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
void TestMtPt(){
	TFile *PtYspectra= new TFile("PtY_spectra.root","RECREATE");
  Double_t mLambda=1.116;
  Int_t totTry=10000.;
  
  TH1F* hmtm0=new TH1F("hmtm0","",100,0.,3.);
  TH1F* hmt=new TH1F("hmt","",100,0.,3.);
  TH1F* hpt=new TH1F("hpt","",100,0.,3.);
	TH1F* hptnotl=new TH1F("hptnotl","",100,0.,3.);
	TH1F* hy=new TH1F("hy","",100,-3.,10.);
	TH1F* hynotl=new TH1F("hynotl","",100,-3.,10.);
	hmtm0->SetDirectory(0);	
	hmt->SetDirectory(0);
	hpt->SetDirectory(0);
	hptnotl->SetDirectory(0);
	hy->SetDirectory(0);
	hynotl->SetDirectory(0);

  // functional form for dN/d(mt-m0) = mt*exponential --- x is (mt-m0)
  TF1* fmtm0=new TF1("fmt","(x+[1])*[0]*exp(-(x/0.301))",0.,3.);
  // functional form for dN/dmt = mt*exponential --- x is mt
  TF1* fmt=new TF1("fmt","x*[0]*exp(-((x-[1])/0.301))",mLambda,3.); //(Check segno + o -)
  // functional form for dN/dpt = pt* exponential --- x is pt
  TF1* fpt=new TF1("fpt","x*[0]*exp(-(TMath::Sqrt(x*x+[1]*[1])-[1])/0.301)",0.,3.);
	TF1* fptnotl=new TF1("fptnotl","x*[0]*exp(-(TMath::Sqrt(x*x+[1]*[1])-[1])/0.303)",0.,3.);
  fmtm0->SetParameter(0,totTry);
  fmt->SetParameter(0,totTry);
  fpt->SetParameter(0,totTry);
	fptnotl->SetParameter(0,totTry);
  fmtm0->FixParameter(1,mLambda);
  fmt->FixParameter(1,mLambda);
  fpt->FixParameter(1,mLambda);
	fptnotl->FixParameter(1,mLambda);
  
  //fmtm0->Draw();

new TCanvas();
//TF1 *f5_lam_rap = new TF1("f5", "[0]*((x>-3 && x<-2)*(8.5*x+25.5) + (x>-2 && x<2)*8.5 + (x>2 && x<3)*(-8.5*x+25.5))", -3, 3);
TF1 *f5_lam_rap = new TF1("f5", "[0]*(((x-2.9)>-3 && (x-2.9)<-2)*(8.5*(x-2.9)+25.5) + ((x-2.9)>-2 && (x-2.9)<2)*8.5 + ((x-2.9)>2 && (x-2.9)<3)*(-8.5*(x-2.9)+25.5))", -3, 10);
f5_lam_rap->SetTitle("Rapidity spectra #Lambda");
f5_lam_rap->GetHistogram()->GetYaxis()->SetTitle("#frac{dN}{dy}");
f5_lam_rap->GetHistogram()->GetXaxis()->SetTitle("y");
f5_lam_rap->Draw();
f5_lam_rap->SetParameter(0,totTry);

//new TCanvas();
TF1 *f6_nlam_rap = new TF1("f6","[0]*(1.2*TMath::Gaus((x-2.9),0,0.95))",-3,10);
f6_nlam_rap->SetTitle("Rapidity spectra #bar{ #Lambda }");
f6_nlam_rap->GetHistogram()->GetYaxis()->SetTitle("#frac{dN}{dy}");
f6_nlam_rap->GetHistogram()->GetXaxis()->SetTitle("y");
//f6_nlam_rap->Draw();
f6_nlam_rap->SetParameter(0,totTry);

  for(Int_t j=0; j<totTry; j++){
    Double_t mtm0=fmtm0->GetRandom();
    Double_t mt=mtm0+mLambda;
    Double_t pt=TMath::Sqrt(mt*mt-mLambda*mLambda);
		Double_t ptnotl=fptnotl->GetRandom();
		Double_t y=f5_lam_rap->GetRandom();
		Double_t ynotl=f6_nlam_rap->GetRandom();
    hmtm0->Fill(mtm0);
    hmt->Fill(mt);
    hpt->Fill(pt);
		hptnotl->Fill(ptnotl);
		hy->Fill(y);
		hynotl->Fill(ynotl);
  }

	hpt->Write();
	hy->Write();
	hptnotl->Write();
	hynotl->Write();
	hmtm0->Fit(fmtm0);
  hmt->Fit(fmt,"","",mLambda,3);
  hpt->Fit(fpt);
	hptnotl->Fit(fptnotl);
	hy->Fit(f5_lam_rap,"","",-3,10);
	hynotl->Fit(f6_nlam_rap,"","",-3,10);
  
	

  TCanvas* ct=new TCanvas("ct","",1500,600);
  ct->Divide(4,1);
  ct->cd(1);
  hmtm0->Draw();
  fmtm0->Draw("same");
  ct->cd(2);
  hmt->Draw();
  fmt->Draw("same");
  ct->cd(3);
  hpt->Draw();
  fpt->Draw("same");
	ct->cd(4);
	hptnotl->Draw();
  fptnotl->Draw("same");

	new TCanvas();
	
  fmtm0->Draw();
  fmt->Draw("same");

	TCanvas* c=new TCanvas("c","",1500,600);
  c->Divide(2,1);
	c->cd(1);
	hy->Draw();
	f5_lam_rap->Draw("same");
	
	c->cd(2);
	hynotl->Draw();
	f6_nlam_rap->Draw("same");

	PtYspectra->Close(); 


  
}
