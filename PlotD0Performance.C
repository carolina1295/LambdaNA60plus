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
#include <TNtuple.h>

void ProjectionforBinMassInv(TH2F* hist2D, char* graphTitle , char* Xaxisname, TCanvas *c0, TFile *f);
void ProjectionforBin(TH2F* hist2D, char* graphTitle , char* Xaxisname, TCanvas *c0, TFile *f);

void PlotD0Performance()
{
  TFile *fin = new TFile("D0-Signal-histos.root");
	TFile *foutD0=new TFile("histogramD0.root","RECREATE");
  // Histograms for inv mass resolution studies:
  //mass of candidates with correct K-pi assignemnt vs pt and y
  TH2F* hMassVsPt=(TH2F*)fin->Get("hMassVsPt");
  TH2F* hMassVsY=(TH2F*)fin->Get("hMassVsY");
  // Histograms for inv mass resolution studies:
  //mass of candidates with swapped K-pi assignemnt vs pt and y
  TH2F* hMassReflVsPt=(TH2F*)fin->Get("hMassReflVsPt");
  TH2F* hMassReflVsY=(TH2F*)fin->Get("hMassReflVsY");

  // Histograms for efficiency calculation (vs pt and y)
  TH1F* PtGen=(TH1F*)fin->Get("PTGen");
  TH1F* hYGen=(TH1F*)fin->Get("hYGen");
  TH1F* hPtReco=(TH1F*)fin->Get("PTAll");
  TH1F* hYReco=(TH1F*)fin->Get("yall");
	// Histograms for efficiency*acceptance calculation 
 	TH1F* hPtGen=(TH1F*)fin->Get("hPtGen");
  TH1F* hPtRecoAll=(TH1F*)fin->Get("hPtRecoAll");
  TH1F* hPtGenRecoAll=(TH1F*)fin->Get("hPtGenRecoAll");
  // Histogram of decay length (x axis = decay length, y axis = pt)
  TH2F* hDecLen=(TH2F*)fin->Get("hDist");
  // Histogram of cosine of pointing angle (x axis = cos(point), y axis = pt)
  TH2F* hCosp=(TH2F*)fin->Get("hCosp");

  // Histograms of decay vertex position (x axis = difference between measured and reconstructed position, y axis = pt)
  TH2F *hResVx = (TH2F*)fin->Get("hResVx");
  TH2F *hResVy = (TH2F*)fin->Get("hResVy");
  TH2F *hResVz = (TH2F*)fin->Get("hResVz");
  // Histograms of decay vertex position (x axis = difference between measured and reconstructed position, y axis = rapidity)
 /* TH2F *hResVxVsY = (TH2F*)fin->Get("hResVxVsY");
  TH2F *hResVyVsY = (TH2F*)fin->Get("hResVyVsY");
  TH2F *hResVzVsY = (TH2F*)fin->Get("hResVzVsY");*/

TFile filin("ntuplaD0.root");
  float nfaketrkKaon, nfaketrkpion,y,xrec,yrec,zrec,secvertgenKaon[3],massinvrec,pResKaon[3],pResPion[3],dca,PrecKaon[3],PrecPion[3];
  TNtuple *variables = (TNtuple*)filin.Get("ntD0;1");
		variables->SetBranchAddress("nfaketrkKaon",&nfaketrkKaon);
    variables->SetBranchAddress("nfaketrkpion",&nfaketrkpion);
    variables->SetBranchAddress("ygen",&y);
    variables->SetBranchAddress("xP",&xrec);
		variables->SetBranchAddress("yP",&yrec);
		variables->SetBranchAddress("zP",&zrec);
		variables->SetBranchAddress("Vxgen",&secvertgenKaon[0]);
		variables->SetBranchAddress("Vygen",&secvertgenKaon[1]);
		variables->SetBranchAddress("Vzgen",&secvertgenKaon[2]);
		variables->SetBranchAddress("massinvrec",&massinvrec);
		variables->SetBranchAddress("diffxprot",&pResKaon[0]);
		variables->SetBranchAddress("diffyprot",&pResKaon[1]);
		variables->SetBranchAddress("diffzprot",&pResKaon[2]);
		variables->SetBranchAddress("diffxpion",&pResPion[0]);
		variables->SetBranchAddress("diffypion",&pResPion[1]);
		variables->SetBranchAddress("diffzpion",&pResPion[2]);
		variables->SetBranchAddress("dca",&dca);
		variables->SetBranchAddress("pxrecprot",&PrecKaon[0]);
		variables->SetBranchAddress("pxrecprot",&PrecKaon[1]);
		variables->SetBranchAddress("pxrecprot",&PrecKaon[2]);
		variables->SetBranchAddress("pxrecpion",&PrecPion[0]);
		variables->SetBranchAddress("pxrecpion",&PrecPion[1]);
		variables->SetBranchAddress("pxrecpion",&PrecPion[2]);

	int events=variables->GetEntries();
	double theta=0;
	TH2F *hMassinrecY = new TH2F("hMassinrecVsY", "Mass vs Y", 3000, 0,3,60, 0, 6);	
	TH2F *hResVxVsY = new TH2F("hResVxVsY", "Res Vx vs Y", 200, -200., 200., 36, 0, 6);
  TH2F *hResVyVsY = new TH2F("hResVyVsY", "Res Vy vs Y", 200, -200., 200., 36, 0, 6);
  TH2F *hResVzVsY = new TH2F("hResVzVsY", "Res Vz vs Y", 200, -2000., 2000., 36, 0, 6);
	TH2F *hResPxVsYKaon = new TH2F("hResPxVsYKaon", "Res Px vs Y", 200, -5., 5., 36, 0, 6);
	TH2F *hResPxVsYpion = new TH2F("hResPxVsYpion", "Res Px vs Y", 200, -5., 5., 36, 0, 6);
	TH2F *hResPyVsYKaon = new TH2F("hResPyVsYKaon", "Res Py vs Y", 200, -5., 5., 36, 0, 6);
	TH2F *hResPyVsYpion = new TH2F("hResPyVsYpion", "Res Py vs Y", 200, -5., 5., 36, 0, 6);
	TH2F *hResPzVsYKaon = new TH2F("hResPzVsYKaon", "Res Pz vs Y", 200, -5., 5., 36, 0, 6);
	TH2F *hResPzVsYpion = new TH2F("hResPzVsYpion", "Res Pz vs Y", 200, -5., 5., 36, 0, 6);
	TH1F *hPxKaon = new TH1F("hPxKaon", " Px Kaonon ", 200, -5., 5.);
	TH1F *hPyKaon = new TH1F("hPyKaon", " Py Kaonon", 200, -5., 5.);
	TH1F *hPzKaon = new TH1F("hPzKaon", " Pz Kaonon ", 200, -5., 5.);
	TH1F *hPxpion = new TH1F("hPxpion", " Px pion", 200, -5., 5.);
	TH1F *hPypion = new TH1F("hPypion", " Py pion", 200, -5., 5.);
	TH1F *hPzpion = new TH1F("hPzpion", " Pz pion", 200, -5., 5.);
	TH1F *htheta = new TH1F("htheta", " Angle between the particles", 200, 0., 7.);
	TH2F *hResPxVsYKaonRel = new TH2F("hResPxVsYKaonRel", "Res Px vs Y Rel Kaon", 200, -5., 5., 36, 0, 6);
	TH2F *hResPxVsYpionRel = new TH2F("hResPxVsYpionRel", "Res Px vs Y Rel Pion", 200, -5., 5., 36, 0, 6);
	TH2F *hResPyVsYKaonRel = new TH2F("hResPyVsYKaonRel", "Res Py vs Y Rel Kaon", 200, -5., 5., 36, 0, 6);
	TH2F *hResPyVsYpionRel = new TH2F("hResPyVsYpionRel", "Res Py vs Y Rel Pion", 200, -5., 5., 36, 0, 6);
	TH2F *hResPzVsYKaonRel = new TH2F("hResPzVsYKaonRel", "Res Pz vs Y Rel Kaon", 200, -5., 5., 36, 0, 6);
	TH2F *hResPzVsYpionRel = new TH2F("hResPzVsYpionRel", "Res Pz vs Y Rel Pion", 200, -5., 5., 36, 0, 6);

	for(int i=0 ; i<events; i++) {
		variables->GetEvent(i);
		double residVx=10000.*(xrec - secvertgenKaon[0]);
    double residVy=10000.*(yrec - secvertgenKaon[1]);
    double residVz=10000.*(zrec - secvertgenKaon[2]);
		double p1abs=sqrt(PrecKaon[0]*PrecKaon[0]+PrecKaon[1]*PrecKaon[1]+PrecKaon[2]*PrecKaon[2]);
		double p2abs=sqrt(PrecPion[0]*PrecPion[0]+PrecPion[1]*PrecPion[1]+PrecPion[2]*PrecPion[2]);
		theta=TMath::ACos((PrecKaon[0]*PrecPion[0]+PrecKaon[1]*PrecPion[1]+PrecKaon[2]*PrecPion[2])/(p1abs*p2abs));
		hMassinrecY->Fill(massinvrec, y);
		hResVxVsY->Fill(residVx, y);
    hResVyVsY->Fill(residVy, y);
    hResVzVsY->Fill(residVz, y);
		hResPxVsYKaon->Fill(pResKaon[0],y);
		hResPxVsYpion->Fill(pResPion[0],y);
		hResPyVsYKaon->Fill(pResKaon[1],y);
		hResPyVsYpion->Fill(pResPion[1],y);
		hResPzVsYKaon->Fill(pResKaon[2],y);
		hResPzVsYpion->Fill(pResPion[2],y);

		hResPxVsYKaonRel->Fill(pResKaon[0]/PrecKaon[0],y);
		hResPxVsYpionRel->Fill(pResPion[0]/PrecPion[0],y);
		hResPyVsYKaonRel->Fill(pResKaon[1]/PrecKaon[1],y);
		hResPyVsYpionRel->Fill(pResPion[1]/PrecPion[1],y);
		hResPzVsYKaonRel->Fill(pResKaon[2]/PrecKaon[2],y);
		hResPzVsYpionRel->Fill(pResPion[2]/PrecPion[2],y);

		hPxKaon->Fill(PrecKaon[0]);
		hPyKaon->Fill(PrecKaon[1]);
		hPzKaon->Fill(PrecKaon[2]);
		hPxpion->Fill(PrecPion[0]);
		hPypion->Fill(PrecPion[1]);
		hPzpion->Fill(PrecPion[2]);
		htheta->Fill(theta);
		}
	foutD0->cd();
	hMassinrecY->Write();
	hResVxVsY->Write();
  	hResVyVsY->Write();
  	hResVzVsY->Write();

	hResPxVsYKaon->Write();
	hResPxVsYpion->Write();
	hResPyVsYKaon->Write();
	hResPyVsYpion->Write();
	hResPzVsYKaon->Write();
	hResPzVsYpion->Write();

	hResPxVsYKaonRel->Write();
	hResPxVsYpionRel->Write();
	hResPyVsYKaonRel->Write();
	hResPyVsYpionRel->Write();
	hResPzVsYKaonRel->Write();
	hResPzVsYpionRel->Write();

	hPxKaon->Write();
	hPyKaon->Write();
	hPzKaon->Write();
	hPxpion->Write();
	hPypion->Write();
	hPzpion->Write();

	htheta->Write();

	//Efficiency*acceptance 
	TH1F *GenRec_ptGen = (TH1F*)hPtGenRecoAll->Clone("GenRec_ptGen");
	GenRec_ptGen->Divide(hPtGenRecoAll,hPtGen,1,1,"B");
	GenRec_ptGen->SetDirectory(0);
	TCanvas *eff=new TCanvas("eff");
	GenRec_ptGen->SetTitleSize(20);
	GenRec_ptGen->SetTitle("Efficiency [#frac{GenRecAll}{PtGen}]");
	GenRec_ptGen->SetMarkerStyle(20);
	GenRec_ptGen->SetMarkerSize(0.8);
	GenRec_ptGen->GetYaxis()->SetTitle("Efficiency #epsilon");
	GenRec_ptGen->Draw("E1");
	TH1F *Rec_ptGen = (TH1F*)hPtRecoAll->Clone("Rec_ptGen");
	Rec_ptGen->Divide(hPtRecoAll,hPtGen,1,1,"B");
	Rec_ptGen->SetDirectory(0);
	TCanvas *eff1=new TCanvas("eff1");
	Rec_ptGen->SetTitle("Efficiency [#frac{RecAll}{PtGen}]");
	Rec_ptGen->SetMarkerStyle(20);
	Rec_ptGen->SetMarkerSize(0.8);
	Rec_ptGen->GetYaxis()->SetTitle("Efficiency #epsilon");
	Rec_ptGen->Draw("E1");

		TCanvas* c2p=new TCanvas("c2p");
	ProjectionforBin(hResPxVsYKaon ,"#D_{0} ResPx vs Y", "Y", c2p, foutD0);
	TCanvas* c3p=new TCanvas("c3p");
	ProjectionforBin(hResPyVsYKaon, "#D_{0} ResPy vs Y", "Y", c3p, foutD0);
	TCanvas* c4p=new TCanvas("c4p");
	ProjectionforBin(hResPzVsYKaon, "#D_{0} ResPz vs Y", "Y", c4p, foutD0);
	TCanvas* c5p=new TCanvas("c5p");
	ProjectionforBin(hResPxVsYpion, "D_{0} ResPx vs Y", "Y", c5p, foutD0);
	TCanvas* c6p=new TCanvas("c6p");
	ProjectionforBin(hResPyVsYpion, "D_{0} ResPy vs Y", "Y", c6p, foutD0);
	TCanvas* c7p=new TCanvas("c7p");
	ProjectionforBin(hResPzVsYpion, "D_{0} ResPz vs Y", "Y", c7p, foutD0);

	TCanvas* c2pRel=new TCanvas("c2pRel");
	ProjectionforBin(hResPxVsYKaonRel ,"#D_{0} ResPx vs Y", "Y", c2pRel, foutD0);
	TCanvas* c3pRel=new TCanvas("c3pRel");
	ProjectionforBin(hResPyVsYKaonRel, "#D_{0} ResPy vs Y", "Y", c3pRel, foutD0);
	TCanvas* c4pRel=new TCanvas("c4pRel");
	ProjectionforBin(hResPzVsYKaonRel, "#D_{0} ResPz vs Y", "Y", c4pRel, foutD0);
	TCanvas* c5pRel=new TCanvas("c5pRel");
	ProjectionforBin(hResPxVsYpionRel, "D_{0} ResPx vs Y", "Y", c5pRel, foutD0);
	TCanvas* c6pRel=new TCanvas("c6pRel");
	ProjectionforBin(hResPyVsYpionRel, "D_{0} ResPy vs Y", "Y", c6pRel, foutD0);
	TCanvas* c7pRel=new TCanvas("c7pRel");
	ProjectionforBin(hResPzVsYpionRel, "D_{0} ResPz vs Y", "Y", c7pRel, foutD0);
	
	

	TCanvas* c0=new TCanvas("c0");
	ProjectionforBinMassInv(hMassVsPt, "D_{0} Inv Mass vs P_{t}", "Pt", c0, foutD0);
	TCanvas* c1=new TCanvas("c1");
	ProjectionforBinMassInv(hMassinrecY, "D_{0} Inv Mass vs Y", "Y", c1, foutD0);
	TCanvas* c2=new TCanvas("c2");
	ProjectionforBin(hResVxVsY, "D_{0} Res Vx vs Y", "Y", c2, foutD0);
	TCanvas* c3=new TCanvas("c3");
	ProjectionforBin(hResVyVsY, "D_{0} Res Vy vs Y", "Y", c3, foutD0);
	TCanvas* c4=new TCanvas("c4");
	ProjectionforBin(hResVzVsY, "D_{0} Res Vz vs Y", "Y", c4, foutD0);
	TCanvas* c5=new TCanvas("c5");
	ProjectionforBin(hResVx, "D_{0} Res Vx vs Pt", "Pt", c5, foutD0);
	TCanvas* c6=new TCanvas("c6");
	ProjectionforBin(hResVy, "D_{0} Res Vy vs Pt", "Pt", c6, foutD0);
	TCanvas* c7=new TCanvas("c7");
	ProjectionforBin(hResVz, "D_{0} Res Vz vs Pt", "Pt", c7, foutD0);

	/*TCanvas* c8=new TCanvas("c8");
	c8->Divide(3,1);
	c8->cd(1);
	TGraphErrors* gr11 = (TGraphErrors*)foutD0->Get("Graph;8");
	TGraphErrors* gr11g = (TGraphErrors*)foutD0->Get("Graph;9");
	gr11g->Draw("A*");
	gr11->Draw("P");
	TLegend* leg1 = new TLegend(0.11,0.7,0.3,0.88);
	leg1->SetHeader("Legend");
	leg1->AddEntry(gr11g,"#sigma","p");
	leg1->AddEntry(gr11,"RMS","p");
	leg1->Draw();
	c8->cd(2);
	TGraphErrors* gr12 = (TGraphErrors*)foutD0->Get("Graph;11");
	TGraphErrors* gr12g = (TGraphErrors*)foutD0->Get("Graph;12");
	gr12g->Draw("A*");
	gr12->Draw("P");
	TLegend* leg2 = new TLegend(0.11,0.7,0.3,0.88);
	leg2->SetHeader("Legend");
	leg2->AddEntry(gr12g,"#sigma","p");
	leg2->AddEntry(gr12,"RMS","p");
	leg2->Draw();
	c8->cd(3);
	TGraphErrors* gr13 = (TGraphErrors*)foutD0->Get("Graph;14");
	TGraphErrors* gr13g = (TGraphErrors*)foutD0->Get("Graph;15");
	gr13g->Draw("A*");
	gr13->Draw("P");
	TLegend* leg3 = new TLegend(0.11,0.7,0.3,0.88);
	leg3->SetHeader("Legend");
	leg3->AddEntry(gr13g,"#sigma","p");
	leg3->AddEntry(gr13,"RMS","p");
	leg3->Draw();

	TCanvas* c9=new TCanvas("c9");
	c9->Divide(3,1);
	c9->cd(1);

	TGraphErrors* gr14 = (TGraphErrors*)foutD0->Get("Graph;17");
	TGraphErrors* gr14g = (TGraphErrors*)foutD0->Get("Graph;18");
	gr14g->Draw("A*");
	gr14->Draw("P");
	TLegend* leg4 = new TLegend(0.11,0.7,0.3,0.88);
	leg4->SetHeader("Legend");
	leg4->AddEntry(gr14g,"#sigma","p");
	leg4->AddEntry(gr14,"RMS","p");
	leg4->Draw();
	c9->cd(2);
	TGraphErrors* gr15 = (TGraphErrors*)foutD0->Get("Graph;20");
	TGraphErrors* gr15g = (TGraphErrors*)foutD0->Get("Graph;21");
	gr15g->Draw("A*");
	gr15->Draw("P");
	TLegend* leg5 = new TLegend(0.11,0.7,0.3,0.88);
	leg5->SetHeader("Legend");
	leg5->AddEntry(gr15g,"#sigma","p");
	leg5->AddEntry(gr15,"RMS","p");
	leg5->Draw();
	c9->cd(3);
	TGraphErrors* gr16 = (TGraphErrors*)foutD0->Get("Graph;23");
	TGraphErrors* gr16g = (TGraphErrors*)foutD0->Get("Graph;24");
	gr16g->Draw("A*");
	gr16->Draw("P");
	TLegend* leg6 = new TLegend(0.11,0.7,0.3,0.88);
	leg6->SetHeader("Legend");
	leg6->AddEntry(gr16g,"#sigma","p");
	leg6->AddEntry(gr16,"RMS","p");
	leg6->Draw();

*/
	
foutD0->Close();

}

void ProjectionforBinMassInv(TH2F* hist2D, char* graphTitle , char* Xaxisname, TCanvas *c0, TFile *f){
	double xaxis[100], mean[100], errmean[100], sigma[100], errsigma[100], sigmagaus[100], errsigmagaus[100];
	int count=0;
  c0->Divide(6,6);
  for(Int_t jj=1; jj<hist2D->GetNbinsY(); jj++){
 		printf("n=%d",hist2D->GetNbinsY());
    TH1D *hInvMass = (TH1D *)hist2D->ProjectionX(Form("%sBin%d",hist2D->GetName(),jj), jj, jj);
		//hInvMass->SetDirectory(0);
    //hInvMass->Draw();
    if(hInvMass->GetSumOfWeights()>20){
			count++;
			c0->cd(count);
			c0->SetTitle(Form("%s",hist2D->GetName()));
			//new TCanvas();
			hInvMass->SetDirectory(0);
			//hInvMass->SetMinimum(0);
			//hInvMass->GetXaxis()->SetRangeUser(1.05, 1.20);
    	hInvMass->Draw();
      Double_t rmsh=hInvMass->GetRMS();
      Double_t meanh=hInvMass->GetMean();
			Double_t Errrmsh=hInvMass->GetRMSError();
      Double_t Errmeanh=hInvMass->GetMeanError();

		TF1* fg=((TF1*)(gROOT->GetFunction("gaus")));
		//	fg->SetParameter(1,1.116);
		//	fg->SetParameter(2,0.0015);
		//	fg->SetParameter(0,hInvMass->GetMaximum());
      hInvMass->Fit(fg,"","",meanh-5*rmsh,meanh+5*rmsh);

      Double_t meang=fg->GetParameter(1);
      Double_t emeang=fg->GetParError(1);
      Double_t sigmag=fg->GetParameter(2);
      Double_t esigma=fg->GetParError(2);
      printf(" interval %f-%f   mean=%f+-%f  sigma=%f+-%f\n",
	    	hist2D->GetYaxis()->GetBinLowEdge(jj),
	    	hist2D->GetYaxis()->GetBinUpEdge(jj),
	    	meanh,Errmeanh,rmsh,Errrmsh);
			xaxis[count-1]=(hist2D->GetYaxis()->GetBinUpEdge(jj));
			mean[count-1]=meanh;
			errmean[count-1]=Errmeanh;
			sigma[count-1]=rmsh;
			errsigma[count-1]=Errrmsh;
			sigmagaus[count-1]=sigmag;
			errsigmagaus[count-1]=esigma;
    }
  }
	printf("count=%d",count); 
	TGraphErrors *gr2=new TGraphErrors(count,xaxis,mean,0,errmean);
		gr2->SetTitle(Form("Mean  %s ",graphTitle));
		gr2->GetXaxis()->SetTitle(Form(" %s ",Xaxisname));
		gr2->SetMarkerColor(2);
		gr2->GetYaxis()->SetTitle("Mean");
		gr2->SetMarkerStyle(20);
		gr2->SetMarkerSize(0.7);
		new TCanvas(Form("Graph vs  %s ",graphTitle),"Mean ",200,10,800,500);
		gr2->Draw();
	TGraphErrors *gr1=new TGraphErrors(count,xaxis,sigma,0,errsigma);
		gr1->SetTitle(Form("RMS  %s ",graphTitle));
		gr1->GetXaxis()->SetTitle(Form(" %s ",Xaxisname));
		gr1->GetYaxis()->SetTitle("RMS");
		gr1->SetMarkerColor(2);
		gr1->SetMarkerStyle(20);
		gr1->SetMarkerSize(0.7);
	TGraphErrors *gr1gaus=new TGraphErrors(count,xaxis,sigmagaus,0,errsigmagaus);
		gr1gaus->SetTitle(Form("#sigma  %s ",graphTitle));
		gr1gaus->GetXaxis()->SetTitle(Form(" %s ",Xaxisname));
		gr1gaus->GetYaxis()->SetTitle("");
		gr1gaus->SetMarkerColor(4);
		gr1gaus->SetMarkerStyle(29);
		gr1gaus->SetMarkerSize(0.7);
		new TCanvas(Form("Graph1 vs  %s ",graphTitle),"RMS ",200,10,800,500);
		gr1->Draw();
	f->cd();
	gr2->Write();
	gr1->Write();
	gr1gaus->Write();
}

void ProjectionforBin(TH2F* hist2D, char* graphTitle , char* Xaxisname, TCanvas *c0, TFile *f){
	double xaxis[100], mean[100], errmean[100], sigma[100], errsigma[100], sigmagaus[100], errsigmagaus[100];
	int count=0, i=0, countGraph=0;
  c0->Divide(6,6);
  for(Int_t jj=1; jj<hist2D->GetNbinsY(); jj++){
 		//printf("n=%d ",hist2D->GetNbinsY());
    TH1D *hInvMass = (TH1D *)hist2D->ProjectionX(Form("%sBin%d",hist2D->GetName(),jj), jj, jj);
    if(hInvMass->GetSumOfWeights()>20){
			count++;
			c0->cd(count);
			c0->SetTitle(Form("%s",hist2D->GetName()));
			//new TCanvas();
			hInvMass->SetDirectory(0);
    	hInvMass->Draw();
      Double_t rmsh=hInvMass->GetRMS();
      Double_t meanh=hInvMass->GetMean();
			Double_t Errrmsh=hInvMass->GetRMSError();
      Double_t Errmeanh=hInvMass->GetMeanError();
			hInvMass->Fit("gaus","","",meanh-5*rmsh,meanh+5*rmsh);
    	TF1* fg=(TF1*)hInvMass->GetListOfFunctions()->FindObject("gaus");
      Double_t meang=fg->GetParameter(1);
      Double_t emeang=fg->GetParError(1);
      Double_t sigmag=fg->GetParameter(2);
      Double_t esigma=fg->GetParError(2);
      printf("%s interval %f-%f   mean=%f+-%f  sigma=%f+-%f , SIGMAGAUS %f, ERRORESIGMA %f\n",Xaxisname,
	    	hist2D->GetYaxis()->GetBinLowEdge(jj),
	    	hist2D->GetYaxis()->GetBinUpEdge(jj),
	    	meanh,Errmeanh,rmsh,Errrmsh,sigmag,esigma);
			Double_t ErrMinSigma=0;
			ErrMinSigma=(sigmag*50)/100;
			if(esigma<ErrMinSigma){
			xaxis[countGraph]=(hist2D->GetYaxis()->GetBinUpEdge(jj));
			mean[countGraph]=meanh;
			errmean[countGraph]=Errmeanh;
			sigma[countGraph]=rmsh;
			errsigma[countGraph]=Errrmsh;
			sigmagaus[countGraph]=sigmag;
			errsigmagaus[countGraph]=esigma;
			countGraph++;
			}
    }
  }
	printf("count=%d, countGraph=%d", count, countGraph); 
	TGraphErrors *gr2=new TGraphErrors(countGraph,xaxis,mean,0,errmean);
		gr2->SetTitle(Form("Mean  %s ",graphTitle));
		gr2->GetXaxis()->SetTitle(Form(" %s ",Xaxisname));
		gr2->SetMarkerColor(2);
		gr2->GetYaxis()->SetTitle("Mean");
		gr2->SetMarkerStyle(20);
		//new TCanvas(Form("Graph vs  %s ",graphTitle),"Mean",200,10,800,500);
		//gr2->Draw();
	TGraphErrors *gr1=new TGraphErrors(countGraph,xaxis,sigma,0,errsigma);
		gr1->SetTitle(Form("RMS  %s ",graphTitle));
		gr1->GetXaxis()->SetTitle(Form(" %s ",Xaxisname));
		gr1->GetYaxis()->SetTitle("RMS");
		gr1->SetMarkerColor(2);
		gr1->SetMarkerStyle(20);
	TGraphErrors *gr1gaus=new TGraphErrors(countGraph,xaxis,sigmagaus,0,errsigmagaus);
		gr1gaus->SetTitle(Form("#sigma  %s ",graphTitle));
		gr1gaus->GetXaxis()->SetTitle(Form(" %s ",Xaxisname));
		gr1gaus->GetYaxis()->SetTitle("#sigma");
		gr1gaus->SetMarkerColor(4);
		gr1gaus->SetMarkerStyle(29);
		//new TCanvas(Form("Graph1 vs  %s ",graphTitle),"RMS",200,10,800,500);
		//gr1->Draw();
	f->cd();
gr2->SetName(Form("Mean%s",hist2D->GetName()));
	gr2->Write();
	gr1->SetName(Form("RMS%s",hist2D->GetName()));
	gr1->Write();
	gr1gaus->SetName(Form("Sigma%s ",hist2D->GetName()));
	gr1gaus->Write();
}
