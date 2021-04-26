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
#include "TPad.h"
#include "TVector3.h"


void Projection(TH2F* h2D,  int range, bool rapidity, const char *name );
void ProjectionResolCorrPy(TH2F* hResol, bool rapid, const char *title);
void ProjectionNfake(TH2F* hNfake, const char *title);
void ProjectionforBin(TH2F* hist2D, const char* graphTitle , const char* Xaxisname, TCanvas *c0, TFile *f);
void ProjectionforBinMassInv(TH2F* hist2D, const char* graphTitle , const char* Xaxisname, TCanvas *c0, TFile *f);
void Armenteros(float mp[3], float mn[3], float mm[3], double v[2]);


void PlotLambdaSum(){

	TFile *hLambda=new TFile("LambdaStandAllEq-Signal-histos.root");
	TFile *hLambdaBar=new TFile("LambdaBarStandAllEq-Signal-histos.root");
/*	TFile *hD0=new TFile("histogramD0.root");*/

	TFile *foutLambda=new TFile("histogram5layerAllEq.root","RECREATE");
	//TFile *fLMoreLayer=new TFile("histogramMoreLayer50.root");
//Distribuzione 3D per efficienza
	TH3F* hYPtLzMC=(TH3F*)hLambda->Get("hYPtLzMC");
	TH3F* hYPtLzGen=(TH3F*)hLambda->Get("hYPtLzGen");
	TH3F* hYPtLzRec=(TH3F*)hLambda->Get("hYPtLzRec");

//Proiezioni degli istogrammi prima nel passaggio del rivelatore
	/*TH2F *hLenghtPtotDist = (TH2F*)hLambda->Get("hLenghtPtotDist");
  Projection(hLenghtPtotDist,3,0,"L vs Ptot");
	TH2F *hctauPtotDist = (TH2F*)hLambda->Get("hctauPtotDist");
  Projection(hctauPtotDist,3,0,"c#tau P");
	TH2F *hLzVsPt = (TH2F*)hLambda->Get("hLzVsPt");
  Projection(hLzVsPt,4,0,"Lz vs Pt");
	TH2F *hLzVsY = (TH2F*)hLambda->Get("hLzVsY");
  Projection(hLzVsY,6,1,"Lz vs Y");
	*/

	//Distribuzioni per risoluzione per projection py dopo il passaggio nel rivelatore senza correzione vertice
/*	TH2F *hResPyVsLzprot = (TH2F*)hLambda->Get("hResPyVsLzprot");
  ProjectionResolCorrPy(hResPyVsLzprot,1);
	TH2F *hResPyVsLzpion = (TH2F*)hLambda->Get("hResPyVsLzpion");
 	ProjectionResolCorrPy(hResPyVsLzpion,0);*/

	//Distribuzioni per risoluzione per projection py dopo il passaggio nel rivelatore DOPO CORREZIONE DEL VERTICE
/*	TH2F *hResPyVsLzprotcorr = (TH2F*)hLambda->Get("hResPyVsLzprotcorr");
  ProjectionResolCorrPy(hResPyVsLzprotcorr,0," ResPy Vs Lz (correction) Proton");
	TH2F *hResPyVsLzpioncorr = (TH2F*)hLambda->Get("hResPyVsLzpioncorr");
 	ProjectionResolCorrPy(hResPyVsLzpioncorr,0," ResPy Vs Lz (correction) Pion");
	TH2F *hResPyVsYprotcorr = (TH2F*)hLambda->Get("hResPyVsYprotcorr");
  ProjectionResolCorrPy(hResPyVsYprotcorr,1,"ResPy Vs Y (correction) Proton");
	TH2F *hResPyVsYpioncorr = (TH2F*)hLambda->Get("hResPyVsYpioncorr");
 	ProjectionResolCorrPy(hResPyVsYpioncorr,1,"ResPy Vs Y (correction) Pion");*/

//Distribuzioni nfake per projection py dopo il passaggio nel rivelatore DOPO CORREZIONE DEL VERTICE
	TH2F *hnfakeVsLzprot = (TH2F*)hLambda->Get("hnfakeVsLzprot");
  ProjectionNfake(hnfakeVsLzprot," Nfake Vs Lz Proton ");
	


/*//Relation between Y and Lz check andamento decrescente dell'efficienza in Lz
	TH2F* hYLz=(TH2F*)hYPtLzRec->Project3D("zx"); 
	hYLz->SetDirectory(0);
	TCanvas *yLz= new TCanvas("yLz");
	hYLz->Draw("COLZ");
	TH1D *projYLz[6];
	double mean[6]={0.},Lz[6]={0.};
	TCanvas *prjy=new TCanvas("prjy");
	prjy->Divide(3,2);
	for(int i=0;i<6;i++){
		projYLz[i]=hYLz->ProjectionX(Form("proj_y%d",i),hYLz->GetYaxis()->FindBin(i+0.01),hYLz->GetYaxis()->FindBin(i+0.99));
		prjy->cd(i+1);
		projYLz[i]->Draw();
		mean[i]=projYLz[i]->GetMean();
		Lz[i]=(i+i+1)/2.;
	}
	TGraph *gry_Lz=new TGraph(6,Lz,mean);
	TCanvas *graph=new TCanvas("graph");
	gry_Lz->SetTitle("Rapidity mean Vs Lz");
	gry_Lz->GetXaxis()->SetTitle("Lz[cm]");
	gry_Lz->GetYaxis()->SetTitle("y");
	gry_Lz->Draw("AL*");
	int minLz=0, maxLz=0;
	TH1F* hYGenRec[6], *hYMCrange[6], *Gen_GenMCYlow5[6];
	TCanvas* c1=new TCanvas("c1");
	c1->Divide(3,2);
	TCanvas *effYlow5=new TCanvas("effYlow5");
	effYlow5->Divide(3,2);
	TCanvas* c=new TCanvas("c");
	c->Divide(3,2);
	for(int j=0;j<6;j++){
		minLz=hYPtLzGen->GetZaxis()->FindBin(j+0.01);
		maxLz=hYPtLzGen->GetZaxis()->FindBin(j+0.99);
		hYPtLzGen->GetZaxis()->SetRange(minLz,maxLz);
		hYGenRec[j]=(TH1F*)hYPtLzGen->Project3D("x");
		hYGenRec[j]->SetName(Form("yGen%d",j));
		c->cd(j+1);
		hYGenRec[j]->SetDirectory(0);
		hYGenRec[j]->SetTitle(Form("Rapidity GenRec with %d<Lz<%d cm",j,j+1));
		hYGenRec[j]->Draw();
		minLz=hYPtLzMC->GetZaxis()->FindBin(j+0.01);
		maxLz=hYPtLzMC->GetZaxis()->FindBin(j+0.99);
		hYPtLzMC->GetZaxis()->SetRange(minLz,maxLz);
		hYMCrange[j]=(TH1F*)hYPtLzMC->Project3D("x");
		hYMCrange[j]->SetName("newname");

		c1->cd(j+1);
		hYMCrange[j]->SetDirectory(0);
		hYMCrange[j]->SetTitle(Form("Rapidity MC with %d<Lz<%d cm",j,j+1));
		hYMCrange[j]->Draw();

		Gen_GenMCYlow5[j] = (TH1F*)hYGenRec[j]->Clone("Gen_GenMCYlow5");
		Gen_GenMCYlow5[j]->SetDirectory(0);
		Gen_GenMCYlow5[j]->Divide(hYGenRec[j],hYMCrange[j],1,1,"B");

		effYlow5->cd(j+1);
		Gen_GenMCYlow5[j]->SetTitle(Form(" Efficiency GenRec over MC with %d<Lz<%d cm",j,j+1));
		Gen_GenMCYlow5[j]->GetXaxis()->SetTitle("y");
		Gen_GenMCYlow5[j]->GetYaxis()->SetTitle("Efficiency #epsilon");
		Gen_GenMCYlow5[j]->Draw("E1");
	}

*/



	//Definizione file dove salvare tutti gli istogrammi nati a partire dalla ntupla

	TH2F* hNfakeVszP = new TH2F("hNfakeVszP", "Nfake vs zP", 6, -0.5, 5.5, 100, 0, 50);
	TH2F *hMassinreczP = new TH2F("hMassinrecVszP", "Invarian Mass vs zP", 3000, 0,3,100, 0, 50);
	TH2F *hMassinrecY = new TH2F("hMassinvVsY", "Invariant Mass vs Y", 3000, 0,3,60, 0, 6);
	TH2F *hMassinrecSwitchY = new TH2F("hMassinrecSwitchVsY", "Mass Switch vs Y", 3000, 0,3,60, 0, 6);
	TH2F *hMassinrecY0fake = new TH2F("hMassinvVsY0fake", "Invariant Mass vs Y 0fake", 3000, 0,3,60, 0, 6);
	TH2F *hMassinrecSwitchY0fake = new TH2F("hMassinrecSwitchVsY0fake", "Mass Switch vs Y 0fake", 3000, 0,3,60, 0, 6);
	/*TH2F *hResVertVsY[3], *hResVertVsLz[3];
	hResVertVsY[0] = new TH2F("hResVxVsY", "Res Vx vs Y", 200, -200., 200., 60, 0, 6);
	hResVertVsY[1] = new TH2F("hResVyVsY", "Res Vy vs Y", 200, -200., 200., 60, 0, 6);
	hResVertVsY[2] = new TH2F("hResVzVsY", "Res Vz vs Y", 200, -3000., 3000., 60, 0, 6);
	hResVertVsLz[0] = new TH2F("hResVxVsLz", "Res Vx vs Lz", 200, -200., 200., 200, 0, 60);
	hResVertVsLz[1] = new TH2F("hResVyVsLz", "Res Vy vs Lz", 200, -200., 200., 200, 0, 60);
	hResVertVsLz[2] = new TH2F("hResVzVsLz", "Res Vz vs Lz", 200, -2000., 2000., 200, 0, 60);*/
	TH2F *hResVxVsY = new TH2F("hResVxVsY", "Res Vx vs Y", 200, -200., 200., 60, 0, 6);
  	TH2F *hResVyVsY = new TH2F("hResVyVsY", "Res Vy vs Y", 200, -200., 200., 60, 0, 6);
  	TH2F *hResVzVsY = new TH2F("hResVzVsY", "Res Vz vs Y", 200, -3000., 3000., 60, 0, 6);
	TH2F *hResVxVsLz = new TH2F("hResVxVsLz", "Res Vx vs Lz", 200, -200., 200., 200, 0, 60);
  	TH2F *hResVyVsLz = new TH2F("hResVyVsLz", "Res Vy vs Lz", 200, -200., 200., 200, 0, 60);
  	TH2F *hResVzVsLz = new TH2F("hResVzVsLz", "Res Vz vs Lz", 200, -2000., 2000., 200, 0, 60);
	TH2F *hResPVsYprot[3], *hResPVsYpion[3];
	hResPVsYprot[0] = new TH2F("hResPxVsYprot", "Res Px vs Y", 200, -0.2, 0.2, 60, 0, 6);
	hResPVsYprot[1] = new TH2F("hResPyVsYprot", "Res Py vs Y", 200, -0.2, 0.2, 60, 0, 6);
	hResPVsYprot[2] = new TH2F("hResPzVsYprot", "Res Pz vs Y", 200, -0.2, 0.2, 60, 0, 6);
	hResPVsYpion[0] = new TH2F("hResPxVsYpion", "Res Px vs Y", 200, -0.08, 0.08, 60, 0, 6);
	hResPVsYpion[1] = new TH2F("hResPyVsYpion", "Res Py vs Y", 200, -0.08, 0.08, 60, 0, 6);
	hResPVsYpion[2] = new TH2F("hResPzVsYpion", "Res Pz vs Y", 200, -0.08, 0.08, 60, 0, 6);
	TH2F *hResPVsYpionRel[3], *hResPVsYprotRel[3];
	hResPVsYpionRel[0] = new TH2F("hResPxVsYpionRel", "Res Px vs Y Relativa", 200, -0.08, 0.08, 60, 0, 6);
	hResPVsYpionRel[1] = new TH2F("hResPyVsYpionRel", "Res Py vs Y Relativa", 200, -0.08, 0.08, 60, 0, 6);
	hResPVsYpionRel[2] = new TH2F("hResPzVsYpionRel", "Res Pz vs Y Relativa", 200, -0.08, 0.08, 60, 0, 6);
	hResPVsYprotRel[0] = new TH2F("hResPxVsYprotRel", "Res Px vs Y Relativa", 200, -0.08, 0.08, 60, 0, 6);
	hResPVsYprotRel[1] = new TH2F("hResPyVsYprotRel", "Res Py vs Y Relativa", 200, -0.08, 0.08, 60, 0, 6);
	hResPVsYprotRel[2] = new TH2F("hResPzVsYprotRel", "Res Pz vs Y Relativa", 200, -0.08, 0.08, 60, 0, 6);
	TH2F *hResPVsLzprot[3], *hResPVsLzpion[3];
	hResPVsLzprot[0] = new TH2F("hResPxVsLzprot", "Res Px vs Lz", 200, -1., 1., 200, 0, 60);
	hResPVsLzprot[1] = new TH2F("hResPyVsLzprot", "Res Py vs Lz", 200, -1., 1., 200, 0, 60);
	hResPVsLzprot[2] = new TH2F("hResPzVsLzprot", "Res Pz vs Lz", 200, -1., 1., 200, 0, 60);
	hResPVsLzpion[0] = new TH2F("hResPxVsLzpion", "Res Px vs Lz", 200, -1., 1., 200, 0, 60);
	hResPVsLzpion[1] = new TH2F("hResPyVsLzpion", "Res Py vs Lz", 200, -1., 1., 200, 0, 60);
	hResPVsLzpion[2] = new TH2F("hResPzVsLzpion", "Res Pz vs Lz", 200, -1., 1., 200, 0, 60);
	TH1F *hPprot[3], *hPpion[3],*hPprotFake[3], *hPpionFake[3];
	hPprot[0] = new TH1F("hPxprot", " Px proton ", 200, -5., 5.);
	hPprot[1] = new TH1F("hPyprot", " Py proton", 200, -2., 2.);
	hPprot[2] = new TH1F("hPzprot", " Pz proton ", 200, -5., 50.);
	hPpion[0] = new TH1F("hPxpion", " Px pion", 200, -5., 5.);
	hPpion[1] = new TH1F("hPypion", " Py pion", 200, -2., 2.);
	hPpion[2] = new TH1F("hPzpion", " Pz pion", 200, -1., 20.);
	hPprotFake[0] = new TH1F("hPxprotFake", " Px proton Fake", 200, -5., 5.);
	hPprotFake[1] = new TH1F("hPyprotFake", " Py proton Fake", 200, -2., 2.);
	hPprotFake[2] = new TH1F("hPzprotFake", " Pz proton Fake", 200, -0., 50.);
	hPpionFake[0] = new TH1F("hPxpionFake", " Px pion Fake", 200, -5., 5.);
	hPpionFake[1] = new TH1F("hPypionFake", " Py pion Fake", 200, -2., 2.);
	hPpionFake[2] = new TH1F("hPzpionFake", " Pz pion Fake", 200, -0., 20.);

	TH2F *hPtprotfakeVsyLambda = new TH2F("hPtprotfakeVsyLambda", " Pt proton fake vs y #Lambda", 60, 0, 6, 200, -1., 10.);
	TH2F *hPtpionfakeVsyLambda = new TH2F("hPtpionfakeVsyLambda", " Pt pion fake vs y #Lambda", 60, 0, 6,200, -1., 10.);

	TH2F *hPzVsYprot = new TH2F("hPzVsYprot", "Pz vs Y proton ", 60, 0, 6, 200, -1, 60.);
	TH2F *hPzVsYpion = new TH2F("hPzVsYpion", "Pz vs Y pion", 60, 0, 6, 200, -1, 60.);
	TH2F *hPzprotVsPzLambda = new TH2F("hPzprotVsPzLambda", "Pz proton vs Pz Lambda  ", 200, -1, 60., 200, -1, 60.);
	TH2F *hPzpionVsPzLambda = new TH2F("hPzpionVsPzLambda", "Pz pion vs Pz Lambda ", 200, -1, 60., 200, -1, 60.);


	TH1F *htheta = new TH1F("htheta", "Angle between the particles Fake", 200, 0., 7.);

	TH2F* hNfakeProtVsChi2 = new TH2F("hNfakeProtVsChi2", "NfakeProt vs Chi2", 6, -0.5, 5.5, 100, 0, 2);
	TH2F* hNfakePionVsChi2 = new TH2F("hNfakePionVsChi2", "NfakePion vs Chi2", 6, -0.5, 5.5, 100, 0, 2);
	TH1F* hMassinrecnfakeprot0 = new TH1F("hMassinrecnfakeprot", "Massinv for 0 fake prot", 3000, 0,3);
	TH1F* hMassinrecnfakeprot1 = new TH1F("hMassinrecnfakeprot1", "Massinv for 1 fake prot", 3000, 0,3);
	TH1F* hMassinrecnfakeprot2 = new TH1F("hMassinrecnfakeprot2", "Massinv for 2 fake prot", 3000, 0,3);
	TH1F* hMassinrecnfakeprot3 = new TH1F("hMassinrecnfakeprot3", "Massinv for 3 fake prot", 3000, 0,3);
	TH1F* hMassinrecnfakeprot4 = new TH1F("hMassinrecnfakeprot4", "Massinv for 4 fake prot", 3000, 0,3);
	TH1F* hMassinrecnfakeprot5= new TH1F("hMassinrecnfakeprot5", "Massinv for 5 fake prot", 3000, 0,3);
	TH1F* hDCAnfakeprot0 = new TH1F("hDCAnfakeprot0", "DCA for 0 fake prot", 100, 0,0.2);
	TH1F* hDCAnfakeprot4 = new TH1F("hDCAnfakeprot4", "DCA for 4 fake prot", 100, 0,0.1);


	TH2F *hLzrecVsLzgen = new TH2F("hLzrecVsLzgen", "Lz rec vs Lz gen", 200, 0, 60, 200, 0, 60);
	TH2F *hDeltaLzVsLzgen = new TH2F("hDeltaLzVsLzgen", "#Delta Lz vs Lz gen", 200, 0, 60, 600, -1, 1);

	TH2F* hArmen= new TH2F("hArmen","Armenteros-Podolansky Plot ",200,-1,1,200,0,1);
	TH2F* hArmenLambda= new TH2F("hArmenLambda","Armenteros Plot for #Lambda ",200,-1,1,200,0,1);
	TH2F* hArmenLambdaFake= new TH2F("hArmenLambdaFake","Armenteros Plot for #Lambda Fake ",200,-1,1,200,0,1);
	TH1F* hMArm = new TH1F("hMArm","Invariant Mass for #alpha=0 #Lambda",200,1,1.8);
	TH1F* hyLambda = new TH1F("hyLambda","Rapidity for #Lambda",200,0,6);
	TH1F* hyLambdaFake = new TH1F("hyLambdaFake","Rapidity for #Lambda with fake",200,0,6);
	TH1F* hptLambda = new TH1F("hptLambda","P_{t} for #Lambda",200,0,6);
	TH1F* hptLambdaFake = new TH1F("hptLambdaFake","P_{t} for #Lambda with fake",200,0,6);
	TH1F* hPtProtFake = new TH1F("hPtProtFake","P_{t} proton fake",200,0,10);
	TH1F* hPtPionFake = new TH1F("hPtPionFake","P_{t} pion fake",200,0,10);


	TH1F* hChi2Prot[2];
	hChi2Prot[0] = new TH1F("hChi2Prot", "Chi2 proton", 3000, 0,3);
	hChi2Prot[1] = new TH1F("hChi2ProtFake", "Chi2 proton ONLY for fake", 3000, 0,3);
	TH1F* hChi2ITSProt[2];
	hChi2ITSProt[0] = new TH1F("hChi2ITSProt", "Chi2ITS proton", 3000, 0,3);
	hChi2ITSProt[1] = new TH1F("hChi2ITSProtFake", "Chi2ITS proton ONLY for fake", 3000, 0,3);
	TH1F* hChi2Pion[2];
	hChi2Pion[0] = new TH1F("hChi2Pion", "Chi2 pion", 3000, 0,3);
	hChi2Pion[1] = new TH1F("hChi2PionFake", "Chi2 pion ONLY for fake", 3000, 0,3);
	TH1F* hChi2ITSPion[2];
	hChi2ITSPion[0] = new TH1F("hChi2ITSPion", "Chi2ITS pion", 3000, 0,3);
	hChi2ITSPion[1] = new TH1F("hChi2ITSPionFake", "Chi2ITS pion ONLY for fake", 3000, 0,3);

	TH1F* hthetaProtFake = new TH1F("hthetaProtFake","Theta Proton fake",300,0,90);
	TH1F* hthetaPionFake = new TH1F("hthetaPionFake","Theta Pion fake",300,0,90);
	TH1F* hthetaProt = new TH1F("hthetaProt","Theta Proton ",300,0,90);
	TH1F* hthetaPion = new TH1F("hthetaPion","Theta Pion ",300,0,90);

	//Apertura file di input con variabili lambda per studio performance
	//Lettura Ntupla 
 	TFile filin("ntuplaStandAllEq.root");
  	float nfaketrkprot, nfaketrkpion,yrap,xrec,yrec,zrec,secvertgenProt[3],massinvrec,pGenProt[3],pGenPion[3],dca,pRecProt[3],pRecPion[3],chi2prot,chi2ITSprot,chi2pion,chi2ITSpion,massinvSwitch,pGenLambda[3];
  	TNtuple *variables = (TNtuple*)filin.Get("nt");
		variables->SetBranchAddress("nfaketrkprot",&nfaketrkprot);
    	variables->SetBranchAddress("nfaketrkpion",&nfaketrkpion);
    	variables->SetBranchAddress("ygen",&yrap);
    	variables->SetBranchAddress("xP",&xrec);
		variables->SetBranchAddress("yP",&yrec);
		variables->SetBranchAddress("zP",&zrec);
		variables->SetBranchAddress("Vxgen",&secvertgenProt[0]);
		variables->SetBranchAddress("Vygen",&secvertgenProt[1]);
		variables->SetBranchAddress("Vzgen",&secvertgenProt[2]);
		variables->SetBranchAddress("massinvrec",&massinvrec);
		variables->SetBranchAddress("pxgenprot",&pGenProt[0]);
		variables->SetBranchAddress("pygenprot",&pGenProt[1]);
		variables->SetBranchAddress("pzgenprot",&pGenProt[2]);
		variables->SetBranchAddress("pxgenpion",&pGenPion[0]);
		variables->SetBranchAddress("pygenpion",&pGenPion[1]);
		variables->SetBranchAddress("pzgenpion",&pGenPion[2]);
		variables->SetBranchAddress("dca",&dca);
		variables->SetBranchAddress("chi2prot",&chi2prot);
		variables->SetBranchAddress("chi2ITSprot",&chi2ITSprot);
		variables->SetBranchAddress("chi2pion",&chi2pion);
		variables->SetBranchAddress("chi2ITSpion",&chi2ITSpion);
		variables->SetBranchAddress("pxrecprot",&pRecProt[0]);
		variables->SetBranchAddress("pyrecprot",&pRecProt[1]);
		variables->SetBranchAddress("pzrecprot",&pRecProt[2]);
		variables->SetBranchAddress("pxrecpion",&pRecPion[0]);
		variables->SetBranchAddress("pyrecpion",&pRecPion[1]);
		variables->SetBranchAddress("pzrecpion",&pRecPion[2]);
		variables->SetBranchAddress("massinvswitch",&massinvSwitch);
		variables->SetBranchAddress("pxlambda",&pGenLambda[0]);
		variables->SetBranchAddress("pylambda",&pGenLambda[1]);
		variables->SetBranchAddress("pzlambda",&pGenLambda[2]);


	int events=variables->GetEntries();
	double thetaProt=0, thetaPion=0;
	for(int i=0 ; i<events; i++) {
		variables->GetEvent(i);
		double residVx=10000.*(xrec - secvertgenProt[0]);
    	double residVy=10000.*(yrec - secvertgenProt[1]);
    	double residVz=10000.*(zrec - secvertgenProt[2]);
		double val[2]={0};
		
		double pTotProt=sqrt((pRecProt[0]*pRecProt[0])+(pRecProt[1]*pRecProt[1])+(pRecProt[2]*pRecProt[2]));
		double pTotPion=sqrt(pRecPion[0]*pRecPion[0]+pRecPion[1]*pRecPion[1]+pRecPion[2]*pRecPion[2]);

		hyLambda->Fill(yrap);
		hptLambda->Fill(sqrt(pGenLambda[0]*pGenLambda[0]+pGenLambda[1]*pGenLambda[1]));
		double p1abs=sqrt(pRecProt[0]*pRecProt[0]+pRecProt[1]*pRecProt[1]+pRecProt[2]*pRecProt[2]);
		thetaProt=TMath::ACos(pRecProt[2]/p1abs);
		hthetaProt->Fill((thetaProt*180)/TMath::Pi());
		double p2abs=sqrt(pRecPion[0]*pRecPion[0]+pRecPion[1]*pRecPion[1]+pRecPion[2]*pRecPion[2]);
		thetaPion=TMath::ACos(pRecPion[2]/p2abs);
		hthetaPion->Fill((thetaPion*180)/TMath::Pi());

		Armenteros(pRecProt, pRecPion, pGenLambda,val);
		hArmenLambda->Fill(val[0],val[1]);
		hArmen->Fill(val[0],val[1]);
		if(zrec<=20.){
			
			hNfakeVszP->Fill(nfaketrkprot, zrec);
			}

	
		hMassinrecY->Fill(massinvrec, yrap);
		hMassinrecSwitchY->Fill(massinvSwitch, yrap);

		hNfakeProtVsChi2->Fill(nfaketrkprot,chi2prot);
		hNfakePionVsChi2->Fill(nfaketrkpion,chi2pion);


			//hLzrecVsLzgen->Fill(secvertgenProt[2],zrec);
			//hDeltaLzVsLzgen->Fill(secvertgenProt[2],zrec-secvertgenProt[2]);
		
		if(nfaketrkprot == 0 && nfaketrkpion == 0) {
			if(zrec<=20.){
				hMassinreczP->Fill(massinvrec, zrec);
				}
			hMassinrecY0fake->Fill(massinvrec, yrap);
			hMassinrecSwitchY0fake->Fill(massinvSwitch, yrap);
			hChi2Prot[0]->Fill(chi2prot);
			hChi2Pion[0]->Fill(chi2pion);
			hChi2ITSProt[0]->Fill(chi2ITSprot);
			hChi2ITSPion[0]->Fill(chi2ITSpion);
			
				/*double p1abs=sqrt(pRecProt[0]*pRecProt[0]+pRecProt[1]*pRecProt[1]+pRecProt[2]*pRecProt[2]);
				double p2abs=sqrt(pRecPion[0]*pRecPion[0]+pRecPion[1]*pRecPion[1]+pRecPion[2]*pRecPion[2]);
				theta=TMath::ACos((pRecProt[0]*pRecPion[0]+pRecProt[1]*pRecPion[1]+pRecProt[2]*pRecPion[2])/(p1abs*p2abs));
				htheta->Fill(theta);*/
			
			/*double p1abs=sqrt(pRecProt[0]*pRecProt[0]+pRecProt[1]*pRecProt[1]+pRecProt[2]*pRecProt[2]);
		thetaProt=TMath::ACos(pRecProt[2]/p1abs);
		hthetaProt->Fill(thetaProt);
		double p2abs=sqrt(pRecPion[0]*pRecPion[0]+pRecPion[1]*pRecPion[1]+pRecPion[2]*pRecPion[2]);
		thetaPion=TMath::ACos(pRecPion[2]/p2abs);
		hthetaPion->Fill(thetaPion);*/

			hPzVsYpion->Fill(yrap,pRecPion[2]);
			hPzVsYprot->Fill(yrap,pRecProt[2]);

			hPzprotVsPzLambda->Fill(pGenLambda[2],pGenProt[2]);
			hPzpionVsPzLambda->Fill(pGenLambda[2],pGenPion[2]);

			for(int i=0;i<3;i++){
				//Resolution on P
				hResPVsYpionRel[i]->Fill((pRecPion[i]-pGenPion[i])/pGenPion[i],yrap);
				hResPVsYprotRel[i]->Fill((pRecProt[i]-pGenProt[i])/pGenProt[i],yrap);
				hResPVsYprot[i]->Fill(pRecProt[i]-pGenProt[i],yrap);
				hResPVsYpion[i]->Fill(pRecPion[i]-pGenPion[i],yrap);
				hResPVsLzprot[i]->Fill(pRecProt[i]-pGenProt[i],zrec);
				hResPVsLzpion[i]->Fill(pRecPion[i]-pGenPion[i],zrec);
				hPprot[i]->Fill(pRecProt[i]);
				hPpion[i]->Fill(pRecPion[i]);
			}

			//Resolution on Vertex
			hResVxVsY->Fill(residVx, yrap);
			hResVyVsY->Fill(residVy, yrap);
			hResVzVsY->Fill(residVz, yrap);
			hResVxVsLz->Fill(residVx, zrec);
			hResVyVsLz->Fill(residVy, zrec);
			hResVzVsLz->Fill(residVz, zrec);

			hLzrecVsLzgen->Fill(secvertgenProt[2],zrec);
			hDeltaLzVsLzgen->Fill(secvertgenProt[2],zrec-secvertgenProt[2]);
		}

		if(nfaketrkprot!=0 && nfaketrkpion!=0){
			hPtprotfakeVsyLambda->Fill(yrap,sqrt((pRecProt[0]*pRecProt[0])+(pRecProt[1]*pRecProt[1])));
			hPtpionfakeVsyLambda->Fill(yrap,sqrt((pRecPion[0]*pRecPion[0])+(pRecPion[1]*pRecPion[1])));
			val[0]=0;
			val[1]=0;
			Armenteros(pRecProt, pRecPion, pGenLambda,val);
			hArmenLambdaFake->Fill(val[0],val[1]);
			hMArm->Fill(massinvrec);
			hyLambdaFake->Fill(yrap);
			hPtProtFake->Fill(sqrt((pRecProt[0]*pRecProt[0])+(pRecProt[1]*pRecProt[1])));
			hPtPionFake->Fill(sqrt((pRecPion[0]*pRecPion[0])+(pRecPion[1]*pRecPion[1])));
			hptLambdaFake->Fill(sqrt(pGenLambda[0]*pGenLambda[0]+pGenLambda[1]*pGenLambda[1]));
			for(int i=0;i<3;i++){
				hPprotFake[i]->Fill(pRecProt[i]);
				hPpionFake[i]->Fill(pRecPion[i]);
			}
			if((val[0]<0.1 && val[0]>-0.1) && (massinvrec>1.26 && massinvrec<1.38)){
				/*printf("#alpha %f , massa invariante %f \n", val[0],massinvrec);
				printf("protone: momRec %f %f %f ;  numfake=%f \n",pRecProt[0],pRecProt[1],pRecProt[2],nfaketrkprot);
				printf("pione: momRec %f %f %f ;  numfake=%f \n",pRecPion[0],pRecPion[1],pRecPion[2],nfaketrkpion);*/
				hChi2Prot[1]->Fill(chi2prot);
				hChi2Pion[1]->Fill(chi2pion);
				hChi2ITSProt[1]->Fill(chi2ITSprot);
				hChi2ITSPion[1]->Fill(chi2ITSpion);
				double p1abs=sqrt(pRecProt[0]*pRecProt[0]+pRecProt[1]*pRecProt[1]+pRecProt[2]*pRecProt[2]);
				thetaProt=TMath::ACos(pRecProt[2]/p1abs);
				hthetaProtFake->Fill((thetaProt*180)/TMath::Pi());
				double p2abs=sqrt(pRecPion[0]*pRecPion[0]+pRecPion[1]*pRecPion[1]+pRecPion[2]*pRecPion[2]);
				thetaPion=TMath::ACos(pRecPion[2]/p2abs);
				hthetaPionFake->Fill((thetaPion*180)/TMath::Pi());

			}
		}
		if(nfaketrkprot == 0) {
			hMassinrecnfakeprot0->Fill(massinvrec);
			hMassinrecnfakeprot0->SetTitle(Form("Invariant Mass for %f faketrk proton",nfaketrkprot));
			hDCAnfakeprot0->Fill(dca);
			hDCAnfakeprot0->SetTitle(Form("DCA for %f faketrk proton",nfaketrkprot));	
		}
		if(nfaketrkprot == 1){
			hMassinrecnfakeprot1->Fill(massinvrec);
			hMassinrecnfakeprot1->SetTitle(Form("Invariant Mass for %f faketrk proton",nfaketrkprot));
		}
		if(nfaketrkprot == 2){
			hMassinrecnfakeprot2->Fill(massinvrec);
			hMassinrecnfakeprot2->SetTitle(Form("Invariant Mass for %f faketrk proton",nfaketrkprot));
		}
		if(nfaketrkprot == 3) {
			hMassinrecnfakeprot3->Fill(massinvrec);
			hMassinrecnfakeprot3->SetTitle(Form("Invariant Mass for %f faketrk proton",nfaketrkprot));
		}			
		if(nfaketrkprot == 4) {
			hMassinrecnfakeprot4->Fill(massinvrec);
			hMassinrecnfakeprot4->SetTitle(Form("Invariant Mass for %f faketrk proton",nfaketrkprot));
			hDCAnfakeprot4->Fill(dca);
			hDCAnfakeprot4->SetTitle(Form("DCA for %f faketrk proton",nfaketrkprot));
		}
		if(nfaketrkprot == 5) {
			hMassinrecnfakeprot5->Fill(massinvrec);
			hMassinrecnfakeprot5->SetTitle(Form("Invariant Mass for %f faketrk proton",nfaketrkprot));
		}
	}
	filin.Close();

	TH2F *hPtprotfakeVsyLambdaBar = new TH2F("hPxprotfakeVsyLambdaBar", " Pt proton fake vs y #bar{#Lambda}", 60, 0, 6, 200, -1., 10.);
	TH2F *hPtpionfakeVsyLambdaBar = new TH2F("hPxpionfakeVsyLambdaBar", " Pt pion fake vs y #bar{#Lambda}", 60, 0, 6,200, -1., 10.);
	TH2F* hArmenLambdaBar= new TH2F("hArmenLambdaBar","Armenteros Plot for #bar{#Lambda}",200,-1,1,200,0,0.8);
	TH2F* hArmenLambdaBarFake= new TH2F("hArmenLambdaBarFake","Armenteros Plot for #bar{#Lambda} Fake",200,-1,1,200,0,0.8);
	TH1F* hMArmBar = new TH1F("hMArmBar","Invariant Mass for #alpha=0 #bar{#Lambda}",200,1,1.8);
	TH1F* hyLambdaBar = new TH1F("hyLambdaBar","Rapidity for #bar{#Lambda}",200,0,6);
	TH1F* hyLambdaBarFake = new TH1F("hyLambdaBarFake","Rapidity for #bar{#Lambda} with fake",200,0,6);
	TH1F* hptLambdaBar = new TH1F("hptLambdaBar","P_{t} for #bar{#Lambda}",200,0,6);
	TH1F* hptLambdaBarFake = new TH1F("hptLambdaBarFake","P_{t} for #bar{#Lambda} with fake",200,0,6);
	TH1F* hPtProtFakeBar = new TH1F("hPtProtFakeBar","P_{t} proton fake",200,0,10);
	TH1F* hPtPionFakeBar = new TH1F("hPtPionFakeBar","P_{t} pion fake",200,0,10);

	TH1F* hPprotFakeBar[3], *hPpionFakeBar[3];
	hPprotFakeBar[0] = new TH1F("hPxprotFakeBar", " Px proton Fake Bar", 200, -5., 5.);
	hPprotFakeBar[1] = new TH1F("hPyprotFakeBar", " Py proton Fake Bar", 200, -2., 2.);
	hPprotFakeBar[2] = new TH1F("hPzprotFakeBar", " Pz proton Fake Bar", 200, 0, 50.);
	hPpionFakeBar[0] = new TH1F("hPxpionFakeBar", " Px pion Fake Bar", 200, -5., 5.);
	hPpionFakeBar[1] = new TH1F("hPypionFakeBar", " Py pion Fake Bar", 200, -2., 2.);
	hPpionFakeBar[2] = new TH1F("hPzpionFakeBar", " Pz pion Fake Bar", 200, 0., 20.);

	TH1F* hthetaProtBarFake = new TH1F("hthetaProtBarFake","Theta Proton fake Bar",200,0,1);
	TH1F* hthetaPionBarFake = new TH1F("hthetaPionBarFake","Theta Pion fake Bar",200,0,1);
	TH1F* hthetaProtBar = new TH1F("hthetaProtBar","Theta Proton Bar ",200,0,1);
	TH1F* hthetaPionBar = new TH1F("hthetaPionBar","Theta Pion Bar",200,0,1);


	TFile filinBar("ntuplaLambdaBarStandAllEq.root");
  	float nfaketrkprotBar, nfaketrkpionBar,yrapBar,xrecBar,yrecBar,zrecBar,secvertgenProtBar[3],massinvrecBar,pGenProtBar[3],pGenPionBar[3],dcaBar,pRecProtBar[3],pRecPionBar[3],chi2protBar,chi2ITSprotBar,chi2pionBar,chi2ITSpionBar,massinvSwitchBar,pGenLambdaBar[3];
  	TNtuple *variablesBar = (TNtuple*)filinBar.Get("ntLambdaBarstand");
		variablesBar->SetBranchAddress("nfaketrkprot",&nfaketrkprotBar);
    	variablesBar->SetBranchAddress("nfaketrkpion",&nfaketrkpionBar);
    	variablesBar->SetBranchAddress("ygen",&yrapBar);
    	variablesBar->SetBranchAddress("xP",&xrecBar);
		variablesBar->SetBranchAddress("yP",&yrecBar);
		variablesBar->SetBranchAddress("zP",&zrecBar);
		variablesBar->SetBranchAddress("Vxgen",&secvertgenProtBar[0]);
		variablesBar->SetBranchAddress("Vygen",&secvertgenProtBar[1]);
		variablesBar->SetBranchAddress("Vzgen",&secvertgenProtBar[2]);
		variablesBar->SetBranchAddress("massinvrec",&massinvrecBar);
		variablesBar->SetBranchAddress("pxgenprot",&pGenProtBar[0]);
		variablesBar->SetBranchAddress("pygenprot",&pGenProtBar[1]);
		variablesBar->SetBranchAddress("pzgenprot",&pGenProtBar[2]);
		variablesBar->SetBranchAddress("pxgenpion",&pGenPionBar[0]);
		variablesBar->SetBranchAddress("pygenpion",&pGenPionBar[1]);
		variablesBar->SetBranchAddress("pzgenpion",&pGenPionBar[2]);
		variablesBar->SetBranchAddress("dca",&dcaBar);
		variablesBar->SetBranchAddress("chi2prot",&chi2protBar);
		variablesBar->SetBranchAddress("chi2ITSprot",&chi2ITSprotBar);
		variablesBar->SetBranchAddress("chi2pion",&chi2pionBar);
		variablesBar->SetBranchAddress("chi2ITSpion",&chi2ITSpionBar);
		variablesBar->SetBranchAddress("pxrecprot",&pRecProtBar[0]);
		variablesBar->SetBranchAddress("pyrecprot",&pRecProtBar[1]);
		variablesBar->SetBranchAddress("pzrecprot",&pRecProtBar[2]);
		variablesBar->SetBranchAddress("pxrecpion",&pRecPionBar[0]);
		variablesBar->SetBranchAddress("pyrecpion",&pRecPionBar[1]);
		variablesBar->SetBranchAddress("pzrecpion",&pRecPionBar[2]);
		variablesBar->SetBranchAddress("massinvswitch",&massinvSwitchBar);
		variablesBar->SetBranchAddress("pxlambda",&pGenLambdaBar[0]);
		variablesBar->SetBranchAddress("pylambda",&pGenLambdaBar[1]);
		variablesBar->SetBranchAddress("pzlambda",&pGenLambdaBar[2]);

	double thetaProtBar=0, thetaPionBar=0;
	int eventsBar=variablesBar->GetEntries();

	for(int i=0 ; i<eventsBar; i++) {
		variablesBar->GetEvent(i);
		double residVxBar=10000.*(xrecBar - secvertgenProtBar[0]);
		double residVyBar=10000.*(yrecBar - secvertgenProtBar[1]);
		double residVzBar=10000.*(zrecBar - secvertgenProtBar[2]);

		hyLambdaBar->Fill(yrapBar);
		hptLambdaBar->Fill(sqrt(pGenLambdaBar[0]*pGenLambdaBar[0]+pGenLambdaBar[1]*pGenLambdaBar[1]));

		double valBar[2]={0};
		Armenteros(pRecPionBar, pRecProtBar, pGenLambdaBar,valBar);
		hArmenLambdaBar->Fill(valBar[0],valBar[1]);
		hArmen->Fill(valBar[0],valBar[1]);
		double pTotProtBar=sqrt(pRecProtBar[0]*pRecProtBar[0]+pRecProtBar[1]*pRecProtBar[1]+pRecProtBar[2]*pRecProtBar[2]);
		double pTotPionBar=sqrt(pRecPionBar[0]*pRecPionBar[0]+pRecPionBar[1]*pRecPionBar[1]+pRecPionBar[2]*pRecPionBar[2]);

		if(zrecBar<=20.){
			hNfakeVszP->Fill(nfaketrkprotBar, zrecBar);
		}

	/*	double p1absBar=sqrt(pRecProtBar[0]*pRecProtBar[0]+pRecProtBar[1]*pRecProtBar[1]+pRecProtBar[2]*pRecProtBar[2]);
		double p2absBar=sqrt(pRecPionBar[0]*pRecPionBar[0]+pRecPionBar[1]*pRecPionBar[1]+pRecPionBar[2]*pRecPionBar[2]);
		theta=TMath::ACos((pRecProtBar[0]*pRecPionBar[0]+pRecProtBar[1]*pRecPionBar[1]+pRecProtBar[2]*pRecPionBar[2])/(p1absBar*p2absBar));

		htheta->Fill(theta);*/
		
		hMassinrecY->Fill(massinvrecBar, yrapBar);
		hMassinrecSwitchY->Fill(massinvSwitchBar, yrapBar);
		hNfakeProtVsChi2->Fill(nfaketrkprotBar,chi2protBar);
		hNfakePionVsChi2->Fill(nfaketrkpionBar,chi2pionBar);

		//hLzrecVsLzgen->Fill(secvertgenProtBar[2],zrecBar);
		//hDeltaLzVsLzgen->Fill(secvertgenProtBar[2],zrecBar-secvertgenProtBar[2]);

		if(nfaketrkprotBar == 0 && nfaketrkpionBar == 0) {
					if(zrecBar<=20.){
			hMassinreczP->Fill(massinvrecBar, zrecBar);
					}
			hMassinrecY0fake->Fill(massinvrecBar, yrapBar);
			hMassinrecSwitchY0fake->Fill(massinvSwitchBar, yrapBar);

			hPzVsYpion->Fill(yrapBar,pRecPionBar[2]);
			hPzVsYprot->Fill(yrapBar,pRecProtBar[2]);

			//hPzprotVsPzLambda->Fill(pGenLambdaBar[2],pGenProtBar[2]);
			//hPzpionVsPzLambda->Fill(pGenLambdaBar[2],pGenPionBar[2]);
				double p1abs=sqrt(pRecProtBar[0]*pRecProtBar[0]+pRecProtBar[1]*pRecProtBar[1]+pRecProtBar[2]*pRecProtBar[2]);
				thetaProtBar=TMath::ACos(pRecProtBar[2]/p1abs);
				hthetaProtBar->Fill(thetaProtBar);
				double p2abs=sqrt(pRecPionBar[0]*pRecPionBar[0]+pRecPionBar[1]*pRecPionBar[1]+pRecPionBar[2]*pRecPionBar[2]);
				thetaPionBar=TMath::ACos(pRecPionBar[2]/p2abs);
				hthetaPionBar->Fill(thetaPionBar);
			for(int i=0;i<3;i++){
				//Resolution on P
				hResPVsYpionRel[i]->Fill((pRecPionBar[i]-pGenPionBar[i])/pGenPionBar[i],yrapBar);
				hResPVsYprotRel[i]->Fill((pRecProtBar[i]-pGenProtBar[i])/pGenProtBar[i],yrapBar);
				hResPVsYprot[i]->Fill(pRecProtBar[i]-pGenProtBar[i],yrapBar);
				hResPVsYpion[i]->Fill(pRecPionBar[i]-pGenPionBar[i],yrapBar);
				hResPVsLzprot[i]->Fill(pRecProtBar[i]-pGenProtBar[i],zrecBar);
				hResPVsLzpion[i]->Fill(pRecPionBar[i]-pGenPionBar[i],zrecBar);
				//hPprot[i]->Fill(pRecProtBar[i]);
				//hPpion[i]->Fill(pRecPionBar[i]);
			}

			//Resolution on vertex
			hResVxVsY->Fill(residVxBar, yrapBar);
			hResVyVsY->Fill(residVyBar, yrapBar);
			hResVzVsY->Fill(residVzBar, yrapBar);
			hResVxVsLz->Fill(residVxBar, zrecBar);
			hResVyVsLz->Fill(residVyBar, zrecBar);
			hResVzVsLz->Fill(residVzBar, zrecBar);
			hLzrecVsLzgen->Fill(secvertgenProtBar[2],zrecBar);
			hDeltaLzVsLzgen->Fill(secvertgenProtBar[2],zrecBar-secvertgenProtBar[2]);
		}

		if(nfaketrkprotBar!=0 && nfaketrkpionBar!=0){
			hPtprotfakeVsyLambdaBar->Fill(yrapBar,sqrt((pRecProtBar[0]*pRecProtBar[0])+(pRecProtBar[1]*pRecProtBar[1])));
			hPtpionfakeVsyLambdaBar->Fill(yrapBar,sqrt((pRecPionBar[0]*pRecPionBar[0])+(pRecPionBar[1]*pRecPionBar[1])));
			valBar[0]=0;
			valBar[1]=0;
			Armenteros(pRecPionBar, pRecProtBar, pGenLambdaBar, valBar);
			hArmenLambdaBarFake->Fill(valBar[0],valBar[1]);
			hMArmBar->Fill(massinvrecBar);
			hyLambdaBarFake->Fill(yrapBar);
			hPtProtFakeBar->Fill(sqrt((pRecProtBar[0]*pRecProtBar[0])+(pRecProtBar[1]*pRecProtBar[1])));
			hPtPionFakeBar->Fill(sqrt((pRecPionBar[0]*pRecPionBar[0])+(pRecPionBar[1]*pRecPionBar[1])));

			hptLambdaBarFake->Fill(sqrt(pGenLambdaBar[0]*pGenLambdaBar[0]+pGenLambdaBar[1]*pGenLambdaBar[1]));
			for(int i=0;i<3;i++){
				hPprotFakeBar[i]->Fill(pRecProtBar[i]);
				hPpionFakeBar[i]->Fill(pRecPionBar[i]);
			}
			if((valBar[0]<0.1 && valBar[0]>-0.1) && (massinvrecBar>1.26 && massinvrecBar<1.38)){
				/*printf("#alpha %f , massa invariante %f \n", val[0],massinvrec);
				printf("protone: momRec %f %f %f ;  numfake=%f \n",pRecProt[0],pRecProt[1],pRecProt[2],nfaketrkprot);
				printf("pione: momRec %f %f %f ;  numfake=%f \n",pRecPion[0],pRecPion[1],pRecPion[2],nfaketrkpion);*/

				double p1abs=sqrt(pRecProtBar[0]*pRecProtBar[0]+pRecProtBar[1]*pRecProtBar[1]+pRecProtBar[2]*pRecProtBar[2]);
				thetaProtBar=TMath::ACos(pRecProtBar[2]/p1abs);
				hthetaProtBarFake->Fill(thetaProtBar);
				double p2abs=sqrt(pRecPionBar[0]*pRecPionBar[0]+pRecPionBar[1]*pRecPionBar[1]+pRecPionBar[2]*pRecPionBar[2]);
				thetaPionBar=TMath::ACos(pRecPionBar[2]/p2abs);
				hthetaPionBarFake->Fill(thetaPionBar);

			}
		}

		if(nfaketrkprotBar == 0) {
			hMassinrecnfakeprot0->Fill(massinvrecBar);
			hMassinrecnfakeprot0->SetTitle(Form("Invariant Mass for %f faketrk proton",nfaketrkprotBar));
			hDCAnfakeprot0->Fill(dcaBar);
			hDCAnfakeprot0->SetTitle(Form("DCA for %f faketrk proton",nfaketrkprotBar));
		}
		if(nfaketrkprotBar == 1){
			hMassinrecnfakeprot1->Fill(massinvrecBar);
			hMassinrecnfakeprot1->SetTitle(Form("Invariant Mass for %f faketrk proton",nfaketrkprotBar));
		}
		if(nfaketrkprotBar == 2){
			hMassinrecnfakeprot2->Fill(massinvrecBar);
			hMassinrecnfakeprot2->SetTitle(Form("Invariant Mass for %f faketrk proton",nfaketrkprotBar));
		}
		if(nfaketrkprotBar == 3) {
			hMassinrecnfakeprot3->Fill(massinvrecBar);
			hMassinrecnfakeprot3->SetTitle(Form("Invariant Mass for %f faketrk proton",nfaketrkprotBar));
		}			
		if(nfaketrkprotBar == 4) {
			hMassinrecnfakeprot4->Fill(massinvrecBar);
			hMassinrecnfakeprot4->SetTitle(Form("Invariant Mass for %f faketrk proton",nfaketrkprotBar));
			hDCAnfakeprot4->Fill(dcaBar);
			hDCAnfakeprot4->SetTitle(Form("DCA for %f faketrk proton",nfaketrkprotBar));
		}
		if(nfaketrkprotBar == 5) {
			hMassinrecnfakeprot5->Fill(massinvrecBar);
			hMassinrecnfakeprot5->SetTitle(Form("Invariant Mass for %f faketrk proton",nfaketrkprotBar));
		}
	}
	
	filinBar.Close();	
	hPzVsYprot->SetDirectory(0);
	hPzVsYpion->SetDirectory(0);
	hPzVsYpion->GetYaxis()->SetTitle("Pz_{rec}(#frac{GeV}{c})");
	hPzVsYpion->GetXaxis()->SetTitle("y");

	TCanvas *Ppion=new TCanvas("Ppion");
	hPzVsYpion->Draw("same");
	TLatex* lat1 = new TLatex();
  	lat1->SetTextFont(42);
  	lat1->SetTextSize(0.04);
  	lat1->SetNDC();
	  lat1->DrawLatex(0.5,0.95,"Pion");
	lat1->DrawLatex(0.5,0.85,"Pz_{rec} vs y");
	lat1->DrawLatex(0.5,0.75,Form("Mean x %0.3f ",hPzVsYpion->GetMean(1)));
	lat1->DrawLatex(0.5,0.65,Form("Mean y %0.2f ",hPzVsYpion->GetMean(2)));
	hPzVsYprot->GetYaxis()->SetTitle("Pz_{rec}(#frac{GeV}{c})");
	hPzVsYprot->GetXaxis()->SetTitle("y");
	TCanvas *Pprot=new TCanvas("Pprot");
	hPzVsYprot->Draw("same");
	TLatex* lat2 = new TLatex();
  lat2->SetTextFont(42);
  lat2->SetTextSize(0.04);
  lat2->SetNDC();
	  lat2->DrawLatex(0.5,0.95,"Proton");
	lat2->DrawLatex(0.5,0.85,"Pz_{rec} vs y");
lat2->DrawLatex(0.5,0.75,Form("Mean x %0.3f ",hPzVsYprot->GetMean(1)));
lat2->DrawLatex(0.5,0.65,Form("Mean y %0.2f ",hPzVsYprot->GetMean(2)));


	hPzpionVsPzLambda->GetYaxis()->SetTitle("Pz_{gen} #pi (#frac{GeV}{c})");
	hPzpionVsPzLambda->GetXaxis()->SetTitle("Pz_{gen} #Lambda(#frac{GeV}{c})");
	hPzprotVsPzLambda->GetYaxis()->SetTitle("Pz_{gen} p (#frac{GeV}{c})");
	hPzprotVsPzLambda->GetXaxis()->SetTitle("Pz_{gen} #Lambda(#frac{GeV}{c})");

	TH1F *hDivProt = (TH1F*)hthetaProtFake->Clone("hDivProt"); 
	hDivProt->Divide(hthetaProtFake,hthetaProt,1,1,"B");
	TH1F *hDivPion = (TH1F*)hthetaPionFake->Clone("hDivPion"); 
	hDivPion->Divide(hthetaPionFake,hthetaPion,1,1,"B");
	new TCanvas();
	hArmen->Draw("same");
	TLatex* lat = new TLatex();
  	lat->SetTextFont(42);
  	lat->SetTextSize(0.04);
  	lat->SetNDC();
	lat->DrawLatex(0.5,0.85,"#Lambda");
lat->DrawLatex(0.5,0.75,"#bar{#Lambda}");

	foutLambda->cd();
	hArmen->Write();
	hNfakeVszP->Write();
	hMassinrecnfakeprot0->Write();
	hMassinrecnfakeprot1->Write();
	hMassinrecnfakeprot2->Write();
	hMassinrecnfakeprot3->Write();
	hMassinrecnfakeprot4->Write();
	hMassinrecnfakeprot5->Write();
	hDCAnfakeprot0->Write();
	hDCAnfakeprot4->Write();
	hMassinrecY->Write();
	hMassinrecSwitchY->Write();
	hMassinrecY0fake->Write();
	hMassinrecSwitchY0fake->Write();
	hDivProt->Write();
	hDivPion->Write();
	hResVxVsY->Write();
  	hResVyVsY->Write();
  	hResVzVsY->Write();
	hResVxVsLz->Write();
  	hResVyVsLz->Write();
  	hResVzVsLz->Write();

	hPzVsYprot->Write();
	hPzVsYpion->Write();

	hPzprotVsPzLambda->Write();
	hPzpionVsPzLambda->Write();

	hMassinreczP->Write();

	for(int i=0;i<3;i++){
		hResPVsYpionRel[i]->Write();
		hResPVsYprotRel[i]->Write();
		hResPVsYprot[i]->Write();
		hResPVsYpion[i]->Write();
		hPprot[i]->Write();
		hPpion[i]->Write();
		hResPVsLzprot[i]->Write();
		hResPVsLzpion[i]->Write();
	}
	for(int i=0;i<3;i++){
		hPprotFake[i]->Write();
		hPpionFake[i]->Write();
		hPprotFakeBar[i]->Write();
		hPpionFakeBar[i]->Write();
	}
	for(int i=0;i<2;i++){
		hChi2Prot[i]->Write();
		hChi2ITSProt[i]->Write();
		hChi2Pion[i]->Write();
		hChi2ITSPion[i]->Write();
	}

	hMArm->Write();
	hMArmBar->Write();
	hPtprotfakeVsyLambda->Write();
	hPtpionfakeVsyLambda->Write();

	hPtprotfakeVsyLambdaBar->Write();
	hPtpionfakeVsyLambdaBar->Write();

	
	hthetaProt->Write();
	hthetaPion->Write();
	hthetaProtFake->Write();
	hthetaPionFake->Write();

	hthetaProtBar->Write();
	hthetaPionBar->Write();
	hthetaProtBarFake->Write();
	hthetaPionBarFake->Write();

	

	
	hNfakeProtVsChi2->Write();
	hNfakePionVsChi2->Write();

	hLzrecVsLzgen->Write();
	hDeltaLzVsLzgen->Write();

	hArmenLambda->Write();
	hArmenLambdaBar->Write();
	hArmenLambdaFake->Write();
	hArmenLambdaBarFake->Write();

	hyLambda->Write();
	hyLambdaFake->Write();
	hptLambda->Write();
	hptLambdaFake->Write();
	hPtProtFake->Write();
	hPtPionFake->Write();


	hyLambdaBar->Write();
	hyLambdaBarFake->Write();
	hptLambdaBar->Write();
	hptLambdaBarFake->Write();
	hPtProtFakeBar->Write();
	hPtPionFakeBar->Write();



	TH1D* hInvMass0fake=(TH1D *)hMassinrecY0fake->ProjectionX();
	hInvMass0fake->SetDirectory(0);
	hInvMass0fake->SetLineColor(2);

	TCanvas *inv=new TCanvas("inv");
	hInvMass0fake->SetTitle("Invariant mass with 0 fake");
	hInvMass0fake->GetXaxis()->SetTitle("m_{p#pi} (#frac{GeV}{c^{2}})");
	hInvMass0fake->GetYaxis()->SetTitle("Counts");
	hInvMass0fake->GetXaxis()->SetRangeUser(1,1.2);
	hInvMass0fake->Draw();
	
	TH1D* hInvMassSwitch0fake=(TH1D *)hMassinrecSwitchY0fake->ProjectionX();
	hInvMassSwitch0fake->SetDirectory(0);
	hInvMassSwitch0fake->SetLineColor(4);
	TCanvas *invSwap=new TCanvas("InvSwap");
	hInvMassSwitch0fake->SetTitle("Invariant mass swap with 0 fake");
	hInvMassSwitch0fake->GetXaxis()->SetTitle("swap m_{p#pi} (#frac{GeV}{c^{2}})");
	hInvMassSwitch0fake->GetYaxis()->SetTitle("Counts");
	hInvMassSwitch0fake->Draw();



	TH1D* MassInvProjswap[27];
	TH1D*  MassInvProj[27];
	TCanvas* cSwap=new TCanvas("cSwap");
	int count=0;
	TCanvas* c=new TCanvas("c");
	c->Divide(6,5);
  	cSwap->Divide(6,5);
  	for(Int_t jj=1, ii=1; jj<hMassinrecSwitchY->GetNbinsY() && ii<hMassinrecY->GetNbinsY(); jj++, ii++){
 		//printf("n=%d ",hist2D->GetNbinsY());
    	TH1D *hInvMassSwitch = (TH1D *)hMassinrecSwitchY->ProjectionX(Form("%sBin%d",hMassinrecSwitchY->GetName(),jj), jj, jj);
		TH1D *hInvMass = (TH1D *)hMassinrecY->ProjectionX(Form("%sBin%d",hMassinrecY->GetName(),ii), ii, ii);
		//hInvMass->SetDirectory(0);
    	//hInvMass->Draw();
		//i=i+2;
    	if(hInvMassSwitch->GetSumOfWeights()>45 && hMassinrecY->GetSumOfWeights()>45){
			MassInvProjswap[count]=hInvMassSwitch;
			MassInvProj[count]=hInvMass;
			count++;
			cSwap->cd(count);
			cSwap->SetTitle(Form("%s",hMassinrecSwitchY->GetName()));
			hInvMassSwitch->SetDirectory(0);
			hInvMassSwitch->SetLineColor(4);
			hInvMassSwitch->Draw();
			c->cd(count);
			hInvMass->SetDirectory(0);
			hInvMass->SetLineColor(2);
    		hInvMassSwitch->Draw("hist");
			hInvMass->Draw("SAME");
		}
	}
	
	TH1F *hsumSwap = (TH1F*)MassInvProjswap[0]->Clone("hsumSwap");
	TH1F *hsum = (TH1F*)MassInvProj[0]->Clone("hsum");
	for(int i=1; i<count;i++){
			hsumSwap->Add(MassInvProjswap[i]);
			hsum->Add(MassInvProj[i]);
		}
	hsumSwap->SetDirectory(0);
	hsum->SetDirectory(0);
	TCanvas *ms= new TCanvas("ms");
	hsumSwap->SetTitle("Invariant mass swap with fake");
	hsumSwap->GetXaxis()->SetTitle("swap m_{p#pi} (#frac{GeV}{c^{2}})");
	hsumSwap->GetYaxis()->SetTitle("Counts");
	hsumSwap->Draw("hist");
	TCanvas *m = new TCanvas("m");
	hsum->SetTitle("Invariant mass with fake");
	hsum->GetXaxis()->SetTitle("m_{p#pi} (#frac{GeV}{c^{2}})");
	hsum->GetYaxis()->SetTitle("Counts");
	hsum->Draw();
	hsum->SetLineColor(2);
	hsumSwap->SetLineColor(4);
	TCanvas *spp= new TCanvas("spp");
	hsumSwap->Draw("hist");

	hsum->Draw("same");


	//Projection for checking Pz resoluztion of pion

	/*int totvalue=7;
	double meanPz[7]={0.}, RMSPz[7]={0.}, meanErrPz[7]={0.}, RMSErrPz[7]={0.}, yRange[7]={0.};
	TH1D* hProjReszY[7];
	double jrange=0.;
	for(int i=0;i<totvalue;i++){
		hProjReszY[i]=hResol->ProjectionX((Form("Res Pz with y range [%f %f]",jrange,jrange+1)), hResol->GetYaxis()->FindBin(jrange+0.01),hResol->GetYaxis()->FindBin(jrange+0.99));
		meanPz[i]=hProjResyLz[i]->GetMean();
		RMSPz[i]=hProjResyLz[i]->GetRMS();
		meanErrPz[i]=hProjResyLz[i]->GetMeanError();
		RMSErrPz[i]=hProjResyLz[i]->GetRMSError();
		yRange[i]=(jrange+jrange+1)/2.;
		hProjReszY[i]->SetDirectory(0);
		new TCanvas();
		hProjReszY[i]->Draw();
		jrange=jrange+1;
		}*
		
	ProjectionResolCorrPy(hPzVsYprot,1,"Pz Prot");
	ProjectionResolCorrPy(hPzVsYpion,1,"Pz Pion");




/*	TCanvas* c2p=new TCanvas("c2p");
	ProjectionforBin(hResPVsYprot[0] ,"#Lambda+#bar{#Lambda} ResPx vs Y prot", "Y", c2p, foutLambda);
	TCanvas* c3p=new TCanvas("c3p");
	ProjectionforBin(hResPVsYprot[1], "#Lambda+#bar{#Lambda} ResPy vs Y prot", "Y", c3p, foutLambda);
	TCanvas* c4p=new TCanvas("c4p");
	ProjectionforBin(hResPVsYprot[2], "#Lambda+#bar{#Lambda} ResPz vs Y prot", "Y", c4p, foutLambda);
	TCanvas* c5p=new TCanvas("c5p");
	ProjectionforBin(hResPVsYpion[0], "#Lambda+#bar{#Lambda} ResPx vs Y pion", "Y", c5p, foutLambda);
	TCanvas* c6p=new TCanvas("c6p");
	ProjectionforBin(hResPVsYpion[1], "#Lambda+#bar{#Lambda} ResPy vs Y pion", "Y", c6p, foutLambda);
	TCanvas* c7p=new TCanvas("c7p");
	ProjectionforBin(hResPVsYpion[2], "#Lambda+#bar{#Lambda} ResPz vs Y pion", "Y", c7p, foutLambda);*/

	TCanvas* c5piRel=new TCanvas("c5piRel");
	ProjectionforBin(hResPVsYpionRel[0], "#Lambda+#bar{#Lambda} pion ResPx vs Y Rel", "Y", c5piRel, foutLambda);
	TCanvas* c6piRel=new TCanvas("c6piRel");
	ProjectionforBin(hResPVsYpionRel[1], "#Lambda+#bar{#Lambda} pion ResPy vs Y Rel", "Y", c6piRel, foutLambda);
	TCanvas* c7piRel=new TCanvas("c7piRel");
	ProjectionforBin(hResPVsYpionRel[2], "#Lambda+#bar{#Lambda} pion ResPz vs Y Rel", "Y", c7piRel, foutLambda);

	TCanvas* c5pRel=new TCanvas("c5pRel");
	ProjectionforBin(hResPVsYprotRel[0], "#Lambda+#bar{#Lambda} prot ResPx vs Y Rel", "Y", c5pRel, foutLambda);
	TCanvas* c6pRel=new TCanvas("c6pRel");
	ProjectionforBin(hResPVsYprotRel[1], "#Lambda+#bar{#Lambda} prot ResPy vs Y Rel", "Y", c6pRel, foutLambda);
	TCanvas* c7pRel=new TCanvas("c7pRel");
	ProjectionforBin(hResPVsYprotRel[2], "#Lambda+#bar{#Lambda} prot ResPz vs Y Rel", "Y", c7pRel, foutLambda);

		TCanvas* c1=new TCanvas("c1");
	ProjectionforBinMassInv(hMassinrecY0fake, "#Lambda+#bar{#Lambda} Inv Mass vs Y", "Y", c1, foutLambda);

	TCanvas* resP=new TCanvas("resP");
	
	resP->Divide(3,1);
	TMultiGraph *mgp[3], *mgpLD[3];
	mgp[0] = new TMultiGraph();
	mgp[1] = new TMultiGraph();
	mgp[2] = new TMultiGraph();
	mgpLD[0] = new TMultiGraph();
	mgpLD[1] = new TMultiGraph();
	mgpLD[2] = new TMultiGraph();
	mgpLD[0]->SetTitle("Relative Resolutions on P_{x} for #Lambda and D_{0} daughters; y ; #sigma");
	mgpLD[1]->SetTitle("Relative Resolutions on P_{y} for #Lambda and D_{0} daughters; y ; #sigma");
	mgpLD[2]->SetTitle("Relative Resolutions on P_{z} for #Lambda and D_{0} daughters; y ; #sigma");
	/*mgp[0]->SetTitle("Px relative resolution; y ; #sigma");
	mgp[1]->SetTitle("Py relative resolution; y ; #sigma");
	mgp[2]->SetTitle("Pz relative resolution; y ; #sigma");*/
	TGraphErrors *pprot[3], *ppion[3], *pkaonD0[3], *ppionD0[3],*pprotML[3], *ppionML[3];
	pprot[0]=(TGraphErrors*)foutLambda->Get("Sigma#Lambda+#bar{#Lambda} prot ResPx vs Y Rel");
	pprot[1]=(TGraphErrors*)foutLambda->Get("Sigma#Lambda+#bar{#Lambda} prot ResPy vs Y Rel");
	pprot[2]=(TGraphErrors*)foutLambda->Get("Sigma#Lambda+#bar{#Lambda} prot ResPz vs Y Rel");
	ppion[0]=(TGraphErrors*)foutLambda->Get("Sigma#Lambda+#bar{#Lambda} pion ResPx vs Y Rel");
	ppion[1]=(TGraphErrors*)foutLambda->Get("Sigma#Lambda+#bar{#Lambda} pion ResPy vs Y Rel");
	ppion[2]=(TGraphErrors*)foutLambda->Get("Sigma#Lambda+#bar{#Lambda} pion ResPz vs Y Rel");
/*	pprotML[0]=(TGraphErrors*)fLMoreLayer->Get("Sigma#Lambda+#bar{#Lambda} prot ResPx vs Y Rel");
	pprotML[1]=(TGraphErrors*)fLMoreLayer->Get("Sigma#Lambda+#bar{#Lambda} prot ResPy vs Y Rel");
	pprotML[2]=(TGraphErrors*)fLMoreLayer->Get("Sigma#Lambda+#bar{#Lambda} prot ResPz vs Y Rel");
	ppionML[0]=(TGraphErrors*)fLMoreLayer->Get("Sigma#Lambda+#bar{#Lambda} pion ResPx vs Y Rel");
	ppionML[1]=(TGraphErrors*)fLMoreLayer->Get("Sigma#Lambda+#bar{#Lambda} pion ResPy vs Y Rel");
	ppionML[2]=(TGraphErrors*)fLMoreLayer->Get("Sigma#Lambda+#bar{#Lambda} pion ResPz vs Y Rel");
	pkaonD0[0]=(TGraphErrors*)hD0->Get("SigmahResPxVsYKaonRel");
	pkaonD0[1]=(TGraphErrors*)hD0->Get("SigmahResPyVsYKaonRel");
	pkaonD0[2]=(TGraphErrors*)hD0->Get("SigmahResPzVsYKaonRel");
	ppionD0[0]=(TGraphErrors*)hD0->Get("SigmahResPxVsYpionRel");
	ppionD0[1]=(TGraphErrors*)hD0->Get("SigmahResPyVsYpionRel");
	ppionD0[2]=(TGraphErrors*)hD0->Get("SigmahResPzVsYpionRel");*/
	for(int i=0;i<3;i++){
		resP->cd(i+1);
		gPad->SetGrid();
		gStyle->SetLineScalePS(1);
		//5LAYER
		pprot[i]->SetMarkerColor(kBlue);
		ppion[i]->SetMarkerColor(kRed-4);
		pprot[i]->SetLineColor(kBlue);
		ppion[i]->SetLineColor(kRed-4);
		pprot[i]->SetMarkerStyle(27);
		ppion[i]->SetMarkerStyle(25);
		pprot[i]->SetMarkerSize(1.2);
		ppion[i]->SetMarkerSize(1);
		//10LAYER
	/*	pprotML[i]->SetMarkerColor(kAzure+10);
		ppionML[i]->SetMarkerColor(93);
		pprotML[i]->SetLineColor(kAzure+10);
		ppionML[i]->SetLineColor(93);
		pprotML[i]->SetMarkerStyle(21);
		ppionML[i]->SetMarkerStyle(21);
		pprotML[i]->SetMarkerSize(0.8);
		ppionML[i]->SetMarkerSize(0.8);*/
		//5LAY
		mgp[i]->Add(pprot[i],"p");
		mgp[i]->Add(ppion[i],"p");
		//10LAY
		//mgp[i]->Add(pprotML[i],"p");
		//mgp[i]->Add(ppionML[i],"p");
		mgp[i]->GetXaxis()->SetTitle("y");
		mgp[i]->GetYaxis()->SetTitle("#sigma");
		mgp[i]->Draw("a");
	}
	resP->cd(1);
	auto legend = new TLegend(0.1,0.7,0.48,0.9);
   	//legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
   	legend->AddEntry(pprot[0],"proton","lep");
	//legend->AddEntry(pprotML[0],"#Lambda proton with 10 layers 50#mum","lep");
  	legend->AddEntry(ppion[0],"pion","lep");
  	//legend->AddEntry(ppionML[0],"#Lambda pion with 10 layers 50#mum","lep");

   	legend->Draw();


	//mg->SetTitle("Different resolutions on P_{x}; y ; #sigma (#frac{Gev}/{c})");
	TCanvas* resPLD=new TCanvas("resPLD");
	resPLD->Divide(3,1);
	//resPLD->SetTitle("pippo");
	//TLatex * tex = new TLatex(0.5,0.8,"My Title"); 
	//tex->SetNDC(); 
	//tex->SetTextSize(0.2); 
	//tex->Draw(); 

/*	for(int i=0;i<3;i++){
		resPLD->cd(i+1);
		gPad->SetGrid();
		
		pprot[i]->SetMarkerColor(kBlue+1);
		ppion[i]->SetMarkerColor(kRed+1);
		pkaonD0[i]->SetMarkerColor(kAzure+7);
		ppionD0[i]->SetMarkerColor(kOrange+7);
		pprot[i]->SetLineColor(kBlue+1);
		ppion[i]->SetLineColor(kRed+1);
		pkaonD0[i]->SetLineColor(kAzure+7);
		ppionD0[i]->SetLineColor(kOrange+7);
		
		pprot[i]->SetMarkerStyle(20);
		ppion[i]->SetMarkerStyle(20);
		pkaonD0[i]->SetMarkerStyle(21);
		ppionD0[i]->SetMarkerStyle(21);
		pprot[i]->SetMarkerSize(0.9);
		ppion[i]->SetMarkerSize(0.9);
		pkaonD0[i]->SetMarkerSize(0.9);
		ppionD0[i]->SetMarkerSize(0.9);
		mgpLD[i]->Add(pprot[i],"p");
		mgpLD[i]->Add(ppion[i],"p");
		mgpLD[i]->Add(pkaonD0[i],"p");
		mgpLD[i]->Add(ppionD0[i],"p");
		mgpLD[i]->GetXaxis()->SetTitle("y");
		mgpLD[i]->GetYaxis()->SetTitle("#sigma");
		mgpLD[i]->Draw("a");
	}
	resPLD->cd(1);
	auto legendD = new TLegend(0.1,0.7,0.48,0.9);
   	//legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
   	legendD->AddEntry(pprot[0],"#Lambda proton","lep");
  	legendD->AddEntry(ppion[0],"#Lambda pion","lep");
   	legendD->AddEntry(pkaonD0[0],"D_{0} kaon","lep");
   	legendD->AddEntry(ppionD0[0],"D_{0} pion","lep");	   
   	legendD->Draw();
	
*/
		
	
	TCanvas* c01=new TCanvas("c01");
	ProjectionforBinMassInv(hMassinreczP, "#Lambda+#bar{#Lambda} Inv Mass vs Lz", "zP", c01, foutLambda);
	
	//Vertex resolution
	TCanvas* c2=new TCanvas("c2");
	ProjectionforBin(hResVxVsY, "#Lambda+#bar{#Lambda} ResX vs Y", "Y", c2, foutLambda);
	TCanvas* c3=new TCanvas("c3");
	ProjectionforBin(hResVyVsY, "#Lambda+#bar{#Lambda} ResY vs Y", "Y", c3, foutLambda);
	TCanvas* c4=new TCanvas("c4");
	ProjectionforBin(hResVzVsY, "#Lambda+#bar{#Lambda} ResZ vs Y", "Y", c4, foutLambda);
	TCanvas* c5=new TCanvas("c5");
	ProjectionforBin(hResVxVsLz, "#Lambda+#bar{#Lambda} ResX vs Lz", "Lz", c5, foutLambda);
	TCanvas* c6=new TCanvas("c6");
	ProjectionforBin(hResVyVsLz, "#Lambda+#bar{#Lambda} ResY vs Lz", "Lz", c6, foutLambda);
	TCanvas* c7=new TCanvas("c7");
	ProjectionforBin(hResVzVsLz, "#Lambda+#bar{#Lambda} ResZ vs Lz", "Lz", c7, foutLambda);

	TGraphErrors *invMassy, *invMassLz, *invMassyML, *invMassLzML;
	invMassy=(TGraphErrors*)foutLambda->Get("SigmahMassinvVsY0fake");
	invMassLz=(TGraphErrors*)foutLambda->Get("SigmahMassinrecVszP");
	//invMassyML=(TGraphErrors*)fLMoreLayer->Get("SigmahMassinvVsY0fake");
	//invMassLzML=(TGraphErrors*)fLMoreLayer->Get("SigmahMassinrecVszP");
	TCanvas *cInvy=new TCanvas("cInvy");
	invMassy->SetMarkerStyle(20);
	invMassy->SetMarkerColor(kBlue-4);
	invMassy->SetMarkerSize(1);
	invMassy->GetXaxis()->SetTitle("y");
	invMassy->GetYaxis()->SetTitle("#sigma_{m} (GeV/c^{2})");
	invMassy->SetTitle("#sigma Inv Mass Vs y");
	//invMassyML->SetMarkerStyle(21);
	//invMassyML->SetMarkerColor(kOrange+7);
	//invMassyML->SetMarkerSize(0.9);
	invMassy->Draw("AP");
	//invMassyML->Draw("P");


	TCanvas *cInvLz=new TCanvas("cInvLz");
	invMassLz->SetMarkerStyle(20);
	//invMassLz->SetFillStyle(3001);
	invMassLz->SetFillStyle(3001);
	invMassLz->SetLineWidth(2);
	invMassLz->SetMarkerColor(kBlue-4);
	invMassLz->SetMarkerSize(1);
	invMassLz->GetXaxis()->SetTitle("Lz_{rec} (cm)");
	invMassLz->GetYaxis()->SetTitle("#sigma_{m} (GeV/c^{2})");
	invMassLz->SetTitle("#sigma Inv Mass Vs Lz_{rec}");
	//invMassLzML->SetMarkerStyle(21);
	//invMassLzML->SetMarkerColor(kOrange+7);
	//invMassLzML->SetMarkerSize(0.9);
	invMassLz->Draw("A2P");
	//invMassLzML->Draw("P");



	TGraphErrors *VertVsy[3], *VertVsLz[3],*VertVsyML[3], *VertVsLzML[3];
	VertVsy[0]=(TGraphErrors*)foutLambda->Get("Sigma#Lambda+#bar{#Lambda} ResX vs Y");
	VertVsy[1]=(TGraphErrors*)foutLambda->Get("Sigma#Lambda+#bar{#Lambda} ResY vs Y");
	VertVsy[2]=(TGraphErrors*)foutLambda->Get("Sigma#Lambda+#bar{#Lambda} ResZ vs Y");
	VertVsLz[0]=(TGraphErrors*)foutLambda->Get("Sigma#Lambda+#bar{#Lambda} ResX vs Lz");
	VertVsLz[1]=(TGraphErrors*)foutLambda->Get("Sigma#Lambda+#bar{#Lambda} ResY vs Lz");
	VertVsLz[2]=(TGraphErrors*)foutLambda->Get("Sigma#Lambda+#bar{#Lambda} ResZ vs Lz");
	/*VertVsyML[0]=(TGraphErrors*)fLMoreLayer->Get("Sigma#Lambda+#bar{#Lambda} ResX vs Y");
	VertVsyML[1]=(TGraphErrors*)fLMoreLayer->Get("Sigma#Lambda+#bar{#Lambda} ResY vs Y");
	VertVsyML[2]=(TGraphErrors*)fLMoreLayer->Get("Sigma#Lambda+#bar{#Lambda} ResZ vs Y");
	VertVsLzML[0]=(TGraphErrors*)fLMoreLayer->Get("Sigma#Lambda+#bar{#Lambda} ResX vs Lz");
	VertVsLzML[1]=(TGraphErrors*)fLMoreLayer->Get("Sigma#Lambda+#bar{#Lambda} ResY vs Lz");
	VertVsLzML[2]=(TGraphErrors*)fLMoreLayer->Get("Sigma#Lambda+#bar{#Lambda} ResZ vs Lz");*/

	TMultiGraph *mgV[3], *mgLz[3];
	mgV[0] = new TMultiGraph();
	mgV[1] = new TMultiGraph();
	mgV[2] = new TMultiGraph();
	mgLz[0] = new TMultiGraph();
	mgLz[1] = new TMultiGraph();
	mgLz[2] = new TMultiGraph();
	/*mgLz[0]->SetTitle("Resolution on V_{x} vs Lz");
	mgLz[1]->SetTitle("Resolution on V_{y} vs Lz");
	mgLz[2]->SetTitle("Resolution on V_{z} vs Lz");
	mgV[0]->SetTitle("Resolution on V_{x} vs y");
	mgV[1]->SetTitle("Resolution on V_{y} vs y");
	mgV[2]->SetTitle("Resolution on V_{z} vs y");*/
	TCanvas *cVy = new TCanvas("cVy");
	cVy->Divide(3,1);
	TCanvas *cVLz = new TCanvas("cVLz");
	cVLz->Divide(3,1);
	/*VertVsy[0]->SetTitle("Resolution on V_{x} vs y");
	VertVsy[1]->SetTitle("Resolution on V_{y} vs y");
	VertVsy[2]->SetTitle("Resolution on V_{z} vs y");

	VertVsLz[0]->SetTitle("Resolution on V_{x} vs Lz");
	VertVsLz[1]->SetTitle("Resolution on V_{y} vs Lz");
	VertVsLz[2]->SetTitle("Resolution on V_{z} vs Lz");*/
	for(int i=0;i<3;i++){
		cVy->cd(i+1);
		gPad->SetGrid();
		VertVsy[i]->SetMarkerStyle(20);
		VertVsy[i]->SetMarkerColor(4);
		VertVsy[i]->SetLineColor(4);
		VertVsy[i]->SetMarkerSize(1.1);
		//VertVsyML[i]->SetMarkerStyle(21);
		//VertVsyML[i]->SetMarkerColor(kOrange+7);
		//VertVsyML[i]->SetMarkerSize(0.9);
		VertVsy[i]->GetXaxis()->SetTitle("y");
		VertVsy[i]->GetXaxis()->SetLimits(1,4.8);
		VertVsy[i]->GetYaxis()->SetTitle("#sigma (#mum) ");
		mgV[i]->Add(VertVsy[i],"p");
		//mgV[i]->Add(VertVsyML[i],"p");
		mgV[i]->GetXaxis()->SetTitle("y");
		mgV[i]->GetXaxis()->SetLimits(1,4.8);
		mgV[i]->GetYaxis()->SetTitle("#sigma (#mum) ");
		mgV[i]->Draw("a");
		cVLz->cd(i+1);
		gPad->SetGrid();
		VertVsLz[i]->SetMarkerStyle(21);
		VertVsLz[i]->SetMarkerColor(kAzure+2);
		VertVsLz[i]->SetLineColor(kAzure+2);
		VertVsLz[i]->SetMarkerSize(1.1);

		//VertVsLzML[i]->SetMarkerStyle(21);
		//VertVsLzML[i]->SetMarkerColor(kOrange+7);
		//VertVsLzML[i]->SetMarkerSize(0.9);
		mgLz[i]->Add(VertVsLz[i],"p");
		//mgLz[i]->Add(VertVsLzML[i],"p");
		mgLz[i]->GetXaxis()->SetTitle("Lz (cm)");
		
		mgLz[i]->GetYaxis()->SetTitle("#sigma (#mum) ");
		mgLz[i]->Draw("a");
	} 
/*	 cVy->cd(1);
	auto legendy = new TLegend(0.1,0.7,0.48,0.9);
   	//legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
   	legendy->AddEntry(VertVsy[0],"Resolution 5 layers","lep");
  	//legendy->AddEntry(VertVsyML[0],"Resolution 10 layers","lep");
	  legendy->Draw();
	  cVLz->cd(1);
	  auto legendLz = new TLegend(0.1,0.7,0.48,0.9);
   	legendLz->AddEntry(VertVsLz[0],"Resolution 5 layers","lep");
   	legendLz->AddEntry(VertVsLzML[0],"Resolution 10 layers","lep");	   
   	legendLz->Draw();
	*/   
	

	
	hLambda->Close();
	hLambdaBar->Close();
/*	hPion->Close();
	hProt->Close();
	hAntiProt->Close();
	hPionpos->Close();
	hD0->Close();*/
	foutLambda->Close();

}

void Armenteros(float mp[3], float mn[3], float mm[3], double v[2]){
	TVector3 vecN(mn[0],mn[1],mn[2]);
  	TVector3 vecP(mp[0],mp[1],mp[2]);
  	TVector3 vecM(mm[0],mm[1],mm[2]);
  
	Double_t costhetaP = (vecP * vecM)/(vecP.Mag() * vecM.Mag());
	Double_t costhetaN = (vecN * vecM)/(vecN.Mag() * vecM.Mag());
	Double_t sin2thetaP = 1 - costhetaP*costhetaP;
	Double_t sinthetaP = 0;
	if(sin2thetaP>0) sinthetaP=TMath::Sqrt(sin2thetaP);

	Double_t alfa = ((vecP.Mag())*costhetaP-(vecN.Mag())*costhetaN)/
		((vecP.Mag())*costhetaP+(vecN.Mag())*costhetaN) ;
	Double_t qt = vecP.Mag()*sinthetaP;

	v[0]=alfa;
	v[1]=qt;
}


void Projection(TH2F* h2D,  int range, bool rapidity,const char *name ){
	double mean[40], RMS[40], meanErr[40], RMSErr[40], Lz[40];
	std::vector<TH1D*>hProject;
	double jrange=0.;
	if(rapidity==1) jrange=1.5;
	for(int i=0;i<range;i++){
		hProject.push_back(h2D->ProjectionY(Form("%s range [%f %f]",name,jrange,jrange+0.5), h2D->GetXaxis()->FindBin(jrange+0.01),h2D->GetXaxis()->FindBin(jrange+0.49)));
		mean[i]=(hProject.at(i)->GetMean());
		RMS[i]=(hProject.at(i)->GetRMS());
		meanErr[i]=(hProject.at(i)->GetMeanError());
		RMSErr[i]=(hProject.at(i)->GetRMSError());
		Lz[i]=((jrange+jrange+0.5)/2.);
		jrange=jrange+0.5;
		new TCanvas();
		hProject.at(i)->SetDirectory(0);
		hProject.at(i)->SetMinimum(0);
		hProject.at(i)->Draw();
	}
	/*if(rapidity==1){
		printf("range  %d",range);
		TGraphErrors* gr = new TGraphErrors(6,Lz,mean,0.,meanErr);
		gr->SetTitle("Mean  vs ?? ");
		new TCanvas();
		gr->Draw("AC*");
	}*/
}



void ProjectionResolCorrPy(TH2F* hResol, bool rapid,const char *title){
	int totvalue=4;
	double meanPy[10]={0.}, RMSPy[10]={0.}, meanErrPy[10]={0.}, RMSErrPy[10]={0.}, Lzrange[10]={0.};
	TH1D* hProjResyLz[10];
	double jrange=1.;

	for(int i=0;i<totvalue;i++){
		hProjResyLz[i]=hResol->ProjectionY((Form("Pz with y range [%f %f]",jrange,jrange+1)), hResol->GetXaxis()->FindBin(jrange+0.01),hResol->GetXaxis()->FindBin(jrange+0.99));
		meanPy[i]=hProjResyLz[i]->GetMean();
		RMSPy[i]=hProjResyLz[i]->GetRMS();
		meanErrPy[i]=hProjResyLz[i]->GetMeanError();
		RMSErrPy[i]=hProjResyLz[i]->GetRMSError();
		Lzrange[i]=(jrange+jrange+1)/2.;
		hProjResyLz[i]->SetDirectory(0);
		hProjResyLz[i]->SetTitle(Form(" %s in y range [%f %f] ",title,jrange,jrange+1));
		//new TCanvas();
		//hProjResyLz[i]->Draw();
		jrange=jrange+1;
		}
	if(rapid==0){
		TGraphErrors* grResy = new TGraphErrors(totvalue,Lzrange,meanPy,0,meanErrPy);
		grResy->SetTitle(Form("Mean %s  ", title));
		//new TCanvas();
		//grResy->Draw("AC*");
		TGraphErrors* grRMSy = new TGraphErrors(totvalue,Lzrange,RMSPy,0,RMSErrPy);
		grRMSy->SetTitle(Form("RMS %s  ",title));
		//new TCanvas();
		//grRMSy->Draw("AC*");
			}
	else{
		TGraphErrors* grResyRapid = new TGraphErrors(totvalue,Lzrange,meanPy,0,meanErrPy);
		grResyRapid->SetTitle(Form("Mean %s ",title));
		//new TCanvas();
		//grResyRapid->Draw("AC*");
		TGraphErrors* grRMSyRapid = new TGraphErrors(totvalue,Lzrange,RMSPy,0,RMSErrPy);
		grRMSyRapid->SetTitle(Form("RMS %s ",title));
		grRMSyRapid->GetXaxis()->SetTitle("y");
		grRMSyRapid->GetYaxis()->SetTitle("RMS (#frac{GeV}{c})");
		//new TCanvas();
		//grRMSyRapid->Draw("AC*");
		}
	}

void ProjectionNfake(TH2F* hNfake,const char *title){
	int totlay=5;
	TH1D* hProjNfakeLz[totlay];
	TH1D* hNfakefract[totlay];
	int layer[6]={0,7,15,20,27,38};
	double counter=0.,conteggibin=0.;
	for(int jlay=0; jlay<totlay; jlay++){
		hProjNfakeLz[jlay]=hNfake->ProjectionX((Form("Nfake vs layer in %d cm ",layer[jlay+1])), hNfake->GetYaxis()->FindBin(layer[jlay]+0.01),hNfake->GetYaxis()->FindBin(layer[jlay+1]-0.01));
		hProjNfakeLz[jlay]->SetDirectory(0);
		//new TCanvas();
		//hProjNfakeLz[jlay]->SetTitle(Form("Nfake vs layer in %d cm \n ",layer[jlay+1]));
		//hProjNfakeLz[jlay]->Draw();
		printf("Nfake vs layer in %d cm \n ",layer[jlay+1]);
		//printf("bincenter %f,bin %d \n",hNfake->GetYaxis()->GetBinCenter((hNfake->GetYaxis()->FindBin(layer[jlay]+0.01))),hNfake->GetYaxis()->FindBin(layer[jlay]+0.01));
		//printf("bincenter %f,bin %d \n",hNfake->GetYaxis()->GetBinCenter((hNfake->GetYaxis()->FindBin(layer[jlay]-0.01))),hNfake->GetYaxis()->FindBin(layer[jlay]-0.01));
		for(int i=1;i<(hProjNfakeLz[jlay]->GetXaxis()->FindBin(6.01));i++){
			counter=counter+(hProjNfakeLz[jlay]->GetBinContent(i));
			conteggibin=hProjNfakeLz[jlay]->GetBinContent(i);
			//printf("frazione di conteggi nel bin %d: %f, conteggi %f  (bin corrispondente %f) \n",i,(conteggibin/hProjNfakeLz[jlay]->GetEntries()),conteggibin,hProjNfakeLz[jlay]->GetXaxis()->GetBinCenter(i));
		}			
	}
}

void ProjectionforBin(TH2F* hist2D,const char* graphTitle ,const char* Xaxisname, TCanvas *c0, TFile *f){
	double xaxis[200], mean[200], errmean[200], sigma[200], errsigma[200], sigmagaus[200], errsigmagaus[200];
	int count=0, i=0, countGraph=0;
  	c0->Divide(10,5);
  	for(Int_t jj=1; jj<hist2D->GetNbinsY(); jj++){
		TH1D *hInvMass = (TH1D *)hist2D->ProjectionX(Form("%sBin%d",hist2D->GetName(),jj), jj, jj);
		if(hInvMass->GetSumOfWeights()>45){
			count++;
			c0->cd(count);
			c0->SetTitle(Form("%s",graphTitle));
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
			ErrMinSigma=(sigmag*30)/100;
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
		//*gr = *gr2;
		//new TCanvas(Form("Graph vs  %s ",graphTitle),"Mean vs L_{z}",200,10,800,500);
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
		gr1gaus->SetMarkerSize(0.9);
		//new TCanvas(Form("Graph1 vs  %s ",graphTitle),"RMS vs L_{z}",200,10,800,500);
		//gr1->Draw();
	f->cd();
	gr2->SetName(Form("Mean%s",hist2D->GetName()));
	gr2->Write();
	gr1->SetName(Form("RMS%s",hist2D->GetName()));
	gr1->Write();
	gr1gaus->SetName(Form("Sigma%s ",graphTitle));
	gr1gaus->Write();
}

void ProjectionforBinMassInv(TH2F* hist2D,const char* graphTitle ,const char* Xaxisname, TCanvas *c0, TFile *f){
	double xaxis[100], mean[100], errmean[100], sigma[100], errsigma[100], sigmagaus[100], errsigmagaus[100];
	int count=0;
  	c0->Divide(6,5);
  	for(Int_t jj=1; jj<hist2D->GetNbinsY(); jj++){
 		printf("n=%d",hist2D->GetNbinsY());
    	TH1D *hInvMass = (TH1D *)hist2D->ProjectionX(Form("%sBin%d",hist2D->GetName(),jj), jj, jj);
    	if(hInvMass->GetSumOfWeights()>45){
			count++;
			c0->cd(count);
			c0->SetTitle(Form("%s",hist2D->GetName()));
			//new TCanvas();
			hInvMass->SetDirectory(0);
			
			hInvMass->GetXaxis()->SetRangeUser(1.05, 1.20);
    		hInvMass->Draw();
      		Double_t rmsh=hInvMass->GetRMS();
      		Double_t meanh=hInvMass->GetMean();
			Double_t Errrmsh=hInvMass->GetRMSError();
      		Double_t Errmeanh=hInvMass->GetMeanError();
			TF1* fg=((TF1*)(gROOT->GetFunction("gaus")));
			fg->SetParameter(1,1.116);
			fg->SetParameter(2,0.0015);
			fg->SetParameter(0,hInvMass->GetMaximum());
      		hInvMass->Fit(fg,"B","",meanh-5*rmsh,meanh+5*rmsh);
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
		//new TCanvas(Form("Graph vs  %s ",graphTitle),Form("Mean  %s ",graphTitle),200,10,800,500);
		//gr2->Draw();
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
		//new TCanvas(Form("Graph1 vs  %s ",graphTitle),Form("RMS  %s ",graphTitle),200,10,800,500);
		//gr1->Draw();
	f->cd();
	gr2->SetName(Form("Mean%s",hist2D->GetName()));
	gr2->Write();
	gr1->SetName(Form("RMS%s",hist2D->GetName()));
	gr1->Write();
	gr1gaus->SetName(Form("Sigma%s ",hist2D->GetName()));
	gr1gaus->Write();

}

