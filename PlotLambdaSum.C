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


void Projection(TH2F* h2D,  int range, bool rapidity, char *name );
void ProjectionResolCorrPy(TH2F* hResol, bool rapid, char *title);
void ProjectionNfake(TH2F* hNfake, char *title);
void ProjectionforBin(TH2F* hist2D, char* graphTitle , char* Xaxisname, TCanvas *c0, TFile *f);
void ProjectionforBinMassInv(TH2F* hist2D, char* graphTitle , char* Xaxisname, TCanvas *c0, TFile *f);
void efficiency(TH1F* hRec, TH1F* hMC, TH1F* hGen, char* xaxis, TFile *f);

void PlotLambdaSum(){

	TFile *hLambda=new TFile("LambdastandconBKG-Signal-histos.root");
	TFile *hLambdaBar=new TFile("LambdaBarstandconBKG-Signal-histos.root");
	TFile *hD0=new TFile("histogramD0.root");
	TFile *hPion=new TFile("PION-Signal-histos.root");
	TFile *hPionpos=new TFile("PIONpos-Signal-histos.root");
	TFile *hProt=new TFile("PROTON-Signal-histos.root");
	TFile *hAntiProt=new TFile("AntiPROTON-Signal-histos.root");

	TFile *foutLambda=new TFile("histogram.root","RECREATE");
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
	


	
/*	//Lambda
	//efficiency for y
	TH1F* hYMC=(TH1F*)hYPtLzMC->ProjectionX();
	TH1F* hYGen=(TH1F*)hYPtLzGen->ProjectionX();
	TH1F* hYRec=(TH1F*)hYPtLzRec->ProjectionX();
	efficiency(hYRec,hYMC,hYGen,"y #Lambda",foutLambda);

	//efficiency for Pt 	
	TH1F* hPtMC=(TH1F*)hYPtLzMC->ProjectionY();
	TH1F* hPtGen=(TH1F*)hYPtLzGen->ProjectionY();
	TH1F* hPtRec=(TH1F*)hYPtLzRec->ProjectionY();
	efficiency(hPtRec,hPtMC,hPtGen,"P_{t} #Lambda",foutLambda);

	//efficiency for Lz
	TH1F* hLzMC=(TH1F*)hYPtLzMC->ProjectionZ();
	TH1F* hLzGen=(TH1F*)hYPtLzGen->ProjectionZ();
	TH1F* hLzRec=(TH1F*)hYPtLzRec->ProjectionZ();
	efficiency(hLzRec,hLzMC,hLzGen,"Lz #Lambda",foutLambda);

*/
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
//Lambda Bar

/*	//Distribuzione 3D per efficienza Lambda Bar
	TH3F* hYPtLzMCBar=(TH3F*)hLambdaBar->Get("hYPtLzMC");
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
	efficiency(hLzRecBar,hLzMCBar,hLzGenBar,"Lz #Lambda Bar",foutLambda);
	
	TCanvas *effLnotL=new TCanvas("effLnotL");
	TH1F* eff[3], *effBar[3];
	eff[0]=(TH1F*)foutLambda->Get("Gen over GenMC for y #Lambda");
	eff[1]=(TH1F*)foutLambda->Get("Gen over GenMC for P_{t} #Lambda");	
	eff[2]=(TH1F*)foutLambda->Get("Gen over GenMC for Lz #Lambda");
	effBar[0]=(TH1F*)foutLambda->Get("Gen over GenMC for y #Lambda Bar");
	effBar[1]=(TH1F*)foutLambda->Get("Gen over GenMC for P_{t} #Lambda Bar");	
	effBar[2]=(TH1F*)foutLambda->Get("Gen over GenMC for Lz #Lambda Bar");

	effLnotL->Divide(3,1);
	for(int i=0;i<3;i++){
		eff[i]->SetDirectory(0);
		effBar[i]->SetDirectory(0);
		eff[i]->SetMarkerStyle(20);
		effBar[i]->SetMarkerStyle(20);
		eff[i]->SetMarkerSize(0.9);
		effBar[i]->SetMarkerSize(0.9);
		eff[i]->SetMarkerColor(1);
		effBar[i]->SetMarkerColor(2);
		
		effLnotL->cd(i+1);
		eff[i]->Draw();
		effBar[i]->Draw("SAME");
		gPad->SetGrid();
	}
	*/



	//Definizione file dove salvare tutti gli istogrammi nati a partire dalla ntupla

	TH2F* hNfakeVszP = new TH2F("hNfakeVszP", "Nfake vs zP", 6, -0.5, 5.5, 100, 0, 50);
	TH2F *hMassinreczP = new TH2F("hMassinrecVszP", "Invarian Mass vs zP", 3000, 0,3,100, 0, 50);
	TH2F *hMassinrecY = new TH2F("hMassinvVsY", "Invariant Mass vs Y", 3000, 0,3,60, 0, 6);
	TH2F *hMassinrecSwitchY = new TH2F("hMassinrecSwitchVsY", "Mass Switch vs Y", 3000, 0,3,60, 0, 6);
	TH2F *hMassinrecY0fake = new TH2F("hMassinvVsY0fake", "Invariant Mass vs Y 0fake", 3000, 0,3,60, 0, 6);
	TH2F *hMassinrecSwitchY0fake = new TH2F("hMassinrecSwitchVsY0fake", "Mass Switch vs Y 0fake", 3000, 0,3,60, 0, 6);

	TH2F *hResVxVsY = new TH2F("hResVxVsY", "Res Vx vs Y", 200, -200., 200., 60, 0, 6);
  	TH2F *hResVyVsY = new TH2F("hResVyVsY", "Res Vy vs Y", 200, -200., 200., 60, 0, 6);
  	TH2F *hResVzVsY = new TH2F("hResVzVsY", "Res Vz vs Y", 200, -3000., 3000., 60, 0, 6);
	TH2F *hResVxVsLz = new TH2F("hResVxVsLz", "Res Vx vs Lz", 200, -200., 200., 200, 0, 60);
  	TH2F *hResVyVsLz = new TH2F("hResVyVsLz", "Res Vy vs Lz", 200, -200., 200., 200, 0, 60);
  	TH2F *hResVzVsLz = new TH2F("hResVzVsLz", "Res Vz vs Lz", 200, -2000., 2000., 200, 0, 60);

	TH2F *hResPxVsYprot = new TH2F("hResPxVsYprot", "Res Px vs Y", 200, -0.2, 0.2, 60, 0, 6);
	TH2F *hResPyVsYprot = new TH2F("hResPyVsYprot", "Res Py vs Y", 200, -0.2, 0.2, 60, 0, 6);
	TH2F *hResPzVsYprot = new TH2F("hResPzVsYprot", "Res Pz vs Y", 200, -0.2, 0.2, 60, 0, 6);
	TH2F *hResPxVsYpion = new TH2F("hResPxVsYpion", "Res Px vs Y", 200, -0.08, 0.08, 60, 0, 6);
	TH2F *hResPyVsYpion = new TH2F("hResPyVsYpion", "Res Py vs Y", 200, -0.08, 0.08, 60, 0, 6);
	TH2F *hResPzVsYpion = new TH2F("hResPzVsYpion", "Res Pz vs Y", 200, -0.08, 0.08, 60, 0, 6);
	TH2F *hResPxVsYpionRel = new TH2F("hResPxVsYpionRel", "Res Px vs Y Relativa", 200, -0.08, 0.08, 60, 0, 6);
	TH2F *hResPyVsYpionRel = new TH2F("hResPyVsYpionRel", "Res Py vs Y Relativa", 200, -0.08, 0.08, 60, 0, 6);
	TH2F *hResPzVsYpionRel = new TH2F("hResPzVsYpionRel", "Res Pz vs Y Relativa", 200, -0.08, 0.08, 60, 0, 6);
	TH2F *hResPxVsYprotRel = new TH2F("hResPxVsYprotRel", "Res Px vs Y Relativa", 200, -0.08, 0.08, 60, 0, 6);
	TH2F *hResPyVsYprotRel = new TH2F("hResPyVsYprotRel", "Res Py vs Y Relativa", 200, -0.08, 0.08, 60, 0, 6);
	TH2F *hResPzVsYprotRel = new TH2F("hResPzVsYprotRel", "Res Pz vs Y Relativa", 200, -0.08, 0.08, 60, 0, 6);
	TH2F *hResPxVsLzprot = new TH2F("hResPxVsLzprot", "Res Px vs Lz", 200, -1., 1., 200, 0, 60);
	TH2F *hResPxVsLzpion = new TH2F("hResPxVsLzpion", "Res Px vs Lz", 200, -1., 1., 200, 0, 60);
	TH2F *hResPyVsLzprot = new TH2F("hResPyVsLzprot", "Res Py vs Lz", 200, -1., 1., 200, 0, 60);
	TH2F *hResPyVsLzpion = new TH2F("hResPyVsLzpion", "Res Py vs Lz", 200, -1., 1., 200, 0, 60);
	TH2F *hResPzVsLzprot = new TH2F("hResPzVsLzprot", "Res Pz vs Lz", 200, -1., 1., 200, 0, 60);
	TH2F *hResPzVsLzpion = new TH2F("hResPzVsLzpion", "Res Pz vs Lz", 200, -1., 1., 200, 0, 60);
	TH1F *hPxprot = new TH1F("hPxprot", " Px proton ", 200, -5., 5.);
	TH1F *hPyprot = new TH1F("hPyprot", " Py proton", 200, -2., 2.);
	TH1F *hPzprot = new TH1F("hPzprot", " Pz proton ", 200, -5., 50.);
	TH1F *hPxpion = new TH1F("hPxpion", " Px pion", 200, -5., 5.);
	TH1F *hPypion = new TH1F("hPypion", " Py pion", 200, -2., 2.);
	TH1F *hPzpion = new TH1F("hPzpion", " Pz pion", 200, -1., 20.);

	TH1F *htheta = new TH1F("htheta", " Angle between the particles", 200, 0., 7.);

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
	//Apertura file di input con variabili lambda per studio performance
	//Lettura Ntupla 
 	TFile filin("ntuplastandconBKG.root");
  	float nfaketrkprot, nfaketrkpion,y,xrec,yrec,zrec,secvertgenProt[3],massinvrec,pGenProt[3],pGenPion[3],dca,pRecProt[3],pRecPion[3],chi2prot,chi2ITSprot,chi2pion,chi2ITSpion,massinvSwitch,pGenLambda[3];
  	TNtuple *variables = (TNtuple*)filin.Get("nt");
		variables->SetBranchAddress("nfaketrkprot",&nfaketrkprot);
    	variables->SetBranchAddress("nfaketrkpion",&nfaketrkpion);
    	variables->SetBranchAddress("ygen",&y);
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
		variables->SetBranchAddress("pzlambda",&pGenLambda[0]);

	int events=variables->GetEntries();
	double theta=0;
	for(int i=0 ; i<events; i++) {
		variables->GetEvent(i);
		double residVx=10000.*(xrec - secvertgenProt[0]);
    	double residVy=10000.*(yrec - secvertgenProt[1]);
    	double residVz=10000.*(zrec - secvertgenProt[2]);
		if(zrec<=20.){
			hMassinreczP->Fill(massinvrec, zrec);
			hNfakeVszP->Fill(nfaketrkprot, zrec);
			}
		double p1abs=sqrt(pRecProt[0]*pRecProt[0]+pRecProt[1]*pRecProt[1]+pRecProt[2]*pRecProt[2]);
		double p2abs=sqrt(pRecPion[0]*pRecPion[0]+pRecPion[1]*pRecPion[1]+pRecPion[2]*pRecPion[2]);
		theta=TMath::ACos((pRecProt[0]*pRecPion[0]+pRecProt[1]*pRecPion[1]+pRecProt[2]*pRecPion[2])/(p1abs*p2abs));

		hMassinrecY->Fill(massinvrec, y);
		hMassinrecSwitchY->Fill(massinvSwitch, y);

		hNfakeProtVsChi2->Fill(nfaketrkprot,chi2prot);
		hNfakePionVsChi2->Fill(nfaketrkpion,chi2pion);

		htheta->Fill(theta);
			//hLzrecVsLzgen->Fill(secvertgenProt[2],zrec);
			//hDeltaLzVsLzgen->Fill(secvertgenProt[2],zrec-secvertgenProt[2]);
		
		if(nfaketrkprot == 0 && nfaketrkpion == 0) {
			hMassinrecY0fake->Fill(massinvrec, y);
			hMassinrecSwitchY0fake->Fill(massinvSwitch, y);
			//Resolution on P
			hResPxVsYpionRel->Fill((pRecPion[0]-pGenPion[0])/pRecPion[0],y);
			hResPyVsYpionRel->Fill((pRecPion[1]-pGenPion[1])/pRecPion[1],y);
			hResPzVsYpionRel->Fill((pRecPion[2]-pGenPion[2])/pRecPion[2],y);

			hResPxVsYprotRel->Fill((pRecProt[0]-pGenProt[0])/pRecProt[0],y);
			hResPyVsYprotRel->Fill((pRecProt[1]-pGenProt[1])/pRecProt[1],y);
			hResPzVsYprotRel->Fill((pRecProt[2]-pGenProt[2])/pRecProt[2],y);

			hResPxVsLzprot->Fill(pRecProt[0]-pGenProt[0],zrec);
			hResPyVsLzprot->Fill(pRecProt[1]-pGenProt[1],zrec);
			hResPzVsLzprot->Fill(pRecProt[2]-pGenProt[2],zrec);

			hResPxVsYprot->Fill(pRecProt[0]-pGenProt[0],y);
			hResPyVsYprot->Fill(pRecProt[1]-pGenProt[1],y);
			hResPzVsYprot->Fill(pRecProt[2]-pGenProt[2],y);

			hResPxVsLzpion->Fill(pRecPion[0]-pGenPion[0],zrec);
			hResPyVsLzpion->Fill(pRecPion[1]-pGenPion[1],zrec);
			hResPzVsLzpion->Fill(pRecPion[2]-pGenPion[2],zrec);

			hPxprot->Fill(pRecProt[0]);
			hPyprot->Fill(pRecProt[1]);
			hPzprot->Fill(pRecProt[2]);
			hPxpion->Fill(pRecPion[0]);
			hPypion->Fill(pRecPion[1]);
			hPzpion->Fill(pRecPion[2]);

			//Resolution on Vertex
			hResVxVsY->Fill(residVx, y);
			hResVyVsY->Fill(residVy, y);
			hResVzVsY->Fill(residVz, y);
			hResVxVsLz->Fill(residVx, zrec);
			hResVyVsLz->Fill(residVy, zrec);
			hResVzVsLz->Fill(residVz, zrec);

			hLzrecVsLzgen->Fill(secvertgenProt[2],zrec);
			hDeltaLzVsLzgen->Fill(secvertgenProt[2],zrec-secvertgenProt[2]);
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

	TFile filinBar("ntuplaLambdaBarstandconBKG.root");
  	float nfaketrkprotBar, nfaketrkpionBar,yBar,xrecBar,yrecBar,zrecBar,secvertgenProtBar[3],massinvrecBar,pGenProtBar[3],pGenPionBar[3],dcaBar,pRecProtBar[3],pRecPionBar[3],chi2protBar,chi2ITSprotBar,chi2pionBar,chi2ITSpionBar,massinvSwitchBar,pGenLambdaBar[3];
  	TNtuple *variablesBar = (TNtuple*)filinBar.Get("ntLambdaBarstand");
		variablesBar->SetBranchAddress("nfaketrkprot",&nfaketrkprotBar);
    	variablesBar->SetBranchAddress("nfaketrkpion",&nfaketrkpionBar);
    	variablesBar->SetBranchAddress("ygen",&yBar);
    	variablesBar->SetBranchAddress("xP",&xrecBar);
		variablesBar->SetBranchAddress("yP",&yrecBar);
		variablesBar->SetBranchAddress("zP",&zrecBar);
		variablesBar->SetBranchAddress("Vxgen",&secvertgenProtBar[0]);
		variablesBar->SetBranchAddress("Vygen",&secvertgenProtBar[1]);
		variablesBar->SetBranchAddress("Vzgen",&secvertgenProtBar[2]);
		variablesBar->SetBranchAddress("massinvrec",&massinvrecBar);
		variables->SetBranchAddress("pxgenprot",&pGenProtBar[0]);
		variables->SetBranchAddress("pygenprot",&pGenProtBar[1]);
		variables->SetBranchAddress("pzgenprot",&pGenProtBar[2]);
		variables->SetBranchAddress("pxgenpion",&pGenPionBar[0]);
		variables->SetBranchAddress("pygenpion",&pGenPionBar[1]);
		variables->SetBranchAddress("pzgenpion",&pGenPionBar[2]);
		variablesBar->SetBranchAddress("dca",&dcaBar);
		variables->SetBranchAddress("chi2prot",&chi2protBar);
		variables->SetBranchAddress("chi2ITSprot",&chi2ITSprotBar);
		variables->SetBranchAddress("chi2pion",&chi2pionBar);
		variables->SetBranchAddress("chi2ITSpion",&chi2ITSpionBar);
		variables->SetBranchAddress("pxrecprot",&pRecProtBar[0]);
		variables->SetBranchAddress("pyrecprot",&pRecProtBar[1]);
		variables->SetBranchAddress("pzrecprot",&pRecProtBar[2]);
		variables->SetBranchAddress("pxrecpion",&pRecPionBar[0]);
		variables->SetBranchAddress("pyrecpion",&pRecPionBar[1]);
		variables->SetBranchAddress("pzrecpion",&pRecPionBar[2]);
		variables->SetBranchAddress("massinvswitch",&massinvSwitchBar);
		variables->SetBranchAddress("pxlambda",&pGenLambdaBar[0]);
		variables->SetBranchAddress("pylambda",&pGenLambdaBar[1]);
		variables->SetBranchAddress("pzlambda",&pGenLambdaBar[0]);

	int eventsBar=variablesBar->GetEntries();

	for(int i=0 ; i<eventsBar; i++) {
		variablesBar->GetEvent(i);
		double residVxBar=10000.*(xrecBar - secvertgenProtBar[0]);
		double residVyBar=10000.*(yrecBar - secvertgenProtBar[1]);
		double residVzBar=10000.*(zrecBar - secvertgenProtBar[2]);
		if(zrecBar<=20.){
			hMassinreczP->Fill(massinvrecBar, zrecBar);
			hNfakeVszP->Fill(nfaketrkprotBar, zrecBar);
		}
		double p1absBar=sqrt(pRecProtBar[0]*pRecProtBar[0]+pRecProtBar[1]*pRecProtBar[1]+pRecProtBar[2]*pRecProtBar[2]);
		double p2absBar=sqrt(pRecPionBar[0]*pRecPionBar[0]+pRecPionBar[1]*pRecPionBar[1]+pRecPionBar[2]*pRecPionBar[2]);
		theta=TMath::ACos((pRecProtBar[0]*pRecPionBar[0]+pRecProtBar[1]*pRecPionBar[1]+pRecProtBar[2]*pRecPionBar[2])/(p1absBar*p2absBar));
		
		hMassinrecY->Fill(massinvrecBar, yBar);
		hMassinrecSwitchY->Fill(massinvSwitchBar, yBar);

		htheta->Fill(theta);
		hNfakeProtVsChi2->Fill(nfaketrkprotBar,chi2protBar);
		hNfakePionVsChi2->Fill(nfaketrkpionBar,chi2pionBar);

		//hLzrecVsLzgen->Fill(secvertgenProtBar[2],zrecBar);
		//hDeltaLzVsLzgen->Fill(secvertgenProtBar[2],zrecBar-secvertgenProtBar[2]);

		if(nfaketrkprotBar == 0 && nfaketrkpionBar == 0) {
			hMassinrecY0fake->Fill(massinvrecBar, yBar);
			hMassinrecSwitchY0fake->Fill(massinvSwitchBar, yBar);

			//Resolution on P
			hResPxVsYpionRel->Fill((pRecPionBar[0]-pGenPionBar[0])/pRecPionBar[0],yBar);
			hResPyVsYpionRel->Fill((pRecPionBar[1]-pGenPionBar[1])/pRecPionBar[1],yBar);
			hResPzVsYpionRel->Fill((pRecPionBar[2]-pGenPionBar[2])/pRecPionBar[2],yBar);

			hResPxVsYprotRel->Fill((pRecProtBar[0]-pGenProtBar[0])/pRecProtBar[0],yBar);
			hResPyVsYprotRel->Fill((pRecProtBar[1]-pGenProtBar[1])/pRecProtBar[1],yBar);
			hResPzVsYprotRel->Fill((pRecProtBar[2]-pGenProtBar[2])/pRecProtBar[2],yBar);

			hResPxVsLzprot->Fill(pRecProtBar[0]-pGenProtBar[0],zrecBar);
			hResPyVsLzprot->Fill(pRecProtBar[1]-pGenProtBar[1],zrecBar);
			hResPzVsLzprot->Fill(pRecProtBar[2]-pGenProtBar[2],zrecBar);
			hResPxVsYprot->Fill(pRecProtBar[0]-pGenProtBar[0],yBar);
			hResPyVsYprot->Fill(pRecProtBar[1]-pGenProtBar[1],yBar);
			hResPzVsYprot->Fill(pRecProtBar[2]-pGenProtBar[2],yBar);

			hResPxVsLzpion->Fill(pRecPionBar[0]-pGenPionBar[0],zrecBar);
			hResPyVsLzpion->Fill(pRecPionBar[1]-pGenPionBar[1],zrecBar);
			hResPzVsLzpion->Fill(pRecPionBar[2]-pGenPionBar[2],zrecBar);
			hResPxVsYpion->Fill(pRecPionBar[0]-pGenPionBar[0],yBar);
			hResPyVsYpion->Fill(pRecPionBar[1]-pGenPionBar[1],yBar);
			hResPzVsYpion->Fill(pRecPionBar[2]-pGenPionBar[2],yBar);
		
			hPxprot->Fill(pRecProtBar[0]);
			hPyprot->Fill(pRecProtBar[1]);
			hPzprot->Fill(pRecProtBar[2]);
			hPxpion->Fill(pRecPionBar[0]);
			hPypion->Fill(pRecPionBar[1]);
			hPzpion->Fill(pRecPionBar[2]);

			//Resolution on vertex
			hResVxVsY->Fill(residVxBar, yBar);
			hResVyVsY->Fill(residVyBar, yBar);
			hResVzVsY->Fill(residVzBar, yBar);
			hResVxVsLz->Fill(residVxBar, zrecBar);
			hResVyVsLz->Fill(residVyBar, zrecBar);
			hResVzVsLz->Fill(residVzBar, zrecBar);
			hLzrecVsLzgen->Fill(secvertgenProtBar[2],zrecBar);
			hDeltaLzVsLzgen->Fill(secvertgenProtBar[2],zrecBar-secvertgenProtBar[2]);

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
	
	foutLambda->cd();
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
	hResVxVsY->Write();
  	hResVyVsY->Write();
  	hResVzVsY->Write();
	hResVxVsLz->Write();
  	hResVyVsLz->Write();
  	hResVzVsLz->Write();
	hResPxVsLzprot->Write();
	hResPyVsLzprot->Write();
	hResPxVsLzpion->Write();
	hResPyVsLzpion->Write();
	hResPzVsLzprot->Write();
	hResPzVsLzpion->Write();

	hResPxVsYprot->Write();
	hResPxVsYpion->Write();
	hResPyVsYprot->Write();
	hResPyVsYpion->Write();
	hResPzVsYprot->Write();
	hResPzVsYpion->Write();
	hMassinreczP->Write();

	hResPxVsYpionRel->Write();
	hResPyVsYpionRel->Write();
	hResPzVsYpionRel->Write();

	hResPxVsYprotRel->Write();
	hResPyVsYprotRel->Write();
	hResPzVsYprotRel->Write();

	hPxprot->Write();
	hPyprot->Write();
	hPzprot->Write();
	hPxpion->Write();
	hPypion->Write();
	hPzpion->Write();
	
	htheta->Write();
	
	hNfakeProtVsChi2->Write();
	hNfakePionVsChi2->Write();

	hLzrecVsLzgen->Write();
	hDeltaLzVsLzgen->Write();

	TH1D* hInvMass0fake=(TH1D *)hMassinrecY0fake->ProjectionX();
	hInvMass0fake->SetDirectory(0);
	//new TCanvas();
	hInvMass0fake->SetLineColor(2);
	
	TH1D* hInvMassSwitch0fake=(TH1D *)hMassinrecSwitchY0fake->ProjectionX();
	hInvMassSwitch0fake->SetDirectory(0);
	
//new TCanvas();
	//hInvMassSwitch0fake->Draw();
//hInvMass0fake->Draw();


	/*TH1D* MassInvProjswap[27];
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
	new TCanvas();
	hsumSwap->Draw();
	new TCanvas();
	hsum->Draw();

	new TCanvas();
	hsumSwap->Draw();
	hsum->Draw("same");
*/
	TCanvas* c1=new TCanvas("c1");
	ProjectionforBinMassInv(hMassinrecY, "#Lambda+#bar{#Lambda} Inv Mass vs Y", "Y", c1, foutLambda);

	TCanvas* c2p=new TCanvas("c2p");
	ProjectionforBin(hResPxVsYprot ,"#Lambda+#bar{#Lambda} ResPx vs Y prot", "Y", c2p, foutLambda);
	TCanvas* c3p=new TCanvas("c3p");
	ProjectionforBin(hResPyVsYprot, "#Lambda+#bar{#Lambda} ResPy vs Y prot", "Y", c3p, foutLambda);
	TCanvas* c4p=new TCanvas("c4p");
	ProjectionforBin(hResPzVsYprot, "#Lambda+#bar{#Lambda} ResPz vs Y prot", "Y", c4p, foutLambda);
	TCanvas* c5p=new TCanvas("c5p");
	ProjectionforBin(hResPxVsYpion, "#Lambda+#bar{#Lambda} ResPx vs Y pion", "Y", c5p, foutLambda);
	TCanvas* c6p=new TCanvas("c6p");
	ProjectionforBin(hResPyVsYpion, "#Lambda+#bar{#Lambda} ResPy vs Y pion", "Y", c6p, foutLambda);
	TCanvas* c7p=new TCanvas("c7p");
	ProjectionforBin(hResPzVsYpion, "#Lambda+#bar{#Lambda} ResPz vs Y pion", "Y", c7p, foutLambda);

	TCanvas* c5piRel=new TCanvas("c5piRel");
	ProjectionforBin(hResPxVsYpionRel, "#Lambda+#bar{#Lambda} pion ResPx vs Y Rel", "Y", c5piRel, foutLambda);
	TCanvas* c6piRel=new TCanvas("c6piRel");
	ProjectionforBin(hResPyVsYpionRel, "#Lambda+#bar{#Lambda} pion ResPy vs Y Rel", "Y", c6piRel, foutLambda);
	TCanvas* c7piRel=new TCanvas("c7piRel");
	ProjectionforBin(hResPzVsYpionRel, "#Lambda+#bar{#Lambda} pion ResPz vs Y Rel", "Y", c7piRel, foutLambda);

	TCanvas* c5pRel=new TCanvas("c5pRel");
	ProjectionforBin(hResPxVsYprotRel, "#Lambda+#bar{#Lambda} prot ResPx vs Y Rel", "Y", c5pRel, foutLambda);
	TCanvas* c6pRel=new TCanvas("c6pRel");
	ProjectionforBin(hResPyVsYprotRel, "#Lambda+#bar{#Lambda} prot ResPy vs Y Rel", "Y", c6pRel, foutLambda);
	TCanvas* c7pRel=new TCanvas("c7pRel");
	ProjectionforBin(hResPzVsYprotRel, "#Lambda+#bar{#Lambda} prot ResPz vs Y Rel", "Y", c7pRel, foutLambda);

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
	mgp[0]->SetTitle("Relative Resolutions on P_{x} for #Lambda daughters ; y ; #sigma");
	mgp[1]->SetTitle("Relative Resolutions on P_{y} for #Lambda daughters; y ; #sigma");
	mgp[2]->SetTitle("Relative Resolutions on P_{z} for #Lambda daughters; y ; #sigma");
	TGraphErrors *pprot[3], *ppion[3], *pkaonD0[3], *ppionD0[3];
	pprot[0]=(TGraphErrors*)foutLambda->Get("Sigma#Lambda+#bar{#Lambda} prot ResPx vs Y Rel");
	pprot[1]=(TGraphErrors*)foutLambda->Get("Sigma#Lambda+#bar{#Lambda} prot ResPy vs Y Rel");
	pprot[2]=(TGraphErrors*)foutLambda->Get("Sigma#Lambda+#bar{#Lambda} prot ResPz vs Y Rel");
	ppion[0]=(TGraphErrors*)foutLambda->Get("Sigma#Lambda+#bar{#Lambda} pion ResPx vs Y Rel");
	ppion[1]=(TGraphErrors*)foutLambda->Get("Sigma#Lambda+#bar{#Lambda} pion ResPy vs Y Rel");
	ppion[2]=(TGraphErrors*)foutLambda->Get("Sigma#Lambda+#bar{#Lambda} pion ResPz vs Y Rel");
	pkaonD0[0]=(TGraphErrors*)hD0->Get("SigmahResPxVsYKaonRel");
	pkaonD0[1]=(TGraphErrors*)hD0->Get("SigmahResPyVsYKaonRel");
	pkaonD0[2]=(TGraphErrors*)hD0->Get("SigmahResPzVsYKaonRel");
	ppionD0[0]=(TGraphErrors*)hD0->Get("SigmahResPxVsYpionRel");
	ppionD0[1]=(TGraphErrors*)hD0->Get("SigmahResPyVsYpionRel");
	ppionD0[2]=(TGraphErrors*)hD0->Get("SigmahResPzVsYpionRel");
	for(int i=0;i<3;i++){
		resP->cd(i+1);
		gPad->SetGrid();
		pprot[i]->SetMarkerColor(1);
		ppion[i]->SetMarkerColor(2);
		pprot[i]->SetMarkerStyle(20);
		ppion[i]->SetMarkerStyle(20);
		pprot[i]->SetMarkerSize(0.9);
		ppion[i]->SetMarkerSize(0.9);
		mgp[i]->Add(pprot[i],"p");
		mgp[i]->Add(ppion[i],"p");
		mgp[i]->GetXaxis()->SetTitle("y");
		mgp[i]->GetYaxis()->SetTitle("#sigma");
		mgp[i]->Draw("a");
	}
	resP->cd(1);
	auto legend = new TLegend(0.1,0.7,0.48,0.9);
   	//legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
   	legend->AddEntry(pprot[0],"#Lambda proton","lep");
  	legend->AddEntry(ppion[0],"#Lambda pion","lep");

   	legend->Draw();


	//mg->SetTitle("Different resolutions on P_{x}; y ; #sigma (#frac{Gev}/{c})");
	TCanvas* resPLD=new TCanvas("resPLD");
	resPLD->Divide(3,1);
	//resPLD->SetTitle("pippo");
	//TLatex * tex = new TLatex(0.5,0.8,"My Title"); 
	//tex->SetNDC(); 
	//tex->SetTextSize(0.2); 
	//tex->Draw(); 

	for(int i=0;i<3;i++){
		resPLD->cd(i+1);
		gPad->SetGrid();
		pprot[i]->SetMarkerColor(1);
		ppion[i]->SetMarkerColor(2);
		pkaonD0[i]->SetMarkerColor(4);
		ppionD0[i]->SetMarkerColor(7);
		pprot[i]->SetMarkerStyle(20);
		ppion[i]->SetMarkerStyle(20);
		pkaonD0[i]->SetMarkerStyle(20);
		ppionD0[i]->SetMarkerStyle(20);
		pprot[i]->SetMarkerSize(0.8);
		ppion[i]->SetMarkerSize(0.8);
		pkaonD0[i]->SetMarkerSize(0.8);
		ppionD0[i]->SetMarkerSize(0.8);
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
	
 
		
	
/*	TH2F* hResPxPionRel=(TH2F*)hPion->Get("hResPxVsYpionRel");
	TH2F* hResPyPionRel=(TH2F*)hPion->Get("hResPyVsYpionRel");
	TH2F* hResPzPionRel=(TH2F*)hPion->Get("hResPzVsYpionRel");

	TCanvas* c5Pion=new TCanvas("c5Pion");
	ProjectionforBin(hResPxPionRel, "#pi - ResPx vs Y relativa", "Y", c5Pion, foutLambda);
	TCanvas* c6Pion=new TCanvas("c6Pion");
	ProjectionforBin(hResPyPionRel, "#pi - ResPy vs Y relativa", "Y", c6Pion, foutLambda);
	TCanvas* c7Pion=new TCanvas("c7Pion");
	ProjectionforBin(hResPzPionRel, "#pi - ResPz vs Y relativa", "Y", c7Pion, foutLambda);

	TH2F* hResPxPionposRel=(TH2F*)hPionpos->Get("hResPxVsYpionRel");
	TH2F* hResPyPionposRel=(TH2F*)hPionpos->Get("hResPyVsYpionRel");
	TH2F* hResPzPionposRel=(TH2F*)hPionpos->Get("hResPzVsYpionRel");

	TCanvas* c5Pionpos=new TCanvas("c5Pionpos");
	ProjectionforBin(hResPxPionposRel, "#pi + ResPx vs Y relativa", "Y", c5Pionpos, foutLambda);
	TCanvas* c6Pionpos=new TCanvas("c6Pionpos");
	ProjectionforBin(hResPyPionposRel, "#pi + ResPy vs Y relativa", "Y", c6Pionpos, foutLambda);
	TCanvas* c7Pionpos=new TCanvas("c7Pionpos");
	ProjectionforBin(hResPzPionposRel, "#pi + ResPz vs Y relativa", "Y", c7Pionpos, foutLambda);

	TH2F* hResPxProtRel=(TH2F*)hProt->Get("hResPxVsYprotRel");
	TH2F* hResPyProtRel=(TH2F*)hProt->Get("hResPyVsYprotRel");
	TH2F* hResPzProtRel=(TH2F*)hProt->Get("hResPzVsYprotRel");

	TCanvas* c5Prot=new TCanvas("c5Prot");
	ProjectionforBin(hResPxProtRel, "proton ResPx vs Y relativa", "Y", c5Prot, foutLambda);
	TCanvas* c6Prot=new TCanvas("c6Prot");
	ProjectionforBin(hResPyProtRel, "proton ResPy vs Y relativa", "Y", c6Prot, foutLambda);
	TCanvas* c7Prot=new TCanvas("c7Prot");
	ProjectionforBin(hResPzProtRel, "proton ResPz vs Y relativa", "Y", c7Prot, foutLambda);

	TH2F* hResPxAntiProtRel=(TH2F*)hAntiProt->Get("hResPxVsYprotRel");
	TH2F* hResPyAntiProtRel=(TH2F*)hAntiProt->Get("hResPyVsYprotRel");
	TH2F* hResPzAntiProtRel=(TH2F*)hAntiProt->Get("hResPzVsYprotRel");	

	TCanvas* c5AntiProt=new TCanvas("c5AntiProt");
	ProjectionforBin(hResPxAntiProtRel, "Antiproton ResPx vs Y relativa", "Y", c5AntiProt, foutLambda);
	TCanvas* c6AntiProt=new TCanvas("c6AntiProt");
	ProjectionforBin(hResPyAntiProtRel, "Antiproton ResPy vs Y relativa", "Y", c6AntiProt, foutLambda);
	TCanvas* c7AntiProt=new TCanvas("c7AntiProt");
	ProjectionforBin(hResPzAntiProtRel, "Antiproton ResPz vs Y relativa", "Y", c7AntiProt, foutLambda);


	
*/

	/*TGraphErrors *pxpionLambda=(TGraphErrors*)foutLambda->Get("Sigma#Lambda+#bar{#Lambda} pion ResPx vs Y Rel");
	TGraphErrors *pypionLambda=(TGraphErrors*)foutLambda->Get("Sigma#Lambda+#bar{#Lambda} pion ResPy vs Y Rel");
	TGraphErrors *pxprotLambda=(TGraphErrors*)foutLambda->Get("Sigma#Lambda+#bar{#Lambda} prot ResPx vs Y Rel");
	TGraphErrors *pyprotLambda=(TGraphErrors*)foutLambda->Get("Sigma#Lambda+#bar{#Lambda} prot ResPy vs Y Rel");
	TGraphErrors *pxpion=(TGraphErrors*)foutLambda->Get("Sigma#pi - ResPx vs Y relativa");
	TGraphErrors *pypion=(TGraphErrors*)foutLambda->Get("Sigma#pi - ResPy vs Y relativa");
	TGraphErrors *pxpionpos=(TGraphErrors*)foutLambda->Get("Sigma#pi + ResPx vs Y relativa");
	TGraphErrors *pypionpos=(TGraphErrors*)foutLambda->Get("Sigma#pi + ResPy vs Y relativa");
	TGraphErrors *pxprot=(TGraphErrors*)foutLambda->Get("Sigmaproton ResPx vs Y relativa");
	TGraphErrors *pyprot=(TGraphErrors*)foutLambda->Get("Sigmaproton ResPy vs Y relativa");
	TGraphErrors *pxantiprot=(TGraphErrors*)foutLambda->Get("SigmaAntiproton ResPx vs Y relativa");
	TGraphErrors *pyantiprot=(TGraphErrors*)foutLambda->Get("SigmaAntiproton ResPy vs Y relativa");
	
	TCanvas *deltaPrelX=new TCanvas("deltaPrelX");
	deltaPrelX->SetGrid();
	TMultiGraph *mg = new TMultiGraph();
	mg->SetTitle("Different resolutions on P_{x}; y ; #sigma (#frac{Gev}/{c})");
	pxpion->SetMarkerColor(1);
	pxpionpos->SetMarkerColor(2);
	pxprot->SetMarkerColor(4);
	pxpionLambda->SetMarkerColor(8);
	pxantiprot->SetMarkerColor(6);
	pxprotLambda->SetMarkerColor(5);
	pxpion->SetMarkerStyle(20);
	pxpionpos->SetMarkerStyle(20);
	pxprot->SetMarkerStyle(20);
	pxpionLambda->SetMarkerStyle(20);
	pxantiprot->SetMarkerStyle(20);
	pxprotLambda->SetMarkerStyle(20);
	pxpion->SetMarkerSize(0.9);
	pxpionpos->SetMarkerSize(0.9);
	pxprot->SetMarkerSize(0.9);
	pxpionLambda->SetMarkerSize(0.9);
	pxantiprot->SetMarkerSize(0.9);
	pxprotLambda->SetMarkerSize(0.9);
	mg->Add(pxpion,"p");
	mg->Add(pxpionpos,"p");
	mg->Add(pxprot,"p");
	mg->Add(pxantiprot,"p");
	mg->Add(pxprotLambda,"p");
	mg->Add(pxpionLambda,"p");
	
	//mg->GetXaxis()->SetTitle("y");
	//mg->GetYaxis()->SetTitle("#sigma");
	mg->Draw("a");

	TCanvas *deltaPrelY=new TCanvas("deltaPrelY");
	deltaPrelY->SetGrid();
	TMultiGraph *mgyrel = new TMultiGraph();
	pypion->SetMarkerColor(1);
	pypionpos->SetMarkerColor(2);
	pyprot->SetMarkerColor(4);
	pyantiprot->SetMarkerColor(6);
	pyprotLambda->SetMarkerColor(5);
	pypionLambda->SetMarkerColor(8);
	pypion->SetMarkerStyle(20);
	pypionpos->SetMarkerStyle(20);
	pyprot->SetMarkerStyle(20);
	pypionLambda->SetMarkerStyle(20);
	pyantiprot->SetMarkerStyle(20);
	pyprotLambda->SetMarkerStyle(20);
	pypion->SetMarkerSize(0.9);
	pypionpos->SetMarkerSize(0.9);
	pyprot->SetMarkerSize(0.9);
	pypionLambda->SetMarkerSize(0.9);
	pyantiprot->SetMarkerSize(0.9);
	pyprotLambda->SetMarkerSize(0.9);
	mgyrel->Add(pypion,"p");
	mgyrel->Add(pypionpos,"p");
	mgyrel->Add(pyprot,"p");
	mgyrel->Add(pyantiprot,"p");
	mgyrel->Add(pyprotLambda,"p");
	mgyrel->Add(pypionLambda,"p");
	mgyrel->GetXaxis()->SetTitle("y");
	mgyrel->GetYaxis()->SetTitle("#sigma");
	mgyrel->Draw("a");*/

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
	
	hLambda->Close();
	hD0->Close();
	foutLambda->Close();

}


void Projection(TH2F* h2D,  int range, bool rapidity, char *name ){
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



void ProjectionResolCorrPy(TH2F* hResol, bool rapid, char *title){
	int totvalue=10;
	double meanPy[10]={0.}, RMSPy[10]={0.}, meanErrPy[10]={0.}, RMSErrPy[10]={0.}, Lzrange[10]={0.};
	TH1D* hProjResyLz[10];
	double jrange=0.;
	if(rapid==1){
		jrange=1;
		totvalue=4;
	}
	for(int i=0;i<totvalue;i++){
		hProjResyLz[i]=hResol->ProjectionY((Form("Res Py Lz range [%f %f]",jrange,jrange+1)), hResol->GetXaxis()->FindBin(jrange+0.01),hResol->GetXaxis()->FindBin(jrange+0.99));
		meanPy[i]=hProjResyLz[i]->GetMean();
		RMSPy[i]=hProjResyLz[i]->GetRMS();
		meanErrPy[i]=hProjResyLz[i]->GetMeanError();
		RMSErrPy[i]=hProjResyLz[i]->GetRMSError();
		Lzrange[i]=(jrange+jrange+1)/2.;
		hProjResyLz[i]->SetDirectory(0);
		new TCanvas();
		hProjResyLz[i]->Draw();
		jrange=jrange+1;
		}
	if(rapid==0){
		TGraphErrors* grResy = new TGraphErrors(totvalue,Lzrange,meanPy,0,meanErrPy);
		grResy->SetTitle(Form("Mean %s  ", title));
		new TCanvas();
		grResy->Draw("AC*");
		TGraphErrors* grRMSy = new TGraphErrors(totvalue,Lzrange,RMSPy,0,RMSErrPy);
		grRMSy->SetTitle(Form("RMS %s  ",title));
		new TCanvas();
		grRMSy->Draw("AC*");
			}
	else{
		TGraphErrors* grResyRapid = new TGraphErrors(totvalue,Lzrange,meanPy,0,meanErrPy);
		grResyRapid->SetTitle(Form("Mean %s ",title));
		//new TCanvas();
		//grResyRapid->Draw("AC*");
		TGraphErrors* grRMSyRapid = new TGraphErrors(totvalue,Lzrange,RMSPy,0,RMSErrPy);
		grRMSyRapid->SetTitle(Form("RMS %s ",title));
		//new TCanvas();
		//grRMSyRapid->Draw("AC*");
		}
	}

void ProjectionNfake(TH2F* hNfake, char *title){
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

void ProjectionforBin(TH2F* hist2D, char* graphTitle , char* Xaxisname, TCanvas *c0, TFile *f){
	double xaxis[100], mean[100], errmean[100], sigma[100], errsigma[100], sigmagaus[100], errsigmagaus[100];
	int count=0, i=0, countGraph=0;
  	c0->Divide(8,5);
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

void ProjectionforBinMassInv(TH2F* hist2D, char* graphTitle , char* Xaxisname, TCanvas *c0, TFile *f){
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
			//hInvMass->SetMinimum(0);
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

void efficiency(TH1F* hRec, TH1F* hMC, TH1F* hGen, char* xaxis, TFile* f){
	TH1F *Rec_Gen = (TH1F*)hRec->Clone("Rec_Gen");
	Rec_Gen->SetDirectory(0);
	Rec_Gen->Divide(hRec,hMC,1,1,"B");
	TCanvas *eff_y=new TCanvas(Form("eff_%s",xaxis));
	Rec_Gen->SetTitle(Form("Rec over GenMC for %s",xaxis));
	Rec_Gen->GetXaxis()->SetTitle(Form("%s",xaxis));
	Rec_Gen->GetYaxis()->SetTitle("Efficiency #epsilon");
	Rec_Gen->Draw("E1");
	TH1F *Gen_GenMC = (TH1F*)hGen->Clone("Gen_GenMC");
	Gen_GenMC->SetDirectory(0);
	Gen_GenMC->Divide(hGen,hMC,1,1,"B");
	TCanvas *effGen=new TCanvas(Form("effGen%s",xaxis));
	Gen_GenMC->SetTitle(Form("Gen over GenMC for %s",xaxis));
	Gen_GenMC->GetXaxis()->SetTitle(Form("%s",xaxis));
	Gen_GenMC->GetYaxis()->SetTitle("Efficiency #epsilon");
	Gen_GenMC->Draw("E1");
	f->cd();
	Gen_GenMC->SetName(Form("Gen over GenMC for %s",xaxis));
	Gen_GenMC->Write();

}
