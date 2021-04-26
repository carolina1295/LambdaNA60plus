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

void Projection(TH2F* h2D,  int range, bool rapidity, char *name );
void ProjectionResol(TH2F* hResol, bool rapid, char *title);
void ProjectionNfake(TH2F* hNfake, char *title);
void ProjectionforBin(TH2F* hist2D, char* graphTitle , char* Xaxisname, TCanvas *c0, TFile *f);
void ProjectionforBinMassInv(TH2F* hist2D, char* graphTitle , char* Xaxisname, TCanvas *c0, TFile *f);

void PlotLambda(){

TFile *hLambda=new TFile("Lambda-Signal-histos.root");
TFile *hD0=new TFile("histogramD0.root");
TFile *fLBar=new TFile("histogramLambdaBar.root");
//TFile *hLambdaPYTHIA=new TFile("Lambda-Signal-histosPYTHIA.root");
	
//Checking pythia and evtgen results
/*	TH1D *hmasstrue = (TH1D*)hLambda->Get("hMassLambdatrue");
	hmasstrue->SetDirectory(0);
	TH1D *hmasstruePYTHIA = (TH1D*)hLambdaPYTHIA->Get("hMassLambdatrue");
	hmasstruePYTHIA->SetDirectory(0);
	TH1D *hmassfalse = (TH1D*)hLambda->Get("hMassLambdafalse");
	hmassfalse->SetDirectory(0);
	TH1D *hmassfalsePYTHIA = (TH1D*)hLambdaPYTHIA->Get("hMassLambdafalse");
	hmassfalsePYTHIA->SetDirectory(0);
	hmasstrue->SetLineColor(2);
  hmasstruePYTHIA->SetLineColor(4);
	new TCanvas();
	hmasstrue->Draw();
	hmasstrue->Draw("SAME");
	hmassfalse->SetLineColor(2);
  hmassfalsePYTHIA->SetLineColor(4);
	new TCanvas();
	hmassfalse->Draw();
	hmassfalse->Draw("SAME");*/

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
  ProjectionResol(hResPyVsLzprot,1);
	TH2F *hResPyVsLzpion = (TH2F*)hLambda->Get("hResPyVsLzpion");
 	ProjectionResol(hResPyVsLzpion,0);*/

	//Distribuzioni per risoluzione per projection py dopo il passaggio nel rivelatore DOPO CORREZIONE DEL VERTICE
	TH2F *hResPyVsLzprotcorr = (TH2F*)hLambda->Get("hResPyVsLzprotcorr");
  ProjectionResol(hResPyVsLzprotcorr,0," ResPy Vs Lz (correction) Proton");
	TH2F *hResPyVsLzpioncorr = (TH2F*)hLambda->Get("hResPyVsLzpioncorr");
 	ProjectionResol(hResPyVsLzpioncorr,0," ResPy Vs Lz (correction) Pion");
	TH2F *hResPyVsYprotcorr = (TH2F*)hLambda->Get("hResPyVsYprotcorr");
  ProjectionResol(hResPyVsYprotcorr,1,"ResPy Vs Y (correction) Proton");
	TH2F *hResPyVsYpioncorr = (TH2F*)hLambda->Get("hResPyVsYpioncorr");
 	ProjectionResol(hResPyVsYpioncorr,1,"ResPy Vs Y (correction) Pion");

//Distribuzioni nfake per projection py dopo il passaggio nel rivelatore DOPO CORREZIONE DEL VERTICE
	TH2F *hnfakeVsLzprot = (TH2F*)hLambda->Get("hnfakeVsLzprot");
  ProjectionNfake(hnfakeVsLzprot," Nfake Vs Lz Proton ");
	
//Apertura file di input con variabili lambda per studio performance

//Lettura Ntupla 
 	TFile filin("ntupla.root");
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
	//Definizione file dove salvare tutti gli istogrammi nati a partire dalla ntupla
	TFile *foutLambda=new TFile("histogramLambda.root","RECREATE");
	//printf("Check events tree(Nfake) = %d and ntupla = %d \n",evNfake,events);
	TH2F* hNfakeVszP = new TH2F("hNfakeVszP", "Nfake vs zP", 6, -0.5, 5.5, 100, 0, 50);
	TH2F *hMassinreczP = new TH2F("hMassinrecVszP", "Mass vs zP", 3000, 0,3,100, 0, 50);
	TH2F *hMassinrecY = new TH2F("hMassinrecVsY", "Mass vs Y", 3000, 0,3,60, 0, 6);	
	TH2F *hResVxVsY = new TH2F("hResVxVsY", "Res Vx vs Y", 200, -200., 200., 60, 0, 6);
  TH2F *hResVyVsY = new TH2F("hResVyVsY", "Res Vy vs Y", 200, -200., 200., 60, 0, 6);
  TH2F *hResVzVsY = new TH2F("hResVzVsY", "Res Vz vs Y", 200, -2000., 2000., 60, 0, 6);
	TH2F *hResVxVsLz = new TH2F("hResVxVsLz", "Res Vx vs Lz", 200, -200., 200., 200, 0, 50);
  TH2F *hResVyVsLz = new TH2F("hResVyVsLz", "Res Vy vs Lz", 200, -200., 200., 200, 0, 50);
  TH2F *hResVzVsLz = new TH2F("hResVzVsLz", "Res Vz vs Lz", 200, -2000., 2000., 200, 0, 50);
	TH2F *hResPxVsYprot = new TH2F("hResPxVsYprot", "Res Px vs Y", 200, -0.2, 0.2, 60, 0, 6);
	TH2F *hResPyVsYprot = new TH2F("hResPyVsYprot", "Res Py vs Y", 200, -0.2, 0.2, 60, 0, 6);
	TH2F *hResPzVsYprot = new TH2F("hResPzVsYprot", "Res Pz vs Y", 200, -0.2, 0.2, 60, 0, 6);
	TH2F *hResPxVsYpion = new TH2F("hResPxVsYpion", "Res Px vs Y", 200, -0.08, 0.08, 60, 0, 6);
	TH2F *hResPyVsYpion = new TH2F("hResPyVsYpion", "Res Py vs Y", 200, -0.08, 0.08, 60, 0, 6);
	TH2F *hResPzVsYpion = new TH2F("hResPzVsYpion", "Res Pz vs Y", 200, -0.08, 0.08, 60, 0, 6);
	TH2F *hResPxVsYpionRel = new TH2F("hResPxVsYpionRel", "Res Px vs Y Relativa", 200, -0.08, 0.08, 60, 0, 6);
	TH2F *hResPyVsYpionRel = new TH2F("hResPyVsYpionRel", "Res Py vs Y Relativa", 200, -0.08, 0.08, 60, 0, 6);
	TH2F *hResPzVsYpionRel = new TH2F("hResPzVsYpionRel", "Res Pz vs Y Relativa", 200, -0.08, 0.08, 60, 0, 6);
	TH2F *hResPxVsLzprot = new TH2F("hResPxVsLzprot", "Res Px vs Lz", 200, -1., 1., 200, 0, 50);
	TH2F *hResPxVsLzpion = new TH2F("hResPxVsLzpion", "Res Px vs Lz", 200, -1., 1., 200, 0, 50);
	TH2F *hResPyVsLzprot = new TH2F("hResPyVsLzprot", "Res Py vs Lz", 200, -1., 1., 200, 0, 50);
	TH2F *hResPyVsLzpion = new TH2F("hResPyVsLzpion", "Res Py vs Lz", 200, -1., 1., 200, 0, 50);
	TH2F *hResPzVsLzprot = new TH2F("hResPzVsLzprot", "Res Pz vs Lz", 200, -1., 1., 200, 0, 50);
	TH2F *hResPzVsLzpion = new TH2F("hResPzVsLzpion", "Res Pz vs Lz", 200, -1., 1., 200, 0, 50);
//	TH1F hmassnfake("hMassinrecnfakeprot", "Massinv for N fake prot", 3000, 0,3);
	TH1F* hMassinrecnfakeprot0 = new TH1F("hMassinrecnfakeprot", "Massinv for 0 fake prot", 3000, 0,3);
	TH1F* hMassinrecnfakeprot1 = new TH1F("hMassinrecnfakeprot1", "Massinv for 1 fake prot", 3000, 0,3);
	TH1F* hMassinrecnfakeprot2 = new TH1F("hMassinrecnfakeprot2", "Massinv for 2 fake prot", 3000, 0,3);
	TH1F* hMassinrecnfakeprot3 = new TH1F("hMassinrecnfakeprot3", "Massinv for 3 fake prot", 3000, 0,3);
	TH1F* hMassinrecnfakeprot4 = new TH1F("hMassinrecnfakeprot4", "Massinv for 4 fake prot", 3000, 0,3);
	TH1F* hMassinrecnfakeprot5= new TH1F("hMassinrecnfakeprot5", "Massinv for 5 fake prot", 3000, 0,3);
	TH1F* hDCAnfakeprot0 = new TH1F("hDCAnfakeprot0", "DCA for 0 fake prot", 100, 0,0.3);
	TH1F* hDCAnfakeprot4 = new TH1F("hDCAnfakeprot4", "DCA for 4 fake prot", 100, 0,0.3);
	for(int i=0 ; i<events; i++) {
		variables->GetEvent(i);
		double residVx=10000.*(xrec - secvertgenProt[0]);
    double residVy=10000.*(yrec - secvertgenProt[1]);
    double residVz=10000.*(zrec - secvertgenProt[2]);
		if(zrec<=20.){
			hMassinreczP->Fill(massinvrec, zrec);
			hNfakeVszP->Fill(nfaketrkprot, zrec);
			}
		hResVxVsY->Fill(residVx, y);
    hResVyVsY->Fill(residVy, y);
    hResVzVsY->Fill(residVz, y);
		hResVxVsLz->Fill(residVx, zrec);
    hResVyVsLz->Fill(residVy, zrec);
    hResVzVsLz->Fill(residVz, zrec);
		hMassinrecY->Fill(massinvrec, y);
		hResPxVsLzprot->Fill(pRecProt[0]-pGenProt[0],zrec);
		hResPxVsYprot->Fill(pRecProt[0]-pGenProt[0],y);
		hResPxVsLzpion->Fill(pRecPion[0]-pGenPion[0],zrec);
		hResPxVsYpion->Fill(pRecPion[0]-pGenPion[0],y);
		hResPyVsLzprot->Fill(pRecProt[1]-pGenProt[1],zrec);
		hResPyVsYprot->Fill(pRecProt[1]-pGenProt[1],y);
		hResPyVsLzpion->Fill(pRecPion[1]-pGenPion[1],zrec);
		hResPyVsYpion->Fill(pRecPion[1]-pGenPion[1],y);
		hResPzVsLzprot->Fill(pRecProt[2]-pGenProt[2],zrec);
		hResPzVsYprot->Fill(pRecProt[2]-pGenProt[2],y);
		hResPzVsLzpion->Fill(pRecPion[2]-pGenPion[2],zrec);
		hResPzVsYpion->Fill(pRecPion[2]-pGenPion[2],y);

		hResPxVsYpionRel->Fill((pRecPion[0]-pGenPion[0])/pRecPion[0],y);
		hResPyVsYpionRel->Fill((pRecPion[1]-pGenPion[1])/pRecPion[1],y);
		hResPzVsYpionRel->Fill((pRecPion[2]-pGenPion[2])/pRecPion[2],y);

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
	hResVxVsY->Write();
  hResVyVsY->Write();
  hResVzVsY->Write();
	hResVxVsLz->Write();
  hResVyVsLz->Write();
  hResVzVsLz->Write();
	hResPxVsLzprot->Write();
	hResPxVsYprot->Write();
	hResPxVsLzpion->Write();
	hResPxVsYpion->Write();
	hResPyVsLzprot->Write();
	hResPyVsYprot->Write();
	hResPyVsLzpion->Write();
	hResPyVsYpion->Write();
	hResPzVsLzprot->Write();
	hResPzVsYprot->Write();
	hResPzVsLzpion->Write();
	hResPzVsYpion->Write();
	hMassinreczP->Write();
	hResPxVsYpionRel->Write();
	hResPyVsYpionRel->Write();
	hResPzVsYpionRel->Write();

	TCanvas* c01=new TCanvas("c01");
	ProjectionforBinMassInv(hMassinreczP, "#Lambda Inv Mass vs Lz", "zP", c01, foutLambda);
	
	TCanvas* c1=new TCanvas("c1");
	ProjectionforBinMassInv(hMassinrecY, "#Lambda Inv Mass vs Y", "Y", c1, foutLambda);

	TCanvas* c2=new TCanvas("c2");
	ProjectionforBin(hResVxVsY, "#Lambda ResX vs Y", "Y", c2, foutLambda);
	TCanvas* c3=new TCanvas("c3");
	ProjectionforBin(hResVyVsY, "#Lambda ResY vs Y", "Y", c3, foutLambda);
	TCanvas* c4=new TCanvas("c4");
	ProjectionforBin(hResVzVsY, "#Lambda ResZ vs Y", "Y", c4, foutLambda);
	TCanvas* c5=new TCanvas("c5");
	ProjectionforBin(hResVxVsLz, "#Lambda ResX vs Lz", "Lz", c5, foutLambda);
	TCanvas* c6=new TCanvas("c6");
	ProjectionforBin(hResVyVsLz, "#Lambda ResY vs Lz", "Lz", c6, foutLambda);
	TCanvas* c7=new TCanvas("c7");
	ProjectionforBin(hResVzVsLz, "#Lambda ResZ vs Lz", "Lz", c7, foutLambda);

	TCanvas* c5pRel=new TCanvas("c5pRel");
	ProjectionforBin(hResPxVsYpionRel, "#Lambda+#bar{#Lambda} ResPx vs Y Rel", "Y", c5pRel, foutLambda);
	TCanvas* c6pRel=new TCanvas("c6pRel");
	ProjectionforBin(hResPyVsYpionRel, "#Lambda+#bar{#Lambda} ResPy vs Y Rel", "Y", c6pRel, foutLambda);
	TCanvas* c7pRel=new TCanvas("c7pRel");
	ProjectionforBin(hResPzVsYpionRel, "#Lambda+#bar{#Lambda} ResPz vs Y Rel", "Y", c7pRel, foutLambda);


	/*TH1F* gr11 = (TH1F*)foutLambda->Get("hRMSvsYhResVxVsY");
	TH1F* gr11g = (TH1F*)fLBar->Get("hRMSvsYhResVxVsY");
	TH1F *RMSlambda_lBarX = (TH1F*)gr11->Clone("RMSlambda_lBarX");
	RMSlambda_lBarX->Divide(gr11,gr11g);
	RMSlambda_lBarX->SetDirectory(0);
	TCanvas *eff=new TCanvas("eff");
	eff->Divide(3,1);
	eff->cd(1);
	RMSlambda_lBarX->SetTitle("RMS X_{rec}-X_{gen} #frac{#Lambda}{#Lambda Bar}");
	RMSlambda_lBarX->Draw();
	TH1F* gr12 = (TH1F*)foutLambda->Get("hRMSvsYhResVyVsY");
	TH1F* gr12g = (TH1F*)fLBar->Get("hRMSvsYhResVyVsY");
	TH1F *RMSlambda_lBarY = (TH1F*)gr11->Clone("RMSlambda_lBarY");
	RMSlambda_lBarY->Divide(gr12,gr12g);
	RMSlambda_lBarY->SetDirectory(0);
	eff->cd(2);
	RMSlambda_lBarY->SetTitle("RMS Y_{rec}-Y_{gen} #frac{#Lambda}{#Lambda Bar}");
	RMSlambda_lBarY->Draw();
	TH1F* gr13 = (TH1F*)foutLambda->Get("hRMSvsYhResVzVsY");
	TH1F* gr13g = (TH1F*)fLBar->Get("hRMSvsYhResVzVsY");
	TH1F *RMSlambda_lBarZ = (TH1F*)gr11->Clone("RMSlambda_lBarZ");
	RMSlambda_lBarZ->Divide(gr13,gr13g);
	RMSlambda_lBarZ->SetDirectory(0);
	eff->cd(3);
	RMSlambda_lBarZ->SetTitle("RMS Z_{rec}-Z_{gen} #frac{#Lambda}{#Lambda Bar}");
	RMSlambda_lBarZ->Draw();*/
	
	TH1F* gr11 = (TH1F*)foutLambda->Get("hsigmavsYhResVxVsY");
	TH1F* gr11g = (TH1F*)fLBar->Get("hsigmavsYhResVxVsY");
	TH1F *RMSlambda_lBarX = (TH1F*)gr11->Clone("RMSlambda_lBarX");
	RMSlambda_lBarX->Divide(gr11,gr11g);
	RMSlambda_lBarX->SetDirectory(0);
	TCanvas *eff=new TCanvas("eff");
	eff->Divide(3,1);
	eff->cd(1);
	RMSlambda_lBarX->SetTitle("#sigma X_{rec}-X_{gen} #frac{#Lambda}{#Lambda Bar}");
	RMSlambda_lBarX->Draw();
	TH1F* gr12 = (TH1F*)foutLambda->Get("hsigmavsYhResVyVsY");
	TH1F* gr12g = (TH1F*)fLBar->Get("hsigmavsYhResVyVsY");
	TH1F *RMSlambda_lBarY = (TH1F*)gr11->Clone("RMSlambda_lBarY");
	RMSlambda_lBarY->Divide(gr12,gr12g);
	RMSlambda_lBarY->SetDirectory(0);
	eff->cd(2);
	RMSlambda_lBarY->SetTitle("#sigma Y_{rec}-Y_{gen} #frac{#Lambda}{#Lambda Bar}");
	RMSlambda_lBarY->Draw();
	TH1F* gr13 = (TH1F*)foutLambda->Get("hsigmavsYhResVzVsY");
	TH1F* gr13g = (TH1F*)fLBar->Get("hsigmavsYhResVzVsY");
	TH1F *RMSlambda_lBarZ = (TH1F*)gr11->Clone("RMSlambda_lBarZ");
	RMSlambda_lBarZ->Divide(gr13,gr13g);
	RMSlambda_lBarZ->SetDirectory(0);
	eff->cd(3);
	RMSlambda_lBarZ->SetTitle("#sigma Z_{rec}-Z_{gen} #frac{#Lambda}{#Lambda Bar}");
	RMSlambda_lBarZ->Draw();

	//Altro metodo per fare il rapporto
	TGraphErrors* grLx = (TGraphErrors*)foutLambda->Get("SigmahResVxVsY");
	TGraphErrors* grLbx = (TGraphErrors*)fLBar->Get("SigmahResVxVsY");
	TGraphErrors* grLy = (TGraphErrors*)foutLambda->Get("SigmahResVyVsY");
	TGraphErrors* grLby = (TGraphErrors*)fLBar->Get("SigmahResVyVsY");
	TGraphErrors* grLz = (TGraphErrors*)foutLambda->Get("SigmahResVzVsY");
	TGraphErrors* grLbz = (TGraphErrors*)fLBar->Get("SigmahResVzVsY");
	printf("numpuntilambda %d e numpuntilambdabar %d \n",grLx->GetN(),grLbx->GetN());
	double ygr[40], sigmagr[40], esigmagr[40];
	int count=0;
	for(int i=0;i<grLx->GetN();i++){
			Double_t yL=grLx->GetPointX(i);
		for(int j=0;j<grLbx->GetN();j++){
			if(yL>(grLbx->GetPointX(j)-0.05) && yL<(grLbx->GetPointX(j)+0.05)){
				printf("ylambda=%f e ylambdabar=%f \n",grLx->GetPointX(i),grLbx->GetPointX(j));
				ygr[count]=grLbx->GetPointX(j);
				sigmagr[count]=(grLx->GetPointY(i)/grLbx->GetPointY(j));
				double a=grLx->GetPointY(i);
				double b=grLbx->GetPointY(j);
				double err_a=grLx->GetErrorY(i);
				double err_b=grLbx->GetErrorY(j);
				esigmagr[count]=sqrt((1/b)*(1/b)*err_a*err_a + (a/(b*b))*(a/(b*b))*err_b*err_b);
				count++;
					}
			}
		}
	TGraphErrors* grsigmaX = new TGraphErrors(count,ygr,sigmagr,0,esigmagr);
	grsigmaX->SetTitle("#sigma X_{rec}-X_{gen} #frac{#Lambda}{#Lambda Bar}");
TCanvas *cR=new TCanvas("cR");
	cR->Divide(3,1);
	cR->cd(1);
	grsigmaX->Draw();
printf("RESY numpuntilambda %d e numpuntilambdabar %d \n",grLy->GetN(),grLby->GetN());
	count=0;
	for(int i=0;i<grLy->GetN();i++){
			Double_t yL=grLy->GetPointX(i);
		for(int j=0;j<grLby->GetN();j++){
			if(yL>(grLby->GetPointX(j)-0.05) && yL<(grLby->GetPointX(j)+0.05)){
				printf("ylambda=%f e ylambdabar=%f \n",grLy->GetPointX(i),grLby->GetPointX(j));
				ygr[count]=grLby->GetPointX(j);
				sigmagr[count]=(grLy->GetPointY(i)/grLby->GetPointY(j));
				double a=grLy->GetPointY(i);
				double b=grLby->GetPointY(j);
				double err_a=grLy->GetErrorY(i);
				double err_b=grLby->GetErrorY(j);
				esigmagr[count]=sqrt((1/b)*(1/b)*err_a*err_a + (a/(b*b))*(a/(b*b))*err_b*err_b);
				count++;
					}
			}
		}
	TGraphErrors* grsigmaY = new TGraphErrors(count,ygr,sigmagr,0,esigmagr);
	grsigmaY->SetTitle("#sigma Y_{rec}-Y_{gen} #frac{#Lambda}{#Lambda Bar}");
cR->cd(2);
	grsigmaY->Draw();

printf("RESZ numpuntilambda %d e numpuntilambdabar %d \n",grLz->GetN(),grLbz->GetN());
	count=0;
	for(int i=0;i<grLz->GetN();i++){
			Double_t yL=grLz->GetPointX(i);
		for(int j=0;j<grLbz->GetN();j++){
			if(yL>(grLbz->GetPointX(j)-0.05) && yL<(grLbz->GetPointX(j)+0.05)){
				printf("ylambda=%f e ylambdabar=%f ",grLz->GetPointX(i),grLbz->GetPointX(j));
				ygr[count]=grLbz->GetPointX(j);
				
				sigmagr[count]=(grLz->GetPointY(i)/grLbz->GetPointY(j));
				double a=grLz->GetPointY(i);
				double b=grLbz->GetPointY(j);
				double err_a=grLz->GetErrorY(i);
				double err_b=grLbz->GetErrorY(j);
				esigmagr[count]=sqrt((1/b)*(1/b)*err_a*err_a + (a/(b*b))*(a/(b*b))*err_b*err_b);
				printf("sigma %f \n",sigmagr[count]);
				count++;
					}
			}
		}
	printf("count %d ",count);
	for(int i=0;i<count;i++){
		printf("y=%f sigma=%f sigmaerr=%f \n",ygr[i],sigmagr[i],esigmagr[i]);	
		}
	TGraphErrors* grsigmaZ = new TGraphErrors(count,ygr,sigmagr,0,esigmagr);
	grsigmaZ->SetTitle("#sigma Z_{rec}-Z_{gen} #frac{#Lambda}{#Lambda Bar}");
	cR->cd(3);
	grsigmaZ->Draw();


/*TCanvas* c8=new TCanvas("c8");
	c8->Divide(3,1);
	c8->cd(1);
	TGraphErrors* gr11 = (TGraphErrors*)foutLambda->Get("RMShResVxVsY");
	TGraphErrors* gr11g = (TGraphErrors*)foutLambda->Get("SigmahResVxVsY");
	gr11g->Draw("A*");
	gr11->Draw("P");
	TLegend* leg1 = new TLegend(0.11,0.7,0.3,0.88);
	leg1->SetHeader("Legend");
	leg1->AddEntry(gr11g,"#sigma","p");
	leg1->AddEntry(gr11,"RMS","p");
	leg1->Draw();
	c8->cd(2);
	TGraphErrors* gr12 = (TGraphErrors*)foutLambda->Get("RMShResVyVsY");
	TGraphErrors* gr12g = (TGraphErrors*)foutLambda->Get("SigmahResVyVsY");
	gr12g->Draw("A*");
	gr12->Draw("P");
	TLegend* leg2 = new TLegend(0.11,0.7,0.3,0.88);
	leg2->SetHeader("Legend");
	leg2->AddEntry(gr12g,"#sigma","p");
	leg2->AddEntry(gr12,"RMS","p");
	leg2->Draw();
	c8->cd(3);
	TGraphErrors* gr13 = (TGraphErrors*)foutLambda->Get("RMShResVzVsY");
	TGraphErrors* gr13g = (TGraphErrors*)foutLambda->Get("SigmahResVzVsY");
	gr13g->Draw("A*");
	gr13->Draw("P");
	TLegend* leg3 = new TLegend(0.11,0.7,0.3,0.88);
	leg3->SetHeader("Legend");
	leg3->AddEntry(gr13g,"#sigma","p");
	leg3->AddEntry(gr13,"RMS","p");
	leg3->Draw();
*/

	hLambda->Close();
	hD0->Close();
	//hfile.Close();
	filin.Close();
	foutLambda->Close();

}


void Projection(TH2F* h2D,  int range, bool rapidity, char *name ){
		std::vector<double> mean, RMS, meanErr, RMSErr, Lz;
		std::vector<TH1D*> hProject;
		double jrange=0.;
		if(rapidity==1) jrange=1.5;
		for(int i=0;i<range;i++){
				hProject.push_back(h2D->ProjectionY(Form("%s range [%f %f]",name,jrange,jrange+0.5), h2D->GetXaxis()->FindBin(jrange+0.01),h2D->GetXaxis()->FindBin(jrange+0.49)));
				mean.push_back(hProject.at(i)->GetMean());
				RMS.push_back(hProject.at(i)->GetRMS());
				meanErr.push_back(hProject.at(i)->GetMeanError());
				RMSErr.push_back(hProject.at(i)->GetRMSError());
				Lz.push_back((jrange+jrange+0.5)/2.);
				jrange=jrange+0.5;
				/*new TCanvas();
				hProject.at(i)->SetDirectory(0);
				hProject.at(i)->SetMinimum(0);
				hProject.at(i)->Draw();*/
			}
		/*if(rapidity==1){
			printf("range  %d",range);
			TGraphErrors* gr = new TGraphErrors(6,Lz,mean,0.,meanErr);
			gr->SetTitle("Mean  vs ?? ");
			new TCanvas();
			gr->Draw("AC*");
			}*/
	}
				
	
		/*for(int i=0;i<range;i++){
				hProject[i]=h2D->ProjectionY(Form("%s range [%f %f]",name,jrange,jrange+0.5), h2D->GetXaxis()->FindBin(jrange+0.01),h2D->GetXaxis()->FindBin(jrange+0.49));
				mean[i]=hProject[i]->GetMean();
				RMS[i]=hProject[i]->GetRMS();
				meanErr[i]=hProject[i]->GetMeanError();
				RMSErr[i]=hProject[i]->GetRMSError();
				Lz[i]=(jrange+jrange+0.5)/2.;
				jrange=jrange+0.5;
				new TCanvas();
				hProject[i]->SetDirectory(0);
				hProject[i]->SetMinimum(0);
				hProject[i]->Draw();
			}
		if(graph==1){
			TGraphErrors* gr = new TGraphErrors(kRange,Lz,mean,0,meanErr);
			gr->SetTitle("Mean  vs ?? ");
			new TCanvas();
			gr->Draw("AC*");
			}
}*/


void ProjectionResol(TH2F* hResol, bool rapid, char *title){
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
		/*hProjResyLz[i]->SetDirectory(0);
		new TCanvas();
		hProjResyLz[i]->Draw();*/
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
  c0->Divide(6,5);
	TH1F* hRMSvsY = new TH1F(Form("hRMSvsY%s",hist2D->GetName()), "RMS", 200, 0,6);
	TH1F* hsigmavsY = new TH1F(Form("hsigmavsY%s",hist2D->GetName()), "SIGMA", 200, 0,6);
  for(Int_t jj=1; jj<hist2D->GetNbinsY(); jj++){
 		
    TH1D *hInvMass = (TH1D *)hist2D->ProjectionX(Form("%sBin%d",hist2D->GetName(),jj), jj, jj);
		//hInvMass->SetDirectory(0);
    //hInvMass->Draw();
		//i=i+2;
    if(hInvMass->GetSumOfWeights()>45){
			count++;
			c0->cd(count);
			c0->SetTitle(Form("#Lambda %s",hist2D->GetName()));
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
			hRMSvsY->Fill(hist2D->GetYaxis()->GetBinUpEdge(jj),rmsh);
			hRMSvsY->SetBinError(hRMSvsY->FindBin(hist2D->GetYaxis()->GetBinUpEdge(jj)),Errrmsh);
			hsigmavsY->Fill(hist2D->GetYaxis()->GetBinUpEdge(jj),sigmag);
			hsigmavsY->SetBinError(hsigmavsY->FindBin(hist2D->GetYaxis()->GetBinUpEdge(jj)),esigma);
			}
    }
  }
	//printf("count=%d, countGraph=%d", count, countGraph); 
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
	gr1gaus->SetName(Form("Sigma%s ",hist2D->GetName()));
	gr1gaus->Write();
	hRMSvsY->Write();
	hsigmavsY->Write();
}

void ProjectionforBinMassInv(TH2F* hist2D, char* graphTitle , char* Xaxisname, TCanvas *c0, TFile *f){
	double xaxis[100], mean[100], errmean[100], sigma[100], errsigma[100], sigmagaus[100], errsigmagaus[100];
	int count=0;
  c0->Divide(4,7);
  for(Int_t jj=1; jj<hist2D->GetNbinsY(); jj++){
 		
    TH1D *hInvMass = (TH1D *)hist2D->ProjectionX(Form("%sBin%d",hist2D->GetName(),jj), jj, jj);
		//hInvMass->SetDirectory(0);
    //hInvMass->Draw();
    if(hInvMass->GetSumOfWeights()>45){
			count++;
			c0->cd(count);
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
	//printf("count=%d",count); 
	TGraphErrors *gr2=new TGraphErrors(count,xaxis,mean,0,errmean);
		gr2->SetTitle(Form("Mean  %s ",graphTitle));
		gr2->GetXaxis()->SetTitle(Form(" %s ",Xaxisname));
		gr2->SetMarkerColor(2);
		gr2->GetYaxis()->SetTitle("Mean");
		gr2->SetMarkerStyle(20);
		gr2->SetMarkerSize(0.7);
		//new TCanvas(Form("Graph vs  %s ",graphTitle),"Mean vs L_{z}",200,10,800,500);
		//gr2->Draw();
	TGraphErrors *gr1=new TGraphErrors(count,xaxis,sigma,0,errsigma);
		gr1->SetTitle(Form("RMS %s ",graphTitle));
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
		//new TCanvas(Form("Graph1 vs  %s ",graphTitle),"RMS vs L_{z}",200,10,800,500);
		//gr1->Draw();
	f->cd();
	gr2->SetName(Form("Mean%s",hist2D->GetName()));
	gr2->Write();
	gr1->SetName(Form("RMS%s",hist2D->GetName()));
	gr1->Write();
	gr1gaus->SetName(Form("Sigma%s ",hist2D->GetName()));
	gr1gaus->Write();
}
