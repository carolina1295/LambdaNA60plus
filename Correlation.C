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
//#include <basic_string>

#define EvBkg 5000
#define EvSig 44.9/1000000

void Correlation(){

    TFile *fCorrel= new TFile("VariablesCorrelations.root","RECREATE");
TFile filin("fntSig1000.root");
float mass, pt, yrap, dist, cosp, d01, d02, d0prod, dca, ptMin, ptMax, Alpha, qt;
TNtuple *variables = (TNtuple*)filin.Get("ntD0cand");
    variables->SetBranchAddress("mass",&mass);
    variables->SetBranchAddress("pt",&pt);
    variables->SetBranchAddress("y",&yrap);
    variables->SetBranchAddress("dist",&dist);
    variables->SetBranchAddress("cosp",&cosp);
    variables->SetBranchAddress("d01",&d01);
    variables->SetBranchAddress("d02",&d02);
    variables->SetBranchAddress("d0prod",&d0prod);
    variables->SetBranchAddress("dca",&dca);
    variables->SetBranchAddress("ptMin",&ptMin);
    variables->SetBranchAddress("ptMax",&ptMax);
    variables->SetBranchAddress("Alpha",&Alpha);
    variables->SetBranchAddress("qt",&qt);

    const char* named02[6]={ "d0prod (cm^{2})", "dist (cm)", "cosp ", "d01 (cm)", "dca (cm)", "ptMax (GeV/c)"};
    const char* named0prod[6]={ "d02 (cm)", "dist (cm)", "cosp ", "d01 (cm)", "dca (cm)", "ptMax (GeV/c)"};
    const char* namedist[6]={ "d02 (cm)", "d0prod (cm^{2})", "cosp ", "d01 (cm)", "dca (cm)", "ptMax (GeV/c)"};

    TH2F* hd02Corr[6], *hd0prodCorr[6], *hdistCorr[6];
    hd02Corr[0] = new TH2F("hd02_0","d0prod vs d02",200,-0.5,0.5,200,-0.05,0.05);
    hd02Corr[1] = new TH2F("hd02_1","dist vs d02",200,-0.5,0.5,200,0,20);
    hd02Corr[2] = new TH2F("hd02_2","cosp vs d02",200,-0.5,0.5,200,-1.2,1.2);
    hd02Corr[3] = new TH2F("hd02_3","d01 vs d02",200,-0.5,0.5,200,-0.5,0.5);
    hd02Corr[4] = new TH2F("hd02_4","dca vs d02",200,-0.5,0.5,200,0,4);
    hd02Corr[5] = new TH2F("hd02_5","ptMax vs d02",200,-0.5,0.5,200,0,20);

   // TH2F* hd02 = new TH2F("hd02_0","d0prod vs d02",200,-0.5,0.5,200,-0.05,0.05);

    hd0prodCorr[0] = new TH2F("hd0prod_0","d02 vs d0prod",200,-0.05,0.05,200,-0.5,0.5);
    hd0prodCorr[1] = new TH2F("hd0prod_1","dist vs d0prod",200,-0.05,0.05,200,0,20);
    hd0prodCorr[2] = new TH2F("hd0prod_2","cosp vs d0prod",200,-0.05,0.05,200,-1.2,1.2);
    hd0prodCorr[3] = new TH2F("hd0prod_3","d01 vs d0prod",200,-0.05,0.05,200,-0.5,0.5);
    hd0prodCorr[4] = new TH2F("hd0prod_4","dca vs d0prod",200,-0.05,0.05,200,0,4);
    hd0prodCorr[5] = new TH2F("hd0prod_5","ptMax vs d0prod",200,-0.05,0.05,200,0,20);

    hdistCorr[0] = new TH2F("hdist_0","d02 vs dist",200,0,20,200,-0.5,0.5);
    hdistCorr[1] = new TH2F("hdist_1","d0prod vs dist",200,0,20,200,-0.05,0.05);
    hdistCorr[2] = new TH2F("hdist_2","cosp vs dist",200,0,20,200,-1.2,1.2);
    hdistCorr[3] = new TH2F("hdist_3","d01 vs dist",200,0,20,200,0,20);
    hdistCorr[4] = new TH2F("hdist_4","dca vs dist",200,0,20,200,0,4);
    hdistCorr[5] = new TH2F("hdist_5","ptMax vs dist",200,0,20,200,0,20);

TH1F* hmassSig=new TH1F("hmassSig","Signal",200,0,1.5);
    int eventsSig=variables->GetEntries();

	for(int i=0 ; i<eventsSig; i++) {
        variables->GetEvent(i);

        hmassSig->Fill(mass);

        hd02Corr[0]->Fill(d02,d0prod);
        hd02Corr[1]->Fill(d02,dist);
        hd02Corr[2]->Fill(d02,cosp);
        hd02Corr[3]->Fill(d02,d01);
        hd02Corr[4]->Fill(d02,dca);
        hd02Corr[5]->Fill(d02,ptMax);

        hd0prodCorr[0]->Fill(d0prod,d02);
        hd0prodCorr[1]->Fill(d0prod,dist);
        hd0prodCorr[2]->Fill(d0prod,cosp);
        hd0prodCorr[3]->Fill(d0prod,d01);
        hd0prodCorr[4]->Fill(d0prod,dca);
        hd0prodCorr[5]->Fill(d0prod,ptMax);

        hdistCorr[0]->Fill(dist,d02);
        hdistCorr[1]->Fill(dist,d0prod);
        hdistCorr[2]->Fill(dist,cosp);
        hdistCorr[3]->Fill(dist,d01);
        hdistCorr[4]->Fill(dist,dca);
        hdistCorr[5]->Fill(dist,ptMax);
    }
    for(int iCorr=0;iCorr<6;iCorr++){
        hd02Corr[iCorr]->GetXaxis()->SetTitle("d02 (cm)");
        hd02Corr[iCorr]->GetYaxis()->SetTitle(Form("%s",named02[iCorr]));
        hd0prodCorr[iCorr]->GetXaxis()->SetTitle("d0prod (cm^{2})");
        hd0prodCorr[iCorr]->GetYaxis()->SetTitle(Form("%s",named0prod[iCorr]));
        hdistCorr[iCorr]->GetXaxis()->SetTitle("dist (cm)");
        hdistCorr[iCorr]->GetYaxis()->SetTitle(Form("%s",namedist[iCorr]));
    }
    fCorrel->cd();
    for(int iCorr=0;iCorr<6;iCorr++){
        hd02Corr[iCorr]->Write();
        hd0prodCorr[iCorr]->Write();
        hdistCorr[iCorr]->Write();
    }
    hmassSig->SetDirectory(0);
    hmassSig->Scale(EvSig);
    new TCanvas();
    hmassSig->Draw();

    filin.Close();

TFile filinBkg("fntBkg1000.root");
  	float massBkg, ptBkg, yrapBkg, distBkg, cospBkg, d01Bkg, d02Bkg, d0prodBkg, dcaBkg, ptMinBkg, ptMaxBkg, AlphaBkg, qtBkg;
  	TNtuple *variablesBkg = (TNtuple*)filinBkg.Get("ntD0cand");
		variablesBkg->SetBranchAddress("mass",&massBkg);
    	variablesBkg->SetBranchAddress("pt",&ptBkg);
    	variablesBkg->SetBranchAddress("y",&yrapBkg);
    	variablesBkg->SetBranchAddress("dist",&distBkg);
		variablesBkg->SetBranchAddress("cosp",&cospBkg);
		variablesBkg->SetBranchAddress("d01",&d01Bkg);
		variablesBkg->SetBranchAddress("d02",&d02Bkg);
		variablesBkg->SetBranchAddress("d0prod",&d0prodBkg);
		variablesBkg->SetBranchAddress("dca",&dcaBkg);
		variablesBkg->SetBranchAddress("ptMin",&ptMinBkg);
		variablesBkg->SetBranchAddress("ptMax",&ptMaxBkg);
        variablesBkg->SetBranchAddress("Alpha",&AlphaBkg);
		variablesBkg->SetBranchAddress("qt",&qtBkg);

    TH2F* hd02CorrBkg[6], *hd0prodCorrBkg[6], *hdistCorrBkg[6];
    hd02CorrBkg[0] = new TH2F("hd02Bkg_0","Bkg d0prod vs d02",200,-0.5,0.5,200,-0.05,0.05);
    hd02CorrBkg[1] = new TH2F("hd02Bkg_1","Bkg dist vs d02",200,-0.5,0.5,200,0,20);
    hd02CorrBkg[2] = new TH2F("hd02Bkg_2","Bkg cosp vs d02",200,-0.5,0.5,200,-1.2,1.2);
    hd02CorrBkg[3] = new TH2F("hd02Bkg_3","Bkg d01 vs d02",200,-0.5,0.5,200,-0.5,0.5);
    hd02CorrBkg[4] = new TH2F("hd02Bkg_4","Bkg dca vs d02",200,-0.5,0.5,200,0,4);
    hd02CorrBkg[5] = new TH2F("hd02_5","Bkg ptMax vs d02",200,-0.5,0.5,200,0,20);

    hd0prodCorrBkg[0] = new TH2F("hd0prodBkg_0","Bkg d02 vs d0prod",200,-0.05,0.05,200,-0.5,0.5);
    hd0prodCorrBkg[1] = new TH2F("hd0prodBkg_1","Bkg dist vs d0prod",200,-0.05,0.05,200,0,20);
    hd0prodCorrBkg[2] = new TH2F("hd0prodBkg_2","Bkg cosp vs d0prod",200,-0.05,0.05,200,-1.2,1.2);
    hd0prodCorrBkg[3] = new TH2F("hd0prodBkg_3","Bkg d01 vs d0prod",200,-0.05,0.05,200,-0.5,0.5);
    hd0prodCorrBkg[4] = new TH2F("hd0prodBkg_4","Bkg dca vs d0prod",200,-0.05,0.05,200,0,4);
    hd0prodCorrBkg[5] = new TH2F("hd0prodBkg_5","Bkg ptMax vs d0prod",200,-0.05,0.05,200,0,20);

    hdistCorrBkg[0] = new TH2F("hdistBkg_0","Bkg d02 vs dist",200,0,20,200,-0.5,0.5);
    hdistCorrBkg[1] = new TH2F("hdistBkg_1","Bkg d0prod vs dist",200,0,20,200,-0.05,0.05);
    hdistCorrBkg[2] = new TH2F("hdistBkg_2","Bkg cosp vs dist",200,0,20,200,-1.2,1.2);
    hdistCorrBkg[3] = new TH2F("hdistBkg_3","Bkg d01 vs dist",200,0,20,200,0,20);
    hdistCorrBkg[4] = new TH2F("hdistBkg_4","Bkg dca vs dist",200,0,20,200,0,4);
    hdistCorrBkg[5] = new TH2F("hdistBkg_5","Bkg ptMax vs dist",200,0,20,200,0,20);
TH1F* hmassBkg=new TH1F("hmassBkg","Background",20000,1,1.5);
    

    int eventsBkg=variablesBkg->GetEntries();
	for(int i=0 ; i<eventsBkg; i++) {
        variablesBkg->GetEvent(i);

        hmassBkg->Fill(massBkg);

        hd02CorrBkg[0]->Fill(d02Bkg,d0prodBkg);
        hd02CorrBkg[1]->Fill(d02Bkg,distBkg);
        hd02CorrBkg[2]->Fill(d02Bkg,d01Bkg);
        hd02CorrBkg[3]->Fill(d02Bkg,cospBkg);
        hd02CorrBkg[4]->Fill(d02Bkg,dcaBkg);
        hd02CorrBkg[5]->Fill(d02Bkg,ptMaxBkg);

        hd0prodCorrBkg[0]->Fill(d0prodBkg,d02Bkg);
        hd0prodCorrBkg[1]->Fill(d0prodBkg,distBkg);
        hd0prodCorrBkg[2]->Fill(d0prodBkg,d01Bkg);
        hd0prodCorrBkg[3]->Fill(d0prodBkg,cospBkg);
        hd0prodCorrBkg[4]->Fill(d0prodBkg,dcaBkg);
        hd0prodCorrBkg[5]->Fill(d0prodBkg,ptMaxBkg);

        hdistCorrBkg[0]->Fill(distBkg,d02Bkg);
        hdistCorrBkg[1]->Fill(distBkg,d0prodBkg);
        hdistCorrBkg[2]->Fill(distBkg,d01Bkg);
        hdistCorrBkg[3]->Fill(distBkg,cospBkg);
        hdistCorrBkg[4]->Fill(distBkg,dcaBkg);
        hdistCorrBkg[5]->Fill(distBkg,ptMaxBkg);
    }

    hmassBkg->SetDirectory(0);
   hmassBkg->Scale(1/EvBkg);
          new TCanvas();
    hmassBkg->Draw();

    TCanvas *c=new TCanvas("c");
            
    hmassBkg->Draw();
    hmassSig->Draw("SAME");


    for(int iCorrBkg=0;iCorrBkg<6;iCorrBkg++){
        hd02CorrBkg[iCorrBkg]->GetXaxis()->SetTitle("d02 (cm)");
        hd02CorrBkg[iCorrBkg]->GetYaxis()->SetTitle(Form("%s",named02[iCorrBkg]));
        hd0prodCorrBkg[iCorrBkg]->GetXaxis()->SetTitle("d0prod (cm)");
        hd0prodCorrBkg[iCorrBkg]->GetYaxis()->SetTitle(Form("%s",named02[iCorrBkg]));
        hdistCorrBkg[iCorrBkg]->GetXaxis()->SetTitle("dist (cm)");
        hdistCorrBkg[iCorrBkg]->GetYaxis()->SetTitle(Form("%s",named02[iCorrBkg]));
    }
    fCorrel->cd();

    for(int iCorr=0;iCorr<6;iCorr++){
        hd02CorrBkg[iCorr]->Write();
        hd0prodCorrBkg[iCorr]->Write();
        hdistCorrBkg[iCorr]->Write();
    }





    fCorrel->Close();

filinBkg.Close();

}





