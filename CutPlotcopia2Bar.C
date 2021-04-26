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
#include "TGraphErrors.h"
#include <string>
#include <vector>
#include <TNtuple.h>
#include "TH3F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TMultiGraph.h"
#include <iostream>
#include <fstream>
#include "TF1.h"
#include <time.h>
//#include <basic_string>
#define DEBUG 0
#define EvBkg 1000
#define EvSig 3.07/1000000
#define BR 0.64
#define NVar 7
#define NStep 21

//void DrawSameFit(TH1F* hBkg, TH1F* hSig, TCanvas *c, const char* sTitle, double val[2]);
void DrawFit( TH1F* hSig, TCanvas *c, const char* sTitle, double val[2]);
void DrawSame(TH1F* hBkg, TH1F* hSig, TCanvas *c, const char* sTitle);
void ReadPreCut(string nameVar[NVar], double cutVar[NVar], bool set[NVar], char cond[NVar]);
void ReadScanCut(string nameVarScan[NVar], double Stepcut[NVar], double min[NVar], double max[NVar], bool setScan[NVar], char condScan[NVar]);
bool AreEqual(bool SetCut[NVar], bool CheckCut[NVar]);
void Plot(double CutSig[NVar][NStep], double CutBkg[NVar][NStep],double Bkg[NVar][NStep], double TotSig, double TotBkg,TFile *f);

void CutPlot(bool ScanX2){
    
    clock_t start=clock();

    int bin=0;
    double norm=1;
    double FitVar[2]={0};

    //PreCut variables
    string CutVariable[NVar];
    double PreCut[NVar];
    bool On_Off[NVar];
    char high_low[NVar];

    //Scan Variables
    string ScanVariable[NVar];
    double StepScan[NVar], ScanMin[NVar], ScanMax[NVar];
    bool On_OffScan[NVar];
    char high_lowScan[NVar];

    TFile *fout= new TFile("ScanX.root","RECREATE");
    ReadPreCut(CutVariable,PreCut,On_Off,high_low);
    ReadScanCut(ScanVariable, StepScan, ScanMin, ScanMax, On_OffScan, high_lowScan);

    for(int i=0;i<NVar;i++){
        cout<<CutVariable[i]<<"\t"<<PreCut[i]<<"\t"<<On_Off[i]<<high_low[i] <<endl;
    }
    for(int i=0;i<NVar;i++){
        cout<<ScanVariable[i]<<"\t"<<StepScan[i]<<"\t"<<"\t"<<ScanMin[i]<<"\t"<<ScanMax[i]<<"\t"<<On_OffScan[i]<<high_lowScan[i] <<endl;
    }

    //Define the Matrix with cut
    double MxCutScan[NVar][NStep];
    for(int iVar=0;iVar<NVar;iVar++){
        double Min=ScanMin[iVar];
        double Max=ScanMax[iVar];
        double Step=StepScan[iVar];
        //printf("Scan %d) ",iVar);
        for (int iStep = 0; iStep < NStep; iStep++){
            MxCutScan[iVar][iStep]=Min;
            //printf("%f \n ", MxCutScan[iVar][iStep]);
            Min=Min+Step;
        }
    }
    double S[NVar][NStep] = {{0},{0}};
    double B[NVar][NStep] = {{0},{0}};
    double Bsig[NVar][NStep] = {{0},{0}};
    //Get positions for the double scan 
    bool GetCut1=0, GetCut2=0;
    int pos1=0, pos2=0, check=0;
    for(int iVar=0;iVar<NVar;iVar++){
        if(On_OffScan[iVar]==1 && GetCut1==0){
            pos1=iVar;
            GetCut1=1;
            check++;
            continue;
        }
        if(On_OffScan[iVar]==1 && GetCut1==1 ){
            pos2=iVar;
            GetCut2=1;
            check++;
        }
        printf("%d) pos1 %d pos2 %d \n",iVar,pos1,pos2);
    }
    int count=0;

    TFile filin("fntSig1000Bar.root");
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

        TH2F* hArmenLambda = new TH2F("hArmenLambda","Armenteros Plot for #Lambda ",200,-1,1,200,0,1);
        TH1F* hmassSig = new TH1F("hmassSig","Invariant Mass for #Lambda",200,1,1.5);
        TH1F* hySig = new TH1F("hySig","Rapidity for #Lambda",200,0,6);
        TH1F* hptSig = new TH1F("hptSig","P_{t} for #Lambda",200,0,6);
        TH1F* hptMinSig = new TH1F("hptMinSig","P_{t} Min for #Lambda",200,0,20);
        TH1F* hptMaxSig = new TH1F("hptMaxSig","P_{t} Max for #Lambda",200,0,80);
        TH1F* hdistSig = new TH1F("hdistSig","dist for #Lambda",200,0,20);
        TH1F* hcospSig = new TH1F("hcospSig","cosp for #Lambda",200,-1.5,1.5);
        TH1F* hd01Sig = new TH1F("hd01Sig","d01 for #Lambda",200,-0.5,0.5);
        TH1F* hd02Sig = new TH1F("hd02Sig","d02 for #Lambda",200,-0.5,0.5);  
        TH1F* hd0prodSig = new TH1F("hd0prodSig","d0prod for #Lambda",200,-0.05,0.05); 
        TH1F* hdcaSig = new TH1F("hdcaSig","dca for #Lambda",200,0,1); 
        
        TH2F* hSigScan = new TH2F("hSigScan","Signal candidates survived",NStep,-0.5,NStep-0.5,NStep,-0.5,NStep-0.5);
                                
        TH1F* hmassSigCut = new TH1F("hmassSigCut","Invariant Mass for #Lambda post PreCut",200,1,1.5);
        
        #if DEBUG
            TH2F* hArmenLambdaCut = new TH2F("hArmenLambdaCut","Armenteros Plot for #Lambda post PreCut",200,-1,1,200,0,1);
            //TH1F* hmassSigCut = new TH1F("hmassSigCut","Invariant Mass for #Lambda post PreCut",200,1,1.5);
            TH1F* hySigCut = new TH1F("hySigCut","Rapidity for #Lambda post PreCut",200,0,6);
            TH1F* hptSigCut = new TH1F("hptSigCut","P_{t} for #Lambda post PreCut",200,0,6);
            TH1F* hptMinSigCut = new TH1F("hptMinSigCut","P_{t} Min for #Lambda post PreCut",200,0,20);
            TH1F* hptMaxSigCut = new TH1F("hptMaxSigCut","P_{t} Max for #Lambda post PreCut",200,0,80);
            TH1F* hdistSigCut = new TH1F("hdistSigCut","dist for #Lambda post PreCut",200,0,20);
            TH1F* hcospSigCut = new TH1F("hcospSigCut","cosp for #Lambda post PreCut",200,-1.5,1.5);
            TH1F* hd01SigCut = new TH1F("hd01SigCut","d01 for #Lambda post PreCut",200,-0.5,0.5);
            TH1F* hd02SigCut = new TH1F("hd02SigCut","d02 for #Lambda post PreCut",200,-0.5,0.5);  
            TH1F* hd0prodSigCut = new TH1F("hd0prodSigCut","d0prod for #Lambda post PreCut",200,-0.05,0.05); 
            TH1F* hdcaSigCut = new TH1F("hdcaSigCut","dca for #Lambda post PreCut",200,0,1);
        #endif
        
        double Scounts[NStep][NStep]={{0},{0}}; //definito per debuggare il doubleScan e riempire th2f significance e s/b
        int scanPos=0;

        int eventsSig=variables->GetEntries();

        //Riempio tutto
        for (int i = 0; i < eventsSig; i++){
            variables->GetEvent(i);
            if(Alpha<-0.4 && qt<0.2){
                hmassSig->Fill(mass);
                hySig->Fill(yrap);
                hptSig->Fill(pt);
                hptMinSig->Fill(ptMin);
                hptMaxSig->Fill(ptMax);
                hdistSig->Fill(dist);  
                hcospSig->Fill(cosp); 
                hd01Sig->Fill(d01);  
                hd02Sig->Fill(d02);  
                hd0prodSig->Fill(d0prod);  
                hdcaSig->Fill(dca);
                hArmenLambda->Fill(Alpha,qt);
            }
        }

        hmassSig->SetDirectory(0);
        hySig->SetDirectory(0);
        hptSig->SetDirectory(0);
        hptMinSig->SetDirectory(0);
        hptMaxSig->SetDirectory(0);
        hdistSig->SetDirectory(0);
        hd01Sig->SetDirectory(0);
        hd02Sig->SetDirectory(0);
        hd0prodSig->SetDirectory(0);
        hcospSig->SetDirectory(0);
        hdcaSig->SetDirectory(0);
        hArmenLambda->SetDirectory(0);

        //Normalization of signal
        hySig->Scale(norm/hySig->Integral());
        hptSig->Scale(norm/hptSig->Integral());
        hptMinSig->Scale(norm/hptMinSig->Integral());
        hptMaxSig->Scale(norm/hptMaxSig->Integral());
        hdistSig->Scale(norm/hdistSig->Integral());
        hd01Sig->Scale(norm/hd01Sig->Integral());
        hd02Sig->Scale(norm/hd02Sig->Integral());
        hcospSig->Scale(norm/hcospSig->Integral());
        hd0prodSig->Scale(norm/hd0prodSig->Integral());
        hdcaSig->Scale(norm/hdcaSig->Integral());
        hmassSig->Scale((norm/1000000)*44.9*BR);
        
        TCanvas *cSig=new TCanvas("cSig");
        hmassSig->SetMarkerStyle(20);
        hmassSig->Draw();

        double Stot = hmassSig->GetEntries();
        
        //Calcolo mean e sigma massa invariante
        TCanvas *cC= new TCanvas("cC");
        DrawFit(hmassSig, cC,"mass",FitVar);
        double meanLow = FitVar[0]-3*FitVar[1];
        double meanHigh = FitVar[0]+3*FitVar[1];
       
        double var[NVar]={0};
        int Sigcandid=0;
        //Apply cut
        for(int i=0 ; i<eventsSig; i++) {
            variables->GetEvent(i);
            bool CutOk[NVar]={0};
            if(Alpha<-0.4 && qt<0.2){
                var[0]=abs(d02);
                var[1]=d0prod;
                var[2]=dist;
                var[3]=cosp;
                var[4]=abs(d01);
                var[5]=dca;
                var[6]=ptMax;

                if(On_Off[0]==1 && var[0]<PreCut[0]) continue;
                if(On_Off[1]==1 && var[1]>PreCut[1]) continue;
                if(On_Off[2]==1 && var[2]<PreCut[2]) continue;
                if(On_Off[3]==1 && var[3]<PreCut[3]) continue;
                if(On_Off[4]==1 && var[4]<PreCut[4]) continue;
                if(On_Off[5]==1 && var[5]<PreCut[5]) continue;
                if(On_Off[6]==1 && var[6]<PreCut[6]) continue;
                
                Sigcandid++;
                //printf("ev %d) %f \n",i,var[0]);
                hmassSigCut->Fill(mass);
                #if DEBUG
                    hySigCut->Fill(yrap);
                    hptSigCut->Fill(pt);
                    hptMinSigCut->Fill(ptMin);
                    hptMaxSigCut->Fill(ptMax);
                    hdistSigCut->Fill(dist);  
                    hcospSigCut->Fill(cosp); 
                    hd01SigCut->Fill(d01);  
                    hd02SigCut->Fill(d02);  
                    hd0prodSigCut->Fill(d0prod);  
                    hdcaSigCut->Fill(dca);
                    hArmenLambdaCut->Fill(Alpha,qt);
                #endif

                // Single Scan 
                for(int iVar=0;iVar<NVar;iVar++){
                    for(int iScan=0;iScan<NStep;iScan++){
                        if (On_OffScan[iVar]==1 && high_lowScan[iVar]=='+'){
                        //printf("Scanmin %f ScanMax %f Step%f \n",ScanMin[iVar],ScanMax[iVar],StepScan[iVar]);
                            //printf("scan %f %f \n",MxCutScan[iVar][iScan],var[iVar]);
                            if(var[iVar]>MxCutScan[iVar][iScan] && mass>=meanLow && mass<=meanHigh){
                                S[iVar][iScan]=(S[iVar][iScan])+1;
                                //printf("S %f \n",S[iVar][iScan]);
                            }
                        }
                        if (On_OffScan[iVar]==1 && high_lowScan[iVar]=='-'){ 
                            //printf("scan %f %f \n",MxCutScan[iVar][iScan],var[iVar]);
                            if(var[iVar]<MxCutScan[iVar][iScan] && mass>=meanLow && mass<=meanHigh){   
                                S[iVar][iScan]=(S[iVar][iScan])+1;
                            }
                        } 
                    }  
                }

                //Double Scan
                if(ScanX2==1){
                    bin=0;
                    double cut1=0, cut2=0;
                    for(int iCut1=0;iCut1<NStep;iCut1++){
                        cut1=MxCutScan[pos1][iCut1];
                        for(int iCut2=0;iCut2<NStep;iCut2++){
                            cut2=MxCutScan[pos2][iCut2];
                            if(high_lowScan[pos1]=='+' && high_lowScan[pos2]=='+'){
                                if(var[pos1]>cut1 && var[pos2]>cut2 && mass>=meanLow && mass<=meanHigh){
                                    //printf("var1 %s var2 %s, cut1 %f cut2 %f \n",ScanVariable[pos1].c_str(),ScanVariable[pos2].c_str(),cut1,cut2);
                                    bin = hSigScan->GetBin(iCut1+1,iCut2+1,0);
                                    hSigScan->AddBinContent(bin);
                                    Scounts[iCut1][iCut2]++;
                                }
                            }
                            if(high_lowScan[pos1]=='+' && high_lowScan[pos2]=='-'){
                                if(var[pos1]>cut1 && var[pos2]<cut2 && mass>=meanLow && mass<=meanHigh){
                                    //printf("var1 %s var2 %s cut1 %f cut2 %f \n",ScanVariable[pos1].c_str(),ScanVariable[pos2].c_str(),cut1,cut2);
                                    bin = hSigScan->GetBin(iCut1+1,iCut2+1,0);
                                    hSigScan->AddBinContent(bin);
                                    Scounts[iCut1][iCut2]++;
                                }
                            }
                            if(high_lowScan[pos1]=='-' && high_lowScan[pos2]=='+'){
                                if(var[pos1]<cut1 && var[pos2]>cut2 && mass>=meanLow && mass<=meanHigh){
                                    //printf("var1 %s var2 %s \n",ScanVariable[pos1].c_str(),ScanVariable[pos2].c_str());
                                    bin = hSigScan->GetBin(iCut1+1,iCut2+1,0);
                                    hSigScan->AddBinContent(bin);
                                    Scounts[iCut1][iCut2]++;
                                }
                            }
                            if(high_lowScan[pos1]=='-' && high_lowScan[pos2]=='-'){
                                if(var[pos1]<cut1 && var[pos2]<cut2 && mass>=meanLow && mass<=meanHigh){
                                    //printf("var1 %s var2 %s \n",ScanVariable[pos1].c_str(),ScanVariable[pos2].c_str());
                                    bin = hSigScan->GetBin(iCut1+1,iCut2+1,0);
                                    hSigScan->AddBinContent(bin);
                                    Scounts[iCut1][iCut2]++;
                                }
                            }
                        }
                    }
                }
            }
        } 

        //printf("Stot %f, %d",Stot,c);

        for(int i=0;i<NStep;i++){
            for(int j=0;j<NStep;j++){
                //printf("Scounts %f, eff %f \n",Scounts[i][j],Scounts[i][j]/Stot);
                bin= hSigScan->GetBin(i+1,j+1,0);
                //printf("S %f \n",hSigScan->GetBinContent(bin));
                
            }
        }

        double Min1, Min2;
        if(ScanX2==1){
            //Candidate sopravvissute
            hSigScan->SetDirectory(0);
            hSigScan->GetXaxis()->SetTitle(Form("%s ",ScanVariable[pos1].c_str()));
            hSigScan->GetYaxis()->SetTitle(Form("%s ",ScanVariable[pos2].c_str()));
            hSigScan->SetEntries(1);
            /*new TCanvas();
            hSigScan->Draw("COLZ");*/

            //calcolo efficienza
            TH2F* hEffScan2Var = (TH2F*)hSigScan->Clone("hEffScan2Var");
            hEffScan2Var->SetDirectory(0);
            hEffScan2Var->Scale(1/Stot);
            hEffScan2Var->SetTitle("Efficiency");
            Min1=ScanMin[pos1];
            Min2=ScanMin[pos2];
            for(int iLab=0;iLab<NStep;iLab++){
                hSigScan->GetXaxis()->SetBinLabel(iLab+1,Form("%f",Min1));
                hSigScan->GetYaxis()->SetBinLabel(iLab+1,Form("%f",Min2));

                hEffScan2Var->GetXaxis()->SetBinLabel(iLab+1,Form("%f",Min1));
                hEffScan2Var->GetYaxis()->SetBinLabel(iLab+1,Form("%f",Min2));

                Min1=Min1+StepScan[pos1];
                Min2=Min2+StepScan[pos2];
            }
            TCanvas *ceff = new TCanvas("ceff");
            hEffScan2Var->Draw("COLZ");
        }
        printf("Stot %d, SXev %f",Sigcandid,Sigcandid*EvSig*BR);
        double SigXevent=(Sigcandid*EvSig*BR);
        printf("pippo %f \n",SigXevent);
        hmassSigCut->SetDirectory(0);
        hmassSigCut->Scale((norm/1000000)*44.9*BR);
        TCanvas *mSig=new TCanvas("mSig");
        //double inteBin=hmassSigCut->Integral();
        TF1 *f1 = new TF1("f1","gaus",1.08,1.2);
        hmassSigCut->Fit(f1,"","",1.08,1.2);
        hmassSigCut->Draw();
    #if DEBUG
        hySigCut->SetDirectory(0);
        hptSigCut->SetDirectory(0);
        hptMinSigCut->SetDirectory(0);
        hptMaxSigCut->SetDirectory(0);
        hdistSigCut->SetDirectory(0);
        hd01SigCut->SetDirectory(0);
        hd02SigCut->SetDirectory(0);
        hd0prodSigCut->SetDirectory(0);
        hcospSigCut->SetDirectory(0);
        hdcaSigCut->SetDirectory(0);
        hArmenLambdaCut->SetDirectory(0);

        hySigCut->Scale(norm/hySigCut->Integral());
        hptSigCut->Scale(norm/hptSigCut->Integral());
        hptMinSigCut->Scale(norm/hptMinSigCut->Integral());
        hptMaxSigCut->Scale(norm/hptMaxSigCut->Integral());
        hdistSigCut->Scale(norm/hdistSigCut->Integral());
        hd01SigCut->Scale(norm/hd01SigCut->Integral());
        hd02SigCut->Scale(norm/hd02SigCut->Integral());
        hcospSigCut->Scale(norm/hcospSigCut->Integral());
        hd0prodSigCut->Scale(norm/hd0prodSigCut->Integral());
        hdcaSigCut->Scale(norm/hdcaSigCut->Integral());

      /*  TCanvas *cy = new TCanvas("cy");
        DrawSame(hySigCut, hySig, cy,"y");
        TCanvas *cpt = new TCanvas("cpt");
        DrawSame(hptSigCut, hptSig, cpt, "pt");
        TCanvas *cptMin = new TCanvas("cptMin");
        DrawSame(hptMinSigCut, hptMinSig, cptMin, "ptMin");
        TCanvas *cptMax = new TCanvas("cptMax");
        DrawSame(hptMaxSigCut, hptMaxSig, cptMax,"ptMax");
        TCanvas *cdist = new TCanvas("cdist");
        DrawSame(hdistSigCut, hdistSig, cdist, "dist");
        TCanvas *ccosp = new TCanvas("ccosp");
        DrawSame(hcospSigCut, hcospSig, ccosp, "cosp");
        TCanvas *cd01 = new TCanvas("cd01");
        DrawSame(hd01SigCut, hd01Sig, cd01, "d01");
        TCanvas *cd02 = new TCanvas("cd02");
        DrawSame(hd02SigCut, hd02Sig, cd02, "d02");
        TCanvas *cd0prod = new TCanvas("cd0prod");
        DrawSame(hd0prodSigCut, hd0prodSig, cd0prod,"d0prod");
        TCanvas *cdca = new TCanvas("cdca");
        DrawSame(hdcaSigCut, hdcaSig, cdca,"dca");*/
    #endif

    filin.Close();

    TFile filinBkg("fntBkg5000Bar.root");
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
    

    TH2F* hArmenBkg= new TH2F("hArmenBkg","Armenteros Plot for Bkg",200,-1,1,200,0,1);
    TH1F* hmassBkg = new TH1F("hmassBkg","Mass for Bkg",200,1.,1.5);
	TH1F* hyBkg = new TH1F("hyBkg","Rapidity for Bkg",200,0,6);
	TH1F* hptBkg = new TH1F("hptBkg","P_{t} for Bkg",200,0,6);
    TH1F* hptMinBkg = new TH1F("hptMinBkg","P_{t} Min for Bkg",200,0,40);
    TH1F* hptMaxBkg = new TH1F("hptMaxBkg","P_{t} Max for Bkg",200,0,60);
    TH1F* hdistBkg = new TH1F("hdistBkg","dist for Bkg",200,0,20);
    TH1F* hcospBkg = new TH1F("hcospBkg","cosp for Bkg",200,-1.5,1.5);
    TH1F* hd01Bkg = new TH1F("hd01Bkg","d01 for Bkg",200,-0.1,0.1);
    TH1F* hd02Bkg = new TH1F("hd02Bkg","d02 for Bkg",200,-0.1,0.1);  
    TH1F* hd0prodBkg = new TH1F("hd0prodBkg","d0prod for Bkg",200,-0.01,0.01); 
    TH1F* hdcaBkg = new TH1F("hdcaBkg","dca for Bkg",200,0,1); 
    TH1F* hAlphaBkg= new TH1F("hAlphaBkg","Alpha for Bkg",200,-1,1);
    TH1F* hqtBkg= new TH1F("hqtBkg","qt for Bkg",200,0,1);

    TH2F* hBkgScan = new TH2F("hBkgScan","Bkg candidates survived",NStep,-0.5,NStep-0.5,NStep,-0.5,NStep-0.5);

    TH1F* hmassBkgCut = new TH1F("hmassBkgCut","Mass for Bkg post PreCut",200,0.75,1.5);

    double Bcounts[NStep][NStep]={{0},{0}};
    #if DEBUG
        TH2F* hArmenBkgCut= new TH2F("hArmenBkgCut","Armenteros Plot for Bkg post PreCut",200,-1,1,200,0,1);
        TH1F* hyBkgCut = new TH1F("hyBkgCut","Rapidity for Bkg post PreCut ",200,0,6);
        TH1F* hptBkgCut = new TH1F("hptBkgCut","P_{t} for Bkg post PreCut",200,0,6);
        TH1F* hptMinBkgCut = new TH1F("hptMinBkgCut","P_{t} Min for Bkg post PreCut",200,0,40);
        TH1F* hptMaxBkgCut = new TH1F("hptMaxBkgCut","P_{t} Max for Bkg post PreCut",200,0,60);
        TH1F* hdistBkgCut = new TH1F("hdistBkgCut","dist for Bkg post PreCut",200,0,20);
        TH1F* hcospBkgCut = new TH1F("hcospBkgCut","cosp for Bkg post PreCut",200,-1.5,1.5);
        TH1F* hd01BkgCut = new TH1F("hd01BkgCut","d01 for Bkg post PreCut",200,-0.1,0.1);
        TH1F* hd02BkgCut = new TH1F("hd02BkgCut","d02 for Bkg post PreCut",200,-0.1,0.1);  
        TH1F* hd0prodBkgCut = new TH1F("hd0prodBkgCut","d0prod for Bkg post PreCut",200,-0.01,0.01); 
        TH1F* hdcaBkgCut = new TH1F("hdcaBkgCut","dca for Bkg post PreCut",200,0,1); 
        TH1F* hAlphaBkgCut= new TH1F("hAlphaBkgCut","Alpha for Bkg post PreCut",200,-1,1);
        TH1F* hqtBkgCut = new TH1F("hqtBkgCut","qt for Bkg post PreCut",200,0,1);
    #endif
    double Bkgentries=0;
    int counts=0;
    int countsCut=0;
    int eventsBkg=variablesBkg->GetEntries();
	for(int i=0 ; i<eventsBkg; i++){
        variablesBkg->GetEvent(i);
        if(AlphaBkg<-0.4 && qtBkg<0.2){
            hArmenBkg->Fill(AlphaBkg,qtBkg);
            hmassBkg->Fill(massBkg);
            if(massBkg>=meanLow && massBkg<=meanHigh) Bkgentries++;
            hyBkg->Fill(yrapBkg);
            hptBkg->Fill(ptBkg);
            hptMinBkg->Fill(ptMinBkg);
            hptMaxBkg->Fill(ptMaxBkg);
            hdistBkg->Fill(distBkg);  
            hcospBkg->Fill(cospBkg); 
            hd01Bkg->Fill(d01Bkg);  
            hd02Bkg->Fill(d02Bkg);  
            hd0prodBkg->Fill(d0prodBkg);  
            hdcaBkg->Fill(dcaBkg);
            hAlphaBkg->Fill(AlphaBkg);
            hqtBkg->Fill(qtBkg);
        }
    }

    hmassBkg->SetDirectory(0);
    hyBkg->SetDirectory(0);
    hptBkg->SetDirectory(0);
    hptMinBkg->SetDirectory(0);
    hptMaxBkg->SetDirectory(0);
    hdistBkg->SetDirectory(0);
    hd01Bkg->SetDirectory(0);
    hd02Bkg->SetDirectory(0);
    hd0prodBkg->SetDirectory(0);
    hdcaBkg->SetDirectory(0);
    hcospBkg->SetDirectory(0);
    hArmenBkg->SetDirectory(0);
    hAlphaBkg->SetDirectory(0);
    hqtBkg->SetDirectory(0);

    double Btot=hmassBkg->GetEntries();
    hmassBkg->Scale(norm/EvBkg);
    TAxis *axis = hmassBkg->GetXaxis();
    int bmin = axis->FindBin(meanLow); 
    int bmax = axis->FindBin(meanHigh); 
    double BInt=hmassBkg->Integral(bmin,bmax);
   /* TCanvas *c1Bkg=new TCanvas("c1Bkg");
    hmassBkg->SetName("BkgPrecut");
    hmassBkg->SetMarkerStyle(22);
    hmassBkg->Draw();*/
   
    hyBkg->Scale(norm/hyBkg->Integral());
    hptBkg->Scale(norm/hptBkg->Integral());
    hptMinBkg->Scale(norm/hptMinBkg->Integral());
    hptMaxBkg->Scale(norm/hptMaxBkg->Integral());
    hdistBkg->Scale(norm/hdistBkg->Integral());
    hd01Bkg->Scale(norm/hd01Bkg->Integral());
    hd02Bkg->Scale(norm/hd02Bkg->Integral());
    hcospBkg->Scale(norm/hcospBkg->Integral());
    hd0prodBkg->Scale(norm/hd0prodBkg->Integral());
    hdcaBkg->Scale(norm/hdcaBkg->Integral());
    hAlphaBkg->Scale(norm/hAlphaBkg->Integral());
    hqtBkg->Scale(norm/hqtBkg->Integral());

    double Bkg3sigma=0;

    TCanvas *NocutMass=new TCanvas("NocutMass");
    NocutMass->SetName("NoCuttt");
    hmassBkg->Draw();
    hmassSig->SetMarkerColor(kOrange+7);
    hmassSig->Draw("Same");
    //Apply Cut Bkg
    for(int i=0 ; i<eventsBkg; i++){
        variablesBkg->GetEvent(i);
        bool CutOkBkg[NVar]={0};
        if(AlphaBkg<-0.4 && qtBkg<0.2){
            var[0]=abs(d02Bkg);
            var[1]=d0prodBkg;
            var[2]=distBkg;
            var[3]=cospBkg;
            var[4]=abs(d01Bkg);
            var[5]=dcaBkg;
            var[6]=ptMaxBkg;

            if(On_Off[0]==1 && var[0]<PreCut[0]) continue;
            if(On_Off[1]==1 && var[1]>PreCut[1]) continue;
            if(On_Off[2]==1 && var[2]<PreCut[2]) continue;
            if(On_Off[3]==1 && var[3]<PreCut[3]) continue;
            if(On_Off[4]==1 && var[4]<PreCut[4]) continue;
            if(On_Off[5]==1 && var[5]<PreCut[5]) continue;
            if(On_Off[6]==1 && var[6]<PreCut[6]) continue;

            if(massBkg>=meanLow && massBkg<=meanHigh){
                Bkg3sigma++;
            }

            hmassBkgCut->Fill(massBkg);
            #if DEBUG
            
                hyBkgCut->Fill(yrapBkg);
                hptBkgCut->Fill(ptBkg);
                hptMinBkgCut->Fill(ptMinBkg);
                hptMaxBkgCut->Fill(ptMaxBkg);
                hdistBkgCut->Fill(distBkg);  
                hcospBkgCut->Fill(cospBkg); 
                hd01BkgCut->Fill(d01Bkg);  
                hd02BkgCut->Fill(d02Bkg);  
                hd0prodBkgCut->Fill(d0prodBkg);  
                hdcaBkgCut->Fill(dcaBkg);
            #endif

            // Single Scan 
           for(int iVar=0;iVar<NVar;iVar++){
                for(int iScan=0;iScan<NStep;iScan++){
                    if (On_OffScan[iVar]==1 && high_lowScan[iVar]=='+'){
                    //printf("Scanmin %f ScanMax %f Step%f \n",ScanMin[iVar],ScanMax[iVar],StepScan[iVar]);
                        //printf("scan BKG %f %f",MxCutScan[iVar][iScan],var[iVar]);
                        if(var[iVar]>MxCutScan[iVar][iScan] && massBkg>=meanLow && massBkg<=meanHigh){
                            B[iVar][iScan]++;
                        }
                    }
                    if (On_OffScan[iVar]==1 && high_lowScan[iVar]=='-'){ 
                        if(var[iVar]<MxCutScan[iVar][iScan] && massBkg>=meanLow && massBkg<=meanHigh){   
                            B[iVar][iScan]++;
                        }
                    } 
                }  
            }

            if(ScanX2==1){
                //Double Scan
                bin=0;
                double cut1=0, cut2=0;
                for(int iCut1=0;iCut1<NStep;iCut1++){
                    cut1=MxCutScan[pos1][iCut1];
                    for(int iCut2=0;iCut2<NStep;iCut2++){
                        cut2=MxCutScan[pos2][iCut2];
                        if(high_lowScan[pos1]=='+' && high_lowScan[pos2]=='+'){
                            if(var[pos1]>cut1 && var[pos2]>cut2 && massBkg>=meanLow && massBkg<=meanHigh){
                                //printf("var1 %s var2 %s, cut1 %f cut2 %f \n",ScanVariable[pos1].c_str(),ScanVariable[pos2].c_str(),cut1,cut2);
                                bin = hBkgScan->GetBin(iCut1+1,iCut2+1,0);
                                hBkgScan->AddBinContent(bin);
                                Bcounts[iCut1][iCut2]++;
                            }
                        }
                        if(high_lowScan[pos1]=='+' && high_lowScan[pos2]=='-'){
                            if(var[pos1]>cut1 && var[pos2]<cut2 && massBkg>=meanLow && massBkg<=meanHigh){
                                //printf("var1 %s var2 %s cut1 %f cut2 %f \n",ScanVariable[pos1].c_str(),ScanVariable[pos2].c_str(),cut1,cut2);
                                bin = hBkgScan->GetBin(iCut1+1,iCut2+1,0);
                                hBkgScan->AddBinContent(bin);
                                Bcounts[iCut1][iCut2]++;
                            }
                        }
                        if(high_lowScan[pos1]=='-' && high_lowScan[pos2]=='+'){
                            if(var[pos1]<cut1 && var[pos2]>cut2 && massBkg>=meanLow && massBkg<=meanHigh){
                                //printf("var1 %s var2 %s \n",ScanVariable[pos1].c_str(),ScanVariable[pos2].c_str());
                                bin = hBkgScan->GetBin(iCut1+1,iCut2+1,0);
                                hBkgScan->AddBinContent(bin);
                                Bcounts[iCut1][iCut2]++;
                            }
                        }
                        if(high_lowScan[pos1]=='-' && high_lowScan[pos2]=='-'){
                            if(var[pos1]<cut1 && var[pos2]<cut2 && massBkg>=meanLow && massBkg<=meanHigh){
                                //printf("var1 %s var2 %s \n",ScanVariable[pos1].c_str(),ScanVariable[pos2].c_str());
                                bin = hBkgScan->GetBin(iCut1+1,iCut2+1,0);
                                hBkgScan->AddBinContent(bin);
                                Bcounts[iCut1][iCut2]++;
                            }
                        }
                    }
                }
            }
        }
    }    
       

  /*  for(int i=0;i<NVar;i++){
        for(int j=0;j<11;j++){
            //S[i][j]=S[i][j]*EvSig;
            //B[i][j]=B[i][j]/EvBkg;
            printf("%d) Bkg cut %f bkg cut/tot %f, bkgSignRange %f \n",i, B[i][j],B[i][j]/Btot, B[i][j]/Bkgentries);
        }
    } */

    //SCANSIONI 2D
    if(ScanX2==1){
        //Candidate sopravvissute
        hBkgScan->SetDirectory(0);
        hBkgScan->GetXaxis()->SetTitle(Form("%s ",ScanVariable[pos1].c_str()));
        hBkgScan->GetYaxis()->SetTitle(Form("%s ",ScanVariable[pos2].c_str()));
        hBkgScan->SetEntries(1);
        /*new TCanvas();
        hBkgScan->Draw("COLZ");*/

        //Frazione di fondo rimasto nel range di massa
        TH2F* hfrBkg = (TH2F*)hBkgScan->Clone("hfrBkg");
        hfrBkg->SetDirectory(0);
        hfrBkg->Scale(1/Bkgentries);
        hfrBkg->SetTitle("Fraction of Bkg");
        Min1=ScanMin[pos1];
        Min2=ScanMin[pos2];
        for(int iLab=0;iLab<NStep;iLab++){
            hBkgScan->GetXaxis()->SetBinLabel(iLab+1,Form("%f",Min1));
            hBkgScan->GetYaxis()->SetBinLabel(iLab+1,Form("%f",Min2));

            hfrBkg->GetXaxis()->SetBinLabel(iLab+1,Form("%f",Min1));
            hfrBkg->GetYaxis()->SetBinLabel(iLab+1,Form("%f",Min2));
            Min1=Min1+StepScan[pos1];
            Min2=Min2+StepScan[pos2];
            }
        TCanvas *cfr=new TCanvas("cfr");
        hfrBkg->SetStats(0);
        hfrBkg->Draw("COLZ");

        TH2F* hS_B = new TH2F("hS_B","Signal over Bkg",NStep,-0.5,NStep-0.5,NStep,-0.5,NStep-0.5);
        TH2F* hSignif = new TH2F("hSignif","Significance",NStep,-0.5,NStep-0.5,NStep,-0.5,NStep-0.5);
        for(int iStep=0;iStep<NStep;iStep++){
            for(int jStep=0;jStep<NStep;jStep++){
                double S=(Scounts[iStep][jStep])*EvSig*BR;//normalizzo al signal
                double B=(Bcounts[iStep][jStep])/EvBkg;//normalizzo al num eventi fondo
                double Signif=((sqrt(5*1e9))*S)/sqrt(S+B);
                hSignif->Fill(iStep,jStep,Signif);

                if(B<1e-6){
                    B=-10;
                    hS_B->Fill(iStep,jStep,S/B);
                }else{
                    hS_B->Fill(iStep,jStep,S/B);
                    
                }
                hS_B->SetMinimum(0);
            }
        }

        hSignif->SetDirectory(0);
        hSignif->SetTitle("Significance (3#sigma)");
        hSignif->GetXaxis()->SetTitle(Form("%s ",ScanVariable[pos1].c_str()));
        hSignif->GetYaxis()->SetTitle(Form("%s ",ScanVariable[pos2].c_str()));

        hS_B->SetDirectory(0);
        hS_B->SetTitle("Signal Over Bkg");
        hS_B->GetXaxis()->SetTitle(Form("%s ",ScanVariable[pos1].c_str()));
        hS_B->GetYaxis()->SetTitle(Form("%s ",ScanVariable[pos2].c_str()));

        Min1=ScanMin[pos1];
        Min2=ScanMin[pos2];
        for(int iLab=0;iLab<NStep;iLab++){
            hSignif->GetXaxis()->SetBinLabel(iLab+1,Form("%f",Min1));
            hSignif->GetYaxis()->SetBinLabel(iLab+1,Form("%f",Min2));

            hS_B->GetXaxis()->SetBinLabel(iLab+1,Form("%f",Min1));
            hS_B->GetYaxis()->SetBinLabel(iLab+1,Form("%f",Min2));

            Min1=Min1+StepScan[pos1];
            Min2=Min2+StepScan[pos2];
            }
        
        TCanvas *cSignif=new TCanvas("cSignif");
        hSignif->SetStats(0);
        hSignif->Draw("COLZ");

        TCanvas *cSB=new TCanvas("cSB");
        hS_B->SetStats(0);
        hS_B->Draw("COLZ");
    }


double Bkg3sigmaXevent=Bkg3sigma/EvBkg;
printf("Counts %f, countsXevent %f \n",Bkg3sigma,Bkg3sigmaXevent);
printf("Significance %f ,S=%f, B= %f \n",SigXevent/(sqrt(SigXevent+Bkg3sigmaXevent)),SigXevent,Bkg3sigmaXevent);
hmassBkgCut->SetDirectory(0);
hmassBkgCut->Scale(norm/EvBkg);
TCanvas *cSigBkgCut = new TCanvas("cSigBkgCut");
hmassBkgCut->SetMarkerColor(kAzure+7);
hmassBkgCut->SetMarkerStyle(22);
hmassBkgCut->SetMarkerSize(1);
hmassBkgCut->Draw();
hmassSigCut->SetMarkerColor(kOrange+7);
hmassSigCut->SetMarkerStyle(20);
hmassSigCut->SetMarkerSize(0.8);
hmassSigCut->Draw("SAME");

TCanvas *cBkg=new TCanvas("cBkg");
TF1 *f = new TF1("f","pol2",1.07,1.155);
f->SetRange(1.07,1.155);

hmassBkgCut->Fit(f,"","",1.07,1.155);
hmassBkgCut->Draw();
double IntParab=f->Integral(1.07,1.155);
double IntMass=f->Integral(meanLow,meanHigh);

printf("par %f, mass %f \n",IntParab,IntMass);
/*double p0 = f->GetParameter(0);
double p1 = f->GetParameter(1);
double p2 = f->GetParameter(2);

TF1 *fprov=new TF1("f","pol2",1.1, 1.155);
fprov->SetParameter(0,p0);
fprov->SetParameter(1,p1);
fprov->SetParameter(2,p2);*/


//SCANSIONE TAGLI PER OTTIMIZZAZIONE 1D
    Plot(S, B, Bsig, Stot, Bkgentries, fout);

    #if DEBUG    
        hyBkgCut->SetDirectory(0);
        hptBkgCut->SetDirectory(0);
        hptMinBkgCut->SetDirectory(0);
        hptMaxBkgCut->SetDirectory(0);
        hdistBkgCut->SetDirectory(0);
        hd01BkgCut->SetDirectory(0);
        hd02BkgCut->SetDirectory(0);
        hd0prodBkgCut->SetDirectory(0);
        hdcaBkgCut->SetDirectory(0);
        hcospBkgCut->SetDirectory(0);
        hArmenBkgCut->SetDirectory(0);
        hAlphaBkgCut->SetDirectory(0);
        hqtBkgCut->SetDirectory(0);

        hyBkgCut->Scale(norm/hyBkgCut->Integral());
        hptBkgCut->Scale(norm/hptBkgCut->Integral());
        hptMinBkgCut->Scale(norm/hptMinBkgCut->Integral());
        hptMaxBkgCut->Scale(norm/hptMaxBkgCut->Integral());
        hdistBkgCut->Scale(norm/hdistBkgCut->Integral());
        hd01BkgCut->Scale(norm/hd01BkgCut->Integral());
        hd02BkgCut->Scale(norm/hd02BkgCut->Integral());
        hcospBkgCut->Scale(norm/hcospBkgCut->Integral());
        hd0prodBkgCut->Scale(norm/hd0prodBkgCut->Integral());
        hdcaBkgCut->Scale(norm/hdcaBkgCut->Integral());
        hAlphaBkgCut->Scale(norm/hAlphaBkgCut->Integral());
        hqtBkgCut->Scale(norm/hqtBkgCut->Integral());

            TCanvas *cut=new TCanvas("cut");
            cut->SetName("PostCut");
    cut->Divide(3,2);

    cut->cd(1);
    hd02BkgCut->SetMarkerColor(kAzure+7);
    hd02BkgCut->SetMarkerStyle(22);
    hd02BkgCut->Draw();
    hd02SigCut->SetMarkerColor(kOrange+7);
    hd02SigCut->SetMarkerStyle(20);
    hd02SigCut->SetMarkerSize(0.8);
    hd02SigCut->SetMarkerColor(kOrange+7);
    hd02SigCut->SetFillColor(0);
    hd02SigCut->SetFillStyle(0);
    hd02SigCut->SetLineWidth(1);
    hd02SigCut->SetLineColor(kOrange+7);
    hd02SigCut->Draw("Same");
    cut->cd(2);
    hd01BkgCut->SetMarkerColor(kAzure+7);
    hd01BkgCut->SetMarkerStyle(22);
    hd01BkgCut->Draw();
    hd01SigCut->SetMarkerColor(kOrange+7);
    hd01SigCut->SetMarkerStyle(20);
    hd01SigCut->SetMarkerSize(0.8);
    hd01SigCut->SetMarkerColor(kOrange+7);
    hd01SigCut->SetFillColor(0);
    hd01SigCut->SetFillStyle(0);
    hd01SigCut->SetLineWidth(1);
    hd01SigCut->SetLineColor(kOrange+7);
    hd01SigCut->Draw("Same");
    cut->cd(3);
    hd0prodBkgCut->SetMarkerColor(kAzure+7);
    hd0prodBkgCut->SetMarkerStyle(22);
    hd0prodBkgCut->Draw();
    hd0prodSigCut->SetMarkerColor(kOrange+7);
    hd0prodSigCut->SetMarkerStyle(20);
    hd0prodSigCut->SetMarkerSize(0.8);
    hd0prodSigCut->SetMarkerColor(kOrange+7);
    hd0prodSigCut->SetFillColor(0);
    hd0prodSigCut->SetFillStyle(0);
    hd0prodSigCut->SetLineWidth(1);
    hd0prodSigCut->SetLineColor(kOrange+7);
    hd0prodSigCut->Draw("Same");
    hd0prodSigCut->Draw("Same");
    cut->cd(4);
    hdistBkgCut->SetMarkerColor(kAzure+7);
    hdistBkgCut->SetMarkerStyle(22);
    hdistBkgCut->Draw();
    hdistSigCut->SetMarkerColor(kOrange+7);
    hdistSigCut->Draw("Same");
    cut->cd(5);
    hcospBkgCut->SetMarkerColor(kAzure+7);
    hcospBkgCut->SetMarkerStyle(22);
    hcospBkgCut->Draw();
    hcospSigCut->SetMarkerColor(kOrange+7);
    hcospSigCut->Draw("Same");
    cut->cd(6);
    hdcaBkgCut->SetMarkerColor(kAzure+7);
    hdcaBkgCut->SetMarkerStyle(22);
    hdcaBkgCut->Draw();
    hdcaSigCut->SetMarkerColor(kOrange+7);
    hdcaSigCut->Draw("Same");


    #endif
    //printf("BtotpostArment= %f \n",hmassBkgCut->GetEntries());
    //printf("Btot= %f, Btoy %f, BcountsPostArm %f \n",Btot,BtotInt,Bkgentries);
    //printf("BcountsPostArm %f \n",Bkgentries);

    //fout->Close();



 /*   TCanvas *Nocut=new TCanvas("Nocut");
    Nocut->Divide(3,2);
    Nocut->cd(1);
    hd02Bkg->SetMarkerColor(kAzure+7);
    hd02Bkg->SetMarkerStyle(22);
    hd02Bkg->Draw();
    hd02Sig->SetMarkerColor(kOrange+7);
    hd02Sig->SetMarkerStyle(20);
    hd02Sig->Draw("Same");
    Nocut->cd(2);
    hd01Bkg->SetMarkerColor(kAzure+7);
    hd01Bkg->SetMarkerStyle(22);
    hd01Bkg->Draw();
    hd01Sig->SetMarkerColor(kOrange+7);
    hd01Sig->SetMarkerStyle(20);
    hd01Sig->Draw("Same");
    Nocut->cd(3);
    hd0prodBkg->SetMarkerColor(kAzure+7);
    hd0prodBkg->SetMarkerStyle(22);
    hd0prodBkg->Draw();
    hd0prodSig->SetMarkerColor(kOrange+7);
    hd0prodSig->SetMarkerStyle(20);
    hd0prodSig->Draw("Same");
   Nocut->cd(4);
    hdistBkg->SetMarkerColor(kAzure+7);
    hdistBkg->SetMarkerStyle(22);
    hdistBkg->Draw();
    hdistSig->SetMarkerColor(kOrange+7);
    hdistSig->SetMarkerStyle(20);
    hdistSig->Draw("Same");
    Nocut->cd(5);
    hcospBkg->SetMarkerColor(kAzure+7);
    hcospBkg->SetMarkerStyle(22);
    hcospBkg->Draw();
    hcospSig->SetMarkerColor(kOrange+7);
    hcospSig->SetMarkerStyle(20);
    hcospSig->Draw("Same");
    Nocut->cd(6);
    hdcaBkg->SetMarkerColor(kAzure+7);
    hdcaBkg->SetMarkerStyle(22);
    hdcaBkg->Draw();
    hdcaSig->SetMarkerColor(kOrange+7);
    hdcaSig->SetMarkerStyle(20);
    hdcaSig->Draw("Same");*/

TH1F* hBkgGen=new TH1F("hBkgGen","Inv mass Bkg",200,1.1,1.16);
int value=0;
printf("Ev fondo %f \n",((double)Bkg3sigmaXevent));
Long64_t NBkg=0;

/*TRandom3 *r3=new TRandom3();
for(int i=0;i<200;i++){
    double IntBin=f->Integral(hBkgGen->GetBinLowEdge(i+1),hBkgGen->GetBinLowEdge(i+2));
    NBkg=(Long64_t)(Bkg3sigmaXevent*(IntBin/IntMass)*5*1e9);
    value=r3->Poisson(NBkg);
    hBkgGen->SetBinContent(i+1,value);
    hBkgGen->SetBinError(i+1,sqrt(value));
}*/

Long64_t BkgGen=5*1e9*Bkg3sigmaXevent*IntParab/IntMass;
double valueB=0;
for(Long64_t i=0;i<BkgGen;i++){
    valueB=f->GetRandom();
    hBkgGen->Fill(valueB);  
}
hBkgGen->SetDirectory(0);
TCanvas *cBkgGen=new TCanvas("cBkgGen");
hBkgGen->Draw();

Long64_t Bcheck=0;
for(int i=hBkgGen->FindBin(meanLow); i<hBkgGen->FindBin(meanHigh);i++){
    Bcheck=Bcheck+(hBkgGen->GetBinContent(i));
}
printf("Bcheck %lld \n",Bcheck);



TH1F* hSigGen=new TH1F("hSigGen","Inv mass Sig",200000,1.1,1.16);
double value1=0;
hSigGen->SetDirectory(0);
Long64_t NSig=5*SigXevent*1e9;
//printf("%d \n",NSig);

for(Long64_t i=0;i<NSig;i++){
    value1=f1->GetRandom();
    hSigGen->Fill(value1);  
}
printf("Sign %lld \n",NSig);
TCanvas *cSigGen=new TCanvas("cSigGen");
hSigGen->SetMarkerStyle(20);
hSigGen->Draw();

Long64_t Scheck=0;
for(int i=hSigGen->FindBin(meanLow); i<hSigGen->FindBin(meanHigh);i++){
    Scheck=Scheck+(hSigGen->GetBinContent(i));
}
printf("Scheck %lld \n",Scheck);
  TList *list = new TList;
  list->Add(hBkgGen);
  list->Add(hSigGen);
TH1F* hSum=(TH1F*)hBkgGen->Clone("hSum");
  hSum->Reset();
  hSum->Merge(list);
  hSum->Draw();

TCanvas *cSum = new TCanvas("cSum");
hSum->SetDirectory(0);
TF1* fSum=new TF1("fSum","gaus+pol2(3)",1.1,1.16);
fSum->SetParLimits(0,1e8,3e9);
fSum->SetParLimits(1,1.114,1.117);
fSum->SetParLimits(2,1e-3,2e-3);
hSum->Fit(fSum);
//TF1* fSumBkg=new TF1("fSumBkg","pol2",1.1,1.16);
//hSum->Fit(fSumBkg,"","same",1.1,1.16);
hSum->SetTitle("Invariant Mass");
hSum->Draw("same");

printf("signif %f",((Scheck)/(sqrt(Scheck+Bcheck))));

  TLatex* lat2 = new TLatex();
  lat2->SetTextFont(42);
  lat2->SetTextSize(0.04);
  lat2->SetNDC();
lat2->DrawLatex(0.5,0.85,"PbPb, #sqrt{#it{s}} = 17.3 GeV, centrality 0-5%");
lat2->DrawLatex(0.5,0.75,Form("#mu %f ",fSum->GetParameter(1)));
lat2->DrawLatex(0.5,0.65,Form("#sigma %f ",fSum->GetParameter(2)));
lat2->DrawLatex(0.5,0.55,Form("Significance(3#sigma) %f #pm %f",((sqrt(5*1e9))*Scheck*1e-9/5)/(sqrt((Scheck*1e-9/5)+(Bcheck*1e-9/5))),sqrt((1e5*Scheck*1e-9/5)/(sqrt((Scheck*1e-9/5)+(Bcheck*1e-9/5))))));
fout->cd();
hSum->Write();


clock_t end=clock();
cout<<"Execution time: "<<((double)(end-start)/CLOCKS_PER_SEC)<<"sec"<<endl;

filinBkg.Close();



}

void DrawSame(TH1F* hBkg, TH1F* hSig, TCanvas *c, const char* sTitle){

    hBkg->SetMarkerStyle(23);
    hBkg->SetMarkerSize(1);
    hBkg->SetMarkerColor(kAzure-3);
    hBkg->SetFillColor(0);
    hBkg->SetFillStyle(0);
    hBkg->SetLineWidth(1);
    hBkg->SetLineColor(kAzure-3);
    hBkg->SetTitle(Form("Overlap of %s for Bkg and Sig",sTitle));
    
    gPad->SetGrid();
    hBkg->Draw();

    hSig->SetMarkerStyle(20);
    hSig->SetMarkerSize(0.8);
    hSig->SetMarkerColor(kOrange+7);
    hSig->SetFillColor(0);
    hSig->SetFillStyle(0);
    hSig->SetLineWidth(1);
    hSig->SetLineColor(kOrange+7);
    gPad->SetGrid();
    hSig->Draw("same");
}

void DrawFit( TH1F* hSig, TCanvas *c, const char* sTitle, double val[2]){


    hSig->SetMarkerStyle(20);
    hSig->SetMarkerColor(kOrange+7);
    hSig->SetFillColor(0);
    hSig->SetFillStyle(0);
    hSig->SetLineWidth(1);
    hSig->SetLineColor(kOrange+7);
    hSig->Draw();

    TF1 *f2 = new TF1("f2","gaus",1.06,1.18);
    hSig->Fit(f2,"","",1.06,1.18);
    double mean=f2->GetParameter(1);
    double sigma=f2->GetParameter(2);
    val[0]=mean;
    val[1]=sigma;
    //printf("fit %f %f \n",mean,sigma);
}

void ReadPreCut(string nameVar[NVar], double cutVar[NVar], bool set[NVar], char cond[NVar]){

    std::ifstream filef("PreCut.txt");
    int i=0;
    while (filef >> nameVar[i] >> cutVar[i] >> set[i] >> cond[i]){
        //std::cout << nameVar[i] <<"\t"<< cutVar[i] <<"\t"<< set[i] <<endl;
        i++;
    }    
}

void ReadScanCut(string nameVarScan[NVar], double Stepcut[NVar], double min[NVar], double max[NVar], bool setScan[NVar], char condScan[NVar]){
        std::ifstream fileScan("ScanCut.txt");
    int i=0;
    while (fileScan >> nameVarScan[i] >> Stepcut[i] >> min[i] >> max[i] >> setScan[i] >> condScan[i]){
        //std::cout << nameVar[i] <<"\t"<< cutVar[i] <<"\t"<< set[i] <<endl;
        i++;
    }  
}

bool AreEqual(bool SetCut[NVar], bool CheckCut[NVar]){
    for(int i=0;i<NVar;i++){
        if(SetCut[i]!=CheckCut[i]){
            return false;
        }
    }
    return true;
}
void Plot(double CutSig[NVar][NStep], double CutBkg[NVar][NStep],double Bkg[NVar][NStep], double TotSig, double TotBkg, TFile *f){

    std::ifstream fileScan("ScanCut.txt");
    string nameVarScan[NVar];
    double StepScan[NVar], ScanMin[NVar], ScanMax[NVar];
    bool On_OffScan[NVar];
    char condScan[NVar];
    int i=0;
    while (fileScan >> nameVarScan[i] >> StepScan[i] >> ScanMin[i] >> ScanMax[i] >> On_OffScan[i] >> condScan[i]){
        i++;
    }
    TH1F* hEffXdiffCut[NVar], *hSB[NVar], *hSignif[NVar], *hBkgRej[NVar];
    TCanvas *c[NVar];

    const char* name[NVar]={ "|d02| (cm)","d0prod (cm^{2})", "dist (cm)", "cosp ", "|d01| (cm)", "dca (cm)", "ptMax (GeV/c)"};
    for(int i=0;i<NVar;i++){

        c[i]=new TCanvas(Form("c%d",i));
        hEffXdiffCut[i]= new TH1F(Form("hEffXdiffCut%d",i),Form("%s scan for efficiency",nameVarScan[i].c_str()),NStep,-0.5,NStep-0.5);
        hSB[i] = new TH1F(Form("hSB%d",i),Form("%s scan for S/B",nameVarScan[i].c_str()),NStep,-0.5,NStep-0.5);
        hSignif[i] = new TH1F(Form("hSignif%d",i),Form("%s scan for Significance (3#sigma)",nameVarScan[i].c_str()),NStep,-0.5,NStep-0.5);
        hBkgRej[i] = new TH1F(Form("hBkgRej%d",i),Form("%s scan for Bkg rejected",nameVarScan[i].c_str()),NStep,-0.5,NStep-0.5);
        }
    f->cd();
    for(int k=0;k<NVar;k++){
        if(On_OffScan[k]==1){
            double Min=ScanMin[k];
            for(int i=0;i<NStep;i++){
                hEffXdiffCut[k]->Fill(i,((CutSig[k][i])/TotSig));
                hBkgRej[k]->Fill(i,((CutBkg[k][i])/(TotBkg)));
                double S =(CutSig[k][i]*(EvSig));
                double B =(CutBkg[k][i]/EvBkg);
                //printf("%d) S %f B %f signif %f \n",k,S,B,S/sqrt(S+B));
                hSignif[k]->Fill(i,(((sqrt(5*1e9))*(S))/sqrt((S)+(B))));
                if(B<1e-6){
                    B=-10;
                    hSB[k]->Fill(i,(S)/(B));
                }else{
                
                hSB[k]->Fill(i,(S)/(B));
                }
                hSB[k]->SetMinimum(0);
                //hSignif[k]->Fill(i,((1e5*(S))/sqrt(S)+(B)));
                //hSB[k]->Fill(i,(S)/(B));
                //printf("%d) %f %f %f \n",i+1,Min,CutBkg[k][i],TotBkg);
                hEffXdiffCut[k]->GetXaxis()->SetBinLabel(i+1,Form("%f",Min));
                hBkgRej[k]->GetXaxis()->SetBinLabel(i+1,Form("%f",Min));
                hSignif[k]->GetXaxis()->SetBinLabel(i+1,Form("%f",Min));
                hSB[k]->GetXaxis()->SetBinLabel(i+1,Form("%f",Min));
                Min=Min+StepScan[k];
                hEffXdiffCut[k]->SetStats(0);
                hBkgRej[k]->SetStats(0);
                hSignif[k]->SetStats(0);
                hSB[k]->SetStats(0);
                

                
            }
        }
    }
    
    for(int i=0;i<NVar;i++){
        if(On_OffScan[i]==1){
            hEffXdiffCut[i]->SetDirectory(0);
            hEffXdiffCut[i]->GetXaxis()->SetTitle(Form("%s",name[i]));
            hEffXdiffCut[i]->GetYaxis()->SetTitle("Efficiency*event");
            hEffXdiffCut[i]->SetMarkerStyle(20);
            hEffXdiffCut[i]->SetMarkerSize(1);
            hEffXdiffCut[i]->SetMarkerColor(kOrange+7);
            hEffXdiffCut[i]->SetLineWidth(0);
            
            c[i]->Divide(2,2);
            c[i]->cd(1);
            gPad->SetGrid();
            hEffXdiffCut[i]->Draw();

            hBkgRej[i]->SetDirectory(0);
            hBkgRej[i]->GetXaxis()->SetTitle(Form("%s",name[i]));
            hBkgRej[i]->GetYaxis()->SetTitle("Fraction of Bkg");
            hBkgRej[i]->SetMarkerStyle(20);
            hBkgRej[i]->SetMarkerSize(1);
            hBkgRej[i]->SetMarkerColor(kOrange+7);
            hBkgRej[i]->SetLineWidth(0);
            c[i]->cd(2);
            gPad->SetGrid();
            hBkgRej[i]->Draw();


            hSignif[i]->SetDirectory(0);
            hSignif[i]->GetXaxis()->SetTitle(Form("%s",name[i]));
            hSignif[i]->GetYaxis()->SetTitle("Significance");
            hSignif[i]->SetMarkerStyle(20);
            hSignif[i]->SetMarkerSize(1);
            hSignif[i]->SetMarkerColor(kOrange+7);
            hSignif[i]->SetLineWidth(0);
           c[i]->cd(3);
            gPad->SetGrid();
            hSignif[i]->Draw();

            hSB[i]->SetDirectory(0);
            hSB[i]->GetXaxis()->SetTitle(Form("%s",name[i]));
            hSB[i]->GetYaxis()->SetTitle("S/B");
            hSB[i]->SetMarkerStyle(20);
            hSB[i]->SetMarkerSize(1);
            hSB[i]->SetMarkerColor(kOrange+7);
            hSB[i]->SetLineWidth(0);
            c[i]->cd(4);
            gPad->SetGrid();
            hSB[i]->Draw();
                hEffXdiffCut[i]->Write();
                hBkgRej[i]->Write();
                hSignif[i]->Write();
                hSB[i]->Write();
        }
    }
}


