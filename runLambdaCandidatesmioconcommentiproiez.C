#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TString.h>
#include <TTree.h>
#include <TArrayF.h>
#include <TMath.h>
#include <TH1F.h>
#include <TNtuple.h>
#include <TFile.h>
#include "KMCProbeFwd.h"
#include "KMCDetectorFwd.h"
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "TRandom.h"
#include "TF1.h"
#include "THnSparse.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "TROOT.h"
#include "GenMUONLMR.h"
#include "TStopwatch.h"
#include "TTree.h"
#include "TParticle.h"
#include "AliAODRecoDecay.h"
#include "AliDecayer.h"
#include "AliDecayerEvtGen.h"
#include "TDatabasePDG.h"
#include "TGraph.h"
#include "TGraphErrors.h"

#endif

// Track Chi2Tot cut
double ChiTot = 1.5;

// settings for signal generation CONTROLLA I VALORI PER SIGNAL GENERATION 
double yminSG = -10.; // min y to generate
double ymaxSG = 10.;  //
double ptminSG = 0.;
double ptmaxSG = 10; //Elena's change, it was 3 GeV/c

double vX = 0, vY = 0, vZ = 0; // event vertex

THnSparseF* CreateSparse();//istogramma multidimensionale
void ComputeVertex(KMCProbeFwd &t0, KMCProbeFwd &t1, Double_t &xV, Double_t &yV, Double_t &zV);//calcolo vertice
Double_t CosPointingAngle(Double_t vprim[3], Double_t vsec[3], TLorentzVector &parent);
TDatime dt;
Double_t ImpParXY(Double_t vprim[3], Double_t vsec[3], TLorentzVector &parent);
Double_t CosThetaStar(TLorentzVector &parent, TLorentzVector &dauk);
double weightmean(TH1D* hproject);


void GenerateD0SignalCandidates(Int_t nevents = 100000, 
				double Eint = 160., 
				const char *setup = "setup-10um-itssa_Eff1.txt", 
				const char *filNamPow="/home/prino/na60plus/POWHEG/pp20/Charm1dot5/pp0_frag-PtSpectra-Boost.root",//dovrò cambiare file, dove ho tutte le info pt,y,...
				const char *privateDecayTable = "/home/prino/na60plus/decaytables/USERTABD0.DEC",//check decay table devo aggiungere quella della lambda
				int optPartAntiPart=3,
				bool writeNtuple = kFALSE, //???????
				bool simulateBg=kTRUE){
  
  // Generate D0->Kpi signals and simulate detector response for decay tracks
  // Need input D0 pt and y ditribution from POWHEG


  int refreshBg = 1000;//cosa è?
  static UInt_t seed = dt.Get();//genero seed sempre diverso tramite date/time con TDatime
  gRandom->SetSeed(seed);
  gSystem->Load("$ALICE_ROOT/lib/libEvtGen.so");//carico libreria.EvtGen is a Monte Carlo event generator that simulates the decays of heavy flavour particles
  gSystem->Load("$ALICE_ROOT/lib/libEvtGenExternal.so");
  gSystem->Load("$ALICE_ROOT/lib/libTEvtGen.so");
  

  //  PYTHIA input -> not used
  // TFile *fin = new TFile("Mergedfkine.root");
  // // TFile *fin = new TFile("fkineNew.root");
  // TH1D *hD0pt = (TH1D *)fin->Get("hD0pt");
  // TH1D *hD0y = (TH1D *)fin->Get("hD0y");
  // hD0y->Rebin(4);
  
  //  POWHEG+PYTHIA input 
  printf("--> pt and y shape of Lambda from %s\n",filNamPow);
  TFile *filPow=new TFile(filNamPow);
  TH1D *hLambdapt = (TH1D*)filPow->Get("hpt");
  TH1D *hLambday = (TH1D*)filPow->Get("hy");
	printf("tutto beneeeeeeee yeeeeeeeeeeeee \n");
  TH1D *hLambdabarpt = (TH1D*)filPow->Get("hptnotl");
  TH1D *hLambdabary = (TH1D*)filPow->Get("hynotl");
	printf("tutto beneeeeeeee NOTLAMBDAAAAAAA \n");
  /*TH3D* h3Dbarpow=(TH3D*)filPow->Get("hptyetam421");
  if(h3Dbarpow){
    TH1D *hD0barpt = (TH1D*)h3Dbarpow->ProjectionX("hD0barpt");
    TH1D *hD0bary = (TH1D*)h3Dbarpow->ProjectionY("hD0bary");*/
  /*if(optPartAntiPart==2){
      hLambdapt=hLambdabarpt;
      hLambday=hLambdabary;
    }*/
/*  if(optPartAntiPart==3){
      hLambdapt->Add(hLambdabarpt);
      hLambday->Add(hLambdabary);
    }else if(optPartAntiPart==2){
      hLambdapt=hLambdabarpt;
      hLambday=hLambdabary;
    }
  }*/
/*  TFile *fntupla = new TFile("ntupla.root", "recreate");
	float variable[16];
	TNtuple *nt = new TNtuple("nt","Variablesforperformance","nfaketrkprot:nfaketrkpion:ygen:xP:yP:zP:Vxgen:Vygen:Vzgen:massinvrec:diffxprot:diffyprot:diffzprot:diffxpion:diffypion:diffzpion");
	
   TFile *hfile = new TFile("Nfake.root","RECREATE");
   TTree *treeNfake = new TTree("treeNfake","Nfake");*/
   
//istogrammi dei prodotti che saranno salvati in decayhistos.root
  TH2F *hptP = new TH2F("hptP", "Protons from Lambda decays", 50,0.,10.,50, 0., 10.);
  TH2F *hptPi = new TH2F("hptPi", "pions from Lambda decays", 50, 0.,10.,50,0., 10.);
  TH1D *hyP = new TH1D("hyP", "y Protons from Lambda decays", 50, 0., 5.);
  TH1D *hyPi = new TH1D("hyPi", "y pions from Lambda decays", 50, 0., 5.);
  TH2F *hyPiP = new TH2F("hyPiP", "y pions vs y Protons from Lambda decays", 50, 0., 5., 50, 0., 5.);

  TFile *fout = new TFile("LambdaNorec-Signal-histos.root", "recreate");
   //int outN = nev/10;
  //if (outN<1) outN=1;

//_____________________________SETUP SPERIMENTALE
  KMCDetectorFwd *det = new KMCDetectorFwd();
  det->ReadSetup(setup, setup);
  det->InitBkg(Eint); //check modificare bkg su altro file
  det->ForceLastActiveLayer(det->GetLastActiveLayerITS()); // will not propagate beyond VT

  det->SetMinITSHits(det->GetNumberOfActiveLayersITS()); //NA60+
  //det->SetMinITSHits(det->GetNumberOfActiveLayersITS()-1); //NA60
  det->SetMinMSHits(0); //NA60+
  //det->SetMinMSHits(det->GetNumberOfActiveLayersMS()-1); //NA60
  det->SetMinTRHits(0);
  //
  // max number of seeds on each layer to propagate (per muon track)
  det->SetMaxSeedToPropagate(3000);
  //
  // set chi2 cuts
  det->SetMaxChi2Cl(10.);  // max track to cluster chi2
  det->SetMaxChi2NDF(3.5); // max total chi2/ndf
  det->SetMaxChi2Vtx(20);  // fiducial cut on chi2 of convergence to vtx  
  // IMPORTANT FOR NON-UNIFORM FIELDS
  det->SetDefStepAir(1);
  det->SetMinP2Propagate(1); //NA60+
  //det->SetMinP2Propagate(2); //NA60
  //
  det->SetIncludeVertex(kFALSE); // count vertex as an extra measured point
  //  det->SetApplyBransonPCorrection();
  det->ImposeVertex(0., 0., 0.);
  det->BookControlHistos();
  //
  
  //Magnetic field and detector parameters_______________mantieni uguale
 //Magnetic field and detector parameters
 TVirtualMagField* fld = TGeoGlobalMagField::Instance()->GetField();
  if (fld->IsA() == MagField::Class()) {
    MagField* mag = (MagField*) fld;
    int BNreg = mag->GetNReg();
    const double *BzMin = mag->GetZMin();
    const double *BzMax = mag->GetZMax();
    const double *BVal;
    printf("*************************************\n");
    printf("number of magnetic field regions = %d\n", BNreg);
    for (int i = 0; i < BNreg; i++){
      BVal = mag->GetBVals(i);
      printf("*** Field region %d ***\n", i);
      if (i == 0){
	printf("Bx = %f B = %f Bz = %f zmin = %f zmax = %f\n", BVal[0], BVal[1], BVal[2], BzMin[i], BzMax[i]);
      }else if (i == 1){
	printf("B = %f Rmin = %f Rmax = %f zmin = %f zmax = %f\n", BVal[0], BVal[1], BVal[2], BzMin[i], BzMax[i]);
      }
    }
  }
  

  
 //printf("pippo");
  
 

  // prepare decays
  TGenPhaseSpace decay;//Utility class to generate n-body event, with constant cross-section (default) or with Fermi energy dependence (opt="Fermi"). The event is generated in the center-of-mass frame, 												but the decay products are finally boosted using the betas of the original particle.

  TLorentzVector parentgen, daugen[2], parent, daurec[2], parentrefl, daurecswapmass[2]; //quadrimpulso ultime due che significato hanno?
  KMCProbeFwd recProbe[2];  
  AliDecayerEvtGen *fDecayer = new AliDecayerEvtGen();
  fDecayer->Init(); //read the default decay table DECAY.DEC and particle table
  bool privTab=kFALSE;
  if (strlen(privateDecayTable)>0){
    if(gSystem->Exec(Form("ls -l %s",privateDecayTable))==0){//???
      fDecayer->SetDecayTablePath((char*)privateDecayTable);
      fDecayer->ReadDecayTable();
      printf("-- Use  decay table from file %s\n",privateDecayTable);
      privTab=kTRUE;
    }
  }
 /* if(!privTab){
    printf("-- Use existing decay modes in aliroot\n");
    fDecayer->SetForceDecay(kHadronicD); 
  }*/
  fDecayer->ForceDecay();

  //TParticleClass: Description of the dynamic properties of a particle. A dynamic particle class created by event generators and used during the propagation in detectors.
  TClonesArray *particles = new TClonesArray("TParticle", 1000);

  TLorentzVector *mom = new TLorentzVector(); //quadrimpulso particella lambda
  
  Int_t pdgParticle = 3122;// define mother particle
  
  TH2F* hYPtGen = new TH2F("hYPtGen", "Y-Pt corr match", 80, 1.0, 5.4, 40, ptminSG, ptmaxSG);// ?
  TH1D* hPtGen = new TH1D("hPtGen", "Pt gen", 40, ptminSG, ptmaxSG);
  TH1D* hYGen = new TH1D("hYGen", "Y full phase space", 80., -3.0, 7.0);
  TH2F* hYPtRecoAll = new TH2F("hYPtRecoAll", "Y-Pt all match", 80, 1.0, 5.4, 40, ptminSG, ptmaxSG);
  TH1D* hPtRecoAll = new TH1D("hPtRecoAll", "Reconstructed Pt all match", 40, ptminSG, ptmaxSG);  
  TH1D* hPtGenRecoAll = new TH1D("hPtGenRecoAll", "Generated Pt all match", 40, ptminSG, ptmaxSG);
  TH2F* hPtRecoVsGenAll = new TH2F("hPtRecoVsGenAll"," ; Generated p_{T} ; Reconstructed p_{T}",40, ptminSG, ptmaxSG,40, ptminSG, ptmaxSG);
  TH2F* hDiffPtRecoGenAll = new TH2F("hDiffPtRecoGenAll"," ; Generated p_{T} ; Reco p_{T} - Gen p_{T}",40, ptminSG, ptmaxSG,100,-0.2,0.2);

  TH1D* hYRecoAll = new TH1D("hYRecoAll", "Reconstructed Y all match", 80., 1., 5.4);
  TH1D* hYGenRecoAll = new TH1D("hYGenRecoAll", "Generated Y all match", 80., 1., 5.4);
  TH2F* hYPtRecoFake = new TH2F("hYPtRecoFake", "Y-Pt fake match", 80, 1.0, 5.4, 40, ptminSG, ptmaxSG);
  TH1D* hPtRecoFake = new TH1D("hPtRecoFake", "Pt fake match", 40, ptminSG, ptmaxSG);
  TH1D* hMassAll = new TH1D("hMassAll", "Mass all match", 200, 0.5, 3.5);//massa
  TH1D* hMassFake = new TH1D("hMassFake", "Mass fake match", 200, 1., 3.5);
  TH1D* hMassRefl = new TH1D("hMassRefl", "Mass reflections", 200, 1., 3.5);

  TH2F *hDistXY = new TH2F("hDistXY", "", 100, 0, 0.1, 30, 0, 3);
  TH2F *hDist = new TH2F("hDist", "", 300, 0, 10, 30, 0, 3);
  TH2F *hDistgenXY = new TH2F("hDistgenXY", "", 100, 0, 0.1, 30, 0, 3);
  TH2F *hDistgen = new TH2F("hDistgen", "", 300, 0, 10, 30, 0, 3);
  TH2F *hCosp = new TH2F("hCosp", "", 300, -1, 1, 30, 0, 3);
  TH2F *hDCA = new TH2F("hDCA", "", 100, 0, 0.1, 30, 0, 3);
  TH2F *hDCAx = new TH2F("hDCAx", "", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hDCAy = new TH2F("hDCAy", "", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hDCAz = new TH2F("hDCAz", "", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hd0XYprod = new TH2F("hd0xyprod", "", 100, -0.01, 0.01, 30, 0, 3);
  TH2F *hd0XY1 = new TH2F("hd0xy1", "", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hd0XY2 = new TH2F("hd0xy2", "", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hMassVsPt = new TH2F("hMassVsPt", "", 200, 1.5, 2.5, 6, 0, 3);
  TH2F *hMassVsY = new TH2F("hMassVsY", "", 200, 1.5, 2.5, 10, 0, 5);
  TH2F *hMassReflVsPt = new TH2F("hMassReflVsPt", "", 200, 1.5, 2.5, 6, 0, 3);
  TH2F *hMassReflVsY = new TH2F("hMassReflVsY", "", 200, 1.5, 2.5, 10, 0, 5);
  
  TH2F *hResVx = new TH2F("hResVx", "", 200, -1000., 1000., 30, 0, 3);
  TH2F *hResVy = new TH2F("hResVy", "", 200, -1000., 1000., 30, 0, 3);
  TH2F *hResVz = new TH2F("hResVz", "", 200, -1000., 1000., 30, 0, 3);
  TH2F *hResPx = new TH2F("hResPx", "", 100, -1, 1, 30, 0, 3); //for Kaons
  TH2F *hResPy = new TH2F("hResPy", "", 100, -1, 1, 30, 0, 3);
  TH2F *hResPz = new TH2F("hResPz", "", 100, -1, 1, 30, 0, 3);
  TH2F *hResVxVsY = new TH2F("hResVxVsY", "", 200, -1000., 1000., 50, 0, 5);
  TH2F *hResVyVsY = new TH2F("hResVyVsY", "", 200, -1000., 1000., 50, 0, 5);
  TH2F *hResVzVsY = new TH2F("hResVzVsY", "", 200, -1000., 1000., 50, 0, 5);
  TH2F *hResPxVsY = new TH2F("hResPxVsY", "", 100, -1, 1, 50, 0, 5); //for Kaons
  TH2F *hResPyVsY = new TH2F("hResPyVsY", "", 100, -1, 1, 50, 0, 5);
  TH2F *hResPzVsY = new TH2F("hResPzVsY", "", 100, -1, 1, 50, 0, 5);
  TH2F *hResDist = new TH2F("hResDist", "", 100, -0.5, 0.5, 30, 0, 3);
  TH2F *hResDistXY = new TH2F("hResDistXY", "", 100, -0.1, 0.1, 30, 0, 3);
  TH1D *hNevents = new TH1D("hNevents", "", 1, 0, 1);//un singolo bin tra zero e 1 io fillo a 0.5 in questo modo sono a centro bin. Numero eventi

//histograms for checking simulation (quindi togli le righe di codice nel passaggio rivelatore)
	TH1D* hMassLambdatrue = new TH1D("hMassLambdatrue", "Mass", 400, 1, 1.2);//massa invariante con massa corretta
	TH1D* hMassLambdafalse = new TH1D("hMassLambdafalse", "Mass", 200, 0.5, 6);//massa invariante con massa switch
	TH1D* hVxGen = new TH1D("hVxGen", "Vx gen", 100, -24.5, 24.);
	TH1D* hVyGen = new TH1D("hVyGen", "Vy gen", 100, -24.5, 24.);
	TH1D* hVzGen = new TH1D("hVzGen", "Vz gen", 100, -24.5, 24.);
	TH1D* hsecVxGen = new TH1D("hsecVertGenX", "SecVerGenX protone", 300, -50., 50.);
	TH1D* hsecVyGen = new TH1D("hsecVertGenY", "SecVerGenY protone", 300, -50., 50.); 
	TH1D* hsecVzGen = new TH1D("hsecVertGenZ", "SecVerGenZ protone", 300, -50., 100.);
	TH1D* hTimeDist = new TH1D("hTimeDist", "Decay lenght c#tau", 500, -0.5, 50);//width=0.1
	TH2F* hctauPtotDist = new TH2F("hctauPtotDist", "c#tau vs p_{tot}", 500, 0. , 50.,1000, 0., 100);//width=0.1
	TH1D* hTime = new TH1D("hTime", "Decay time #tau", 500, -0.5, 50);//width=0.1
	TH1D* hLenght = new TH1D("hLenght", "L (no boost)", 500, -0.5, 50);//width=0.1
	TH2F* hLenghtPtotDist = new TH2F("hLenghtPtotDist", "L vs p_{tot}", 500, 0. , 50.,1000, 0., 100);//width=0.1
	TH2F* hLzVsPt = new TH2F("hLzVsPt", "Lz vs p_{t}", 500, 0 , 50.,1000, 0, 100);//width=0.1
	TH2F* hLzVsY = new TH2F("hLzVsY", "Lz vs Y", 200, 0 , 20,1000, 0, 100); //width=0.1
	TH2F* hLzVsYlay1 = new TH2F("hLzVsYlay1", "Lz vs Y first layer", 300, 0. , 50.,300, 0.5, 50);

//histograms for checking reconstruction (quindi tieni le righe di codice nel passaggio rivelatore)
	TH1D* hPrecgenXprot = new TH1D("hPrecgenXprot", "P_{Xrec}-P_{Xgen} proton", 300, -1, 1);//width=0.01
	TH1D* hPrecgenYprot = new TH1D("hPrecgenYprot", "P_{Yrec}-P_{Ygen} proton", 300, -1, 1);//width=0.01
	TH1D* hPrecgenZprot = new TH1D("hPrecgenZprot", "P_{Zrec}-P_{Zgen} proton", 300, -1, 1);//width=0.01
	TH1D* hPrecgenXpion = new TH1D("hPrecgenXpion", "P_{Xrec}-P_{Xgen} pion", 300, -1, 1);//width=0.01
	TH1D* hPrecgenYpion = new TH1D("hPrecgenYpion", "P_{Yrec}-P_{Ygen} pion", 300, -1, 1);//width=0.01
	TH1D* hPrecgenZpion = new TH1D("hPrecgenZpion", "P_{Zrec}-P_{Zgen} pion", 300, -1, 1);//width=0.01
	TH2F* hResPxVsLzprot = new TH2F("hResPxVsLzprot", "P_{Xrec}-P_{Xgen}  vs Lz proton",10000, 0., 100, 300, -1 , 1);//width=0.01
	TH2F* hResPyVsLzprot = new TH2F("hResPyVsLzprot", "P_{Yrec}-P_{Ygen}  vs Lz proton",10000, 0., 100, 300, -1 , 1);
	TH2F* hResPzVsLzprot = new TH2F("hResPzVsLzprot", "P_{Zrec}-P_{Zgen}  vs Lz proton",10000, 0., 100, 300, -1 , 1);
	TH2F* hResPxVsLzpion = new TH2F("hResPxVsLzpion", "P_{Xrec}-P_{Xgen}  vs Lz pion",10000, 0., 100, 300, -1 , 1);//width=0.01
	TH2F* hResPyVsLzpion = new TH2F("hResPyVsLzpion", "P_{Yrec}-P_{Ygen}  vs Lz pion",10000, 0., 100, 300, -1 , 1);
	TH2F* hResPzVsLzpion = new TH2F("hResPzVsLzpion", "P_{Zrec}-P_{Zgen}  vs Lz pion",10000, 0., 100, 300, -1 , 1);
	TH2F* hResPxVsYprot = new TH2F("hResPxVsYprot", "P_{Xrec}-P_{Xgen}  vs Y proton",10000, 0., 100, 300, -1 , 1);//width=0.01
	TH2F* hResPyVsYprot = new TH2F("hResPyVsYprot", "P_{Yrec}-P_{Ygen}  vs Y proton",10000, 0., 100, 300, -1 , 1);
	TH2F* hResPzVsYprot = new TH2F("hResPzVsYprot", "P_{Zrec}-P_{Zgen}  vs Y proton",10000, 0., 100, 300, -1 , 1);
	TH2F* hResPxVsYpion = new TH2F("hResPxVsYpion", "P_{Xrec}-P_{Xgen}  vs Y pion",10000, 0., 100, 300, -1 , 1);//width=0.01
	TH2F* hResPyVsYpion = new TH2F("hResPyVsYpion", "P_{Yrec}-P_{Ygen}  vs Y pion",10000, 0., 100, 300, -1 , 1);
	TH2F* hResPzVsYpion = new TH2F("hResPzVsYpion", "P_{Zrec}-P_{Zgen}  vs Y pion",10000, 0., 100, 300, -1 , 1);
//correction Py
	TH2F* hResPyVsLzprotcorr = new TH2F("hResPyVsLzprotcorr", "P_{Yrec}-P_{Ygen}  vs Lz proton",10000, 0., 100, 300, -1 , 1);
	TH2F* hResPyVsLzpioncorr = new TH2F("hResPyVsLzpioncorr", "P_{Yrec}-P_{Ygen}  vs Lz pion",10000, 0., 100, 300, -1 , 1);
	TH2F* hResPyVsYprotcorr = new TH2F("hResPyVsYprotcorr", "P_{Yrec}-P_{Ygen}  vs Y proton",10000, 0., 100, 300, -1 , 1);
	TH2F* hResPyVsYpioncorr = new TH2F("hResPyVsYpioncorr", "P_{Yrec}-P_{Ygen}  vs Y pion",10000, 0., 100, 300, -1 , 1);
//Nfake
	TH2F* hnfakeVsLzprot = new TH2F("hnfakeVsLzprot", "L_{z} Proton vs Nfake ", 6, -0.5 , 5.5,500, 0, 50);//width=0.1Lz
	TH2F* hResPyVsnfakeprot= new TH2F("hResPyVsnfakeprot", " P_{Yrec}-P_{Ygen} vs Nfake Proton", 6, -0.5 , 5.5,300, -1 , 1);//width=0.01Lz
	TH1F* hResPycorrection= new TH1F("hResPycorrection", " P_{Ynotcorrect}-P_{Ycorrect} ",300, -1 , 1);

//invariant mass reconstruction
	TH1D* hMassreconstr = new TH1D("hMassreconstr", "Mass", 3000, 0,3);
	TH2F* hMassinrecVsNfake = new TH2F("hMassinrecVsNfake", "Nfake vs Mass", 6, -0.5 , 5.5,3000, 0,3);
	TH2F* hMassinrecVszP = new TH2F("hMassinrecVszP", "Mass vs zP", 500, 0, 50,3000, 0,3);
  
	TH2F *hd0 = new TH2F("hd0", "", 100, 0, 0.1, 30, 0, 3);
  THnSparseF *hsp = CreateSparse();
  const Int_t nDim=static_cast<const Int_t>(hsp->GetNdimensions());
  Double_t arrsp[nDim];

  TFile *fnt = 0x0;//INDIRIZZO??
  TNtuple *ntD0cand = 0x0;
  if (writeNtuple)//????
    {
      fnt = new TFile("fntSig.root", "recreate");
      ntD0cand = new TNtuple("ntD0cand", "ntD0cand", "mass:pt:y:dist:cosp:d01:d02:d0prod:dca:ptMin:ptMax", 32000);
    }
  Float_t arrnt[11];
	float nfaketrk = 0, nfaketrkprot = 0, nfaketrkpion = 0;
	
//GENERAZIONE
  for (Int_t iev = 0; iev < nevents; iev++){
    hNevents->Fill(0.5);
    Double_t vprim[3] = {0, 0, 0};//vertice primario
    if(iev%1000==0) printf(" ***************  ev = %d \n", iev);
    int nrec = 0;
    int nfake = 0;
		nfaketrk = 0;
		nfaketrkprot = 0;
		nfaketrkpion = 0;
    double pxyz[3], pProtRec[3], pPionRec[3], pProtGen[3], pPionGen[3];
		double y = 0;
    
    if (simulateBg && (iev%refreshBg)==0) det->GenBgEvent(0.,0.,0.);//bkg c'è sempre
    Double_t ptGenD = hLambdapt->GetRandom(); // get Lambda distribution from file
    Double_t yGenD = hLambday->GetRandom();
    Double_t phi = gRandom->Rndm() * 2 * TMath::Pi();
    Double_t pxGenD = ptGenD * TMath::Cos(phi);
    Double_t pyGenD = ptGenD * TMath::Sin(phi);
    
    Double_t mass = TDatabasePDG::Instance()->GetParticle(pdgParticle)->Mass();
    Double_t mt = TMath::Sqrt(ptGenD * ptGenD + mass * mass);
    Double_t pzGenD = mt * TMath::SinH(yGenD);
    Double_t en = mt * TMath::CosH(yGenD);

    Double_t massProton = TDatabasePDG::Instance()->GetParticle(2212)->Mass();
    Double_t massPion = TDatabasePDG::Instance()->GetParticle(-211)->Mass();
	//printf("massproton= %f, masspion= %f", massProton, massPion);
    mom->SetPxPyPzE(pxGenD, pyGenD, pzGenD, en); //setto i valori del quadrimpulso della particella lambda generata
		//printf("Momenti gen: px= %f, py= %f, pz= %f \n", mom->Px(), mom->Py(), mom->Pz());
		hVxGen->Fill(vprim[0]);
		hVyGen->Fill(vprim[1]);
		hVzGen->Fill(vprim[2]);
		//printf("Vertici= %f %f %f \n",vprim[0],vprim[1],vprim[2]);
//DEFINISCO IL NUMERO DI PRODOTTI
    Int_t np;
    do{
      fDecayer->Decay(pdgParticle, mom);
			//Decay a particle----input: pdg code and momentum of the particle to be decayed, all informations about decay products are stored in fEvtstdhep(non so cosa sia) 
      np = fDecayer->ImportParticles(particles);
			//fDecayer->ImportParticles() -----Input: pointer to a TClonesArray - Output(Int_t): number of decay products. Put all the informations about the decay products in the TClonesArray particles
    } while (np < 0);
    
		double diff;
    Int_t arrpdgdau[2];
    Double_t ptK = -999.;//CAMBIA VARIABILI PRODOTTI
    Double_t ptPi = -999.;
    Double_t yK=-999.;
    Double_t yPi = -999.;
    Int_t icount = 0;
    Double_t secvertgenK[3]={0.,0.,0.};
    Double_t secvertgenPi[3]={0.,0.,0.};
		Double_t ptot = 0, massal = 0, p = 0, Respxyzprot[3] = {0.,0.,0.}, Respxyzpion[3] = {0.,0.,0.}, PyCorrectProt = 0, PyCorrectPion = 0;
 //printf("evento= %d, numero prodotti:%d \n",iev,np);  
    // loop on decay products
    for (int i = 0; i < np; i++) { 
      TParticle *iparticle1 = (TParticle *)particles->At(i);
      Int_t kf = TMath::Abs(iparticle1->GetPdgCode());
      vX = iparticle1->Vx();
      vY = iparticle1->Vy();
      vZ = iparticle1->Vz();
			
      if (kf == pdgParticle){				
			//Lambda particle

				hYGen->Fill(iparticle1->Y());
				hPtGen->Fill(iparticle1->Pt());
				hYPtGen->Fill(iparticle1->Y(), iparticle1->Pt());
				ptot=iparticle1->P();
				massal=iparticle1->GetMass();
				p=iparticle1->Pt();
				y=iparticle1->Y();
				//printf("mother part = %d, code=%d, pt=%f, y=%f \n", i, kf, iparticle1->Pt(), iparticle1->Y());
      }	
      if (kf == 2212 || kf == 211){
//devo mettere il segno positivo perchè tanto ho fatto abs e perchè na60plus non riconosce carica
	// daughters that can be reconstructed: p and pi
			// printf("protone o pione \n");
				Double_t e = iparticle1->Energy();
				Double_t px = iparticle1->Px();
				Double_t py = iparticle1->Py();
				Double_t pz = iparticle1->Pz();	    
				TLorentzVector *pDecDau = new TLorentzVector(0., 0., 0., 0.);//quadrimpulso figli
				pDecDau->SetXYZM(iparticle1->Px(), iparticle1->Py(), iparticle1->Pz(), iparticle1->GetMass());
				//printf("pdg code= %d \n",kf);
				Int_t crg=1;//setta il segno della carica
				if(iparticle1->GetPdgCode()<0) crg=-1;//se il valore del codice è negativo allora la carica è negativa
		
				//togli per fare controlli iniziali	di massa inv e vertici
		/*		if (!det->SolveSingleTrack(pDecDau->Pt(), pDecDau->Rapidity(), pDecDau->Phi(), iparticle1->GetMass(), crg, vX, vY, vZ, 0, 1, 99)) continue;
				KMCProbeFwd *trw = det->GetLayer(0)->GetWinnerMCTrack();
				if (!trw) continue;
				if (trw->GetNormChi2(kTRUE) > ChiTot) continue;

				nrec++;
				nfake += trw->GetNFakeITSHits();
				trw->GetPXYZ(pxyz);
				nfaketrk = trw->GetNFakeITSHits();*/

				//if(iev%100==0) printf("SENZA correzione P = %f %f %f \n",pxyz[0],pxyz[1],pxyz[2] );

				if (kf == 2212){
					//printf("PROTONE \n");
					// Proton daughter
				//	nfaketrkprot = trw->GetNFakeITSHits();
					ptK = iparticle1->Pt();
					yK = iparticle1->Y();
					hptP->Fill(ptGenD,ptK);
					hyP->Fill(iparticle1->Y());
					secvertgenK[0] = iparticle1->Vx();
					secvertgenK[1] = iparticle1->Vy();
					secvertgenK[2] = iparticle1->Vz();
					daugen[0].SetXYZM(iparticle1->Px(), iparticle1->Py(), iparticle1->Pz(), iparticle1->GetMass());
					//printf("impulsi= %f ; %f ; %f ; massa %f \n",daugen[0].X(),daugen[0].Y(), daugen[0].Z(),daugen[0].M());
					//daurec[0].SetXYZM(pxyz[0], pxyz[1], pxyz[2], iparticle1->GetMass());
					//printf("pxrec %f pxgen %f \n", pxyz[0],iparticle1->Px());
					pProtGen[0]=iparticle1->Px();
					pProtGen[1]=iparticle1->Py();
					pProtGen[2]=iparticle1->Pz();
		/*			hPrecgenXprot->Fill(pxyz[0]-iparticle1->Px());
					hPrecgenYprot->Fill(pxyz[1]-iparticle1->Py());
					hPrecgenZprot->Fill(pxyz[2]-iparticle1->Pz());
					Respxyzprot[0]=(pxyz[0]-iparticle1->Px());
					Respxyzprot[1]=(pxyz[1]-iparticle1->Py());
					Respxyzprot[2]=(pxyz[2]-iparticle1->Pz());
					hResPxVsLzprot->Fill(secvertgenK[2],Respxyzprot[0]);
					hResPyVsLzprot->Fill(secvertgenK[2],Respxyzprot[1]);
					hResPzVsLzprot->Fill(secvertgenK[2],Respxyzprot[2]);
					hnfakeVsLzprot->Fill(nfaketrk,secvertgenK[2]);*/
					//hResPyVsnfakeprot->Fill(nfaketrk,Respxyzprot[1]);
			/*		diff=pxyz[1];
					trw->PropagateToZBxByBz(secvertgenK[2]);
					trw->GetPXYZ(pxyz);
					daurec[0].SetXYZM(pxyz[0], pxyz[1], pxyz[2], iparticle1->GetMass());
					//if(iev%100==0) printf(" correzione P = %f %f %f \n",pxyz[0],pxyz[1],pxyz[2] );
					printf("p %f %f %f \n",pxyz[0],pxyz[1],pxyz[2]);
					hResPycorrection->Fill(diff-pxyz[1]);
					PyCorrectProt=pxyz[1]-iparticle1->Py();
					hResPyVsLzprotcorr->Fill(secvertgenK[2],PyCorrectProt);
					hResPyVsnfakeprot->Fill(nfaketrk,PyCorrectProt);
					daurecswapmass[0].SetXYZM(pxyz[0], pxyz[1], pxyz[2],massPion);
					recProbe[0] = *trw; */

				}else if (kf == 211){
					// Pion daughter
					//printf("PIONEEEEEE \n");
			//		nfaketrkpion = trw->GetNFakeITSHits();
					ptPi = iparticle1->Pt();
					yPi = iparticle1->Y();
					hptPi->Fill(ptGenD,ptPi);
					hyPi->Fill(iparticle1->Y());
					secvertgenPi[0] = iparticle1->Vx();
					secvertgenPi[1] = iparticle1->Vy();
					secvertgenPi[2] = iparticle1->Vz();
					daugen[1].SetXYZM(iparticle1->Px(), iparticle1->Py(), iparticle1->Pz(), iparticle1->GetMass());
					//printf("impulsi= %f ; %f ; %f ; massa %f \n",daugen[1].X(),daugen[1].Y(), daugen[1].Z(),daugen[1].M());
				//	daurec[1].SetXYZM(pxyz[0], pxyz[1], pxyz[2],  iparticle1->GetMass());
					pPionGen[0]=iparticle1->Px();
					pPionGen[1]=iparticle1->Py();
					pPionGen[2]=iparticle1->Pz();
				/*	hPrecgenXpion->Fill(pxyz[0]-iparticle1->Px());
					hPrecgenYpion->Fill(pxyz[1]-iparticle1->Py());
					hPrecgenZpion->Fill(pxyz[2]-iparticle1->Pz());
					Respxyzpion[0]=(pxyz[0]-iparticle1->Px());
					Respxyzpion[1]=(pxyz[1]-iparticle1->Py());
					Respxyzpion[2]=(pxyz[2]-iparticle1->Pz());
					hResPxVsLzpion->Fill(secvertgenPi[2],Respxyzpion[0]);
					hResPyVsLzpion->Fill(secvertgenPi[2],Respxyzpion[1]);
					hResPzVsLzpion->Fill(secvertgenPi[2],Respxyzpion[2]);
					trw->PropagateToZBxByBz(secvertgenPi[2]);
					trw->GetPXYZ(pxyz);
					PyCorrectPion=pxyz[1]-iparticle1->Py();
					hResPyVsLzpioncorr->Fill(secvertgenPi[2],PyCorrectPion);
					daurecswapmass[1].SetXYZM(pxyz[0], pxyz[1], pxyz[2],massProton);
					recProbe[1] = *trw;*/
					}
				}
			}
		
	//Check Lambda invariant mass 
		double E1= sqrt(daugen[0].X() * daugen[0].X() + daugen[0].Y() * daugen[0].Y() + daugen[0].Z() * daugen[0].Z() + daugen[0].M() * daugen[0].M());
		double E2= sqrt(daugen[1].X() * daugen[1].X() + daugen[1].Y() * daugen[1].Y() + daugen[1].Z() * daugen[1].Z() + daugen[1].M() * daugen[1].M());
		double p1= sqrt(daugen[0].X() * daugen[0].X() + daugen[0].Y() * daugen[0].Y() + daugen[0].Z() * daugen[0].Z());
		double p2= sqrt(daugen[1].X() * daugen[1].X() + daugen[1].Y() * daugen[1].Y() + daugen[1].Z() * daugen[1].Z());
		double massainv= sqrt(daugen[0].M() * daugen[0].M() + daugen[1].M() * daugen[1].M() + 2*(E1*E2 - daugen[0].X() * daugen[1].X() - daugen[0].Y() * daugen[1].Y() - daugen[0].Z() * daugen[1].Z()));

	//Check Lambda invariant mass, switching masses of proton and pion
		double E1false= sqrt(daugen[0].X() * daugen[0].X() + daugen[0].Y() * daugen[0].Y() + daugen[0].Z() * daugen[0].Z() + daugen[1].M() * daugen[1].M());
		double E2false= sqrt(daugen[1].X() * daugen[1].X() + daugen[1].Y() * daugen[1].Y() + daugen[1].Z() * daugen[1].Z() + daugen[0].M() * daugen[0].M());
		double massainvfal=sqrt(daugen[1].M() * daugen[1].M() + daugen[0].M() * daugen[0].M() + 2*(E1false*E2false - daugen[0].X() * daugen[1].X() - daugen[0].Y() * daugen[1].Y() - daugen[0].Z() * daugen[1].Z()));

	//Decay length evaluation ct=L/gammabeta
		double L = sqrt(secvertgenK[0]*secvertgenK[0] + secvertgenK[1]*secvertgenK[1] + secvertgenK[2]*secvertgenK[2]);
		double gammabeta=ptot/massal;
		double ctau=L/(gammabeta);
		double tau=ctau/3.;
		//printf("momento totale %f \n",ptot);
		//printf("Decay Lenght=%f \n",ctau);
		
		//printf("massa invar=%f \n", massainv);
    if (ptK > 0 && ptPi > 0) hyPiP->Fill(yPi, yK);
//vertice gen e massa invariante 
	hMassLambdatrue->Fill(massainv);
		hMassLambdafalse->Fill(massainvfal);
		hsecVxGen->Fill(secvertgenK[0]);
		hsecVyGen->Fill(secvertgenK[1]);
		hsecVzGen->Fill(secvertgenK[2]);
		hTimeDist->Fill(ctau);
		hTime->Fill(tau);
		hLenght->Fill(L);
		hctauPtotDist->Fill(ptot,ctau);
		hLenghtPtotDist->Fill(ptot,L);
		hLzVsPt->Fill(p,secvertgenK[2]);
		hLzVsY->Fill(y,secvertgenK[2]);

    //if (nrec < 2) continue; //togli per fare controlli iniziali
	 	hResPxVsYprot->Fill(y,Respxyzprot[0]);
		hResPyVsYprot->Fill(y,Respxyzprot[1]);
		hResPzVsYprot->Fill(y,Respxyzprot[2]);
		hResPxVsYpion->Fill(y,Respxyzpion[0]);
		hResPyVsYpion->Fill(y,Respxyzpion[1]);
		hResPzVsYpion->Fill(y,Respxyzpion[2]);
		
		hResPyVsYprotcorr->Fill(y,PyCorrectProt);
		hResPyVsYpioncorr->Fill(y,PyCorrectPion);
    recProbe[0].PropagateToDCA(&recProbe[1]);
    
    parent = daurec[0];
    parent += daurec[1];// equivalent to parent=parent+daurec[]

    parentgen = daugen[0];
    parentgen += daugen[1];
    parentrefl = daurecswapmass[0];
    parentrefl += daurecswapmass[1];
      
    Double_t  ptRecD=parent.Pt();
    Double_t  massRecD=parent.M();

    Double_t  massRecReflD=parentrefl.M();
    Double_t yRecD = 0.5 * TMath::Log((parent.E() + parent.Pz()) / (parent.E() - parent.Pz()));
    hYPtRecoAll->Fill(yRecD, ptRecD);
    hPtRecoAll->Fill(ptRecD);
    hPtGenRecoAll->Fill(ptGenD);
    hPtRecoVsGenAll->Fill(ptGenD,ptRecD);
    hDiffPtRecoGenAll->Fill(ptGenD,(ptRecD-ptGenD));
    hYRecoAll->Fill(yRecD);
    hYGenRecoAll->Fill(yGenD);
    hMassAll->Fill(massRecD);
    hMassRefl->Fill(massRecReflD);
    if (nfake > 0){
      hYPtRecoFake->Fill(yRecD, ptRecD);
      hPtRecoFake->Fill(ptRecD);
      hMassFake->Fill(massRecD);
    }
    hMassVsPt->Fill(massRecD,ptRecD);
    hMassVsY->Fill(massRecD,yRecD);
    hMassReflVsPt->Fill(massRecReflD,ptRecD);
    hMassReflVsY->Fill(massRecReflD,yRecD);
   
    Float_t d1 = recProbe[1].GetX() - recProbe[0].GetX();
    Float_t d2 = recProbe[1].GetY() - recProbe[0].GetY();
    Float_t d3 = recProbe[1].GetZ() - recProbe[0].GetZ();
    Float_t dca = sqrt(d1 * d1 + d2 * d2 + d3 * d3);
    // printf(" DCA = %f\n", sqrt(d1 * d1 + d2 * d2 + d3 * d3));
    hDCA->Fill(dca, ptRecD);
    hDCAx->Fill(d1, ptRecD);
    hDCAy->Fill(d2, ptRecD);
    hDCAz->Fill(d3, ptRecD);
    
    //      Double_t xP = (recProbe[1].GetX() + recProbe[0].GetX()) / 2.;
    //      Double_t yP = (recProbe[1].GetY() + recProbe[0].GetY()) / 2.;
    //      Double_t zP = (recProbe[1].GetZ() + recProbe[0].GetZ()) / 2.;
    
    Double_t xP, yP, zP;
    ComputeVertex(recProbe[0],recProbe[1],xP,yP,zP);
    Double_t residVx=10000.*(xP - secvertgenK[0]);
    Double_t residVy=10000.*(yP - secvertgenK[1]);
    Double_t residVz=10000.*(zP - secvertgenK[2]);
    hResVx->Fill(residVx, ptRecD);
    hResVy->Fill(residVy, ptRecD);
    hResVz->Fill(residVz, ptRecD);
    hResVxVsY->Fill(residVx, yRecD);
    hResVyVsY->Fill(residVy, yRecD);
    hResVzVsY->Fill(residVz, yRecD);

	//PROPAGO Zrec protone
/*		printf("pippoo \n");
		recProbe[0].PropagateToZBxByBz(zP);
		recProbe[0].GetPXYZ(pxyz);
		printf("p REC PROTONE %f %f %f \n",pxyz[0],pxyz[1],pxyz[2]);
		pProtRec[0]=pxyz[0];
		pProtRec[1]=pxyz[1];
		pProtRec[2]=pxyz[2];

	//PROPAGO Zrec pione
		recProbe[1].PropagateToZBxByBz(zP);
		recProbe[1].GetPXYZ(pxyz);
		printf("p REC PIONE %f %f %f \n",pxyz[0],pxyz[1],pxyz[2]);
		pPionRec[0]=pxyz[0];
		pPionRec[1]=pxyz[1];
		pPionRec[2]=pxyz[2];

	//Calcolo massa invariante ricostruita
		double E1rec= sqrt(pProtRec[0] * pProtRec[0] + pProtRec[1] * pProtRec[1] + pProtRec[2] * pProtRec[2] + daugen[0].M() * daugen[0].M());
		double E2rec= sqrt(pPionRec[0] * pPionRec[0] + pPionRec[1] * pPionRec[1] + pPionRec[2] * pPionRec[2] + daugen[1].M() * daugen[1].M());
		double p1rec= sqrt(pProtRec[0] * pProtRec[0] + pProtRec[1] * pProtRec[1] + pProtRec[2] * pProtRec[2]);
		double p2rec= sqrt(pPionRec[0] * pPionRec[0] + pPionRec[1] * pPionRec[1] + pPionRec[2] * pPionRec[2]);
		double massainvrec= sqrt(daugen[0].M() * daugen[0].M() + daugen[1].M() * daugen[1].M() + 2*(E1rec*E2rec - pProtRec[0] * pPionRec[0]- pProtRec[1] * pPionRec[1] - pProtRec[2] * pPionRec[2]));
		
		hMassreconstr->Fill(massainvrec);
		hMassinrecVsNfake->Fill(nfaketrk,massainvrec);
		hMassinrecVszP->Fill(zP,massainvrec);
		variable[0]=nfaketrkprot;
		variable[1]=nfaketrkpion;
		variable[2]=(float)y;
		variable[3]=(float)xP;
		variable[4]=(float)yP;
		variable[5]=(float)zP;
		variable[6]=(float)secvertgenK[0];
		variable[7]=(float)secvertgenK[1];
		variable[8]=(float)secvertgenK[2];
		variable[9]=(float)massainvrec;
		variable[10]=(float)(pProtRec[0]-pProtGen[0]);
		variable[11]=(float)(pProtRec[1]-pProtGen[1]);
		variable[12]=(float)(pProtRec[2]-pProtGen[2]);
		variable[13]=(float)(pPionRec[0]-pPionGen[0]);
		variable[14]=(float)(pPionRec[1]-pPionGen[1]);
		variable[15]=(float)(pPionRec[2]-pPionGen[2]);
	//Riempimento Ntupla
		//nt->Fill(y,xP,yP,zP,secvertgenK[0],secvertgenK[1],secvertgenK[2],massainvrec,pProtRec[0]-pProtGen[0],pProtRec[1]-pProtGen[1],pProtRec[2]-pProtGen[2],pPionRec[0]-pPionGen[0],pPionRec[1]-pPionGen[1],pPionRec[2]-pPionGen[2]);
		nt->Fill(variable);
		treeNfake->Fill();
    */
    hResPx->Fill(daurec[0].Px() - daugen[0].Px(), ptRecD);
    hResPy->Fill(daurec[0].Py() - daugen[0].Py(), ptRecD);
    hResPz->Fill(daurec[0].Pz() - daugen[0].Pz(), ptRecD);
    hResPx->Fill(daurec[1].Px() - daugen[1].Px(), ptRecD);
    hResPy->Fill(daurec[1].Py() - daugen[1].Py(), ptRecD);
    hResPz->Fill(daurec[1].Pz() - daugen[1].Pz(), ptRecD);

    hResPxVsY->Fill(daurec[0].Px() - daugen[0].Px(), yRecD);
    hResPyVsY->Fill(daurec[0].Py() - daugen[0].Py(), yRecD);
    hResPzVsY->Fill(daurec[0].Pz() - daugen[0].Pz(), yRecD);
    hResPxVsY->Fill(daurec[1].Px() - daugen[1].Px(), yRecD);
    hResPyVsY->Fill(daurec[1].Py() - daugen[1].Py(), yRecD);
    hResPzVsY->Fill(daurec[1].Pz() - daugen[1].Pz(), yRecD);
    
    // cout << "secvert generated Pion: " << secvertgenPi[0] << "  " << secvertgenPi[1] << "  " << secvertgenPi[2] << endl;
    // cout << "Reco Vert  Pion: " << xP << "  " << yP << "  " << zP << endl;
    
    Float_t dist = TMath::Sqrt(xP * xP + yP * yP + zP * zP);
    Float_t distXY = TMath::Sqrt(xP * xP + yP * yP);
    Float_t distgen = TMath::Sqrt(secvertgenPi[0] * secvertgenPi[0] + secvertgenPi[1] * secvertgenPi[1] + secvertgenPi[2] * secvertgenPi[2]);
    Float_t distgenXY = TMath::Sqrt(secvertgenPi[0] * secvertgenPi[0] + secvertgenPi[1] * secvertgenPi[1]);
    // printf("dist = %f , distXY=%f , dx=%f, dy=%f, dz=%f z1=%f, z2=%f \n", dist, distXY, xP, yP, zP, recProbe[0].GetZ(), recProbe[1].GetZ());
    // printf("distgen = %f , distgenXY=%f \n", distgen, distgenXY);
    
    Double_t vsec[3] = {xP, yP, zP};
    Double_t cosp = CosPointingAngle(vprim, vsec, parent);
    Double_t cts = CosThetaStar(parent,daurec[0]);
    Double_t ipD = ImpParXY(vprim, vsec, parent);
    hCosp->Fill(cosp, ptRecD);
    // printf(" ***** ***** cos point = %f \n", cosp);
    //if (cosp < -0.98)
    //    printf("SMALL COSPOINT");
    
    hResDist->Fill(dist - distgen, ptRecD);
    hResDistXY->Fill(distXY - distgenXY, ptRecD);
    
    //recProbe[0].PropagateToDCA(&recProbe[1]);
    
    // hYPtAll->Fill(parent.Y(), ptRecD);
    // hPtAll->Fill(ptRecD);
    hDistXY->Fill(distXY, ptRecD);
    hDist->Fill(dist, ptRecD);
    hDistgenXY->Fill(distgenXY, ptRecD);
    hDistgen->Fill(distgen, ptRecD);
      
    //AliExternalTrackParam *track1 = (AliExternalTrackParam *)recProbe[0].GetTrack();
    //AliExternalTrackParam *track2 = (AliExternalTrackParam *)recProbe[1].GetTrack();
    recProbe[0].PropagateToZBxByBz(0);
    Double_t d0x1 = recProbe[0].GetX();
    Double_t d0y1 = recProbe[0].GetY();
    Double_t d0xy1 = TMath::Sqrt(d0x1 * d0x1 + d0y1 * d0y1);
    if (d0x1 < 0)
      d0xy1 *= -1;
    
    recProbe[1].PropagateToZBxByBz(0);
    Double_t d0x2 = recProbe[1].GetX();
    Double_t d0y2 = recProbe[1].GetY();
    Double_t d0xy2 = TMath::Sqrt(d0x2 * d0x2 + d0y2 * d0y2);
    if (d0x2 < 0)
      d0xy2 *= -1;
    
    // printf("d0xy1 = %f, d0xy2 = %f \n", d0xy1, d0xy2);
    
    hd0XYprod->Fill(d0xy1 * d0xy2, ptRecD);
    hd0XY1->Fill(d0xy1, ptRecD);
    hd0XY2->Fill(d0xy2, ptRecD);
      
    arrsp[0] = massRecD;
    arrsp[1] = ptRecD;
    arrsp[2] = yRecD;
    arrsp[3] = dist;
    arrsp[4] = cosp;
    arrsp[5] = TMath::Min(TMath::Abs(d0xy1),TMath::Abs(d0xy2));
    arrsp[6] = d0xy1 * d0xy2;
    arrsp[7] = dca;
    arrsp[8] = TMath::Min(recProbe[0].GetTrack()->Pt(),recProbe[1].GetTrack()->Pt());
    arrsp[9] = TMath::Abs(ipD);	    
    arrsp[10] = cts;      
    hsp->Fill(arrsp);
    
    if (ntD0cand){
      arrnt[0] = massRecD;
      arrnt[1] = ptRecD;
      arrnt[2] = yRecD;
      arrnt[3] = dist;
      arrnt[4] = cosp;
      arrnt[5] = d0xy1;
      arrnt[6] = d0xy2;
      arrnt[7] = d0xy1 * d0xy2;
      arrnt[8] = dca;
      arrnt[9] = TMath::Min(recProbe[0].GetTrack()->Pt(),recProbe[1].GetTrack()->Pt());
      arrnt[10] = TMath::Max(recProbe[0].GetTrack()->Pt(),recProbe[1].GetTrack()->Pt());
      ntD0cand->Fill(arrnt);
    }
  } //event loop
  
  hMassAll->SetLineColor(kBlue);
  hMassAll->Draw();
  hMassAll->SetMinimum(0.1);
  hMassFake->SetLineColor(kRed);
  hMassFake->Draw("same");


//CI SONO ANCORA LE PROIEZIONI E I GRAFICI FATTI PRIMA DI FARLI NELLA MACRO PLOTLAMBDA.C
 //CHECKING SIMULATION
//Projection L vs Ptot
	TH1D* hProjLP0 = hLenghtPtotDist->ProjectionY("L vs Ptot range [0 0.5] GeV/c",hLenghtPtotDist->GetXaxis()->FindBin(0.01),hLenghtPtotDist->GetXaxis()->FindBin(0.49));
	TH1D* hProjLP1 = hLenghtPtotDist->ProjectionY("L vs Ptot range [0.5 1] GeV/c",hLenghtPtotDist->GetXaxis()->FindBin(0.51),hLenghtPtotDist->GetXaxis()->FindBin(0.99));
	TH1D* hProjLP2 = hLenghtPtotDist->ProjectionY("L vs Ptot range [1 1.5] GeV/c",hLenghtPtotDist->GetXaxis()->FindBin(1.01),hLenghtPtotDist->GetXaxis()->FindBin(1.49));
	printf("bincenterPTOTDIST %f,bin %d \n",hLenghtPtotDist->GetXaxis()->GetBinCenter((hLenghtPtotDist->GetXaxis()->FindBin(0.01))),hLenghtPtotDist->GetXaxis()->FindBin(0.01));
	printf("bincenterPTOTDIST %f,bin %d \n",hLenghtPtotDist->GetXaxis()->GetBinCenter((hLenghtPtotDist->GetXaxis()->FindBin(0.49))),hLenghtPtotDist->GetXaxis()->FindBin(0.49));
	printf("bincenterPTOTDIST %f,bin %d \n",hLenghtPtotDist->GetXaxis()->GetBinCenter((hLenghtPtotDist->GetXaxis()->FindBin(0.51))),hLenghtPtotDist->GetXaxis()->FindBin(0.51));
//Projection ctau vs Ptot
	TH1D* hProjctauP0 = hctauPtotDist->ProjectionY("c#tau P range [0 0.5] GeV/c",hctauPtotDist->GetXaxis()->FindBin(0.01),hctauPtotDist->GetXaxis()->FindBin(0.49));
	TH1D* hProjctauP1 = hctauPtotDist->ProjectionY("c#tau P range [0.5 1] GeV/c",hctauPtotDist->GetXaxis()->FindBin(0.51),hctauPtotDist->GetXaxis()->FindBin(0.99));
	TH1D* hProjctauP2 = hctauPtotDist->ProjectionY("c#tau P range [1 1.5] GeV/c",hctauPtotDist->GetXaxis()->FindBin(1.01),hctauPtotDist->GetXaxis()->FindBin(1.49));
//Projection Lz vs Pt
	TH1D* hProjLzPt0 = hLzVsPt->ProjectionY("Lz vs Pt range [0 0.5] GeV/c",hLzVsPt->GetXaxis()->FindBin(0.01),hLzVsPt->GetXaxis()->FindBin(0.49));
	TH1D* hProjLzPt1 = hLzVsPt->ProjectionY("Lz vs Pt range [0.5 1] GeV/c",hLzVsPt->GetXaxis()->FindBin(0.51),hLzVsPt->GetXaxis()->FindBin(0.99));
	TH1D* hProjLzPt2 = hLzVsPt->ProjectionY("Lz vs Pt range [1 1.5] GeV/c",hLzVsPt->GetXaxis()->FindBin(1.01),hLzVsPt->GetXaxis()->FindBin(1.49));
	TH1D* hProjLzPt3 = hLzVsPt->ProjectionY("Lz vs Pt range [1.5 2] GeV/c",hLzVsPt->GetXaxis()->FindBin(1.51),hLzVsPt->GetXaxis()->FindBin(1.99));
//projection Lz vs Y	
	TH1D* hProjLzY0 = hLzVsY->ProjectionY("Lz vs Y range [1.5 2] ",hLzVsY->GetXaxis()->FindBin(1.51),hLzVsY->GetXaxis()->FindBin(1.99));
	TH1D* hProjLzY1 = hLzVsY->ProjectionY("Lz vs Y range [2 2.5] ",hLzVsY->GetXaxis()->FindBin(2.01),hLzVsY->GetXaxis()->FindBin(2.49));
	TH1D* hProjLzY2 = hLzVsY->ProjectionY("Lz vs Y range [2.5 3] ",hLzVsY->GetXaxis()->FindBin(2.51),hLzVsY->GetXaxis()->FindBin(2.99));
	TH1D* hProjLzY3 = hLzVsY->ProjectionY("Lz vs Y range [3 3.5] ",hLzVsY->GetXaxis()->FindBin(3.01),hLzVsY->GetXaxis()->FindBin(3.49));
	TH1D* hProjLzY4 = hLzVsY->ProjectionY("Lz vs Y range [3.5 4] ",hLzVsY->GetXaxis()->FindBin(3.51),hLzVsY->GetXaxis()->FindBin(3.99));
	TH1D* hProjLzY5 = hLzVsY->ProjectionY("Lz vs Y range [4 4.5] ",hLzVsY->GetXaxis()->FindBin(4.01),hLzVsY->GetXaxis()->FindBin(4.49));
//printf("bincenter %f,bin %d \n",hProjLzY0->GetXaxis()->GetBinCenter((hProjLzY0->GetXaxis()->FindBin(6.99))),hProjLzY0->GetXaxis()->FindBin(6.99));
printf("bincenter %f,bin %d \n",hLzVsY->GetXaxis()->GetBinCenter((hLzVsY->GetXaxis()->FindBin(1.51))),hLzVsY->GetXaxis()->FindBin(1.51));
printf("bincenter %f,bin %d \n",hLzVsY->GetXaxis()->GetBinCenter((hLzVsY->GetXaxis()->FindBin(1.99))),hLzVsY->GetXaxis()->FindBin(1.99));
printf("bincenter %f,bin %d \n",hLzVsY->GetXaxis()->GetBinCenter((hLzVsY->GetXaxis()->FindBin(3.51))),hLzVsY->GetXaxis()->FindBin(3.51));
printf("bincenter %f,bin %d \n",hLzVsY->GetXaxis()->GetBinCenter((hLzVsY->GetXaxis()->FindBin(3.49))),hLzVsY->GetXaxis()->FindBin(3.49));
//printf("bincenter %f,bin %d \n",hLzVsY->GetXaxis()->GetBinCenter((hLzVsY->GetXaxis()->FindBin(4.01))),hLzVsY->GetXaxis()->FindBin(4.01));
//printf("bincenter %f,bin %d \n",hLzVsY->GetXaxis()->GetBinCenter((hLzVsY->GetXaxis()->FindBin(4.51))),hLzVsY->GetXaxis()->FindBin(4.51));
	
	double mean[6],rangey[6],weightvalue[6];
	double counter=0.;

	for(int i=1; i <= hProjLzPt0->FindBin(7); i++){
			counter=counter+(hProjLzPt0->GetBinContent(i));

	}
	printf("Lz vs Pt range [0 0.5] GeV/c : Entries %f, Counter %f, frazione %f \n:",hProjLzPt0->GetEntries(),counter,counter/(hProjLzPt0->GetEntries()));
	 counter=0.;
	for(int i=1; i <= hProjLzPt1->FindBin(7); i++){
			counter=counter+(hProjLzPt1->GetBinContent(i));
	}
	printf("Lz vs Pt range [0.5 1] GeV/c : Entries %f, Counter %f, frazione %f \n:",hProjLzPt1->GetEntries(),counter,counter/(hProjLzPt1->GetEntries()));
	counter=0.;

	for(int i=1; i <= hProjLzPt2->FindBin(7); i++){
			counter=counter+(hProjLzPt2->GetBinContent(i));
	}
	printf("Lz vs Pt range [1 1.5] GeV/c : Entries %f, Counter %f, frazione %f \n",hProjLzPt2->GetEntries(),counter,counter/(hProjLzPt2->GetEntries()));
double countery=0;
double meany=0;

printf("_____________________________________________\n ");

	printf("Lz vs Y range [1.5 2 ] : ");
	weightvalue[0]=weightmean(hProjLzY0);
	printf("Lz vs Y range [2 2.5 ] : ");
	weightvalue[1]=weightmean(hProjLzY1);
	printf("Lz vs Y range [2.5 3 ] : ");
	weightvalue[2]=weightmean(hProjLzY2); 
	printf("Lz vs Y range [3 3.5 ] : ");
	weightvalue[3]=weightmean(hProjLzY3);
	printf("Lz vs Y range [3.5 4 ] : ");
	weightvalue[4]=weightmean(hProjLzY4);
	printf("Lz vs Y range [4 4.5 ] : ");
	weightvalue[5]=weightmean(hProjLzY5);
	rangey[0]=(1.5+2)/2.;
	rangey[1]=(2+2.5)/2.;
	rangey[2]=(2.5+3)/2.;
	rangey[3]=(3+3.5)/2.;
	rangey[4]=(3.5+4)/2.;
	rangey[5]=(4+4.5)/2.;

	mean[0]=hProjLzY0->GetMean();
	mean[1]=hProjLzY1->GetMean();
	mean[2]=hProjLzY2->GetMean();
	mean[3]=hProjLzY3->GetMean();
	mean[4]=hProjLzY4->GetMean();
	mean[5]=hProjLzY5->GetMean();
	

	printf("Valori medi Lz \n");
	printf("Lz vs Y range [1.5 2] valor medio= %f \n",hProjLzY0->GetMean()); 
	printf("Lz vs Y range [2 2.5] valor medio= %f \n",hProjLzY1->GetMean());
	printf("Lz vs Y range [2.5 3] valor medio= %f \n",hProjLzY2->GetMean());
	printf("Lz vs Y range [3 3.5] valor medio= %f \n",hProjLzY3->GetMean());
	printf("Lz vs Y range [3.5 4] valor medio= %f \n",hProjLzY4->GetMean());
	printf("Lz vs Y range [4 4.5] valor medio= %f \n",hProjLzY5->GetMean());
	TGraph* gr = new TGraph(6,rangey,mean);
  gr->SetTitle("Mean value of Lz Vs rapidity");
  gr->GetXaxis()->SetTitle("y");
  gr->GetYaxis()->SetTitle("Lz(cm)");
   gr->Draw("AC*");
	TGraph* gr1 = new TGraph(6,rangey,weightvalue);
  gr1->SetTitle("Fraction of #Lambda with Lz<7cm Vs rapidity");
  gr1->GetXaxis()->SetTitle("y");
  gr1->GetYaxis()->SetTitle("fraction of #Lambda with Lz < 7cm");
  gr1->Draw("AC*");
	

//Checking reconstruction	
/*	double Resyprot[6]={0.} , RMSyprot[6]={0.} , Lzrange[6]={0.};
	TH1D* hProjResyLzprot0 = hResPyVsLzprot->ProjectionY("Res Py Lz range [0 1] ",hResPyVsLzprot->GetXaxis()->FindBin(0.01),hResPyVsLzprot->GetXaxis()->FindBin(1.01));
	Resyprot[0]=hProjResyLzprot0->GetMean();
	RMSyprot[0]=hProjResyLzprot0->GetRMS();
	Lzrange[0]=0.5;
	TH1D* hProjResyLzprot1 = hResPyVsLzprot->ProjectionY("Res Py Lz range [1 2] ",hResPyVsLzprot->GetXaxis()->FindBin(1.01),hResPyVsLzprot->GetXaxis()->FindBin(2.01));
	Resyprot[1]=hProjResyLzprot1->GetMean();
	RMSyprot[1]=hProjResyLzprot1->GetRMS();
	Lzrange[1]=1.5;
	TH1D* hProjResyLzprot2 = hResPyVsLzprot->ProjectionY("Res Py Lz range [2 3] ",hResPyVsLzprot->GetXaxis()->FindBin(2.01),hResPyVsLzprot->GetXaxis()->FindBin(3.01));
	Resyprot[2]=hProjResyLzprot2->GetMean();
	RMSyprot[2]=hProjResyLzprot2->GetRMS();
	Lzrange[2]=2.5;
	TH1D* hProjResyLzprot3 = hResPyVsLzprot->ProjectionY("Res Py Lz range [3 4] ",hResPyVsLzprot->GetXaxis()->FindBin(3.01),hResPyVsLzprot->GetXaxis()->FindBin(4.01));
	Resyprot[3]=hProjResyLzprot3->GetMean();
	RMSyprot[3]=hProjResyLzprot3->GetRMS();
	Lzrange[3]=3.5;
	TH1D* hProjResyLzprot4 = hResPyVsLzprot->ProjectionY("Res Py Lz range [4 5] ",hResPyVsLzprot->GetXaxis()->FindBin(4.01),hResPyVsLzprot->GetXaxis()->FindBin(5.01));
	Resyprot[4]=hProjResyLzprot4->GetMean();
	RMSyprot[4]=hProjResyLzprot4->GetRMS();
	Lzrange[4]=4.5;
	TH1D* hProjResyLzprot5 = hResPyVsLzprot->ProjectionY("Res Py Lz range [5 6] ",hResPyVsLzprot->GetXaxis()->FindBin(5.01),hResPyVsLzprot->GetXaxis()->FindBin(6.01));
	Resyprot[5]=hProjResyLzprot5->GetMean();
	RMSyprot[5]=hProjResyLzprot5->GetRMS();
	Lzrange[5]=5.5;

	TGraph* grResyProt = new TGraph(6,Lzrange,Resyprot);
	grResyProt->SetTitle("Mean Py vs Lz Proton"); 	
  grResyProt->Draw("AC*");
	TGraph* grRMSyProt = new TGraph(6,Lzrange,RMSyprot);
	grRMSyProt->SetTitle("RMS Py vs Lz Proton");
  grRMSyProt->Draw("AC*");

	double Resypion[6]={0.} , RMSypion[6]={0.};
	TH1D* hProjResyLzpion0 = hResPyVsLzpion->ProjectionY("Res Py Lz range [0 1] ",hResPyVsLzpion->GetXaxis()->FindBin(0.01),hResPyVsLzpion->GetXaxis()->FindBin(1.01));
	Resypion[0]=hProjResyLzpion0->GetMean();
	RMSypion[0]=hProjResyLzpion0->GetRMS();
	TH1D* hProjResyLzpion1 = hResPyVsLzpion->ProjectionY("Res Py Lz range [1 2] ",hResPyVsLzpion->GetXaxis()->FindBin(1.01),hResPyVsLzpion->GetXaxis()->FindBin(2.01));
	Resypion[1]=hProjResyLzpion1->GetMean();
	RMSypion[1]=hProjResyLzpion1->GetRMS();
	TH1D* hProjResyLzpion2 = hResPyVsLzpion->ProjectionY("Res Py Lz range [2 3] ",hResPyVsLzpion->GetXaxis()->FindBin(2.01),hResPyVsLzpion->GetXaxis()->FindBin(3.01));
	Resypion[2]=hProjResyLzpion2->GetMean();
	RMSypion[2]=hProjResyLzpion2->GetRMS();
	TH1D* hProjResyLzpion3 = hResPyVsLzpion->ProjectionY("Res Py Lz range [3 4] ",hResPyVsLzpion->GetXaxis()->FindBin(3.01),hResPyVsLzpion->GetXaxis()->FindBin(4.01));
	Resypion[3]=hProjResyLzpion3->GetMean();
	RMSypion[3]=hProjResyLzpion3->GetRMS();
	TH1D* hProjResyLzpion4 = hResPyVsLzpion->ProjectionY("Res Py Lz range [4 5] ",hResPyVsLzpion->GetXaxis()->FindBin(4.01),hResPyVsLzpion->GetXaxis()->FindBin(5.01));
	Resypion[4]=hProjResyLzpion4->GetMean();
	RMSypion[4]=hProjResyLzpion4->GetRMS();
	TH1D* hProjResyLzpion5 = hResPyVsLzpion->ProjectionY("Res Py Lz range [5 6] ",hResPyVsLzpion->GetXaxis()->FindBin(5.01),hResPyVsLzpion->GetXaxis()->FindBin(6.01));
	Resypion[5]=hProjResyLzpion5->GetMean();
	RMSypion[5]=hProjResyLzpion5->GetRMS();
	
	TGraph* grResyPion = new TGraph(6,Lzrange,Resypion);
	grResyPion->SetTitle("Mean Py vs Lz Pion");
  grResyPion->Draw("AC*");
	TGraph* grRMSyPion = new TGraph(6,Lzrange,RMSypion);
	grRMSyPion->SetTitle("RMS Py vs Lz Pion");
  grRMSyPion->Draw("AC*");*/

  fout->cd();
	hMassinrecVszP->Write();
	hMassinrecVsNfake->Write();
	hMassreconstr->Write();
//histograms for checking reconstruction (con il tracciatore)
/*	hResPyVsLzprotcorr->Write();
	hResPyVsLzpioncorr->Write();
	hResPyVsYprotcorr->Write();
	hResPyVsYpioncorr->Write();
	hResPycorrection->Write();
	hnfakeVsLzprot->Write();
	hResPyVsnfakeprot->Write();
	hPrecgenXprot->Write();	
	hPrecgenYprot->Write();
	hPrecgenZprot->Write();
	hPrecgenXpion->Write();	
	hPrecgenYpion->Write();
	hPrecgenZpion->Write();
	hResPxVsYprot->Write();
	hResPyVsYprot->Write();
	hResPzVsYprot->Write();
	hResPxVsYpion->Write();
	hResPyVsYpion->Write();
	hResPzVsYpion->Write();
	hResPxVsLzprot->Write();
	hResPyVsLzprot->Write();
	hResPzVsLzprot->Write();
	hResPxVsLzpion->Write();
	hResPyVsLzpion->Write();
	hResPzVsLzpion->Write();
	*/
//histograms for checking simulation (senza il tracciatore)
	gr1->Write();
	gr->Write();
	hProjLzY0->Write();
	hProjLzY1->Write();
	hProjLzY2->Write();
	hProjLzY3->Write();
	hProjLzY4->Write();
	hProjLzY5->Write();
	hLzVsPt->Write();
	hLzVsY->Write();
	hLenght->Write();
	hLenghtPtotDist->Write();
	hTime->Write();
	hctauPtotDist->Write();
	hTimeDist->Write();
	hsecVxGen->Write();
	hsecVyGen->Write(); 
	hsecVzGen->Write(); 
	hVxGen->Write();
	hVyGen->Write();
	hVzGen->Write();
	hMassLambdatrue->Write();
	hMassLambdafalse->Write();
//histograms
  hMassAll->Write();
  hMassFake->Write();
  hMassRefl->Write();
  hMassVsPt->Write();
  hMassVsY->Write();
  hMassReflVsPt->Write();
  hMassReflVsY->Write();
  hYPtGen->Write();
  hPtGen->Write();
  hYGen->Write();
  hYPtRecoAll->Write();  
  hYPtRecoFake->Write();
  hPtRecoAll->Write();
  hPtGenRecoAll->Write();
  hPtRecoVsGenAll->Write();
  hDiffPtRecoGenAll->Write();
  hYRecoAll->Write();
  hYGenRecoAll->Write();
  hPtRecoFake->Write();
  hDistXY->Write();
  hDist->Write();
  hDistgenXY->Write();
  hDistgen->Write();
  hCosp->Write();
  hDCA->Write();
  hDCAx->Write();
  hDCAy->Write();
  hDCAz->Write();
  hResVx->Write();
  hResVy->Write();
  hResVz->Write();
  hResPx->Write();
  hResPy->Write();
  hResPz->Write();
  hResVxVsY->Write();
  hResVyVsY->Write();
  hResVzVsY->Write();
  hResPxVsY->Write();
  hResPyVsY->Write();
  hResPzVsY->Write();
  hResDist->Write();
  hResDistXY->Write();
  hd0XYprod->Write();
  hd0XY1->Write();
  hd0XY2->Write();
  hNevents->Write();
  hsp->Write();
  if (ntD0cand){
    fnt->cd();
    ntD0cand->Write();
    fnt->Close();
  }
  // TCanvas *ccdau = new TCanvas();
  // ccdau->Divide(3, 2);
  // ccdau->cd(1)->SetLogy();
  // hPtGen->Draw();
  // ccdau->cd(2)->SetLogy();
  // hptK->Draw();
  // ccdau->cd(3)->SetLogy();
  // hptPi->Draw();
  // ccdau->cd(4)->SetLogy();
  // hYGen->Draw();
  // ccdau->cd(5)->SetLogy();
  // hyK->Draw();
  // ccdau->cd(6)->SetLogy();
  // hyPi->Draw();
  
  TFile fout2("DecayHistos.root", "RECREATE");
  // TFile fout2("DecayHistostest.root", "RECREATE");
  hPtGen->Write();
  hptP->Write();
  hptPi->Write();
  hyP->Write();
  hyPi->Write();
  hYGen->Write();
  hyPiP->Write();
  fout2.Close();

  fout->Close();

//	fntupla->Write();
//	fntupla->Close();
//	hfile->Write();
//	hfile->Close();
}

double weightmean(TH1D* hproject){
int counter=0;
double wmean=0.,frac=0.;
for(int i=1; i <= hproject->FindBin(6.99); i++){
			counter=counter+(hproject->GetBinContent(i));
			wmean=wmean+((hproject->GetXaxis()->GetBinCenter(i))*(hproject->GetBinContent(i)));
			//printf("tot conteggi in bin %d = %f \n",i,countery);
			//printf("tot conteggi in bin %d = %f \n",i,hProjLzY0->GetBinContent(i));
	}
	printf(" frazione di Lz<7 = %f e media pesata dei valori= %f \n",counter/(hproject->GetEntries()),wmean/counter);

return counter/(hproject->GetEntries()) ;
}


void MakeD0CombinBkgCandidates(const char* trackTreeFile="treeBkgEvents.root",
			       Int_t nevents = 999999, 
			       Int_t writeNtuple = kFALSE){

  // Read the TTree of tracks produced with runBkgVT.C
  // Create D0 combinatorial background candidates (= OS pairs of tracks)
  // Store in THnSparse and (optionally) TNtuple

  TFile *filetree = new TFile(trackTreeFile);
  TTree *tree = (TTree *)filetree->Get("tree");
  TClonesArray *arr = 0;
  tree->SetBranchAddress("tracks", &arr);
  Int_t entries = tree->GetEntries();
  printf("Number of events in tree = %d\n",entries);
  if(nevents>entries) nevents=entries;
  else printf(" --> Analysis performed on first %d events\n",nevents);

  TDatime dt;
  static UInt_t seed = dt.Get();
  gRandom->SetSeed(seed);

  // define mother particle
  Int_t pdgParticle = 421;
  
  TFile *fout = new TFile("D0-Bkg-histos.root", "recreate");
  TH1D* hPtRecoAll = new TH1D("hPtRecoAll", "Pt all match", 50, 0., 5.);
  TH2F* hYPtRecoAll = new TH2F("hYPtRecoAll", "Y-Pt all match", 40, 1., 5., 50, 0., 5.);
  TH1D* hMassAll = new TH1D("hMassAll", "Mass all match", 250, 0., 2.5);

  TH2F *hDistXY = new TH2F("hDistXY", "", 100, 0, 0.1, 30, 0, 3);
  TH2F *hDist = new TH2F("hDist", "", 300, 0, 10, 30, 0, 3);
  TH2F *hDistgenXY = new TH2F("hDistgenXY", "", 100, 0, 0.1, 30, 0, 3);
  TH2F *hDistgen = new TH2F("hDistgen", "", 300, 0, 10, 30, 0, 3);
  TH2F *hDCA = new TH2F("hDCA", "", 100, 0, 0.1, 30, 0, 3);
  TH2F *hDCAx = new TH2F("hDCAx", "", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hDCAy = new TH2F("hDCAy", "", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hDCAz = new TH2F("hDCAz", "", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hd0XYprod = new TH2F("hd0xyprod", "", 100, -0.01, 0.01, 30, 0, 3);
  TH2F *hd0XY1 = new TH2F("hd0xy1", "", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hd0XY2 = new TH2F("hd0xy2", "", 100, -0.1, 0.1, 30, 0, 3);
  
  TH2F *hResVx = new TH2F("hResVx", "", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hResVy = new TH2F("hResVy", "", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hResVz = new TH2F("hResVz", "", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hResPx = new TH2F("hResPx", "", 100, -1, 1, 30, 0, 3); //for Kaons
  TH2F *hResPy = new TH2F("hResPy", "", 100, -1, 1, 30, 0, 3);
  TH2F *hResPz = new TH2F("hResPz", "", 100, -1, 1, 30, 0, 3);
  TH2F *hResDist = new TH2F("hResDist", "", 100, -0.5, 0.5, 30, 0, 3);
  TH2F *hResDistXY = new TH2F("hResDistXY", "", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hCosp = new TH2F("hCosp", "", 100, -1, 1, 30, 0, 3);
  TH2F *hCosThStVsMass = new TH2F("hCosThStVsMass", "", 50, 1.5, 2.5, 40, -1, 1);
  
  TH2F *hd0 = new TH2F("hd0", "", 100, 0, 0.1, 30, 0, 3);
  
  TH1D *hcand = new TH1D("hcand", "", 1000, 0, 500000);
  TH1D *hcandpeak = new TH1D("hcandpeak", "", 500, 0, 15000);
  TH1D *hNevents = new TH1D("hNevents", "", 1, 0, 1);
    
  THnSparseF *hsp = CreateSparse();
  const Int_t nDim=static_cast<const Int_t>(hsp->GetNdimensions());
  Double_t arrsp[nDim];

  TFile *fnt = 0x0;
  TNtuple *ntD0cand = 0x0;
  Float_t arrnt[11];
  if (writeNtuple){
    fnt = new TFile("fntBkg.root", "recreate");
    ntD0cand = new TNtuple("ntD0cand", "ntD0cand", "mass:pt:y:dist:cosp:d01:d02:d0prod:dca:ptMin:ptMax", 32000);
  }

  KMCProbeFwd recProbe[2];
  TLorentzVector parent, daurec[2];
  Double_t massD0 = 1.864;
  for (Int_t iev = 0; iev < nevents; iev++){
    hNevents->Fill(0.5);
    Double_t vprim[3] = {0, 0, 0};
    Double_t countCandInPeak = 0;
    Double_t countCand = 0;
    tree->GetEvent(iev);
    Int_t arrentr = arr->GetEntriesFast();
    
    for (Int_t itr = 0; itr < arrentr; itr++){
      KMCProbeFwd *tr1 = (KMCProbeFwd *)arr->At(itr);
      // cout << "tr P=" << tr1->GetP() << endl;
      Float_t ch1 = tr1->GetCharge();
      for (Int_t itr2 = itr; itr2 < arrentr; itr2++){
	KMCProbeFwd *tr2 = (KMCProbeFwd *)arr->At(itr2);
	Float_t ch2 = tr2->GetCharge();
	if (ch1 * ch2 > 0) continue;
	if (ch1 < 0){ //convention: first track negative
	  recProbe[0] = *tr1;
	  recProbe[1] = *tr2;
	}else if (ch2 < 0){
	  recProbe[0] = *tr2;
	  recProbe[1] = *tr1;
	}
	recProbe[0].PropagateToDCA(&recProbe[1]);
	Double_t pxyz[3];
	recProbe[0].GetPXYZ(pxyz);
	
	Double_t pxyz2[3];
	recProbe[1].GetPXYZ(pxyz2);
	
	for(Int_t iMassHyp=0; iMassHyp<2; iMassHyp++){
	  // mass hypothesis: Kpi, piK
	  Int_t iKaon=-1;
	  if(iMassHyp==0){
	    daurec[0].SetXYZM(pxyz[0], pxyz[1], pxyz[2], KMCDetectorFwd::kMassK);
	    daurec[1].SetXYZM(pxyz2[0], pxyz2[1], pxyz2[2], KMCDetectorFwd::kMassPi);
	    iKaon=0;
	  }else{
	    daurec[0].SetXYZM(pxyz[0], pxyz[1], pxyz[2], KMCDetectorFwd::kMassPi);
	    daurec[1].SetXYZM(pxyz2[0], pxyz2[1], pxyz2[2], KMCDetectorFwd::kMassK);
	    iKaon=1;
	  }
	  parent = daurec[0];
	  parent += daurec[1];
	  countCand++;
	  Float_t ptD=parent.Pt();
	  Float_t invMassD=parent.M();
	  Float_t yD = 0.5 * TMath::Log((parent.E() + parent.Pz()) / (parent.E() - parent.Pz()));
	  hYPtRecoAll->Fill(yD, ptD);
	  hPtRecoAll->Fill(ptD);
	  hMassAll->Fill(invMassD);
	  if(invMassD>1.65  && invMassD<2.15){
	    // range to fill histos
	    if(invMassD>1.805 && invMassD<1.925) countCandInPeak++;
	    Float_t d1 = recProbe[1].GetX() - recProbe[0].GetX();
	    Float_t d2 = recProbe[1].GetY() - recProbe[0].GetY();
	    Float_t d3 = recProbe[1].GetZ() - recProbe[0].GetZ();
	    Float_t dca = sqrt(d1 * d1 + d2 * d2 + d3 * d3);
	    
	    //printf(" DCA = %f\n", sqrt(d1 * d1 + d2 * d2 + d3 * d3));
	    hDCA->Fill(dca, ptD);
	    hDCAx->Fill(d1, ptD);
	    hDCAy->Fill(d2, ptD);
	    hDCAz->Fill(d3, ptD);
	    // Float_t xP = (recProbe[1].GetX() + recProbe[0].GetX()) / 2.;
	    // Float_t yP = (recProbe[1].GetY() + recProbe[0].GetY()) / 2.;
	    // Float_t zP = (recProbe[1].GetZ() + recProbe[0].GetZ()) / 2.;
	    Double_t xP, yP, zP;
	    ComputeVertex(recProbe[0],recProbe[1],xP,yP,zP);
	    Float_t dist = TMath::Sqrt(xP * xP + yP * yP + zP * zP);
	    Float_t distXY = TMath::Sqrt(xP * xP + yP * yP);
	    Double_t vsec[3] = {xP, yP, zP};
	    Double_t cosp = CosPointingAngle(vprim, vsec, parent);
	    Double_t cts = CosThetaStar(parent,daurec[iKaon]);
	    Double_t ipD = ImpParXY(vprim, vsec, parent);
	    hCosp->Fill(cosp, ptD);
	    hCosThStVsMass->Fill(invMassD,cts);
	    //printf(" ***** ***** cos point = %f \n", cosp);	    
	    hDistXY->Fill(distXY, ptD);
	    hDist->Fill(dist, ptD);
	    
	    recProbe[0].PropagateToZBxByBz(0);
	    Double_t d0x1 = recProbe[0].GetX();
	    Double_t d0y1 = recProbe[0].GetY();
	    Double_t d0xy1 = TMath::Sqrt(d0x1 * d0x1 + d0y1 * d0y1);
	    if (d0x1 < 0) d0xy1 *= -1;
	    
	    recProbe[1].PropagateToZBxByBz(0);
	    Double_t d0x2 = recProbe[1].GetX();
	    Double_t d0y2 = recProbe[1].GetY();
	    Double_t d0xy2 = TMath::Sqrt(d0x2 * d0x2 + d0y2 * d0y2);
	    if (d0x2 < 0) d0xy2 *= -1;
	      
	    //printf("d0xy1 = %f, d0xy2 = %f \n", d0xy1, d0xy2);
	    
	    hd0XYprod->Fill(d0xy1 * d0xy2, ptD);
	    hd0XY1->Fill(d0xy1, ptD);
	    hd0XY2->Fill(d0xy2, ptD);
	    arrsp[0] = invMassD;
	    arrsp[1] = ptD;
	    arrsp[2] = yD;
	    arrsp[3] = dist;
	    arrsp[4] = cosp;
	    arrsp[5] = TMath::Min(TMath::Abs(d0xy1),TMath::Abs(d0xy2));
	    arrsp[6] = d0xy1 * d0xy2;
	    arrsp[7] = dca;
	    arrsp[8] = TMath::Min(recProbe[0].GetTrack()->Pt(),recProbe[1].GetTrack()->Pt());
	    arrsp[9] = TMath::Abs(ipD);	    
	    arrsp[10] = cts;
	    hsp->Fill(arrsp);
	    
	    if (ntD0cand){
	      arrnt[0] = invMassD;
	      arrnt[1] = ptD;
	      arrnt[2] = yD;
	      arrnt[3] = dist;
	      arrnt[4] = cosp;
	      arrnt[5] = d0xy1;
	      arrnt[6] = d0xy2;
	      arrnt[7] = d0xy1 * d0xy2;
	      arrnt[8] = dca;
	      arrnt[9] = TMath::Min(recProbe[0].GetTrack()->Pt(),recProbe[1].GetTrack()->Pt());
	      arrnt[10] = TMath::Max(recProbe[0].GetTrack()->Pt(),recProbe[1].GetTrack()->Pt());
	      ntD0cand->Fill(arrnt);
	    }
	    
	  } // check on inv mass
	} // loop on mass hypothesis
      } // loop on first track
    } // loop on second track
    
    hcand->Fill(countCand);
    hcandpeak->Fill(countCandInPeak);
    printf(" --> Event %d, tot D0 candidates = %.0f  in peak = %.0f\n",iev,countCand,countCandInPeak);
  }
  
  fout->cd();
  hNevents->Write();
  hcand->Write();
  hcandpeak->Write();  
  hMassAll->Write();
  hYPtRecoAll->Write();
  hPtRecoAll->Write();
  hDistXY->Write();
  hDist->Write();
  hDCA->Write();
  hDCAx->Write();
  hDCAy->Write();
  hDCAz->Write();
  hCosp->Write();
  hCosThStVsMass->Write();
  hd0XYprod->Write();
  hd0XY1->Write();
  hd0XY2->Write();
  hsp->Write();
  fout->Close();
  if (ntD0cand){
    fnt->cd();
    ntD0cand->Write();
    fnt->Close();
  }
}




Double_t CosPointingAngle(Double_t vprim[3], Double_t vsec[3], TLorentzVector &parent)
{

    // /// XY
    // /// Cosine of pointing angle in transverse plane assuming it is produced
    // /// at "point"

    // TVector3 momXY(parent.Px(), parent.Py(), 0.);
    // TVector3 flineXY(vsec[0] - vprim[0],
    //                  vsec[1] - vprim[1],
    //                  0.);

    // Double_t ptot2 = momXY.Mag2() * flineXY.Mag2();
    // if (ptot2 <= 0)
    // {
    //     return 0.0;
    // }
    // else
    // {
    //     Double_t cos = momXY.Dot(flineXY) / TMath::Sqrt(ptot2);
    //     if (cos > 1.0)
    //         cos = 1.0;
    //     if (cos < -1.0)
    //         cos = -1.0;
    //     return cos;
    // }
    /// Cosine of pointing angle in space assuming it is produced at "point"

    TVector3 mom(parent.Px(), parent.Py(), parent.Pz());
    TVector3 fline(vsec[0] - vprim[0],
                   vsec[1] - vprim[1],
                   vsec[2] - vprim[2]);

    Double_t ptot2 = mom.Mag2() * fline.Mag2();
    if (ptot2 <= 0)
    {
        return 0.0;
    }
    else
    {
        Double_t cos = mom.Dot(fline) / TMath::Sqrt(ptot2);
        if (cos > 1.0)
            cos = 1.0;
        if (cos < -1.0)
            cos = -1.0;
        return cos;
    }
}

Double_t ImpParXY(Double_t vprim[3], Double_t vsec[3], TLorentzVector &parent){

  Double_t k = -(vsec[0]-vprim[0])*parent.Px()-(vsec[1]-vprim[1])*parent.Py();
  k /= parent.Pt()*parent.Pt();
  Double_t dx = vsec[0]-vprim[0]+k*parent.Px();
  Double_t dy = vsec[1]-vprim[1]+k*parent.Py();
  Double_t absImpPar = TMath::Sqrt(dx*dx+dy*dy);
  TVector3 mom(parent.Px(), parent.Py(), parent.Pz());
  TVector3 fline(vsec[0] - vprim[0],
		 vsec[1] - vprim[1],
		 vsec[2] - vprim[2]);
  TVector3 cross = mom.Cross(fline);
  return (cross.Z()>0. ? absImpPar : -absImpPar);
}

Double_t CosThetaStar(TLorentzVector &parent, TLorentzVector &dauk) {

  Double_t massMoth = TDatabasePDG::Instance()->GetParticle(421)->Mass();
  Double_t massK = TDatabasePDG::Instance()->GetParticle(321)->Mass();
  Double_t massPi = TDatabasePDG::Instance()->GetParticle(211)->Mass();

  Double_t pStar = TMath::Sqrt((massMoth*massMoth-massK*massK-massPi*massPi)*(massMoth*massMoth-massK*massK-massPi*massPi)-4.*massK*massK*massPi*massPi)/(2.*massMoth);

  Double_t pMoth=parent.P();
  Double_t e=TMath::Sqrt(massMoth*massMoth+pMoth*pMoth);
  Double_t beta = pMoth/e;
  Double_t gamma = e/massMoth;
  TVector3 momDau(dauk.Px(),dauk.Py(),dauk.Pz());
  TVector3 momMoth(parent.Px(),parent.Py(),parent.Pz());
  Double_t qlProng=momDau.Dot(momMoth)/momMoth.Mag();
  Double_t cts = (qlProng/gamma-beta*TMath::Sqrt(pStar*pStar+massK*massK))/pStar;

  return cts;
}

THnSparseF* CreateSparse(){
  const Int_t nAxes=11;
  TString axTit[nAxes]={"Inv. mass (GeV/c^{2})","p_{T} (GeV/c)","y",
			"Dec Len (cm)","cos(#vartheta_{p})",
			"d_0^{min} (cm)",
			"d_0*d_0 (cm^{2})","DCA",
			"p_{T}^{min} (GeV/c)",
			"d_0^{D} (cm)","cos(#theta*)"};
  Int_t bins[nAxes] =   {100,   5,  40, 30,  20,   10,   10,      12,      8,  16,  10}; 
  Double_t min[nAxes] = {1.65,  0., 1., 0., 0.98, 0.,   -0.0006, 0.0,   0.,  0.,  -1.};
  Double_t max[nAxes] = {2.15,  5., 5., 0.3, 1.,   0.05, 0.,      0.03,  4.,  0.04, 1.};  
  THnSparseF *hsp = new THnSparseF("hsp", "hsp", nAxes, bins, min, max);
  for(Int_t iax=0; iax<nAxes; iax++) hsp->GetAxis(iax)->SetTitle(axTit[iax].Data());
  return hsp;
}


void ComputeVertex(KMCProbeFwd &t0, KMCProbeFwd &t1, Double_t &xV, Double_t &yV, Double_t &zV){
  // Based on AliESDv0

  //Trivial estimation of the vertex parameters
  Double_t tmp[3];
  t0.GetXYZ(tmp);
  Double_t  x1=tmp[0],  y1=tmp[1],  z1=tmp[2];
  const Double_t ss=0.0005*0.0005;//a kind of a residual misalignment precision
  Double_t sx1=t0.GetSigmaX2()+ss, sy1=t0.GetSigmaY2()+ss;
  t1.GetXYZ(tmp);
  Double_t  x2=tmp[0],  y2=tmp[1],  z2=tmp[2];
  Double_t sx2=t1.GetSigmaX2()+ss, sy2=t1.GetSigmaY2()+ss; 
    
  Double_t wx1=sx2/(sx1+sx2), wx2=1.- wx1;
  Double_t wy1=sy2/(sy1+sy2), wy2=1.- wy1;
  Double_t wz1=0.5, wz2=1.- wz1;
  xV=wx1*x1 + wx2*x2; 
  yV=wy1*y1 + wy2*y2; 
  zV=wz1*z1 + wz2*z2;
  return;
}
