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
				bool writeNtuple = kFALSE, 
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
 //TFile *fntupla = new TFile("ntupla.root", "recreate");
	float variable[31];
	//TNtuple *nt = new	TNtuple("nt","Variablesforperformance","nfaketrkprot:nfaketrkpion:ygen:xP:yP:zP:Vxgen:Vygen:Vzgen:massinvrec:pxgenprot:pygenprot:pzgenprot:pxgenpion:pygenpion:pzgenpion:dca:chi2prot:chi2ITSprot:chi2pion:chi2ITSpion:pxrecprot:pyrecprot:pzrecprot:pxrecpion:pyrecpion:pzrecpion:massinvswitch:pxlambda:pylambda:pzlambda");

  
   
//istogrammi dei prodotti che saranno salvati in decayhistos.root
  TH2F *hptP = new TH2F("hptP", "Protons from Lambda decays", 50,0.,10.,50, 0., 10.);
  TH2F *hptPi = new TH2F("hptPi", "pions from Lambda decays", 50, 0.,10.,50,0., 10.);
  TH1D *hyP = new TH1D("hyP", "y Protons from Lambda decays", 50, 0., 5.);
  TH1D *hyPi = new TH1D("hyPi", "y pions from Lambda decays", 50, 0., 5.);
  TH2F *hyPiP = new TH2F("hyPiP", "y pions vs y Protons from Lambda decays", 50, 0., 5., 50, 0., 5.);
 
  
  TFile *fout = new TFile("PROTON-Signal-histos.root", "recreate");
  
  //int outN = nev/10;
  //if (outN<1) outN=1;

//_____________________________SETUP SPERIMENTALE
  KMCDetectorFwd *det = new KMCDetectorFwd();
  det->ReadSetup(setup, setup);
  det->InitBkg(Eint); //check modificare bkg su altro file
  det->ForceLastActiveLayer(det->GetLastActiveLayerITS()); // will not propagate beyond VT
  det->SetMinITSHits(7);
  //det->SetMinITSHits(det->GetNumberOfActiveLayersITS()); //NA60+
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


  // prepare decays
  //TGenPhaseSpace decay;//Utility class to generate n-body event, with constant cross-section (default) or with Fermi energy dependence (opt="Fermi"). The event is generated in the center-of-mass frame, 												but the decay products are finally boosted using the betas of the original particle.

  TLorentzVector parentgen, daugen[2], parent, daurec[2], parentrefl, daurecswapmass[2]; 
  KMCProbeFwd recProbe[2];  

  //TParticleClass: Description of the dynamic properties of a particle. A dynamic particle class created by event generators and used during the propagation in detectors.
  TClonesArray *particles = new TClonesArray("TParticle", 1000);

  TLorentzVector *mom = new TLorentzVector(); //quadrimpulso particella lambda
  
  Int_t pdgParticle = 2212;// define PROTON particle
  
	TH3F* hYPtLzGen=new TH3F("hYPtLzGen", "Y-Pt-Lz corr match Gen", 80, 0, 6, 40, ptminSG, ptmaxSG,100,0,50);
	TH3F* hYPtLzRec=new TH3F("hYPtLzRec", "Y-Pt-Lz corr match Rec", 80, 0, 6, 40, ptminSG, ptmaxSG,100,0,50);
	TH3F* hYPtLzMC=new TH3F("hYPtLzMC", "Y-Pt-Lz corr match MC", 80, 0, 6, 40, ptminSG, ptmaxSG,100,0,50);
  TH2F* hYPtGen = new TH2F("hYPtGen", "Y-Pt corr match", 80, 1.0, 5.4, 40, ptminSG, ptmaxSG);
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

  TH2F *hResVxVsY = new TH2F("hResVxVsY", "", 200, -1000., 1000., 50, 0, 5);
  TH2F *hResVyVsY = new TH2F("hResVyVsY", "", 200, -1000., 1000., 50, 0, 5);
  TH2F *hResVzVsY = new TH2F("hResVzVsY", "", 200, -1000., 1000., 50, 0, 5);

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


//histograms for checking reconstruction (quindi tieni le righe di codice nel passaggio rivelatore)
	TH1D* hPrecXprot = new TH1D("hPrecXprot", "P_{Xrec} prot", 300, -5, 5);//width=0.01
	TH1D* hPrecYprot = new TH1D("hPrecYprot", "P_{Yrec} prot", 300, -2, 2);//width=0.01
	TH1D* hPrecZprot = new TH1D("hPrecZprot", "P_{Zrec} prot", 300, 0, 50);//width=0.01
	TH1D* hPrecgenXprot = new TH1D("hPrecgenXprot", "P_{Xrec}-P_{Xgen} prot", 300, -0.5, 0.5);//width=0.01
	TH1D* hPrecgenYprot = new TH1D("hPrecgenYprot", "P_{Yrec}-P_{Ygen} prot", 300, -0.5, 0.5);//width=0.01
	TH1D* hPrecgenZprot = new TH1D("hPrecgenZprot", "P_{Zrec}-P_{Zgen} prot", 300, -0.5, 0.5);//width=0.01
	TH2F* hResPxVsYprot = new TH2F("hResPxVsYprot", "P_{Xrec}-P_{Xgen}  vs Y prot", 300, -0.5 , 0.5,60, 0., 6);//width=0.01
	TH2F* hResPyVsYprot = new TH2F("hResPyVsYprot", "P_{Yrec}-P_{Ygen}  vs Y prot", 300, -0.5 , 0.5,60, 0., 6);
	TH2F* hResPzVsYprot = new TH2F("hResPzVsYprot", "P_{Zrec}-P_{Zgen}  vs Y prot", 300, -0.5 , 0.5,60, 0., 6);
	TH2F* hResPxVsYprotRel = new TH2F("hResPxVsYprotRel", "P_{Xrec}-P_{Xgen}  vs Y prot relativa", 300, -0.5 , 0.5,60, 0., 6);//width=0.01
	TH2F* hResPyVsYprotRel = new TH2F("hResPyVsYprotRel", "P_{Yrec}-P_{Ygen}  vs Y prot relativa", 300, -0.5 , 0.5,60, 0., 6);
	TH2F* hResPzVsYprotRel = new TH2F("hResPzVsYprotRel", "P_{Zrec}-P_{Zgen}  vs Y prot relativa", 300, -0.5 , 0.5,60, 0., 6);
//Nfake
	TH2F* hnfakeVsLzprot = new TH2F("hnfakeVsLzprot", "L_{z} Proton vs Nfake ", 6, -0.5 , 5.5,500, 0, 50);//width=0.1Lz
	TH2F* hResPyVsnfakeprot= new TH2F("hResPyVsnfakeprot", " P_{Yrec}-P_{Ygen} vs Nfake Proton", 6, -0.5 , 5.5,300, -1 , 1);//width=0.01Lz
	TH1F* hResPycorrection= new TH1F("hResPycorrection", " P_{Ynotcorrect}-P_{Ycorrect} ",300, -1 , 1);

//efficiency
  TH1F* hEfficiency = new TH1F("hEfficiency","Efficiency #pi+",32,6.5,38.5);
  
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
  	Double_t vprim[3];
  //for(int i=0;i<30;i++){
    	float nfaketrk = 0, nfaketrkprot = 0, nfaketrkpion = 0;
	  int count=0;
    double countrec=0;
     vprim[0] = 1;
     vprim[1] = 1; 
     vprim[2] = 10;
   
	  //GENERAZIONE
    for (Int_t iev = 0; iev < nevents; iev++){
      hNevents->Fill(0.5);
      //Double_t vprim[3] = {1, 1, 8};//vertice primario
      if(iev%10000==0) printf(" ***************  ev = %d \n", iev);
      int nrec = 0;
      int nfake = 0;
      nfaketrk = 0;
      nfaketrkprot = 0;
      nfaketrkpion = 0;
      double pxyz[3], pProtRec[3], pPionRec[3], pProtGen[3], pPionGen[3];
      
      if (simulateBg && (iev%refreshBg)==0) det->GenBgEvent(0.,0.,0.);//bkg c'è sempre
      Double_t ptGenD = 0.5; // get Lambda distribution from file
      Double_t yGenD = 3;
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
      
      double diff;
      Int_t arrpdgdau[2];
      Double_t ptK = -999.;
      Double_t ptPi = -999.;
      Double_t yK=-999.;
      Double_t yPi = -999.;
      Int_t icount = 0;
      Double_t secvertgenK[3]={0.,0.,0.};
      Double_t secvertgenPi[3]={0.,0.,0.};
      Double_t ptot = 0, massal = 0, p = 0, Respxyzprot[3] = {0.,0.,0.};
      //printf("evento= %d, numero prodotti:%d \n",iev,np);  
      hYGen->Fill(yGenD);
      hPtGen->Fill(ptGenD);
      hYPtGen->Fill(yGenD, ptGenD);

      TLorentzVector *pDecDau = new TLorentzVector(0., 0., 0., 0.);//quadrimpulso figli
      pDecDau->SetXYZM(pxGenD, pyGenD, pzGenD, mass);

      Int_t crg=1;//setta il segno della carica
      count++;
      if (det->SolveSingleTrack(pDecDau->Pt(), pDecDau->Rapidity(), pDecDau->Phi(), mass, crg, vprim[0], vprim[1], vprim[2], 0, 1, 99)){
        KMCProbeFwd *trw = det->GetLayer(0)->GetWinnerMCTrack();
        if (trw){
          //print getwinnermc e chi2
          if (trw->GetNormChi2(kTRUE) < ChiTot){
            nrec++;
            nfake += trw->GetNFakeITSHits();
            trw->GetPXYZ(pxyz);
            nfaketrk = trw->GetNFakeITSHits();
            printf("evento %d Macro: numhits %d, num fake %f \n",iev,trw->GetNITSHits(),nfaketrk);
						printf("Macro: getInnerLayerCheck %d \n",trw->GetInnerLayerChecked());
            printf("Bit patt :");
            for(unsigned int mask = 0x80000000; mask; mask >>= 1){
            printf("%d", !!(mask & trw->GetHitsPatt()));
            }
            printf("\n");
            /*printf("pxrec %f pxgen %f \n", pxyz[0],pxGenD);
            printf("pyrec %f pygen %f \n", pxyz[1],pyGenD);
            printf("pzrec %f pzgen %f \n", pxyz[2],pzGenD);*/
            daugen[0].SetXYZM(pxGenD, pyGenD, pzGenD, mass);
            daurec[0].SetXYZM(pxyz[0], pxyz[1], pxyz[2], mass);
            recProbe[0] = *trw; 
            //printf("z %f \n",recProbe[0].GetX());
            hPrecXprot->Fill(pxyz[0]);
            hPrecYprot->Fill(pxyz[1]);
            hPrecZprot->Fill(pxyz[2]);	
            hPrecgenXprot->Fill(pxyz[0]-pxGenD);
            hPrecgenYprot->Fill(pxyz[1]-pyGenD);
            hPrecgenZprot->Fill(pxyz[2]-pzGenD);
            Respxyzprot[0]=(pxyz[0]-pxGenD);
            Respxyzprot[1]=(pxyz[1]-pyGenD);
            Respxyzprot[2]=(pxyz[2]-pzGenD);
            if(nfaketrk == 0){
              hResPxVsYprot->Fill(Respxyzprot[0],yGenD);
              hResPyVsYprot->Fill(Respxyzprot[1],yGenD);
              hResPzVsYprot->Fill(Respxyzprot[2],yGenD);
              hResPxVsYprotRel->Fill(Respxyzprot[0]/pxyz[0],yGenD);
              hResPyVsYprotRel->Fill(Respxyzprot[1]/pxyz[1],yGenD);
              hResPzVsYprotRel->Fill(Respxyzprot[2]/pxyz[2],yGenD);
            }
            countrec++;
          }	
        }
      }
    } //event loop

    hEfficiency->Fill(25,countrec/nevents);
    hEfficiency->SetMarkerStyle(20);
    hEfficiency->SetLineColor(kBlue);
  //}
  /*hMassAll->SetLineColor(kBlue);
  hMassAll->Draw();
  hMassAll->SetMinimum(0.1);
  hMassFake->SetLineColor(kRed);
  hMassFake->Draw("same");
	*/

  fout->cd();
  hEfficiency->Write();
	hYPtLzMC->Write();
	hYPtLzGen->Write();
	hYPtLzRec->Write();
	//histograms for checking reconstruction (con il tracciatore)

	hResPyVsnfakeprot->Write();
	hPrecXprot->Write();
	hPrecYprot->Write();
	hPrecZprot->Write();
	hPrecgenXprot->Write();	
	hPrecgenYprot->Write();
	hPrecgenZprot->Write();
	hResPxVsYprot->Write();
	hResPyVsYprot->Write();
	hResPzVsYprot->Write();
	hResPxVsYprotRel->Write();
	hResPyVsYprotRel->Write();
	hResPzVsYprotRel->Write();
	
	
  TFile fout2("DecayHistosAntiProton.root", "RECREATE");
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

return wmean/counter ;
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
