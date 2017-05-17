#include "utility.h"

#define CHECK_EVERY_N 100000

#define PFMET_CUT 0.0 // global PFMET cut (both SR and CR): met > this
#define TKMET_CUT 0.0 // global TKMET cut (both SR and CR): met > this
#define QCD_INVERTED_DXY_THR 0.00 // cut dxy > this  //0.005

#define PFMET_SR 0.0 // pfmet > this value in SR
#define PFMET_CR 0.0 // pfmet < this value in CR

using namespace std;

static Int_t nMtBins = 120;
static Double_t mtMin = 20;
static Double_t mtMax = 140;
static Int_t nMt2over4Bins = 300;
static Double_t mt2over4Min = 100;
static Double_t mt2over4Max = 4900;
static string recoilCorrectionsPathToFile = "/afs/cern.ch/user/m/mciprian/www/wmass/analysisPlots_8TeV/metResoResp_mz70to110_eleTrigID_noHLTmatch/";
static Bool_t correctRecoilResponse = false; 

//=============================================================

void fillHistograms(const string& inputDIR = "./", const string& outputDIR = "./", 
		    const Sample& sample = Sample::zjets, 
		    TFile* outputFile = NULL, 
		    const Bool_t QCD_enriched_region = false, const Bool_t isMuon = false) {

  cout << endl;
  cout << "================================================" << endl;
  cout << endl;

  if (outputFile == NULL) {
    cout << "Error: file is NULL. please check. Exit ..." << endl;
    exit(EXIT_FAILURE);
  }

  outputFile->cd();

  // define electron ID at 8 TeV from --> https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaCutBasedIdentification
  // use use tight trigger ID: for quantities not mentioned in this ID, use those from lepton loose ID
  // will use PF isolation, we don't have detector based ID in ntuples yet
  // electronID tightEleID_EB_8TeV(0.004, 0.03, 0.01, 0.12, 0.02, 0.1, 0.05);
  // electronID tightEleID_EE_8TeV(0.005, 0.02, 0.03, 0.10, 0.02, 0.1, 0.05); 
  electronID tightEleID_EB_8TeV(0.007, 0.15, 0.01, 0.12, 0.02, 0.2, 0.05, 1, 1);
  electronID tightEleID_EE_8TeV(0.009, 0.10, 0.03, 0.10, 0.02, 0.2, 0.05, 1, 1); 

  // triggering MVA ID --> https://twiki.cern.ch/twiki/bin/viewauth/CMS/MultivariateElectronIdentification#Triggering_MVA
  electronTriggeringMVAID innerEB_TrigMVAID_highPt(0.1, 0.15, 0, 1);
  electronTriggeringMVAID outerEB_TrigMVAID_highPt(0.85, 0.15, 0, 1);
  electronTriggeringMVAID EE_TrigMVAID_highPt(0.92, 0.15, 0, 1);

  // electronTriggeringMVAID innerEB_TrigMVAID_lowPt(0.0, 0.15, 0, 1);
  // electronTriggeringMVAID outerEB_TrigMVAID_lowPt(0.10, 0.15, 0, 1);
  // electronTriggeringMVAID EE_TrigMVAID_lowPt(0.62, 0.15, 0, 1);
    
  electronID* eleID = NULL; // choose in loop for each event between EB or EE id
  electronTriggeringMVAID* eleTrigMVAID = NULL; //

  Double_t eleIso03thr_QCD = 0.1; 
  Double_t eleIso04thr_QCD = 0.15;
  Double_t muIso04thr = 0.12;       // 0.15,  0.2,   0.125
  Double_t muIso04thr_QCD = 0.12;       // 0.15,  e.g., sig region for iso < 0.15 and QCD region for iso > 0.2

  // be consistent with threshold when assigning histogram boundaries

  Int_t nBins_lepIso03 = 25;  // 20
  Double_t min_lepIso03 = 0.0; // 0.0
  Double_t max_lepIso03 = 0.1; // 0.060
  Int_t nBins_lepIso04 = 48;  // 30
  Double_t min_lepIso04 = 0.0; // 0.0
  Double_t max_lepIso04 = 0.12; // 0.150

  if (not isMuon) {
    nBins_lepIso04 = 60;  // 35,  30
    min_lepIso04 = 0.0; // 0.15, 0.2
    max_lepIso04 = 0.15; // 0.500
  }

  if (QCD_enriched_region) {
    nBins_lepIso03 = 180;  // 90, 85
    min_lepIso03 = 0.1; // 0.05, 0.075
    max_lepIso03 = 1.0; // 0.50
    nBins_lepIso04 = 88;  // 35,  30
    min_lepIso04 = 0.12; // 0.15, 0.2
    max_lepIso04 = 1.00; // 0.500
    if (not isMuon) {
      nBins_lepIso04 = 85;  // 35,  30
      min_lepIso04 = 0.15; // 0.15, 0.2
      max_lepIso04 = 1.000; // 0.500
    }
  }

  Int_t chargedLeptonFlavour = isMuon ? 13 : 11;
  Double_t lepEtaThreshold = isMuon ? 2.1 : 1.479;  // EB only for electron for now
  Int_t lepTightIdThreshold = isMuon ? 1 : 3;

  TDirectory *dirSample = NULL;
  string sampleDir = getStringFromEnumSample(sample).c_str();
  cout << "Sample --> " << sampleDir << endl;
  if (outputFile->GetKey(sampleDir.c_str())) dirSample = outputFile->GetDirectory(sampleDir.c_str());
  else dirSample = outputFile->mkdir(sampleDir.c_str());
  // if (sample == Sample::data_doubleEG)
  //   dirSample = outputFile->mkdir("data_doubleEG");
  // else if (sample == Sample::data_singleEG)
  //   dirSample = outputFile->mkdir("data_singleEG");
  // if (sample == Sample::data_doubleMu)
  //   dirSample = outputFile->mkdir("data_doubleMu");
  // else if (sample == Sample::data_singleMu)
  //   dirSample = outputFile->mkdir("data_singleMu");
  // else if (sample == Sample::zjets) 
  //   dirSample = outputFile->mkdir("zjets");
  // else if (sample == Sample::wjets)
  //   dirSample = outputFile->mkdir("wjets");
  // else if (sample == Sample::qcd_mu)
  //   dirSample = outputFile->mkdir("qcd_mu");
  // else if (sample == Sample::qcd_ele)
  //   dirSample = outputFile->mkdir("qcd_ele");
  // else if (sample == Sample::top)
  //   dirSample = outputFile->mkdir("top");
  // else if (sample == Sample::diboson)
  //   dirSample = outputFile->mkdir("diboson");

  dirSample->cd();

  gROOT->SetBatch(kTRUE);
  TH1::SetDefaultSumw2(); //all the following histograms will automatically call TH1::Sumw2()  

  cout << endl;

  TChain* chain = new TChain("treeProducerWMassEle");  
  TChain* friendChain = new TChain("mjvars/t");  // leave as NULL if you don't use friend trees
  TChain* SfFriendChain = NULL;
  if (sampleDir.find("data") == string::npos) SfFriendChain = new TChain("sf/t");  // leave as NULL if you don't use friend trees

  vector<Double_t> genwgtVec;
  buildChain8TeV(chain, genwgtVec, inputDIR, sample, friendChain, SfFriendChain); 

  // change directory again, when building chain something was messed up
  dirSample->cd();
  //  cout << "check" << endl;

  TTreeReader reader (chain);

  TTreeReaderValue<Int_t> isData  (reader,"isData");

  // reco met
  TTreeReaderValue<Float_t> tkmet    (reader,"tkmet_pt");
  TTreeReaderValue<Float_t> tkmet_phi(reader,"tkmet_phi");
  TTreeReaderValue<Float_t> pfmet    (reader,"met_pt");
  TTreeReaderValue<Float_t> pfmet_phi(reader,"met_phi");

  // LepGood for muons, LepCorr for electrons
  string lepVarToUse = isMuon ? "nLepGood" : "nLepCorr";
  TTreeReaderValue<Int_t> nlep  (reader,lepVarToUse.c_str());
  lepVarToUse = isMuon ? "LepGood_pt" : "LepCorr_pt";
  TTreeReaderArray<Float_t> lep_pt (reader,lepVarToUse.c_str());
  lepVarToUse = isMuon ? "LepGood_eta" : "LepCorr_eta";
  TTreeReaderArray<Float_t> lep_eta (reader,lepVarToUse.c_str());

  // LepGood branch
  TTreeReaderArray<Int_t> lep_pdgId (reader,"LepGood_pdgId");
  TTreeReaderArray<Float_t> lep_phi (reader,"LepGood_phi");
  TTreeReaderArray<Float_t> lep_mass (reader,"LepGood_mass");
  TTreeReaderArray<Float_t> lep_relIso03 (reader,"LepGood_relIso03");
  TTreeReaderArray<Float_t> lep_relIso04 (reader,"LepGood_relIso04");
  TTreeReaderArray<Int_t> lep_tightId (reader,"LepGood_tightId");  // must make sure it is tuned on 8 TeV (it might be based on 13 TeV) 

  // for electronID
  TTreeReaderArray<Float_t> lep_detaIn (reader,"LepGood_detaIn");
  TTreeReaderArray<Float_t> lep_dphiIn (reader,"LepGood_dphiIn");
  TTreeReaderArray<Float_t> lep_sigmaIetaIeta (reader,"LepGood_sigmaIetaIeta");
  TTreeReaderArray<Float_t> lep_HoE (reader,"LepGood_hcalOverEcal");
  TTreeReaderArray<Float_t> lep_dxy (reader,"LepGood_dxy");
  TTreeReaderArray<Float_t> lep_dz (reader,"LepGood_dz");
  TTreeReaderArray<Int_t> lep_lostHits (reader,"LepGood_lostHits");
  TTreeReaderArray<Int_t> lep_convVeto (reader,"LepGood_convVeto");
  TTreeReaderArray<Float_t> lep_ecalEnergy (reader,"LepGood_ecalEnergy");
  TTreeReaderArray<Float_t> lep_eSuperClusterOverP (reader,"LepGood_eSuperClusterOverP");
  // 1/E - 1/p = 1/ecalEnergy - eSuperClusterOverP/ecalEnergy == (1 - eSuperClusterOverP)/ecalEnergy;

  // related to electron triggering MVA ID
  //TTreeReaderArray<Float_t> lep_sip3ds (reader,"LepGood_sip3d");
  TTreeReaderArray<Float_t> lep_mvaId (reader,"LepGood_mvaId");  // value of non triggering MVA ID
  TTreeReaderArray<Float_t> lep_mvaIdTrig (reader,"LepGood_mvaIdTrig"); // value of triggering MVA ID

  TTreeReaderArray<Int_t> lep_eleMVAId (reader,"LepGood_eleMVAId"); // 
  TTreeReaderArray<Int_t> lep_convVetoFull (reader,"LepGood_convVetoFull"); // like conVetoFull, but also includes no missing hits
  TTreeReaderArray<Float_t> lep_eleMVAPreselId (reader,"LepGood_eleMVAPreselId");

  // other electron related branch
  TTreeReaderArray<Float_t> lep_etaSc (reader,"LepGood_scEta");
  TTreeReaderArray<Float_t> lep_r9 (reader,"LepGood_r9");
  
  // MC reweight
  // must activate it only for non data sample
  ///////////////////////////////
  // use a dummy variable in trees to activate TTreeReaderValue even for data
  string dummybranch = "rho";
  if (sampleDir.find("data") == string::npos) dummybranch = "puWeight";  // PU weight. Can use this in tree because it was computed with all 8 TeV dataset
  TTreeReaderValue<Float_t> puw (reader,dummybranch.c_str());  // if running on data, create it as a branch existing in tree (it won't be used)
  dummybranch = "m2l";
  if (sampleDir.find("data") == string::npos) dummybranch = "SF_LepTight_1l";
  TTreeReaderValue<Float_t> lepEfficiency (reader,dummybranch.c_str());

  // Match of lepton with trigger object 
  TTreeReaderValue<Int_t> HLT_SingleMu (reader,"HLT_SingleMu");
  TTreeReaderValue<Int_t> HLT_SingleEl (reader,"HLT_SingleEl");
  //TTreeReaderArray<Float_t> lep_trgMatch (reader,"LepGood_trgMatch");

  //cout << "check 2" << endl;

  // TH1
  TH1D* hmT = new TH1D("hmT","",nMtBins, mtMin, mtMax);
  TH1D* hmT2over4 = new TH1D("hmT2over4","",nMt2over4Bins, mt2over4Min, mt2over4Max);
  TH1D* hpfmet = new TH1D("hpfmet","",60,0,300);
  TH1D* htkmet = new TH1D("htkmet","",60,0,300);
  TH1D* hlep1pt = new TH1D("hlep1pt","",28,24,80);  //60, 0, 300
  TH1D* hlep1pt2 = new TH1D("hlep1pt2","",125,500,3000);
  TH1D* hlep1eta = new TH1D("hlep1eta","",48,-2.4,2.4);  //60, 0, 300
  TH1D* hlep2pt = new TH1D("hlep2pt","",40,0,80);
  TH1D* hbosonpt = new TH1D("hbosonpt","",30,0,30);
  //TH1D* hbosoneta = new TH1D("hbosoneta","",100,-5,5);
  TH1D* hlep1relIso03 = new TH1D("hlep1relIso03","",nBins_lepIso03,min_lepIso03,max_lepIso03);
  TH1D* hlep1relIso04 = new TH1D("hlep1relIso04","",nBins_lepIso04,min_lepIso04,max_lepIso04);
  TH1D* hrecoil = new TH1D("hrecoil","",30,0,30);
  TH1D* hdxy = new TH1D("hdxy","",40,0,0.1);
  TH1D* hdphiLepMet = new TH1D("hdphiLepMet","",32,0,3.2);

  // TH2
  TH2D* h2_mT_lep1pt = new TH2D("h2_mT_lep1pt","",nMtBins, mtMin, mtMax, 28,24,80);
  TH2D* h2_mT_lep1eta = new TH2D("h2_mT_lep1eta","",nMtBins, mtMin, mtMax, 30, -1.5, 1.5);
  TH2D* h2_mT_pfmet = new TH2D("h2_mT_pfmet","",nMtBins, mtMin, mtMax, 60,0,300);
  TH2D* h2_mT_tkmet = new TH2D("h2_mT_tkmet","",nMtBins, mtMin, mtMax, 60,0,300);
  TH2D* h2_mT_bosonPt = new TH2D("h2_mT_bosonPt","",nMtBins, mtMin, mtMax, 30,0,30);
  TH2D* h2_mT_lep1relIso03 = new TH2D("h2_mT_lep1relIso03","",nMtBins, mtMin, mtMax, nBins_lepIso03,min_lepIso03,max_lepIso03);
  TH2D* h2_mT_lep1relIso04 = new TH2D("h2_mT_lep1relIso04","",nMtBins, mtMin, mtMax, nBins_lepIso04,min_lepIso04,max_lepIso04);

  TH2D* h2_lep1pt_lep1eta = new TH2D("h2_lep1pt_lep1eta","",28,24,80, 30, -1.5, 1.5);
  TH2D* h2_lep1pt_pfmet = new TH2D("h2_lep1pt_pfmet","",28,24,80, 60,0,300);
  TH2D* h2_lep1pt_tkmet = new TH2D("h2_lep1pt_tkmet","",28,24,80, 60,0,300);
  TH2D* h2_lep1pt_bosonPt = new TH2D("h2_lep1pt_bosonPt","",28,24,80, 30,0,30);
  TH2D* h2_lep1pt_lep1relIso03 = new TH2D("h2_lep1pt_lep1relIso03","",28,24,80, nBins_lepIso03,min_lepIso03,max_lepIso03);
  TH2D* h2_lep1pt_lep1relIso04 = new TH2D("h2_lep1pt_lep1relIso04","",28,24,80, nBins_lepIso04,min_lepIso04,max_lepIso04);

  TH2D* h2_pfmet_tkmet = new TH2D("h2_pfmet_tkmet","",60,0,300, 60,0,300);

  TH1D* hlep1relIso03_noIsoCut = new TH1D("hlep1relIso03_noIsoCut","",100,0.0,0.5);
  TH2D* h2_mT_lep1relIso03_noIsoCut = new TH2D("h2_mT_lep1relIso03_noIsoCut","",nMtBins, mtMin, mtMax,100,0.0,0.5);
  TH2D* h2_lep1pt_lep1relIso03_noIsoCut = new TH2D("h2_lep1pt_lep1relIso03_noIsoCut","",28,24,80, 100,0.0,0.5);
  TH1D* hlep1relIso04_noIsoCut = new TH1D("hlep1relIso04_noIsoCut","",100,0.0,0.5);
  TH2D* h2_mT_lep1relIso04_noIsoCut = new TH2D("h2_mT_lep1relIso04_noIsoCut","",nMtBins, mtMin, mtMax,100,0.0,0.5);
  TH2D* h2_lep1pt_lep1relIso04_noIsoCut = new TH2D("h2_lep1pt_lep1relIso04_noIsoCut","",28,24,80, 100,0.0,0.5);

  ////////////////////////////////
  // electron specific histograms
  ////////////////////////////////
  TH1D* hlep1sigIetaIeta = NULL;
  TH1D* hlep1r9 = NULL;
  TH2D* h2_mT_lep1sigIetaIeta = NULL;
  TH2D* h2_mT_lep1r9 = NULL;
  TH2D* h2_lep1pt_lep1sigIetaIeta = NULL;
  TH2D* h2_lep1pt_lep1r9 = NULL;
  TH1D* hdetaIn = NULL;
  TH1D* hdphiIn = NULL;
  // TH1D* hdetaIn_noCut = NULL;
  // TH1D* hdphiIn_noCut = NULL;
  // TH1D* hdxy_noCut = NULL;
  // TH2D* h2_mT_detaIn_noCut = NULL;
  // TH2D* h2_mT_dphiIn_noCut = NULL;
  // TH2D* h2_mT_dxy_noCut = NULL;
  // TH2D* h2_lep1pt_detaIn_noCut = NULL;
  // TH2D* h2_lep1pt_dphiIn_noCut = NULL;
  // TH2D* h2_lep1pt_dxy_noCut = NULL;
  if (not isMuon) {
    hlep1sigIetaIeta = new TH1D("hlep1sigIetaIeta","",25, 0.0, 0.025);
    hlep1r9 = new TH1D("hlep1r9","", 55, 0.0, 1.1);
    h2_mT_lep1sigIetaIeta = new TH2D("h2_mT_lep1sigIetaIeta","",nMtBins, mtMin, mtMax, 25, 0.0, 0.025);
    h2_mT_lep1r9 = new TH2D("h2_mT_lep1r9","",nMtBins, mtMin, mtMax, 55, 0.0, 1.1);
    h2_lep1pt_lep1sigIetaIeta = new TH2D("h2_lep1pt_lep1sigIetaIeta","",28,24,80, 25, 0.0, 0.025);
    h2_lep1pt_lep1r9 = new TH2D("h2_lep1pt_lep1r9","",28,24,80, 55, 0.0, 1.1);
    hdetaIn = new TH1D("hdetaIn","",300,0.0,0.300);
    hdphiIn = new TH1D("hdphiIn","",60,0.0,0.30);
    // hdetaIn_noCut = new TH1D("hdetaIn_noCut","",300,0.0,0.300);
    // hdphiIn_noCut = new TH1D("hdphiIn_noCut","",60,0.0,0.30);
    // hdxy_noCut = new TH1D("hdxy_noCut","",20,0.0,0.05);
    // h2_mT_detaIn_noCut = new TH2D("h2_mT_detaIn_noCut","",nMtBins, mtMin, mtMax, 300,0.0,0.300);
    // h2_mT_dphiIn_noCut = new TH2D("h2_mT_dphiIn_noCut","",nMtBins, mtMin, mtMax, 60,0.0,0.30);
    // h2_mT_dxy_noCut = new TH2D("h2_mT_dxy_noCut","",nMtBins, mtMin, mtMax, 20,0.0,0.05);
    // h2_lep1pt_detaIn_noCut = new TH2D("h2_lep1pt_detaIn_noCut","",28,24,80, 300,0.0,0.300);
    // h2_lep1pt_dphiIn_noCut = new TH2D("h2_lep1pt_dphiIn_noCut","",28,24,80, 60,0.0,0.30);
    // h2_lep1pt_dxy_noCut = new TH2D("h2_lep1pt_dxy_noCut","",28,24,80, 20,0.0,0.05);
  }
  ///////////////////////////////////

  /////////////////////////////////
  // only for positive lepton
  TH1D* hmT_plus = new TH1D("hmT_plus","",nMtBins, mtMin, mtMax);
  TH1D* hmT2over4_plus = new TH1D("hmT2over4_plus","",nMt2over4Bins, mt2over4Min, mt2over4Max);
  TH1D* hpfmet_plus = new TH1D("hpfmet_plus","",60,0,300);
  TH1D* htkmet_plus = new TH1D("htkmet_plus","",60,0,300);
  TH1D* hlep1pt_plus = new TH1D("hlep1pt_plus","",28,24,80);  //60, 0, 300
  TH1D* hlep1pt2_plus = new TH1D("hlep1pt2_plus","",125,500,3000);
  TH1D* hlep1eta_plus = new TH1D("hlep1eta_plus","",48,-2.4,2.4);  //60, 0, 300
  TH1D* hlep2pt_plus = new TH1D("hlep2pt_plus","",35,10,80);
  TH1D* hbosonpt_plus = new TH1D("hbosonpt_plus","",30,0,30);
  //TH1D* hbosoneta_plus = new TH1D("hbosoneta_plus","",100,-5,5);
  TH1D* hlep1sigIetaIeta_plus = new TH1D("hlep1sigIetaIeta_plus","",25, 0.0, 0.025);
  TH1D* hlep1relIso03_plus = new TH1D("hlep1relIso03_plus","",nBins_lepIso03,min_lepIso03,max_lepIso03);
  TH1D* hlep1relIso04_plus = new TH1D("hlep1relIso04_plus","",nBins_lepIso04,min_lepIso04,max_lepIso04);
  TH1D* hlep1r9_plus = new TH1D("hlep1r9_plus","", 55, 0.0, 1.1);
  TH1D* hrecoil_plus = new TH1D("hrecoil_plus","",30,0,30);
  TH1D* hdxy_plus = new TH1D("hdxy_plus","",20,0,0.1);
  TH1D* hdphiLepMet_plus = new TH1D("hdphiLepMet_plus","",32,0,3.2);

  TH2D* h2_mT_lep1pt_plus = new TH2D("h2_mT_lep1pt_plus","",nMtBins, mtMin, mtMax, 28,24,80);
  TH2D* h2_mT_lep1eta_plus = new TH2D("h2_mT_lep1eta_plus","",nMtBins, mtMin, mtMax, 30, -1.5, 1.5);
  TH2D* h2_mT_lep1sigIetaIeta_plus = new TH2D("h2_mT_lep1sigIetaIeta_plus","",nMtBins, mtMin, mtMax, 25, 0.0, 0.025);
  TH2D* h2_mT_lep1r9_plus = new TH2D("h2_mT_lep1r9_plus","",nMtBins, mtMin, mtMax, 55, 0.0, 1.1);
  TH2D* h2_mT_pfmet_plus = new TH2D("h2_mT_pfmet_plus","",nMtBins, mtMin, mtMax, 60,0,300);
  TH2D* h2_mT_tkmet_plus = new TH2D("h2_mT_tkmet_plus","",nMtBins, mtMin, mtMax, 60,0,300);
  TH2D* h2_mT_bosonPt_plus = new TH2D("h2_mT_bosonPt_plus","",nMtBins, mtMin, mtMax, 30,0,30);
  TH2D* h2_mT_lep1relIso03_plus = new TH2D("h2_mT_lep1relIso03_plus","",nMtBins, mtMin, mtMax, nBins_lepIso03,min_lepIso03,max_lepIso03);
  TH2D* h2_mT_lep1relIso04_plus = new TH2D("h2_mT_lep1relIso04_plus","",nMtBins, mtMin, mtMax, nBins_lepIso04,min_lepIso04,max_lepIso04);

  TH2D* h2_lep1pt_lep1eta_plus = new TH2D("h2_lep1pt_lep1eta_plus","",28,24,80, 30, -1.5, 1.5);
  TH2D* h2_lep1pt_lep1sigIetaIeta_plus = new TH2D("h2_lep1pt_lep1sigIetaIeta_plus","",28,24,80, 25, 0.0, 0.025);
  TH2D* h2_lep1pt_lep1r9_plus = new TH2D("h2_lep1pt_lep1r9_plus","",28,24,80, 55, 0.0, 1.1);
  TH2D* h2_lep1pt_pfmet_plus = new TH2D("h2_lep1pt_pfmet_plus","",28,24,80, 60,0,300);
  TH2D* h2_lep1pt_tkmet_plus = new TH2D("h2_lep1pt_tkmet_plus","",28,24,80, 60,0,300);
  TH2D* h2_lep1pt_bosonPt_plus = new TH2D("h2_lep1pt_bosonPt_plus","",28,24,80, 30,0,30);
  TH2D* h2_lep1pt_lep1relIso03_plus = new TH2D("h2_lep1pt_lep1relIso03_plus","",28,24,80, nBins_lepIso03,min_lepIso03,max_lepIso03);
  TH2D* h2_lep1pt_lep1relIso04_plus = new TH2D("h2_lep1pt_lep1relIso04_plus","",28,24,80, nBins_lepIso04,min_lepIso04,max_lepIso04);

  TH1D* hlep1relIso03_noIsoCut_plus = new TH1D("hlep1relIso03_noIsoCut_plus","",100,0.0,0.5);
  TH2D* h2_mT_lep1relIso03_noIsoCut_plus = new TH2D("h2_mT_lep1relIso03_noIsoCut_plus","",nMtBins, mtMin, mtMax,100,0.0,0.5);
  TH2D* h2_lep1pt_lep1relIso03_noIsoCut_plus = new TH2D("h2_lep1pt_lep1relIso03_noIsoCut_plus","",28,24,80, 100,0.0,0.5);
  TH1D* hlep1relIso04_noIsoCut_plus = new TH1D("hlep1relIso04_noIsoCut_plus","",100,0.0,0.5);
  TH2D* h2_mT_lep1relIso04_noIsoCut_plus = new TH2D("h2_mT_lep1relIso04_noIsoCut_plus","",nMtBins, mtMin, mtMax,100,0.0,0.5);
  TH2D* h2_lep1pt_lep1relIso04_noIsoCut_plus = new TH2D("h2_lep1pt_lep1relIso04_noIsoCut_plus","",28,24,80, 100,0.0,0.5);

  ///////////////////////////////
  // negative lepton
  TH1D* hmT_minus = new TH1D("hmT_minus","",nMtBins, mtMin, mtMax);
  TH1D* hmT2over4_minus = new TH1D("hmT2over4_minus","",nMt2over4Bins, mt2over4Min, mt2over4Max);
  TH1D* hpfmet_minus = new TH1D("hpfmet_minus","",60,0,300);
  TH1D* htkmet_minus = new TH1D("htkmet_minus","",60,0,300);
  TH1D* hlep1pt_minus = new TH1D("hlep1pt_minus","",28,24,80);  //60, 0, 300
  TH1D* hlep1pt2_minus = new TH1D("hlep1pt2_minus","",125,500,3000);
  TH1D* hlep1eta_minus = new TH1D("hlep1eta_minus","",48,-2.4,2.4);  //60, 0, 300
  TH1D* hlep2pt_minus = new TH1D("hlep2pt_minus","",35,10,80);
  TH1D* hbosonpt_minus = new TH1D("hbosonpt_minus","",30,0,30);
  //TH1D* hbosoneta_minus = new TH1D("hbosoneta_minus","",100,-5,5);
  TH1D* hlep1sigIetaIeta_minus = new TH1D("hlep1sigIetaIeta_minus","",25, 0.0, 0.025);
  TH1D* hlep1relIso03_minus = new TH1D("hlep1relIso03_minus","",nBins_lepIso03,min_lepIso03,max_lepIso03);
  TH1D* hlep1relIso04_minus = new TH1D("hlep1relIso04_minus","",nBins_lepIso04,min_lepIso04,max_lepIso04);
  TH1D* hlep1r9_minus = new TH1D("hlep1r9_minus","", 55, 0.0, 1.1);
  TH1D* hrecoil_minus = new TH1D("hrecoil_minus","",30,0,30);
  TH1D* hdxy_minus = new TH1D("hdxy_minus","",20,0,0.1);
  TH1D* hdphiLepMet_minus = new TH1D("hdphiLepMet_minus","",32,0,3.2);

  TH2D* h2_mT_lep1pt_minus = new TH2D("h2_mT_lep1pt_minus","",nMtBins, mtMin, mtMax, 28,24,80);
  TH2D* h2_mT_lep1eta_minus = new TH2D("h2_mT_lep1eta_minus","",nMtBins, mtMin, mtMax, 30, -1.5, 1.5);
  TH2D* h2_mT_lep1sigIetaIeta_minus = new TH2D("h2_mT_lep1sigIetaIeta_minus","",nMtBins, mtMin, mtMax, 25, 0.0, 0.025);
  TH2D* h2_mT_lep1r9_minus = new TH2D("h2_mT_lep1r9_minus","",nMtBins, mtMin, mtMax, 55, 0.0, 1.1);
  TH2D* h2_mT_pfmet_minus = new TH2D("h2_mT_pfmet_minus","",nMtBins, mtMin, mtMax, 60,0,300);
  TH2D* h2_mT_tkmet_minus = new TH2D("h2_mT_tkmet_minus","",nMtBins, mtMin, mtMax, 60,0,300);
  TH2D* h2_mT_bosonPt_minus = new TH2D("h2_mT_bosonPt_minus","",nMtBins, mtMin, mtMax, 30,0,30);
  TH2D* h2_mT_lep1relIso03_minus = new TH2D("h2_mT_lep1relIso03_minus","",nMtBins, mtMin, mtMax, nBins_lepIso03,min_lepIso03,max_lepIso03);
  TH2D* h2_mT_lep1relIso04_minus = new TH2D("h2_mT_lep1relIso04_minus","",nMtBins, mtMin, mtMax, nBins_lepIso04,min_lepIso04,max_lepIso04);

  TH2D* h2_lep1pt_lep1eta_minus = new TH2D("h2_lep1pt_lep1eta_minus","",28,24,80, 30, -1.5, 1.5);
  TH2D* h2_lep1pt_lep1sigIetaIeta_minus = new TH2D("h2_lep1pt_lep1sigIetaIeta_minus","",28,24,80, 25, 0.0, 0.025);
  TH2D* h2_lep1pt_lep1r9_minus = new TH2D("h2_lep1pt_lep1r9_minus","",28,24,80, 55, 0.0, 1.1);
  TH2D* h2_lep1pt_pfmet_minus = new TH2D("h2_lep1pt_pfmet_minus","",28,24,80, 60,0,300);
  TH2D* h2_lep1pt_tkmet_minus = new TH2D("h2_lep1pt_tkmet_minus","",28,24,80, 60,0,300);
  TH2D* h2_lep1pt_bosonPt_minus = new TH2D("h2_lep1pt_bosonPt_minus","",28,24,80, 30,0,30);
  TH2D* h2_lep1pt_lep1relIso03_minus = new TH2D("h2_lep1pt_lep1relIso03_minus","",28,24,80, nBins_lepIso03,min_lepIso03,max_lepIso03);
  TH2D* h2_lep1pt_lep1relIso04_minus = new TH2D("h2_lep1pt_lep1relIso04_minus","",28,24,80, nBins_lepIso04,min_lepIso04,max_lepIso04);

  TH1D* hlep1relIso03_noIsoCut_minus = new TH1D("hlep1relIso03_noIsoCut_minus","",100,0.0,0.5);
  TH2D* h2_mT_lep1relIso03_noIsoCut_minus = new TH2D("h2_mT_lep1relIso03_noIsoCut_minus","",nMtBins, mtMin, mtMax,100,0.0,0.5);
  TH2D* h2_lep1pt_lep1relIso03_noIsoCut_minus = new TH2D("h2_lep1pt_lep1relIso03_noIsoCut_minus","",28,24,80, 100,0.0,0.5);
  TH1D* hlep1relIso04_noIsoCut_minus = new TH1D("hlep1relIso04_noIsoCut_minus","",100,0.0,0.5);
  TH2D* h2_mT_lep1relIso04_noIsoCut_minus = new TH2D("h2_mT_lep1relIso04_noIsoCut_minus","",nMtBins, mtMin, mtMax,100,0.0,0.5);
  TH2D* h2_lep1pt_lep1relIso04_noIsoCut_minus = new TH2D("h2_lep1pt_lep1relIso04_noIsoCut_minus","",28,24,80, 100,0.0,0.5);

  ////////////////////
  ////////////////////
  // get recoil correction histogram as a function of boson pT  
  ////////////////////
  string recoilCorrFileName = "";
  TFile* recoilCorrFile = NULL;
  TH1D* hRecoilCorrectionTkmet = NULL;
  TH1D* hRecoilCorrectionPfmet = NULL;

  if (correctRecoilResponse) {

    recoilCorrFileName = recoilCorrectionsPathToFile;
    if (isMuon) recoilCorrFileName += "Zmumu/";
    else recoilCorrFileName += "Zee/";
    recoilCorrFileName += "wmass_resoresphists.root";

    recoilCorrFile = new TFile(recoilCorrFileName.c_str(),"READ");
    if (!recoilCorrFile || recoilCorrFile->IsZombie()) {
      cout << "Error: file with recoil corrections not opened. Exit" << endl;
      exit(EXIT_FAILURE);
    }

    if (sampleDir.find("data") != string::npos) hRecoilCorrectionTkmet = (TH1D*) getHistCloneFromFile(recoilCorrFile, "hRecoilResponseCorrFactor_tkmet_data", ""); 
    else hRecoilCorrectionTkmet = (TH1D*) getHistCloneFromFile(recoilCorrFile, "hRecoilResponseCorrFactor_tkmet_zjets", "");

    if (sampleDir.find("data") != string::npos) hRecoilCorrectionPfmet= (TH1D*) getHistCloneFromFile(recoilCorrFile, "hRecoilResponseCorrFactor_pfmet_data", ""); 
    else hRecoilCorrectionPfmet = (TH1D*) getHistCloneFromFile(recoilCorrFile, "hRecoilResponseCorrFactor_pfmet_zjets", "");

    checkNotNullPtr(hRecoilCorrectionTkmet,"hRecoilCorrectionTkmet");
    checkNotNullPtr(hRecoilCorrectionPfmet,"hRecoilCorrectionPfmet");
    
  }
  ////////////////////
  ////////////////////

  // start event loop                                                                                                                     
                    
  long int nTotal = chain->GetEntries();
  long int nEvents = 0;

  Bool_t negativeLeptonHasPassedSelection = false;
  Bool_t positiveLeptonHasPassedSelection = false;

  // flags for electrons to save whether some selection is passed and allow to fill some histograms before only one of these selection is applied
  // Bool_t passDphiSel = false;
  // Bool_t passDetaSel = false;
  // Bool_t passDxySel = false;
  Bool_t passIsoSel = false;

  Double_t wgt = 1.0;
  Double_t mT = 0.0;
  Double_t mT2over4 = 0.0;
  Double_t lep1pt = 0.0;
  Double_t lep1pt2 = 0.0;
  Double_t dphiLepMet = -999.0;
  ////////////////////////////////////////////
  // to get correct weight depending on sample in chain
  string currentFile = "";
  Int_t ifile = 0;
  ////////////////////

  while(reader.Next()){
  
    cout.flush();
    if(nEvents % CHECK_EVERY_N == 0) cout<<"\r"<<"Analyzing events "<<double(nEvents)/nTotal*100<<" % ";
    //cout << "entry : " << nEvents << endl;
    nEvents++;

    // to debug
    // if (nEvents == 100) break;

    if(dynamic_cast<TChain*>(reader.GetTree())->GetFile()->GetName() != currentFile and currentFile != ""){ 
      currentFile = dynamic_cast<TChain*>(reader.GetTree())->GetFile()->GetName();                   
      ifile ++;                                                                                      
    } else if(dynamic_cast<TChain*>(reader.GetTree())->GetFile()->GetName() != currentFile and currentFile == ""){
      currentFile = dynamic_cast<TChain*>(reader.GetTree())->GetFile()->GetName();                         
    }      

    negativeLeptonHasPassedSelection = false;
    positiveLeptonHasPassedSelection = false;
    // passDphiSel = false;
    // passDetaSel = false;
    // passDxySel = false;
    passIsoSel = false;
  
    // selection
    // if (*nlep != 1) continue;    // 1 leptons                                                  
    if (*nlep < 1) continue;    // at least 1 lepton, but use leading on to cut                                                   
    if (fabs(lep_pdgId[0]) != chargedLeptonFlavour) continue;  // electrons                                                     
    if (lep_pt[0] < 30.0) continue;
    if (fabs(lep_eta[0]) > lepEtaThreshold) continue; 
    //if (lep_trgMatch[0] < 0.5) continue;  // avoid using it for now, it could be ill-defined

    // remove Z candidate when nlep > 1
    Bool_t isCandidateZJetsEvent = false;
    if (*nlep > 1) {
      for (Int_t ilep = 1; ilep < *nlep && !isCandidateZJetsEvent; ilep++) {
	if ((lep_pdgId[0] + lep_pdgId[ilep]) == 0) {
	  TLorentzVector lep1, zboson;
	  lep1.SetPtEtaPhiM(lep_pt[0],lep_eta[0],lep_phi[0], lep_mass[0]) ;
	  zboson.SetPtEtaPhiM(lep_pt[ilep],lep_eta[ilep],lep_phi[ilep], lep_mass[ilep]);
	  zboson += lep1;
	  Double_t zmass = zboson.M();
	  if (zmass > 80.0 && zmass < 100.0) isCandidateZJetsEvent = true;
	}
      }
    }
    if (isCandidateZJetsEvent) continue;

    // tight ID
    // for electron don't cut now on quantities to be inverted to select QCD

    if (isMuon) {

      if (*HLT_SingleMu == 0) continue;
      if ( lep_tightId[0] < lepTightIdThreshold ) continue;   // tight ID

    } else {

      if (*HLT_SingleEl == 0) continue;

      // if (fabs(lep_etaSc[0]) < 1.479) eleID = &tightEleID_EB_8TeV;
      // else                            eleID = &tightEleID_EE_8TeV;
      // if (lep_sigmaIetaIeta[0] > eleID->sigmaIetaIeta) continue;
      // if (lep_HoE[0] > eleID->HoE) continue;
      // if (lep_dz[0] > eleID->dz) continue;
      // if ((1 - lep_eSuperClusterOverP[0])/lep_ecalEnergy[0] > eleID->invE_minus_invP) continue;
      // if (lep_lostHits[0] > eleID->missingHits) continue;
      // if (lep_convVeto[0] == 0) continue;

      // triggering ID
      // if (fabs(lep_eta[0]) < 0.8) eleTrigMVAID = &innerEB_TrigMVAID_highPt;
      // else if (fabs(lep_eta[0]) < 1.479) eleTrigMVAID = &outerEB_TrigMVAID_highPt;
      // else if (fabs(lep_eta[0]) < 2.5) eleTrigMVAID = &EE_TrigMVAID_highPt;
      // if (lep_mvaIdTrig[0] < eleTrigMVAID->mva) continue;
      // if (lep_lostHits[0] > eleTrigMVAID->maxMissingHits) continue;
      // if (lep_convVeto[0] < eleTrigMVAID->convVeto) continue;

      // implementing electron triggering ID
      if (lep_eleMVAPreselId[0] < 0.5) continue;      // trigger level ID 1 if ok, but since it is float, distinguish from 0 by asking > 0.5)
      if (lep_eleMVAId[0] < 2) continue;      
      if (lep_convVetoFull[0] != 1) continue; // it includes both vertex fit probability and missing hits selections
      // relative isolation cut applied afterwards, to allow for inversion depending on the region


    }

    if (lep_pdgId[0] > 0) negativeLeptonHasPassedSelection = true;
    else                  positiveLeptonHasPassedSelection = true;

    if (*isData == 1) wgt = 1.0;
    else wgt = intLumi * genwgtVec[ifile] * *puw * *lepEfficiency; 

    TLorentzVector lep1Reco;
    TVector2 lep1Reco2D, tkmet2D, pfmet2D, metReco, recoilPfmet, recoilTkmet, recoilReco, bosonReco;

    lep1Reco.SetPtEtaPhiM(lep_pt[0],lep_eta[0],lep_phi[0],lep_mass[0]); 
    lep1pt = lep1Reco.Pt();
    lep1Reco2D = lep1Reco.Vect().XYvector();
    tkmet2D.SetMagPhi(*tkmet,*tkmet_phi);
    recoilTkmet = -1.0 * (lep1Reco2D + tkmet2D);
    pfmet2D.SetMagPhi(*pfmet,*pfmet_phi);
    recoilPfmet = -1.0 * (lep1Reco2D + pfmet2D);
    if (useTrackMet) {
      metReco = tkmet2D;
      recoilReco = recoilTkmet;
    } else {
      metReco = pfmet2D;
      recoilReco = recoilPfmet;
    }
 
    ///////////////////////////////
    if (correctRecoilResponse) {

      // to get the "true" pT(W) bin, we use pT(W) computed with PFmet, which should be a more accurate estimate of pT(W) with respect to using trackMET
      // the problem is that pT(W) is computed using MET, which is not corrected yet
      // for W+jets MC we could use the gen level pT, but this cannot be used on data or other MC samples 

      Double_t wptUsingPFmet = (lep1Reco2D + pfmet2D).Mod();
      TH1D* hRecoilCorrection = useTrackMet ? hRecoilCorrectionTkmet : hRecoilCorrectionPfmet;
      Int_t recoilBin = hRecoilCorrection->GetXaxis()->FindFixBin(wptUsingPFmet);

      if (recoilBin <= 0) recoilBin = 1;
      else if (recoilBin > hRecoilCorrection->GetNbinsX()) recoilBin = hRecoilCorrection->GetNbinsX();

      recoilTkmet *= hRecoilCorrectionTkmet->GetBinContent(recoilBin);
      recoilPfmet *= hRecoilCorrectionPfmet->GetBinContent(recoilBin);
      tkmet2D = -1.0 * (recoilTkmet + lep1Reco2D);
      pfmet2D = -1.0 * (recoilPfmet + lep1Reco2D);

      // recoilReco *= hRecoilCorrection->GetBinContent(recoilBin);
      // metReco = -1.0 * (recoilReco + lep1Reco2D); 
      // bosonReco = lep1Reco2D + metReco;

      if (useTrackMet) {
	metReco = tkmet2D;
      } else {
	metReco = pfmet2D;
      }

    }
    ////////////////////////////
    bosonReco = lep1Reco2D + metReco;
    recoilReco = -1.0 * bosonReco;

    if (bosonReco.Mod() > 30) continue;
    if (pfmet2D.Mod() < PFMET_CUT) continue;    
    if (tkmet2D.Mod() < TKMET_CUT) continue;    

    mT = sqrt( 2. * lep1Reco2D.Mod() * metReco.Mod() * ( 1. - cos(lep1Reco2D.DeltaPhi(metReco)) ) );
    if (mT < mtMin) continue;
    //if (mT > mtMax) continue;    
    mT2over4 = mT*mT/4.0;
    lep1pt2 = lep1pt * lep1pt; 

    dphiLepMet = metReco.DeltaPhi(lep1Reco2D);
    //////////////////////////////////////////////
    // apply some selections but do not use continue. Just fill some bool to keep track of passed selections
    ////////////////////////////////////////////
    if (isMuon) {

      if (QCD_enriched_region) {

	if (pfmet2D.Mod() > PFMET_CR) continue;
	if (lep_dxy[0] < QCD_INVERTED_DXY_THR) continue; // apply dxy cut to select QCD in muon region
	if (lep_relIso04[0] > muIso04thr_QCD) passIsoSel = true; 

      } else {

	if (lep_relIso04[0] < muIso04thr) passIsoSel = true; 
	if (pfmet2D.Mod() < PFMET_SR) continue;

      }
      
    } else {

      if (QCD_enriched_region) {

	if (lep_relIso04[0] > eleIso04thr_QCD) passIsoSel = true;

      	// if (lep_detaIn[0] > eleID->dEtaIn) passDetaSel = true;
      	// if (lep_dphiIn[0] > eleID->dPhiIn) passDphiSel = true;
      	// if (lep_dxy[0] > eleID->dxy) passDxySel = true;
      	// if (lep_relIso03[0] > leptonIso03threshold) passIsoSel = true;  

      } else {

	//if (lep_relIso04[0] < eleTrigMVAID->relPFiso) passIsoSel = true;	
	if (lep_relIso04[0] < 0.15) passIsoSel = true;

      	// if (lep_detaIn[0] < eleID->dEtaIn) passDetaSel = true;
      	// if (lep_dphiIn[0] < eleID->dPhiIn) passDphiSel = true;
      	// if (lep_dxy[0] < eleID->dxy) passDxySel = true;
      	// if (lep_relIso03[0] < leptonIso03threshold) passIsoSel = true;  
      }

    }
    ////////////////////////////////

    /////////////////////////////////////
    // iso histo before iso cut
    ////////////////////////////////////////////

    // fill histograms for isolation before cutting on it
    fillTH1(hlep1relIso03_noIsoCut, lep_relIso03[0], wgt);
    fillTH2(h2_mT_lep1relIso03_noIsoCut, mT, lep_relIso03[0], wgt);
    fillTH2(h2_lep1pt_lep1relIso03_noIsoCut, lep1pt, lep_relIso03[0], wgt);
    fillTH1(hlep1relIso04_noIsoCut, lep_relIso04[0], wgt);
    fillTH2(h2_mT_lep1relIso04_noIsoCut, mT, lep_relIso04[0], wgt);
    fillTH2(h2_lep1pt_lep1relIso04_noIsoCut, lep1pt, lep_relIso04[0], wgt);

    if (negativeLeptonHasPassedSelection) {
      fillTH1(hlep1relIso03_noIsoCut_minus, lep_relIso03[0], wgt);
      fillTH2(h2_mT_lep1relIso03_noIsoCut_minus, mT, lep_relIso03[0], wgt);
      fillTH2(h2_lep1pt_lep1relIso03_noIsoCut_minus, lep1pt, lep_relIso03[0], wgt);
      fillTH1(hlep1relIso04_noIsoCut_minus, lep_relIso04[0], wgt);
      fillTH2(h2_mT_lep1relIso04_noIsoCut_minus, mT, lep_relIso04[0], wgt);
      fillTH2(h2_lep1pt_lep1relIso04_noIsoCut_minus, lep1pt, lep_relIso04[0], wgt);
    }

    if (positiveLeptonHasPassedSelection) {
      fillTH1(hlep1relIso03_noIsoCut_plus, lep_relIso03[0], wgt);
      fillTH2(h2_mT_lep1relIso03_noIsoCut_plus, mT, lep_relIso03[0], wgt);
      fillTH2(h2_lep1pt_lep1relIso03_noIsoCut_plus, lep1pt, lep_relIso03[0], wgt);
      fillTH1(hlep1relIso04_noIsoCut_plus, lep_relIso04[0], wgt);
      fillTH2(h2_mT_lep1relIso04_noIsoCut_plus, mT, lep_relIso04[0], wgt);
      fillTH2(h2_lep1pt_lep1relIso04_noIsoCut_plus, lep1pt, lep_relIso04[0], wgt);
    }

    /////////////////////////////////////

    /////////////////////////////////////
    // electron cut to be inverted in QCD region, fill histogram before each of those cuts
    ////////////////////////////////////////////
    // if ((not isMuon) && passIsoSel) {

    //   if (passDphiSel && passDxySel) {
    // 	fillTH1(hdetaIn_noCut, lep_detaIn[0], wgt);
    // 	fillTH2(h2_mT_detaIn_noCut, mT, lep_detaIn[0], wgt);
    // 	fillTH2(h2_lep1pt_detaIn_noCut, lep1pt, lep_detaIn[0], wgt);
    //   }
    //   if (passDetaSel && passDxySel) {
    // 	fillTH1(hdphiIn_noCut, lep_dphiIn[0], wgt);
    // 	fillTH2(h2_mT_dphiIn_noCut, mT, lep_dphiIn[0], wgt);
    // 	fillTH2(h2_lep1pt_dphiIn_noCut, lep1pt, lep_dphiIn[0], wgt);
    //   }
    //   if (passDetaSel && passDphiSel) {
    // 	fillTH1(hdxy_noCut, lep_dxy[0], wgt);
    // 	fillTH2(h2_mT_dxy_noCut, mT, lep_dxy[0], wgt);
    // 	fillTH2(h2_lep1pt_dxy_noCut, lep1pt, lep_dxy[0], wgt);
    //   }

    // }
    ////////////////////////////////////////////

    ////////////////////////////
    // LEP ISO CUT AND OTHERS TO SEPARATE QCD REGION
    /////////////////////////////
    if (not passIsoSel) continue;
    //cout << "passed iso sel" << endl;
 
    // if (isMuon) {
    //   if (not passIsoSel) continue;
    // }
    // else {
    //   if (QCD_enriched_region) {
    // 	if (!(passDetaSel || passDphiSel || passDxySel || passIsoSel)) continue;
    //   } else 
    // 	if (!(passDetaSel && passDphiSel && passDxySel && passIsoSel)) continue;
    // }

    fillTH1(hmT,(Double_t) mT, wgt);
    fillTH1(hmT2over4,(Double_t) mT2over4, wgt);
    fillTH1(hpfmet,(Double_t) pfmet2D.Mod(), wgt);
    fillTH1(htkmet,(Double_t) tkmet2D.Mod(), wgt);
    fillTH1(hlep1pt,(Double_t) lep1pt, wgt);
    fillTH1(hlep1pt2,(Double_t) lep1pt2, wgt);
    fillTH1(hlep1eta,(Double_t) lep1Reco.Eta(), wgt);
    fillTH1(hlep2pt,(Double_t) metReco.Mod(), wgt);
    fillTH1(hbosonpt,(Double_t) bosonReco.Mod(), wgt);
    //fillTH1(hbosoneta,(Double_t) bosonReco.Eta(), wgt);
    fillTH1(hlep1relIso03,(Double_t) lep_relIso03[0], wgt);
    fillTH1(hlep1relIso04,(Double_t) lep_relIso04[0], wgt);
    fillTH1(hrecoil, recoilReco.Mod(), wgt);
    fillTH1(hdxy,lep_dxy[0], wgt);
    fillTH1(hdphiLepMet, dphiLepMet, wgt);

    fillTH2(h2_mT_lep1pt, mT, lep_pt[0], wgt);
    fillTH2(h2_mT_lep1eta, mT, lep_eta[0], wgt);
    fillTH2(h2_mT_pfmet, mT, pfmet2D.Mod(), wgt);
    fillTH2(h2_mT_tkmet, mT, tkmet2D.Mod(), wgt);
    fillTH2(h2_mT_bosonPt, mT, bosonReco.Mod(), wgt);
    fillTH2(h2_mT_lep1relIso03, mT, lep_relIso03[0], wgt);
    fillTH2(h2_mT_lep1relIso04, mT, lep_relIso04[0], wgt);

    fillTH2(h2_lep1pt_lep1eta, lep_pt[0], lep_eta[0], wgt);
    fillTH2(h2_lep1pt_pfmet, lep_pt[0], pfmet2D.Mod(), wgt);
    fillTH2(h2_lep1pt_tkmet, lep_pt[0], tkmet2D.Mod(), wgt);
    fillTH2(h2_lep1pt_bosonPt, lep_pt[0], bosonReco.Mod(), wgt);
    fillTH2(h2_lep1pt_lep1relIso03, lep_pt[0], lep_relIso03[0], wgt);
    fillTH2(h2_lep1pt_lep1relIso04, lep_pt[0], lep_relIso04[0], wgt);

    fillTH2(h2_pfmet_tkmet, pfmet2D.Mod(), tkmet2D.Mod(), wgt); 

    if (not isMuon) {
      fillTH1(hlep1sigIetaIeta,(Double_t) lep_sigmaIetaIeta[0], wgt);
      fillTH1(hlep1r9,(Double_t) lep_r9[0], wgt);
      fillTH1(hdetaIn,(Double_t) lep_detaIn[0], wgt);    
      fillTH1(hdphiIn,(Double_t) lep_dphiIn[0], wgt);
      fillTH2(h2_mT_lep1sigIetaIeta, mT, lep_sigmaIetaIeta[0], wgt);
      fillTH2(h2_mT_lep1r9, mT, lep_r9[0], wgt);
      fillTH2(h2_lep1pt_lep1sigIetaIeta, lep_pt[0], lep_sigmaIetaIeta[0], wgt);
      fillTH2(h2_lep1pt_lep1r9, lep_pt[0], lep_r9[0], wgt); 
    }      

    if (negativeLeptonHasPassedSelection) {
      fillTH1(hmT_minus,(Double_t) mT, wgt);
      fillTH1(hmT2over4_minus,(Double_t) mT2over4, wgt);
      fillTH1(hpfmet_minus,(Double_t) pfmet2D.Mod(), wgt);
      fillTH1(htkmet_minus,(Double_t) tkmet2D.Mod(), wgt);
      fillTH1(hlep1pt_minus,(Double_t) lep1pt, wgt);
      fillTH1(hlep1pt2_minus,(Double_t) lep1pt2, wgt);
      fillTH1(hlep1eta_minus,(Double_t) lep1Reco.Eta(), wgt);
      fillTH1(hlep2pt_minus,(Double_t) metReco.Mod(), wgt);
      fillTH1(hbosonpt_minus,(Double_t) bosonReco.Mod(), wgt);
      //fillTH1(hbosoneta_minus,(Double_t) bosonReco.Eta(), wgt);
      fillTH1(hlep1relIso03_minus,(Double_t) lep_relIso03[0], wgt);
      fillTH1(hlep1relIso04_minus,(Double_t) lep_relIso04[0], wgt);
      fillTH1(hrecoil_minus, recoilReco.Mod(), wgt);
      fillTH1(hdxy_minus,lep_dxy[0], wgt);
      fillTH1(hdphiLepMet_minus, dphiLepMet, wgt);
    
      fillTH2(h2_mT_lep1pt_minus, mT, lep_pt[0], wgt);
      fillTH2(h2_mT_lep1eta_minus, mT, lep_eta[0], wgt);
      fillTH2(h2_mT_pfmet_minus, mT, pfmet2D.Mod(), wgt);
      fillTH2(h2_mT_tkmet_minus, mT, tkmet2D.Mod(), wgt);
      fillTH2(h2_mT_bosonPt_minus, mT, bosonReco.Mod(), wgt);
      fillTH2(h2_mT_lep1relIso03_minus, mT, lep_relIso03[0], wgt);
      fillTH2(h2_mT_lep1relIso04_minus, mT, lep_relIso04[0], wgt);

      fillTH2(h2_lep1pt_lep1eta_minus, lep_pt[0], lep_eta[0], wgt);
      fillTH2(h2_lep1pt_pfmet_minus, lep_pt[0], pfmet2D.Mod(), wgt);
      fillTH2(h2_lep1pt_tkmet_minus, lep_pt[0], tkmet2D.Mod(), wgt);
      fillTH2(h2_lep1pt_bosonPt_minus, lep_pt[0], bosonReco.Mod(), wgt);
      fillTH2(h2_lep1pt_lep1relIso03_minus, lep_pt[0], lep_relIso03[0], wgt);
      fillTH2(h2_lep1pt_lep1relIso04_minus, lep_pt[0], lep_relIso04[0], wgt);

      if (not isMuon) {
	fillTH1(hlep1sigIetaIeta_minus,(Double_t) lep_sigmaIetaIeta[0], wgt);
	fillTH1(hlep1r9_minus,(Double_t) lep_r9[0], wgt);
	fillTH2(h2_mT_lep1sigIetaIeta_minus, mT, lep_sigmaIetaIeta[0], wgt);
	fillTH2(h2_mT_lep1r9_minus, mT, lep_r9[0], wgt);
	fillTH2(h2_lep1pt_lep1sigIetaIeta_minus, lep_pt[0], lep_sigmaIetaIeta[0], wgt);
	fillTH2(h2_lep1pt_lep1r9_minus, lep_pt[0], lep_r9[0], wgt);
      }

    }

    if (positiveLeptonHasPassedSelection) {
      fillTH1(hmT_plus,(Double_t) mT, wgt);
      fillTH1(hmT2over4_plus,(Double_t) mT2over4, wgt);
      fillTH1(hpfmet_plus,(Double_t) pfmet2D.Mod(), wgt);
      fillTH1(htkmet_plus,(Double_t) tkmet2D.Mod(), wgt);
      fillTH1(hlep1pt_plus,(Double_t) lep1pt, wgt);
      fillTH1(hlep1pt2_plus,(Double_t) lep1pt2, wgt);
      fillTH1(hlep1eta_plus,(Double_t) lep1Reco.Eta(), wgt);
      fillTH1(hlep2pt_plus,(Double_t) metReco.Mod(), wgt);
      fillTH1(hbosonpt_plus,(Double_t) bosonReco.Mod(), wgt);
      //fillTH1(hbosoneta_plus,(Double_t) bosonReco.Eta(), wgt);
      fillTH1(hlep1relIso03_plus,(Double_t) lep_relIso03[0], wgt);
      fillTH1(hlep1relIso04_plus,(Double_t) lep_relIso04[0], wgt);
      fillTH1(hrecoil_plus, recoilReco.Mod(), wgt);
      fillTH1(hdxy_plus,lep_dxy[0], wgt);
      fillTH1(hdphiLepMet_plus, dphiLepMet, wgt);

      fillTH2(h2_mT_lep1pt_plus, mT, lep_pt[0], wgt);
      fillTH2(h2_mT_lep1eta_plus, mT, lep_eta[0], wgt);
      fillTH2(h2_mT_pfmet_plus, mT, pfmet2D.Mod(), wgt);
      fillTH2(h2_mT_tkmet_plus, mT, tkmet2D.Mod(), wgt);
      fillTH2(h2_mT_bosonPt_plus, mT, bosonReco.Mod(), wgt);
      fillTH2(h2_mT_lep1relIso03_plus, mT, lep_relIso03[0], wgt);
      fillTH2(h2_mT_lep1relIso04_plus, mT, lep_relIso04[0], wgt);

      fillTH2(h2_lep1pt_lep1eta_plus, lep_pt[0], lep_eta[0], wgt);
      fillTH2(h2_lep1pt_pfmet_plus, lep_pt[0], pfmet2D.Mod(), wgt);
      fillTH2(h2_lep1pt_tkmet_plus, lep_pt[0], tkmet2D.Mod(), wgt);
      fillTH2(h2_lep1pt_bosonPt_plus, lep_pt[0], bosonReco.Mod(), wgt);
      fillTH2(h2_lep1pt_lep1relIso03_plus, lep_pt[0], lep_relIso03[0], wgt);
      fillTH2(h2_lep1pt_lep1relIso04_plus, lep_pt[0], lep_relIso04[0], wgt);

      if (not isMuon) {
	fillTH1(hlep1r9_plus,(Double_t) lep_r9[0], wgt);
	fillTH1(hlep1sigIetaIeta_plus,(Double_t) lep_sigmaIetaIeta[0], wgt);
	fillTH2(h2_mT_lep1sigIetaIeta_plus, mT, lep_sigmaIetaIeta[0], wgt);
	fillTH2(h2_mT_lep1r9_plus, mT, lep_r9[0], wgt);
	fillTH2(h2_lep1pt_lep1sigIetaIeta_plus, lep_pt[0], lep_sigmaIetaIeta[0], wgt);
	fillTH2(h2_lep1pt_lep1r9_plus, lep_pt[0], lep_r9[0], wgt);
      }

    }

  }

  cout << endl;
  cout << "Writing on output file" << endl;
  // delete hRecoilCorrectionTkmet;
  // delete hRecoilCorrectionPfmet;
  delete recoilCorrFile;

  // if the file is opened in UPDATE mode, the following should overwrite an object if its key inside the file already exists
  // this needs to be tested
  outputFile->Write(0,TObject::kOverwrite);

  ////////////////////////////////////////////
  // tmp plot to be removed to adjust settings in CMS_lumi     
  TH2D* h2tmp = new TH2D("h2tmp","",2,0,2,2,0,2);
  h2tmp->Fill(0.8,0.8,4);
  h2tmp->Fill(0.5,1.2,2);
  h2tmp->Fill(1.2,1.1,5);
  h2tmp->Fill(1.8,0.5,3);
  drawCorrelationPlot(h2tmp, "variable 1", "variable 2", "tmpToBeRemoved", "tmp object", outputDIR);
  system(("rm " + outputDIR + "*tmpToBeRemoved*").c_str());
  delete h2tmp;
  ////////////////////////////////////////////

  string corr_outputDIR = "";

  ////////////////////////////////////////////
  // plot correlation for combined charge
  ////////////////////////////////////////////
  // if (isMuon && (sample == Sample::wjets || sample == Sample::qcd_mu || sample == Sample::data_singleMu) ) {

  //   corr_outputDIR = outputDIR + "combined/correlation/" + sampleDir + "/";
  //   createPlotDirAndCopyPhp(corr_outputDIR);

  //   drawCorrelationPlot(h2_mT_lep1relIso04_noIsoCut, "W(#mu#nu) transverse mass [GeV]", "muon isolation (relIso04) [GeV]", "correlation_"+sampleDir+"_mT_lep1relIso04_noIsoCut",getTexLabel(sampleDir),corr_outputDIR,2);  
  //   drawCorrelationPlot(h2_mT_lep1relIso04, "W(#mu#nu) transverse mass [GeV]", "muon isolation (relIso04) [GeV]", "correlation_"+sampleDir+"_mT_lep1relIso04",getTexLabel(sampleDir),corr_outputDIR);  
  //   drawCorrelationPlot(h2_mT_lep1pt, "W(#mu#nu) transverse mass [GeV]", "muon p_{T} [GeV]", "correlation_"+sampleDir+"_mT_lep1pt",getTexLabel(sampleDir),corr_outputDIR);      
  //   drawCorrelationPlot(h2_mT_pfmet, "W(#mu#nu) transverse mass [GeV]", "PF E_{T}^{miss} [GeV]", "correlation_"+sampleDir+"_mT_pfmet",getTexLabel(sampleDir),corr_outputDIR);  
  //   drawCorrelationPlot(h2_mT_tkmet, "W(#mu#nu) transverse mass [GeV]", "tracker E_{T}^{miss} [GeV]", "correlation_"+sampleDir+"_mT_tkmet",getTexLabel(sampleDir),corr_outputDIR);  
  //   drawCorrelationPlot(h2_mT_bosonPt, "W(#mu#nu) transverse mass [GeV]", "boson p_{T} [GeV]", "correlation_"+sampleDir+"_mT_bosonPt",getTexLabel(sampleDir),corr_outputDIR);  
  //   drawCorrelationPlot(h2_lep1pt_lep1relIso04_noIsoCut, "muon p_{T} [GeV]", "muon isolation (relIso04) [GeV]", "correlation_"+sampleDir+"_lep1pt_lep1relIso04_noIsoCut",getTexLabel(sampleDir),corr_outputDIR,2);  
  //   drawCorrelationPlot(h2_lep1pt_lep1relIso04, "muon p_{T} [GeV]", "muon isolation (relIso04) [GeV]", "correlation_"+sampleDir+"_lep1pt_lep1relIso04",getTexLabel(sampleDir),corr_outputDIR);  
  //   drawCorrelationPlot(h2_lep1pt_pfmet, "muon p_{T} [GeV]", "PF E_{T}^{miss} [GeV]", "correlation_"+sampleDir+"_lep1pt_pfmet",getTexLabel(sampleDir),corr_outputDIR);  
  //   drawCorrelationPlot(h2_lep1pt_tkmet, "muon p_{T} [GeV]", "tracker E_{T}^{miss} [GeV]", "correlation_"+sampleDir+"_lep1pt_tkmet",getTexLabel(sampleDir),corr_outputDIR); 
  //   drawCorrelationPlot(h2_lep1pt_bosonPt, "muon p_{T} [GeV]", "boson p_{T} [GeV]", "correlation_"+sampleDir+"_lep1pt_bosonPt",getTexLabel(sampleDir),corr_outputDIR);  
  //   drawCorrelationPlot(h2_pfmet_tkmet, "PF E_{T}^{miss} [GeV]", "tracker E_{T}^{miss} [GeV]", "correlation_"+sampleDir+"_pfmet_tkmet",getTexLabel(sampleDir),corr_outputDIR);  


  // } else if ((not isMuon) && (sample == Sample::wjets || sample == Sample::qcd_ele || sample == Sample::data_singleEG) ) {

  //   corr_outputDIR = outputDIR + "combined/correlation/" + sampleDir + "/";
  //   createPlotDirAndCopyPhp(corr_outputDIR);

  //   drawCorrelationPlot(h2_mT_lep1relIso03_noIsoCut, "W(e#nu) transverse mass [GeV]", "electron isolation (relIso03) [GeV]", "correlation_"+sampleDir+"_mT_lep1relIso03_noIsoCut",getTexLabel(sampleDir),corr_outputDIR,2);  
  //   drawCorrelationPlot(h2_mT_lep1relIso03, "W(e#nu) transverse mass [GeV]", "electron isolation (relIso03) [GeV]", "correlation_"+sampleDir+"_mT_lep1relIso03",getTexLabel(sampleDir),corr_outputDIR);  
  //   drawCorrelationPlot(h2_mT_lep1pt, "W(e#nu) transverse mass [GeV]", "electron p_{T} [GeV]", "correlation_"+sampleDir+"_mT_lep1pt",getTexLabel(sampleDir),corr_outputDIR);      
  //   drawCorrelationPlot(h2_mT_pfmet, "W(e#nu) transverse mass [GeV]", "PF E_{T}^{miss} [GeV]", "correlation_"+sampleDir+"_mT_pfmet",getTexLabel(sampleDir),corr_outputDIR);  
  //   drawCorrelationPlot(h2_mT_tkmet, "W(e#nu) transverse mass [GeV]", "tracker E_{T}^{miss} [GeV]", "correlation_"+sampleDir+"_mT_tkmet",getTexLabel(sampleDir),corr_outputDIR);  
  //   drawCorrelationPlot(h2_mT_bosonPt, "W(e#nu) transverse mass [GeV]", "boson p_{T} [GeV]", "correlation_"+sampleDir+"_mT_bosonPt",getTexLabel(sampleDir),corr_outputDIR); 
  //   // drawCorrelationPlot(h2_mT_detaIn_noCut, "W(e#nu) transverse mass [GeV]", "electron #Delta#eta(track,SC)", "correlation_"+sampleDir+"_mT_detaIn_noCut",getTexLabel(sampleDir),corr_outputDIR);
  //   // drawCorrelationPlot(h2_mT_dphiIn_noCut, "W(e#nu) transverse mass [GeV]", "electron #Delta#phi(track,SC)", "correlation_"+sampleDir+"_mT_dphiIn_noCut",getTexLabel(sampleDir),corr_outputDIR);    
  //   // drawCorrelationPlot(h2_mT_dxy_noCut, "W(e#nu) transverse mass [GeV]", "electron track #Deltaxy", "correlation_"+sampleDir+"_mT_dxy_noCut",getTexLabel(sampleDir),corr_outputDIR);      

  //   drawCorrelationPlot(h2_lep1pt_lep1relIso04_noIsoCut, "electron p_{T} [GeV]", "electron isolation (relIso04) [GeV]", "correlation_"+sampleDir+"_lep1pt_lep1relIso04_noIsoCut",getTexLabel(sampleDir),corr_outputDIR,2);  
  //   drawCorrelationPlot(h2_lep1pt_lep1relIso04, "electron p_{T} [GeV]", "electron isolation (relIso04) [GeV]", "correlation_"+sampleDir+"_lep1pt_lep1relIso04",getTexLabel(sampleDir),corr_outputDIR);  
  //   drawCorrelationPlot(h2_lep1pt_pfmet, "electron p_{T} [GeV]", "PF E_{T}^{miss} [GeV]", "correlation_"+sampleDir+"_lep1pt_pfmet",getTexLabel(sampleDir),corr_outputDIR);  
  //   drawCorrelationPlot(h2_lep1pt_tkmet, "electron p_{T} [GeV]", "tracker E_{T}^{miss} [GeV]", "correlation_"+sampleDir+"_lep1pt_tkmet",getTexLabel(sampleDir),corr_outputDIR); 
  //   drawCorrelationPlot(h2_lep1pt_bosonPt, "electron p_{T} [GeV]", "boson p_{T} [GeV]", "correlation_"+sampleDir+"_lep1pt_bosonPt",getTexLabel(sampleDir),corr_outputDIR); 
  //   // drawCorrelationPlot(h2_lep1pt_detaIn_noCut, "electron p_{T} [GeV]", "electron #Delta#eta(track,SC)", "correlation_"+sampleDir+"_lep1pt_detaIn_noCut",getTexLabel(sampleDir),corr_outputDIR);      
  //   // drawCorrelationPlot(h2_lep1pt_dphiIn_noCut, "electron p_{T} [GeV]", "electron #Delta#eta(track,SC)", "correlation_"+sampleDir+"_lep1pt_detaIn_noCut",getTexLabel(sampleDir),corr_outputDIR);      
  //   // drawCorrelationPlot(h2_lep1pt_dxy_noCut, "electron p_{T} [GeV]", "electron track #Deltaxy", "correlation_"+sampleDir+"_lep1pt_dxy_noCut",getTexLabel(sampleDir),corr_outputDIR);      
  //   drawCorrelationPlot(h2_pfmet_tkmet, "PF E_{T}^{miss} [GeV]", "tracker E_{T}^{miss} [GeV]", "correlation_"+sampleDir+"_pfmet_tkmet",getTexLabel(sampleDir),corr_outputDIR);

  // }

  cout << endl;

}


//=============================================================

void makeWBosonVariableHistograms8TeV(const string& inputDIR = "./", const string& outputDIR = "./", 
				      const string& outfileName = "wmass_varhists.root", 
				      const Bool_t QCD_enriched_region = false, 
				      const Bool_t isMuon = false) {

  if (outputDIR!= "./") system(("mkdir -p " + outputDIR).c_str());
  system(("cp "+ PhpToCopy + " " + outputDIR).c_str());

  gROOT->SetBatch(kTRUE);
  TH1::SetDefaultSumw2(); //all the following histograms will automatically call TH1::Sumw2()  

  cout << endl;

  //TFile* outputFile = new TFile((outputDIR + outfileName).c_str(),"RECREATE");
  TFile* outputFile = new TFile((outputDIR + outfileName).c_str(),"UPDATE");
  if (!outputFile || outputFile->IsZombie()) {
    cout << "Error: file not opened. Exit" << endl;
    exit(EXIT_FAILURE);
  }

  outputFile->cd();

  if (isMuon) {
    fillHistograms(inputDIR, outputDIR, Sample::data_singleMu, outputFile, QCD_enriched_region, isMuon);
    fillHistograms(inputDIR, outputDIR, Sample::qcd_mu, outputFile, QCD_enriched_region, isMuon);
  } else {
    fillHistograms(inputDIR, outputDIR, Sample::data_singleEG, outputFile, QCD_enriched_region, isMuon);
    // fillHistograms(inputDIR, outputDIR, Sample::qcd_ele, outputFile, QCD_enriched_region, isMuon);
  }
  fillHistograms(inputDIR, outputDIR, Sample::wjets, outputFile, QCD_enriched_region, isMuon);
  fillHistograms(inputDIR, outputDIR, Sample::zjets, outputFile, QCD_enriched_region, isMuon);
  fillHistograms(inputDIR, outputDIR, Sample::top, outputFile, QCD_enriched_region, isMuon);
  fillHistograms(inputDIR, outputDIR, Sample::diboson, outputFile, QCD_enriched_region, isMuon);
  
  outputFile->Close();
  delete outputFile;

  system(("cp makeWBosonVariableHistograms8TeV.C " + outputDIR).c_str());
  
}
