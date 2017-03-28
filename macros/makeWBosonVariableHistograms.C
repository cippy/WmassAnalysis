#include "utility.h"

using namespace std;

//static Double_t intLumi = 36.4;

static Int_t nMtBins = 120;
static Double_t mtMin = 20;
static Double_t mtMax = 140;

// static Bool_t useTrackMet = true;
// static Bool_t useAbsIso = true;  // rel iso or abs iso to cut (abs iso is rel iso times lepton pT) 

//=============================================================

void fillHistograms(const string& inputDIR = "./", const string& outputDIR = "./", 
		    const Sample& sample = Sample::zjets, 
		    TFile* outputFile = NULL, 
		    const Bool_t QCD_enriched_region = false, const Bool_t isMuon = false) {

  if (outputFile == NULL) {
    cout << "Error: file is NULL. please check. Exit ..." << endl;
    exit(EXIT_FAILURE);
  }

  outputFile->cd();

  cout << "useAbsIso = " << useAbsIso << endl;

  Double_t eleIso03thrEB = 0.0588; // 0.0588, 0.075, 0.05
  Double_t eleIso03thrEE = 0.0571; // 0.0571, 0.075, 0.05
  Double_t muIso04thr = 0.15;       // 0.15,  0.2,   0.125
  Double_t relToAbsIsoFactor = 40; // 40
  Double_t eleIso03thr_QCD = 0.0588;
  Double_t muIso04thr_QCD = 0.15;       // 0.15,  e.g., sig region for iso < 0.15 and QCD region for iso > 0.2

  Double_t miniRelIsoThr = 0.08;
  Double_t miniRelIsoThr_QCD = 0.08;
 
  // be consistent with threshold when assigning histogram boundaries

  Int_t nBins_lepIso03 = 20;  // 20
  Double_t min_lepIso03 = 0.0; // 0.0
  Double_t max_lepIso03 = 0.060; // 0.060
  Int_t nBins_lepIso04 = 30;  // 30
  Double_t min_lepIso04 = 0.0; // 0.0
  Double_t max_lepIso04 = 0.15; // 0.150

  Int_t nBins_miniRelIso = 16;  // 30
  Double_t min_miniRelIso = 0.0; // 0.0
  Double_t max_miniRelIso = 0.08; // 0.150

  if (QCD_enriched_region) {
    nBins_lepIso03 = 90;  // 90, 85
    min_lepIso03 = 0.05; // 0.05, 0.075
    max_lepIso03 = 0.5; // 0.50
    nBins_lepIso04 = 35;  // 35,  30
    min_lepIso04 = 0.15; // 0.15, 0.2
    max_lepIso04 = 0.500; // 0.500
    nBins_miniRelIso = 34;  // 35,  30
    min_miniRelIso = 0.08; // 0.15, 0.2
    max_miniRelIso = 0.25; // 0.500
  }

  if (useAbsIso) {
    // absIso = relIso * 40 GeV (40 GeV is the average pT for a lepton from W, thresholds to be used for selection must be taken accordingly)
    min_lepIso03 *= relToAbsIsoFactor; // 0.0
    max_lepIso03 *= relToAbsIsoFactor; // 0.060
    min_lepIso04 *= relToAbsIsoFactor; // 0.0
    max_lepIso04 *= relToAbsIsoFactor; // 0.150

  }

  Int_t chargedLeptonFlavour = isMuon ? 13 : 11;
  Double_t lepEtaThreshold = isMuon ? 2.4 : 1.479;  // EB only for electron for now
  Int_t lepTightIdThreshold = isMuon ? 1 : 3;

  TDirectory *dirSample = NULL;
  if (sample == Sample::data_doubleEG)
    dirSample = outputFile->mkdir("data_doubleEG");
  else if (sample == Sample::data_singleEG)
    dirSample = outputFile->mkdir("data_singleEG");
  if (sample == Sample::data_doubleMu)
    dirSample = outputFile->mkdir("data_doubleMu");
  else if (sample == Sample::data_singleMu)
    dirSample = outputFile->mkdir("data_singleMu");
  else if (sample == Sample::zjets) 
    dirSample = outputFile->mkdir("zjets");
  else if (sample == Sample::wjets)
    dirSample = outputFile->mkdir("wjets");
  else if (sample == Sample::qcd_mu)
    dirSample = outputFile->mkdir("qcd_mu");
  else if (sample == Sample::qcd_ele)
    dirSample = outputFile->mkdir("qcd_ele");
  else if (sample == Sample::top)
    dirSample = outputFile->mkdir("top");
  else if (sample == Sample::diboson)
    dirSample = outputFile->mkdir("diboson");

  dirSample->cd();

  gROOT->SetBatch(kTRUE);
  TH1::SetDefaultSumw2(); //all the following histograms will automatically call TH1::Sumw2()  

  cout << endl;

  TChain* chain = new TChain("tree");
  TChain* chFriend = new TChain("mjvars/t");

  vector<Double_t> genwgtVec;
  buildChain(chain, chFriend, inputDIR, sample);

  TTreeReader reader (chain);

  // reco met
  TTreeReaderValue<Float_t> tkmet    (reader,"met_trkPt");
  TTreeReaderValue<Float_t> tkmet_phi(reader,"met_trkPhi");
  TTreeReaderValue<Float_t> pfmet    (reader,"met_pt");
  TTreeReaderValue<Float_t> pfmet_phi(reader,"met_phi");

  // lepGood branch
  TTreeReaderValue<Int_t> nlep  (reader,"nLepGood");
  TTreeReaderArray<Int_t> lep_pdgId (reader,"LepGood_pdgId");
  TTreeReaderArray<Float_t> lep_pt (reader,"LepGood_pt");
  TTreeReaderArray<Float_t> lep_eta (reader,"LepGood_eta");
  TTreeReaderArray<Float_t> lep_phi (reader,"LepGood_phi");
  TTreeReaderArray<Float_t> lep_mass (reader,"LepGood_mass");
  TTreeReaderArray<Float_t> lep_relIso03 (reader,"LepGood_relIso03");
  TTreeReaderArray<Float_t> lep_relIso04 (reader,"LepGood_relIso04");
  TTreeReaderArray<Int_t> lep_tightId (reader,"LepGood_tightId");
  TTreeReaderArray<Float_t> lep_etaSc (reader,"LepGood_etaSc");
  TTreeReaderArray<Float_t> lep_r9 (reader,"LepGood_full5x5_r9");
  TTreeReaderArray<Float_t> lep_sigmaIetaIeta (reader,"LepGood_full5x5_sigmaIetaIeta");
  TTreeReaderArray<Int_t> lep_hltSafeID (reader,"LepGood_hltId");

  TTreeReaderArray<Float_t> lep_miniRelIso(reader,"LepGood_miniRelIso");

  // MC reweight
  TTreeReaderValue<Float_t> weight(reader,"weight");
  TTreeReaderValue<Float_t> puw(reader,"puw");

  // TH1
  TH1D* hmT = new TH1D("hmT","",nMtBins, mtMin, mtMax);
  TH1D* hpfmet = new TH1D("hpfmet","",60,0,300);
  TH1D* htkmet = new TH1D("htkmet","",60,0,300);
  TH1D* hlep1pt = new TH1D("hlep1pt","",28,24,80);  //60, 0, 300
  TH1D* hlep2pt = new TH1D("hlep2pt","",40,0,80);
  TH1D* hbosonpt = new TH1D("hbosonpt","",30,0,30);
  TH1D* hbosonpt_wlike = new TH1D("hbosonpt_wlike","",30,0,30);
  TH1D* hbosoneta = new TH1D("hbosoneta","",100,-5,5);
  TH1D* hlep1sigIetaIeta = new TH1D("hlep1sigIetaIeta","",25, 0.0, 0.025);
  TH1D* hlep1relIso03 = new TH1D("hlep1relIso03","",nBins_lepIso03,min_lepIso03,max_lepIso03);
  TH1D* hlep1relIso04 = new TH1D("hlep1relIso04","",nBins_lepIso04,min_lepIso04,max_lepIso04);
  TH1D* hlep1miniRelIso = new TH1D("hlep1miniRelIso","",nBins_miniRelIso, min_miniRelIso, max_miniRelIso);
  TH1D* hlep1r9 = new TH1D("hlep1r9","", 55, 0.0, 1.1);
  TH1D* hrecoil = new TH1D("hrecoil","",30,0,30);

  TH2D* h2_mT_lep1pt = new TH2D("h2_mT_lep1pt","",nMtBins, mtMin, mtMax, 28,24,80);
  TH2D* h2_mT_lep1eta = new TH2D("h2_mT_lep1eta","",nMtBins, mtMin, mtMax, 30, -1.5, 1.5);
  TH2D* h2_mT_lep1sigIetaIeta = new TH2D("h2_mT_lep1sigIetaIeta","",nMtBins, mtMin, mtMax, 25, 0.0, 0.025);
  TH2D* h2_mT_lep1r9 = new TH2D("h2_mT_lep1r9","",nMtBins, mtMin, mtMax, 55, 0.0, 1.1);
  TH2D* h2_mT_pfmet = new TH2D("h2_mT_pfmet","",nMtBins, mtMin, mtMax, 60,0,300);
  TH2D* h2_mT_tkmet = new TH2D("h2_mT_tkmet","",nMtBins, mtMin, mtMax, 60,0,300);
  TH2D* h2_mT_bosonPt = new TH2D("h2_mT_bosonPt","",nMtBins, mtMin, mtMax, 30,0,30);
  TH2D* h2_mT_lep1relIso03 = new TH2D("h2_mT_lep1relIso03","",nMtBins, mtMin, mtMax, nBins_lepIso03,min_lepIso03,max_lepIso03);
  TH2D* h2_mT_lep1relIso04 = new TH2D("h2_mT_lep1relIso04","",nMtBins, mtMin, mtMax, nBins_lepIso04,min_lepIso04,max_lepIso04);

  TH2D* h2_lep1pt_lep1eta = new TH2D("h2_lep1pt_lep1eta","",28,24,80, 30, -1.5, 1.5);
  TH2D* h2_lep1pt_lep1sigIetaIeta = new TH2D("h2_lep1pt_lep1sigIetaIeta","",28,24,80, 25, 0.0, 0.025);
  TH2D* h2_lep1pt_lep1r9 = new TH2D("h2_lep1pt_lep1r9","",28,24,80, 55, 0.0, 1.1);
  TH2D* h2_lep1pt_pfmet = new TH2D("h2_lep1pt_pfmet","",28,24,80, 60,0,300);
  TH2D* h2_lep1pt_tkmet = new TH2D("h2_lep1pt_tkmet","",28,24,80, 60,0,300);
  TH2D* h2_lep1pt_bosonPt = new TH2D("h2_lep1pt_bosonPt","",28,24,80, 30,0,30);
  TH2D* h2_lep1pt_lep1relIso03 = new TH2D("h2_lep1pt_lep1relIso03","",28,24,80, nBins_lepIso03,min_lepIso03,max_lepIso03);
  TH2D* h2_lep1pt_lep1relIso04 = new TH2D("h2_lep1pt_lep1relIso04","",28,24,80, nBins_lepIso04,min_lepIso04,max_lepIso04);

  TH1D* hlep1relIso03_noIsoCut = new TH1D("hlep1relIso03_noIsoCut","",100,0.0,0.5);
  TH2D* h2_mT_lep1relIso03_noIsoCut = new TH2D("h2_mT_lep1relIso03_noIsoCut","",nMtBins, mtMin, mtMax,100,0.0,0.5);
  TH2D* h2_lep1pt_lep1relIso03_noIsoCut = new TH2D("h2_lep1pt_lep1relIso03_noIsoCut","",28,24,80, 100,0.0,0.5);
  TH1D* hlep1relIso04_noIsoCut = new TH1D("hlep1relIso04_noIsoCut","",100,0.0,0.5);
  TH2D* h2_mT_lep1relIso04_noIsoCut = new TH2D("h2_mT_lep1relIso04_noIsoCut","",nMtBins, mtMin, mtMax,100,0.0,0.5);
  TH2D* h2_lep1pt_lep1relIso04_noIsoCut = new TH2D("h2_lep1pt_lep1relIso04_noIsoCut","",28,24,80, 100,0.0,0.5);

  TH1D* hlep1miniRelIso_noIsoCut = new TH1D("hlep1miniRelIso_noIsoCut","",50,0.0,0.25);
  TH2D* h2_lep1pt_lep1miniRelIso_noIsoCut = new TH2D("h2_lep1pt_lep1miniRelIso_noIsoCut","",28,24,80, 50, 0.0, 0.25);
  TH2D* h2_mT_lep1miniRelIso_noIsoCut = new TH2D("h2_mT_lep1miniRelIso_noIsoCut","",nMtBins, mtMin, mtMax, 50, 0.0, 0.25);


  /////////////////////////////////
  // only for positive lepton
  TH1D* hmT_plus = new TH1D("hmT_plus","",nMtBins, mtMin, mtMax);
  TH1D* hpfmet_plus = new TH1D("hpfmet_plus","",60,0,300);
  TH1D* htkmet_plus = new TH1D("htkmet_plus","",60,0,300);
  TH1D* hlep1pt_plus = new TH1D("hlep1pt_plus","",28,24,80);  //60, 0, 300
  TH1D* hlep2pt_plus = new TH1D("hlep2pt_plus","",35,10,80);
  TH1D* hbosonpt_plus = new TH1D("hbosonpt_plus","",30,0,30);
  TH1D* hbosonpt_wlike_plus = new TH1D("hbosonpt_wlike_plus","",30,0,30);
  TH1D* hbosoneta_plus = new TH1D("hbosoneta_plus","",100,-5,5);
  TH1D* hlep1sigIetaIeta_plus = new TH1D("hlep1sigIetaIeta_plus","",25, 0.0, 0.025);
  TH1D* hlep1relIso03_plus = new TH1D("hlep1relIso03_plus","",nBins_lepIso03,min_lepIso03,max_lepIso03);
  TH1D* hlep1relIso04_plus = new TH1D("hlep1relIso04_plus","",nBins_lepIso04,min_lepIso04,max_lepIso04);
  TH1D* hlep1miniRelIso_plus = new TH1D("hlep1miniRelIso_plus","",nBins_miniRelIso, min_miniRelIso, max_miniRelIso);
  TH1D* hlep1r9_plus = new TH1D("hlep1r9_plus","", 55, 0.0, 1.1);
  TH1D* hrecoil_plus = new TH1D("hrecoil_plus","",30,0,30);

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
  TH1D* hpfmet_minus = new TH1D("hpfmet_minus","",60,0,300);
  TH1D* htkmet_minus = new TH1D("htkmet_minus","",60,0,300);
  TH1D* hlep1pt_minus = new TH1D("hlep1pt_minus","",28,24,80);  //60, 0, 300
  TH1D* hlep2pt_minus = new TH1D("hlep2pt_minus","",35,10,80);
  TH1D* hbosonpt_minus = new TH1D("hbosonpt_minus","",30,0,30);
  TH1D* hbosonpt_wlike_minus = new TH1D("hbosonpt_wlike_minus","",30,0,30);
  TH1D* hbosoneta_minus = new TH1D("hbosoneta_minus","",100,-5,5);
  TH1D* hlep1sigIetaIeta_minus = new TH1D("hlep1sigIetaIeta_minus","",25, 0.0, 0.025);
  TH1D* hlep1relIso03_minus = new TH1D("hlep1relIso03_minus","",nBins_lepIso03,min_lepIso03,max_lepIso03);
  TH1D* hlep1relIso04_minus = new TH1D("hlep1relIso04_minus","",nBins_lepIso04,min_lepIso04,max_lepIso04);
  TH1D* hlep1miniRelIso_minus = new TH1D("hlep1miniRelIso_minus","",nBins_miniRelIso, min_miniRelIso, max_miniRelIso);
  TH1D* hlep1r9_minus = new TH1D("hlep1r9_minus","", 55, 0.0, 1.1);
  TH1D* hrecoil_minus = new TH1D("hrecoil_minus","",30,0,30);

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

  // start event loop                                                                                                                                                       
  long int nTotal = chain->GetEntries();
  long int nEvents = 0;

  Double_t metToUse_pt = 0.0;
  Double_t metToUse_phi = 0.0;

  Bool_t negativeLeptonHasPassedSelection = false;
  Bool_t positiveLeptonHasPassedSelection = false;
  
  Double_t wgt = 1.0;

  while(reader.Next()){
  
    cout.flush();
    if(nEvents % 100000 == 0) cout<<"\r"<<"Analyzing events "<<double(nEvents)/nTotal*100<<" % ";
    nEvents++;

    if (useTrackMet) {
      metToUse_pt = *tkmet;
      metToUse_phi = *tkmet_phi;
    } else {
      metToUse_pt = *pfmet;
      metToUse_phi = *pfmet_phi;
    }

    negativeLeptonHasPassedSelection = false;
    positiveLeptonHasPassedSelection = false;

    // selection
    if (*nlep != 1) continue;    // 1 leptons                                                  
    if (fabs(lep_pdgId[0]) != chargedLeptonFlavour ) continue;  // electrons                                                     
    if ( lep_pt[0] < 30.0) continue;
    if ( fabs(lep_eta[0]) > lepEtaThreshold ) continue; 
    if ( lep_tightId[0] < lepTightIdThreshold ) continue;   // tight ID
    if (not isMuon) {
      if ( lep_hltSafeID[0] != 1 ) continue;  //HLT safe ID for electrons
    }
    //if ( *pfmet < 30) continue;    

    if (lep_pdgId[0] > 0) negativeLeptonHasPassedSelection = true;
    else                                   positiveLeptonHasPassedSelection = true;

    if (sample == Sample::data_doubleEG) wgt = 1.0;
    else wgt = intLumi * *puw * *weight;

    TLorentzVector lep1Reco, lep2Reco, bosonReco, metReco, metWlikeReco, bosonRecoWlike;

    lep1Reco.SetPtEtaPhiM(lep_pt[0],lep_eta[0],lep_phi[0],lep_mass[0]); 
    metReco.SetPtEtaPhiM(metToUse_pt,0,metToUse_phi,0);
    metWlikeReco = metReco;
    lep2Reco = metWlikeReco;

    Double_t leptonIso03toUse = useAbsIso ? (lep_relIso03[0] * lep1Reco.Pt()) : lep_relIso03[0];
    Double_t leptonIso04toUse = useAbsIso ? (lep_relIso04[0] * lep1Reco.Pt()) : lep_relIso04[0];
    Double_t leptonIso03threshold = (fabs(lep_etaSc[0]) < 1.479) ? eleIso03thrEB : eleIso03thrEE;
    Double_t leptonIso04threshold = muIso04thr;
    Double_t leptonIso03threshold_QCD = eleIso03thr_QCD;
    Double_t leptonIso04threshold_QCD = muIso04thr_QCD;
    if (useAbsIso) {
      // multiply by 40 (average lepton pT in GeV, when lepton comes from a W decay almost at rest)
      leptonIso03threshold = leptonIso03threshold * relToAbsIsoFactor;
      leptonIso04threshold = leptonIso04threshold * relToAbsIsoFactor;
      leptonIso03threshold_QCD = leptonIso03threshold_QCD * relToAbsIsoFactor;
      leptonIso04threshold_QCD = leptonIso04threshold_QCD * relToAbsIsoFactor;
    }
    
    bosonReco = lep1Reco + metWlikeReco;
    bosonRecoWlike = lep1Reco + metWlikeReco;
    if (bosonReco.Pt() > 30) continue;

    Double_t mT = sqrt(2. * lep1Reco.Pt() * metWlikeReco.Pt() * (1. - cos(lep1Reco.DeltaPhi(metWlikeReco)) ) );
    if (mT < mtMin) continue;
    //if (mT > mtMax) continue;    

    TVector2 recoilReco = metWlikeReco.Vect().XYvector() + lep1Reco.Vect().XYvector();
    //   if (recoilReco.Mod() > 30) continue;

    // fill histograms for isolation before cutting on it
    fillTH1(hlep1relIso03_noIsoCut, leptonIso03toUse, wgt);
    fillTH2(h2_mT_lep1relIso03_noIsoCut, mT, leptonIso03toUse, wgt);
    fillTH2(h2_mT_lep1relIso03_noIsoCut, lep1Reco.Pt(), leptonIso03toUse, wgt);
    fillTH1(hlep1relIso04_noIsoCut, leptonIso04toUse, wgt);
    fillTH2(h2_mT_lep1relIso04_noIsoCut, mT, leptonIso04toUse, wgt);
    fillTH2(h2_lep1pt_lep1relIso04_noIsoCut, lep1Reco.Pt(), leptonIso04toUse, wgt);

    fillTH1(hlep1miniRelIso_noIsoCut, lep_miniRelIso[0], wgt);
    fillTH2(h2_mT_lep1miniRelIso_noIsoCut, mT, lep_miniRelIso[0], wgt);
    fillTH2(h2_lep1pt_lep1miniRelIso_noIsoCut, lep1Reco.Pt(), lep_miniRelIso[0], wgt);

    if (negativeLeptonHasPassedSelection) {
      fillTH1(hlep1relIso03_noIsoCut_minus, leptonIso03toUse, wgt);
      fillTH2(h2_mT_lep1relIso03_noIsoCut_minus, mT, leptonIso03toUse, wgt);
      fillTH2(h2_lep1pt_lep1relIso03_noIsoCut_minus, lep1Reco.Pt(), leptonIso03toUse, wgt);
      fillTH1(hlep1relIso04_noIsoCut_minus, leptonIso04toUse, wgt);
      fillTH2(h2_mT_lep1relIso04_noIsoCut_minus, mT, leptonIso04toUse, wgt);
      fillTH2(h2_lep1pt_lep1relIso04_noIsoCut_minus, lep1Reco.Pt(), leptonIso04toUse, wgt);
    }

    if (positiveLeptonHasPassedSelection) {
      fillTH1(hlep1relIso03_noIsoCut_plus, leptonIso03toUse, wgt);
      fillTH2(h2_mT_lep1relIso03_noIsoCut_plus, mT, leptonIso03toUse, wgt);
      fillTH2(h2_lep1pt_lep1relIso03_noIsoCut_plus, lep1Reco.Pt(), leptonIso03toUse, wgt);
      fillTH1(hlep1relIso04_noIsoCut_plus, leptonIso04toUse, wgt);
      fillTH2(h2_mT_lep1relIso04_noIsoCut_plus, mT, leptonIso04toUse, wgt);
      fillTH2(h2_lep1pt_lep1relIso04_noIsoCut_plus, lep1Reco.Pt(), leptonIso04toUse, wgt);
    }

    ////////////////////////////
    // LEP ISO CUT
    /////////////////////////////
    if (isMuon) {
      
      if (QCD_enriched_region) {
	//if ( leptonIso04toUse < leptonIso04threshold_QCD) continue;  // tight iso first lep                 
	if ( lep_miniRelIso[0] < miniRelIsoThr_QCD) continue;  // tight iso first lep                 
      } else {
	//if ( leptonIso04toUse > leptonIso04threshold) continue;  // tight iso first lep                                  
	if ( lep_miniRelIso[0] > miniRelIsoThr) continue;  // tight iso first lep          	
      }
      
    } else {
      
      if (QCD_enriched_region) {	  
	if ( leptonIso03toUse < leptonIso03threshold_QCD) continue;  // tight iso first lep
      } else {
	if ( leptonIso03toUse > leptonIso03threshold) continue;  // tight iso first lep
      }
      
    }
    
    fillTH1(hmT,(Double_t) mT, wgt);
    fillTH1(hpfmet,(Double_t) *pfmet, wgt);
    fillTH1(htkmet,(Double_t) *tkmet, wgt);
    fillTH1(hlep1pt,(Double_t) lep1Reco.Pt(), wgt);
    fillTH1(hlep2pt,(Double_t) metWlikeReco.Pt(), wgt);
    fillTH1(hbosonpt,(Double_t) bosonReco.Pt(), wgt);
    fillTH1(hbosonpt_wlike,(Double_t) bosonRecoWlike.Pt(), wgt);
    fillTH1(hbosoneta,(Double_t) bosonReco.Eta(), wgt);
    fillTH1(hlep1sigIetaIeta,(Double_t) lep_sigmaIetaIeta[0], wgt);
    fillTH1(hlep1relIso03,(Double_t) leptonIso03toUse, wgt);
    fillTH1(hlep1relIso04,(Double_t) leptonIso04toUse, wgt);
    fillTH1(hlep1miniRelIso,(Double_t) lep_miniRelIso[0], wgt);
    fillTH1(hlep1r9,(Double_t) lep_r9[0], wgt);
    fillTH1(hrecoil, recoilReco.Mod(), wgt);

    fillTH2(h2_mT_lep1pt, mT, lep_pt[0], wgt);
    fillTH2(h2_mT_lep1eta, mT, lep_eta[0], wgt);
    fillTH2(h2_mT_lep1sigIetaIeta, mT, lep_sigmaIetaIeta[0], wgt);
    fillTH2(h2_mT_lep1r9, mT, lep_r9[0], wgt);
    fillTH2(h2_mT_pfmet, mT, *pfmet, wgt);
    fillTH2(h2_mT_tkmet, mT, *tkmet, wgt);
    fillTH2(h2_mT_bosonPt, mT, bosonReco.Pt(), wgt);
    fillTH2(h2_mT_lep1relIso03, mT, leptonIso03toUse, wgt);
    fillTH2(h2_mT_lep1relIso04, mT, leptonIso04toUse, wgt);

    fillTH2(h2_lep1pt_lep1eta, lep_pt[0], lep_eta[0], wgt);
    fillTH2(h2_lep1pt_lep1sigIetaIeta, lep_pt[0], lep_sigmaIetaIeta[0], wgt);
    fillTH2(h2_lep1pt_lep1r9, lep_pt[0], lep_r9[0], wgt);
    fillTH2(h2_lep1pt_pfmet, lep_pt[0], *pfmet, wgt);
    fillTH2(h2_lep1pt_tkmet, lep_pt[0], *tkmet, wgt);
    fillTH2(h2_lep1pt_bosonPt, lep_pt[0], bosonReco.Pt(), wgt);
    fillTH2(h2_lep1pt_lep1relIso03, lep_pt[0], leptonIso03toUse, wgt);
    fillTH2(h2_lep1pt_lep1relIso04, lep_pt[0], leptonIso04toUse, wgt);

    if (negativeLeptonHasPassedSelection) {
      fillTH1(hmT_minus,(Double_t) mT, wgt);
      fillTH1(hpfmet_minus,(Double_t) *pfmet, wgt);
      fillTH1(htkmet_minus,(Double_t) *tkmet, wgt);
      fillTH1(hlep1pt_minus,(Double_t) lep1Reco.Pt(), wgt);
      fillTH1(hlep2pt_minus,(Double_t) metWlikeReco.Pt(), wgt);
      fillTH1(hbosonpt_minus,(Double_t) bosonReco.Pt(), wgt);
      fillTH1(hbosonpt_wlike_minus,(Double_t) bosonRecoWlike.Pt(), wgt);
      fillTH1(hbosoneta_minus,(Double_t) bosonReco.Eta(), wgt);
      fillTH1(hlep1sigIetaIeta_minus,(Double_t) lep_sigmaIetaIeta[0], wgt);
      fillTH1(hlep1relIso03_minus,(Double_t) leptonIso03toUse, wgt);
      fillTH1(hlep1relIso04_minus,(Double_t) leptonIso04toUse, wgt);
      fillTH1(hlep1miniRelIso_minus,(Double_t) lep_miniRelIso[0], wgt);
      fillTH1(hlep1r9_minus,(Double_t) lep_r9[0], wgt);
      fillTH1(hrecoil_minus, recoilReco.Mod(), wgt);

      fillTH2(h2_mT_lep1pt_minus, mT, lep_pt[0], wgt);
      fillTH2(h2_mT_lep1eta_minus, mT, lep_eta[0], wgt);
      fillTH2(h2_mT_lep1sigIetaIeta_minus, mT, lep_sigmaIetaIeta[0], wgt);
      fillTH2(h2_mT_lep1r9_minus, mT, lep_r9[0], wgt);
      fillTH2(h2_mT_pfmet_minus, mT, *pfmet, wgt);
      fillTH2(h2_mT_tkmet_minus, mT, *tkmet, wgt);
      fillTH2(h2_mT_bosonPt_minus, mT, bosonReco.Pt(), wgt);
      fillTH2(h2_mT_lep1relIso03_minus, mT, leptonIso03toUse, wgt);
      fillTH2(h2_mT_lep1relIso04_minus, mT, leptonIso04toUse, wgt);

      fillTH2(h2_lep1pt_lep1eta_minus, lep_pt[0], lep_eta[0], wgt);
      fillTH2(h2_lep1pt_lep1sigIetaIeta_minus, lep_pt[0], lep_sigmaIetaIeta[0], wgt);
      fillTH2(h2_lep1pt_lep1r9_minus, lep_pt[0], lep_r9[0], wgt);
      fillTH2(h2_lep1pt_pfmet_minus, lep_pt[0], *pfmet, wgt);
      fillTH2(h2_lep1pt_tkmet_minus, lep_pt[0], *tkmet, wgt);
      fillTH2(h2_lep1pt_bosonPt_minus, lep_pt[0], bosonReco.Pt(), wgt);
      fillTH2(h2_lep1pt_lep1relIso03_minus, lep_pt[0], leptonIso03toUse, wgt);
      fillTH2(h2_lep1pt_lep1relIso04_minus, lep_pt[0], leptonIso04toUse, wgt);

    }

    if (positiveLeptonHasPassedSelection) {
      fillTH1(hmT_plus,(Double_t) mT, wgt);
      fillTH1(hpfmet_plus,(Double_t) *pfmet, wgt);
      fillTH1(htkmet_plus,(Double_t) *tkmet, wgt);
      fillTH1(hlep1pt_plus,(Double_t) lep1Reco.Pt(), wgt);
      fillTH1(hlep2pt_plus,(Double_t) metWlikeReco.Pt(), wgt);
      fillTH1(hbosonpt_plus,(Double_t) bosonReco.Pt(), wgt);
      fillTH1(hbosonpt_wlike_plus,(Double_t) bosonRecoWlike.Pt(), wgt);
      fillTH1(hbosoneta_plus,(Double_t) bosonReco.Eta(), wgt);
      fillTH1(hlep1sigIetaIeta_plus,(Double_t) lep_sigmaIetaIeta[0], wgt);
      fillTH1(hlep1relIso03_plus,(Double_t) leptonIso03toUse, wgt);
      fillTH1(hlep1relIso04_plus,(Double_t) leptonIso04toUse, wgt);
      fillTH1(hlep1miniRelIso_plus,(Double_t) lep_miniRelIso[0], wgt);
      fillTH1(hlep1r9_plus,(Double_t) lep_r9[0], wgt);
      fillTH1(hrecoil_plus, recoilReco.Mod(), wgt);

      fillTH2(h2_mT_lep1pt_plus, mT, lep_pt[0], wgt);
      fillTH2(h2_mT_lep1eta_plus, mT, lep_eta[0], wgt);
      fillTH2(h2_mT_lep1sigIetaIeta_plus, mT, lep_sigmaIetaIeta[0], wgt);
      fillTH2(h2_mT_lep1r9_plus, mT, lep_r9[0], wgt);
      fillTH2(h2_mT_pfmet_plus, mT, *pfmet, wgt);
      fillTH2(h2_mT_tkmet_plus, mT, *tkmet, wgt);
      fillTH2(h2_mT_bosonPt_plus, mT, bosonReco.Pt(), wgt);
      fillTH2(h2_mT_lep1relIso03_plus, mT, leptonIso03toUse, wgt);
      fillTH2(h2_mT_lep1relIso04_plus, mT, leptonIso04toUse, wgt);

      fillTH2(h2_lep1pt_lep1eta_plus, lep_pt[0], lep_eta[0], wgt);
      fillTH2(h2_lep1pt_lep1sigIetaIeta_plus, lep_pt[0], lep_sigmaIetaIeta[0], wgt);
      fillTH2(h2_lep1pt_lep1r9_plus, lep_pt[0], lep_r9[0], wgt);
      fillTH2(h2_lep1pt_pfmet_plus, lep_pt[0], *pfmet, wgt);
      fillTH2(h2_lep1pt_tkmet_plus, lep_pt[0], *tkmet, wgt);
      fillTH2(h2_lep1pt_bosonPt_plus, lep_pt[0], bosonReco.Pt(), wgt);
      fillTH2(h2_lep1pt_lep1relIso03_plus, lep_pt[0], leptonIso03toUse, wgt);
      fillTH2(h2_lep1pt_lep1relIso04_plus, lep_pt[0], leptonIso04toUse, wgt);

    }



  }

  cout << endl;
  cout << "Writing on output file" << endl;
  outputFile->Write();

  cout << endl;

}


//=============================================================

void makeWBosonVariableHistograms(const string& inputDIR = "./", const string& outputDIR = "./", 
				  const string& outfileName = "wmass_varhists.root", 
				  const Bool_t QCD_enriched_region = false, 
				  const Bool_t isMuon = false) {

  if (outputDIR!= "./") system(("mkdir -p " + outputDIR).c_str());
  system(("cp "+ PhpToCopy + " " + outputDIR).c_str());

  gROOT->SetBatch(kTRUE);
  TH1::SetDefaultSumw2(); //all the following histograms will automatically call TH1::Sumw2()  

  cout << endl;

  TFile* outputFile = new TFile((outputDIR + outfileName).c_str(),"RECREATE");
  if (!outputFile || outputFile->IsZombie()) {
    cout << "Error: file not opened. Exit" << endl;
    exit(EXIT_FAILURE);
  }

  outputFile->cd();

  fillHistograms(inputDIR, outputDIR, Sample::wjets, outputFile, QCD_enriched_region, isMuon);
  if (isMuon) fillHistograms(inputDIR, outputDIR, Sample::qcd_mu, outputFile, QCD_enriched_region, isMuon);
  else fillHistograms(inputDIR, outputDIR, Sample::qcd_ele, outputFile, QCD_enriched_region, isMuon);
  fillHistograms(inputDIR, outputDIR, Sample::zjets, outputFile, QCD_enriched_region, isMuon);

  outputFile->Close();

  
}
