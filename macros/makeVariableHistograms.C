#include "utility.h"

using namespace std;

static Double_t intLumi = 36.4;

static Int_t nZmassBins = 61;
static Double_t zMassMin = 59.5;
static Double_t zMassMax = 120.5;

static Int_t nMtBins = 120;
static Double_t mtMin = 20;
static Double_t mtMax = 140;

static Bool_t useTrackMet = true;

//=============================================================

void fillHistograms(const string& inputDIR = "./", const string& outputDIR = "./", 
		    const Sample& sample = Sample::zjets, 
		    TFile* outputFile = NULL, 
		    const Bool_t QCD_enriched_region = false, const Bool_t isWregion = false, const Bool_t isPositiveLepton = false, const Bool_t isMuon = false) {

  if (outputFile == NULL) {
    cout << "Error: file is NULL. please check. Exit ..." << endl;
    exit(EXIT_FAILURE);
  }

  outputFile->cd();

  Int_t nBins_lepIso03 = 20;  // 20
  Double_t min_lepIso03 = 0.0; // 0.0
  Double_t max_lepIso03 = 0.060; // 0.060
  Int_t nBins_lepIso04 = 30;  // 20
  Double_t min_lepIso04 = 0.0; // 0.0
  Double_t max_lepIso04 = 0.150; // 0.060

  if (QCD_enriched_region) {
    nBins_lepIso03 = 45;  // 20
    min_lepIso03 = 0.050; // 0.0
    max_lepIso03 = 0.500; // 0.060
    nBins_lepIso04 = 35;  // 30
    min_lepIso04 = 0.150; // 0.0
    max_lepIso04 = 0.500; // 0.150
  }

  Int_t chargedLeptonFlavour = isMuon ? 13 : 11;
  Double_t lepEtaThreshold = isMuon ? 2.4 : 1.479;  // EB only for electron for now
  Int_t lepTightIdThreshold = isMuon ? 1 : 3;

  // the following candidateLep_index is an index used to select the lepton in the LepGood array. This is used as the charged lepton candidate from W or Z(W-like).
  // for W, it is automatically set to 0 (first one, which should be the only one), while otherLep_index is not used
  // for Z, it is the lepton used to mimic the charged lepton in the W-like analysis, otherLep_index is the other one, with looser selection
  Int_t candidateLep_index = -1;
  Int_t otherLep_index = -1;

  TDirectory *dirSample = NULL;
  if (sample == Sample::data_doubleEG)
    dirSample = outputFile->mkdir("data_doubleEG");
  else if (sample == Sample::data_singleEG)
    dirSample = outputFile->mkdir("data_singleEG");
  else if (sample == Sample::zjets) 
    dirSample = outputFile->mkdir("zjets");
  else if (sample == Sample::wjets)
    dirSample = outputFile->mkdir("wjets");
  else if (sample == Sample::qcd_mu)
    dirSample = outputFile->mkdir("qcd_mu");
  else if (sample == Sample::qcd_ele)
    dirSample = outputFile->mkdir("qcd_ele");

  dirSample->cd();

  gROOT->SetBatch(kTRUE);
  TH1::SetDefaultSumw2(); //all the following histograms will automatically call TH1::Sumw2()  

  cout << endl;

  TChain* chain = new TChain("tree");
  TChain* chFriend = new TChain("mjvars/t");

  buildChain(chain, chFriend, inputDIR, sample); 

  TTreeReader reader (chain);

  // reco met
  TTreeReaderValue<Float_t> tkmet    (reader,"met_trkPt");
  TTreeReaderValue<Float_t> tkmet_phi(reader,"met_trkPhi");
  TTreeReaderValue<Float_t> pfmet    (reader,"met_pt");
  TTreeReaderValue<Float_t> pfmet_phi(reader,"met_phi");
  //TTreeReaderValue<Float_t> pfmet_mass(reader,"met_mass");

  // gen met
  // TTreeReaderValue<Float_t> genpfmet    (reader,"met_genPt");
  // TTreeReaderValue<Float_t> genpfmet_phi(reader,"met_genPhi");

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

  // mZ
  TTreeReaderValue<Float_t> mZ1 (reader,"mZ1");

  // gen particles
  // TTreeReaderValue<Int_t> nGenPart(reader,"nGenPart");
  // TTreeReaderArray<Int_t> GenPart_pdgId(reader,"GenPart_pdgId");
  // TTreeReaderArray<Int_t> GenPart_motherId(reader,"GenPart_motherId");
  // TTreeReaderArray<Int_t> GenPart_motherIndex(reader,"GenPart_motherIndex");
  // TTreeReaderArray<Float_t> GenPart_pt(reader,"GenPart_pt");
  // TTreeReaderArray<Float_t> GenPart_eta(reader,"GenPart_eta");
  // TTreeReaderArray<Float_t> GenPart_phi(reader,"GenPart_phi");
  // TTreeReaderArray<Float_t> GenPart_mass(reader,"GenPart_mass");

  // MC reweight
  TTreeReaderValue<Float_t> weight(reader,"weight");
  TTreeReaderValue<Float_t> puw(reader,"puw");

  // TH1
  TH1D* hmT = new TH1D("hmT","",nMtBins, mtMin, mtMax);
  TH1D* hmZ = new TH1D("hmZ","",nZmassBins, zMassMin, zMassMax);
  TH1D* hpfmet = new TH1D("hpfmet","",60,0,300);
  TH1D* htkmet = new TH1D("htkmet","",60,0,300);
  TH1D* hlep1pt = new TH1D("hlep1pt","",28,24,80);  //60, 0, 300
  TH1D* hlep2pt = new TH1D("hlep2pt","",35,10,80);
  TH1D* hbosonpt = new TH1D("hbosonpt","",30,0,30);
  TH1D* hbosonpt_wlike = new TH1D("hbosonpt_wlike","",30,0,30);
  TH1D* hbosoneta = new TH1D("hbosoneta","",100,-5,5);
  TH1D* hlep1sigIetaIeta = new TH1D("hlep1sigIetaIeta","",25, 0.0, 0.025);
  TH1D* hlep1relIso03 = new TH1D("hlep1relIso03","",nBins_lepIso03,min_lepIso03,max_lepIso03);
  TH1D* hlep1relIso04 = new TH1D("hlep1relIso04","",nBins_lepIso04,min_lepIso04,max_lepIso04);
  TH1D* hlep1r9 = new TH1D("hlep1r9","", 55, 0.0, 1.1);
  TH1D* hdilepInvMass = new TH1D("hdilepInvMass","", 71, 49.5, 120.5);  // same as mZ1 for Z region, and approximately W mass in W region (charged lepton and neutrino)
  TH1D* hrecoil = new TH1D("hrecoil","",30,0,30);

  TH2D* h2_mT_mZ = new TH2D("h2_mT_mZ","",nMtBins, mtMin, mtMax, nZmassBins, zMassMin, zMassMax);
  TH2D* h2_mT_lep1pt = new TH2D("h2_mT_lep1pt","",nMtBins, mtMin, mtMax, 28,24,80);
  TH2D* h2_mT_lep1eta = new TH2D("h2_mT_lep1eta","",nMtBins, mtMin, mtMax, 30, -1.5, 1.5);
  TH2D* h2_mT_lep1sigIetaIeta = new TH2D("h2_mT_lep1sigIetaIeta","",nMtBins, mtMin, mtMax, 25, 0.0, 0.025);
  TH2D* h2_mT_lep1r9 = new TH2D("h2_mT_lep1r9","",nMtBins, mtMin, mtMax, 55, 0.0, 1.1);
  TH2D* h2_mT_pfmet = new TH2D("h2_mT_pfmet","",nMtBins, mtMin, mtMax, 60,0,300);
  TH2D* h2_mT_tkmet = new TH2D("h2_mT_tkmet","",nMtBins, mtMin, mtMax, 60,0,300);
  TH2D* h2_mT_bosonPt = new TH2D("h2_mT_bosonPt","",nMtBins, mtMin, mtMax, 30,0,30);
  TH2D* h2_mT_lep1relIso03 = new TH2D("h2_mT_lep1relIso03","",nMtBins, mtMin, mtMax, nBins_lepIso03,min_lepIso03,max_lepIso03);
  TH2D* h2_mT_lep1relIso04 = new TH2D("h2_mT_lep1relIso04","",nMtBins, mtMin, mtMax, nBins_lepIso04,min_lepIso04,max_lepIso04);

  TH2D* h2_lep1pt_mZ = new TH2D("h2_lep1pt_mZ","",28,24,80, nZmassBins, zMassMin, zMassMax); 
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

  /////////////////////////////////
  // only for positive lepton
  TH1D* hmT_plus = new TH1D("hmT_plus","",nMtBins, mtMin, mtMax);
  TH1D* hmZ_plus = new TH1D("hmZ_plus","",nZmassBins, zMassMin, zMassMax);
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
  TH1D* hlep1r9_plus = new TH1D("hlep1r9_plus","", 55, 0.0, 1.1);
  TH1D* hdilepInvMass_plus = new TH1D("hdilepInvMass_plus","", 71, 49.5, 120.5);  // same as mZ1 for Z region, and approximately W mass in W region (charged lepton and neutrino)
  TH1D* hrecoil_plus = new TH1D("hrecoil_plus","",30,0,30);

  TH2D* h2_mT_mZ_plus = new TH2D("h2_mT_mZ_plus","",nMtBins, mtMin, mtMax, nZmassBins, zMassMin, zMassMax);
  TH2D* h2_mT_lep1pt_plus = new TH2D("h2_mT_lep1pt_plus","",nMtBins, mtMin, mtMax, 28,24,80);
  TH2D* h2_mT_lep1eta_plus = new TH2D("h2_mT_lep1eta_plus","",nMtBins, mtMin, mtMax, 30, -1.5, 1.5);
  TH2D* h2_mT_lep1sigIetaIeta_plus = new TH2D("h2_mT_lep1sigIetaIeta_plus","",nMtBins, mtMin, mtMax, 25, 0.0, 0.025);
  TH2D* h2_mT_lep1r9_plus = new TH2D("h2_mT_lep1r9_plus","",nMtBins, mtMin, mtMax, 55, 0.0, 1.1);
  TH2D* h2_mT_pfmet_plus = new TH2D("h2_mT_pfmet_plus","",nMtBins, mtMin, mtMax, 60,0,300);
  TH2D* h2_mT_tkmet_plus = new TH2D("h2_mT_tkmet_plus","",nMtBins, mtMin, mtMax, 60,0,300);
  TH2D* h2_mT_bosonPt_plus = new TH2D("h2_mT_bosonPt_plus","",nMtBins, mtMin, mtMax, 30,0,30);
  TH2D* h2_mT_lep1relIso03_plus = new TH2D("h2_mT_lep1relIso03_plus","",nMtBins, mtMin, mtMax, nBins_lepIso03,min_lepIso03,max_lepIso03);
  TH2D* h2_mT_lep1relIso04_plus = new TH2D("h2_mT_lep1relIso04_plus","",nMtBins, mtMin, mtMax, nBins_lepIso04,min_lepIso04,max_lepIso04);

  TH2D* h2_lep1pt_mZ_plus = new TH2D("h2_lep1pt_mZ_plus","",28,24,80, nZmassBins, zMassMin, zMassMax); 
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
  TH1D* hmZ_minus = new TH1D("hmZ_minus","",nZmassBins, zMassMin, zMassMax);
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
  TH1D* hlep1r9_minus = new TH1D("hlep1r9_minus","", 55, 0.0, 1.1);
  TH1D* hdilepInvMass_minus = new TH1D("hdilepInvMass_minus","", 71, 49.5, 120.5);  // same as mZ1 for Z region, and approximately W mass in W region (charged lepton and neutrino)
  TH1D* hrecoil_minus = new TH1D("hrecoil_minus","",30,0,30);

  TH2D* h2_mT_mZ_minus = new TH2D("h2_mT_mZ_minus","",nMtBins, mtMin, mtMax, nZmassBins, zMassMin, zMassMax);
  TH2D* h2_mT_lep1pt_minus = new TH2D("h2_mT_lep1pt_minus","",nMtBins, mtMin, mtMax, 28,24,80);
  TH2D* h2_mT_lep1eta_minus = new TH2D("h2_mT_lep1eta_minus","",nMtBins, mtMin, mtMax, 30, -1.5, 1.5);
  TH2D* h2_mT_lep1sigIetaIeta_minus = new TH2D("h2_mT_lep1sigIetaIeta_minus","",nMtBins, mtMin, mtMax, 25, 0.0, 0.025);
  TH2D* h2_mT_lep1r9_minus = new TH2D("h2_mT_lep1r9_minus","",nMtBins, mtMin, mtMax, 55, 0.0, 1.1);
  TH2D* h2_mT_pfmet_minus = new TH2D("h2_mT_pfmet_minus","",nMtBins, mtMin, mtMax, 60,0,300);
  TH2D* h2_mT_tkmet_minus = new TH2D("h2_mT_tkmet_minus","",nMtBins, mtMin, mtMax, 60,0,300);
  TH2D* h2_mT_bosonPt_minus = new TH2D("h2_mT_bosonPt_minus","",nMtBins, mtMin, mtMax, 30,0,30);
  TH2D* h2_mT_lep1relIso03_minus = new TH2D("h2_mT_lep1relIso03_minus","",nMtBins, mtMin, mtMax, nBins_lepIso03,min_lepIso03,max_lepIso03);
  TH2D* h2_mT_lep1relIso04_minus = new TH2D("h2_mT_lep1relIso04_minus","",nMtBins, mtMin, mtMax, nBins_lepIso04,min_lepIso04,max_lepIso04);

  TH2D* h2_lep1pt_mZ_minus = new TH2D("h2_lep1pt_mZ_minus","",28,24,80, nZmassBins, zMassMin, zMassMax); 
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
  
  // for the Z region to separate positive and negative W-like Z
  Bool_t negativeLeptonHasPassedSelection = false;
  Bool_t positiveLeptonHasPassedSelection = false;
  Int_t negLepIndex = -1;
  Int_t posLepIndex = -1;


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

    if (isWregion) {

      candidateLep_index = 0;
      // selection
      if (*nlep != 1) continue;    // 1 leptons                                                  
      if (fabs(lep_pdgId[candidateLep_index]) != chargedLeptonFlavour ) continue;  // electrons                                                     
      // if (isPositiveLepton) {
      // 	if (lep_pdgId[candidateLep_index] > 0) continue;   // remember, negative charged leptons have positive pdgID
      // } else {
      // 	if (lep_pdgId[candidateLep_index] < 0) continue;   // remember, negative charged leptons have positive pdgID
      // }
      if ( lep_pt[candidateLep_index] < 30.0) continue;
      if ( fabs(lep_eta[candidateLep_index]) > lepEtaThreshold ) continue; 
      if ( lep_tightId[candidateLep_index] < lepTightIdThreshold ) continue;   // tight ID
      if (not isMuon) {
	if ( lep_hltSafeID[candidateLep_index] != 1 ) continue;  //HLT safe ID for electrons
      }
      //if ( *pfmet < 30) continue;

      if ( lep_pdgId[candidateLep_index] > 0 ) negativeLeptonHasPassedSelection = true;
      if ( lep_pdgId[candidateLep_index] < 0 ) positiveLeptonHasPassedSelection = true;  	


    } else {

      // selection
      // here put only everything that is symmetric between first and second lepton
      if (*nlep != 2) continue;    // 2 leptons                                                  
      if (fabs(lep_pdgId[0]) != chargedLeptonFlavour || fabs(lep_pdgId[1]) != chargedLeptonFlavour) continue;  // electrons
      if ( (lep_pdgId[0] + lep_pdgId[1]) != 0) continue;   // opposite sign                      
      if ( *mZ1 < zMassMin || *mZ1 > zMassMax) continue;
      if ( lep_pt[0] < 30.0 || lep_pt[0] < 30.0 ) continue;
      if ( fabs(lep_eta[0]) > lepEtaThreshold || fabs(lep_eta[1]) > lepEtaThreshold ) continue;
      if (not isMuon) {
	if ( lep_hltSafeID[0] != 1 || lep_hltSafeID[1] != 1 ) continue;  //HLT safe ID for electrons
      }

      if ( lep_pdgId[0] > 0 ) {
	negLepIndex = 0;
	posLepIndex = 1;
      }
      else {
	posLepIndex = 0;
	negLepIndex = 1;
      }

      if (lep_tightId[negLepIndex] >=  lepTightIdThreshold) negativeLeptonHasPassedSelection = true;
      if (lep_tightId[posLepIndex] >=  lepTightIdThreshold) positiveLeptonHasPassedSelection = true;
      
      // electron or muon specific selections                                                          
      // if (isPositiveLepton) {
      // 	if (lep_pdgId[0] < 0) {
      // 	  candidateLep_index = 0;
      // 	  otherLep_index = 1;
      // 	} else {
      // 	  candidateLep_index = 1;
      // 	  otherLep_index = 0;
      // 	} 
      // }
      // if ( lep_pt[candidateLep_index] < 30.0 || lep_pt[otherLep_index] < 10.0) continue;

      // // electron or muon specific selections                                                          
      // if (isMuon) {
      //   if ( lep_tightId[candidateLep_index] < 1 ) continue;   // tight ID
      // } else {
      // 	if (lep_tightId[candidateLep_index] < 3 ) continue;   // tight ID
      // }

    }
    

    if (sample == Sample::data_doubleEG) wgt = 1.0;
    else wgt = intLumi * *puw * *weight;

    // convention:
    // lep1 is the one with candidateLep_index
    // lep2 is the other one, if any (this means only for Z region). For the W region, lep2 is given by MET

    TLorentzVector lep1Reco, lep2Reco, bosonReco, metReco, metWlikeReco, bosonRecoWlike;
    Double_t mT = -1.0;
    Double_t dilepInvMass -1.0;
    TVector2 recoilReco;

    if (isWregion) {

      metReco.SetPtEtaPhiM(metToUse_pt,0,metToUse_phi,0);
      metWlikeReco = metReco;
      lep1Reco.SetPtEtaPhiM(lep_pt[candidateLep_index],lep_eta[candidateLep_index],lep_phi[candidateLep_index],lep_mass[candidateLep_index]); 
      lep2Reco = metWlikeReco;

      dilepInvMass = (lep1Reco + lep2Reco).M();

      bosonReco = lep1Reco + lep2Reco;
      bosonRecoWlike = lep1Reco + metWlikeReco;
      if (bosonRecoWlike.Pt() > 30) continue;

      mT = sqrt(2. * lep1Reco.Pt() * metWlikeReco.Pt() * (1. - cos(lep1Reco.DeltaPhi(metWlikeReco)) ) );
      if (mT < mtMin) continue;
      //if (mT > mtMax) continue;    

      recoilReco = metWlikeReco.Vect().XYvector() + lep1Reco.Vect().XYvector();;
      // if (isWregion) {
      //   //recoilReco += lep1Reco.Vect().XYvector();
      //   if (recoilReco.Mod() > 30) continue;
      // } else {
      //   //recoilReco += bosonReco.Vect().XYvector();
      //   if (recoilReco.Mod() > 15) continue;
      // }


      // fill histograms for isolation before cutting on it
      fillTH1(hlep1relIso03_noIsoCut, lep_relIso03[candidateLep_index], wgt);
      fillTH2(h2_mT_lep1relIso03_noIsoCut, mT, lep_relIso03[candidateLep_index], wgt);
      fillTH2(h2_mT_lep1relIso03_noIsoCut, lep1Reco.Pt(), lep_relIso03[candidateLep_index], wgt);
      fillTH1(hlep1relIso04_noIsoCut, lep_relIso04[candidateLep_index], wgt);
      fillTH2(h2_mT_lep1relIso04_noIsoCut, mT, lep_relIso04[candidateLep_index], wgt);
      fillTH2(h2_mT_lep1relIso04_noIsoCut, lep1Reco.Pt(), lep_relIso04[candidateLep_index], wgt);


      ////////////////////////////
      // LEP ISO CUT
      /////////////////////////////
      if (isWregion) {

	if (isMuon) {

	  if (QCD_enriched_region) {
	    if ( lep_relIso04[candidateLep_index] < 0.15) continue;  // tight iso first lep                 
	  } else {
	    if ( lep_relIso04[candidateLep_index] > 0.15) continue;  // tight iso first lep                                  
	  }

	} else {

	  Double_t eleIsoThreshold = (fabs(lep_etaSc[candidateLep_index]) < 1.479) ? 0.0588 : 0.0571;
	  if (QCD_enriched_region) {	  
	    if ( lep_relIso03[candidateLep_index] < eleIsoThreshold) continue;  // tight iso first lep
	  } else {
	    if ( lep_relIso03[candidateLep_index] > eleIsoThreshold) continue;  // tight iso first lep
	  }

	}
     
      } else {

	if (isMuon) {

	  if (QCD_enriched_region) {
	    if ( lep_relIso04[candidateLep_index] < 0.15) continue;  // tight iso first lep                 
	    if ( lep_relIso04[otherLep_index] < 0.15) continue;  // tight iso second lep                 
	  } else {
	    if ( lep_relIso04[candidateLep_index] > 0.15) continue;  // tight iso first lep                                  
	    if ( lep_relIso04[otherLep_index] > 0.15) continue;  // tight iso first lep                 
	  }

	} else {

	  Double_t eleIsoThreshold = (fabs(lep_etaSc[candidateLep_index]) < 1.479) ? 0.0588 : 0.0571;
	  if (QCD_enriched_region) {	  
	    if ( lep_relIso03[candidateLep_index] < eleIsoThreshold) continue;  // tight iso first lep
	    if ( lep_relIso03[otherLep_index] < 0.175 ) continue;               // loose iso second lep
	  } else {
	    if ( lep_relIso03[candidateLep_index] > eleIsoThreshold) continue;  // tight iso first lep
	    if ( lep_relIso03[otherLep_index] > 0.175 ) continue;               // loose iso second lep
	  }

	}

      }


      //========================================================

    }

    metReco.SetPtEtaPhiM(metToUse_pt,0,metToUse_phi,0);
    metWlikeReco = metReco;

    if (negativeLeptonHasPassedSelection) {

      if (isWregion) {
	lep1Reco.SetPtEtaPhiM(lep_pt[candidateLep_index],lep_eta[candidateLep_index],lep_phi[candidateLep_index],lep_mass[candidateLep_index]); 
	lep2Reco = metWlikeReco;
      }
      else {
	
	lep2Reco.SetPtEtaPhiM(lep_pt[otherLep_index],lep_eta[otherLep_index],lep_phi[otherLep_index],lep_mass[otherLep_index]);
	metWlikeReco += lep2Reco;
      }


    }

    if (positiveLeptonHasPassedSelection) {

    } 

    if (negativeLeptonHasPassedSelection || positiveLeptonHasPassedSelection) {


    }

    dilepInvMass = (lep1Reco + lep2Reco).M();

    bosonReco = lep1Reco + lep2Reco;
    bosonRecoWlike = lep1Reco + metWlikeReco;
    if (bosonRecoWlike.Pt() > 30) continue;

    mT = sqrt(2. * lep1Reco.Pt() * metWlikeReco.Pt() * (1. - cos(lep1Reco.DeltaPhi(metWlikeReco)) ) );
    if (mT < mtMin) continue;
    //if (mT > mtMax) continue;    

    recoilReco = metWlikeReco.Vect().XYvector() + lep1Reco.Vect().XYvector();;
    // if (isWregion) {
    //   //recoilReco += lep1Reco.Vect().XYvector();
    //   if (recoilReco.Mod() > 30) continue;
    // } else {
    //   //recoilReco += bosonReco.Vect().XYvector();
    //   if (recoilReco.Mod() > 15) continue;
    // }


    // fill histograms for isolation before cutting on it
    fillTH1(hlep1relIso03_noIsoCut, lep_relIso03[candidateLep_index], wgt);
    fillTH2(h2_mT_lep1relIso03_noIsoCut, mT, lep_relIso03[candidateLep_index], wgt);
    fillTH2(h2_mT_lep1relIso03_noIsoCut, lep1Reco.Pt(), lep_relIso03[candidateLep_index], wgt);
    fillTH1(hlep1relIso04_noIsoCut, lep_relIso04[candidateLep_index], wgt);
    fillTH2(h2_mT_lep1relIso04_noIsoCut, mT, lep_relIso04[candidateLep_index], wgt);
    fillTH2(h2_mT_lep1relIso04_noIsoCut, lep1Reco.Pt(), lep_relIso04[candidateLep_index], wgt);


    ////////////////////////////
    // LEP ISO CUT
    /////////////////////////////
    if (isWregion) {

      if (isMuon) {

	if (QCD_enriched_region) {
          if ( lep_relIso04[candidateLep_index] < 0.15) continue;  // tight iso first lep                 
        } else {
          if ( lep_relIso04[candidateLep_index] > 0.15) continue;  // tight iso first lep                                  
        }

      } else {

	Double_t eleIsoThreshold = (fabs(lep_etaSc[candidateLep_index]) < 1.479) ? 0.0588 : 0.0571;
	if (QCD_enriched_region) {	  
	  if ( lep_relIso03[candidateLep_index] < eleIsoThreshold) continue;  // tight iso first lep
	} else {
	  if ( lep_relIso03[candidateLep_index] > eleIsoThreshold) continue;  // tight iso first lep
	}

      }
     
    } else {

      if (isMuon) {

	if (QCD_enriched_region) {
          if ( lep_relIso04[candidateLep_index] < 0.15) continue;  // tight iso first lep                 
          if ( lep_relIso04[otherLep_index] < 0.15) continue;  // tight iso second lep                 
        } else {
          if ( lep_relIso04[candidateLep_index] > 0.15) continue;  // tight iso first lep                                  
          if ( lep_relIso04[otherLep_index] > 0.15) continue;  // tight iso first lep                 
        }

      } else {

	Double_t eleIsoThreshold = (fabs(lep_etaSc[candidateLep_index]) < 1.479) ? 0.0588 : 0.0571;
	if (QCD_enriched_region) {	  
	  if ( lep_relIso03[candidateLep_index] < eleIsoThreshold) continue;  // tight iso first lep
	  if ( lep_relIso03[otherLep_index] < 0.175 ) continue;               // loose iso second lep
	} else {
	  if ( lep_relIso03[candidateLep_index] > eleIsoThreshold) continue;  // tight iso first lep
	  if ( lep_relIso03[otherLep_index] > 0.175 ) continue;               // loose iso second lep
	}

      }

    }
    
    fillTH1(hmT,(Double_t) mT, wgt);
    fillTH1(hmZ,(Double_t) *mZ1, wgt);
    fillTH1(hpfmet,(Double_t) *pfmet, wgt);
    fillTH1(htkmet,(Double_t) *tkmet, wgt);
    fillTH1(hlep1pt,(Double_t) lep1Reco.Pt(), wgt);
    fillTH1(hlep2pt,(Double_t) metWlikeReco.Pt(), wgt);
    fillTH1(hbosonpt,(Double_t) bosonReco.Pt(), wgt);
    fillTH1(hbosonpt_wlike,(Double_t) bosonRecoWlike.Pt(), wgt);
    fillTH1(hbosoneta,(Double_t) bosonReco.Eta(), wgt);
    fillTH1(hlep1sigIetaIeta,(Double_t) lep_sigmaIetaIeta[candidateLep_index], wgt);
    fillTH1(hlep1relIso03,(Double_t) lep_relIso03[candidateLep_index], wgt);
    fillTH1(hlep1relIso04,(Double_t) lep_relIso04[candidateLep_index], wgt);
    fillTH1(hlep1r9,(Double_t) lep_r9[candidateLep_index], wgt);
    fillTH1(hdilepInvMass, dilepInvMass, wgt);
    fillTH1(hrecoil, recoilReco.Mod(), wgt);

    fillTH2(h2_mT_mZ, mT, *mZ1, wgt);
    fillTH2(h2_mT_lep1pt, mT, lep_pt[candidateLep_index], wgt);
    fillTH2(h2_mT_lep1eta, mT, lep_eta[candidateLep_index], wgt);
    fillTH2(h2_mT_lep1sigIetaIeta, mT, lep_sigmaIetaIeta[candidateLep_index], wgt);
    fillTH2(h2_mT_lep1r9, mT, lep_r9[candidateLep_index], wgt);
    fillTH2(h2_mT_pfmet, mT, *pfmet, wgt);
    fillTH2(h2_mT_tkmet, mT, *tkmet, wgt);
    fillTH2(h2_mT_bosonPt, mT, bosonReco.Pt(), wgt);
    fillTH2(h2_mT_lep1relIso03, mT, lep_relIso03[candidateLep_index], wgt);
    fillTH2(h2_mT_lep1relIso04, mT, lep_relIso04[candidateLep_index], wgt);

    fillTH2(h2_lep1pt_mZ, lep_pt[candidateLep_index], *mZ1, wgt);
    fillTH2(h2_lep1pt_lep1eta, lep_pt[candidateLep_index], lep_eta[candidateLep_index], wgt);
    fillTH2(h2_lep1pt_lep1sigIetaIeta, lep_pt[candidateLep_index], lep_sigmaIetaIeta[candidateLep_index], wgt);
    fillTH2(h2_lep1pt_lep1r9, lep_pt[candidateLep_index], lep_r9[candidateLep_index], wgt);
    fillTH2(h2_lep1pt_pfmet, lep_pt[candidateLep_index], *pfmet, wgt);
    fillTH2(h2_lep1pt_tkmet, lep_pt[candidateLep_index], *tkmet, wgt);
    fillTH2(h2_lep1pt_bosonPt, lep_pt[candidateLep_index], bosonReco.Pt(), wgt);
    fillTH2(h2_lep1pt_lep1relIso03, lep_pt[candidateLep_index], lep_relIso03[candidateLep_index], wgt);
    fillTH2(h2_lep1pt_lep1relIso04, lep_pt[candidateLep_index], lep_relIso04[candidateLep_index], wgt);

  }

  // if (sample == Sample::data_doubleEG) { 

  //   drawCorrelationPlot(h2_mT_mZ, "transverse mass(e^{+}e^{-}) [GeV]", "invariant mass(e^{+}e^{-}) [GeV]", "correlation_dataDoubleEG_mT_mZ","data DoubleEG",outputDIR);
  //   drawCorrelationPlot(h2_mT_lep1pt, "transverse mass(e^{+}e^{-}) [GeV]", "leading electron p_{T} [GeV]", "correlation_dataDoubleEG_mT_lep1pt","data DoubleEG",outputDIR);
  //   drawCorrelationPlot(h2_mT_lep1eta, "transverse mass(e^{+}e^{-}) [GeV]", "leading electron #eta", "correlation_dataDoubleEG_mT_lep1eta","data DoubleEG",outputDIR);
  //   drawCorrelationPlot(h2_mT_lep1sigIetaIeta, "transverse mass(e^{+}e^{-}) [GeV]", "leading electron #sigma_{i#etai#eta}", "correlation_dataDoubleEG_mT_lep1sigIetaIeta","data DoubleEG",outputDIR);
  //   drawCorrelationPlot(h2_mT_lep1r9, "transverse mass(e^{+}e^{-}) [GeV]", "leading electron R9", "correlation_dataDoubleEG_mT_lep1r9","data DoubleEG",outputDIR);
  //   drawCorrelationPlot(h2_mT_pfmet, "transverse mass(e^{+}e^{-}) [GeV]", "PF E_{T}^{miss} [GeV]", "correlation_dataDoubleEG_mT_pfmet","data DoubleEG",outputDIR);
  //   drawCorrelationPlot(h2_mT_tkmet, "transverse mass(e^{+}e^{-}) [GeV]", "Tracker E_{T}^{miss} [GeV]", "correlation_dataDoubleEG_mT_tkmet","data DoubleEG",outputDIR);
  //   drawCorrelationPlot(h2_mT_bosonPt, "transverse mass(e^{+}e^{-}) [GeV]", "boson p_{T} [GeV]", "correlation_dataDoubleEG_mT_bosonPt","data DoubleEG",outputDIR);
  //   drawCorrelationPlot(h2_mT_lep1relIso03, "transverse mass(e^{+}e^{-}) [GeV]", "leading electron isolation (relIso03)", "correlation_dataDoubleEG_mT_lep1relIso03","data DoubleEG",outputDIR);
  //   drawCorrelationPlot(h2_mT_lep1relIso04, "transverse mass(e^{+}e^{-}) [GeV]", "leading electron isolation (relIso04)", "correlation_dataDoubleEG_mT_lep1relIso04","data DoubleEG",outputDIR);

  //   // cout << "Check" << endl;

  //   drawCorrelationPlot(h2_lep1pt_mZ, "leading electron p_{T} [GeV]", "invariant mass(e^{+}e^{-}) [GeV]", "correlation_dataDoubleEG_lep1pt_mZ","data DoubleEG",outputDIR);
  //   drawCorrelationPlot(h2_lep1pt_lep1eta, "leading electron p_{T} [GeV]", "leading electron #eta", "correlation_dataDoubleEG_lep1pt_lep1eta","data DoubleEG",outputDIR);
  //   drawCorrelationPlot(h2_lep1pt_lep1sigIetaIeta, "leading electron p_{T} [GeV]", "leading electron #sigma_{i#etai#eta}", "correlation_dataDoubleEG_lep1pt_lep1sigIetaIeta","data DoubleEG",outputDIR);
  //   drawCorrelationPlot(h2_lep1pt_lep1r9, "leading electron p_{T} [GeV]", "leading electron R9", "correlation_dataDoubleEG_lep1pt_lep1r9","data DoubleEG",outputDIR);
  //   drawCorrelationPlot(h2_lep1pt_pfmet, "leading electron p_{T} [GeV]", "PF E_{T}^{miss} [GeV]", "correlation_dataDoubleEG_lep1pt_pfmet","data DoubleEG",outputDIR);
  //   drawCorrelationPlot(h2_lep1pt_tkmet, "leading electron p_{T} [GeV]", "Tracker E_{T}^{miss} [GeV]", "correlation_dataDoubleEG_lep1pt_tkmet","data DoubleEG",outputDIR);
  //   drawCorrelationPlot(h2_lep1pt_bosonPt, "leading electron p_{T} [GeV]", "boson p_{T} [GeV]", "correlation_dataDoubleEG_lep1pt_bosonPt","data DoubleEG",outputDIR);
  //   drawCorrelationPlot(h2_lep1pt_lep1relIso03, "leading electron p_{T} [GeV]", "leading electron isolation (relIso03)", "correlation_dataDoubleEG_lep1pt_lep1relIso03","data DoubleEG",outputDIR);

  // }

  cout << endl;
  cout << "Writing on output file" << endl;
  outputFile->Write();

  cout << endl;

}


//=============================================================

void makeVariableHistograms(const string& inputDIR = "./", const string& outputDIR = "./", 
			    const string& outfileName = "wmass_varhists.root", 
			    const Bool_t QCD_enriched_region = false, const Bool_t isWregion = false, 
			    const Bool_t isPositiveLepton = false, const Bool_t isMuon = false) {

  if (outputDIR!= "./") system(("mkdir -p "+outputDIR).c_str());

  gROOT->SetBatch(kTRUE);
  TH1::SetDefaultSumw2(); //all the following histograms will automatically call TH1::Sumw2()  

  cout << endl;

  TFile* outputFile = new TFile((outputDIR + outfileName).c_str(),"RECREATE");
  if (!outputFile || outputFile->IsZombie()) {
    cout << "Error: file not opened. Exit" << endl;
    exit(EXIT_FAILURE);
  }

  outputFile->cd();

  if (isWregion) fillHistograms(inputDIR, outputDIR, Sample::wjets, outputFile, QCD_enriched_region, isWregion, isPositiveLepton, isMuon);
  else {
    fillHistograms(inputDIR, outputDIR, Sample::data_doubleEG, outputFile, QCD_enriched_region, isWregion, isPositiveLepton, isMuon);
    fillHistograms(inputDIR, outputDIR, Sample::zjets, outputFile, QCD_enriched_region, isWregion, isPositiveLepton, isMuon);
  }
  fillHistograms(inputDIR, outputDIR, Sample::qcd_ele, outputFile, QCD_enriched_region, isWregion, isPositiveLepton, isMuon);

  outputFile->Close();

  
}
