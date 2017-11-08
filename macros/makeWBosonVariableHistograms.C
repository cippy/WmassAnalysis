#include "utility.h"

#define CHECK_EVERY_N 50000
#define N_MAX_ENTRIES_PER_SAMPLE -1 // for tests, use number <= 0 to use all events in each sample

#define PFMET_CUT 0.0 // global PFMET cut (both SR and CR): met > this
#define TKMET_CUT 0.0 // global TKMET cut (both SR and CR): met > this
#define QCD_INVERTED_DXY_THR 0.00 // cut dxy > this  //0.005

#define PFMET_MU_SR 20.0 // pfmet > this value in SR
#define PFMET_MU_CR 30.0 // pfmet < this value in CR

#define PFMET_ELE_SR 20.0 // pfmet > this value in SR
#define PFMET_ELE_CR 30.0 // pfmet < this value in CR

#define WPT_MAX_CUT 40.0
#define MT_MIN_CUT 40.0
#define MT_MAX_CUT 110.0

using namespace std;

static Int_t nMtBins = 140;
static Double_t mtMin = 0;
static Double_t mtMax = 140;
static Int_t nMt2over4Bins = 300;
static Double_t mt2over4Min = 0;
static Double_t mt2over4Max = 4800;
static string recoilCorrectionsPathToFile = "/afs/cern.ch/user/m/mciprian/www/wmass/analysisPlots_8TeV/metResoResp_23Aug2017/";
static Bool_t correctRecoilResponse = false; 
static Bool_t removeZtaggedEvents = false;

// update names to 13 TeV samples
static string eleFakeRateFileName = "/afs/cern.ch/work/m/mciprian/w_mass_analysis/CMSSW_5_3_22_patch1/src/CMGTools/WMass/data/fakerate/FR_data_el_mvatrg.root";
static string eleFakeRateHistoName = "FR_FullSel_el_data_comb";
static string muFakeRateFileName = "/afs/cern.ch/work/m/mciprian/w_mass_analysis/CMSSW_5_3_22_patch1/src/CMGTools/WMass/data/fakerate/FR_data_mu_mvatrg.root";
static string muFakeRateHistoName = "FR_FullSel_mu_data_comb";

static vector<Double_t> muEtaBinEdges_double = {0.0, 0.8, 1.6, 2.1};
static vector<Double_t> eleEtaBinEdges_double = {0.0, 1.0, 1.479, 2.5};


//=============================================================

void fillHistograms(const string& inputDIR = "./", const string& outputDIR = "./", 
		    const Sample& sample = Sample::zjets, 
		    TFile* outputFile = NULL, 
		    const Bool_t QCD_enriched_region = false, const Bool_t isMuon = false) {

  gROOT->SetBatch(kTRUE);
  TH1::SetDefaultSumw2(); //all the following histograms will automatically call TH1::Sumw2()  

  cout << endl;
  cout << "================================================" << endl;
  cout << endl;

  if (outputFile == NULL) {
    cout << "Error: file is NULL. please check. Exit ..." << endl;
    exit(EXIT_FAILURE);
  }

  outputFile->cd();

  //Double_t eleIso03thr_QCD = 0.1;   // iso < this value
  Double_t eleIso04thr = 0.15;  // set in loop depending on EB or EE // iso < this value
  Double_t eleIso04thr_QCD = 0.1; // = 0.05;  // iso > this value
  Double_t muIso04thr = 0.12;       // 0.15,  0.2,   0.125
  Double_t muIso04thr_QCD = 0.20;       // 0.15,  e.g., sig region for iso < 0.15 and QCD region for iso > 0.2

  // dxy, dz, relIso04, convVeto, maxMissingHits, id; 
  eleIdWorkingPoint tightID_EB(0.05,0.1,0.0588,1,1,3);                                                                                       
  eleIdWorkingPoint tightID_EE(0.1, 0.2,0.0571,1,1,3);                                                                                       
  eleIdWorkingPoint* tightID = NULL;
 
  // be consistent with threshold when assigning histogram boundaries
  // let's use a single binning, then do the plot only in the range.

  Int_t nBins_lepIso04 = 250;  // 30
  Double_t min_lepIso04 = 0.0; // 0.0
  Double_t max_lepIso04 = 0.500; // 0.150

  Int_t chargedLeptonFlavour = isMuon ? 13 : 11;
  Double_t lepEtaMaxThreshold = isMuon ? 2.4 : 2.5; // mu: 2.1 max, ele: 1.479 (EB max), 2.4 (EE max) // EB only for electron for now
  Double_t lepEtaMinThreshold = isMuon ? -0.001 : -0.001;
  Int_t lepTightIdThreshold = isMuon ? 1 : 3;  // for electrons 2 is medium cut based id

  TDirectory *dirSample = NULL;
  string sampleDir = getStringFromEnumSample(sample).c_str();
  cout << "Sample --> " << sampleDir << endl;
  if (outputFile->GetKey(sampleDir.c_str())) dirSample = outputFile->GetDirectory(sampleDir.c_str());
  else dirSample = outputFile->mkdir(sampleDir.c_str());
  dirSample->cd();

  cout << endl;

  //string treeName = "tree";

  TChain* chain = new TChain("tree");
  // INFO: the new friend trees at 13 TeV are inside a "tree_Friend_<sampleName>.root" file, and the tree's name is "Friends"
  // friend trees are still located in a directory called "friends" with respect to base trees
  TChain* friendChain = NULL;
  if (use8TeVSample) friendChain = new TChain("mjvars/t");  // leave as NULL if you don't use friend trees
  else {
    if (sampleDir.find("data") == string::npos && sampleDir.find("fake") == string::npos)
      friendChain = new TChain("Friends");      
  }
  //TChain* friendChain = NULL;  // leave as NULL if you don't use friend trees
  TChain* SfFriendChain = NULL;
  //if (sampleDir.find("data") == string::npos && sampleDir.find("fake") == string::npos) SfFriendChain = new TChain("sf/t");  // leave as NULL if you don't use friend trees

  vector<Double_t> genwgtVec;
  buildChain(chain, genwgtVec, use8TeVSample, inputDIR, sample, friendChain, SfFriendChain); 

  // change directory again, when building chain something was messed up
  dirSample->cd();
  //  cout << "check" << endl;

  TTreeReader reader (chain);

  TTreeReaderValue<Int_t> isData  (reader,"isData");

  TTreeReaderValue<Int_t> nVert  (reader,"nVert");
  TTreeReaderValue<Float_t> rho  (reader,"rho");

  // trigger
  string trigger = isMuon ? "HLT_SingleMu" : "HLT_SingleEl"; 
  TTreeReaderValue<Int_t> HLT_Single(reader,trigger.c_str());

  // reco met
  TTreeReaderValue<Float_t> tkmet    (reader,"met_trkPt");
  TTreeReaderValue<Float_t> tkmet_phi(reader,"met_trkPhi");
  TTreeReaderValue<Float_t> pfmet    (reader,"met_pt");
  TTreeReaderValue<Float_t> pfmet_phi(reader,"met_phi");

  TTreeReaderValue<Int_t> nJet  (reader,"nJet");
  TTreeReaderArray<Float_t> Jet_pt(reader,"Jet_pt");
  TTreeReaderArray<Float_t> Jet_phi(reader,"Jet_phi");

  // LepGood for muons, LepCorr for electrons
  string lepVarToUse = isMuon ? "nLepGood" : "nLepCorr";
  //TTreeReaderArray<Int_t> nlep (reader,lepVarToUse.c_str());
  TTreeReaderValue<Int_t> nlep  (reader,"nLepGood");
  lepVarToUse = isMuon ? "LepGood_pt" : "LepCorr_pt";
  TTreeReaderArray<Float_t> lep_pt (reader,"LepGood_pt");
  lepVarToUse = isMuon ? "LepGood_eta" : "LepCorr_eta";
  TTreeReaderArray<Float_t> lep_eta (reader,"LepGood_eta");

  // LepGood branch
  TTreeReaderArray<Int_t> lep_pdgId (reader,"LepGood_pdgId");
  TTreeReaderArray<Float_t> lep_phi (reader,"LepGood_phi");
  TTreeReaderArray<Float_t> lep_mass (reader,"LepGood_mass");
  TTreeReaderArray<Float_t> lep_relIso04 (reader,"LepGood_relIso04");
  TTreeReaderArray<Int_t> lep_tightId (reader,"LepGood_tightId");

  // for electronID
  TTreeReaderArray<Float_t> lep_etaSc (reader,"LepGood_etaSc");    
  TTreeReaderArray<Int_t> lep_hltId (reader,"LepGood_hltId");  
  TTreeReaderArray<Int_t> lep_convVeto (reader,"LepGood_convVeto"); // 
  TTreeReaderArray<Int_t> lep_lostHits (reader,"LepGood_lostHits"); // 
  TTreeReaderArray<Float_t> lep_dxy (reader,"LepGood_dxy");
  TTreeReaderArray<Float_t> lep_dz (reader,"LepGood_dz");

  // TTreeReaderArray<Float_t> lep_detaIn (reader,"LepGood_dEtaScTrkIn");
  // TTreeReaderArray<Float_t> lep_dphiIn (reader,"LepGood_dPhiScTrkIn");
  // TTreeReaderArray<Float_t> lep_sigmaIetaIeta (reader,"LepGood_sigmaIetaIeta");
  // TTreeReaderArray<Float_t> lep_HoE (reader,"LepGood_hcalOverEcal");
  // TTreeReaderArray<Float_t> lep_ecalEnergy (reader,"LepGood_ecalEnergy");
  // TTreeReaderArray<Float_t> lep_eSuperClusterOverP (reader,"LepGood_eSuperClusterOverP");
  // 1/E - 1/p = 1/ecalEnergy - eSuperClusterOverP/ecalEnergy == (1 - eSuperClusterOverP)/ecalEnergy;


  // electron MVA ID at 13 TeV (2016)
  //TTreeReaderArray<Float_t> lep_eleMVAId (reader,"LepGood_mvaIdSpring16GP"); // 


  // MC reweight
  // must activate it only for non data sample
  ///////////////////////////////
  // use a dummy variable in trees to activate TTreeReaderValue even for data
  TTreeReaderValue<Float_t> *xsec = NULL;
  TTreeReaderValue<Float_t> *genWeight = NULL;
  TTreeReaderArray<Float_t> *lep1_effSF = NULL;
  if (sampleDir.find("data") == string::npos && sampleDir.find("fake") == string::npos) {
    xsec = new TTreeReaderValue<Float_t>(reader,"xsec");
    genWeight = new TTreeReaderValue<Float_t>(reader,"genWeight");
    lep1_effSF = new TTreeReaderArray<Float_t>(reader,"LepGood_effSF");
  }


  // Gen particles
  Int_t WJetsGenDaughterPdgId = 0;  // avoid using the boolean above, if this integer is > 0 then we are considering a wjets sample
  if (sampleDir.find("wenujets") != string::npos) {
    WJetsGenDaughterPdgId = 11;
  } else if (sampleDir.find("wmunujets") != string::npos) {
    WJetsGenDaughterPdgId = 13;
  } else if (sampleDir.find("wtaunujets") != string::npos) {
    WJetsGenDaughterPdgId = 15;
  }
 

  TTreeReaderValue<Int_t> *nGenPart = NULL;
  TTreeReaderArray<Int_t> *GenPart_pdgId = NULL;
  TTreeReaderArray<Int_t> *GenPart_motherId = NULL;
  if (WJetsGenDaughterPdgId > 0 ) {
    nGenPart = new TTreeReaderValue<Int_t>(reader,"nGenPart"); 
    GenPart_pdgId = new TTreeReaderArray<Int_t>(reader,"GenPart_pdgId");
    GenPart_motherId = new TTreeReaderArray<Int_t>(reader,"GenPart_motherId");
  }

  // (2) refers to 2 bins of eta, EB and EE

  // TH1
  vector< vector<TH1D*> > hmT(2);
  vector< vector<TH1D*> > hmT2over4(2);
  vector< vector<TH1D*> > hpfmet(2);
  vector< vector<TH1D*> > htkmet(2);
  vector< vector<TH1D*> > hlep1pt(2);
  vector< vector<TH1D*> > hlep1pt2(2);
  vector< vector<TH1D*> > hlep2pt(2);
  vector< vector<TH1D*> > hbosonpt(2);
  vector< vector<TH1D*> > hlep1relIso04(2);
  vector< vector<TH1D*> > hlep1relIso04_noCut(2);
  vector< vector<TH1D*> > hrecoil(2);
  vector< vector<TH1D*> > hdxy(2);
  vector< vector<TH1D*> > hdphiLepMet(2);
  vector< vector<TH1D*> > hnJet(2);

  // histograms for which we do not divide in EB or EE
  vector<TH1D*> hlep1eta;
  vector<TH1D*> hnVert;
  vector<TH1D*> hrho; 

    // TH2
  vector< vector<TH2D*> > h2_mT_lep1pt(2);
  vector< vector<TH2D*> > h2_mT_pfmet(2);
  vector< vector<TH2D*> > h2_mT_tkmet(2);
  vector< vector<TH2D*> > h2_mT_bosonPt(2);
  vector< vector<TH2D*> > h2_mT_lep1relIso04(2);
  vector< vector<TH2D*> > h2_lep1pt_pfmet(2);
  vector< vector<TH2D*> > h2_lep1pt_tkmet(2);
  vector< vector<TH2D*> > h2_lep1pt_bosonPt(2);
  vector< vector<TH2D*> > h2_lep1pt_lep1relIso04(2);
  vector< vector<TH2D*> > h2_pfmet_tkmet(2);

  vector<TH2D*> h2_mT_lep1eta;
  vector<TH2D*> h2_mT_nVert;
  vector<TH2D*> h2_lep1pt_lep1eta;
  vector<TH2D*> h2_lep1pt_nVert;

  // histograms binned in eta
  vector<Double_t> etaBinEdges_double;
  vector<string> etaBinEdges_str;
  if (isMuon) {
    for (UInt_t i = 0; i < muEtaBinEdges_double.size(); i++) {
      etaBinEdges_str.push_back( getStringFromDouble(muEtaBinEdges_double[i]) );
      etaBinEdges_double.push_back(muEtaBinEdges_double[i]);
    }
  } else {
    for (UInt_t i = 0; i < eleEtaBinEdges_double.size(); i++) {
      etaBinEdges_str.push_back( getStringFromDouble(eleEtaBinEdges_double[i]) );
      etaBinEdges_double.push_back(eleEtaBinEdges_double[i]);
    }
  }

  // here I distinguish based on charge
  vector< vector<TH1D*> > hmT_lepEtaBin(2);
  vector< vector<TH1D*> > hlep1pt_lepEtaBin(2);

  vector<string> charge = {"pos", "neg"};
  vector<string> subDetId = {"EB", "EE"};

  // TH1
  for (Int_t icharge = 0; icharge < 2 ; icharge++) {

    // histograms divided in EB or EE

    for (Int_t iEB0orEE1 = 0; iEB0orEE1 < 2; iEB0orEE1++) {

      // TH1
      hmT[iEB0orEE1].push_back( new TH1D(Form("hmT_%s_%s",subDetId[iEB0orEE1].c_str(),charge[icharge].c_str()),"",nMtBins, mtMin, mtMax) );
      hmT2over4[iEB0orEE1].push_back( new TH1D(Form("hmT2over4_%s_%s",subDetId[iEB0orEE1].c_str(),charge[icharge].c_str()),"",nMt2over4Bins, mt2over4Min, mt2over4Max) );
      hpfmet[iEB0orEE1].push_back( new TH1D(Form("hpfmet_%s_%s",subDetId[iEB0orEE1].c_str(),charge[icharge].c_str()),"",120,0,300) );
      htkmet[iEB0orEE1].push_back( new TH1D(Form("htkmet_%s_%s",subDetId[iEB0orEE1].c_str(),charge[icharge].c_str()),"",120,0,300) );
      hlep1pt[iEB0orEE1].push_back( new TH1D(Form("hlep1pt_%s_%s",subDetId[iEB0orEE1].c_str(),charge[icharge].c_str()),"",112,24,80) );  //60, 0, 300
      hlep1pt2[iEB0orEE1].push_back( new TH1D(Form("hlep1pt2_%s_%s",subDetId[iEB0orEE1].c_str(),charge[icharge].c_str()),"",125,500,3000) );
      hlep2pt[iEB0orEE1].push_back( new TH1D(Form("hlep2pt_%s_%s",subDetId[iEB0orEE1].c_str(),charge[icharge].c_str()),"",40,0,80) );
      hbosonpt[iEB0orEE1].push_back( new TH1D(Form("hbosonpt_%s_%s",subDetId[iEB0orEE1].c_str(),charge[icharge].c_str()),"",40,0,40) );
      hlep1relIso04[iEB0orEE1].push_back( new TH1D(Form("hlep1relIso04_%s_%s",subDetId[iEB0orEE1].c_str(),charge[icharge].c_str()),"",nBins_lepIso04,min_lepIso04,max_lepIso04) );
      hlep1relIso04_noCut[iEB0orEE1].push_back( new TH1D(Form("hlep1relIso04_noCut_%s_%s",subDetId[iEB0orEE1].c_str(),charge[icharge].c_str()),"",nBins_lepIso04,min_lepIso04,max_lepIso04) );
      hrecoil[iEB0orEE1].push_back( new TH1D(Form("hrecoil_%s_%s",subDetId[iEB0orEE1].c_str(),charge[icharge].c_str()),"",40,0,40) );
      hdxy[iEB0orEE1].push_back( new TH1D(Form("hdxy_%s_%s",subDetId[iEB0orEE1].c_str(),charge[icharge].c_str()),"",40,0,0.1) );
      hdphiLepMet[iEB0orEE1].push_back( new TH1D(Form("hdphiLepMet_%s_%s",subDetId[iEB0orEE1].c_str(),charge[icharge].c_str()),"",32,0,3.2) );
      hnJet[iEB0orEE1].push_back( new TH1D(Form("hnJet_%s_%s",subDetId[iEB0orEE1].c_str(),charge[icharge].c_str()),"",6,-0.5,5.5) );

      // TH2
      h2_mT_lep1pt[iEB0orEE1].push_back( new TH2D(Form("h2_mT_lep1pt_%s_%s",subDetId[iEB0orEE1].c_str(),charge[icharge].c_str()),"",nMtBins, mtMin, mtMax, 112,24,80) );
      h2_mT_pfmet[iEB0orEE1].push_back( new TH2D(Form("h2_mT_pfmet_%s_%s",subDetId[iEB0orEE1].c_str(),charge[icharge].c_str()),"",nMtBins, mtMin, mtMax, 120,0,300) );
      h2_mT_tkmet[iEB0orEE1].push_back( new TH2D(Form("h2_mT_tkmet_%s_%s",subDetId[iEB0orEE1].c_str(),charge[icharge].c_str()),"",nMtBins, mtMin, mtMax, 120,0,300) );
      h2_mT_bosonPt[iEB0orEE1].push_back( new TH2D(Form("h2_mT_bosonPt_%s_%s",subDetId[iEB0orEE1].c_str(),charge[icharge].c_str()),"",nMtBins, mtMin, mtMax, 40,0,40) );
      h2_mT_lep1relIso04[iEB0orEE1].push_back( new TH2D(Form("h2_mT_lep1relIso04_%s_%s",subDetId[iEB0orEE1].c_str(),charge[icharge].c_str()),"",nMtBins, mtMin, mtMax, nBins_lepIso04,min_lepIso04,max_lepIso04) );
      h2_lep1pt_pfmet[iEB0orEE1].push_back( new TH2D(Form("h2_lep1pt_pfmet_%s_%s",subDetId[iEB0orEE1].c_str(),charge[icharge].c_str()),"",112,24,80, 120,0,300) );
      h2_lep1pt_tkmet[iEB0orEE1].push_back( new TH2D(Form("h2_lep1pt_tkmet_%s_%s",subDetId[iEB0orEE1].c_str(),charge[icharge].c_str()),"",112,24,80, 120,0,300) );
      h2_lep1pt_bosonPt[iEB0orEE1].push_back( new TH2D(Form("h2_lep1pt_bosonPt_%s_%s",subDetId[iEB0orEE1].c_str(),charge[icharge].c_str()),"",112,24,80, 40,0,40) );
      h2_lep1pt_lep1relIso04[iEB0orEE1].push_back( new TH2D(Form("h2_lep1pt_lep1relIso04_%s_%s",subDetId[iEB0orEE1].c_str(),charge[icharge].c_str()),"",112,24,80, nBins_lepIso04,min_lepIso04,max_lepIso04) );
      h2_pfmet_tkmet[iEB0orEE1].push_back( new TH2D(Form("h2_pfmet_tkmet_%s_%s",subDetId[iEB0orEE1].c_str(),charge[icharge].c_str()),"",120,0,300, 120,0,300) );

    }

    // TH1
    hlep1eta.push_back( new TH1D(Form("hlep1eta_%s",charge[icharge].c_str()),"",50,-2.5,2.5) );  //60, 0, 300
    hnVert.push_back( new TH1D(Form("hnVert_%s",charge[icharge].c_str()),"",50,0.5,50.5) );
    hrho.push_back( new TH1D(Form("hrho_%s",charge[icharge].c_str()),"",100,0,50) );

    for (UInt_t ibin = 0; ibin < etaBinEdges_str.size() -1; ibin++) {
      hmT_lepEtaBin[icharge].push_back(new TH1D(Form("hmT_etaBin%sTo%s_%s",
						     etaBinEdges_str[ibin].c_str(),
						     etaBinEdges_str[ibin+1].c_str(),
						     charge[icharge].c_str()),
						"",nMtBins, mtMin, mtMax));
      hlep1pt_lepEtaBin[icharge].push_back(new TH1D(Form("hlep1pt_etaBin%sTo%s_%s",
							 etaBinEdges_str[ibin].c_str(),
							 etaBinEdges_str[ibin+1].c_str(),
							 charge[icharge].c_str()),
						    "", 112,24,80));
    }

    //TH2
    h2_lep1pt_lep1eta.push_back( new TH2D(Form("h2_lep1pt_lep1eta_%s",charge[icharge].c_str()),"",112,24,80, 50, -2.5, 2.5) );
    h2_lep1pt_nVert.push_back( new TH2D(Form("h2_lep1pt_nVert_%s",charge[icharge].c_str()),"",112,24,80, 50,0.5,50.5) );
    h2_mT_lep1eta.push_back( new TH2D(Form("h2_mT_lep1eta_%s",charge[icharge].c_str()),"",nMtBins, mtMin, mtMax, 50, -2.5, 2.5) );
    h2_mT_nVert.push_back( new TH2D(Form("h2_mT_nVert_%s",charge[icharge].c_str()),"",nMtBins, mtMin, mtMax, 50, 0.5, 50.5) );

  }


  // histograms with mass weights
  // vector<TH1D*> hmT_mw;
  // vector<TH1D*> hlep1pt_mw;
  //  if (sampleDir.find("data") == string::npos && sampleDir.find("fake") == string::npos) {
    //////////////////////////////////////////////////////
    // must read the tree at least once to get nWMassSteps
    //////////////////////////////////////////////////////
    // TTreeReader reader_tmp (chain);
    // TTreeReaderValue<Int_t> nWMassSteps_tmp(reader_tmp, "nWMassSteps"); // number of weights, not used, value chosen by user and put in nMassWeights
    // long int nEvents_tmp = 0;
    // while(reader_tmp.Next() && nEvents_tmp == 0){
    //   nMassWeights = *nWMassSteps_tmp;
    //   nEvents_tmp++;
    // }
    // for (Int_t im = 0; im < nMassWeights; im++) {
    //   hmT_mw.push_back(new TH1D(Form("hmT_mw_%d",im),"",nMtBins, mtMin, mtMax));
    //   hlep1pt_mw.push_back(new TH1D(Form("hlep1pt_mw_%d",im),"",112,24,80));
    // }
  // }
  


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


  //////////////////////
  // use fake rate for electrons to estimate QCD
  //////////////////////
  TFile* fakeRateFile = NULL;
  TH2F* h2fakeRate = NULL;

  if (useFakeRateForElectron && (sampleDir.find("qcd_ele_fake") != string::npos)) {

    fakeRateFile = new TFile(eleFakeRateFileName.c_str(),"READ");

    if (!fakeRateFile || fakeRateFile->IsZombie()) {
      cout << "Error: file with fake rate not opened. Exit" << endl;
      exit(EXIT_FAILURE);
    }
    h2fakeRate = (TH2F*) getHist2CloneFromFile(fakeRateFile, eleFakeRateHistoName.c_str(), ""); 
    
    checkNotNullPtr(h2fakeRate,"h2fakeRate");

  }

  //////////////////////
  // use fake rate for muons to estimate QCD
  //////////////////////
  if (useFakeRateForMuon && (sampleDir.find("qcd_mu_fake") != string::npos)) {

    fakeRateFile = new TFile(muFakeRateFileName.c_str(),"READ");

    if (!fakeRateFile || fakeRateFile->IsZombie()) {
      cout << "Error: file with fake rate not opened. Exit" << endl;
      exit(EXIT_FAILURE);
    }
    h2fakeRate = (TH2F*) getHist2CloneFromFile(fakeRateFile, muFakeRateHistoName.c_str(), ""); 
    
    checkNotNullPtr(h2fakeRate,"h2fakeRate");

  }

  Bool_t passLooseAndNotTightSeleForFR = false;
  Bool_t passLooseAndNotTightSeleForFR_noIso = false;

  ////////////////////
  ////////////////////

  // start event loop                                                                                                                     
                    
  long int nTotal = chain->GetEntries();
  long int nEvents = 0;
  long int nEventsInSample = 0; // events processed for each sample

  Int_t EB0orEE1 = -1;

  Bool_t negativeLeptonHasPassedSelection = false;
  Bool_t positiveLeptonHasPassedSelection = false;

  // flags for electrons to save whether some selection is passed and allow to fill some histograms before only one of these selection is applied
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


  // cout << "CHECK 7" << endl;


  while(reader.Next()){
  
    cout.flush();
    if(nEvents % CHECK_EVERY_N == 0) cout<<"\r"<<"Analyzing events "<<double(nEvents)/nTotal*100<<" % ";
    //cout << "entry : " << nEvents << endl;
    nEvents++;
    nEventsInSample++;

    // to debug
    // if (WJetsGenDaughterPdgId == 0) break;

    if(dynamic_cast<TChain*>(reader.GetTree())->GetFile()->GetName() != currentFile and currentFile != ""){ 
      currentFile = dynamic_cast<TChain*>(reader.GetTree())->GetFile()->GetName();                   
      ifile ++;                                                                                      
      nEventsInSample = 1; // reset nEvents when sub sample is changed (useful with N_MAX_ENTRIES_PER_SAMPLE for debugging)
    } else if(dynamic_cast<TChain*>(reader.GetTree())->GetFile()->GetName() != currentFile and currentFile == ""){
      currentFile = dynamic_cast<TChain*>(reader.GetTree())->GetFile()->GetName();                         
    }      
    if (N_MAX_ENTRIES_PER_SAMPLE > 0 && nEventsInSample > N_MAX_ENTRIES_PER_SAMPLE) continue;

    //checkPoint(currentFile,200,nEvents);

    negativeLeptonHasPassedSelection = false;
    positiveLeptonHasPassedSelection = false;
    // passDxySel = false;
    passIsoSel = false;
    passLooseAndNotTightSeleForFR = false;
    passLooseAndNotTightSeleForFR_noIso = false;
  
    // selection that is common to everything
    // selection
    if (not *HLT_Single) continue; // Single* trigger (SingleMu or SingleEG)
    if (*nlep != 1) continue;    // 1 leptons                                                  
    if (fabs(lep_pdgId[0]) != chargedLeptonFlavour) continue;  // electrons                                                     
    if (lep_pt[0] < 30.0) continue;
    if (fabs(lep_eta[0]) < lepEtaMinThreshold || fabs(lep_eta[0]) > lepEtaMaxThreshold) continue; 

    if (WJetsGenDaughterPdgId > 0) {
      Bool_t GenPartNotFound = true;
      for (Int_t igen = 0; GenPartNotFound && igen < **nGenPart; igen++) {      
    	if (fabs((*GenPart_pdgId)[igen]) == WJetsGenDaughterPdgId && fabs((*GenPart_motherId)[igen]) == 24) GenPartNotFound = false;
    	//if (not GenPartNotFound) cout << "GenPart_pdgId " << (*GenPart_pdgId)[igen] << endl;
      }
      if (GenPartNotFound) continue;
    }

    // remove Z candidate when nlep > 1
    // if (removeZtaggedEvents) {
    //   Bool_t isCandidateZJetsEvent = false;
    //   if (*nlep > 1) {
    // 	for (Int_t ilep = 1; ilep < *nlep && !isCandidateZJetsEvent; ilep++) {
    // 	  if ((lep_pdgId[0] + lep_pdgId[ilep]) == 0) {
    // 	    TLorentzVector lep1, zboson;
    // 	    lep1.SetPtEtaPhiM(lep_pt[0],lep_eta[0],lep_phi[0], lep_mass[0]) ;
    // 	    zboson.SetPtEtaPhiM(lep_pt[ilep],lep_eta[ilep],lep_phi[ilep], lep_mass[ilep]);
    // 	    zboson += lep1;
    // 	    Double_t zmass = zboson.M();
    // 	    if (zmass > 80.0 && zmass < 100.0) isCandidateZJetsEvent = true;
    // 	  }
    // 	}
    //   }
    //   if (isCandidateZJetsEvent) continue;
    // }

    if (isMuon) {

      if ( lep_tightId[0] < lepTightIdThreshold ) continue;   // tight ID
      if (fabs(lep_eta[0]) < 1.2) EB0orEE1 = 0;
      else EB0orEE1 = 1;

    } else {

      if (lep_hltId[0] < 1) continue;

      if (fabs(lep_etaSc[0]) < 1.479) {
	tightID = &tightID_EB;
      } else {
	tightID = &tightID_EB;
      }	

      if (lep_tightId[0]  < tightID->cutBasedId() || 
	  lep_lostHits[0] > tightID->maxMissingHits() || 
	  fabs(lep_dz[0]) > tightID->dz() || 
	  fabs(lep_dxy[0]) > tightID->dxy() || 
	  lep_convVeto[0] != tightID->convVeto()
	  ) 
	{ 
	  passLooseAndNotTightSeleForFR_noIso = true;
	}
      if (lep_relIso04[0] < tightID->relIso04()) passIsoSel = true;
      if (passLooseAndNotTightSeleForFR_noIso || not passIsoSel) passLooseAndNotTightSeleForFR = true;
      if (not useFakeRateForElectron and passLooseAndNotTightSeleForFR_noIso) continue; // if not using fake rate, require tight selection (without isolation here)

      if (fabs(lep_eta[0]) < 1.479) EB0orEE1 = 0;
      else EB0orEE1 = 1;

    }

    Int_t chargeIndex = 0; // 0 for positive, 1 for negative
    if (lep_pdgId[0] > 0) {
      negativeLeptonHasPassedSelection = true;
      chargeIndex = 1;
    } else {
      positiveLeptonHasPassedSelection = true;
      chargeIndex = 0;
    }


    // hardcoded value for Z cross-section, because it was updated wrt our ntuple content
    if (*isData == 1) wgt = 1.0;
    else if (sampleDir.find("zjets") != string::npos) wgt = 1000.0 * intLumi * (*lep1_effSF)[0] * **genWeight * 5765.4 * genwgtVec[ifile]; 
    else  wgt = 1000.0 * intLumi * (*lep1_effSF)[0] * **genWeight * **xsec * genwgtVec[ifile]; 

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

    if (bosonReco.Mod() > WPT_MAX_CUT) continue;
    if (pfmet2D.Mod() < PFMET_CUT) continue;    
    if (tkmet2D.Mod() < TKMET_CUT) continue;    

    mT = sqrt( 2. * lep1Reco2D.Mod() * metReco.Mod() * ( 1. - cos(lep1Reco2D.DeltaPhi(metReco)) ) );
    if (mT < MT_MIN_CUT) continue;
    if (mT > MT_MAX_CUT) continue;
    //if (mT > mtMax) continue;    
    mT2over4 = mT*mT/4.0;
    lep1pt2 = lep1pt * lep1pt; 

    dphiLepMet = fabs(metReco.DeltaPhi(lep1Reco2D));
    //////////////////////////////////////////////
    // apply some selections but do not use continue. Just fill some bool to keep track of passed selections
    ////////////////////////////////////////////
    if (isMuon) {

      if (QCD_enriched_region) {

	if (PFMET_MU_CR > 0.0 && pfmet2D.Mod() > PFMET_MU_CR) continue;
	if (fabs(lep_dxy[0]) < QCD_INVERTED_DXY_THR) continue; // apply dxy cut to select QCD in muon region
	if (lep_relIso04[0] > muIso04thr_QCD) passIsoSel = true; 

      } else {

	if (useFakeRateForMuon && (sampleDir.find("qcd_mu_fake") != string::npos)) {

	  if (lep_relIso04[0] < muIso04thr) {
	    wgt = 0.0;
	  } else {
	    // FIXME
	    // avoid overflow on x axis (which has pT <= 100.0), otherwise I don't know what happens
	    // same for y axis with eta
	    Double_t ptToFindBin = (lep_pt[0] > 100.0) ? 99.0 : lep_pt[0]; 
	    Double_t etaToFindBin = (fabs(lep_eta[0]) > 2.4) ? 2.395 : fabs(lep_eta[0]); 
	    wgt *= h2fakeRate->GetBinContent( h2fakeRate->FindBin(ptToFindBin,etaToFindBin) ); // search for global bin of the TH2
	  }
	  passIsoSel = true; 
	} else {
	  if (lep_relIso04[0] < muIso04thr) passIsoSel = true; 
	  if (fabs(lep_dxy[0]) > 0.02) continue; 
	  if (pfmet2D.Mod() < PFMET_MU_SR) continue;
	}

      }
      
    } else {

      if (QCD_enriched_region) {

	if (PFMET_ELE_CR > 0.0 && pfmet2D.Mod() > PFMET_ELE_CR) continue;
	if (lep_relIso04[0] > eleIso04thr_QCD) passIsoSel = true;


      } else {

	// when using fake rate, must use events that pass loose selection, but not tight one
	if (useFakeRateForElectron && (sampleDir.find("qcd_ele_fake") != string::npos)) {

	  if (not passLooseAndNotTightSeleForFR) {
	    wgt = 0.0;
	  } else {
	    // FIXME
	    // avoid overflow on x axis (which has pT <= 100.0), otherwise I don't know what happens
	    // same for y axis with eta
	    Double_t ptToFindBin = (lep1pt > 100.0) ? 99.0 : lep1pt; 
	    Double_t etaToFindBin = (fabs(lep_eta[0]) > 2.5) ? 2.495 : fabs(lep_eta[0]); 
	    wgt *= h2fakeRate->GetBinContent( h2fakeRate->FindBin(ptToFindBin,etaToFindBin) ); // search for global bin of the TH2
	  }
	  // when using fake rate, consider the event as if it had passed the tight selection, therefore also isolation
	  passLooseAndNotTightSeleForFR = false;
	  passLooseAndNotTightSeleForFR_noIso = false;
	  passIsoSel = true;  
	}
       
	if (pfmet2D.Mod() < PFMET_ELE_SR) continue;

      }

    }
    /////////////////////////////////////

    fillTH1(hlep1relIso04_noCut[EB0orEE1][chargeIndex],(Double_t) lep_relIso04[0], wgt);

    // now we can apply the isolation cut as well
    if (not passIsoSel) continue;

    fillTH1(hmT[EB0orEE1][chargeIndex],(Double_t) mT, wgt);
    fillTH1(hmT2over4[EB0orEE1][chargeIndex],(Double_t) mT2over4, wgt);
    fillTH1(hpfmet[EB0orEE1][chargeIndex],(Double_t) pfmet2D.Mod(), wgt);
    fillTH1(htkmet[EB0orEE1][chargeIndex],(Double_t) tkmet2D.Mod(), wgt);
    fillTH1(hlep1pt[EB0orEE1][chargeIndex],(Double_t) lep1pt, wgt);
    fillTH1(hlep1pt2[EB0orEE1][chargeIndex],(Double_t) lep1pt2, wgt);
    fillTH1(hlep2pt[EB0orEE1][chargeIndex],(Double_t) metReco.Mod(), wgt);
    fillTH1(hbosonpt[EB0orEE1][chargeIndex],(Double_t) bosonReco.Mod(), wgt);
    fillTH1(hlep1relIso04[EB0orEE1][chargeIndex],(Double_t) lep_relIso04[0], wgt);
    fillTH1(hrecoil[EB0orEE1][chargeIndex], recoilReco.Mod(), wgt);
    fillTH1(hdxy[EB0orEE1][chargeIndex],fabs(lep_dxy[0]), wgt);
    fillTH1(hdphiLepMet[EB0orEE1][chargeIndex], dphiLepMet, wgt);
    fillTH1(hnJet[EB0orEE1][chargeIndex],*nJet, wgt);

    fillTH1(hlep1eta[chargeIndex],(Double_t) lep1Reco.Eta(), wgt);
    fillTH1(hnVert[chargeIndex], *nVert, wgt);
    fillTH1(hrho[chargeIndex], *rho, wgt);

    Bool_t etaBinFound = false;
    for (UInt_t etaBinEdge = 1; etaBinEdge < etaBinEdges_double.size() && not etaBinFound; etaBinEdge++) {
      if (fabs(lep1Reco.Eta()) < etaBinEdges_double[etaBinEdge]) {
	fillTH1(hmT_lepEtaBin[chargeIndex][etaBinEdge-1],(Double_t) mT, wgt);
	fillTH1(hlep1pt_lepEtaBin[chargeIndex][etaBinEdge-1],(Double_t) lep1pt, wgt);
	etaBinFound = true;
      }
    }

    // if (WJetsGenDaughterPdgId > 0) {
    //   for (Int_t im = 0; im < nMassWeights; im++) {
    // 	fillTH1(hmT_mw[im], (Double_t) mT, wgt* (*mwWeight)[im]);
    // 	fillTH1(hlep1pt_mw[im], lep1pt, wgt* (*mwWeight)[im]);
    //   }
    // }      

    fillTH2(h2_mT_lep1pt[EB0orEE1][chargeIndex], mT, lep1pt, wgt);
    fillTH2(h2_mT_pfmet[EB0orEE1][chargeIndex], mT, pfmet2D.Mod(), wgt);
    fillTH2(h2_mT_tkmet[EB0orEE1][chargeIndex], mT, tkmet2D.Mod(), wgt);
    fillTH2(h2_mT_bosonPt[EB0orEE1][chargeIndex], mT, bosonReco.Mod(), wgt);
    fillTH2(h2_mT_lep1relIso04[EB0orEE1][chargeIndex], mT, lep_relIso04[0], wgt);
    fillTH2(h2_lep1pt_pfmet[EB0orEE1][chargeIndex], lep1pt, pfmet2D.Mod(), wgt);
    fillTH2(h2_lep1pt_tkmet[EB0orEE1][chargeIndex], lep1pt, tkmet2D.Mod(), wgt);
    fillTH2(h2_lep1pt_bosonPt[EB0orEE1][chargeIndex], lep1pt, bosonReco.Mod(), wgt);
    fillTH2(h2_lep1pt_lep1relIso04[EB0orEE1][chargeIndex], lep1pt, lep_relIso04[0], wgt);
    fillTH2(h2_pfmet_tkmet[EB0orEE1][chargeIndex], pfmet2D.Mod(), tkmet2D.Mod(), wgt); 

    fillTH2(h2_lep1pt_lep1eta[chargeIndex], lep1pt, lep_eta[0], wgt);
    fillTH2(h2_lep1pt_nVert[chargeIndex], lep1pt, *nVert, wgt);
    fillTH2(h2_mT_lep1eta[chargeIndex], mT, lep_eta[0], wgt);
    fillTH2(h2_mT_nVert[chargeIndex], mT, *nVert, wgt);


  }

  cout << endl;
  cout << "Writing on output file" << endl;
  delete recoilCorrFile;
  delete fakeRateFile;

  // if the file is opened in UPDATE mode, the following should overwrite an object if its key inside the file already exists
  // this needs to be tested
  outputFile->Write(0,TObject::kOverwrite);

  cout << endl;

}


//=============================================================

void makeWBosonVariableHistograms(const string& inputDIR = "./", const string& outputDIR = "./", 
				  const string& outfileName = "wmass_varhists.root", 
				  const Bool_t QCD_enriched_region = false, 
				  const Bool_t isMuon = false) {

  createPlotDirAndCopyPhp(outputDIR);

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
    // do both fake rate and MC for QCD (former only if required, because it might not be implemented)
    if (useFakeRateForMuon && not QCD_enriched_region) fillHistograms(inputDIR, outputDIR, Sample::qcd_mu_fake, outputFile, QCD_enriched_region, isMuon);
    fillHistograms(inputDIR, outputDIR, Sample::qcd_mu, outputFile, QCD_enriched_region, isMuon);
    fillHistograms(inputDIR, outputDIR, Sample::wmunujets, outputFile, QCD_enriched_region, isMuon);
    fillHistograms(inputDIR, outputDIR, Sample::wenujets, outputFile, QCD_enriched_region, isMuon);
  } else {
    fillHistograms(inputDIR, outputDIR, Sample::data_singleEG, outputFile, QCD_enriched_region, isMuon);
    if (useFakeRateForElectron && not QCD_enriched_region) 
      fillHistograms(inputDIR, outputDIR, Sample::qcd_ele_fake, outputFile, QCD_enriched_region, isMuon); // use fake rate only in SR
    fillHistograms(inputDIR, outputDIR, Sample::qcd_ele, outputFile, QCD_enriched_region, isMuon);
    fillHistograms(inputDIR, outputDIR, Sample::wenujets, outputFile, QCD_enriched_region, isMuon);
    fillHistograms(inputDIR, outputDIR, Sample::wmunujets, outputFile, QCD_enriched_region, isMuon);
  }
  fillHistograms(inputDIR, outputDIR, Sample::wtaunujets, outputFile, QCD_enriched_region, isMuon);
  fillHistograms(inputDIR, outputDIR, Sample::zjets, outputFile, QCD_enriched_region, isMuon);
  fillHistograms(inputDIR, outputDIR, Sample::top, outputFile, QCD_enriched_region, isMuon);
  fillHistograms(inputDIR, outputDIR, Sample::diboson, outputFile, QCD_enriched_region, isMuon);
  
  outputFile->Close();
  delete outputFile;

  system(("cp makeWBosonVariableHistograms.C " + outputDIR).c_str());
  
}
