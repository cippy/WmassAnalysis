#include "utility.h"

using namespace std;

// define sample
//static vector<string> sampleDir = {data_doubleEG, data_singleEG, wjets, zjets, qcd, top, diboson};
//static vector<string> sampleDir = {"data_doubleEG", "zjets"};
//static Double_t intLumi = 36.4;

void makeDataMCPlots(const string& outputDIR_tmp = "./", 
		     const string& inputFileName_tmp = "wmass_varhists.root", 
		     const Bool_t isWregion = false, const Bool_t isMuon = false, 
		     const Int_t separatePositiveAndNegative = 0) {

  // separatePositiveAndNegative = 0, 1, 2 to do plot for combined +ve and -ve leptons, +ve lepton, -ve lepton
  string subdir = "combined/";
  if (separatePositiveAndNegative == 1) subdir = "positive/";
  else if (separatePositiveAndNegative == 2) subdir = "negative/";

  string outputDIR = outputDIR_tmp + subdir;
  string inputFileName = outputDIR_tmp + inputFileName_tmp;

  if (outputDIR != "./") system(("mkdir -p " + outputDIR).c_str());
  system(("cp "+ PhpToCopy + " " + outputDIR).c_str());

  gROOT->SetBatch(kTRUE);
  TH1::SetDefaultSumw2(); //all the following histograms will automatically call TH1::Sumw2()                  
  cout << endl;

  TFile* inputFile = new TFile(inputFileName.c_str(),"READ");
  if (!inputFile || inputFile->IsZombie()) {
    cout << "Error: file not opened. Exit" << endl;
    exit(EXIT_FAILURE);
  }

  inputFile->cd();

  string plotNameID = "";
  if (isWregion) {
    if (isMuon) plotNameID = "wmn";
    else plotNameID = "wen"; 
  } else {
    if (isMuon) plotNameID = "zmm";
    else plotNameID = "zee"; 
  }
  string lepton = isMuon ? "muon" : "electron";

  vector<Int_t> rebinFactor;

  vector<string> hvarName;
  hvarName.push_back("hbosoneta");      rebinFactor.push_back(2);
  hvarName.push_back("hlep1eta");      rebinFactor.push_back(1);
  if (separatePositiveAndNegative == 0) hvarName.push_back("hlep1miniRelIso_noIsoCut");      rebinFactor.push_back(1);
  hvarName.push_back("hpfmet");      rebinFactor.push_back(2);
  hvarName.push_back("htkmet");      rebinFactor.push_back(2);
  hvarName.push_back("hlep1pt");      rebinFactor.push_back(1);
  hvarName.push_back("hlep2pt");      rebinFactor.push_back(1);
  hvarName.push_back("hbosonpt");      rebinFactor.push_back(1);
  hvarName.push_back("hbosonpt_wlike");      rebinFactor.push_back(1);
  hvarName.push_back("hmT");      rebinFactor.push_back(2);
  hvarName.push_back("hrecoil");      rebinFactor.push_back(1);
  if (isMuon) {
    hvarName.push_back("hlep1relIso04");      rebinFactor.push_back(1);   
  } else {
    hvarName.push_back("hlep1sigIetaIeta");      rebinFactor.push_back(1);
    hvarName.push_back("hlep1relIso03");      rebinFactor.push_back(2);
    hvarName.push_back("hlep1r9");      rebinFactor.push_back(1);     
  }
  if (not isWregion) { hvarName.push_back("hmZ");  rebinFactor.push_back(1); }

  for (UInt_t i = 0; i < hvarName.size(); i++) {
    if (separatePositiveAndNegative == 1) hvarName[i] = hvarName[i] + "_plus";
    if (separatePositiveAndNegative == 2) hvarName[i] = hvarName[i] + "_minus";
  }     


  vector<string> xAxisTitle;  // use "title::lowrange::uprange to pass lower and upper value for range
  xAxisTitle.push_back("boson #eta");
  xAxisTitle.push_back((lepton + " #eta").c_str() );
  if (separatePositiveAndNegative == 0) xAxisTitle.push_back((lepton + " isolation (miniRelIso)").c_str());
  xAxisTitle.push_back("PF E_{T}^{miss} [GeV]");
  xAxisTitle.push_back("tracker E_{T}^{miss} [GeV]");
  xAxisTitle.push_back((lepton + " p_{T} [GeV]").c_str());
  if (isWregion) {
    if (useTrackMet) xAxisTitle.push_back("tracker E_{T}^{miss} (neutrino p_{T}) [GeV]");
    else xAxisTitle.push_back("PF E_{T}^{miss} (neutrino p_{T}) [GeV]");
  }
  else xAxisTitle.push_back("E_{T}^{miss} W-like (E_{T,no-lep2}^{miss}) [GeV]");
  xAxisTitle.push_back("boson p_{T} [GeV]");
  xAxisTitle.push_back("boson p_{T} (W-like) [GeV]");
  if (isWregion) xAxisTitle.push_back("W transverse mass [GeV]");
  else xAxisTitle.push_back("Z transverse mass (W-like) [GeV]");
  xAxisTitle.push_back("recoil [GeV]");
  if (isMuon) { 
    xAxisTitle.push_back((lepton + " isolation (relIso04)").c_str());
  } else {
    xAxisTitle.push_back((lepton + " #sigma_{i#etai#eta}").c_str());
    xAxisTitle.push_back((lepton + " isolation (relIso03)").c_str());
    xAxisTitle.push_back((lepton + " R9").c_str());
  }
  if (not isWregion) {
    if (isMuon) xAxisTitle.push_back("invariant mass(#mu^{+}#mu^{-}) [GeV]");
    else xAxisTitle.push_back("invariant mass(e^{+}e^{-}) [GeV]");
  }

  vector<string> canvasTitle;
  canvasTitle.push_back("etaBoson");
  canvasTitle.push_back("etaLep1");
  if (separatePositiveAndNegative == 0) canvasTitle.push_back("lep1miniRelIso");
  canvasTitle.push_back("pfmet");
  canvasTitle.push_back("tkmet");
  canvasTitle.push_back("ptLep1");
  canvasTitle.push_back("ptLep2");
  canvasTitle.push_back("ptBoson");
  canvasTitle.push_back("ptBosonWlike");
  canvasTitle.push_back("mtBosonWlike");
  canvasTitle.push_back("recoil");
  if (isMuon) {
    canvasTitle.push_back("lep1relIso04");
  } else {
    canvasTitle.push_back("lep1sigmaIetaIeta");
    canvasTitle.push_back("lep1relIso03");
    canvasTitle.push_back("lep1r9");
  }
  if (not isWregion) {
    canvasTitle.push_back("invMass");
  }

  //////////////////////////
  // see if using absolute or relative iso
  // add suffix about the region in canvas title
  for (UInt_t i = 0; i < canvasTitle.size(); i++) {

    if (useAbsIso) {
      if (canvasTitle[i] == "lep1relIso03") {
	canvasTitle[i] = "lep1absIso03";
	xAxisTitle[i] = lepton + " isolation (absIso03) [GeV]";
      }
      else if (canvasTitle[i] == "lep1relIso04") {
	canvasTitle[i] = "lep1absIso04";
	xAxisTitle[i] = lepton + " isolation (absIso04) [GeV]";
      }
    }

    canvasTitle[i] = canvasTitle[i] + "_" + plotNameID;
  }

  TH1D* hdata = NULL;
  TH1D* hzjets = NULL;
  TH1D* hqcd = NULL;
  TH1D* hwjets = NULL;

  for (UInt_t i = 0; i < hvarName.size(); i++) {

    //if ( i > 0 ) break;  for tests with only first entry

    vector<TH1*> stackElementMC;  // first element is the one on top of the stack
    vector<string> stackLegendMC;

    if (isWregion) {

      hwjets = (TH1D*) getHistCloneFromFile(inputFile, hvarName[i], "wjets");
      if (isMuon) hqcd   = (TH1D*) getHistCloneFromFile(inputFile, hvarName[i], "qcd_mu");
      else hqcd   = (TH1D*) getHistCloneFromFile(inputFile, hvarName[i], "qcd_ele");    
      hzjets = (TH1D*) getHistCloneFromFile(inputFile, hvarName[i], "zjets");

      checkNotNullPtr(hwjets,"hwjets");
      checkNotNullPtr(hqcd,"hqcd");
      checkNotNullPtr(hzjets,"hzjets");

      stackElementMC.push_back(hwjets);
      stackElementMC.push_back(hqcd);
      stackElementMC.push_back(hzjets);

      stackLegendMC.push_back("W(l#nu)+jets");
      stackLegendMC.push_back("QCD");
      stackLegendMC.push_back("Z(ll)+jets");

    } else {

      hdata  = (TH1D*) getHistCloneFromFile(inputFile, hvarName[i], "data_doubleEG");
      hzjets = (TH1D*) getHistCloneFromFile(inputFile, hvarName[i], "zjets");
      if (isMuon) hqcd   = (TH1D*) getHistCloneFromFile(inputFile, hvarName[i], "qcd_mu");
      else hqcd   = (TH1D*) getHistCloneFromFile(inputFile, hvarName[i], "qcd_ele");    

      checkNotNullPtr(hdata,"hdata");
      checkNotNullPtr(hzjets,"hzjets");
      checkNotNullPtr(hqcd,"hqcd");

      stackElementMC.push_back(hzjets);
      stackElementMC.push_back(hqcd);

      stackLegendMC.push_back("Z(ll)+jets");
      stackLegendMC.push_back("QCD");

    }

    if (hdata == NULL) drawTH1dataMCstack(NULL, stackElementMC, xAxisTitle[i], "Events", canvasTitle[i], outputDIR, "pseudoData", stackLegendMC, "data/MC", intLumi, rebinFactor[i]);
    else drawTH1dataMCstack(hdata, stackElementMC, xAxisTitle[i], "Events", canvasTitle[i], outputDIR, "data", stackLegendMC, "data/MC", intLumi, rebinFactor[i]);

  }
  // here plot
  // drawTH1pair(TH1* h1, TH1* h2,
  // 	      const string& xAxisName = "", const string& yAxisName = "Events", const string& canvasName = "default",
  // 	      const string& outputDIR = "./",
  // 	      const string& legEntry1 = "data", const string& legEntry2 = "MC", const string& ratioPadYaxisName = "data/MC", const Double_t lumi = -1.0)


  inputFile->Close();

  //system(("cp /afs/cern.ch/user/m/mciprian/www/index.php " + outputDIR).c_str());

  cout << endl;

}
