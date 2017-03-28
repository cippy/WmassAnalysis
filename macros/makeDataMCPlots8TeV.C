#include "utility.h"

using namespace std;

// define sample
//static vector<string> sampleDir = {data_doubleEG, data_singleEG, wjets, zjets, qcd, top, diboson};
//static vector<string> sampleDir = {"data_doubleEG", "zjets"};
//static Double_t intLumi = 36.4;

void makeDataMCPlots8TeV(const string& outputDIR_tmp = "./", 
			 const string& inputFileName_tmp = "wmass_varhists.root", 
			 const Bool_t isWregion = false, 
			 const Bool_t isMuon = false, 
			 const Int_t separatePositiveAndNegative = 0,
			 const Bool_t noPlotData = false
			 ) 
{
  
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
  //hvarName.push_back("hbosoneta");      rebinFactor.push_back(2);
  hvarName.push_back("hlep1eta");      rebinFactor.push_back(1);
  hvarName.push_back("hpfmet");      rebinFactor.push_back(2);
  hvarName.push_back("htkmet");      rebinFactor.push_back(2);
  hvarName.push_back("hlep1pt");      rebinFactor.push_back(1);
  hvarName.push_back("hlep2pt");      rebinFactor.push_back(1);
  hvarName.push_back("hbosonpt");      rebinFactor.push_back(1);
  hvarName.push_back("hbosonpt_wlike");      rebinFactor.push_back(1);
  hvarName.push_back("hmT");      rebinFactor.push_back(2);
  hvarName.push_back("hrecoil");      rebinFactor.push_back(1);
  hvarName.push_back("hdxy");      rebinFactor.push_back(1);
  if (isMuon) {
    hvarName.push_back("hlep1relIso04");      rebinFactor.push_back(1);   
    hvarName.push_back("hlep1relIso04_noIsoCut");      rebinFactor.push_back(1);   
  } else {
    hvarName.push_back("hlep1sigIetaIeta");      rebinFactor.push_back(1);
    hvarName.push_back("hlep1relIso03");      rebinFactor.push_back(2);
    hvarName.push_back("hlep1r9");      rebinFactor.push_back(1);     
    hvarName.push_back("hdetaIn_noCut");      rebinFactor.push_back(1);     
    hvarName.push_back("hdphiIn_noCut");      rebinFactor.push_back(1);     
    hvarName.push_back("hdxy_noCut");      rebinFactor.push_back(1);     
  }
  if (not isWregion) { hvarName.push_back("hmZ");  rebinFactor.push_back(1); }

  for (UInt_t i = 0; i < hvarName.size(); i++) {
    if (separatePositiveAndNegative == 1) hvarName[i] = hvarName[i] + "_plus";
    if (separatePositiveAndNegative == 2) hvarName[i] = hvarName[i] + "_minus";
  }     


  vector<string> xAxisTitle;  // use "title::lowrange::uprange to pass lower and upper value for range
  //xAxisTitle.push_back("boson #eta");
  xAxisTitle.push_back((lepton + " #eta").c_str() );
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
  xAxisTitle.push_back("track #Deltaxy [cm]");
  if (isMuon) { 
    xAxisTitle.push_back((lepton + " isolation (relIso04)").c_str());
    xAxisTitle.push_back((lepton + " isolation (relIso04)").c_str());
  } else {
    xAxisTitle.push_back((lepton + " #sigma_{i#etai#eta}").c_str());
    xAxisTitle.push_back((lepton + " isolation (relIso03)").c_str());
    xAxisTitle.push_back((lepton + " R9").c_str());
    xAxisTitle.push_back((lepton + " #Delta#eta(track,SC)").c_str());
    xAxisTitle.push_back((lepton + " #Delta#eta(track,SC)").c_str());
    xAxisTitle.push_back((lepton + " #Deltaxy [cm]").c_str());
  }
  if (not isWregion) {
    if (isMuon) xAxisTitle.push_back("invariant mass(#mu^{+}#mu^{-}) [GeV]");
    else xAxisTitle.push_back("invariant mass(e^{+}e^{-}) [GeV]");
  }

  vector<string> canvasTitle;
  //canvasTitle.push_back("etaBoson");
  canvasTitle.push_back("etaLep1");
  canvasTitle.push_back("pfmet");
  canvasTitle.push_back("tkmet");
  canvasTitle.push_back("ptLep1");
  canvasTitle.push_back("ptLep2");
  canvasTitle.push_back("ptBoson");
  canvasTitle.push_back("ptBosonWlike");
  canvasTitle.push_back("mtBosonWlike");
  canvasTitle.push_back("recoil");
  canvasTitle.push_back("dxy");
  if (isMuon) {
    canvasTitle.push_back("lep1relIso04");
    canvasTitle.push_back("lep1relIso04_noIsoCut");
  } else {
    canvasTitle.push_back("lep1sigmaIetaIeta");
    canvasTitle.push_back("lep1relIso03");
    canvasTitle.push_back("lep1r9");
    canvasTitle.push_back("lep1detaIn_noCut");
    canvasTitle.push_back("lep1dphiIn_noCut");
    canvasTitle.push_back("lep1dxy_noCut");
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
  TH1D* htop = NULL;
  TH1D* hdiboson = NULL;

  // tmp plot to be removed to adjust settings in CMS_lumi
  TH1D* htmp1 = new TH1D("htmp1","",1,0,1);
  TH1D* htmp2 = new TH1D("htmp2","",1,0,1);
  htmp1->Fill(0.5);
  htmp2->Fill(0.5);
  vector<TH1*> htmpVec; htmpVec.push_back(htmp2);
  drawTH1dataMCstack(htmp1, htmpVec, "variable", "Events", "tmpToBeRemoved", outputDIR);
  system(("rm " + outputDIR + "*tmpToBeRemoved*").c_str());

  for (UInt_t i = 0; i < hvarName.size(); i++) {

    //if ( i > 0 ) break;  for tests with only first entry

    vector<TH1*> stackElementMC;  // first element is the one on top of the stack
    vector<string> stackLegendMC;

    if (isWregion) {

      if (isMuon) {
	hdata  = (TH1D*) getHistCloneFromFile(inputFile, hvarName[i], "data_singleMu");
	hwjets = (TH1D*) getHistCloneFromFile(inputFile, hvarName[i], "wjets");
	hqcd   = (TH1D*) getHistCloneFromFile(inputFile, hvarName[i], "qcd_mu");
	hzjets = (TH1D*) getHistCloneFromFile(inputFile, hvarName[i], "zjets");
	htop   = (TH1D*) getHistCloneFromFile(inputFile, hvarName[i], "top");
	hdiboson = (TH1D*) getHistCloneFromFile(inputFile, hvarName[i], "diboson");

	checkNotNullPtr(hdata,"hdata");
	checkNotNullPtr(hwjets,"hwjets");
	checkNotNullPtr(hqcd,"hqcd");
	checkNotNullPtr(hzjets,"hzjets");
	checkNotNullPtr(htop,"htop");
	checkNotNullPtr(hdiboson,"hdiboson");

	// cout << hqcd->Integral() << endl;	
	// hqcd->Scale(0.00037);
	// cout << hqcd->Integral() << endl;	

	stackElementMC.push_back(hwjets);
	stackElementMC.push_back(hqcd);
	stackElementMC.push_back(hzjets);
	stackElementMC.push_back(htop);
	stackElementMC.push_back(hdiboson);
	
	stackLegendMC.push_back("W(l#nu)+jets");
	stackLegendMC.push_back("QCD");
	stackLegendMC.push_back("Z(ll)+jets");
	stackLegendMC.push_back("T#bar{T}, single top");
	stackLegendMC.push_back("WW, WZ");

      } else {
	hdata  = (TH1D*) getHistCloneFromFile(inputFile, hvarName[i], "data_singleEG");
	hwjets = (TH1D*) getHistCloneFromFile(inputFile, hvarName[i], "wjets");
	//hqcd   = (TH1D*) getHistCloneFromFile(inputFile, hvarName[i], "qcd_ele");    
	hzjets = (TH1D*) getHistCloneFromFile(inputFile, hvarName[i], "zjets");

	checkNotNullPtr(hdata,"hdata");
	checkNotNullPtr(hwjets,"hwjets");
	//checkNotNullPtr(hqcd,"hqcd");
	checkNotNullPtr(hzjets,"hzjets");
	
	stackElementMC.push_back(hwjets);
	//stackElementMC.push_back(hqcd);
	stackElementMC.push_back(hzjets);
	
	stackLegendMC.push_back("W(l#nu)+jets");
	//stackLegendMC.push_back("QCD");
	stackLegendMC.push_back("Z(ll)+jets");
      }      


    } else {

      hzjets = (TH1D*) getHistCloneFromFile(inputFile, hvarName[i], "zjets");
      if (isMuon) hqcd = (TH1D*) getHistCloneFromFile(inputFile, hvarName[i], "qcd_mu");
      else hqcd = (TH1D*) getHistCloneFromFile(inputFile, hvarName[i], "qcd_ele");    

      checkNotNullPtr(hdata,"hdata");
      checkNotNullPtr(hzjets,"hzjets");
      checkNotNullPtr(hqcd,"hqcd");

      stackElementMC.push_back(hzjets);
      stackElementMC.push_back(hqcd);

      stackLegendMC.push_back("Z(ll)+jets");
      stackLegendMC.push_back("QCD");

    }


    if (hdata == NULL || noPlotData) drawTH1MCstack(stackElementMC, xAxisTitle[i], "Events", canvasTitle[i], outputDIR, stackLegendMC, intLumi, rebinFactor[i]);
    else drawTH1dataMCstack(hdata, stackElementMC, xAxisTitle[i], "Events", canvasTitle[i], outputDIR, "data", stackLegendMC, "data/MC", intLumi, rebinFactor[i]);

  }

  inputFile->Close();

  cout << endl;

}
