#include "utility.h"

using namespace std;

// define sample
//static vector<string> sampleDir = {data_doubleEG, data_singleEG, wjets, zjets, qcd, top, diboson};
//static vector<string> sampleDir = {"data_doubleEG", "zjets"};
//static Double_t intLumi = 36.4;

//=============================================


class plotManager {

public:
  plotManager(const string & histName,
	      const string & xAxisName,
	      const string & canvasName,
	      const int    & rebinFactor
	      ) 
  {
    histName_    = histName;
    xAxisName_   = xAxisName;
    canvasName_  = canvasName;
    rebinFactor_ = rebinFactor;
  };

  ~plotManager() {};

  string getHistName()    const { return histName_;    };
  string getXaxisName()   const { return xAxisName_;   };
  string getCanvasName()  const { return canvasName_;  };
  int    getRebinFactor() const { return rebinFactor_; };

  void setHistName   (const string & histName  )  { histName_    = histName;    };
  void setXaxisName  (const string & xAxisName )  { xAxisName_   = xAxisName;   };
  void setCanvasName (const string & canvasName)  { canvasName_  = canvasName;  };  
  void setRebinFactor(const int    & rebinFactor) { rebinFactor_ = rebinFactor; };

private:

  string histName_;
  string xAxisName_;
  string canvasName_;    
  int    rebinFactor_;

};

//=============================================

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

  vector<plotManager> myPlot;
  myPlot.push_back(plotManager("hlep1eta",(lepton + " #eta").c_str(),"etaLep1",1));
  myPlot.push_back(plotManager("hpfmet","PF E_{T}^{miss} [GeV]","pfmet",2));
  myPlot.push_back(plotManager("htkmet","tracker E_{T}^{miss} [GeV]","tkmet",2));
  myPlot.push_back(plotManager("hlep1pt",(lepton + " p_{T} [GeV]").c_str(),"ptLep1",1));
  myPlot.push_back(plotManager("hlep1pt2",(lepton + " p_{T}^{2} [GeV^{2}]::850,3100").c_str(),"ptLep1square",2));
  myPlot.push_back(plotManager("hlep2pt","tracker E_{T}^{miss} (neutrino p_{T}) [GeV]","ptLep2",1));
  myPlot.push_back(plotManager("hbosonpt","W p_{T} [GeV]","ptBoson",1));
  myPlot.push_back(plotManager("hmT","W transverse mass [GeV]","mtBosonWlike",2));
  myPlot.push_back(plotManager("hmT2over4","W m_{T}^{2}/4 [GeV^{2}]","mt2over4",4));
  myPlot.push_back(plotManager("hrecoil","recoil [GeV]","recoil",1));
  myPlot.push_back(plotManager("hdxy","track #Deltaxy [cm]","dxy",1));
  myPlot.push_back(plotManager("hdphiLepMet",("#Delta#phi(" + lepton + ",E_{T}^{miss})::1.5,3.2").c_str(),"dphiLepMet",1));
  myPlot.push_back(plotManager("hlep1relIso04",(lepton + " isolation (relIso04)").c_str(),"lep1relIso04",1));
  myPlot.push_back(plotManager("hlep1relIso04_noIsoCut",(lepton + " isolation (relIso04)").c_str(),"lep1relIso04_noIsoCut",1));
  if (not isMuon && separatePositiveAndNegative == 0) {
    myPlot.push_back(plotManager("hlep1sigIetaIeta",(lepton + " #sigma_{i#etai#eta}").c_str(),"lep1sigmaIetaIeta",1));
    myPlot.push_back(plotManager("hlep1relIso03",(lepton + " isolation (relIso03)").c_str(),"lep1relIso03",1));
    myPlot.push_back(plotManager("hlep1relIso03_noIsoCut",(lepton + " isolation (relIso03)::0.0,0.5").c_str(),"lep1relIso03_noIsoCut",1));
    myPlot.push_back(plotManager("hlep1r9",(lepton + " R9").c_str(),"lep1r9",1));
    myPlot.push_back(plotManager("hdetaIn",(lepton + " #Delta#eta(track,SC)").c_str(),"lep1detaIn",1));
    myPlot.push_back(plotManager("hdphiIn",(lepton + " #Delta#phi(track,SC)").c_str(),"lep1dphiIn",1));
    // myPlot.push_back(plotManager("hdetaIn_noCut",(lepton + " #Delta#eta(track,SC)").c_str(),"lep1detaIn_noCut",1));
    // myPlot.push_back(plotManager("hdphiIn_noCut",(lepton + " #Delta#phi(track,SC)").c_str(),"lep1dphiIn_noCut",1));
    // myPlot.push_back(plotManager("hdetaIn_noCutDphiDeta",(lepton + " #Delta#eta(track,SC)").c_str(),"lep1detaIn_noCutDphiDeta",1));
    // myPlot.push_back(plotManager("hdphiIn_noCutDphiDeta",(lepton + " #Delta#phi(track,SC)").c_str(),"lep1dphiIn_noCutDphiDeta",1));
    // myPlot.push_back(plotManager("hpfmet_noCutDphiDeta","PF E_{T}^{miss} [GeV]","pfmet_noCutDphiDeta",2));
    // myPlot.push_back(plotManager("hmT_noCutDphiDeta","W transverse mass [GeV]","mtBosonWlike_noCutDphiDeta",2));
    // myPlot.push_back(plotManager("hlep1pt_noCutDphiDeta",(lepton + " p_{T} [GeV]").c_str(),"ptLep1_noCutDphiDeta",1));
    // myPlot.push_back(plotManager("hmT_noCutIsoDphiDeta","W transverse mass [GeV]","mtBosonWlike_noCutIsoDphiDeta",2));
    // myPlot.push_back(plotManager("hlep1pt_noCutIsoDphiDeta",(lepton + " p_{T} [GeV]").c_str(),"ptLep1_noCutIsoDphiDeta",1));
  }

  TH1D* hdata = NULL;
  TH1D* hzjets = NULL;
  TH1D* hqcd = NULL;
  TH1D* hwjets = NULL;
  TH1D* htop = NULL;
  TH1D* hdiboson = NULL;
  TH1D* hwenujets = NULL;
  TH1D* hwmunujets = NULL;
  TH1D* hwtaunujets = NULL;

  // tmp plot to be removed to adjust settings in CMS_lumi
  TH1D* htmp1 = new TH1D("htmp1","",1,0,1);
  TH1D* htmp2 = new TH1D("htmp2","",1,0,1);
  htmp1->Fill(0.5);
  htmp2->Fill(0.5);
  vector<TH1*> htmpVec; htmpVec.push_back(htmp2);
  drawTH1dataMCstack(htmp1, htmpVec, "variable", "Events", "tmpToBeRemoved", outputDIR);
  system(("rm " + outputDIR + "*tmpToBeRemoved*").c_str());

  for (UInt_t i = 0; i < myPlot.size(); i++) {

    //if ( i > 0 ) break;  for tests with only first entry

    if (separatePositiveAndNegative == 1) myPlot[i].setHistName(myPlot[i].getHistName() + "_plus");
    if (separatePositiveAndNegative == 2) myPlot[i].setHistName(myPlot[i].getHistName() + "_plus");
    //////////////////////////
    // add suffix about the region in canvas title
    myPlot[i].setCanvasName(myPlot[i].getCanvasName() + "_" + plotNameID);


    vector<TH1*> stackElementMC;  // first element is the one on top of the stack
    vector<string> stackLegendMC;

    if (isWregion) {

      if (isMuon) {
	hdata  = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName(), "data_singleMu");
	hwmunujets = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName(), "wmunujets");
	if (useFakeRateForMuon) hqcd   = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName(), "qcd_mu_fake");
	else hqcd   = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName(), "qcd_mu");
	hwtaunujets = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName(), "wtaunujets");
	hzjets = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName(), "zjets");
	htop   = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName(), "top");
	hdiboson = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName(), "diboson");
	hwenujets = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName(), "wenujets");

	checkNotNullPtr(hdata,"hdata");
	checkNotNullPtr(hwmunujets,"hwmunujets");
	checkNotNullPtr(hqcd,"hqcd");
	checkNotNullPtr(hwtaunujets,"hwtaunujets");
	checkNotNullPtr(hzjets,"hzjets");
	checkNotNullPtr(htop,"htop");
	checkNotNullPtr(hdiboson,"hdiboson");
	checkNotNullPtr(hwenujets,"hwenujets"); // negligible, merge with Wtaunu

	if ( not hwtaunujets->Add(hwenujets) ) {
	  cout << "#### Error in makeDataMCPlots8TeV() function: failed to merge Wenu histogram in Wtaunu. Exit ..." << endl;
	  exit(EXIT_FAILURE);
	}

	// cout << hqcd->Integral() << endl;	
	// hqcd->Scale(0.00037);
	// cout << hqcd->Integral() << endl;	

	stackElementMC.push_back(hwmunujets);
	stackElementMC.push_back(hqcd);
	stackElementMC.push_back(hwtaunujets);
	stackElementMC.push_back(hzjets);
	stackElementMC.push_back(htop);
	stackElementMC.push_back(hdiboson);
	//stackElementMC.push_back(hwenujets);	

	stackLegendMC.push_back("W(#mu#nu)+jets");
	if (useFakeRateForMuon) stackLegendMC.push_back("QCD fake rate");
	else stackLegendMC.push_back("QCD");
	stackLegendMC.push_back("W(#tau#nu)+jets");
	stackLegendMC.push_back("Z(ll)+jets");
	stackLegendMC.push_back("t#bar{t}, single top");
	stackLegendMC.push_back("WW, WZ");
	//stackLegendMC.push_back("W(e#nu)+jets");

      } else {

	hdata  = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName(), "data_singleEG");
	hwenujets = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName(), "wenujets");
	hwtaunujets = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName(), "wtaunujets");
	//hqcd   = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName(), "qcd_ele");
	if (useFakeRateForElectron) hqcd   = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName(), "qcd_ele_fake");
	hzjets = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName(), "zjets");
	htop   = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName(), "top");
	hdiboson = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName(), "diboson");
	hwmunujets = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName(), "wmunujets");

	checkNotNullPtr(hdata,"hdata");
	checkNotNullPtr(hwenujets,"hwenujets");
	checkNotNullPtr(hwtaunujets,"hwtaunujets");
	//checkNotNullPtr(hqcd,"hqcd");
	if (useFakeRateForElectron) checkNotNullPtr(hqcd,"hqcd");
	checkNotNullPtr(hzjets,"hzjets");
	checkNotNullPtr(htop,"htop");
	checkNotNullPtr(hdiboson,"hdiboson");
	checkNotNullPtr(hwmunujets,"hwjmunuets"); // negligible, merge with Wtaunu

	if ( not hwtaunujets->Add(hwmunujets) ) {
	  cout << "#### Error in makeDataMCPlots8TeV() function: failed to merge Wmunu histogram in Wtaunu. Exit ..." << endl;
	  exit(EXIT_FAILURE);
	}


	stackElementMC.push_back(hwenujets);
	stackElementMC.push_back(hwtaunujets);
	//stackElementMC.push_back(hqcd);
	if (useFakeRateForElectron) stackElementMC.push_back(hqcd);
	stackElementMC.push_back(hzjets);
	stackElementMC.push_back(htop);
	stackElementMC.push_back(hdiboson);	
	// stackElementMC.push_back(hwmunujets);

	stackLegendMC.push_back("W(e#nu)+jets");
	stackLegendMC.push_back("W(#tau#nu)+jets");
	//stackLegendMC.push_back("QCD");
	if (useFakeRateForElectron) stackLegendMC.push_back("QCD fake rate");
	stackLegendMC.push_back("Z(ll)+jets");
	stackLegendMC.push_back("t#bar{t}, single top");
	stackLegendMC.push_back("WW, WZ");
	// stackLegendMC.push_back("W(#mu#nu)+jets");

      }      


    } else {

      hzjets = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName(), "zjets");
      if (isMuon) hqcd = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName(), "qcd_mu");
      else hqcd = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName(), "qcd_ele");    

      checkNotNullPtr(hdata,"hdata");
      checkNotNullPtr(hzjets,"hzjets");
      checkNotNullPtr(hqcd,"hqcd");

      stackElementMC.push_back(hzjets);
      stackElementMC.push_back(hqcd);

      stackLegendMC.push_back("Z(ll)+jets");
      stackLegendMC.push_back("QCD");

    }


    if (hdata == NULL || noPlotData) drawTH1MCstack(stackElementMC, myPlot[i].getXaxisName(), "Events", myPlot[i].getCanvasName(), outputDIR, stackLegendMC, intLumi, myPlot[i].getRebinFactor());
    else drawTH1dataMCstack(hdata, stackElementMC, myPlot[i].getXaxisName(), "Events", myPlot[i].getCanvasName(), outputDIR, "data", stackLegendMC, "data/MC", intLumi, myPlot[i].getRebinFactor());

  }

  inputFile->Close();

  cout << endl;

}
