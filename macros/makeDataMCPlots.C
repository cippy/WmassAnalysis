#include "utility.h"

using namespace std;



//=============================================

void makeDataMCPlots(const string& outputDIR_tmp = "./", 
		     const string& inputFileName_tmp = "wmass_varhists.root", 
		     const Bool_t isWregion = true, 
		     const Bool_t isMuon = false, 
		     const Int_t plot_all0_pos1_neg2_comb3 = 0,
		     const Int_t plot_all0_EB1_EE2 = 0,
		     const Bool_t noPlotData = false
		     ) 
{
  
  // plot_all0_pos1_neg2_comb3 = 1, 2, 3 to do plots for +ve lepton, -ve lepton, +ve and -ve summed (0 will result in 3, i.e. combination of the two charges)
  // 0 is never used actually, this macro is called by a script with plot_all0_pos1_neg2_comb3 = 1, 2 or 3
  // However, passing 0 will be equivalent to 3
  string subdir = "";
  if (plot_all0_EB1_EE2 == 1) subdir = "barrel/";
  else if (plot_all0_EB1_EE2 == 2) subdir = "endcap/";

  if (plot_all0_pos1_neg2_comb3 == 1) subdir += "positive/";
  else if (plot_all0_pos1_neg2_comb3 == 2) subdir += "negative/";
  else subdir += "comb/";

  string outputDIR = outputDIR_tmp + subdir;
  string inputFileName = outputDIR_tmp + inputFileName_tmp;

  createPlotDirAndCopyPhp(outputDIR);

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
  myPlot.push_back(plotManager("hnVert","number of vertices::0,40","nVert",1));
  myPlot.push_back(plotManager("hrho","#rho::0,35","rho",1));

  myPlot.push_back(plotManager("hpfmet","PF E_{T}^{miss} [GeV]::20,250","pfmet",1));
  myPlot.push_back(plotManager("htkmet","tracker E_{T}^{miss} [GeV]::0,100","tkmet",2));
  myPlot.push_back(plotManager("hlep1pt",(lepton + " p_{T} [GeV]::25,75").c_str(),"ptLep1",2));
  myPlot.push_back(plotManager("hlep1pt2",(lepton + " p_{T}^{2} [GeV^{2}]::700,3100").c_str(),"ptLep1square",2));
  myPlot.push_back(plotManager("hlep2pt","tracker E_{T}^{miss} (neutrino p_{T}) [GeV]","ptLep2",1));
  myPlot.push_back(plotManager("hbosonpt","W p_{T} [GeV]","ptBoson",1));
  myPlot.push_back(plotManager("hmT","W transverse mass [GeV]::40,110","mtBosonWlike",2));
  myPlot.push_back(plotManager("hmT2over4","W m_{T}^{2}/4 [GeV^{2}]","mt2over4",4));
  myPlot.push_back(plotManager("hrecoil","recoil [GeV]","recoil",1));
  myPlot.push_back(plotManager("hdxy","track #Deltaxy [cm]::0,0.05","dxy",1));
  myPlot.push_back(plotManager("hdphiLepMet",("#Delta#phi(" + lepton + ",E_{T}^{miss})::1.5,3.2").c_str(),"dphiLepMet",1));
  myPlot.push_back(plotManager("hnJet","number of jets (not cleaned)","njet",1));
  myPlot.push_back(plotManager("hlep1relIso04",(lepton + " isolation (relIso04)::0,0.1").c_str(),"lep1relIso04",1));
  myPlot.push_back(plotManager("hlep1relIso04_noCut",(lepton + " isolation (relIso04)").c_str(),"lep1relIso04_noCut",1));
  if (not isMuon && plot_all0_pos1_neg2_comb3 == 0) {
    // myPlot.push_back(plotManager("hlep1sigIetaIeta",(lepton + " #sigma_{i#etai#eta}").c_str(),"lep1sigmaIetaIeta",1));
    // myPlot.push_back(plotManager("hlep1relIso03",(lepton + " isolation (relIso03)").c_str(),"lep1relIso03",1));
    // myPlot.push_back(plotManager("hlep1relIso03_noIsoCut",(lepton + " isolation (relIso03)::0.0,0.5").c_str(),"lep1relIso03_noIsoCut",1));
    // myPlot.push_back(plotManager("hlep1r9",(lepton + " R9").c_str(),"lep1r9",1));
    // myPlot.push_back(plotManager("hdetaIn",(lepton + " #Delta#eta(track,SC)").c_str(),"lep1detaIn",1));
    // myPlot.push_back(plotManager("hdphiIn",(lepton + " #Delta#phi(track,SC)").c_str(),"lep1dphiIn",1));
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

  // used to get histograms of one charge and some to the othe to make the combination
  TH1D* hdata_otherCharge = NULL;
  TH1D* hzjets_otherCharge = NULL;
  TH1D* hqcd_otherCharge = NULL;
  TH1D* hwjets_otherCharge = NULL;
  TH1D* htop_otherCharge = NULL;
  TH1D* hdiboson_otherCharge = NULL;
  TH1D* hwenujets_otherCharge = NULL;
  TH1D* hwmunujets_otherCharge = NULL;
  TH1D* hwtaunujets_otherCharge = NULL;


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

    // FIXME: find a better way to exclude some histograms
    if (myPlot[i].getHistName() != "hlep1eta" && myPlot[i].getHistName() != "hnVert" && myPlot[i].getHistName() != "hrho") {
      if (plot_all0_EB1_EE2 == 1) myPlot[i].setHistName(myPlot[i].getHistName() + "_EB");
      else if (plot_all0_EB1_EE2 == 2) myPlot[i].setHistName(myPlot[i].getHistName() + "_EE");  
      if (plot_all0_EB1_EE2 == 1) myPlot[i].setCanvasName(myPlot[i].getCanvasName() + "_EB");
      else if (plot_all0_EB1_EE2 == 2) myPlot[i].setCanvasName(myPlot[i].getCanvasName() + "_EE");
    }
    if (plot_all0_pos1_neg2_comb3 == 1) myPlot[i].setHistName(myPlot[i].getHistName() + "_pos");
    else if (plot_all0_pos1_neg2_comb3 == 2) myPlot[i].setHistName(myPlot[i].getHistName() + "_neg");
    if (plot_all0_pos1_neg2_comb3 == 1) myPlot[i].setCanvasName(myPlot[i].getCanvasName() + "_pos");
    else if (plot_all0_pos1_neg2_comb3 == 2) myPlot[i].setCanvasName(myPlot[i].getCanvasName() + "_neg");
    //////////////////////////
    // add suffix about the region in canvas title
    myPlot[i].setCanvasName(myPlot[i].getCanvasName() + "_" + plotNameID);


    vector<TH1*> stackElementMC;  // first element is the one on top of the stack
    vector<string> stackLegendMC;

    if (isWregion) {

      if (isMuon) {

	if (plot_all0_pos1_neg2_comb3 > 0 && plot_all0_pos1_neg2_comb3 < 3) {

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

	} else {
	    
	    hdata  = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName()+"_pos", "data_singleMu");
	    hwmunujets = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName()+"_pos", "wmunujets");
	    if (useFakeRateForMuon) hqcd   = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName()+"_pos", "qcd_mu_fake");
	    else hqcd   = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName()+"_pos", "qcd_mu");
	    hwtaunujets = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName()+"_pos", "wtaunujets");
	    hzjets = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName()+"_pos", "zjets");
	    htop   = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName()+"_pos", "top");
	    hdiboson = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName()+"_pos", "diboson");
	    hwenujets = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName()+"_pos", "wenujets");

	    checkNotNullPtr(hdata,"hdata");
	    checkNotNullPtr(hwmunujets,"hwmunujets");
	    checkNotNullPtr(hqcd,"hqcd");
	    checkNotNullPtr(hwtaunujets,"hwtaunujets");
	    checkNotNullPtr(hzjets,"hzjets");
	    checkNotNullPtr(htop,"htop");
	    checkNotNullPtr(hdiboson,"hdiboson");
	    checkNotNullPtr(hwenujets,"hwenujets"); // negligible, merge with Wtaunu

	    hdata_otherCharge = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName()+"_neg", "data_singleMu");
	    hwmunujets_otherCharge = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName()+"_neg", "wmunujets");
	    if (useFakeRateForMuon) hqcd_otherCharge = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName()+"_neg", "qcd_mu_fake");
	    else hqcd_otherCharge = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName()+"_neg", "qcd_mu");
	    hwtaunujets_otherCharge = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName()+"_neg", "wtaunujets");
	    hzjets_otherCharge = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName()+"_neg", "zjets");
	    htop_otherCharge = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName()+"_neg", "top");
	    hdiboson_otherCharge = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName()+"_neg", "diboson");
	    hwenujets_otherCharge = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName()+"_neg", "wenujets");

	    checkNotNullPtr(hdata_otherCharge,"hdata_otherCharge");
	    checkNotNullPtr(hwmunujets_otherCharge,"hwmunujets_otherCharge");
	    checkNotNullPtr(hqcd_otherCharge,"hqcd_otherCharge");
	    checkNotNullPtr(hwtaunujets_otherCharge,"hwtaunujets_otherCharge");
	    checkNotNullPtr(hzjets_otherCharge,"hzjets_otherCharge");
	    checkNotNullPtr(htop_otherCharge,"htop_otherCharge");
	    checkNotNullPtr(hdiboson_otherCharge,"hdiboson_otherCharge");
	    checkNotNullPtr(hwenujets_otherCharge,"hwenujets_otherCharge"); // negligible, merge with Wtaunu

	    hdata->Add(hdata_otherCharge);
	    hzjets->Add(hzjets_otherCharge);
	    hqcd->Add(hqcd_otherCharge);
	    hwjets->Add(hwjets_otherCharge);
	    htop->Add(htop_otherCharge);
	    hdiboson->Add(hdiboson_otherCharge);
	    hwenujets->Add(hwenujets_otherCharge);
	    hwmunujets->Add(hwmunujets_otherCharge);
	    hwtaunujets->Add(hwtaunujets_otherCharge);

	}

	if ( not hwtaunujets->Add(hwenujets) ) {
	  cout << "#### Error in makeDataMCPlots() function: failed to merge Wenu histogram in Wtaunu. Exit ..." << endl;
	  exit(EXIT_FAILURE);
	}

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
	stackLegendMC.push_back("WW, WZ, ZZ");
	//stackLegendMC.push_back("W(e#nu)+jets");

      } else {

	if (plot_all0_pos1_neg2_comb3 > 0 && plot_all0_pos1_neg2_comb3 < 3) {

	  hdata  = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName(), "data_singleEG");
	  hwenujets = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName(), "wenujets");
	  hwtaunujets = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName(), "wtaunujets");
	  if (useFakeRateForElectron) hqcd = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName(), "qcd_ele_fake");
	  else hqcd = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName(), "qcd_ele");
	  hzjets = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName(), "zjets");
	  htop   = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName(), "top");
	  hdiboson = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName(), "diboson");
	  hwmunujets = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName(), "wmunujets");

	  checkNotNullPtr(hdata,"hdata");
	  checkNotNullPtr(hwenujets,"hwenujets");
	  checkNotNullPtr(hwtaunujets,"hwtaunujets");
	  if (useFakeRateForElectron) checkNotNullPtr(hqcd,"hqcd");
	  else checkNotNullPtr(hqcd,"hqcd");
	  checkNotNullPtr(hzjets,"hzjets");
	  checkNotNullPtr(htop,"htop");
	  checkNotNullPtr(hdiboson,"hdiboson");
	  checkNotNullPtr(hwmunujets,"hwjmunuets"); // negligible, merge with Wtaunu

	} else {

	  // get positive, then negative
	  hdata  = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName()+"_pos", "data_singleEG");
	  hwenujets = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName()+"_pos", "wenujets");
	  hwtaunujets = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName()+"_pos", "wtaunujets");
	  if (useFakeRateForElectron) hqcd   = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName()+"_pos", "qcd_ele_fake");
	  else hqcd   = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName()+"_pos", "qcd_ele");
	  hzjets = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName()+"_pos", "zjets");
	  htop   = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName()+"_pos", "top");
	  hdiboson = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName()+"_pos", "diboson");
	  hwmunujets = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName()+"_pos", "wmunujets");

	  checkNotNullPtr(hdata,"hdata");
	  checkNotNullPtr(hwenujets,"hwenujets"); 
	  checkNotNullPtr(hwtaunujets,"hwtaunujets");
	  checkNotNullPtr(hqcd,"hqcd");
	  checkNotNullPtr(hzjets,"hzjets");
	  checkNotNullPtr(htop,"htop");
	  checkNotNullPtr(hdiboson,"hdiboson");
	  checkNotNullPtr(hwmunujets,"hwmunujets"); // negligible, merge with Wtaunu

	  hdata_otherCharge = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName()+"_neg", "data_singleEG");
	  hwmunujets_otherCharge = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName()+"_neg", "wmunujets");
	  if (useFakeRateForElectron) hqcd_otherCharge = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName()+"_neg", "qcd_ele_fake");
	  else hqcd_otherCharge = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName()+"_neg", "qcd_ele");
	  hwtaunujets_otherCharge = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName()+"_neg", "wtaunujets");
	  hzjets_otherCharge = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName()+"_neg", "zjets");
	  htop_otherCharge = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName()+"_neg", "top");
	  hdiboson_otherCharge = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName()+"_neg", "diboson");
	  hwenujets_otherCharge = (TH1D*) getHistCloneFromFile(inputFile, myPlot[i].getHistName()+"_neg", "wenujets");

	  checkNotNullPtr(hdata_otherCharge,"hdata_otherCharge");
	  checkNotNullPtr(hwenujets_otherCharge,"hwenujets_otherCharge"); 
	  checkNotNullPtr(hwtaunujets_otherCharge,"hwtaunujets_otherCharge");
	  checkNotNullPtr(hqcd_otherCharge,"hqcd_otherCharge");
	  checkNotNullPtr(hzjets_otherCharge,"hzjets_otherCharge");
	  checkNotNullPtr(htop_otherCharge,"htop_otherCharge");
	  checkNotNullPtr(hdiboson_otherCharge,"hdiboson_otherCharge");
	  checkNotNullPtr(hwmunujets_otherCharge,"hwmunujets_otherCharge"); // negligible, merge with Wtaunu


	  hdata->Add(hdata_otherCharge);
	  hwenujets->Add(hwenujets_otherCharge);
	  hzjets->Add(hzjets_otherCharge);
	  hqcd->Add(hqcd_otherCharge);
	  htop->Add(htop_otherCharge);
	  hdiboson->Add(hdiboson_otherCharge);
	  hwmunujets->Add(hwmunujets_otherCharge);
	  hwtaunujets->Add(hwtaunujets_otherCharge);

	}

	if ( not hwtaunujets->Add(hwmunujets) ) {
	  cout << "#### Error in makeDataMCPlots() function: failed to merge Wmunu histogram in Wtaunu. Exit ..." << endl;
	  exit(EXIT_FAILURE);
	}

	stackElementMC.push_back(hwenujets);
	stackElementMC.push_back(hwtaunujets);
	if (useFakeRateForElectron) stackElementMC.push_back(hqcd);
	else stackElementMC.push_back(hqcd);
	stackElementMC.push_back(hzjets);
	stackElementMC.push_back(htop);
	stackElementMC.push_back(hdiboson);	
	// stackElementMC.push_back(hwmunujets);

	stackLegendMC.push_back("W(e#nu)+jets");
	stackLegendMC.push_back("W(#tau#nu)+jets");
	if (useFakeRateForElectron) stackLegendMC.push_back("QCD fake rate");
	else stackLegendMC.push_back("QCD");
	stackLegendMC.push_back("Z(ll)+jets");
	stackLegendMC.push_back("t#bar{t}, single top");
	stackLegendMC.push_back("WW, WZ, ZZ");
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

    cout << "check 2" << endl;

    if (hdata == NULL || noPlotData) drawTH1MCstack(stackElementMC, myPlot[i].getXaxisName(), "Events", myPlot[i].getCanvasName(), outputDIR, stackLegendMC, intLumi, myPlot[i].getRebinFactor());
    else drawTH1dataMCstack(hdata, stackElementMC, myPlot[i].getXaxisName(), "Events", myPlot[i].getCanvasName(), outputDIR, "data", stackLegendMC, "data/MC", intLumi, myPlot[i].getRebinFactor());

  }

  inputFile->Close();

  cout << endl;

}
