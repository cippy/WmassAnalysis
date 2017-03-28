#include "utility.h"

using namespace std;

void makeQCDstudy8TeV(const string& outputDIR_tmp = "./",
		      const string& inputFileName_sigReg = "wmass_varhists.root", // give absolute path
		      const string& inputFileName_bkgReg = "wmass_varhists.root", // give absolute path
		      const Bool_t isWregion = true, 
		      const Bool_t isMuon = false 
		      //const Int_t separatePositiveAndNegative = 0
		      ) {

  string outputDIR = outputDIR_tmp;
  if (isMuon) outputDIR += "wmunu/";
  else outputDIR += "wenu/";

  if (outputDIR != "./") system(("mkdir -p " + outputDIR).c_str());
  system(("cp "+ PhpToCopy + " " + outputDIR).c_str());

  gROOT->SetBatch(kTRUE);
  TH1::SetDefaultSumw2(); //all the following histograms will automatically call TH1::Sumw2()                                        

  cout << endl;

  TH1D* hvar_sigReg = NULL;
  TH1D* hvar_bkgReg = NULL;

  TFile* inputFile_sigReg = new TFile(inputFileName_sigReg.c_str(),"READ");
  if (!inputFile_sigReg || inputFile_sigReg->IsZombie()) {
    cout << "Error: signal file not opened. Exit" << endl;
    exit(EXIT_FAILURE);
  }

  TFile* inputFile_bkgReg = new TFile(inputFileName_bkgReg.c_str(),"READ");
  if (!inputFile_bkgReg || inputFile_bkgReg->IsZombie()) {
    cout << "Error: background file not opened. Exit" << endl;
    exit(EXIT_FAILURE);
  }

  ///////////////////////////////
  //  compare QCD in signal and background region

  for (Int_t icharge = 0; icharge < 3; icharge++) {

    string hist_suffix = "";
    string xAxisName_prefix = "";
    string canvas_suffix = "_combined";
    if (icharge == 1) {
      hist_suffix = "_plus";
      xAxisName_prefix = "positive ";
      canvas_suffix = "_posLep";
    }
    else if (icharge == 2) {
      hist_suffix = "_minus";
      xAxisName_prefix = "negative ";
      canvas_suffix = "_negLep";
    }
 
    vector<string> hvarName;
    vector<Int_t> rebinFactor;
    vector<string> xAxisname;
    vector<string> canvasTitle;

    hvarName.push_back("hmT");  
    hvarName.push_back("hlep1pt");

    if (isMuon) {

      rebinFactor.push_back(3); 
      rebinFactor.push_back(1);

      xAxisname.push_back("W(#mu#nu) transverse mass [GeV]");
      xAxisname.push_back("muon p_{T} [GeV]");

      canvasTitle.push_back("qcd_mT_wmunu");
      canvasTitle.push_back("qcd_chargLepPt_wmunu");

    } else {

      rebinFactor.push_back(6); 
      rebinFactor.push_back(2);

      xAxisname.push_back("W(e#nu) transverse mass [GeV]");
      xAxisname.push_back("electron p_{T} [GeV]");

      canvasTitle.push_back("qcd_mT_wenu");
      canvasTitle.push_back("qcd_chargLepPt_wenu");

    }

    for (UInt_t i = 0; i < hvarName.size(); i++) {
      hvarName[i] = hvarName[i] + hist_suffix;
      xAxisname[i] = xAxisName_prefix + xAxisname[i];
      canvasTitle[i] = canvasTitle[i] + canvas_suffix;
    }
    
    TH1D* htmp1 = new TH1D("htmp1","",1,0,1);
    TH1D* htmp2 = new TH1D("htmp2","",1,0,1);
    htmp1->Fill(0.5);
    htmp2->Fill(0.5);
    drawTH1pair(htmp1, htmp2, "variable", "Events", "tmpToBeRemoved", outputDIR);
    system(("rm " + outputDIR + "*tmpToBeRemoved*").c_str());

    for (UInt_t i = 0; i < hvarName.size(); i++) {
    
      inputFile_sigReg->cd();
      if (isMuon) hvar_sigReg   = (TH1D*) getHistCloneFromFile(inputFile_sigReg, hvarName[i], "qcd_mu");
      else        hvar_sigReg   = (TH1D*) getHistCloneFromFile(inputFile_sigReg, hvarName[i], "qcd_ele");
      checkNotNullPtr(hvar_sigReg,"hvar_sigReg");

      inputFile_bkgReg->cd();
      if (isMuon) hvar_bkgReg   = (TH1D*) getHistCloneFromFile(inputFile_bkgReg, hvarName[i], "qcd_mu");
      else        hvar_bkgReg   = (TH1D*) getHistCloneFromFile(inputFile_bkgReg, hvarName[i], "qcd_ele");
      checkNotNullPtr(hvar_bkgReg,"hvar_bkgReg");

      drawTH1pair(hvar_sigReg, hvar_bkgReg, xAxisname[i], "a.u.", canvasTitle[i], outputDIR, "QCD signal iso", "QCD inverted iso", "sig/bkg",-1,rebinFactor[i]);

    }

  }

  
  /////////////////////////////////
  // compare W+ and W- distributions in signal or background region separately

  TH1D* hvar_plus = NULL;
  TH1D* hvar_minus = NULL;
  vector<string> hvarName;
  vector<Int_t> rebinFactor;
  vector<string> xAxisname;
  vector<string> canvasTitle;

  hvarName.push_back("hmT");  
  hvarName.push_back("hlep1pt");

  if (isMuon) {

    rebinFactor.push_back(3); 
    rebinFactor.push_back(1); 

    xAxisname.push_back("W(#mu#nu) transverse mass [GeV]");
    xAxisname.push_back("muon p_{T} [GeV]");

    canvasTitle.push_back("qcd_mT_wmunu_comparePlusMinus");
    canvasTitle.push_back("qcd_chargLepPt_wmunu_comparePlusMinus");

  } else {

    rebinFactor.push_back(6); 
    rebinFactor.push_back(2); 

    xAxisname.push_back("W(e#nu) transverse mass [GeV]");
    xAxisname.push_back("electron p_{T} [GeV]");

    canvasTitle.push_back("qcd_mT_wenu_comparePlusMinus");
    canvasTitle.push_back("qcd_chargLepPt_wenu_comparePlusMinus");

  }

  for (UInt_t i = 0; i < hvarName.size(); i++) { 

    // signal region
    inputFile_sigReg->cd();
    if (isMuon) {
      hvar_plus   = (TH1D*) getHistCloneFromFile(inputFile_sigReg, hvarName[i]+"_plus", "qcd_mu");
      hvar_minus   = (TH1D*) getHistCloneFromFile(inputFile_sigReg, hvarName[i]+"_minus", "qcd_mu");
    } else {
      hvar_plus   = (TH1D*) getHistCloneFromFile(inputFile_sigReg, hvarName[i]+"_plus", "qcd_ele");
      hvar_minus   = (TH1D*) getHistCloneFromFile(inputFile_sigReg, hvarName[i]+"_minus", "qcd_ele");
    }
    checkNotNullPtr(hvar_plus,"hvar_plus");
    checkNotNullPtr(hvar_minus,"hvar_minus");
    
    drawTH1pair(hvar_plus, hvar_minus, xAxisname[i], "a.u.", canvasTitle[i]+"_sigRegion", outputDIR, "QCD sig W+", "QCD sig W-", "W+ / W-",-1,rebinFactor[i]);    
 
    // background region
    inputFile_bkgReg->cd();
    if (isMuon) {
      hvar_plus   = (TH1D*) getHistCloneFromFile(inputFile_bkgReg, hvarName[i]+"_plus", "qcd_mu");
      hvar_minus   = (TH1D*) getHistCloneFromFile(inputFile_bkgReg, hvarName[i]+"_minus", "qcd_mu");
    } else {
      hvar_plus   = (TH1D*) getHistCloneFromFile(inputFile_bkgReg, hvarName[i]+"_plus", "qcd_ele");
      hvar_minus   = (TH1D*) getHistCloneFromFile(inputFile_bkgReg, hvarName[i]+"_minus", "qcd_ele");
    }
    checkNotNullPtr(hvar_plus,"hvar_plus");
    checkNotNullPtr(hvar_minus,"hvar_minus");

    drawTH1pair(hvar_plus, hvar_minus, xAxisname[i], "a.u.", canvasTitle[i]+"_bkgRegion", outputDIR, "QCD bkg W+", "QCD bkg W-", "W+ / W-",-1,rebinFactor[i]); 

  }

}
