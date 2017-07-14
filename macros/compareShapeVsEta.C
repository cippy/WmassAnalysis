#include "utility.h"

using namespace std;

Bool_t histBinnedInEtaAreInSameFile = true;

void doCompareQCDshape(const string& outputDIR_tmp = "./",
		       const string& canvasName_tmp = "shapeCompVsEta",
		       const Bool_t isMuon = false,
		       const string& hvarName_base = "hmT",
		       const Bool_t scaleHistograms = true,
		       const Bool_t QCD_enriched_region = false,
		       const Bool_t useDataOnlyForQCD_CR = false,
		       const Bool_t useQCDMonteCarloForMuonsInCR = false,
		       const Int_t plot_comb0_pos1_neg2 = 0
		       ) {
  
  string outputDIR = outputDIR_tmp + "shape_study/";
  if (isMuon) outputDIR += "wmunu/";
  else outputDIR += "wenu/";
  outputDIR += "shapeVsEta/";
  if (plot_comb0_pos1_neg2 == 0) outputDIR += "combined/";
  else if (plot_comb0_pos1_neg2 == 1) outputDIR += "positive/";
  else if (plot_comb0_pos1_neg2 == 2) outputDIR += "negative/";

  if (outputDIR != "./") system(("mkdir -p " + outputDIR).c_str());
  system(("cp "+ PhpToCopy + " " + outputDIR).c_str());

  string canvasName = canvasName_tmp;

  if (hvarName_base == "hmT") canvasName += "_mT";
  else if (hvarName_base == "hlep1pt") canvasName += "_lep1pt";

  if (not scaleHistograms) canvasName += "_realNorm";

  gROOT->SetBatch(kTRUE);
  TH1::SetDefaultSumw2(); //all the following histograms will automatically call TH1::Sumw2()                                        

  cout << endl;

  vector<TH1D*> hvar;
  string hvarName = hvarName_base;

  vector<string> inputFileName;

  string subDir = "";
  if (QCD_enriched_region) {
    if (isMuon) subDir = "Wmunu_QCD_CR/";
    else subDir = "Wenu_QCD_CR/";
  } else {
    if (isMuon) subDir = "Wmunu/";
    else subDir = "Wenu/";
  }

  if (not histBinnedInEtaAreInSameFile) {
    inputFileName.push_back("/afs/cern.ch/user/m/mciprian/www/wmass/analysisPlots_8TeV/varHist_treesV2_eleSkimTree_20may2017_1848_separateWjets/" + subDir + "wmass_varhists.root");
    inputFileName.push_back("/afs/cern.ch/user/m/mciprian/www/wmass/analysisPlots_8TeV/varHist_treesV2_eleSkimTree_20may2017_1848_separateWjets_erMu0p8Ele1p0/" + subDir + "wmass_varhists.root");
    inputFileName.push_back("/afs/cern.ch/user/m/mciprian/www/wmass/analysisPlots_8TeV/varHist_treesV2_eleSkimTree_erMu0p8to1p6Ele1p0to1p479/" + subDir + "wmass_varhists.root");
    inputFileName.push_back("/afs/cern.ch/user/m/mciprian/www/wmass/analysisPlots_8TeV/varHist_treesV2_eleSkimTree_erMu1p6to2p1Ele1p479to2p4/" + subDir + "wmass_varhists.root");
  } else {
    // use n equal file, jus to keep the framework the same
    inputFileName.push_back(outputDIR_tmp + subDir + "wmass_varhists.root");
    inputFileName.push_back(outputDIR_tmp + subDir + "wmass_varhists.root");
    inputFileName.push_back(outputDIR_tmp + subDir + "wmass_varhists.root");
    inputFileName.push_back(outputDIR_tmp + subDir + "wmass_varhists.root");
  }

  // used when histBinnedInEtaAreInSameFile = true
  vector<string> histNameEtaBin_suffix;
  histNameEtaBin_suffix.push_back("");  // leave first empty for unbinned histogram

  vector<string> legEntry;
  legEntry.push_back("no #eta cut");
  if (isMuon) {
    legEntry.push_back("|#eta| < 0.8");
    legEntry.push_back("0.8 < |#eta| < 1.6");
    legEntry.push_back("1.6 < |#eta| < 2.1");
    if (histBinnedInEtaAreInSameFile) {
      histNameEtaBin_suffix.push_back("_etaBin0p0To0p8");
      histNameEtaBin_suffix.push_back("_etaBin0p8To1p6");
      histNameEtaBin_suffix.push_back("_etaBin1p6To2p1");
    }
  } else {
    legEntry.push_back("|#eta| < 1.0");
    legEntry.push_back("1.0 < |#eta| < 1.479");
    legEntry.push_back("1.479 < |#eta| < 2.5");
    if (histBinnedInEtaAreInSameFile) {
      histNameEtaBin_suffix.push_back("_etaBin0p0To1p0");
      histNameEtaBin_suffix.push_back("_etaBin1p0To1p479");
      histNameEtaBin_suffix.push_back("_etaBin1p479To2p5");
    }
  }

  if (plot_comb0_pos1_neg2 == 1) {
    for (UInt_t i = 0; i < histNameEtaBin_suffix.size(); i++) {
      histNameEtaBin_suffix[i] += "_plus"; 
    }
  } else if (plot_comb0_pos1_neg2 == 2) {
    for (UInt_t i = 0; i < histNameEtaBin_suffix.size(); i++) {
      histNameEtaBin_suffix[i] += "_minus";
    }
  }


  // W in SR from MC
  if (not QCD_enriched_region) {

    for (UInt_t i = 0; i < inputFileName.size(); i++) {
    
      TFile* inputFile = new TFile(inputFileName[i].c_str(),"READ");
      if (!inputFile || inputFile->IsZombie()) {
	cout << "Error: file not opened. Exit" << endl;
	exit(EXIT_FAILURE);
      }
      if (histBinnedInEtaAreInSameFile) hvarName = hvarName_base + histNameEtaBin_suffix[i];

      if (isMuon) hvar.push_back( (TH1D*) getHistCloneFromFile(inputFile, hvarName, "wmunujets") );
      else        hvar.push_back( (TH1D*) getHistCloneFromFile(inputFile, hvarName, "wenujets") );
      checkNotNullPtr(hvar[i],"hvar[i]");
      //if (hvarName_base == "hmT") myRebinHisto(hvar.back(),1);
      if (scaleHistograms) hvar.back()->Scale(1./hvar.back()->Integral());

    }

  } else {

    if (useDataOnlyForQCD_CR) {

      for (UInt_t i = 0; i < inputFileName.size(); i++) {
      
	TFile* inputFile = new TFile(inputFileName[i].c_str(),"READ");
	if (!inputFile || inputFile->IsZombie()) {
	  cout << "Error: file not opened. Exit" << endl;
	  exit(EXIT_FAILURE);
	}

	if (histBinnedInEtaAreInSameFile) hvarName = hvarName_base + histNameEtaBin_suffix[i];
	if (isMuon) hvar.push_back( (TH1D*) getHistCloneFromFile(inputFile, hvarName, "data_singleMu") );
	else        hvar.push_back( (TH1D*) getHistCloneFromFile(inputFile, hvarName, "data_singleEG") );
	checkNotNullPtr(hvar.back(),"hvar.back()");
	if (hvarName_base == "hmT") {
	  if (isMuon) myRebinHisto(hvar.back(),2);
	  else myRebinHisto(hvar.back(),6);
	}
	if (scaleHistograms) hvar.back()->Scale(1./hvar.back()->Integral());

      }

    } else if (isMuon && useQCDMonteCarloForMuonsInCR) {

      // QCD in CR obtained from MonteCarlo
      
      for (UInt_t i = 0; i < inputFileName.size(); i++) {
	  
	TFile* inputFile = new TFile(inputFileName[i].c_str(),"READ");
	if (!inputFile || inputFile->IsZombie()) {
	  cout << "Error: file not opened. Exit" << endl;
	  exit(EXIT_FAILURE);
	}
	if (histBinnedInEtaAreInSameFile) hvarName = hvarName_base + histNameEtaBin_suffix[i];
	hvar.push_back( (TH1D*) getHistCloneFromFile(inputFile, hvarName, "qcd_mu") );
	checkNotNullPtr(hvar.back(),"hvar.back()");
	
	if (hvarName_base == "hmT") myRebinHisto(hvar.back(),2);
	if (scaleHistograms) hvar.back()->Scale(1./hvar.back()->Integral());
	
      }
      
    } else {

      // QCD in CR obtained as data-MC (MC does not include QCD). 
      // This assumes that any discrepancy between data and MC can be attributed to QCD MC, which could be false

      for (UInt_t i = 0; i < inputFileName.size(); i++) {
      
	TFile* inputFile = new TFile(inputFileName[i].c_str(),"READ");
	if (!inputFile || inputFile->IsZombie()) {
	  cout << "Error: file not opened. Exit" << endl;
	  exit(EXIT_FAILURE);
	}

	TH1D* hdata = NULL;
	TH1D* hallMCnoQCD = NULL;
	TH1D* hdataMinusAllMCnoQCD = NULL;
	TH1D* hmc = NULL;

	if (histBinnedInEtaAreInSameFile) hvarName = hvarName_base + histNameEtaBin_suffix[i];
	if (isMuon) hdata = (TH1D*) getHistCloneFromFile(inputFile, hvarName, "data_singleMu");
	else        hdata = (TH1D*) getHistCloneFromFile(inputFile, hvarName, "data_singleEG");
	checkNotNullPtr(hdata,"hdata");
	if (hvarName_base == "hmT") {
	  if (isMuon) myRebinHisto(hdata,2);
	  else myRebinHisto(hdata,6);
	}

	vector<string> listAllMC = {"zjets", "wtaunujets", "wmunujets", "wenujets", "top", "diboson"};
	hallMCnoQCD = (TH1D*) hdata->Clone("hallMCnoQCD");
	hallMCnoQCD->Reset("ICES");

	for (UInt_t imc = 0; imc < listAllMC.size(); imc++) {

	  inputFile->cd();
	  hmc = (TH1D*) getHistCloneFromFile(inputFile, hvarName, listAllMC[imc]);
	  checkNotNullPtr(hmc,"hmc");
	  if (hvarName_base == "hmT") {
	    if (isMuon) myRebinHisto(hmc,2);
	    else myRebinHisto(hmc,6);
	  }
	  hallMCnoQCD->Add(hmc);     

	}

	hdataMinusAllMCnoQCD = (TH1D*) hdata->Clone("hdataMinusAllMCnoQCD");
	hdataMinusAllMCnoQCD->Add(hallMCnoQCD,-1.0);
	if (scaleHistograms) hdataMinusAllMCnoQCD->Scale(1./hdataMinusAllMCnoQCD->Integral());
      
	hvar.push_back( (TH1D*) hdataMinusAllMCnoQCD->Clone(Form("dataMinusMC_%d",i)) );
	
      }

    }

  }

  TCanvas* cToRemove = new TCanvas("cToRemove","",700,700);
  TH1D* htmp1 = new TH1D("htmp1","",1,0,1);
  cToRemove->cd();
  htmp1->Fill(0.5);
  htmp1->Draw("HE");
  CMS_lumi(cToRemove,"",true,false);
  setTDRStyle();
  // cToRemove->SaveAs("/afs/cern.ch/user/m/mciprian/www/test_plot/test_roofit_toBeRemoved.png");
  // system("rm /afs/cern.ch/user/m/mciprian/www/test_plot/test_roofit_toBeRemoved.png");
  delete htmp1;
  delete cToRemove;

  TCanvas* canvas = new TCanvas("canvas","",600,700);
  canvas->cd();
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->cd();
  // canvas->SetBottomMargin(0.3);
  canvas->SetRightMargin(0.06);

  hvar[0]->Draw("hist");
  if (hvarName_base == "hmT") {
    if (isMuon) hvar[0]->GetXaxis()->SetTitle("W(#mu#nu) transverse mass [GeV]");
    else hvar[0]->GetXaxis()->SetTitle("W(e#nu) transverse mass [GeV]");  
  } else if (hvarName_base == "hlep1pt") {
    if (isMuon) hvar[0]->GetXaxis()->SetTitle("muon p_{T} [GeV]");
    else hvar[0]->GetXaxis()->SetTitle("electron p_{T} [GeV]");  
  }

  if (scaleHistograms) {
    CMS_lumi(canvas,"",true,false);
    hvar[0]->GetYaxis()->SetTitle("a.u.");
  }
  else {
    CMS_lumi(canvas,"",false,false);
    hvar[0]->GetYaxis()->SetTitle("events");
  }
  setTDRStyle();

  Int_t colorList[] = {kBlack, kBlue, kRed, kGreen+2, kOrange+1, kCyan+2, kViolet, kCyan};

  TLegend leg (0.5,0.7,0.9,0.9);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);

  Double_t maxY = -1;
  for (UInt_t i = 0; i < inputFileName.size(); i++) {
    hvar[i]->SetStats(0);
    hvar[i]->SetLineColor(colorList[i]);
    hvar[i]->SetLineWidth(2);
    hvar[i]->Draw("hist same");
    leg.AddEntry(hvar[i],legEntry[i].c_str(),"L");
    if (maxY < hvar[i]->GetBinContent(hvar[i]->GetMaximumBin())) maxY = hvar[i]->GetBinContent(hvar[i]->GetMaximumBin());

  }
  hvar[0]->SetMaximum(1.25 * maxY);

  leg.Draw("same");
  canvas->RedrawAxis("sameaxis");

  canvas->SaveAs((outputDIR + canvasName + ".png").c_str());
  canvas->SaveAs((outputDIR + canvasName + ".pdf").c_str());

  delete canvas;

}


void runCompareShapeVsEta(const string& outputDIR_tmp = "./",
			  const string& canvasName_base = "shapeCompVsEta",
			  const Bool_t isMuon = false,
			  const Bool_t scaleHistograms = true,
			  const Bool_t QCD_enriched_region = false,
			  const Bool_t useDataOnlyForQCD_CR = false,
			  const Bool_t useQCDMonteCarloForMuonsInCR = false,
			  const Int_t plot_comb0_pos1_neg2 = 0
			  ) {
  
  if (useQCDMonteCarloForMuonsInCR) {
    if (not isMuon || not QCD_enriched_region) {
      cout << endl;
      cout << "### WARNING: useQCDMonteCarloForMuonsInCR option is supported only for muons in QCD CR. Exit!" << endl;
      cout << endl;
      exit(EXIT_FAILURE);
    }
  }

  string canvasName = canvasName_base;
  if (QCD_enriched_region) {
    if (useDataOnlyForQCD_CR) canvasName = "QCDfromData_CR_" + canvasName_base;
    else canvasName = "QCDfromDataMinusMC_CR_" + canvasName_base;
    if (isMuon && useQCDMonteCarloForMuonsInCR) canvasName = "QCDfromMC_CR_" + canvasName_base;
  } else {
    canvasName = "WfromMC_SR_" + canvasName_base;
  }


  doCompareQCDshape(outputDIR_tmp, canvasName, isMuon, "hmT", scaleHistograms, QCD_enriched_region, useDataOnlyForQCD_CR, useQCDMonteCarloForMuonsInCR, plot_comb0_pos1_neg2);
  doCompareQCDshape(outputDIR_tmp, canvasName, isMuon, "hlep1pt", scaleHistograms, QCD_enriched_region, useDataOnlyForQCD_CR, useQCDMonteCarloForMuonsInCR, plot_comb0_pos1_neg2);
  
}


//========================


void compareShapeVsEta(const string& outputDIR_tmp = "./",
		       const string& canvasName_base = "shapeCompVsEta",
		       const Bool_t isMuon = false
		       ) {

  for (Int_t i = 0; i < 3; i++) {

    runCompareShapeVsEta(outputDIR_tmp, canvasName_base, isMuon, true, false, false, false, i); // SR MC scaledHist
    runCompareShapeVsEta(outputDIR_tmp, canvasName_base, isMuon, false, false, false, false, i); // SR MC notScaledHist 

    runCompareShapeVsEta(outputDIR_tmp, canvasName_base, isMuon, true, true, false, false, i); // CR data-MC scaledHist 
    runCompareShapeVsEta(outputDIR_tmp, canvasName_base, isMuon, false, true, false, false, i); // CR data-MC notScaledHist 

    runCompareShapeVsEta(outputDIR_tmp, canvasName_base, isMuon, true, true, true, false, i); // CR data only scaledHist 
    runCompareShapeVsEta(outputDIR_tmp, canvasName_base, isMuon, false, true, true, false, i); // CR data only notScaledHist 

    if (isMuon) {
      runCompareShapeVsEta(outputDIR_tmp, canvasName_base, isMuon, true, true, false, true, i); // CR QCD MC scaledHist 
      runCompareShapeVsEta(outputDIR_tmp, canvasName_base, isMuon, false, true, false, true, i); // CR QCD MC notScaledHist 
    }

  }

}
