#include "utility.h"

using namespace std;

void doCompareShape(const string& outputDIR_tmp = "./",
		    const string& canvasName_tmp = "QCDshapeComp",
		    const Bool_t isMuon = false,
		    const Int_t tkmet0_pfmet1 = 0,
		    const string& hvarName = "hmT",
		    const Bool_t scaleHistograms = true,
		    const Bool_t QCD_enriched_region = false
		    ) {
  
  if (tkmet0_pfmet1 == 0) return;

  string outputDIR = outputDIR_tmp;
  if (canvasName_tmp.find("WJETS") != string::npos) outputDIR += "WJets_study/";
  else if (canvasName_tmp.find("QCD") != string::npos) outputDIR += "QCD_study/";

  if (isMuon) outputDIR += "wmunu/";
  else outputDIR += "wenu/";

  outputDIR += "shapeComparison/";

  if (QCD_enriched_region) outputDIR += "bkg_region/";
  else outputDIR += "sig_region/";

  if (outputDIR != "./") system(("mkdir -p " + outputDIR).c_str());
  system(("cp "+ PhpToCopy + " " + outputDIR).c_str());

  string canvasName = canvasName_tmp;

  if (hvarName == "hmT") canvasName += "_mT";
  else if (hvarName == "hlep1pt") canvasName += "_lep1pt";

  if (not scaleHistograms) canvasName += "_realNorm";

  if (tkmet0_pfmet1 == 0) canvasName += "_tkmet";
  else canvasName += "_pfmet";

  gROOT->SetBatch(kTRUE);
  TH1::SetDefaultSumw2(); //all the following histograms will automatically call TH1::Sumw2()                                        

  cout << endl;

  vector<TH1*> hvar;
  vector<string> inputFileName;

  // for QCD in electron SR, get QCD as data-MC(noQCD)
  TH1D* htmp = NULL;
  TH1D* hdata_eleSR = NULL; 
  TH1D* hallMCnoQCD_eleSR = NULL;

  string subDir = "";
  if (QCD_enriched_region) {
    if (isMuon) subDir = "Wmunu_QCD_CR/";
    else subDir = "Wenu_QCD_CR/";
  } else {
    if (isMuon) subDir = "Wmunu/";
    else subDir = "Wenu/";
  }

  inputFileName.push_back("/afs/cern.ch/user/m/mciprian/www/wmass/analysisPlots_8TeV/varHist_03July_massWeightOnW_noPFmetCut/" + subDir + "wmass_varhists.root");

  if (tkmet0_pfmet1 == 0) {
    inputFileName.push_back("/afs/cern.ch/user/m/mciprian/www/wmass/analysisPlots_8TeV/variableHistograms_new_v3_nlep1p_HLTbit_eta2p1_tkmet10/" + subDir + "wmass_varhists.root");
    inputFileName.push_back("/afs/cern.ch/user/m/mciprian/www/wmass/analysisPlots_8TeV/variableHistograms_new_v3_nlep1p_HLTbit_eta2p1_tkmet20/" + subDir + "wmass_varhists.root");
    inputFileName.push_back("/afs/cern.ch/user/m/mciprian/www/wmass/analysisPlots_8TeV/variableHistograms_new_v3_nlep1p_HLTbit_eta2p1_tkmet30/" + subDir + "wmass_varhists.root");
  } else {
    inputFileName.push_back("/afs/cern.ch/user/m/mciprian/www/wmass/analysisPlots_8TeV/varHist_03July_massWeightOnW_PFmetCut_SR10_CR20_/" + subDir + "wmass_varhists.root");
    inputFileName.push_back("/afs/cern.ch/user/m/mciprian/www/wmass/analysisPlots_8TeV/varHist_03July_massWeightOnW_PFmetCut_SR20_CR30_/" + subDir + "wmass_varhists.root");
  }


  vector<string> legEntry;
  legEntry.push_back("no MET cut");
  if (tkmet0_pfmet1 == 0) {
    legEntry.push_back("tracker E_{T}^{miss} > 10");
    legEntry.push_back("tracker E_{T}^{miss} > 20");
    legEntry.push_back("tracker E_{T}^{miss} > 30");
  } else {
    if (QCD_enriched_region) {
      legEntry.push_back("PF E_{T}^{miss} < 20");
      legEntry.push_back("PF E_{T}^{miss} < 30");
    } else {
      legEntry.push_back("PF E_{T}^{miss} > 10");
      legEntry.push_back("PF E_{T}^{miss} > 20");
    }
  }

  for (UInt_t i = 0; i < inputFileName.size(); i++) {
    
    TFile* inputFile = new TFile(inputFileName[i].c_str(),"READ");
    if (!inputFile || inputFile->IsZombie()) {
      cout << "Error: signal file not opened. Exit" << endl;
      exit(EXIT_FAILURE);
    }

    if (canvasName_tmp.find("WJETS") != string::npos) {

      if (isMuon) hvar.push_back( (TH1D*) getHistCloneFromFile(inputFile, hvarName, "wmunujets") );
      else hvar.push_back( (TH1D*) getHistCloneFromFile(inputFile, hvarName, "wenujets") );
      checkNotNullPtr(hvar[i],"hvar[i]");

      if (hvarName == "hmT" && QCD_enriched_region) {
	if (isMuon) myRebinHisto(hvar.back(),2);
	else myRebinHisto(hvar.back(),3);
      }
      if (hvarName == "hlep1pt" && QCD_enriched_region) {
	if (isMuon) myRebinHisto(hvar.back(),2);
	else myRebinHisto(hvar.back(),3);
      }

    } else if (canvasName_tmp.find("QCD") != string::npos) {

      if (isMuon) hvar.push_back( (TH1D*) getHistCloneFromFile(inputFile, hvarName, "qcd_mu") );
      //else hvar.push_back( (TH1D*) getHistCloneFromFile(inputFile, hvarName, "qcd_ele") );
      else {
	if (QCD_enriched_region) {
	  cout << "Warning: QCD shape for electrons in CR currently obtained from data without subtracting other backgrounds." << endl;
	  hvar.push_back( (TH1D*) getHistCloneFromFile(inputFile, hvarName, "data_singleEG") );
	} else {
	  //cout << "Warning: QCD shape for electrons in SR currently not implemented. Skipping this round!" << endl;
	  //return;
	  // get data
	  cout << "Warning: QCD shape for electrons in SR currently obtained from data after subtracting other MC backgrounds." << endl;
	  hdata_eleSR = (TH1D*) getHistCloneFromFile(inputFile, hvarName, "data_singleEG");
	  checkNotNullPtr(hdata_eleSR,"hdata_eleSR");
	  // copy structure of data histogram, then reset bin content
	  hallMCnoQCD_eleSR = (TH1D*) hdata_eleSR->Clone("hallMCnoQCD_eleSR");
	  hallMCnoQCD_eleSR->Reset("ICES");
	  // get all MC except QCD
	  vector<string> listAllMCnoQCD = {"wenujets", "wmunujets", "wtaunujets", "zjets", "top", "diboson"};
	  for (UInt_t imc = 0; imc < listAllMCnoQCD.size(); imc++) {
	    htmp = (TH1D*) getHistCloneFromFile(inputFile, hvarName, listAllMCnoQCD[imc]);
	    checkNotNullPtr(htmp,"htmp");
	    hallMCnoQCD_eleSR->Add(htmp);
	  }
	  hvar.push_back(hdata_eleSR);
	  hvar[i]->Add(hallMCnoQCD_eleSR, -1.0);

	}  
      }
      checkNotNullPtr(hvar[i],"hvar[i]");

      if (hvarName == "hmT") {
	if (QCD_enriched_region) {
	  if (isMuon) myRebinHisto(hvar.back(),4);
	  else myRebinHisto(hvar.back(),4);
	} else {
	  if (isMuon) myRebinHisto(hvar.back(),3);
          else myRebinHisto(hvar.back(),3);
	}
      }
      if (hvarName == "hlep1pt" && QCD_enriched_region) {
	if (isMuon) myRebinHisto(hvar.back(),2);
	else myRebinHisto(hvar.back(),2);
      }

    }

    if (scaleHistograms) hvar.back()->Scale(1./hvar.back()->Integral());

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

  if (hvarName == "hmT") {
    if (isMuon) hvar[0]->GetXaxis()->SetTitle("W(#mu#nu) transverse mass [GeV]");
    else hvar[0]->GetXaxis()->SetTitle("W(e#nu) transverse mass [GeV]");  
  } else if (hvarName == "hlep1pt") {
    if (isMuon) hvar[0]->GetXaxis()->SetTitle("muon p_{T} [GeV]");
    else hvar[0]->GetXaxis()->SetTitle("electron p_{T} [GeV]");  
  }
  if (scaleHistograms) {
    hvar[0]->GetYaxis()->SetTitle("a.u.");
  }
  else {
    hvar[0]->GetYaxis()->SetTitle("events");
  }

  string xAxisNames(hvar[0]->GetXaxis()->GetTitle());
  string yAxisNames(hvar[0]->GetYaxis()->GetTitle());

  if (scaleHistograms) draw_nTH1(hvar,xAxisNames, yAxisNames, canvasName, outputDIR, legEntry, "with/no cut",-1.0, 1, false, true);
  else draw_nTH1(hvar,xAxisNames, yAxisNames, canvasName, outputDIR, legEntry, "with/no cut",-1.0, 1, false, false);

  // TCanvas* canvas = new TCanvas("canvas","",600,700);
  // canvas->cd();
  // canvas->SetTickx(1);
  // canvas->SetTicky(1);
  // canvas->cd();
  // // canvas->SetBottomMargin(0.3);
  // canvas->SetRightMargin(0.06);

  // hvar[0]->Draw("hist");
  // if (hvarName == "hmT") {
  //   if (isMuon) hvar[0]->GetXaxis()->SetTitle("W(#mu#nu) transverse mass [GeV]");
  //   else hvar[0]->GetXaxis()->SetTitle("W(e#nu) transverse mass [GeV]");  
  // } else if (hvarName == "hlep1pt") {
  //   if (isMuon) hvar[0]->GetXaxis()->SetTitle("muon p_{T} [GeV]");
  //   else hvar[0]->GetXaxis()->SetTitle("electron p_{T} [GeV]");  
  // }

  // if (scaleHistograms) {
  //   CMS_lumi(canvas,"",true,false);
  //   hvar[0]->GetYaxis()->SetTitle("a.u.");
  // }
  // else {
  //   CMS_lumi(canvas,"",false,false);
  //   hvar[0]->GetYaxis()->SetTitle("events");
  // }
  // setTDRStyle();

  // Int_t colorList[] = {kBlack, kBlue, kRed, kGreen+2, kOrange+1, kCyan+2, kViolet, kCyan};

  // TLegend leg (0.5,0.7,0.9,0.9);
  // leg.SetFillColor(0);
  // leg.SetFillStyle(0);
  // leg.SetBorderSize(0);

  // Double_t maxY = -1;
  // for (UInt_t i = 0; i < inputFileName.size(); i++) {
  //   hvar[i]->SetStats(0);
  //   hvar[i]->SetLineColor(colorList[i]);
  //   hvar[i]->SetLineWidth(2);
  //   hvar[i]->Draw("hist same");
  //   leg.AddEntry(hvar[i],legEntry[i].c_str(),"L");
  //   if (maxY < hvar[i]->GetBinContent(hvar[i]->GetMaximumBin())) maxY = hvar[i]->GetBinContent(hvar[i]->GetMaximumBin());

  // }
  // hvar[0]->SetMaximum(1.25 * maxY);

  // leg.Draw("same");
  // canvas->RedrawAxis("sameaxis");

  // canvas->SaveAs((outputDIR + canvasName + ".png").c_str());
  // canvas->SaveAs((outputDIR + canvasName + ".pdf").c_str());

}


void compareShape(const string& outputDIR_tmp = "./",
		  const string& canvasName = "ShapeComp",
		  const Bool_t isMuon = false
		  ) {
  
  // signal region
  doCompareShape(outputDIR_tmp, canvasName, isMuon, 0, "hmT", true, false);
  doCompareShape(outputDIR_tmp, canvasName, isMuon, 0, "hmT", false, false);
  doCompareShape(outputDIR_tmp, canvasName, isMuon, 1, "hmT", true, false);
  doCompareShape(outputDIR_tmp, canvasName, isMuon, 1, "hmT", false, false);
  doCompareShape(outputDIR_tmp, canvasName, isMuon, 0, "hlep1pt", true, false);
  doCompareShape(outputDIR_tmp, canvasName, isMuon, 0, "hlep1pt", false, false);
  doCompareShape(outputDIR_tmp, canvasName, isMuon, 1, "hlep1pt", true, false);
  doCompareShape(outputDIR_tmp, canvasName, isMuon, 1, "hlep1pt", false, false);

  // background region
  doCompareShape(outputDIR_tmp, canvasName, isMuon, 0, "hmT", true, true);
  doCompareShape(outputDIR_tmp, canvasName, isMuon, 0, "hmT", false, true);
  doCompareShape(outputDIR_tmp, canvasName, isMuon, 1, "hmT", true, true);
  doCompareShape(outputDIR_tmp, canvasName, isMuon, 1, "hmT", false, true);
  doCompareShape(outputDIR_tmp, canvasName, isMuon, 0, "hlep1pt", true, true);
  doCompareShape(outputDIR_tmp, canvasName, isMuon, 0, "hlep1pt", false, true);
  doCompareShape(outputDIR_tmp, canvasName, isMuon, 1, "hlep1pt", true, true);
  doCompareShape(outputDIR_tmp, canvasName, isMuon, 1, "hlep1pt", false, true);


}



void compareShapeVsMet(const string& outputDIR_tmp = "./",
		       const Bool_t isMuon = false
		       ) {

  compareShape(outputDIR_tmp, "QCD_shapeComp", isMuon);
  //compareShape(outputDIR_tmp, "WJETS_shapeComp", isMuon);
  
}
