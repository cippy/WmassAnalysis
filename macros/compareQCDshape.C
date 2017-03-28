#include "utility.h"

using namespace std;

void doCompareQCDshape(const string& outputDIR_tmp = "./",
		       const string& canvasName_tmp = "QCDshapeComp",
		       const Bool_t isMuon = false,
		       const Int_t tkmet0_pfmet1 = 0,
		       const string& hvarName = "hmT",
		       const Bool_t scaleHistograms = true,
		       const Bool_t QCD_enriched_region = false
		       ) {
  
  string outputDIR = outputDIR_tmp + "QCD_study/";
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

  vector<TH1D*> hvar;

  vector<string> inputFileName;

  string subDir = "";
  if (QCD_enriched_region) {
    if (isMuon) subDir = "Wmunu_invertIsoCut/";
    else subDir = "Wenu_invertIsoCut/";
  } else {
    if (isMuon) subDir = "Wmunu/";
    else subDir = "Wenu/";
  }

  inputFileName.push_back("/afs/cern.ch/user/m/mciprian/www/wmass/analysisPlots_8TeV/variableHistograms_new_v3_nlep1p_HLTbit_eta2p1/" + subDir + "wmass_varhists.root");
  if (tkmet0_pfmet1 == 0) {
    inputFileName.push_back("/afs/cern.ch/user/m/mciprian/www/wmass/analysisPlots_8TeV/variableHistograms_new_v3_nlep1p_HLTbit_eta2p1_tkmet10/" + subDir + "wmass_varhists.root");
    inputFileName.push_back("/afs/cern.ch/user/m/mciprian/www/wmass/analysisPlots_8TeV/variableHistograms_new_v3_nlep1p_HLTbit_eta2p1_tkmet20/" + subDir + "wmass_varhists.root");
    inputFileName.push_back("/afs/cern.ch/user/m/mciprian/www/wmass/analysisPlots_8TeV/variableHistograms_new_v3_nlep1p_HLTbit_eta2p1_tkmet30/" + subDir + "wmass_varhists.root");
  } else {
    inputFileName.push_back("/afs/cern.ch/user/m/mciprian/www/wmass/analysisPlots_8TeV/variableHistograms_new_v3_nlep1p_HLTbit_eta2p1_pfmet10/" + subDir + "wmass_varhists.root");
    inputFileName.push_back("/afs/cern.ch/user/m/mciprian/www/wmass/analysisPlots_8TeV/variableHistograms_new_v3_nlep1p_HLTbit_eta2p1_pfmet20/" + subDir + "wmass_varhists.root");
    inputFileName.push_back("/afs/cern.ch/user/m/mciprian/www/wmass/analysisPlots_8TeV/variableHistograms_new_v3_nlep1p_HLTbit_eta2p1_pfmet30/" + subDir + "wmass_varhists.root");
  }


  vector<string> legEntry;
  legEntry.push_back("no MET cut");
  if (tkmet0_pfmet1 == 0) {
    legEntry.push_back("tracker E_{T}^{miss} > 10");
    legEntry.push_back("tracker E_{T}^{miss} > 20");
    legEntry.push_back("tracker E_{T}^{miss} > 30");
  } else {
    legEntry.push_back("PF E_{T}^{miss} > 10");
    legEntry.push_back("PF E_{T}^{miss} > 20");
    legEntry.push_back("PF E_{T}^{miss} > 30");
  }

  for (UInt_t i = 0; i < inputFileName.size(); i++) {
    
    TFile* inputFile = new TFile(inputFileName[i].c_str(),"READ");
    if (!inputFile || inputFile->IsZombie()) {
      cout << "Error: signal file not opened. Exit" << endl;
      exit(EXIT_FAILURE);
    }
    if (isMuon) hvar.push_back( (TH1D*) getHistCloneFromFile(inputFile, hvarName, "qcd_mu") );
    else        hvar.push_back( (TH1D*) getHistCloneFromFile(inputFile, hvarName, "qcd_ele") );
    checkNotNullPtr(hvar[i],"hvar[i]");

    if (hvarName == "hmT") myRebinHisto(hvar.back(),3);
    if (scaleHistograms) hvar.back()->Scale(1./hvar.back()->Integral());

  }

  TCanvas* cToRemove = new TCanvas("cToRemove","",700,700);
  TH1D* htmp1 = new TH1D("htmp1","",1,0,1);
  cToRemove->cd();
  htmp1->Fill(0.5);
  htmp1->Draw("HE");
  CMS_lumi(cToRemove,"",true,false);
  setTDRStyle();
  cToRemove->SaveAs("/afs/cern.ch/user/m/mciprian/www/test_plot/test_roofit_toBeRemoved.png");
  system("rm /afs/cern.ch/user/m/mciprian/www/test_plot/test_roofit_toBeRemoved.png");
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
  if (hvarName == "hmT") {
    if (isMuon) hvar[0]->GetXaxis()->SetTitle("W(#mu#nu) transverse mass [GeV]");
    else hvar[0]->GetXaxis()->SetTitle("W(e#nu) transverse mass [GeV]");  
  } else if (hvarName == "hlep1pt") {
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

}


void compareQCDshape(const string& outputDIR_tmp = "./",
		     const string& canvasName = "QCDshapeComp",
		     const Bool_t isMuon = false
		     ) {
  
  // signal region
  doCompareQCDshape(outputDIR_tmp, canvasName, isMuon, 0, "hmT", true, false);
  doCompareQCDshape(outputDIR_tmp, canvasName, isMuon, 0, "hmT", false, false);
  doCompareQCDshape(outputDIR_tmp, canvasName, isMuon, 1, "hmT", true, false);
  doCompareQCDshape(outputDIR_tmp, canvasName, isMuon, 1, "hmT", false, false);
  doCompareQCDshape(outputDIR_tmp, canvasName, isMuon, 0, "hlep1pt", true, false);
  doCompareQCDshape(outputDIR_tmp, canvasName, isMuon, 0, "hlep1pt", false, false);
  doCompareQCDshape(outputDIR_tmp, canvasName, isMuon, 1, "hlep1pt", true, false);
  doCompareQCDshape(outputDIR_tmp, canvasName, isMuon, 1, "hlep1pt", false, false);

  // background region
  doCompareQCDshape(outputDIR_tmp, canvasName, isMuon, 0, "hmT", true, true);
  doCompareQCDshape(outputDIR_tmp, canvasName, isMuon, 0, "hmT", false, true);
  doCompareQCDshape(outputDIR_tmp, canvasName, isMuon, 1, "hmT", true, true);
  doCompareQCDshape(outputDIR_tmp, canvasName, isMuon, 1, "hmT", false, true);
  doCompareQCDshape(outputDIR_tmp, canvasName, isMuon, 0, "hlep1pt", true, true);
  doCompareQCDshape(outputDIR_tmp, canvasName, isMuon, 0, "hlep1pt", false, true);
  doCompareQCDshape(outputDIR_tmp, canvasName, isMuon, 1, "hlep1pt", true, true);
  doCompareQCDshape(outputDIR_tmp, canvasName, isMuon, 1, "hlep1pt", false, true);


}
