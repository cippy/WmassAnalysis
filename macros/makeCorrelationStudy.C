#include "utility.h"

using namespace std;

string getTexLabel(const string& sampleDir = "wjets") {

  string label = "";
  if (sampleDir == "qcd_ele") label = "QCD (electron enriched)";
  else if (sampleDir == "qcd_mu") label = "QCD (muon enriched)";
  else if (sampleDir == "wjets") label = "W(l#nu)+jets";
  else if (sampleDir == "zjets") label = "Z(ll)+jets";
  else if (sampleDir == "data_doubleEG") label = "data DoubleEG";
  else if (sampleDir == "data_singleEG") label = "data singleEG";
  else {
    cout << "Error in getTexLabel(): directory name " << sampleDir << " unknown, please check. Exit" << endl;
    exit(EXIT_FAILURE);
  }

  return label;

}

void makeCorrelationStudy(const string& outputDIR_tmp = "./", 
			  const string& inputFileName_tmp = "wmass_varhists.root", // file in outputDIR_tmp
			  const string& sampleDir = "wjets",   // select folder in input file where histograms will be taken
			  const Bool_t isWregion = false, 
			  const Bool_t isMuon = false, 
			  const Int_t separatePositiveAndNegative = 0  // can be 0,1,2 for combined, positive or negative
			  ) {

  string hist_suffix = "";

  string subdir = "combined/";
  if (separatePositiveAndNegative == 1) {
    subdir = "positive/";
    hist_suffix = "_plus";
  } else if (separatePositiveAndNegative == 2) {
    subdir = "negative/";
    hist_suffix = "_minus";
  }

  string outputDIR = outputDIR_tmp + subdir + "correlation/" + sampleDir + "/";
  string inputFileName = outputDIR_tmp + inputFileName_tmp;

  system(("mkdir -p " + outputDIR).c_str());
  system(("cp "+ PhpToCopy + " " + outputDIR).c_str());

  string sampleDirTex = getTexLabel(sampleDir);

  gROOT->SetBatch(kTRUE);
  TH1::SetDefaultSumw2(); //all the following histograms will automatically call TH1::Sumw2()                           
  cout << endl;

  TFile* inputFile = new TFile(inputFileName.c_str(),"READ");
  if (!inputFile || inputFile->IsZombie()) {
    cout << "Error: file not opened. Exit" << endl;
    exit(EXIT_FAILURE);
  }

  inputFile->cd();

  vector<string> hvarName;
  vector<string> xAxisname;
  vector<string> yAxisname;

  if (isMuon) {

    hvarName.push_back("h2_mT_lep1relIso04_noIsoCut");
    hvarName.push_back("h2_mT_lep1miniRelIso_noIsoCut");
    hvarName.push_back("h2_mT_tkmet");
    hvarName.push_back("h2_mT_pfmet");
    hvarName.push_back("h2_mT_bosonPt");
    hvarName.push_back("h2_lep1pt_lep1relIso04_noIsoCut");
    hvarName.push_back("h2_lep1pt_lep1miniRelIso_noIsoCut");
    hvarName.push_back("h2_lep1pt_tkmet");
    hvarName.push_back("h2_lep1pt_pfmet");
    hvarName.push_back("h2_lep1pt_bosonPt");

    xAxisname.push_back("W(#mu#nu) transverse mass [GeV]");
    xAxisname.push_back("W(#mu#nu) transverse mass [GeV]");
    xAxisname.push_back("W(#mu#nu) transverse mass [GeV]");
    xAxisname.push_back("W(#mu#nu) transverse mass [GeV]");
    xAxisname.push_back("W(#mu#nu) transverse mass [GeV]");
    xAxisname.push_back("muon p_{T} [GeV]");
    xAxisname.push_back("muon p_{T} [GeV]");
    xAxisname.push_back("muon p_{T} [GeV]");
    xAxisname.push_back("muon p_{T} [GeV]");
    xAxisname.push_back("muon p_{T} [GeV]");

    yAxisname.push_back("muon isolation (relIso04)");
    yAxisname.push_back("muon isolation (miniRelIso)");
    yAxisname.push_back("tracker E_{T}^{miss} [GeV]");
    yAxisname.push_back("PF E_{T}^{miss} [GeV]");
    yAxisname.push_back("boson p_{T} [GeV]");
    yAxisname.push_back("muon isolation (relIso04)");
    yAxisname.push_back("muon isolation (miniRelIso)");
    yAxisname.push_back("tracker E_{T}^{miss} [GeV]");
    yAxisname.push_back("PF E_{T}^{miss} [GeV]");
    yAxisname.push_back("boson p_{T} [GeV]");

  } else {

    hvarName.push_back("h2_mT_lep1relIso03_noIsoCut");
    hvarName.push_back("h2_mT_tkmet");
    hvarName.push_back("h2_mT_pfmet");
    hvarName.push_back("h2_mT_bosonPt");
    hvarName.push_back("h2_lep1pt_lep1relIso03_noIsoCut");
    hvarName.push_back("h2_lep1pt_tkmet");
    hvarName.push_back("h2_lep1pt_pfmet");
    hvarName.push_back("h2_lep1pt_bosonPt");

    xAxisname.push_back("W(e#nu) transverse mass [GeV]");
    xAxisname.push_back("W(e#nu) transverse mass [GeV]");
    xAxisname.push_back("W(e#nu) transverse mass [GeV]");
    xAxisname.push_back("W(e#nu) transverse mass [GeV]");
    xAxisname.push_back("electron p_{T} [GeV]");
    xAxisname.push_back("electron p_{T} [GeV]");
    xAxisname.push_back("electron p_{T} [GeV]");
    xAxisname.push_back("electron p_{T} [GeV]");

    yAxisname.push_back("electron isolation (relIso03)");
    yAxisname.push_back("tracker E_{T}^{miss} [GeV]");
    yAxisname.push_back("PF E_{T}^{miss} [GeV]");
    yAxisname.push_back("boson p_{T} [GeV]");
    yAxisname.push_back("electron isolation (relIso03)");
    yAxisname.push_back("tracker E_{T}^{miss} [GeV]");
    yAxisname.push_back("PF E_{T}^{miss} [GeV]");
    yAxisname.push_back("boson p_{T} [GeV]");

  }


  TH2D* h2tmp = new TH2D("h2tmp","",2,0,2,2,0,2);
  h2tmp->Fill(0.8,0.8,4);
  h2tmp->Fill(0.5,1.2,2);
  h2tmp->Fill(1.2,1.1,5);
  h2tmp->Fill(1.8,0.5,3);
  drawCorrelationPlot(h2tmp, "variable 1", "variable 2", "tmpToBeRemoved", "tmp object", outputDIR);
  system(("rm " + outputDIR + "*tmpToBeRemoved*").c_str());

  for (UInt_t i = 0; i < hvarName.size(); i++) {
    
    hvarName[i] = hvarName[i] + hist_suffix;

    TH2D* h2var = NULL;
    h2var = (TH2D*) getHist2CloneFromFile(inputFile, hvarName[i], sampleDir);
    checkNotNullPtr(h2var,"h2var");
    
    string prunedhvarName = "";
    prunedhvarName.assign(hvarName[i],2,string::npos);
    
    drawCorrelationPlot(h2var, xAxisname[i], yAxisname[i], "correlation_"+sampleDir+prunedhvarName,getTexLabel(sampleDir),outputDIR);
    
  }

  /*
  drawCorrelationPlot(h2_mT_mZ, "transverse mass(e^{+}e^{-}) [GeV]", "invariant mass(e^{+}e^{-}) [GeV]", "correlation_dataDoubleEG_mT_mZ","data DoubleEG",outputDIR);
  drawCorrelationPlot(h2_mT_lep1pt, "transverse mass(e^{+}e^{-}) [GeV]", "leading electron p_{T} [GeV]", "correlation_dataDoubleEG_mT_lep1pt","data DoubleEG",outputDIR);
  drawCorrelationPlot(h2_mT_lep1eta, "transverse mass(e^{+}e^{-}) [GeV]", "leading electron #eta", "correlation_dataDoubleEG_mT_lep1eta","data DoubleEG",outputDIR);
  drawCorrelationPlot(h2_mT_lep1sigIetaIeta, "transverse mass(e^{+}e^{-}) [GeV]", "leading electron #sigma_{i#etai#eta}", "correlation_dataDoubleEG_mT_lep1sigIetaIeta","d\
ata DoubleEG",outputDIR);
  drawCorrelationPlot(h2_mT_lep1r9, "transverse mass(e^{+}e^{-}) [GeV]", "leading electron R9", "correlation_dataDoubleEG_mT_lep1r9","data DoubleEG",outputDIR);
  drawCorrelationPlot(h2_mT_pfmet, "transverse mass(e^{+}e^{-}) [GeV]", "PF E_{T}^{miss} [GeV]", "correlation_dataDoubleEG_mT_pfmet","data DoubleEG",outputDIR);
  drawCorrelationPlot(h2_mT_tkmet, "transverse mass(e^{+}e^{-}) [GeV]", "Tracker E_{T}^{miss} [GeV]", "correlation_dataDoubleEG_mT_tkmet","data DoubleEG",outputDIR);
  drawCorrelationPlot(h2_mT_bosonPt, "transverse mass(e^{+}e^{-}) [GeV]", "boson p_{T} [GeV]", "correlation_dataDoubleEG_mT_bosonPt","data DoubleEG",outputDIR);
  drawCorrelationPlot(h2_mT_lep1relIso03, "transverse mass(e^{+}e^{-}) [GeV]", "leading electron isolation (relIso03)", "correlation_dataDoubleEG_mT_lep1relIso03","data D\
oubleEG",outputDIR);
  drawCorrelationPlot(h2_mT_lep1relIso04, "transverse mass(e^{+}e^{-}) [GeV]", "leading electron isolation (relIso04)", "correlation_dataDoubleEG_mT_lep1relIso04","data D\
oubleEG",outputDIR);

  // cout << "Check" << endl;                                                                                                                                              

  drawCorrelationPlot(h2_lep1pt_mZ, "leading electron p_{T} [GeV]", "invariant mass(e^{+}e^{-}) [GeV]", "correlation_dataDoubleEG_lep1pt_mZ","data DoubleEG",outputDIR);
  drawCorrelationPlot(h2_lep1pt_lep1eta, "leading electron p_{T} [GeV]", "leading electron #eta", "correlation_dataDoubleEG_lep1pt_lep1eta","data DoubleEG",outputDIR);
  drawCorrelationPlot(h2_lep1pt_lep1sigIetaIeta, "leading electron p_{T} [GeV]", "leading electron #sigma_{i#etai#eta}", "correlation_dataDoubleEG_lep1pt_lep1sigIetaIeta"\
		      ,"data DoubleEG",outputDIR);
  drawCorrelationPlot(h2_lep1pt_lep1r9, "leading electron p_{T} [GeV]", "leading electron R9", "correlation_dataDoubleEG_lep1pt_lep1r9","data DoubleEG",outputDIR);
  drawCorrelationPlot(h2_lep1pt_pfmet, "leading electron p_{T} [GeV]", "PF E_{T}^{miss} [GeV]", "correlation_dataDoubleEG_lep1pt_pfmet","data DoubleEG",outputDIR);
  drawCorrelationPlot(h2_lep1pt_tkmet, "leading electron p_{T} [GeV]", "Tracker E_{T}^{miss} [GeV]", "correlation_dataDoubleEG_lep1pt_tkmet","data DoubleEG",outputDIR);
  drawCorrelationPlot(h2_lep1pt_bosonPt, "leading electron p_{T} [GeV]", "boson p_{T} [GeV]", "correlation_dataDoubleEG_lep1pt_bosonPt","data DoubleEG",outputDIR);
  drawCorrelationPlot(h2_lep1pt_lep1relIso03, "leading electron p_{T} [GeV]", "leading electron isolation (relIso03)", "correlation_dataDoubleEG_lep1pt_lep1relIso03","dat\
a DoubleEG",outputDIR);

  */


}
