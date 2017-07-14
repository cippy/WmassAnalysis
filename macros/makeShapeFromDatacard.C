#include "utility.h"

using namespace std;

//============================================================


//============================================================

void searchMinMaxIntegralHistogram(string& minHistName, string& maxHistName, TFile* inputFile = NULL, const string& histNamePattern = "") {

  gROOT->SetBatch(kTRUE);
  TH1::SetDefaultSumw2(); //all the following histograms will automatically call TH1::Sumw2()                                                          

  if (inputFile == NULL) {
    cout << "Error in searchMinMaxIntegralHistogram(): inputFile is NULL. please check. Exit ..." << endl;
    exit(EXIT_FAILURE);
  }

  inputFile->cd();

  vector<TH1D*> histFromFile;
  TH1D* hvar = NULL;

  TList* list = inputFile->GetListOfKeys() ;
  if (!list) { 
    cout << "Error in searchMinMaxIntegralHistogram(). No keys found in file" << endl;  
    exit(EXIT_FAILURE); 
  }
  TIter next(list) ;
  TKey* key ;
  TObject* obj ;    
  while ( (key = (TKey*)next()) ) {
    obj = key->ReadObj() ;
    if ( obj->InheritsFrom("TH1")) {
  	string hname(obj->GetName());
	if (histNamePattern == "" || hname.find(histNamePattern) != string::npos) {
	  histFromFile.push_back( (TH1D*) getHistCloneFromFile(inputFile, hname, ""));
	  checkNotNullPtr(histFromFile.back(),"histFromFile.back(): " + hname);
	  //cout << histFromFile.back()->GetName() << endl;
	}
    }
  }
  
  Double_t minIntegral = 1e34;
  Double_t maxIntegral = -1.;
  // Int_t indexMin = -1;
  // Int_t indexMax = -1;

  for (UInt_t i = 0; i < histFromFile.size(); i++) {

    if ( minIntegral > histFromFile[i]->Integral() ) {
      minIntegral = histFromFile[i]->Integral();
      minHistName = string(histFromFile[i]->GetName());
      //indexMin = i;
    }

    if (maxIntegral < histFromFile[i]->Integral()) {
      maxIntegral = histFromFile[i]->Integral();
      maxHistName = string(histFromFile[i]->GetName());
      //indexMax = i;
    }

  }

  // cout << "Histo with max integral --> " << maxHistName << " : " << maxIntegral << endl;
  // cout << "Histo with min integral --> " << minHistName << " : " << minIntegral << endl;

}

//============================================================


TH1* runShapeFromDatacard(//TH1* htoreturn = NULL,
			  const string& inputDIR = "./",
			  const string& inputFileName = "./", // nominal mass
			  const string& outputDIR = "./",
			  const string& canvasName = "",
			  const string& hnominalName = "x_W",   // name of nominal histogram
			  const string& hvarNamePattern = "x_W_CMS_We_pdf",  // patter to select variated histogram
			  const string& legendSyst = "pdf"   // name of syst changed up or down, pdf, recoil, electron scale, ...
			  ) {

  cout << endl;
  //cout << "================================================" << endl;
  cout << endl;
  
  gROOT->SetBatch(kTRUE);
  TH1::SetDefaultSumw2(); //all the following histograms will automatically call TH1::Sumw2()                                        
  TH1::AddDirectory(kFALSE);

  TFile* inputFile = new TFile((inputDIR+inputFileName).c_str(),"READ");
  if (!inputFile || inputFile->IsZombie()) {
    cout << "Error: file not opened. Exit" << endl;
    exit(EXIT_FAILURE);
  }
  inputFile->cd();

  TH1D* hmin = NULL;
  TH1D* hmax = NULL;
  TH1D* hnominal = NULL;

  //string hnominalName = "x_W";
  string minHistName = "";
  string maxHistName = "";

  searchMinMaxIntegralHistogram(minHistName, maxHistName, inputFile, hvarNamePattern);

  hmin = (TH1D*) getHistCloneFromFile(inputFile, minHistName, "");
  hmax = (TH1D*) getHistCloneFromFile(inputFile, maxHistName, "");
  checkNotNullPtr(hmin,"hmin");
  checkNotNullPtr(hmax,"hmax");
  // cout << "Histo with max integral --> " << maxHistName << " : " << hmax->Integral() << endl;
  // cout << "Histo with min integral --> " << minHistName << " : " << hmin->Integral() << endl;

  // TH1D htmp( *((TH1D*) getHistCloneFromFile(inputFile, hnominalName, "")) );
  // hnominal = new TH1D(htmp);
  hnominal = (TH1D*) getHistCloneFromFile(inputFile, hnominalName, "");
  checkNotNullPtr(hnominal,"hnominal");
  //TH1D* htmp = (TH1D*) getHistCloneFromFile(inputFile, hnominalName, "");
  // TH1D* htmp = new TH1D(hnominal->GetName(),"",hnominal->GetNbinsX(),hnominal->GetBinLowEdge(1), hnominal->GetBinLowEdge(hnominal->GetNbinsX()+1));

  // htoreturn = new TH1D();
  // for (Int_t i = 1; i <= hnominal->GetNbinsX(); i++) {
  //   htoreturn->SetBinContent(i,hnominal->GetBinContent(i));
  //   htoreturn->SetBinError(i,hnominal->GetBinError(i));
  // }

  vector<TH1*> vecHist1d;
  vector<string> vecLegEntry;
  vecHist1d.push_back(hnominal); vecLegEntry.push_back("nominal");
  vecHist1d.push_back(hmin);     (minHistName.find("Down") == string::npos) ? vecLegEntry.push_back(legendSyst + " Up") : vecLegEntry.push_back(legendSyst + " Down");
  vecHist1d.push_back(hmax);     (maxHistName.find("Down") == string::npos) ? vecLegEntry.push_back(legendSyst + " Up") : vecLegEntry.push_back(legendSyst + " Down");

  //=================================================
  // as usual, must call draw_* function at least once, otherwise plots are screwed up. This depends on the settings in CMS_Lumi.h, used by the function
  TH1D* htmp1 = new TH1D("htmp1","",1,0,1);
  TH1D* htmp2 = new TH1D("htmp2","",1,0,1);
  htmp1->Fill(0.5);
  htmp2->Fill(0.5);
  vector<TH1*> htmpVec; htmpVec.push_back(htmp1); htmpVec.push_back(htmp2);
  vector<string> legtmp; legtmp.push_back("htmp1"); legtmp.push_back("htmp2");
  draw_nTH1(htmpVec, "variable", "Events", "tmpToBeRemoved", "/afs/cern.ch/user/m/mciprian/www/test_plot/", legtmp, "var/nominal", -1.0, 1, false);
  //system("rm /afs/cern.ch/user/m/mciprian/www/test_plot/*tmpToBeRemoved*");
  //=================================================

  draw_nTH1(vecHist1d, "electron p_{T} [GeV]", "Events", canvasName, outputDIR, vecLegEntry, "var/nominal", -1.0, 1, false, true);

  cout << endl;
  //  cout << endl;
  inputFile->Close();
  delete inputFile;

  // checkNotNullPtr(htmp,"htmp");
  // return htmp;

  return hnominal;

}


//=============================================

void callShapeFromDatacard(const string& inputFileDIR = "./",
			   const string& inputFileName = "./",  // nominal mass
			   const string& inputFileNameMassUp = "./",  // mass Up
			   const string& inputFileNameMassDown = "./", // mass Down
			   const string& outputDIR = "./",
			   const Int_t massShift = 0,
			   Double_t *chiSq_massUpDown = NULL,
			   Int_t *chiSq_massUpDown_NDF = NULL
			   )
{

  // test outputDIR --> /afs/cern.ch/user/m/mciprian/www/test_plot/

  // https://github.com/emanueledimarco/cmg-cmssw/blob/wmass_53X/CMGTools/WMass/python/tools/eventVars_wmass.py#L41
  //  self.wmass_steps = [x for x in range(0,10,2)] + [x for x in range(10,25,5)] + [x for x in range(25,55,10)] + [x for x in range(55,141,20)]

  // name from file name subtracting last part
  string canvasNameNominal = inputFileName.substr(0, inputFileName.find(".input.root"));
  string canvasNameMassUp = inputFileNameMassUp.substr(0, inputFileNameMassUp.find(".input.root"));
  string canvasNameMassDown = inputFileNameMassDown.substr(0, inputFileNameMassDown.find(".input.root"));

  createPlotDirAndCopyPhp(outputDIR);

  // W pdf
  TH1D* hnominal = (TH1D*) runShapeFromDatacard(inputFileDIR, inputFileName, outputDIR, canvasNameNominal+"_Wpdf", "x_W","x_W_CMS_We_pdf","pdf");
  TH1D* hMassUp = (TH1D*) runShapeFromDatacard(inputFileDIR, inputFileNameMassUp, outputDIR, canvasNameMassUp+"_Wpdf", "x_W","x_W_CMS_We_pdf","pdf");
  TH1D* hMassDown = (TH1D*) runShapeFromDatacard(inputFileDIR, inputFileNameMassDown, outputDIR, canvasNameMassDown+"_Wpdf", "x_W","x_W_CMS_We_pdf","pdf");
  // TH1D* hnominal = NULL;
  // TH1D* hMassUp = NULL;
  // TH1D* hMassDown = NULL;
  // runShapeFromDatacard(hnominal, inputFileDIR, inputFileName, outputDIR, canvasNameNominal+"_Wpdf", "x_W","x_W_CMS_We_pdf","pdf");
  // runShapeFromDatacard(hMassUp, inputFileDIR, inputFileNameMassUp, outputDIR, canvasNameMassUp+"_Wpdf", "x_W","x_W_CMS_We_pdf","pdf");
  // runShapeFromDatacard(hMassDown, inputFileDIR, inputFileNameMassDown, outputDIR, canvasNameMassDown+"_Wpdf", "x_W","x_W_CMS_We_pdf","pdf");
  checkNotNullPtr(hnominal,"hnominal");
  checkNotNullPtr(hMassUp,"hMassUp");
  checkNotNullPtr(hMassDown,"hMassDown");

  // W recoil scale
  runShapeFromDatacard(inputFileDIR, inputFileName, outputDIR, canvasNameNominal+"_Welescale", "x_W","x_W_CMS_We_elescale","e scale");
  runShapeFromDatacard(inputFileDIR, inputFileNameMassUp, outputDIR, canvasNameMassUp+"_Welescale", "x_W","x_W_CMS_We_elescale","e scale");
  runShapeFromDatacard(inputFileDIR, inputFileNameMassDown, outputDIR, canvasNameMassDown+"_Welescale", "x_W","x_W_CMS_We_elescale","e scale");
  // W electron scale
  runShapeFromDatacard(inputFileDIR, inputFileName, outputDIR, canvasNameNominal+"_Wrecoil", "x_W","x_W_CMS_We_recoil","recoil");
  runShapeFromDatacard(inputFileDIR, inputFileNameMassUp, outputDIR, canvasNameMassUp+"_Wrecoil", "x_W","x_W_CMS_We_recoil","recoil");
  runShapeFromDatacard(inputFileDIR, inputFileNameMassDown, outputDIR, canvasNameMassDown+"_Wrecoil", "x_W","x_W_CMS_We_recoil","recoil");
  // QCD fake rate: norm
  runShapeFromDatacard(inputFileDIR, inputFileName, outputDIR, canvasNameNominal+"_QCD_FR_Norm", "x_data_fakes", "x_data_fakes_CMS_We_FRe_norm", "QCD FR norm");
  runShapeFromDatacard(inputFileDIR, inputFileNameMassUp, outputDIR, canvasNameMassUp+"_QCD_FR_Norm", "x_data_fakes", "x_data_fakes_CMS_We_FRe_norm", "QCD fake FR norm");
  runShapeFromDatacard(inputFileDIR, inputFileNameMassDown, outputDIR, canvasNameMassDown+"_QCD_FR_Norm", "x_data_fakes", "x_data_fakes_CMS_We_FRe_norm", "QCD FR norm");

  // QCD fake rate: pt
  runShapeFromDatacard(inputFileDIR, inputFileName, outputDIR, canvasNameNominal+"_QCD_FR_pt", "x_data_fakes", "x_data_fakes_CMS_We_FRe_pt", "QCD FR pt");
  runShapeFromDatacard(inputFileDIR, inputFileNameMassUp, outputDIR, canvasNameMassUp+"_QCD_FR_pt", "x_data_fakes", "x_data_fakes_CMS_We_FRe_pt", "QCD fake FR pt");
  runShapeFromDatacard(inputFileDIR, inputFileNameMassDown, outputDIR, canvasNameMassDown+"_QCD_FR_pt", "x_data_fakes", "x_data_fakes_CMS_We_FRe_pt", "QCD FR pt");

  vector<TH1*> histList; histList.push_back(hnominal); histList.push_back(hMassUp); histList.push_back(hMassDown);
  vector<string> legEntryList; legEntryList.push_back("M_{0}"); legEntryList.push_back(Form("M_{0} + %d MeV",massShift)); legEntryList.push_back(Form("M_{0} - %d MeV",massShift));

  // for (UInt_t i = 0; i < histList.size(); i++) {
  //   cout << "histList[i]->Integral(): " << histList[i]->Integral() << endl;
  // }

  if (inputFileName.find("_pos_") != string::npos) {
    draw_nTH1(histList, "electron p_{T} [GeV]", "Events", Form("comparisonMassVariation_Wplus_%dMeV",massShift), outputDIR, legEntryList, "var/nominal", -1.0, 1, false, true);
  } else {
    draw_nTH1(histList, "electron p_{T} [GeV]", "Events", Form("comparisonMassVariation_Wminus_%dMeV",massShift), outputDIR, legEntryList, "var/nominal", -1.0, 1, false, true);
  }

  Double_t myChiSquare_mw_up = 0.0;
  Double_t myChiSquare_mw_down = 0.0;
  for (Int_t i = 1; i <= hnominal->GetNbinsX(); i++) {
    myChiSquare_mw_up   += pow( (hnominal->GetBinContent(i) - hMassUp->GetBinContent(i))/hnominal->GetBinError(i), 2);
    myChiSquare_mw_down += pow( (hnominal->GetBinContent(i) - hMassUp->GetBinContent(i))/hnominal->GetBinError(i), 2);
  }

  chiSq_massUpDown[0] = myChiSquare_mw_up;
  chiSq_massUpDown[1] = myChiSquare_mw_down;
  chiSq_massUpDown_NDF[0] = hnominal->GetNbinsX();
  chiSq_massUpDown_NDF[1] = hnominal->GetNbinsX();

  //cout << endl;
  cout << endl;

  // TFile* inputFileMassUp = new TFile(inputFileNameMassUp.c_str(),"READ");
  // if (!inputFileMassUp || inputFileMassUp->IsZombie()) {
  //   cout << "Error: file not opened. Exit" << endl;
  //   exit(EXIT_FAILURE);
  // }


  // TFile* inputFileMassDown = new TFile(inputFileNameMassDown.c_str(),"READ");
  // if (!inputFileMassDown || inputFileMassDown->IsZombie()) {
  //   cout << "Error: file not opened. Exit" << endl;
  //   exit(EXIT_FAILURE);
  // }

  // inputFileMassUp->cd();
  // TH1D* hnominalMassUp = NULL;
  // hnominalMassUp = (TH1D*) getHistCloneFromFile(inputFile, , "");
  // checkNotNullPtr(hnominalMassUp,"hnominalMassUp");


}

//=============================================================================

void makeShapeFromDatacard(const string& inputFileDIR_base = "./",
			   const string& outputDIR_base = "./"
			   )
{

  createPlotDirAndCopyPhp(outputDIR_base);

  vector<Double_t> etaBinEdges;
  etaBinEdges.push_back(0.0);
  etaBinEdges.push_back(1.0);
  etaBinEdges.push_back(1.5);
  etaBinEdges.push_back(2.5);
  Int_t nEtaBins = (Int_t) (etaBinEdges.size() -1);  // 3 bins numbered as 0, 1, 2 


  vector<string> etaBinEdges_str;
  for (UInt_t ietabin = 0; ietabin < etaBinEdges.size(); ietabin++) {
    etaBinEdges_str.push_back(getStringFromDouble(etaBinEdges[ietabin]));
  }    

  //==============   WARNING  =================
  // must be consistent with mass values assigned here: 
  // /afs/cern.ch/work/m/mciprian/w_mass_analysis/CMSSW_5_3_22_patch1/src/CMGTools/WMass/python/tools/eventVars_wmass.py

  UInt_t howManyMassIndex = 2;  // hardcoded temporary number to decide how many index around the nominal must be considered 
  // usually they goes from 0 to N, where N is an even number, but we can have central value X=11 and we may want to consider only 2 variations, from X-2 to X+2
  // in this case, set howManyMassIndex = 2
  vector<Bool_t> indexIsGood;  // store if index must be considered

  // shift in MeV
  vector<Int_t> tmp_wMassShift;
  for(Int_t i = 0; i < 24; i += 2) {
    tmp_wMassShift.push_back(i);
  }
  for(Int_t i = 24; i < 54; i += 10) {
    tmp_wMassShift.push_back(i);
  }
  for(Int_t i = 54; i < 141; i += 20) {
    tmp_wMassShift.push_back(i);
  }

  for (UInt_t i = 0; i < tmp_wMassShift.size(); i++) {
    if (i <= howManyMassIndex) indexIsGood.push_back(true);
    else indexIsGood.push_back(false);
  }
  
  vector<Int_t> wMassShift;
  UInt_t tmp_wMassShift_size = tmp_wMassShift.size();
  for(UInt_t i = 0; i < tmp_wMassShift_size; i++) {
    wMassShift.push_back(-1 * tmp_wMassShift[tmp_wMassShift_size -1 - i]);
  }
  // don't double count 0
  for(UInt_t i = 1; i < tmp_wMassShift_size; i++) {  
    wMassShift.push_back(tmp_wMassShift[i]);
  }

  vector<string> chargeID_str; chargeID_str.push_back("pos"); chargeID_str.push_back("neg");

  if (wMassShift.size()%2 == 0) {
    cout << "Warning, should have an odd number of mass variations (including the nominal), while wMassShift.size() = " << wMassShift.size() << endl;
    exit(EXIT_FAILURE);
  } 
  Int_t nominalMass_index = (wMassShift.size() -1 ) / 2;
  stringstream ss_mwNominal;  ss_mwNominal << nominalMass_index;

  map<string, Double_t> myChiSquare_mass_charge_eta;
  map<string, Int_t> myChiSquare_mass_charge_eta_NDF;
  map<string, Double_t> myChiSquare_mass;
  map<string, Int_t> myChiSquare_mass_NDF;

  for (UInt_t imw = 0; imw < wMassShift.size(); imw++) {

    stringstream ss_mw;  ss_mw   << imw;
    
    for (UInt_t icharge = 0; icharge < chargeID_str.size(); icharge++) {    
      
      for (Int_t ietabin = 0; ietabin < nEtaBins; ietabin++) {

	string etaLow = etaBinEdges_str[ietabin];
	string etaUp = etaBinEdges_str[ietabin+1];

	myChiSquare_mass_charge_eta["wenu_mass" + ss_mw.str()     + "_" + chargeID_str[icharge] + "_eta_" + etaLow + "_" + etaUp] = 0.0;
	myChiSquare_mass_charge_eta_NDF["wenu_mass" + ss_mw.str()     + "_" + chargeID_str[icharge] + "_eta_" + etaLow + "_" + etaUp] = 0;

      }

    }

    myChiSquare_mass[ss_mw.str()]   = 0.0;	
    myChiSquare_mass_NDF[ss_mw.str()]   = 0;	


  }


  for (Int_t imw = 1; imw <= nominalMass_index && indexIsGood[imw]; imw++) {

    stringstream ss_mw_nom;  ss_mw_nom   << nominalMass_index;
    stringstream ss_mw_up;   ss_mw_up   << (nominalMass_index + imw);
    stringstream ss_mw_down; ss_mw_down << (nominalMass_index - imw);
    
    for (UInt_t icharge = 0; icharge < chargeID_str.size(); icharge++) {    
      
      for (Int_t ietabin = 0; ietabin < nEtaBins; ietabin++) {
      //for (UInt_t ietabin = 2; ietabin < 3; ietabin++) {

	//if ((nominalMass_index - imw) == 7 && icharge == 0 && ietabin == 1) continue;     

	string etaLow = etaBinEdges_str[ietabin];
	string etaUp = etaBinEdges_str[ietabin+1];

	string inputFileDIR = inputFileDIR_base + "eta_" + etaLow + "_" + etaUp + "/";
	string outputDIR    = outputDIR_base    + "eta_" + etaLow + "_" + etaUp + "/";
	createPlotDirAndCopyPhp(outputDIR);

	string inputFileName         = "wenu_mass" + ss_mw_nom.str()  + "_" + chargeID_str[icharge] + "_eta_" + etaLow + "_" + etaUp + ".input.root";
	string inputFileNameMassUp   = "wenu_mass" + ss_mw_up.str()   + "_" + chargeID_str[icharge] + "_eta_" + etaLow + "_" + etaUp + ".input.root";
	string inputFileNameMassDown = "wenu_mass" + ss_mw_down.str() + "_" + chargeID_str[icharge] + "_eta_" + etaLow + "_" + etaUp + ".input.root";

	Double_t chiSq_massUpDown[2] = {0.0, 0.0};
	Int_t chiSq_massUpDown_NDF[2] = {0, 0};

	callShapeFromDatacard(inputFileDIR, 
			      inputFileName, inputFileNameMassUp, inputFileNameMassDown, 
			      outputDIR, 
			      wMassShift[nominalMass_index + imw], 
			      chiSq_massUpDown,
			      chiSq_massUpDown_NDF);

	myChiSquare_mass_charge_eta["wenu_mass" + ss_mw_up.str()     + "_" + chargeID_str[icharge] + "_eta_" + etaLow + "_" + etaUp ] = chiSq_massUpDown[0];
	myChiSquare_mass_charge_eta["wenu_mass" + ss_mw_down.str()   + "_" + chargeID_str[icharge] + "_eta_" + etaLow + "_" + etaUp ] = chiSq_massUpDown[1];
	myChiSquare_mass_charge_eta_NDF["wenu_mass" + ss_mw_up.str()     + "_" + chargeID_str[icharge] + "_eta_" + etaLow + "_" + etaUp ]   = chiSq_massUpDown_NDF[0];	
	myChiSquare_mass_charge_eta_NDF["wenu_mass" + ss_mw_down.str()   + "_" + chargeID_str[icharge] + "_eta_" + etaLow + "_" + etaUp ]   = chiSq_massUpDown_NDF[1];	
    
	myChiSquare_mass[ss_mw_up.str()]   += chiSq_massUpDown[0]; 
	myChiSquare_mass[ss_mw_down.str()] += chiSq_massUpDown[1];
	myChiSquare_mass_NDF[ss_mw_up.str()]   += chiSq_massUpDown_NDF[0];	
	myChiSquare_mass_NDF[ss_mw_down.str()] += chiSq_massUpDown_NDF[1];	

      }

    }

  }

  TGraphErrors * gChiSquareVsMassShift = new TGraphErrors();
  
  for (UInt_t imw = 0; imw < wMassShift.size(); imw++) {

    stringstream ss_mw;  ss_mw << imw;
  
    gChiSquareVsMassShift->SetPoint(imw, wMassShift[imw], myChiSquare_mass[ss_mw.str()]/((Double_t)myChiSquare_mass_NDF[ss_mw.str()]));
    //gChiSquareVsMassShift->SetPointError(imw, wMassShift[imw], myChiSquare_mass[ss_mw.str()]/((Double_t)myChiSquare_mass_NDF[ss_mw.str()]));

  }

  Double_t Chi2Sigma = sqrt(2.0/(Double_t)myChiSquare_mass_NDF["0"]);
  cout <<"Chi2 NDF = " << myChiSquare_mass_NDF["0"] << endl;
  cout <<"Chi2Sigma = " << Chi2Sigma << endl;
  TH1F* horizLine = new TH1F("horizLine","",400,-200,200);
  for (Int_t i = 1; i <= horizLine->GetNbinsX(); i++) {
    horizLine->SetBinContent(i, 1.0);
    horizLine->SetBinError(i, Chi2Sigma);
  }
  horizLine->SetStats(0);
  horizLine->SetLineColor(kBlack);
  horizLine->SetLineWidth(2);
  horizLine->SetFillColor(kGray);

  TCanvas* c = new TCanvas("c","");
  gChiSquareVsMassShift->Draw("AL*");
  horizLine->Draw("E2same");
  //gChiSquareVsMassShift->Draw("AL* same");
  c->RedrawAxis("sameaxis");
  gChiSquareVsMassShift->GetXaxis()->SetTitle("m_{W}(var) - m_{W}(nominal) [MeV]");
  gChiSquareVsMassShift->GetYaxis()->SetTitle("#chi^{2}");
  gChiSquareVsMassShift->GetYaxis()->SetRangeUser(0.0, 1.1);
  TF1* line = new TF1("horiz_line","1",horizLine->GetXaxis()->GetBinLowEdge(1),horizLine->GetXaxis()->GetBinLowEdge(horizLine->GetNbinsX()+1));
  line->SetLineColor(kRed);
  line->SetLineWidth(2);
  line->Draw("Lsame");
  c->SaveAs((outputDIR_base + "chiSquareVsMassShift.png").c_str()); 
  c->SaveAs((outputDIR_base + "chiSquareVsMassShift.pdf").c_str()); 
  delete c;

  cout << endl;
  cout << "The end" << endl;
  cout << endl;

}
