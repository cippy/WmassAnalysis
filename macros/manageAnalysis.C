#include "TSystem.h"
#include "TROOT.h"

#include <iostream>
#include <string>
#include <cstdlib> //as stdlib.h    
#include <cstdio>

using namespace std;

// usage :
// root -l -q -b 'manageAnalysis.C++("argument1","argument2",...)' &> logfile.log  

// &> logfile.log is useful if you want to save output in a log for debugging
// You call this macro with the same arguments as makeVariableHistograms.C, but using only strings (to build command in an easier way)
// e.g., if you pass a number or a bool, you write these as if they were strings

void manageAnalysis(const Int_t skip_no0_Mu1_ele2 = 0, const Int_t skip_no0_sel1_plot2 = 0) {

  string outputPath = "/afs/cern.ch/user/m/mciprian/www/wmass/analysisPlots/variableHistograms_miniIso/";
    
  string inputDIR  = "/u2/emanuele/TREES_1LEP_80X_V1/";
  string outfileName = "wmass_varhists.root";

  string outputDIR = "";
  string QCD_enriched_region = "";
  string isMuon = "";

  string isWregion = "";

  string command = "";
  string plotCommand = "";
  string qcdStudyCommand = "";
  //string correlationCommand = "";

  string infile_sigReg = "";   
  string infile_bkgReg = ""; 

  gROOT->ProcessLine(".L makeWBosonVariableHistograms.C++");
  gROOT->ProcessLine(".L makeDataMCPlots.C++");
  gROOT->ProcessLine(".L makeQCDstudy.C++");
  //gROOT->ProcessLine(".L makeCorrelationStudy.C++");

  cout << "===========================" << endl;
  cout << " BEGIN!" << endl;
  cout << "===========================" << endl;

  ///////////////////////////////////
  // W region
  isWregion = "true";
  cout << endl;
  cout << "W REGION" << endl;
  cout << endl;

  if (skip_no0_Mu1_ele2 != 2) {

    /////////////////////
    // Wenu
    outputDIR = outputPath + "Wenu/";
    QCD_enriched_region = "false";
    isMuon = "false";
    command = "makeWBosonVariableHistograms(\"" + inputDIR + "\",\"" + outputDIR + "\",\"" + outfileName + "\"," + QCD_enriched_region + "," + isMuon + ")";
    if (skip_no0_sel1_plot2 != 1) {
      cout << command << endl;
      gROOT->ProcessLine(command.c_str());
    }

    //plot with option 0,1,2 to do combined, positive or negative W
    if (skip_no0_sel1_plot2 != 2) {
      for (Int_t i = 0; i < 3; i++) {
	plotCommand = "makeDataMCPlots(\"" + outputDIR + "\",\"" + outfileName + "\"," + isWregion + "," + isMuon + "," + string(Form("%d",i)) + ")";
	cout << plotCommand << endl;
	gROOT->ProcessLine(plotCommand.c_str());
      }
    }
    cout << endl;

    //correlationCommand = "makeCorrelationStudy(\"" + outputDIR + "\",\"" + outputDIR + "\",\"" + outfileName + "\"," + QCD_enriched_region + "," + isMuon + ")";

    //////////////////////
    // Wenu inverted iso cut
    outputDIR = outputPath + "Wenu_invertIsoCut/";
    QCD_enriched_region = "true";
    isMuon = "false";
    command = "makeWBosonVariableHistograms(\"" + inputDIR + "\",\"" + outputDIR + "\",\"" + outfileName + "\"," + QCD_enriched_region + "," + isMuon + ")";
    if (skip_no0_sel1_plot2 != 1) {
      cout << command << endl;
      gROOT->ProcessLine(command.c_str());
    }

    //plot with option 0,1,2 to do combined, positive or negative W
    if (skip_no0_sel1_plot2 != 2) {
      for (Int_t i = 0; i < 3; i++) {
	plotCommand = "makeDataMCPlots(\"" + outputDIR + "\",\"" + outfileName + "\"," + isWregion + "," + isMuon + "," + string(Form("%d",i)) + ")";
	cout << plotCommand << endl;
	gROOT->ProcessLine(plotCommand.c_str());
      }
    }
    cout << endl;

    //////////////////////////////
    // QCD study W->enu
    outputDIR = outputPath + "QCD_study/";
    infile_sigReg = outputPath + "Wenu/wmass_varhists.root"; 
    infile_bkgReg = outputPath + "Wenu_invertIsoCut/wmass_varhists.root";
    isMuon = "false";
    qcdStudyCommand = "makeQCDstudy(\"" + outputDIR + "\",\"" + infile_sigReg + "\",\"" + infile_bkgReg + "\"," + isWregion + "," + isMuon + ")";
    gROOT->ProcessLine(qcdStudyCommand.c_str());


  }

  if (skip_no0_Mu1_ele2 != 1) {

    /////////////////////
    // Wmunu
    outputDIR = outputPath + "Wmunu/";
    QCD_enriched_region = "false";
    isMuon = "true";
    command = "makeWBosonVariableHistograms(\"" + inputDIR + "\",\"" + outputDIR + "\",\"" + outfileName + "\"," + QCD_enriched_region + "," + isMuon + ")";
    if (skip_no0_sel1_plot2 != 1) {
      cout << command << endl;
      gROOT->ProcessLine(command.c_str());
    }

    //plot with option 0,1,2 to do combined, positive or negative W
    if (skip_no0_sel1_plot2 != 2) {
      for (Int_t i = 0; i < 3; i++) {
	plotCommand = "makeDataMCPlots(\"" + outputDIR + "\",\"" + outfileName + "\"," + isWregion + "," + isMuon + "," + string(Form("%d",i)) + ")";
	cout << plotCommand << endl;
	gROOT->ProcessLine(plotCommand.c_str());
      }
    }
    cout << endl;

    //////////////////////
    // Wmunu inverted iso cut
    outputDIR = outputPath + "Wmunu_invertIsoCut/";
    QCD_enriched_region = "true";
    isMuon = "true";
    command = "makeWBosonVariableHistograms(\"" + inputDIR + "\",\"" + outputDIR + "\",\"" + outfileName + "\"," + QCD_enriched_region + "," + isMuon + ")";
    if (skip_no0_sel1_plot2 != 2) {
      cout << command << endl;
      gROOT->ProcessLine(command.c_str());
    }  

    //plot with option 0,1,2 to do combined, positive or negative W
    if (skip_no0_sel1_plot2 != 2) {
      for (Int_t i = 0; i < 3; i++) {
	plotCommand = "makeDataMCPlots(\"" + outputDIR + "\",\"" + outfileName + "\"," + isWregion + "," + isMuon + "," + string(Form("%d",i)) + ")";
	cout << plotCommand << endl;
	gROOT->ProcessLine(plotCommand.c_str());
      }
    }
    cout << endl;

    //////////////////////////////
    // QCD study W->munu
    outputDIR = outputPath + "QCD_study/";
    infile_sigReg = outputPath + "Wmunu/wmass_varhists.root"; 
    infile_bkgReg = outputPath + "Wmunu_invertIsoCut/wmass_varhists.root";
    isMuon = "true";
    qcdStudyCommand = "makeQCDstudy(\"" + outputDIR + "\",\"" + infile_sigReg + "\",\"" + infile_bkgReg + "\"," + isWregion + "," + isMuon + ")";
    gROOT->ProcessLine(qcdStudyCommand.c_str());

  }


  ///////////////////////////////////
  // Z region
  isWregion = "false";
  cout << endl;
  cout << "Z REGION" << endl;
  cout << endl;


  cout << "===========================" << endl;
  cout << " THE END!" << endl;
  cout << "===========================" << endl;



}		    
