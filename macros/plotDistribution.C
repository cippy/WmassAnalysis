#include <stdio.h>
#include <stdlib.h>
#include <cstdlib> //as stdlib.h
#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>      // std::istringstream ; to read array of numbers from a line in a file
#include <string>
#include <vector>
#include <iomanip> //for input/output manipulators
//ROOT header files
#include <TAxis.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TMatrixDSym.h>
#include <TPad.h>
#include <TPaveStats.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TVector2.h>
#include <TVector3.h>
#include <TVirtualFitter.h>

#include <RooStats/RooStatsUtils.h>
#include <RooStats/HLFactory.h>
#include <RooRealVar.h>
#include <RooArgSet.h>
#include <RooDataHist.h>
#include <RooFitResult.h>
#include <RooWorkspace.h>
#include <RooPlot.h>
#include <RooHistPdf.h>
#include <RooAddPdf.h>
#include <RooAbsPdf.h>
#include <RooRandom.h>

using namespace std;

void plotDistribution(const string& inputDIR = "./", const string& outputDIR = "./") {

  if (outputDIR!= "./") system(("mkdir -p "+outputDIR).c_str());

}
