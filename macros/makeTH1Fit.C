#include "utility.h"

#include <RooCBShape.h>
#include <RooGaussian.h>
#include <RooGlobalFunc.h>

using namespace std;

static Bool_t debug = true;
static Bool_t fitWithGaussian = true;
static Double_t minFitGaus = 35.0;
static Double_t maxFitGaus = 80.0;
static Double_t minFitCB = 20.0;
static Double_t maxFitCB  = 130.0;

///////////////////////////
// NOTE: for the cystal ball, alpha > 0 (< 0) --> non gaussian tail is on the left (right)


//============================================================                                                                                                              
void makeRooFit(TH1* hist) {


  // Int_t binStartFit = hist->GetBin(minFit);
  // Int_t binEndFit = hist->GetBin(maxFit);

  // get binning from histograms
  vector<double> varBinning;
  for (int i = 1; i < hist->GetNbinsX()+2; i++) {                                                                                                       
    varBinning.push_back(hist->GetXaxis()->GetBinLowEdge(i));
  }

  RooBinning *binning = new RooBinning(varBinning.size()-1,&varBinning[0],"varBinning");                                      
  RooRealVar observable ("mT","transverse mass [GeV]",varBinning.front(),varBinning.back());
  observable.setBinning(*binning);
  if(debug) observable.Print();
    
  // data-histogram                                                                                                                                        
                                          
  RooArgList observable_list (observable);
  RooDataHist qcdhist("qcdhist","",observable_list, hist);
  if (debug) qcdhist.Print();

  RooHistPdf qcdTemplatePdf ("qcdTemplatePdf","",observable_list,qcdhist);
  if(debug) qcdTemplatePdf.Print();

  RooRealVar  mean_sig   ("mean_sig","",60.0, 45.0, 65.0);
  RooRealVar var_sig ("var_sig","" ,10.0, 1.0, 15.0); 
  RooRealVar  alpha_sig  ("alpha_sig","",-1,-20,-0.1);
  RooRealVar n_sig ("n_sig","",10,2,50);

  RooCBShape cb_sig ("cb_sig","",observable,mean_sig,var_sig,alpha_sig,n_sig);
  RooGaussian gaus_sig("gaus_sig","",observable,mean_sig,var_sig);

  RooPlot* frame = observable.frame(RooFit::Title("transverse mass [GeV]")) ;
  qcdhist.plotOn(frame) ;
  RooFitResult* rgaus = gaus_sig.fitTo(qcdhist,RooFit::SumW2Error(kTRUE),RooFit::Save(),RooFit::Range(minFitGaus,maxFitGaus));
  RooFitResult* rcb = cb_sig.fitTo(qcdhist,RooFit::SumW2Error(kTRUE),RooFit::Save(),RooFit::Range(minFitCB,maxFitCB));

  if (debug) {
    cout << "###### Printing fit result for gaussian" << endl;
    rgaus->Print();
    cout << "###### Printing fit result for crystal ball" << endl;
    rcb->Print();
  }

  gaus_sig.plotOn(frame,RooFit::LineColor(kRed));
  cb_sig.plotOn(frame,RooFit::LineColor(kBlue));

  // RooAbsReal* nll = NULL;
  // nll = cb_sig.createNLL(qcdhist,RooFit::Extended(kTRUE),RooFit::Verbose(-1));
  // if (debug) nll->Print();

  // RooMinimizer mfit(*nll);
  // if(not debug){
  //   mfit.setVerbose(kFALSE);
  //   mfit.setPrintLevel(-1);
  //   mfit.setVerbose(kFALSE);
  //   mfit.setPrintLevel(-1);
  // }

  // cout <<"######### Minimize " << endl;
  // mfit.minimize("Minuit2","minimize");
  // cout << "######### Minimize hesse" << endl;
  // mfit.minimize("Minuit2","hesse");
  // // cout << "######### Estimate minos errors for all parameters" << endl;
  // // mfit.minos(RooArgSet(signalNorm,backgroundNorm));
  // RooFitResult* fitResult = mfit.save("fitResult");
  // fitResult->Print();

  // TH1D* fit_hist = (TH1D*) qcdTemplatePdf.createHistogram("fit_hist",observable,RooFit::Binning(*binning));

  TCanvas* cToRemove = new TCanvas("cToRemove","",700,700);
  TH1D* htmp1 = new TH1D("htmp1","",1,0,1);
  cToRemove->cd();
  htmp1->Fill(0.5);
  htmp1->Draw("HE");
  CMS_lumi(cToRemove,"",false,false);
  setTDRStyle();
  cToRemove->SaveAs("/afs/cern.ch/user/m/mciprian/www/test_plot/test_roofit_toBeRemoved.png");
  system("rm /afs/cern.ch/user/m/mciprian/www/test_plot/test_roofit_toBeRemoved.png");
							   
  TCanvas* canvas = new TCanvas("canvas","",700,700);
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->SetRightMargin(0.06);

  // hist->Draw("HE");
  // fit_hist->Draw("HE SAME");
  frame->Draw();

  CMS_lumi(canvas,"",true,false);
  setTDRStyle();

  //  canvas->SaveAs("/afs/cern.ch/user/m/mciprian/www/test_plot/rayleigh.png");
  canvas->SaveAs("/afs/cern.ch/user/m/mciprian/www/test_plot/test_roofit_bkg.png");

  cout << "CB fit parameters" << endl;
  RooRealVar* mean_fitresult = (RooRealVar*) rcb->floatParsFinal().find("mean_sig");
  RooRealVar* sigma_fitresult = (RooRealVar*) rcb->floatParsFinal().find("var_sig");
  mean_fitresult->Print();
  sigma_fitresult->Print();


}



//============================================================                                                                                                              


void makeTH1Fit(const string& inputFileName_sigReg = "wmass_varhists.root",		  
		const Bool_t isMuon = false
		) {

  gStyle->SetOptFit();

  gROOT->SetBatch(kTRUE);
  TH1::SetDefaultSumw2(); //all the following histograms will automatically call TH1::Sumw2()                                        

  cout << endl;

  TH1D* hvar_sigReg = NULL;
  
  TFile* inputFile_sigReg = new TFile(inputFileName_sigReg.c_str(),"READ");
  if (!inputFile_sigReg || inputFile_sigReg->IsZombie()) {
    cout << "Error: signal file not opened. Exit" << endl;
    exit(EXIT_FAILURE);
  }
  
  inputFile_sigReg->cd();
  if (isMuon) hvar_sigReg   = (TH1D*) getHistCloneFromFile(inputFile_sigReg, "hmT", "qcd_mu");
  else        hvar_sigReg   = (TH1D*) getHistCloneFromFile(inputFile_sigReg, "hmT", "qcd_ele");
  checkNotNullPtr(hvar_sigReg,"hvar_sigReg");

  myRebinHisto(hvar_sigReg,3);
  hvar_sigReg->Scale(1./hvar_sigReg->Integral());
  makeRooFit(hvar_sigReg);

  cout << endl;

}
