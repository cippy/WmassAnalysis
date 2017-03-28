#include "utility.h"

using namespace std;

static UInt_t max_event_debug = 500000;
static Bool_t debug = false;

static UInt_t check_every_N = 100000;

static Double_t ptWmax = 30.0;
static vector<Double_t> wptBins = {0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 26.0, 30.0};

//================================================

class recoilResoResp {
public:
  recoilResoResp(const Double_t & ptMin, const Double_t & ptMax, TH1D* histo_perp, TH1D* histo_par, TH1D* histo_par_plusZpt):
    ptMin(ptMin),
    ptMax(ptMax),
    uperp(histo_perp),
    upar(histo_par),
    upar_plusZpt(histo_par_plusZpt)
  {
    ptMean = 0;
  }
  recoilResoResp():
    ptMin(0),
    ptMax(0),
    uperp(NULL),
    upar(NULL),
    upar_plusZpt(NULL){
    ptMean = 0;
  }

    ~recoilResoResp(){};

  Double_t ptMin;
  Double_t ptMax;
  Double_t ptMean;
  TH1D* uperp;
  TH1D* upar;
  TH1D* upar_plusZpt;

};


//================================================

void printSummary(const string& text = "", const vector<recoilResoResp> recoil = {}) {
  
  cout << endl;
  cout << "============================" << endl;
  cout << "SUMMARY: " << text << endl;
  cout << "============================" << endl;
  cout << "bin    pT(Z)    ";
  cout << "integral(uL)    mean(uL)    sigma(uL)    ";    
  cout << "integral(u||)    mean(u||)    sigma(u||)    ";
  cout << "integral(u||+pT(Z))    mean(u||+pT(Z))    sigma(u||+pT(Z))" << endl;

  for(UInt_t ibin = 0; ibin < recoil.size(); ibin++){

    cout << ibin << "    " << recoil.at(ibin).ptMean << "    ";
    cout << recoil.at(ibin).uperp->Integral() << "    ";
    cout << recoil.at(ibin).uperp->GetMean() << "    " << recoil.at(ibin).uperp->GetStdDev() << "    ";
    cout << recoil.at(ibin).upar->Integral() << "    ";
    cout << recoil.at(ibin).upar->GetMean() << "    " << recoil.at(ibin).upar->GetStdDev() << "    ";
    cout << recoil.at(ibin).upar_plusZpt->Integral() << "    ";
    cout << recoil.at(ibin).upar_plusZpt->GetMean() << "    " << recoil.at(ibin).upar_plusZpt->GetStdDev() << "    ";
    cout << endl;

  }

  cout << endl;

}

//================================================


void fillHistograms(vector<recoilResoResp> & recoilTkmet, 
		    vector<recoilResoResp> & recoilPfmet,
		    const string& inputDIR = "./",		    
		    const string& outputDIR = "./",
		    TFile* outputFile = NULL,
		    const Sample& sample = Sample::zjets,
		    const Bool_t isMuon = false
		    ) {


  gROOT->SetBatch(kTRUE);
  TH1::SetDefaultSumw2(); //all the following histograms will automatically call TH1::Sumw2()                                        

  cout << endl;
  cout << "================================================" << endl;
  cout << endl;

  TDirectory *dirSample = NULL;
  string sampleDir = getStringFromEnumSample(sample).c_str();
  cout << "Sample --> " << sampleDir << endl;
  if (outputFile->GetKey(sampleDir.c_str())) dirSample = outputFile->GetDirectory(sampleDir.c_str());
  else dirSample = outputFile->mkdir(sampleDir.c_str());
  dirSample->cd();

  // define electron ID at 8 TeV from --> https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaCutBasedIdentification   
  // use use tight trigger ID: for quantities not mentioned in this ID, use those from lepton loose ID                              
  // will use PF isolation, we don't have detector based ID in ntuples yet                                                   
  // electronID tightEleID_EB_8TeV(0.004, 0.03, 0.01, 0.12, 0.02, 0.1, 0.05);                                   
  // electronID tightEleID_EE_8TeV(0.005, 0.02, 0.03, 0.10, 0.02, 0.1, 0.05);                                                     
  electronID tightEleID_EB_8TeV(0.007, 0.15, 0.01, 0.12, 0.02, 0.2, 0.05, 1, 1);
  electronID tightEleID_EE_8TeV(0.009, 0.10, 0.03, 0.10, 0.02, 0.2, 0.05, 1, 1);
  electronID* eleID = NULL; // choose in loop for each event between EB or EE id    

  Int_t chargedLeptonFlavour = isMuon ? 13 : 11;
  Double_t lepEtaThreshold = isMuon ? 2.1 : 1.479;  // EB only for electron for now                                   
  Int_t lepTightIdThreshold = isMuon ? 1 : 3;
  Double_t looseIsoThreshold = isMuon ? 0.2 : 0.2;
  Double_t tightIsoThreshold = isMuon ? 0.12 : 0.15;
  
  cout << endl;

  TChain* chain = new TChain("treeProducerWMassEle");
  vector<Double_t> genwgtVec;
  buildChain8TeV(chain, genwgtVec, inputDIR, sample);

  TTreeReader reader (chain);

  TTreeReaderValue<Int_t> isData  (reader,"isData");

  // reco met
  TTreeReaderValue<Float_t> tkmet    (reader,"tkmet_pt");
  TTreeReaderValue<Float_t> tkmet_phi(reader,"tkmet_phi");
  TTreeReaderValue<Float_t> pfmet    (reader,"met_pt");
  TTreeReaderValue<Float_t> pfmet_phi(reader,"met_phi");

  // lepGood branch
  TTreeReaderValue<Int_t> nlep  (reader,"nLepGood");
  TTreeReaderArray<Int_t> lep_pdgId (reader,"LepGood_pdgId");
  TTreeReaderArray<Float_t> lep_pt (reader,"LepGood_pt");
  TTreeReaderArray<Float_t> lep_eta (reader,"LepGood_eta");
  TTreeReaderArray<Float_t> lep_phi (reader,"LepGood_phi");
  TTreeReaderArray<Float_t> lep_mass (reader,"LepGood_mass");
  TTreeReaderArray<Float_t> lep_relIso03 (reader,"LepGood_relIso03");
  TTreeReaderArray<Float_t> lep_relIso04 (reader,"LepGood_relIso04");
  TTreeReaderArray<Int_t> lep_tightId (reader,"LepGood_tightId");  // must make sure it is tuned on 8 TeV (it might be based on 13 TeV) 

  // for electronID
  TTreeReaderArray<Float_t> lep_detaIn (reader,"LepGood_detaIn");
  TTreeReaderArray<Float_t> lep_dphiIn (reader,"LepGood_dphiIn");
  TTreeReaderArray<Float_t> lep_sigmaIetaIeta (reader,"LepGood_sigmaIetaIeta");
  TTreeReaderArray<Float_t> lep_HoE (reader,"LepGood_hcalOverEcal");
  TTreeReaderArray<Float_t> lep_dxy (reader,"LepGood_dxy");
  TTreeReaderArray<Float_t> lep_dz (reader,"LepGood_dz");
  TTreeReaderArray<Int_t> lep_lostHits (reader,"LepGood_lostHits");
  TTreeReaderArray<Int_t> lep_convVeto (reader,"LepGood_convVetoFull");
  TTreeReaderArray<Float_t> lep_ecalEnergy (reader,"LepGood_ecalEnergy");
  TTreeReaderArray<Float_t> lep_eSuperClusterOverP (reader,"LepGood_eSuperClusterOverP");
  // 1/E - 1/p = 1/ecalEnergy - eSuperClusterOverP/ecalEnergy == (1 - eSuperClusterOverP)/ecalEnergy;

  // other electron related branches
  TTreeReaderArray<Float_t> lep_etaSc (reader,"LepGood_scEta");
  TTreeReaderArray<Float_t> lep_r9 (reader,"LepGood_r9");
  
  // MC reweight
  // must activate it only for non data sample
  ///////////////////////////////
  string dummybranch = "rho";
  if (sampleDir.find("data") == string::npos) dummybranch = "puWeight";  // PU weight. Can use this in tree because it was computed with all 8 TeV dataset
  TTreeReaderValue<Float_t> puw (reader,dummybranch.c_str());  // if running on data, create it as a branch existing in tree (it won't be used)
  dummybranch = "m2l";
  if (sampleDir.find("data") == string::npos) dummybranch = "LepEff_1lep";  
  TTreeReaderValue<Float_t> lepEfficiency (reader,dummybranch.c_str());

  // Match of lepton with trigger object 
  TTreeReaderValue<Int_t> HLT_SingleMu (reader,"HLT_SingleMu");
  TTreeReaderValue<Int_t> HLT_SingleEl (reader,"HLT_SingleEl");
  TTreeReaderArray<Float_t> lep_trgMatch (reader,"LepGood_trgMatch");

  long int nTotal = chain->GetEntries();
  long int nEvents = 0;
  
  Double_t wgt = 1.0;
  Double_t sumwgt = 0.0;

  TH1D *hbosonpt = new TH1D("hbosonpt","",30,0,30);
  TH1D *hlep1pt = new TH1D("hlep1pt","", 70, 0.0, 70.0);
  TH1D *hlep1eta = new TH1D("hlep1eta","", 48, -2.4, 2.4);
  TH1D *htkmet = new TH1D("htkmet","", 50, 0.0, 250);
  TH1D *hpfmet = new TH1D("hpfmet","", 50, 0.0, 250);
  TH1D *hlep1relIso03 = new TH1D("hlep1relIso03","", 45, 0.0, tightIsoThreshold);
  TH1D *hlep1relIso04 = new TH1D("hlep1relIso04","", 45, 0.0, tightIsoThreshold);

  TH1D* hchargLepPtMean = new TH1D("hchargLepPtMean","",wptBins.size()-1,wptBins.data());

  ////////////////////////////////////////////                                                         
  // to get correct weight depending on sample in chain                                        
  string currentFile = "";
  Int_t ifile = 0;
  ////////////////////            

  while(reader.Next()){

    cout.flush();
    if(nEvents % check_every_N == 0) cout<<"\r"<<"Analyzing events "<<double(nEvents)/nTotal*100<<" % ";
    nEvents++;
    if (debug && nEvents > max_event_debug) break; // to exit at some point, when debugging

    if(dynamic_cast<TChain*>(reader.GetTree())->GetFile()->GetName() != currentFile and currentFile != ""){
      currentFile = dynamic_cast<TChain*>(reader.GetTree())->GetFile()->GetName();
      ifile ++;
    } else if(dynamic_cast<TChain*>(reader.GetTree())->GetFile()->GetName() != currentFile and currentFile == ""){
      currentFile = dynamic_cast<TChain*>(reader.GetTree())->GetFile()->GetName();
    }

    // selection
    if (*nlep != 1) continue;    // 2 leptons
    if (fabs(lep_pdgId[0]) != chargedLeptonFlavour) continue;  
    if (lep_pt[0] < 30.0) continue;
    if (fabs(lep_eta[0]) > lepEtaThreshold) continue; 
    if (lep_trgMatch[0] < 0.5) continue;

    // tight ID on leading lepton
    if (isMuon) {

      if (*HLT_SingleMu == 0) continue;
      // ID + isolation
      if (lep_tightId[0] < lepTightIdThreshold) continue;   // tight ID
      if (lep_relIso04[0] > tightIsoThreshold) continue; 

    } else {

      if (*HLT_SingleEl == 0) continue;

      // ID + isolation
      if (fabs(lep_etaSc[0]) < 1.479) eleID = &tightEleID_EB_8TeV;
      else                            eleID = &tightEleID_EE_8TeV;

      // cutting on R9 produces large disagreement between data and MC. It could be because R9 is badly reproduced in MC
      // if (lep_r9[0] < 0.94) continue;
      // if (lep_r9[1] < 0.94) continue;

      // ID lep1
      if (lep_sigmaIetaIeta[0] > eleID->sigmaIetaIeta) continue;
      if (lep_HoE[0] > eleID->HoE) continue;
      if (lep_dz[0] > eleID->dz) continue;
      if ((1 - lep_eSuperClusterOverP[0])/lep_ecalEnergy[0] > eleID->invE_minus_invP) continue;
      if (lep_lostHits[0] > eleID->missingHits) continue;
      if (lep_convVeto[0] == 0) continue;
      if (lep_detaIn[0] > eleID->dEtaIn) continue;                                                                           
      if (lep_dphiIn[0] > eleID->dPhiIn) continue;                                                                             
      if (lep_dxy[0] > eleID->dxy) continue;                                             
      // iso 
      if (lep_relIso03[0] > tightIsoThreshold) continue; 

    }

    if (*isData == 1) wgt = 1.0;
    else wgt = intLumi * genwgtVec[ifile] * *puw * *lepEfficiency; 

    TLorentzVector lep1Reco;

    lep1Reco.SetPtEtaPhiM(lep_pt[0],lep_eta[0],lep_phi[0],lep_mass[0]); 
    Double_t chargLepPt = lep1Reco.Pt();

    // cout << " lep 1: " << lep1Reco.Pt() << " " << lep1Reco.Eta() << " " << lep1Reco.Phi() << " " << lep1Reco.M() << endl;
    // cout << " lep 2: " << lep2Reco.Pt() << " " << lep2Reco.Eta() << " " << lep2Reco.Phi() << " " << lep2Reco.M() << endl;
    // cout << " Z    : " << zReco.Pt() << " " << zReco.Eta() << " " << zReco.Phi() << " " << zReco.M() << endl;
    // cout << "mZ1 = " << *mZ1 << endl;


    TVector2 tkmetReco, pfmetReco, wReco2D, lep1Reco2D;
    tkmetReco.SetMagPhi(*tkmet,*tkmet_phi);
    pfmetReco.SetMagPhi(*pfmet,*pfmet_phi);
    // use PF MET to estimate W pt, since PF MET should be the "real" MET in the event
    lep1Reco2D = lep1Reco.Vect().XYvector();
    wReco2D = lep1Reco2D + pfmetReco;
    Double_t wRecoPt = wReco2D.Mod();
    if (wRecoPt > ptWmax) continue;
    TVector2 recoilRecoTk = -1 * (lep1Reco2D + tkmetReco);
    TVector2 recoilRecoPf = -1 * (lep1Reco2D + pfmetReco);

    UInt_t ptbin = 0;
    if (recoilTkmet.size() != 0) {
      
      for (UInt_t bin = 0; bin <= recoilTkmet.size()-1; bin++) {
	if (wRecoPt >= recoilTkmet.at(bin).ptMin and wRecoPt < recoilTkmet.at(bin).ptMax) ptbin = bin;
      }
      if (wRecoPt >  recoilTkmet.back().ptMax) ptbin = recoilTkmet.size()-1;

    }

    fillTH1(hbosonpt, wRecoPt, wgt);
    fillTH1(hlep1eta, lep1Reco.Eta(), wgt);
    fillTH1(hlep1pt, lep1Reco.Pt(), wgt);
    fillTH1(htkmet, *tkmet, wgt);
    fillTH1(hpfmet, *pfmet, wgt);
    fillTH1(hlep1relIso03, lep_relIso03[0], wgt);
    fillTH1(hlep1relIso04, lep_relIso04[0], wgt);

    /////////////////////
    // We compute recoil projection along the charged lepton (not along the W boson, as we would do for the Z)
    //////////////////////
    // tracker MET
    Double_t recoil_par = recoilRecoTk.Mod() * cos(recoilRecoTk.DeltaPhi(lep1Reco2D));
    Double_t recoil_perp = recoilRecoTk.Mod() * sin(recoilRecoTk.DeltaPhi(lep1Reco2D));
    Double_t recoil_par_plusZpt = recoil_par + zRecoPt;

    recoilTkmet.at(ptbin).upar_plusZpt->Fill(recoil_par_plusZpt,wgt);
    recoilTkmet.at(ptbin).uperp->Fill(recoil_perp,wgt);
    recoilTkmet.at(ptbin).upar->Fill(recoil_par,wgt);
    recoilTkmet.at(ptbin).ptMean += (zRecoPt * wgt);

    //////////////////////
    // PF MET
    recoil_par = recoilRecoPf.Mod() * cos(recoilRecoPf.DeltaPhi(zReco2D));
    recoil_perp = recoilRecoPf.Mod() * sin(recoilRecoPf.DeltaPhi(zReco2D));
    recoil_par_plusZpt = recoil_par + zRecoPt;

    recoilPfmet.at(ptbin).upar_plusZpt->Fill(recoil_par_plusZpt,wgt);
    recoilPfmet.at(ptbin).upar->Fill(recoil_par,wgt);
    recoilPfmet.at(ptbin).uperp->Fill(recoil_perp,wgt);
    recoilPfmet.at(ptbin).ptMean += (zRecoPt * wgt);

    ///////////////////

    sumwgt += wgt;

  }

  cout << endl;
  cout << "sum(wgt) = " << sumwgt << endl;

  Double_t mean = 0.0;

  for(UInt_t ibin = 0; ibin < recoilTkmet.size(); ibin++){
    if (recoilTkmet.at(ibin).uperp->Integral() > 0) {
      mean = recoilTkmet.at(ibin).ptMean / recoilTkmet.at(ibin).uperp->Integral(0,recoilTkmet.at(ibin).uperp->GetNbinsX()+1);
      recoilTkmet.at(ibin).ptMean = mean;
    } else {
      recoilTkmet.at(ibin).ptMean = (recoilTkmet.at(ibin).ptMax - recoilTkmet.at(ibin).ptMin) / 2.;
    }
    hzptMean->SetBinContent(ibin, recoilTkmet.at(ibin).ptMean);
    recoilTkmet.at(ibin).uperp->Write(0,TObject::kOverwrite);
    recoilTkmet.at(ibin).upar->Write(0,TObject::kOverwrite);
    recoilTkmet.at(ibin).upar_plusZpt->Write(0,TObject::kOverwrite);
  }

  for(UInt_t ibin = 0; ibin < recoilPfmet.size(); ibin++){
    if (recoilTkmet.at(ibin).uperp->Integral() > 0) {
      mean = recoilPfmet.at(ibin).ptMean / recoilPfmet.at(ibin).uperp->Integral(0,recoilPfmet.at(ibin).uperp->GetNbinsX()+1);
      recoilPfmet.at(ibin).ptMean = mean;
    } else {
      recoilPfmet.at(ibin).ptMean = (recoilPfmet.at(ibin).ptMax - recoilPfmet.at(ibin).ptMin) / 2.;
    }
    recoilPfmet.at(ibin).uperp->Write(0,TObject::kOverwrite);
    recoilPfmet.at(ibin).upar->Write(0,TObject::kOverwrite);
    recoilPfmet.at(ibin).upar_plusZpt->Write(0,TObject::kOverwrite);
  }

  cout << endl;
  cout << "Writing on output file" << endl;
  // if the file is opened in UPDATE mode, the following should overwrite an object if its key inside the file already exists
  // this needs to be tested                                     
  outputFile->Write(0,TObject::kOverwrite);

  // TCanvas* canvas = new TCanvas("canvas","",700,700);
  // canvas->cd();
  // hbosonpt->Draw("HE");
  // hbosonpt->GetXaxis()->SetTitle("Z p_{T} [GeV]");
  // hbosonpt->GetYaxis()->SetTitle("Events");
  // canvas->SaveAs((outputDIR + "ptZ_" + sampleDir + ".png" ).c_str());  
  // delete canvas;

  // canvas = new TCanvas("canvas","",700,700);
  // canvas->cd();
  // hzmass->Draw("HE");
  // hzmass->GetXaxis()->SetTitle("Z p_{T} [GeV]");
  // hzmass->GetYaxis()->SetTitle("Events");
  // canvas->SaveAs((outputDIR + "massZ_" + sampleDir + ".png" ).c_str());  
  // delete canvas;

  delete hbosonpt;
  delete hzmass;
  delete hlep1pt;
  delete hlep2pt;
  delete hlep1eta;
  delete hlep2eta;
  delete htkmet;
  delete hpfmet;
  delete hlep1relIso03;
  delete hlep1relIso04;
  delete hlep2relIso03;
  delete hlep2relIso04;

}

//================================================

void makeRecoilResoRespAna(const string& inputDIR = "./", const string& outputDIR = "./", const string& outfileName = "./",
			   const Bool_t isMuon = false
			   ) {

  createPlotDirAndCopyPhp(outputDIR);

  gROOT->SetBatch(kTRUE);
  TH1::SetDefaultSumw2(); //all the following histograms will automatically call TH1::Sumw2() 

  cout << endl;

  vector<recoilResoResp> recoilTkmet_zjets;
  vector<recoilResoResp> recoilPfmet_zjets;

  vector<recoilResoResp> recoilTkmet_data;
  vector<recoilResoResp> recoilPfmet_data;

  for(size_t ibin = 0; ibin < wptBins.size()-1; ibin++){

    // first TH1D is u_perp, second is u_par + ZpT
    
    recoilTkmet_zjets.push_back(recoilResoResp(wptBins.at(ibin),
					       wptBins.at(ibin+1),
					       new TH1D(Form("uperp_tkmet_zjets_pt_%d_%d",int(wptBins.at(ibin)),int(wptBins.at(ibin+1))),"",111,-50.5,60.5),
					       new TH1D(Form("upar_tkmet_zjets_pt_%d_%d",int(wptBins.at(ibin)),int(wptBins.at(ibin+1))),"",111,-50.5,60.5),
					       new TH1D(Form("uparPlusZpt_tkmet_zjets_pt_%d_%d",int(wptBins.at(ibin)),int(wptBins.at(ibin+1))),"",111,-50.5,60.5)
					       )
				); 
    
    recoilPfmet_zjets.push_back(recoilResoResp(wptBins.at(ibin),
					       wptBins.at(ibin+1),
					       new TH1D(Form("uperp_pfmet_zjets_pt_%d_%d",int(wptBins.at(ibin)),int(wptBins.at(ibin+1))),"",111,-50.5,60.5),
					       new TH1D(Form("upar_pfmet_zjets_pt_%d_%d",int(wptBins.at(ibin)),int(wptBins.at(ibin+1))),"",111,-50.5,60.5),
					       new TH1D(Form("uparPlusZpt_pfmet_zjets_pt_%d_%d",int(wptBins.at(ibin)),int(wptBins.at(ibin+1))),"",111,-50.5,60.5)
					       )
				); 

    recoilTkmet_data.push_back(recoilResoResp(wptBins.at(ibin),
					      wptBins.at(ibin+1),
					      new TH1D(Form("uperp_tkmet_data_pt_%d_%d",int(wptBins.at(ibin)),int(wptBins.at(ibin+1))),"",111,-50.5,60.5),
					      new TH1D(Form("upar_tkmet_data_pt_%d_%d",int(wptBins.at(ibin)),int(wptBins.at(ibin+1))),"",111,-50.5,60.5),
					      new TH1D(Form("uparPlusZpt_tkmet_data_pt_%d_%d",int(wptBins.at(ibin)),int(wptBins.at(ibin+1))),"",111,-50.5,60.5)
					      )
			       ); 
    
    recoilPfmet_data.push_back(recoilResoResp(wptBins.at(ibin),
					      wptBins.at(ibin+1),
					      new TH1D(Form("uperp_pfmet_data_pt_%d_%d",int(wptBins.at(ibin)),int(wptBins.at(ibin+1))),"",111,-50.5,60.5),
					      new TH1D(Form("upar_pfmet_data_pt_%d_%d",int(wptBins.at(ibin)),int(wptBins.at(ibin+1))),"",111,-50.5,60.5),
					      new TH1D(Form("uparPlusZpt_pfmet_data_pt_%d_%d",int(wptBins.at(ibin)),int(wptBins.at(ibin+1))),"",111,-50.5,60.5)
					      )
			       ); 

  }

  TFile* outputFile = new TFile((outputDIR + outfileName).c_str(),"UPDATE");
  if (!outputFile || outputFile->IsZombie()) {
    cout << "Error: file not opened. Exit" << endl;
    exit(EXIT_FAILURE);
  }

  outputFile->cd();

  fillHistograms(recoilTkmet_zjets, recoilPfmet_zjets, inputDIR, outputDIR, outputFile, Sample::zjets, isMuon); 
  if (isMuon)  fillHistograms(recoilTkmet_data, recoilPfmet_data, inputDIR, outputDIR, outputFile, Sample::data_doubleMu, isMuon); 
  else         fillHistograms(recoilTkmet_data, recoilPfmet_data, inputDIR, outputDIR, outputFile, Sample::data_doubleEG, isMuon); 

  outputFile->Close();

  printSummary("tracker MET --- MC",recoilTkmet_zjets);
  printSummary("PF MET --- MC",recoilPfmet_zjets);
  printSummary("tracker MET --- MC",recoilTkmet_data);
  printSummary("PF MET --- MC",recoilPfmet_data);


}


//================================================

void plotDistributions(const string& outputDIR = "./", 
		       const string& inputFileName_tmp = "./", 
		       const Bool_t isMuon = false
		       ) 
{

  string inputFileName = outputDIR + inputFileName_tmp;

  TFile* inputFile = new TFile(inputFileName.c_str(),"READ");
  if (!inputFile || inputFile->IsZombie()) {
    cout << "Error: file not opened. Exit" << endl;
    exit(EXIT_FAILURE);
  }

  inputFile->cd();

  string lepton = isMuon ? "muon" : "electron";

  vector<Int_t> rebinFactor;
  vector<string> hvarName;
  hvarName.push_back("hlep1eta");      rebinFactor.push_back(1);
  hvarName.push_back("hlep2eta");      rebinFactor.push_back(1);
  hvarName.push_back("hpfmet");      rebinFactor.push_back(2);
  hvarName.push_back("htkmet");      rebinFactor.push_back(2);
  hvarName.push_back("hlep1pt");      rebinFactor.push_back(1);
  hvarName.push_back("hlep2pt");      rebinFactor.push_back(1);
  hvarName.push_back("hzmass");      rebinFactor.push_back(1);
  hvarName.push_back("hbosonpt");      rebinFactor.push_back(1);
  if (isMuon) {
    hvarName.push_back("hlep1relIso04");      rebinFactor.push_back(1);
    hvarName.push_back("hlep2relIso04");      rebinFactor.push_back(1);
  } else {
    hvarName.push_back("hlep1relIso03");      rebinFactor.push_back(1);
    hvarName.push_back("hlep2relIso03");      rebinFactor.push_back(1);
  }

  vector<string> xAxisTitle;  // use "title::lowrange::uprange to pass lower and upper value for range                            
  xAxisTitle.push_back(("leading " + lepton + " #eta").c_str() );
  xAxisTitle.push_back(("trailing " + lepton + " #eta").c_str() );
  xAxisTitle.push_back("PF E_{T}^{miss} [GeV]");
  xAxisTitle.push_back("tracker E_{T}^{miss} [GeV]");
  xAxisTitle.push_back(("leading " + lepton + " p_{T} [GeV]").c_str());
  xAxisTitle.push_back(("trailing " + lepton + " p_{T} [GeV]").c_str());
  if (isMuon) {
    xAxisTitle.push_back("Z(#mu#mu) mass [GeV]");
    xAxisTitle.push_back("Z(#mu#mu) p_{T} [GeV]");
    xAxisTitle.push_back(("leading " + lepton + " isolation (#DeltaR = 0.4)").c_str());
    xAxisTitle.push_back(("trailing " + lepton + " isolation (#DeltaR = 0.4)").c_str());
  } else {
    xAxisTitle.push_back("Z(ee) mass [GeV]");
    xAxisTitle.push_back("Z(ee) p_{T} [GeV]");
    xAxisTitle.push_back(("leading " + lepton + " isolation (#DeltaR = 0.3)").c_str());
    xAxisTitle.push_back(("trailing " + lepton + " isolation (#DeltaR = 0.3)").c_str());
  }

  vector<string> canvasTitle;
  canvasTitle.push_back("etaLep1");
  canvasTitle.push_back("etaLep2");
  canvasTitle.push_back("pfmet");
  canvasTitle.push_back("tkmet");
  canvasTitle.push_back("ptLep1");
  canvasTitle.push_back("ptLep2");
  canvasTitle.push_back("Zmass");
  canvasTitle.push_back("Zpt");
  if (isMuon) {
    canvasTitle.push_back("lep1relIso04");
    canvasTitle.push_back("lep2relIso04");
  } else {
    canvasTitle.push_back("lep1relIso03");
    canvasTitle.push_back("lep2relIso03"); 
  }

  TH1D* hdata = NULL;
  TH1D* hzjets = NULL;

  // tmp plot to be removed to adjust settings in CMS_lumi                                                                                        
  TH1D* htmp1 = new TH1D("htmp1","",1,0,1);
  TH1D* htmp2 = new TH1D("htmp2","",1,0,1);
  htmp1->Fill(0.5);
  htmp2->Fill(0.5);
  vector<TH1*> htmpVec; htmpVec.push_back(htmp2);
  drawTH1dataMCstack(htmp1, htmpVec, "variable", "Events", "tmpToBeRemoved", outputDIR);
  system(("rm " + outputDIR + "*tmpToBeRemoved*").c_str());

  for (UInt_t i = 0; i < hvarName.size(); i++) {

    vector<TH1*> stackElementMC;  // first element is the one on top of the stack                   
    vector<string> stackLegendMC;
    
    if (isMuon) {
      hdata  = (TH1D*) getHistCloneFromFile(inputFile, hvarName[i], "data_doubleMu");
    } else {
      hdata  = (TH1D*) getHistCloneFromFile(inputFile, hvarName[i], "data_doubleEG");
    }
    hzjets = (TH1D*) getHistCloneFromFile(inputFile, hvarName[i], "zjets");

    checkNotNullPtr(hdata,"hdata");
    checkNotNullPtr(hzjets,"hzjets");

    stackElementMC.push_back(hzjets);

    if (isMuon) stackLegendMC.push_back("Z(#mu#mu)+jets");
    else stackLegendMC.push_back("Z(ee)+jets");

    drawTH1dataMCstack(hdata, stackElementMC, xAxisTitle[i], "Events", canvasTitle[i], outputDIR, "data", stackLegendMC, "data/MC", intLumi, rebinFactor[i]);

  }

  inputFile->Close();
  cout << endl;


}


//====================================================================

void fillTGraphAndDrawRecoil(vector<TH1D*>& hvec_recoil, 
			     TGraphErrors* graph = NULL,
			     TFile* inputFile = NULL, 
			     const string& outputDIR = "",
			     const Bool_t isMuon = false, 
			     const string& hvarName_base = "", 
			     const string& xAxisName = "",
			     const string& sampleLegEntry = "data",
			     const Bool_t saveCanvas = true
			     )

{
  
  for (UInt_t ptbin = 0; ptbin < wptBins.size()-1; ptbin++) {
    
    string hvarName = hvarName_base;
    hvarName += string(Form("_pt_%d_%d", (int) wptBins[ptbin], (int) wptBins[ptbin+1]));
    if (sampleLegEntry == "data") {
      if (isMuon) hvec_recoil.push_back((TH1D*) getHistCloneFromFile(inputFile, hvarName, "data_doubleMu"));
      else hvec_recoil.push_back((TH1D*) getHistCloneFromFile(inputFile, hvarName, "data_doubleEG"));
    } else {
      hvec_recoil.push_back((TH1D*) getHistCloneFromFile(inputFile, hvarName, sampleLegEntry.c_str()));
    }
    checkNotNullPtr(hvec_recoil.back(),Form("%s",hvarName.c_str()));
    
    TFitResultPtr fitRes;
    fitRes = drawTH1(hvec_recoil.back(),
		     xAxisName.c_str(),"Events", Form("%s_ptZ_%1.0f_%1.0f",hvarName_base.c_str(),wptBins[ptbin],wptBins[ptbin+1]),
		     outputDIR,
		     Form("%s::%1.0f < p_{T}^{Z} [GeV] < %1.0f",sampleLegEntry.c_str(),wptBins[ptbin],wptBins[ptbin+1]),
		     intLumi, 1, false, saveCanvas);
    
    // fitRes->Parameter(2); // CB mean                                                                                  
    // fitRes->ParError(2); // CB mean error                                                                                             
    // fitRes->Parameter(3); // CB sigma                                                                                  
    // fitRes->ParError(3); // CB sigma error                                                                                             
    
    graph->SetPoint(ptbin,(wptBins[ptbin+1]+wptBins[ptbin])/2.0,fitRes->Parameter(2));
    graph->SetPointError(ptbin,0.0,fitRes->ParError(2));
    
  }


}


//====================================================================

void plotRecoilResoResp(const string& outputDIR_tmp = "./", 
			const string& inputFileName_tmp = "./", 
			const Bool_t isMuon = false
			) 
{

  string outputDIR = outputDIR_tmp + "recoilPlots/";

  createPlotDirAndCopyPhp(outputDIR);

  gROOT->SetBatch(kTRUE);
  TH1::SetDefaultSumw2(); //all the following histograms will automatically call TH1::Sumw2() 

  string inputFileName = outputDIR_tmp + inputFileName_tmp;

  TFile* inputFile = new TFile(inputFileName.c_str(),"UPDATE");
  if (!inputFile || inputFile->IsZombie()) {
    cout << "Error: file not opened. Exit" << endl;
    exit(EXIT_FAILURE);
  }
  
  inputFile->cd();

  TH1D* htmp1 = new TH1D("htmp1","",1,0,1);
  htmp1->Fill(0.5);  
  drawTH1(htmp1,"x","y","tmpToBeRemoved",outputDIR,"leg");
  system(("rm " + outputDIR + "*tmpToBeRemoved*").c_str());
  delete htmp1;

  TH1D* hzptMean_data = NULL;
  if (isMuon) hzptMean_data = (TH1D*) getHistCloneFromFile(inputFile, "hzptMean", "data_doubleMu");
  else hzptMean_data = (TH1D*) getHistCloneFromFile(inputFile, "hzptMean", "data_doubleEG");
  checkNotNullPtr(hzptMean_data,"hzptMean_data");

  TH1D* hzptMean_zjets = NULL;
  hzptMean_zjets = (TH1D*) getHistCloneFromFile(inputFile, "hzptMean", "zjets");
  checkNotNullPtr(hzptMean_zjets,"hzptMean_zjets");

  string hvarName_base = "";
  //////////////////////////////////
  // data: PF MET
  ///////////////////////////////////
  hvarName_base = "uperp_pfmet_data";
  vector<TH1D *> hvec_uperp_pfmet_data;
  TGraphErrors *graph_uperp_pfmet_data = new TGraphErrors();

  fillTGraphAndDrawRecoil(hvec_uperp_pfmet_data, graph_uperp_pfmet_data, 
			  inputFile, outputDIR, isMuon, 
			  hvarName_base, "u_{#perp}  [GeV]", "data", false);

  hvarName_base = "uparPlusZpt_pfmet_data";
  vector<TH1D *> hvec_uparPlusZpt_pfmet_data;
  TGraphErrors *graph_uparPlusZpt_pfmet_data = new TGraphErrors();

  fillTGraphAndDrawRecoil(hvec_uparPlusZpt_pfmet_data, graph_uparPlusZpt_pfmet_data, 
			  inputFile, outputDIR, isMuon, 
			  hvarName_base, "u_{||} + p_{T}^{Z}  [GeV]", "data", true);

  hvarName_base = "upar_pfmet_data";
  vector<TH1D *> hvec_upar_pfmet_data;
  TGraphErrors *graph_upar_pfmet_data = new TGraphErrors();

  fillTGraphAndDrawRecoil(hvec_upar_pfmet_data, graph_upar_pfmet_data, 
			  inputFile, outputDIR, isMuon, 
			  hvarName_base, "u_{||}  [GeV]", "data", false);
  //////////////////////////////////////////

  //////////////////////////////////
  // data: TK MET
  ///////////////////////////////////
  hvarName_base = "uperp_tkmet_data";
  vector<TH1D *> hvec_uperp_tkmet_data;
  TGraphErrors *graph_uperp_tkmet_data = new TGraphErrors();

  fillTGraphAndDrawRecoil(hvec_uperp_tkmet_data, graph_uperp_tkmet_data, 
			  inputFile, outputDIR, isMuon, 
			  hvarName_base, "u_{#perp}  [GeV]", "data", false);

  hvarName_base = "uparPlusZpt_tkmet_data";
  vector<TH1D *> hvec_uparPlusZpt_tkmet_data;
  TGraphErrors *graph_uparPlusZpt_tkmet_data = new TGraphErrors();

  fillTGraphAndDrawRecoil(hvec_uparPlusZpt_tkmet_data, graph_uparPlusZpt_tkmet_data, 
			  inputFile, outputDIR, isMuon, 
			  hvarName_base, "u_{||} + p_{T}^{Z}  [GeV]", "data", true);

  hvarName_base = "upar_tkmet_data";
  vector<TH1D *> hvec_upar_tkmet_data;
  TGraphErrors *graph_upar_tkmet_data = new TGraphErrors();

  fillTGraphAndDrawRecoil(hvec_upar_tkmet_data, graph_upar_tkmet_data, 
			  inputFile, outputDIR, isMuon, 
			  hvarName_base, "u_{||}  [GeV]", "data", false);
  //////////////////////////////////////////

  //////////////////////////////////
  // zjets: PF MET
  ///////////////////////////////////
  hvarName_base = "uperp_pfmet_zjets";
  vector<TH1D *> hvec_uperp_pfmet_zjets;
  TGraphErrors *graph_uperp_pfmet_zjets = new TGraphErrors();

  fillTGraphAndDrawRecoil(hvec_uperp_pfmet_zjets, graph_uperp_pfmet_zjets, 
			  inputFile, outputDIR, isMuon, 
			  hvarName_base, "u_{#perp}  [GeV]", "zjets", false);

  hvarName_base = "uparPlusZpt_pfmet_zjets";
  vector<TH1D *> hvec_uparPlusZpt_pfmet_zjets;
  TGraphErrors *graph_uparPlusZpt_pfmet_zjets = new TGraphErrors();

  fillTGraphAndDrawRecoil(hvec_uparPlusZpt_pfmet_zjets, graph_uparPlusZpt_pfmet_zjets, 
			  inputFile, outputDIR, isMuon, 
			  hvarName_base, "u_{||} + p_{T}^{Z}  [GeV]", "zjets", true);

  hvarName_base = "upar_pfmet_zjets";
  vector<TH1D *> hvec_upar_pfmet_zjets;
  TGraphErrors *graph_upar_pfmet_zjets = new TGraphErrors();

  fillTGraphAndDrawRecoil(hvec_upar_pfmet_zjets, graph_upar_pfmet_zjets, 
			  inputFile, outputDIR, isMuon, 
			  hvarName_base, "u_{||}  [GeV]", "zjets", false);
  //////////////////////////////////////////

  //////////////////////////////////
  // zjets: TK MET
  ///////////////////////////////////
  hvarName_base = "uperp_tkmet_zjets";
  vector<TH1D *> hvec_uperp_tkmet_zjets;
  TGraphErrors *graph_uperp_tkmet_zjets = new TGraphErrors();

  fillTGraphAndDrawRecoil(hvec_uperp_tkmet_zjets, graph_uperp_tkmet_zjets, 
			  inputFile, outputDIR, isMuon, 
			  hvarName_base, "u_{#perp}  [GeV]", "zjets", false);

  hvarName_base = "uparPlusZpt_tkmet_zjets";
  vector<TH1D *> hvec_uparPlusZpt_tkmet_zjets;
  TGraphErrors *graph_uparPlusZpt_tkmet_zjets = new TGraphErrors();

  fillTGraphAndDrawRecoil(hvec_uparPlusZpt_tkmet_zjets, graph_uparPlusZpt_tkmet_zjets, 
			  inputFile, outputDIR, isMuon, 
			  hvarName_base, "u_{||} + p_{T}^{Z}  [GeV]", "zjets", true);

  hvarName_base = "upar_tkmet_zjets";
  vector<TH1D *> hvec_upar_tkmet_zjets;
  TGraphErrors *graph_upar_tkmet_zjets = new TGraphErrors();

  fillTGraphAndDrawRecoil(hvec_upar_tkmet_zjets, graph_upar_tkmet_zjets, 
			  inputFile, outputDIR, isMuon, 
			  hvarName_base, "u_{||}  [GeV]", "zjets", false);
  //////////////////////////////////////////

  TCanvas* canvas2 = new TCanvas("canvas2","",600,650);
  canvas2->cd();
  canvas2->SetTickx(1);
  canvas2->SetTicky(1);
  canvas2->cd();
  canvas2->SetRightMargin(0.06);

  graph_uparPlusZpt_pfmet_data->SetLineColor(kBlack);
  graph_uparPlusZpt_pfmet_data->SetMarkerColor(kBlack);
  graph_uparPlusZpt_pfmet_data->SetMarkerStyle(20);
  graph_uparPlusZpt_pfmet_data->SetMarkerSize(1);
  graph_uparPlusZpt_pfmet_data->GetXaxis()->SetTitle("Z p_{T} [GeV]");
  graph_uparPlusZpt_pfmet_data->GetYaxis()->SetTitle("u_{||} + p_{T}^{Z}  [GeV]");
  graph_uparPlusZpt_pfmet_data->GetYaxis()->SetRangeUser(-10.0,30.0);
  graph_uparPlusZpt_pfmet_data->Draw("AP");

  graph_uparPlusZpt_pfmet_zjets->SetLineColor(kRed);
  graph_uparPlusZpt_pfmet_zjets->SetMarkerColor(kRed);
  graph_uparPlusZpt_pfmet_zjets->SetMarkerStyle(24);
  graph_uparPlusZpt_pfmet_zjets->SetMarkerSize(1);
  graph_uparPlusZpt_pfmet_zjets->Draw("Psame");

  graph_uparPlusZpt_tkmet_data->SetLineColor(kBlue);
  graph_uparPlusZpt_tkmet_data->SetMarkerColor(kBlue);
  graph_uparPlusZpt_tkmet_data->SetMarkerStyle(24);
  graph_uparPlusZpt_tkmet_data->SetMarkerSize(1);
  graph_uparPlusZpt_tkmet_data->Draw("Psame");

  graph_uparPlusZpt_tkmet_zjets->SetLineColor(kGreen+2);
  graph_uparPlusZpt_tkmet_zjets->SetMarkerColor(kGreen+2);
  graph_uparPlusZpt_tkmet_zjets->SetMarkerStyle(22);
  graph_uparPlusZpt_tkmet_zjets->SetMarkerSize(1);
  graph_uparPlusZpt_tkmet_zjets->Draw("Psame");

  TLegend leg (0.15,0.6,0.53,0.9);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  leg.AddEntry(graph_uparPlusZpt_pfmet_data,"data: PF MET","PLE");
  leg.AddEntry(graph_uparPlusZpt_tkmet_data,"data: track MET","PLE");
  if (isMuon) {
    leg.AddEntry(graph_uparPlusZpt_pfmet_zjets,"Z(#mu#mu)+jets: PF MET","PLE");
    leg.AddEntry(graph_uparPlusZpt_tkmet_zjets,"Z(#mu#mu)+jets: track MET","PLE");
  } else {
    leg.AddEntry(graph_uparPlusZpt_pfmet_zjets,"Z(ee)+jets: PF MET","PLE");
    leg.AddEntry(graph_uparPlusZpt_tkmet_zjets,"Z(ee)+jets: track MET","PLE");
  }
  leg.Draw("same");
  CMS_lumi(canvas2,Form("%.1f",intLumi),true,false);

  canvas2->SaveAs((outputDIR+"graph_uparPlusZpt_vs_Zpt.png").c_str(),"png");
  canvas2->SaveAs((outputDIR+"graph_uparPlusZpt_vs_Zpt.pdf").c_str(),"pdf");

  // graph->Write(("graph_" + hvarName_base).c_str(), TObject::kOverwrite);
  // TCanvas *c = new TCanvas("c","",700,700);
  // graph->SetLineColor(kBlack);
  // graph->SetMarkerColor(kBlack);
  // graph->SetMarkerStyle(20);
  // graph->SetMarkerSize(1);
  // graph->Draw("AP");
  // c->SaveAs("/afs/cern.ch/user/m/mciprian/www/test_plot/graph.png");
  // delete c;
    
  inputFile->Close();
  delete inputFile;

}

//====================================================================



void makeRecoilResoResp(const string& inputDIR = "./", const string& outputDIR_tmp = "./", const string& outfileName = "wmass_resoresphists.root",
			const Bool_t isMuon = false, 
			const UInt_t skip_no0_sel1_plot2 = 0
			) {

  string outputDIR = outputDIR_tmp;
  outputDIR += isMuon ? "Zmumu/" : "Zee/";

  createPlotDirAndCopyPhp(outputDIR);

  gROOT->SetBatch(kTRUE);
  TH1::SetDefaultSumw2(); //all the following histograms will automatically call TH1::Sumw2() 

  cout << endl;

  if (skip_no0_sel1_plot2 != 1) makeRecoilResoRespAna(inputDIR, outputDIR, outfileName, isMuon);
  if (skip_no0_sel1_plot2 != 2) {
    //plotDistributions(outputDIR, outfileName, isMuon);
    plotRecoilResoResp(outputDIR, outfileName, isMuon);
  }

}
