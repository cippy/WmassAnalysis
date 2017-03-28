#include "utility.h"

using namespace std;

static Int_t nZmassBins = 61;
static Double_t zMassMin = 59.5;
static Double_t zMassMax = 120.5;

static Int_t nMtBins = 120;
static Double_t mtMin = 20;
static Double_t mtMax = 140;

static Double_t DR_match = 0.1;

static Bool_t useTrackMet = true;

//=============================================================


void fillRecoGenHistograms(const string& inputDIR = "./", const string& outputDIR = "./", const Sample& sample = Sample::zjets, const Bool_t isWregion = false) {

  gROOT->SetBatch(kTRUE);
  TH1::SetDefaultSumw2(); //all the following histograms will automatically call TH1::Sumw2()  

  cout << endl;

  TChain*chain = new TChain("tree");
  TChain* chFriend = new TChain("mjvars/t");

  //string treePath = treeLocation;
  //string sampleName = "WJetsToLNu";
  //string sampleName = "DYJetsToLL";
  //buildChain(chain, treePath, sampleName); 
  buildChain(chain, chFriend, inputDIR, sample); 

  TTreeReader reader (chain);

  //TTreeReaderValue<Double_t> evtHasGoodVtx(reader,"evtHasGoodVtx");

  // reco met
  TTreeReaderValue<Float_t> tkmet    (reader,"met_trkPt");
  TTreeReaderValue<Float_t> tkmet_phi(reader,"met_trkPhi");
  TTreeReaderValue<Float_t> pfmet    (reader,"met_pt");
  TTreeReaderValue<Float_t> pfmet_phi(reader,"met_phi");
  //TTreeReaderValue<Float_t> pfmet_mass(reader,"met_mass");

  // gen met
  TTreeReaderValue<Float_t> genpfmet    (reader,"met_genPt");
  TTreeReaderValue<Float_t> genpfmet_phi(reader,"met_genPhi");

  // lepGood branch
  TTreeReaderValue<Int_t> nlep  (reader,"nLepGood");
  TTreeReaderArray<Int_t> lep_pdgId (reader,"LepGood_pdgId");
  TTreeReaderArray<Float_t> lep_pt (reader,"LepGood_pt");
  TTreeReaderArray<Float_t> lep_eta (reader,"LepGood_eta");
  TTreeReaderArray<Float_t> lep_phi (reader,"LepGood_phi");
  TTreeReaderArray<Float_t> lep_mass (reader,"LepGood_mass");
  TTreeReaderArray<Float_t> lep_relIso03 (reader,"LepGood_relIso03");
  TTreeReaderArray<Int_t> lep_tightId (reader,"LepGood_tightId");
  TTreeReaderArray<Float_t> lep_etaSc (reader,"LepGood_etaSc");
  TTreeReaderArray<Float_t> lep_r9 (reader,"LepGood_full5x5_r9");
  TTreeReaderArray<Float_t> lep_sigmaIetaIeta (reader,"LepGood_full5x5_sigmaIetaIeta");
  TTreeReaderArray<Int_t> lep_hltSafeID (reader,"LepGood_hltId");

  // mZ
  TTreeReaderValue<Float_t> mZ1 (reader,"mZ1");

  // gen particles
  TTreeReaderValue<Int_t> nGenPart(reader,"nGenPart");
  TTreeReaderArray<Int_t> GenPart_pdgId(reader,"GenPart_pdgId");
  TTreeReaderArray<Int_t> GenPart_motherId(reader,"GenPart_motherId");
  TTreeReaderArray<Int_t> GenPart_motherIndex(reader,"GenPart_motherIndex");
  TTreeReaderArray<Float_t> GenPart_pt(reader,"GenPart_pt");
  TTreeReaderArray<Float_t> GenPart_eta(reader,"GenPart_eta");
  TTreeReaderArray<Float_t> GenPart_phi(reader,"GenPart_phi");
  TTreeReaderArray<Float_t> GenPart_mass(reader,"GenPart_mass");

  // MC reweight
  TTreeReaderValue<Float_t> weight(reader,"weight");
  TTreeReaderValue<Float_t> puw(reader,"puw");

  Int_t WZ_firstGenPartIndex = -1; // used to find gen particles
  Int_t WZ_secondGenPartIndex = -1; // used to find gen particles
  Int_t WZ_bosonGenPartIndex = -1;

  // TH1
  TH1D* hmT = new TH1D("hmT","",nMtBins, mtMin, mtMax);
  TH1D* hmZ = new TH1D("hmZ","",nZmassBins, zMassMin, zMassMax);   // 61, 59.5, 120.5   31, 74.5, 105.5
  TH1D* hpfmet = new TH1D("hpfmet","",30,0,300);
  TH1D* hlep1pt = new TH1D("hlep1pt","",28,24,80);  //60, 0, 300
  TH1D* hlep2pt = new TH1D("hlep2pt","",35,10,80);
  TH1D* hbosonpt = new TH1D("hbosonpt","",30,0,30);
  TH1D* hbosoneta = new TH1D("hbosoneta","",100,-5,5);

  TH1D* hgen_mT = new TH1D("hgen_mT","",nMtBins, mtMin, mtMax);
  TH1D* hgen_mZ = new TH1D("hgen_mZ","",nZmassBins, zMassMin, zMassMax);
  TH1D* hgen_pfmet = new TH1D("hgen_pfmet","",30,0,300);
  TH1D* hgen_lep1pt = new TH1D("hgen_lep1pt","",28,24,80);
  TH1D* hgen_lep2pt = new TH1D("hgen_lep2pt","",35,10,80);
  TH1D* hgen_bosonpt = new TH1D("hgen_bosonpt","",30,0,30);
  TH1D* hgen_bosoneta = new TH1D("hgen_bosoneta","",100,-5,5);

  TH2D* h2_mT_mZ = new TH2D("h2_mT_mZ","",nMtBins, mtMin, mtMax, nZmassBins, zMassMin, zMassMax);
  TH2D* h2_mT_lep1pt = new TH2D("h2_mT_lep1pt","",nMtBins, mtMin, mtMax, 28,24,80);
  TH2D* h2_mT_lep1eta = new TH2D("h2_mT_lep1eta","",nMtBins, mtMin, mtMax, 30, -1.5, 1.5);
  TH2D* h2_mT_lep1sigIetaIeta = new TH2D("h2_mT_lep1sigIetaIeta","",nMtBins, mtMin, mtMax, 25, 0.0, 0.025);
  TH2D* h2_mT_lep1r9 = new TH2D("h2_mT_lep1r9","",nMtBins, mtMin, mtMax, 55, 0.0, 1.1);
  TH2D* h2_mT_pfmet = new TH2D("h2_mT_pfmet","",nMtBins, mtMin, mtMax, 30,0,300);
  TH2D* h2_mT_tkmet = new TH2D("h2_mT_tkmet","",nMtBins, mtMin, mtMax, 30,0,300);
  TH2D* h2_mT_bosonPt = new TH2D("h2_mT_bosonPt","",nMtBins, mtMin, mtMax, 30,0,30);
  TH2D* h2_mT_lep1relIso03 = new TH2D("h2_mT_lep1relIso03","",nMtBins, mtMin, mtMax, 20,0,0.060);

  TH2D* h2_lep1pt_mZ = new TH2D("h2_lep1pt_mZ","",28,24,80, nZmassBins, zMassMin, zMassMax); 
  TH2D* h2_lep1pt_lep1eta = new TH2D("h2_lep1pt_lep1eta","",28,24,80, 30, -1.5, 1.5);
  TH2D* h2_lep1pt_lep1sigIetaIeta = new TH2D("h2_lep1pt_lep1sigIetaIeta","",28,24,80, 25, 0.0, 0.025);
  TH2D* h2_lep1pt_lep1r9 = new TH2D("h2_lep1pt_lep1r9","",28,24,80, 55, 0.0, 1.1);
  TH2D* h2_lep1pt_pfmet = new TH2D("h2_lep1pt_pfmet","",28,24,80, 30,0,300);
  TH2D* h2_lep1pt_tkmet = new TH2D("h2_lep1pt_tkmet","",28,24,80, 30,0,300);
  TH2D* h2_lep1pt_bosonPt = new TH2D("h2_lep1pt_bosonPt","",28,24,80, 30,0,30);
  TH2D* h2_lep1pt_lep1relIso03 = new TH2D("h2_lep1pt_lep1relIso03","",28,24,80, 20,0,0.060);

  // start event loop                                                                                                                                                       
  long int nTotal = chain->GetEntries();
  long int nEvents = 0;
  
  Double_t wgt = 1.0;
  Double_t evtCounter = 0.0;

  Double_t metToUse_pt = 0.0;
  Double_t metToUse_phi = 0.0;

  while(reader.Next()){
  
    cout.flush();
    if(nEvents % 100000 == 0) cout<<"\r"<<"Analyzing events "<<double(nEvents)/nTotal*100<<" % ";
    nEvents++;

    if (useTrackMet) {
      metToUse_pt = *pfmet;
      metToUse_phi = *pfmet_phi;
    } else {
      metToUse_pt = *tkmet;
      metToUse_phi = *tfmet_phi;
    }

    if (isWregion) {

      // selection
      if (*nlep != 1) continue;    // 1 leptons                                                  
      if (fabs(lep_pdgId[0]) != 11 ) continue;  // electrons                                                     
      //if ( lep_pt[0] < 30.0) continue;
      //if ( lep_tightId[0] < 3 ) continue;   // tight ID
      if ( fabs(lep_eta[0]) > 1.479 ) continue;  // EB only for now
      //if ( !(fabs(lep_pdgId[0]) == 11 && lep_hltSafeID[0] == 1) ) continue;  //HLT safe ID for electrons
      //if ( metToUse_pt < 30) continue;
      //if ( lep_relIso03[0] > (fabs(lep_etaSc[0]) < 1.479 ? 0.0588 : 0.0571)) continue;  // tight iso first lep

    } else {

      // selection
      if (*nlep != 2) continue;    // 2 leptons                                                  
      if (fabs(lep_pdgId[0]) != 11 || fabs(lep_pdgId[1]) != 11) continue;  // electrons                                                     
      if ( (lep_pdgId[0] + lep_pdgId[1]) != 0) continue;   // opposite sign                                              
      if ( *mZ1 < zMassMin || *mZ1 > zMassMax) continue;
      //if ( lep_pt[0] < 30.0) continue;
      //if ( lep_tightId[0] < 3 ) continue;   // tight ID
      //if ( lep_pt[1] < 10.0 ) continue;
      if ( fabs(lep_eta[0]) > 1.479 || fabs(lep_eta[1]) > 1.479 ) continue;  // EB only for now
      // if ( !( (fabs(lep_pdgId[0]) == 11 && lep_hltSafeID[0] == 1) && (fabs(lep_pdgId[1]) == 11 && lep_hltSafeID[1] == 1) ) ) continue;  //HLT safe ID for electrons
      // if ( lep_relIso03[0] > (fabs(lep_etaSc[0]) < 1.479 ? 0.0588 : 0.0571)) continue;  // tight iso first lep
      // if ( lep_relIso03[1] > 0.175 ) continue;                          // loose iso second lep

    }
    
    wgt = *puw * *weight;

    TLorentzVector lep1Reco, lep2Reco, bosonReco, metReco, metWlikeReco;
    lep1Reco.SetPtEtaPhiM(lep_pt[0],lep_eta[0],lep_phi[0],lep_mass[0]); 
    metReco.SetPtEtaPhiM(metToUse_pt,0,metToUse_phi,0);
    metWlikeReco = metReco;

    if (isWregion) lep2Reco = metWlikeReco;  // this would be the neutrino for the W phase space, btw pfmet_mass is always 0 or < 10^-6 
    else {
      lep2Reco.SetPtEtaPhiM(lep_pt[1],lep_eta[1],lep_phi[1],lep_mass[1]);
      metWlikeReco += lep2Reco;
    }

    bosonReco = lep1Reco + lep2Reco;
    if (bosonReco.Pt() > 30) continue;

    Double_t mT = sqrt(2. * lep1Reco.Pt() * metWlikeReco.Pt() * (1. - cos(lep1Reco.DeltaPhi(metWlikeReco)) ) );
    if (mT < mtMin) continue;
    //if (mT > mtMax) continue;

    TVector2 recoilReco;
    recoilReco.SetMagPhi((Double_t) *pfmet, (Double_t) *pfmet_phi);
    if (isWregion) {
      recoilReco += lep1Reco.Vect().XYvector();
      if (recoilReco.Mod() > 30) continue;
    } else {
      recoilReco += bosonReco.Vect().XYvector();
      if (recoilReco.Mod() > 15) continue;
    }

    // for W, the first gen lep returned is the charged lepton
    if (isWregion) 
      WZ_firstGenPartIndex = findGenPartFromWLNu(*nGenPart, GenPart_pdgId, GenPart_motherId, LepFlavour::electron, WZ_secondGenPartIndex, WZ_bosonGenPartIndex);
    else 
      WZ_firstGenPartIndex = findGenPartFromZLL(*nGenPart, GenPart_pdgId, GenPart_motherId, LepFlavour::electron, WZ_secondGenPartIndex, WZ_bosonGenPartIndex);
	
    if (WZ_firstGenPartIndex < 0 || WZ_secondGenPartIndex < 0) continue;

    TLorentzVector lep1Gen, lep2Gen;
    // here lep1 and lep2 are ordered with respect to the position in the genPart array, which is not ordered by pt
    // in the following, we define lep1gen as the one matching lep1reco, using DR < 0.3 for the matching
    lep1Gen.SetPtEtaPhiM(GenPart_pt[WZ_firstGenPartIndex],GenPart_eta[WZ_firstGenPartIndex],GenPart_phi[WZ_firstGenPartIndex],GenPart_mass[WZ_firstGenPartIndex]);
    lep2Gen.SetPtEtaPhiM(GenPart_pt[WZ_secondGenPartIndex],GenPart_eta[WZ_secondGenPartIndex],GenPart_phi[WZ_secondGenPartIndex],GenPart_mass[WZ_secondGenPartIndex]);
    TLorentzVector bosonGen;
    bosonGen.SetPtEtaPhiM(GenPart_pt[WZ_bosonGenPartIndex],GenPart_eta[WZ_bosonGenPartIndex],GenPart_phi[WZ_bosonGenPartIndex],GenPart_mass[WZ_bosonGenPartIndex]);

    // if lep2Gen is closer to lep1Reco than lep1Gen, switch lep1Gen and lep2Gen (so that lep1 and lep2 naming scheme is consistent)
    // switch WZ_secondGenPartIndex and WZ_firstGenPartIndex as well
    // this is only for Z, because for W the first lepton is the charged one, so we only require a DR match for it

    if (isWregion) {

      if (lep1Gen.DrEtaPhi(lep1Reco) > DR_match) continue;

    } else {

      if ( lep1Gen.DrEtaPhi(lep1Reco) > lep2Gen.DrEtaPhi(lep1Reco)) {
    	lep2Gen.SetPtEtaPhiM(GenPart_pt[WZ_firstGenPartIndex],GenPart_eta[WZ_firstGenPartIndex],GenPart_phi[WZ_firstGenPartIndex],GenPart_mass[WZ_firstGenPartIndex]);
    	lep1Gen.SetPtEtaPhiM(GenPart_pt[WZ_secondGenPartIndex],GenPart_eta[WZ_secondGenPartIndex],GenPart_phi[WZ_secondGenPartIndex],GenPart_mass[WZ_secondGenPartIndex]);
    	Int_t tmpIndex = WZ_secondGenPartIndex;
    	WZ_secondGenPartIndex = WZ_firstGenPartIndex;
    	WZ_firstGenPartIndex = tmpIndex;
      } 
      if (lep1Gen.DrEtaPhi(lep1Reco) > DR_match || lep2Gen.DrEtaPhi(lep2Reco) > DR_match) continue;

    }

    // reco TH1
    fillTH1(hmT,(Double_t) mT, wgt);
    fillTH1(hmZ,(Double_t) *mZ1, wgt);
    fillTH1(hpfmet,(Double_t) *pfmet, wgt);
    fillTH1(hlep1pt,(Double_t) lep1Reco.Pt(), wgt);
    fillTH1(hlep2pt,(Double_t) lep2Reco.Pt(), wgt);
    fillTH1(hbosonpt,(Double_t) bosonReco.Pt(), wgt);
    fillTH1(hbosoneta,(Double_t) bosonReco.Eta(), wgt);

    // reco TH2
    fillTH2(h2_mT_mZ, mT, *mZ1, wgt);
    fillTH2(h2_mT_lep1pt, mT, lep_pt[0], wgt);
    fillTH2(h2_mT_lep1eta, mT, lep_eta[0], wgt);
    fillTH2(h2_mT_lep1sigIetaIeta, mT, lep_sigmaIetaIeta[0], wgt);
    fillTH2(h2_mT_lep1r9, mT, lep_r9[0], wgt);
    fillTH2(h2_mT_pfmet, mT, *pfmet, wgt);
    fillTH2(h2_mT_tkmet, mT, *tkmet, wgt);
    fillTH2(h2_mT_bosonPt, mT, bosonReco.Pt(), wgt);
    fillTH2(h2_mT_lep1relIso03, mT, lep_relIso03[0], wgt);

    fillTH2(h2_lep1pt_mZ, lep_pt[0], *mZ1, wgt);
    fillTH2(h2_lep1pt_lep1eta, lep_pt[0], lep_eta[0], wgt);
    fillTH2(h2_lep1pt_lep1sigIetaIeta, lep_pt[0], lep_sigmaIetaIeta[0], wgt);
    fillTH2(h2_lep1pt_lep1r9, lep_pt[0], lep_r9[0], wgt);
    fillTH2(h2_lep1pt_pfmet, lep_pt[0], *pfmet, wgt);
    fillTH2(h2_lep1pt_tkmet, lep_pt[0], *tkmet, wgt);
    fillTH2(h2_lep1pt_bosonPt, lep_pt[0], bosonReco.Pt(), wgt);
    fillTH2(h2_lep1pt_lep1relIso03, lep_pt[0], lep_relIso03[0], wgt);

    mT = sqrt(2. * lep1Gen.Pt() * lep2Gen.Pt() * (1. - cos(lep1Gen.DeltaPhi(lep2Gen)) ) );

    fillTH1(hgen_mT,(Double_t) mT, wgt);
    fillTH1(hgen_mZ,(Double_t) bosonGen.M(), wgt);
    fillTH1(hgen_pfmet,(Double_t) *genpfmet, wgt);
    fillTH1(hgen_lep1pt,(Double_t) lep1Gen.Pt(), wgt);
    fillTH1(hgen_lep2pt,(Double_t) lep2Gen.Pt(), wgt);
    fillTH1(hgen_bosonpt,(Double_t) bosonGen.Pt(), wgt);
    fillTH1(hgen_bosoneta,(Double_t) bosonGen.Eta(), wgt);

    evtCounter += wgt;

  }

  // drawPlot1D((TH1D*)hmZ->Clone(), hmZ, "invariant mass(e^{+}e^{-}) [GeV]", "massZee", outputDIR, "Events");
  // drawPlot1D((TH1D*)hmT->Clone(), hmT, "transverse mass(e^{+}e^{-}) [GeV]", "transverseMassZee", outputDIR, "Events");

  cout << endl;
  cout << "evtCounter = " << evtCounter << endl;

  if (isWregion) {
    
    drawTH1pair(hbosoneta, hgen_bosoneta, "boson #eta", "Events", "bosonEtaRecoGen", outputDIR, "Reco W(l#nu)+jets", "Gen W(l#nu)+jets", "Reco/Gen");
    drawTH1pair(hmT, hgen_mT, "W transverse mass [GeV]", "Events", "transverseMassWenuRecoGen", outputDIR, "Reco W(l#nu)+jets", "Gen W(l#nu)+jets", "Reco/Gen");
    drawTH1pair(hpfmet, hgen_pfmet, "PF E_{T}^{miss} [GeV]", "Events", "pfmetRecoGen", outputDIR, "Reco W(l#nu)+jets", "Gen W(l#nu)+jets", "Reco/Gen");
    drawTH1pair(hlep1pt, hgen_lep1pt, "leading electron p_{T} [GeV]", "Events", "ptCharLepRecoGen", outputDIR, "Reco W(l#nu)+jets", "Gen W(l#nu)+jets", "Reco/Gen");
    drawTH1pair(hlep2pt, hgen_lep2pt, "neutrino p_{T} (PF E_{T}^{miss}) [GeV]", "Events", "ptNuRecoGen", outputDIR, "Reco W(l#nu)+jets", "Gen W(l#nu)+jets", "Reco/Gen");
    drawTH1pair(hbosonpt, hgen_bosonpt, "boson p_{T} [GeV]", "Events", "bosonPtRecoGen", outputDIR, "Reco W(l#nu)+jets", "Gen W(l#nu)+jets", "Reco/Gen");

    drawCorrelationPlot(h2_mT_lep1pt, "W transverse mass [GeV]", "leading electron p_{T} [GeV]", "correlation_recoMC_mT_lep1pt","reco MC",outputDIR);
    drawCorrelationPlot(h2_mT_lep1eta, "W transverse mass [GeV]", "leading electron #eta", "correlation_recoMC_mT_lep1eta","reco MC",outputDIR);
    drawCorrelationPlot(h2_mT_lep1sigIetaIeta, "W transverse mass [GeV]", "leading electron #sigma_{i#etai#eta}", "correlation_recoMC_mT_lep1sigIetaIeta","reco MC",outputDIR);
    drawCorrelationPlot(h2_mT_lep1r9, "W transverse mass [GeV]", "leading electron R9", "correlation_recoMC_mT_lep1r9","reco MC",outputDIR);
    drawCorrelationPlot(h2_mT_pfmet, "W transverse mass [GeV]", "PF E_{T}^{miss} [GeV]", "correlation_recoMC_mT_pfmet","reco MC",outputDIR);
    drawCorrelationPlot(h2_mT_tkmet, "W transverse mass [GeV]", "Tracker E_{T}^{miss} [GeV]", "correlation_recoMC_mT_tkmet","reco MC",outputDIR);
    drawCorrelationPlot(h2_mT_lep1relIso03, "W transverse mass [GeV]", "leading electron isolation (relIso03)", "correlation_recoMC_mT_lep1relIso03","reco MC",outputDIR);
    drawCorrelationPlot(h2_mT_bosonPt, "W transverse mass [GeV]", "boson p_{T} [GeV]", "correlation_recoMC_mT_bosonPt","reco MC",outputDIR);

    drawCorrelationPlot(h2_lep1pt_lep1eta, "leading electron p_{T} [GeV]", "leading electron #eta", "correlation_recoMC_lep1pt_lep1eta","reco MC",outputDIR);
    drawCorrelationPlot(h2_lep1pt_lep1sigIetaIeta, "leading electron p_{T} [GeV]", "leading electron #sigma_{i#etai#eta}", "correlation_recoMC_lep1pt_lep1sigIetaIeta","reco MC",outputDIR);
    drawCorrelationPlot(h2_lep1pt_lep1r9, "leading electron p_{T} [GeV]", "leading electron R9", "correlation_recoMC_lep1pt_lep1r9","reco MC",outputDIR);
    drawCorrelationPlot(h2_lep1pt_pfmet, "leading electron p_{T} [GeV]", "PF E_{T}^{miss} [GeV]", "correlation_recoMC_lep1pt_pfmet","reco MC",outputDIR);
    drawCorrelationPlot(h2_lep1pt_tkmet, "leading electron p_{T} [GeV]", "Tracker E_{T}^{miss} [GeV]", "correlation_recoMC_lep1pt_tkmet","reco MC",outputDIR);
    drawCorrelationPlot(h2_lep1pt_lep1relIso03, "leading electron p_{T} [GeV]", "leading electron isolation (relIso03)", "correlation_recoMC_lep1pt_lep1relIso03","reco MC",outputDIR);
    drawCorrelationPlot(h2_lep1pt_bosonPt, "leading electron p_{T} [GeV]", "boson p_{T} [GeV]", "correlation_recoMC_lep1pt_bosonPt","reco MC",outputDIR);


  }

  else {
    drawTH1pair(hbosoneta, hgen_bosoneta, "boson #eta", "Events", "bosonEtaRecoGen", outputDIR, "Reco Z(ll)+jets", "Gen Z(ll)+jets", "Reco/Gen");
    drawTH1pair(hmZ, hgen_mZ, "invariant mass(e^{+}e^{-}) [GeV]", "Events", "massZeeRecoGen", outputDIR, "Reco Z(ll)+jets", "Gen Z(ll)+jets", "Reco/Gen");
    drawTH1pair(hmT, hgen_mT, "transverse mass(e^{+}e^{-}) [GeV]", "Events", "transverseMassZeeRecoGen", outputDIR, "Reco Z(ll)+jets", "Gen Z(ll)+jets", "Reco/Gen");
    drawTH1pair(hpfmet, hgen_pfmet, "PF E_{T}^{miss} [GeV]", "Events", "pfmetRecoGen", outputDIR, "Reco Z(ll)+jets", "Gen Z(ll)+jets", "Reco/Gen");
    drawTH1pair(hlep1pt, hgen_lep1pt, "leading electron p_{T} [GeV]", "Events", "lep1ptRecoGen", outputDIR, "Reco Z(ll)+jets", "Gen Z(ll)+jets", "Reco/Gen");
    drawTH1pair(hlep2pt, hgen_lep2pt, "trailing electron p_{T} [GeV]", "Events", "lep2ptRecoGen", outputDIR, "Reco Z(ll)+jets", "Gen Z(ll)+jets", "Reco/Gen");
    drawTH1pair(hlep2pt, hgen_lep2pt, "trailing electron p_{T} [GeV]", "Events", "lep2ptRecoGen", outputDIR, "Reco Z(ll)+jets", "Gen Z(ll)+jets", "Reco/Gen");
    drawTH1pair(hbosonpt, hgen_bosonpt, "boson p_{T} [GeV]", "Events", "bosonPtRecoGen", outputDIR, "Reco Z(ll)+jets", "Gen Z(ll)+jets", "Reco/Gen");

    drawCorrelationPlot(h2_mT_mZ, "transverse mass(e^{+}e^{-}) [GeV]", "invariant mass(e^{+}e^{-}) [GeV]", "correlation_recoMC_mT_mZ","reco MC",outputDIR);
    drawCorrelationPlot(h2_mT_lep1pt, "transverse mass(e^{+}e^{-}) [GeV]", "leading electron p_{T} [GeV]", "correlation_recoMC_mT_lep1pt","reco MC",outputDIR);
    drawCorrelationPlot(h2_mT_lep1eta, "transverse mass(e^{+}e^{-}) [GeV]", "leading electron #eta", "correlation_recoMC_mT_lep1eta","reco MC",outputDIR);
    drawCorrelationPlot(h2_mT_lep1sigIetaIeta, "transverse mass(e^{+}e^{-}) [GeV]", "leading electron #sigma_{i#etai#eta}", "correlation_recoMC_mT_lep1sigIetaIeta","reco MC",outputDIR);
    drawCorrelationPlot(h2_mT_lep1r9, "transverse mass(e^{+}e^{-}) [GeV]", "leading electron R9", "correlation_recoMC_mT_lep1r9","reco MC",outputDIR);
    drawCorrelationPlot(h2_mT_pfmet, "transverse mass(e^{+}e^{-}) [GeV]", "PF E_{T}^{miss} [GeV]", "correlation_recoMC_mT_pfmet","reco MC",outputDIR);
    drawCorrelationPlot(h2_mT_tkmet, "transverse mass(e^{+}e^{-}) [GeV]", "Tracker E_{T}^{miss} [GeV]", "correlation_recoMC_mT_tkmet","reco MC",outputDIR);
    drawCorrelationPlot(h2_mT_lep1relIso03, "transverse mass(e^{+}e^{-}) [GeV]", "leading electron isolation (relIso03)", "correlation_recoMC_mT_lep1relIso03","reco MC",outputDIR);
    drawCorrelationPlot(h2_mT_bosonPt, "transverse mass(e^{+}e^{-}) [GeV]", "boson p_{T} [GeV]", "correlation_recoMC_mT_bosonPt","reco MC",outputDIR);

    drawCorrelationPlot(h2_lep1pt_mZ, "leading electron p_{T} [GeV]", "invariant mass(e^{+}e^{-}) [GeV]", "correlation_recoMC_lep1pt_mZ","reco MC",outputDIR);
    drawCorrelationPlot(h2_lep1pt_lep1eta, "leading electron p_{T} [GeV]", "leading electron #eta", "correlation_recoMC_lep1pt_lep1eta","reco MC",outputDIR);
    drawCorrelationPlot(h2_lep1pt_lep1sigIetaIeta, "leading electron p_{T} [GeV]", "leading electron #sigma_{i#etai#eta}", "correlation_recoMC_lep1pt_lep1sigIetaIeta","reco MC",outputDIR);
    drawCorrelationPlot(h2_lep1pt_lep1r9, "leading electron p_{T} [GeV]", "leading electron R9", "correlation_recoMC_lep1pt_lep1r9","reco MC",outputDIR);
    drawCorrelationPlot(h2_lep1pt_pfmet, "leading electron p_{T} [GeV]", "PF E_{T}^{miss} [GeV]", "correlation_recoMC_lep1pt_pfmet","reco MC",outputDIR);
    drawCorrelationPlot(h2_lep1pt_tkmet, "leading electron p_{T} [GeV]", "Tracker E_{T}^{miss} [GeV]", "correlation_recoMC_lep1pt_tkmet","reco MC",outputDIR);
    drawCorrelationPlot(h2_lep1pt_lep1relIso03, "leading electron p_{T} [GeV]", "leading electron isolation (relIso03)", "correlation_recoMC_lep1pt_lep1relIso03","reco MC",outputDIR);
    drawCorrelationPlot(h2_lep1pt_bosonPt, "leading electron p_{T} [GeV]", "boson p_{T} [GeV]", "correlation_recoMC_lep1pt_bosonPt","reco MC",outputDIR);

  }

  cout << endl;

}



//=====================================================================

// useW = false to do Z

void makeRecoGenHistograms(const string& inputDIR = "./", const string& outputDIR = "./", const Bool_t& isWregion = false) {

  if (outputDIR!= "./") system(("mkdir -p "+outputDIR).c_str());

  gROOT->SetBatch(kTRUE);
  TH1::SetDefaultSumw2(); //all the following histograms will automatically call TH1::Sumw2()  

  cout << endl;

  if (isWregion) fillRecoGenHistograms(inputDIR, outputDIR, Sample::wjets, isWregion);
  else      fillRecoGenHistograms(inputDIR, outputDIR, Sample::zjets, isWregion);

}
