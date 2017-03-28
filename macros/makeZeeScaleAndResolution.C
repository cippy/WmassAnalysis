//#include "makeZeeScaleAndResolutionUtils.h"
#include "utility.h"
#include <fstream>


// use leading lepton (tag) to divide in categories: its pT must be > 25 GeV
static vector<Double_t> ptBins = {25, 30, 35, 40, 45, 50, 55, 60};
static vector<Double_t> etaBins = {0.0, 0.8, 1.479};
static Int_t nPtBins = (Int_t) ptBins.size() -1;
static Int_t nEtaBins = (Int_t) etaBins.size() -1;

static Int_t nBinsMassZee = 61;
static Double_t massZeeMin = 59.5;
static Double_t massZeeMax = 120.5;

static Double_t intLumi = 36.4;

//static vector<float> nBinsMassZee = 


//======================================================

class zeeScale {
public:
  zeeScale(const Double_t& ptMin = 0, const Double_t& ptMax = 1, const Double_t& etaMin = 0, const Double_t& etaMax = 1, TH1D* histo = NULL):
    ptMin(ptMin),
    ptMax(ptMax),
    etaMin(etaMin),
    etaMax(etaMax),
    zeeHisto(histo){
    ptMean = 0;
  }
    ~zeeScale() {};

  float ptMin;
  float ptMax;
  float ptMean;

  float etaMin;
  float etaMax;

  TH1D* zeeHisto;

  void print() const {
    cout << "histo --> " << zeeHisto->GetName() << "\t pt --> [" << ptMin << "," << ptMax << "]\t eta --> [" << etaMin << "," << etaMax << "]" << endl; 
  }
};

//======================================================


string getStringFromDouble(const Double_t& num = 1.0, const Double_t epsilon = 0.00001) {

  Int_t i = (Int_t) num;
  // stringstream ss; 
  // ss<<i; 
  // string numStr = ss.str();
  Int_t int_decim = (Int_t) (1000 * (num - (Double_t) i + epsilon));

  if      (int_decim%1000 == 0) return string(Form("%dp%d",i,int_decim/1000));
  else if (int_decim%100 == 0)  return string(Form("%dp%d",i,int_decim/100));
  else if (int_decim%10 == 0)   return string(Form("%dp%d",i,int_decim/10));
  else                          return string(Form("%dp%d",i,int_decim));

}

//======================================================

Int_t getHistIndexByPtEta(const Double_t& pt = -1.0, const Double_t& abseta = -1.0) {

  Bool_t ptFound = false;
  Bool_t etaFound = false;

  Int_t ipt = 0;
  Int_t ieta = 0;

  if (pt < ptBins.front() || pt > ptBins.back()) return -1;
  if (abseta < etaBins.front() || abseta > etaBins.back()) return -2;

  while ( !ptFound && (ipt < nPtBins) ) {

    if (pt <= ptBins[ipt+1]) ptFound = true;
    ipt++;
    
  }

  while ( !etaFound && (ieta < nEtaBins) ) {

    if (abseta <= etaBins[ieta+1]) etaFound = true;
    ieta++;
    
  }

  // subtract 1 from ipt and ieta so that bin 0 of vector (first vector bin) is returned when the value is in the first pt and eta bin
  return (ipt -1) * nEtaBins + (ieta -1);   

} 

//======================================================

void fillHistograms(vector<zeeScale>& vecHist, const string& inputDIR = "./", const string& outputDIR = "./", const Sample& sample = Sample::zjets) {

  gROOT->SetBatch(kTRUE);
  TH1::SetDefaultSumw2(); //all the following histograms will automatically call TH1::Sumw2()                               

  cout << endl;

  TChain* chain = new TChain("tree");
  TChain* chFriend = new TChain("mjvars/t");

  buildChain(chain, chFriend, inputDIR, sample); 

  TTreeReader reader (chain);

  // reco met
  TTreeReaderValue<Float_t> tkmet    (reader,"met_trkPt");
  TTreeReaderValue<Float_t> tkmet_phi(reader,"met_trkPhi");
  TTreeReaderValue<Float_t> pfmet    (reader,"met_pt");
  TTreeReaderValue<Float_t> pfmet_phi(reader,"met_phi");

  // gen met
  // TTreeReaderValue<Float_t> genpfmet    (reader,"met_genPt");
  // TTreeReaderValue<Float_t> genpfmet_phi(reader,"met_genPhi");

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

  // mZ
  TTreeReaderValue<Float_t> mZ1 (reader,"mZ1");

  // gen particles
  // TTreeReaderValue<Int_t> nGenPart(reader,"nGenPart");
  // TTreeReaderArray<Int_t> GenPart_pdgId(reader,"GenPart_pdgId");
  // TTreeReaderArray<Int_t> GenPart_motherId(reader,"GenPart_motherId");
  // TTreeReaderArray<Int_t> GenPart_motherIndex(reader,"GenPart_motherIndex");
  // TTreeReaderArray<Float_t> GenPart_pt(reader,"GenPart_pt");
  // TTreeReaderArray<Float_t> GenPart_eta(reader,"GenPart_eta");
  // TTreeReaderArray<Float_t> GenPart_phi(reader,"GenPart_phi");
  // TTreeReaderArray<Float_t> GenPart_mass(reader,"GenPart_mass");

  // MC reweight
  TTreeReaderValue<Float_t> weight(reader,"weight");
  TTreeReaderValue<Float_t> puw(reader,"puw");


  // start event loop                                                                                                                                                       
  long int nTotal = chain->GetEntries();
  long int nEvents = 0;
  
  Double_t wgt = 1.0;

  while(reader.Next()){
  
    cout.flush();
    if(nEvents % 100000 == 0) cout<<"\r"<<"Analyzing events "<<double(nEvents)/nTotal*100<<" % ";
    nEvents++;

    // selection
    if ( *mZ1 < 60 || *mZ1 > 120) continue;
    if (*nlep != 2) continue;    // 2 leptons
    if (fabs(lep_pdgId[0]) != 11 || fabs(lep_pdgId[1]) != 11) continue;  // electrons
    if ( (lep_pdgId[0] + lep_pdgId[1]) != 0) continue;   // same sign
    if ( lep_pt[0] < 25.0 || lep_tightId[0] < 3 || lep_relIso03[0] > 0.0588 ||  lep_relIso03[0] < 0) continue;
    if ( lep_pt[1] < 10.0 || lep_relIso03[1] > 0.175  || lep_relIso03[1]  < 0.0) continue;
    if ( fabs(lep_eta[0]) > 1.479 || fabs(lep_eta[1]) > 1.479 ) continue;

    if (sample == Sample::data_doubleEG) wgt = 1.0;
    else wgt = intLumi * *puw * *weight;

    TLorentzVector lep1Reco, lep2Reco, bosonreco;
    lep1Reco.SetPtEtaPhiM(lep_pt[0],lep_eta[0],lep_phi[0],lep_mass[0]); 
    lep2Reco.SetPtEtaPhiM(lep_pt[1],lep_eta[1],lep_phi[1],lep_mass[1]); 
    bosonreco = lep1Reco + lep2Reco;

    if (bosonreco.Pt() > 30) continue;

    //Double_t mT = sqrt(2. * lep_pt[0] * lep_pt[1] * (1. - cos(lep1Reco.DeltaPhi(lep2Reco)) ) );

    //if (mT < 60 || mT > 100) continue;    

    Int_t vecBin = getHistIndexByPtEta(lep_pt[0], fabs(lep_eta[0]));
    if (vecBin < 0) continue;

    fillTH1(vecHist.at(vecBin).zeeHisto,(Double_t) *mZ1, wgt);

  }

  cout << endl;

}


//=============================================================

void makeZeeScaleAndResolution(const string& inputDIR = "./", const string& outputDIR = "./", const Sample& sample = Sample::zjets) {

  if (outputDIR!= "./") system(("mkdir -p "+outputDIR).c_str());
  
  gROOT->SetBatch(kTRUE);
  TH1::SetDefaultSumw2(); //all the following histograms will automatically call TH1::Sumw2()                               
  cout << endl;

  // string etaStr = "";
  // etaStr = getStringFromDouble(1.47956632);
  // cout <<  etaStr << endl;
  // etaStr = getStringFromDouble(1.479);
  // cout <<  etaStr << endl;
  // etaStr = getStringFromDouble(1.4);
  // cout <<  etaStr << endl;
  // etaStr = getStringFromDouble(0.80);
  // cout <<  etaStr << endl;

  vector<zeeScale> dataHist;
  vector<zeeScale> mcHist;

  Int_t ptBinsSize = (Int_t) ptBins.size() -1;
  Int_t etaBinsSize = (Int_t) etaBins.size() -1;

  for (Int_t ptbin = 0; ptbin < ptBinsSize; ptbin++) {
    for (Int_t etabin = 0; etabin < etaBinsSize; etabin++) {

      dataHist.push_back(zeeScale(ptBins[ptbin],ptBins[ptbin+1],
				  etaBins[etabin],etaBins[etabin+1],
				  new TH1D(Form("dataHist_pt_%d_%d_eta_%s_%s",
						Int_t(ptBins[ptbin]), Int_t(ptBins[ptbin+1]),
						getStringFromDouble(etaBins[etabin]).c_str(), getStringFromDouble(etaBins[etabin+1]).c_str()),"",
					   nBinsMassZee, massZeeMin, massZeeMax)
				  )
			 );
      dataHist.back().print();
      mcHist.push_back(zeeScale(ptBins[ptbin],ptBins[ptbin+1],
				  etaBins[etabin],etaBins[etabin+1],
				  new TH1D(Form("mcHist_pt_%d_%d_eta_%s_%s",
						Int_t(ptBins[ptbin]), Int_t(ptBins[ptbin+1]),
						getStringFromDouble(etaBins[etabin]).c_str(), getStringFromDouble(etaBins[etabin+1]).c_str()),"",
					   nBinsMassZee, massZeeMin, massZeeMax)
				  )
			 );
      mcHist.back().print();

    }
  }

  // cout << "dataHist.size() = " << dataHist.size() << endl;
  // Int_t bin = -8;
  // bin = getHistIndexByPtEta(37, 0.5);
  // cout << "bin = " << bin << endl; 
  // bin = getHistIndexByPtEta(10, 0.5);
  // cout << "bin = " << bin << endl; 
  // bin = getHistIndexByPtEta(30, 2.5);
  // cout << "bin = " << bin << endl; 
  // bin = getHistIndexByPtEta(27, 0.5);
  // cout << "bin = " << bin << endl; 
  // bin = getHistIndexByPtEta(57, 1.4);
  // cout << "bin = " << bin << endl; 

  fillHistograms(dataHist, inputDIR, outputDIR, Sample::data_doubleEG);
  fillHistograms(mcHist, inputDIR, outputDIR, Sample::zjets);

  for (Int_t ptbin = 0; ptbin < ptBinsSize; ptbin++) {
    for (Int_t etabin = 0; etabin < etaBinsSize; etabin++) {
      
      string canvasTitle = string(Form("ZeeMass_pt_%d_%d_eta_%s_%s",
				       Int_t(ptBins[ptbin]), Int_t(ptBins[ptbin+1]),
				       getStringFromDouble(etaBins[etabin]).c_str(), getStringFromDouble(etaBins[etabin+1]).c_str()));

      drawTH1pair(dataHist.at(ptbin * etaBinsSize + etabin).zeeHisto, mcHist.at(ptbin * etaBinsSize + etabin).zeeHisto, 
		  "invariant mass(e^{+}e^{-}) [GeV]", "Events", canvasTitle, outputDIR, "data", "Z(ll)+jets", "data/MC", intLumi);

    }
  }

  cout << endl;

}
