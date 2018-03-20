#ifndef centerOfMassEnergy_h
#define centerOfMassEnergy_h

static Bool_t use8TeVSample = false;
static Int_t centerOfMassEnergy = use8TeVSample ? 8 : 13;
static Double_t intLumi = use8TeVSample ? 19.7 : 19.3; // 35.9 or 36.4 with 2016 G and H included
// following are more for the maco that loops on events
static Bool_t useTrackMet = true;
static Bool_t useAbsIso = false;  // rel iso or abs iso to cut (abs iso is rel iso times lepton pT)          
static bool useFakeRateForElectron = false;
static bool useFakeRateForMuon = false;


#endif
