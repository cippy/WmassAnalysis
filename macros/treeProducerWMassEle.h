//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Mar  4 14:59:37 2017 by ROOT version 6.06/01
// from TTree treeProducerWMassEle/treeProducerWMassEle
// found on file: /u2/emanuele/TREES_1LEP_53X_V1/DYJetsM50/treeProducerWMassEle/treeProducerWMassEle_tree.root
//////////////////////////////////////////////////////////

#ifndef treeProducerWMassEle_h
#define treeProducerWMassEle_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class treeProducerWMassEle {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           run;
   Int_t           lumi;
   Int_t           evt;
   Int_t           isData;
   Int_t           HLT_SingleMu;
   Int_t           HLT_SingleEl;
   Int_t           HLT_MuEG;
   Int_t           HLT_TripleEl;
   Int_t           HLT_DoubleEl;
   Int_t           HLT_DoubleMu;
   Float_t         puWeight;
   Float_t         pdfWeight_CT10[53];
   Float_t         pdfWeight_MSTW2008lo68cl[41];
   Float_t         pdfWeight_NNPDF21_100[101];
   Float_t         rho;
   Int_t           nVert;
   Int_t           GenHeaviestQCDFlavour;
   Float_t         LepEff_1lep;
   Float_t         LepEff_2lep;
   Float_t         mZ1;
   Float_t         m2l;
   Float_t         tkmet_pt;
   Float_t         tkmet_eta;
   Float_t         tkmet_phi;
   Float_t         tkmet_mass;
   Float_t         pucmet_pt;
   Float_t         pucmet_eta;
   Float_t         pucmet_phi;
   Float_t         pucmet_mass;
   Float_t         met_pt;
   Float_t         met_eta;
   Float_t         met_phi;
   Float_t         met_mass;
   Float_t         met_sumEt;
   Float_t         met_genPt;
   Float_t         met_genPhi;
   Float_t         met_genEta;
   Float_t         metraw_pt;
   Float_t         metraw_eta;
   Float_t         metraw_phi;
   Float_t         metraw_mass;
   Float_t         metraw_sumEt;
   Float_t         metraw_genPt;
   Float_t         metraw_genPhi;
   Float_t         metraw_genEta;
   Float_t         pumet_pt;
   Float_t         pumet_eta;
   Float_t         pumet_phi;
   Float_t         pumet_mass;
   Float_t         metNoPU_pt;
   Float_t         metNoPU_eta;
   Float_t         metNoPU_phi;
   Float_t         metNoPU_mass;
   Int_t           nLepOther;
   Float_t         LepOther_pt[7];   //[nLepOther]
   Float_t         LepOther_eta[7];   //[nLepOther]
   Float_t         LepOther_phi[7];   //[nLepOther]
   Float_t         LepOther_mass[7];   //[nLepOther]
   Int_t           LepOther_pdgId[7];   //[nLepOther]
   Int_t           LepOther_charge[7];   //[nLepOther]
   Float_t         LepOther_dxy[7];   //[nLepOther]
   Float_t         LepOther_dz[7];   //[nLepOther]
   Float_t         LepOther_edxy[7];   //[nLepOther]
   Float_t         LepOther_edz[7];   //[nLepOther]
   Float_t         LepOther_ip3d[7];   //[nLepOther]
   Float_t         LepOther_sip3d[7];   //[nLepOther]
   Int_t           LepOther_tightId[7];   //[nLepOther]
   Int_t           LepOther_convVeto[7];   //[nLepOther]
   Int_t           LepOther_lostHits[7];   //[nLepOther]
   Int_t           LepOther_looseIdSusy[7];   //[nLepOther]
   Float_t         LepOther_relIso03[7];   //[nLepOther]
   Float_t         LepOther_relIso04[7];   //[nLepOther]
   Float_t         LepOther_chargedHadRelIso03[7];   //[nLepOther]
   Float_t         LepOther_chargedHadRelIso04[7];   //[nLepOther]
   Int_t           LepOther_convVetoFull[7];   //[nLepOther]
   Int_t           LepOther_eleCutId[7];   //[nLepOther]
   Int_t           LepOther_eleMVAId[7];   //[nLepOther]
   Int_t           LepOther_tightCharge[7];   //[nLepOther]
   Float_t         LepOther_mvaId[7];   //[nLepOther]
   Float_t         LepOther_mvaIdTrig[7];   //[nLepOther]
   Float_t         LepOther_nStations[7];   //[nLepOther]
   Float_t         LepOther_trkKink[7];   //[nLepOther]
   Float_t         LepOther_caloCompatibility[7];   //[nLepOther]
   Float_t         LepOther_globalTrackChi2[7];   //[nLepOther]
   Int_t           LepOther_trackerLayers[7];   //[nLepOther]
   Int_t           LepOther_pixelLayers[7];   //[nLepOther]
   Float_t         LepOther_mvaTTH[7];   //[nLepOther]
   Float_t         LepOther_jetPtRatio[7];   //[nLepOther]
   Float_t         LepOther_jetBTagCSV[7];   //[nLepOther]
   Float_t         LepOther_jetDR[7];   //[nLepOther]
   Int_t           LepOther_mcMatchId[7];   //[nLepOther]
   Int_t           LepOther_mcMatchAny[7];   //[nLepOther]
   Int_t           LepOther_mcMatchAny2[7];   //[nLepOther]
   Int_t           LepOther_mcMatchTau[7];   //[nLepOther]
   Int_t           LepOther_softMuID[7];   //[nLepOther]
   Float_t         LepOther_trgMatch[7];   //[nLepOther]
   Float_t         LepOther_scEta[7];   //[nLepOther]
   Float_t         LepOther_r9[7];   //[nLepOther]
   Float_t         LepOther_classification[7];   //[nLepOther]
   Float_t         LepOther_detaIn[7];   //[nLepOther]
   Float_t         LepOther_dphiIn[7];   //[nLepOther]
   Float_t         LepOther_sigmaIetaIeta[7];   //[nLepOther]
   Float_t         LepOther_sigmaIphiIphi[7];   //[nLepOther]
   Float_t         LepOther_hcalOverEcal[7];   //[nLepOther]
   Float_t         LepOther_correctedEcalEnergyError[7];   //[nLepOther]
   Float_t         LepOther_regressionEnergyError[7];   //[nLepOther]
   Float_t         LepOther_ecalEnergy[7];   //[nLepOther]
   Float_t         LepOther_superCluster_rawEnergy[7];   //[nLepOther]
   Float_t         LepOther_superCluster_preshowerEnergy[7];   //[nLepOther]
   Float_t         LepOther_superCluster_energy[7];   //[nLepOther]
   Float_t         LepOther_eSuperClusterOverP[7];   //[nLepOther]
   Int_t           nLepGood;
   Float_t         LepGood_pt[5];   //[nLepGood]
   Float_t         LepGood_eta[5];   //[nLepGood]
   Float_t         LepGood_phi[5];   //[nLepGood]
   Float_t         LepGood_mass[5];   //[nLepGood]
   Int_t           LepGood_pdgId[5];   //[nLepGood]
   Int_t           LepGood_charge[5];   //[nLepGood]
   Float_t         LepGood_dxy[5];   //[nLepGood]
   Float_t         LepGood_dz[5];   //[nLepGood]
   Float_t         LepGood_edxy[5];   //[nLepGood]
   Float_t         LepGood_edz[5];   //[nLepGood]
   Float_t         LepGood_ip3d[5];   //[nLepGood]
   Float_t         LepGood_sip3d[5];   //[nLepGood]
   Int_t           LepGood_tightId[5];   //[nLepGood]
   Int_t           LepGood_convVeto[5];   //[nLepGood]
   Int_t           LepGood_lostHits[5];   //[nLepGood]
   Int_t           LepGood_looseIdSusy[5];   //[nLepGood]
   Float_t         LepGood_relIso03[5];   //[nLepGood]
   Float_t         LepGood_relIso04[5];   //[nLepGood]
   Float_t         LepGood_chargedHadRelIso03[5];   //[nLepGood]
   Float_t         LepGood_chargedHadRelIso04[5];   //[nLepGood]
   Int_t           LepGood_convVetoFull[5];   //[nLepGood]
   Int_t           LepGood_eleCutId[5];   //[nLepGood]
   Int_t           LepGood_eleMVAId[5];   //[nLepGood]
   Int_t           LepGood_tightCharge[5];   //[nLepGood]
   Float_t         LepGood_mvaId[5];   //[nLepGood]
   Float_t         LepGood_mvaIdTrig[5];   //[nLepGood]
   Float_t         LepGood_nStations[5];   //[nLepGood]
   Float_t         LepGood_trkKink[5];   //[nLepGood]
   Float_t         LepGood_caloCompatibility[5];   //[nLepGood]
   Float_t         LepGood_globalTrackChi2[5];   //[nLepGood]
   Int_t           LepGood_trackerLayers[5];   //[nLepGood]
   Int_t           LepGood_pixelLayers[5];   //[nLepGood]
   Float_t         LepGood_mvaTTH[5];   //[nLepGood]
   Float_t         LepGood_jetPtRatio[5];   //[nLepGood]
   Float_t         LepGood_jetBTagCSV[5];   //[nLepGood]
   Float_t         LepGood_jetDR[5];   //[nLepGood]
   Int_t           LepGood_mcMatchId[5];   //[nLepGood]
   Int_t           LepGood_mcMatchAny[5];   //[nLepGood]
   Int_t           LepGood_mcMatchAny2[5];   //[nLepGood]
   Int_t           LepGood_mcMatchTau[5];   //[nLepGood]
   Int_t           LepGood_softMuID[5];   //[nLepGood]
   Float_t         LepGood_trgMatch[5];   //[nLepGood]
   Float_t         LepGood_scEta[5];   //[nLepGood]
   Float_t         LepGood_r9[5];   //[nLepGood]
   Float_t         LepGood_classification[5];   //[nLepGood]
   Float_t         LepGood_detaIn[5];   //[nLepGood]
   Float_t         LepGood_dphiIn[5];   //[nLepGood]
   Float_t         LepGood_sigmaIetaIeta[5];   //[nLepGood]
   Float_t         LepGood_sigmaIphiIphi[5];   //[nLepGood]
   Float_t         LepGood_hcalOverEcal[5];   //[nLepGood]
   Float_t         LepGood_correctedEcalEnergyError[5];   //[nLepGood]
   Float_t         LepGood_regressionEnergyError[5];   //[nLepGood]
   Float_t         LepGood_ecalEnergy[5];   //[nLepGood]
   Float_t         LepGood_superCluster_rawEnergy[5];   //[nLepGood]
   Float_t         LepGood_superCluster_preshowerEnergy[5];   //[nLepGood]
   Float_t         LepGood_superCluster_energy[5];   //[nLepGood]
   Float_t         LepGood_eSuperClusterOverP[5];   //[nLepGood]
   Int_t           nJetFwd;
   Float_t         JetFwd_pt[4];   //[nJetFwd]
   Float_t         JetFwd_eta[4];   //[nJetFwd]
   Float_t         JetFwd_phi[4];   //[nJetFwd]
   Float_t         JetFwd_mass[4];   //[nJetFwd]
   Float_t         JetFwd_btagCSV[4];   //[nJetFwd]
   Float_t         JetFwd_rawPt[4];   //[nJetFwd]
   Float_t         JetFwd_mcPt[4];   //[nJetFwd]
   Int_t           JetFwd_mcFlavour[4];   //[nJetFwd]
   Float_t         JetFwd_quarkGluonID[4];   //[nJetFwd]
   Int_t           JetFwd_mcMatchId[4];   //[nJetFwd]
   Int_t           JetFwd_mcMatchFlav[4];   //[nJetFwd]
   Int_t           JetFwd_puId[4];   //[nJetFwd]
   Float_t         JetFwd_area[4];   //[nJetFwd]
   Int_t           JetFwd_id[4];   //[nJetFwd]
   Float_t         JetFwd_CHEF[4];   //[nJetFwd]
   Float_t         JetFwd_NHEF[4];   //[nJetFwd]
   Float_t         JetFwd_PHEF[4];   //[nJetFwd]
   Float_t         JetFwd_MUEF[4];   //[nJetFwd]
   Float_t         JetFwd_ELEF[4];   //[nJetFwd]
   Int_t           nJet;
   Float_t         Jet_pt[8];   //[nJet]
   Float_t         Jet_eta[8];   //[nJet]
   Float_t         Jet_phi[8];   //[nJet]
   Float_t         Jet_mass[8];   //[nJet]
   Float_t         Jet_btagCSV[8];   //[nJet]
   Float_t         Jet_rawPt[8];   //[nJet]
   Float_t         Jet_mcPt[8];   //[nJet]
   Int_t           Jet_mcFlavour[8];   //[nJet]
   Float_t         Jet_quarkGluonID[8];   //[nJet]
   Int_t           Jet_mcMatchId[8];   //[nJet]
   Int_t           Jet_mcMatchFlav[8];   //[nJet]
   Int_t           Jet_puId[8];   //[nJet]
   Float_t         Jet_area[8];   //[nJet]
   Int_t           Jet_id[8];   //[nJet]
   Float_t         Jet_CHEF[8];   //[nJet]
   Float_t         Jet_NHEF[8];   //[nJet]
   Float_t         Jet_PHEF[8];   //[nJet]
   Float_t         Jet_MUEF[8];   //[nJet]
   Float_t         Jet_ELEF[8];   //[nJet]
   Int_t           nGenLep;
   Float_t         GenLep_pt[2];   //[nGenLep]
   Float_t         GenLep_eta[2];   //[nGenLep]
   Float_t         GenLep_phi[2];   //[nGenLep]
   Float_t         GenLep_mass[2];   //[nGenLep]
   Int_t           GenLep_pdgId[2];   //[nGenLep]
   Float_t         GenLep_charge[2];   //[nGenLep]
   Int_t           GenLep_sourceId[2];   //[nGenLep]
   Int_t           nGenP6StatusThree;
   Float_t         GenP6StatusThree_pt[13];   //[nGenP6StatusThree]
   Float_t         GenP6StatusThree_eta[13];   //[nGenP6StatusThree]
   Float_t         GenP6StatusThree_phi[13];   //[nGenP6StatusThree]
   Float_t         GenP6StatusThree_mass[13];   //[nGenP6StatusThree]
   Int_t           GenP6StatusThree_pdgId[13];   //[nGenP6StatusThree]
   Float_t         GenP6StatusThree_charge[13];   //[nGenP6StatusThree]
   Int_t           GenP6StatusThree_motherId[13];   //[nGenP6StatusThree]
   Int_t           GenP6StatusThree_grandmaId[13];   //[nGenP6StatusThree]

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_evt;   //!
   TBranch        *b_isData;   //!
   TBranch        *b_HLT_SingleMu;   //!
   TBranch        *b_HLT_SingleEl;   //!
   TBranch        *b_HLT_MuEG;   //!
   TBranch        *b_HLT_TripleEl;   //!
   TBranch        *b_HLT_DoubleEl;   //!
   TBranch        *b_HLT_DoubleMu;   //!
   TBranch        *b_puWeight;   //!
   TBranch        *b_pdfWeight_CT10;   //!
   TBranch        *b_pdfWeight_MSTW2008lo68cl;   //!
   TBranch        *b_pdfWeight_NNPDF21_100;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_nVert;   //!
   TBranch        *b_GenHeaviestQCDFlavour;   //!
   TBranch        *b_LepEff_1lep;   //!
   TBranch        *b_LepEff_2lep;   //!
   TBranch        *b_mZ1;   //!
   TBranch        *b_m2l;   //!
   TBranch        *b_tkmet_pt;   //!
   TBranch        *b_tkmet_eta;   //!
   TBranch        *b_tkmet_phi;   //!
   TBranch        *b_tkmet_mass;   //!
   TBranch        *b_pucmet_pt;   //!
   TBranch        *b_pucmet_eta;   //!
   TBranch        *b_pucmet_phi;   //!
   TBranch        *b_pucmet_mass;   //!
   TBranch        *b_met_pt;   //!
   TBranch        *b_met_eta;   //!
   TBranch        *b_met_phi;   //!
   TBranch        *b_met_mass;   //!
   TBranch        *b_met_sumEt;   //!
   TBranch        *b_met_genPt;   //!
   TBranch        *b_met_genPhi;   //!
   TBranch        *b_met_genEta;   //!
   TBranch        *b_metraw_pt;   //!
   TBranch        *b_metraw_eta;   //!
   TBranch        *b_metraw_phi;   //!
   TBranch        *b_metraw_mass;   //!
   TBranch        *b_metraw_sumEt;   //!
   TBranch        *b_metraw_genPt;   //!
   TBranch        *b_metraw_genPhi;   //!
   TBranch        *b_metraw_genEta;   //!
   TBranch        *b_pumet_pt;   //!
   TBranch        *b_pumet_eta;   //!
   TBranch        *b_pumet_phi;   //!
   TBranch        *b_pumet_mass;   //!
   TBranch        *b_metNoPU_pt;   //!
   TBranch        *b_metNoPU_eta;   //!
   TBranch        *b_metNoPU_phi;   //!
   TBranch        *b_metNoPU_mass;   //!
   TBranch        *b_nLepOther;   //!
   TBranch        *b_LepOther_pt;   //!
   TBranch        *b_LepOther_eta;   //!
   TBranch        *b_LepOther_phi;   //!
   TBranch        *b_LepOther_mass;   //!
   TBranch        *b_LepOther_pdgId;   //!
   TBranch        *b_LepOther_charge;   //!
   TBranch        *b_LepOther_dxy;   //!
   TBranch        *b_LepOther_dz;   //!
   TBranch        *b_LepOther_edxy;   //!
   TBranch        *b_LepOther_edz;   //!
   TBranch        *b_LepOther_ip3d;   //!
   TBranch        *b_LepOther_sip3d;   //!
   TBranch        *b_LepOther_tightId;   //!
   TBranch        *b_LepOther_convVeto;   //!
   TBranch        *b_LepOther_lostHits;   //!
   TBranch        *b_LepOther_looseIdSusy;   //!
   TBranch        *b_LepOther_relIso03;   //!
   TBranch        *b_LepOther_relIso04;   //!
   TBranch        *b_LepOther_chargedHadRelIso03;   //!
   TBranch        *b_LepOther_chargedHadRelIso04;   //!
   TBranch        *b_LepOther_convVetoFull;   //!
   TBranch        *b_LepOther_eleCutId;   //!
   TBranch        *b_LepOther_eleMVAId;   //!
   TBranch        *b_LepOther_tightCharge;   //!
   TBranch        *b_LepOther_mvaId;   //!
   TBranch        *b_LepOther_mvaIdTrig;   //!
   TBranch        *b_LepOther_nStations;   //!
   TBranch        *b_LepOther_trkKink;   //!
   TBranch        *b_LepOther_caloCompatibility;   //!
   TBranch        *b_LepOther_globalTrackChi2;   //!
   TBranch        *b_LepOther_trackerLayers;   //!
   TBranch        *b_LepOther_pixelLayers;   //!
   TBranch        *b_LepOther_mvaTTH;   //!
   TBranch        *b_LepOther_jetPtRatio;   //!
   TBranch        *b_LepOther_jetBTagCSV;   //!
   TBranch        *b_LepOther_jetDR;   //!
   TBranch        *b_LepOther_mcMatchId;   //!
   TBranch        *b_LepOther_mcMatchAny;   //!
   TBranch        *b_LepOther_mcMatchAny2;   //!
   TBranch        *b_LepOther_mcMatchTau;   //!
   TBranch        *b_LepOther_softMuID;   //!
   TBranch        *b_LepOther_trgMatch;   //!
   TBranch        *b_LepOther_scEta;   //!
   TBranch        *b_LepOther_r9;   //!
   TBranch        *b_LepOther_classification;   //!
   TBranch        *b_LepOther_detaIn;   //!
   TBranch        *b_LepOther_dphiIn;   //!
   TBranch        *b_LepOther_sigmaIetaIeta;   //!
   TBranch        *b_LepOther_sigmaIphiIphi;   //!
   TBranch        *b_LepOther_hcalOverEcal;   //!
   TBranch        *b_LepOther_correctedEcalEnergyError;   //!
   TBranch        *b_LepOther_regressionEnergyError;   //!
   TBranch        *b_LepOther_ecalEnergy;   //!
   TBranch        *b_LepOther_superCluster_rawEnergy;   //!
   TBranch        *b_LepOther_superCluster_preshowerEnergy;   //!
   TBranch        *b_LepOther_superCluster_energy;   //!
   TBranch        *b_LepOther_eSuperClusterOverP;   //!
   TBranch        *b_nLepGood;   //!
   TBranch        *b_LepGood_pt;   //!
   TBranch        *b_LepGood_eta;   //!
   TBranch        *b_LepGood_phi;   //!
   TBranch        *b_LepGood_mass;   //!
   TBranch        *b_LepGood_pdgId;   //!
   TBranch        *b_LepGood_charge;   //!
   TBranch        *b_LepGood_dxy;   //!
   TBranch        *b_LepGood_dz;   //!
   TBranch        *b_LepGood_edxy;   //!
   TBranch        *b_LepGood_edz;   //!
   TBranch        *b_LepGood_ip3d;   //!
   TBranch        *b_LepGood_sip3d;   //!
   TBranch        *b_LepGood_tightId;   //!
   TBranch        *b_LepGood_convVeto;   //!
   TBranch        *b_LepGood_lostHits;   //!
   TBranch        *b_LepGood_looseIdSusy;   //!
   TBranch        *b_LepGood_relIso03;   //!
   TBranch        *b_LepGood_relIso04;   //!
   TBranch        *b_LepGood_chargedHadRelIso03;   //!
   TBranch        *b_LepGood_chargedHadRelIso04;   //!
   TBranch        *b_LepGood_convVetoFull;   //!
   TBranch        *b_LepGood_eleCutId;   //!
   TBranch        *b_LepGood_eleMVAId;   //!
   TBranch        *b_LepGood_tightCharge;   //!
   TBranch        *b_LepGood_mvaId;   //!
   TBranch        *b_LepGood_mvaIdTrig;   //!
   TBranch        *b_LepGood_nStations;   //!
   TBranch        *b_LepGood_trkKink;   //!
   TBranch        *b_LepGood_caloCompatibility;   //!
   TBranch        *b_LepGood_globalTrackChi2;   //!
   TBranch        *b_LepGood_trackerLayers;   //!
   TBranch        *b_LepGood_pixelLayers;   //!
   TBranch        *b_LepGood_mvaTTH;   //!
   TBranch        *b_LepGood_jetPtRatio;   //!
   TBranch        *b_LepGood_jetBTagCSV;   //!
   TBranch        *b_LepGood_jetDR;   //!
   TBranch        *b_LepGood_mcMatchId;   //!
   TBranch        *b_LepGood_mcMatchAny;   //!
   TBranch        *b_LepGood_mcMatchAny2;   //!
   TBranch        *b_LepGood_mcMatchTau;   //!
   TBranch        *b_LepGood_softMuID;   //!
   TBranch        *b_LepGood_trgMatch;   //!
   TBranch        *b_LepGood_scEta;   //!
   TBranch        *b_LepGood_r9;   //!
   TBranch        *b_LepGood_classification;   //!
   TBranch        *b_LepGood_detaIn;   //!
   TBranch        *b_LepGood_dphiIn;   //!
   TBranch        *b_LepGood_sigmaIetaIeta;   //!
   TBranch        *b_LepGood_sigmaIphiIphi;   //!
   TBranch        *b_LepGood_hcalOverEcal;   //!
   TBranch        *b_LepGood_correctedEcalEnergyError;   //!
   TBranch        *b_LepGood_regressionEnergyError;   //!
   TBranch        *b_LepGood_ecalEnergy;   //!
   TBranch        *b_LepGood_superCluster_rawEnergy;   //!
   TBranch        *b_LepGood_superCluster_preshowerEnergy;   //!
   TBranch        *b_LepGood_superCluster_energy;   //!
   TBranch        *b_LepGood_eSuperClusterOverP;   //!
   TBranch        *b_nJetFwd;   //!
   TBranch        *b_JetFwd_pt;   //!
   TBranch        *b_JetFwd_eta;   //!
   TBranch        *b_JetFwd_phi;   //!
   TBranch        *b_JetFwd_mass;   //!
   TBranch        *b_JetFwd_btagCSV;   //!
   TBranch        *b_JetFwd_rawPt;   //!
   TBranch        *b_JetFwd_mcPt;   //!
   TBranch        *b_JetFwd_mcFlavour;   //!
   TBranch        *b_JetFwd_quarkGluonID;   //!
   TBranch        *b_JetFwd_mcMatchId;   //!
   TBranch        *b_JetFwd_mcMatchFlav;   //!
   TBranch        *b_JetFwd_puId;   //!
   TBranch        *b_JetFwd_area;   //!
   TBranch        *b_JetFwd_id;   //!
   TBranch        *b_JetFwd_CHEF;   //!
   TBranch        *b_JetFwd_NHEF;   //!
   TBranch        *b_JetFwd_PHEF;   //!
   TBranch        *b_JetFwd_MUEF;   //!
   TBranch        *b_JetFwd_ELEF;   //!
   TBranch        *b_nJet;   //!
   TBranch        *b_Jet_pt;   //!
   TBranch        *b_Jet_eta;   //!
   TBranch        *b_Jet_phi;   //!
   TBranch        *b_Jet_mass;   //!
   TBranch        *b_Jet_btagCSV;   //!
   TBranch        *b_Jet_rawPt;   //!
   TBranch        *b_Jet_mcPt;   //!
   TBranch        *b_Jet_mcFlavour;   //!
   TBranch        *b_Jet_quarkGluonID;   //!
   TBranch        *b_Jet_mcMatchId;   //!
   TBranch        *b_Jet_mcMatchFlav;   //!
   TBranch        *b_Jet_puId;   //!
   TBranch        *b_Jet_area;   //!
   TBranch        *b_Jet_id;   //!
   TBranch        *b_Jet_CHEF;   //!
   TBranch        *b_Jet_NHEF;   //!
   TBranch        *b_Jet_PHEF;   //!
   TBranch        *b_Jet_MUEF;   //!
   TBranch        *b_Jet_ELEF;   //!
   TBranch        *b_nGenLep;   //!
   TBranch        *b_GenLep_pt;   //!
   TBranch        *b_GenLep_eta;   //!
   TBranch        *b_GenLep_phi;   //!
   TBranch        *b_GenLep_mass;   //!
   TBranch        *b_GenLep_pdgId;   //!
   TBranch        *b_GenLep_charge;   //!
   TBranch        *b_GenLep_sourceId;   //!
   TBranch        *b_nGenP6StatusThree;   //!
   TBranch        *b_GenP6StatusThree_pt;   //!
   TBranch        *b_GenP6StatusThree_eta;   //!
   TBranch        *b_GenP6StatusThree_phi;   //!
   TBranch        *b_GenP6StatusThree_mass;   //!
   TBranch        *b_GenP6StatusThree_pdgId;   //!
   TBranch        *b_GenP6StatusThree_charge;   //!
   TBranch        *b_GenP6StatusThree_motherId;   //!
   TBranch        *b_GenP6StatusThree_grandmaId;   //!

   treeProducerWMassEle(TTree *tree=0);
   virtual ~treeProducerWMassEle();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef treeProducerWMassEle_cxx
treeProducerWMassEle::treeProducerWMassEle(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/u2/emanuele/TREES_1LEP_53X_V1/DYJetsM50/treeProducerWMassEle/treeProducerWMassEle_tree.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/u2/emanuele/TREES_1LEP_53X_V1/DYJetsM50/treeProducerWMassEle/treeProducerWMassEle_tree.root");
      }
      f->GetObject("treeProducerWMassEle",tree);

   }
   Init(tree);
}

treeProducerWMassEle::~treeProducerWMassEle()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t treeProducerWMassEle::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t treeProducerWMassEle::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void treeProducerWMassEle::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("evt", &evt, &b_evt);
   fChain->SetBranchAddress("isData", &isData, &b_isData);
   fChain->SetBranchAddress("HLT_SingleMu", &HLT_SingleMu, &b_HLT_SingleMu);
   fChain->SetBranchAddress("HLT_SingleEl", &HLT_SingleEl, &b_HLT_SingleEl);
   fChain->SetBranchAddress("HLT_MuEG", &HLT_MuEG, &b_HLT_MuEG);
   fChain->SetBranchAddress("HLT_TripleEl", &HLT_TripleEl, &b_HLT_TripleEl);
   fChain->SetBranchAddress("HLT_DoubleEl", &HLT_DoubleEl, &b_HLT_DoubleEl);
   fChain->SetBranchAddress("HLT_DoubleMu", &HLT_DoubleMu, &b_HLT_DoubleMu);
   fChain->SetBranchAddress("puWeight", &puWeight, &b_puWeight);
   fChain->SetBranchAddress("pdfWeight_CT10", pdfWeight_CT10, &b_pdfWeight_CT10);
   fChain->SetBranchAddress("pdfWeight_MSTW2008lo68cl", pdfWeight_MSTW2008lo68cl, &b_pdfWeight_MSTW2008lo68cl);
   fChain->SetBranchAddress("pdfWeight_NNPDF21_100", pdfWeight_NNPDF21_100, &b_pdfWeight_NNPDF21_100);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("nVert", &nVert, &b_nVert);
   fChain->SetBranchAddress("GenHeaviestQCDFlavour", &GenHeaviestQCDFlavour, &b_GenHeaviestQCDFlavour);
   fChain->SetBranchAddress("LepEff_1lep", &LepEff_1lep, &b_LepEff_1lep);
   fChain->SetBranchAddress("LepEff_2lep", &LepEff_2lep, &b_LepEff_2lep);
   fChain->SetBranchAddress("mZ1", &mZ1, &b_mZ1);
   fChain->SetBranchAddress("m2l", &m2l, &b_m2l);
   fChain->SetBranchAddress("tkmet_pt", &tkmet_pt, &b_tkmet_pt);
   fChain->SetBranchAddress("tkmet_eta", &tkmet_eta, &b_tkmet_eta);
   fChain->SetBranchAddress("tkmet_phi", &tkmet_phi, &b_tkmet_phi);
   fChain->SetBranchAddress("tkmet_mass", &tkmet_mass, &b_tkmet_mass);
   fChain->SetBranchAddress("pucmet_pt", &pucmet_pt, &b_pucmet_pt);
   fChain->SetBranchAddress("pucmet_eta", &pucmet_eta, &b_pucmet_eta);
   fChain->SetBranchAddress("pucmet_phi", &pucmet_phi, &b_pucmet_phi);
   fChain->SetBranchAddress("pucmet_mass", &pucmet_mass, &b_pucmet_mass);
   fChain->SetBranchAddress("met_pt", &met_pt, &b_met_pt);
   fChain->SetBranchAddress("met_eta", &met_eta, &b_met_eta);
   fChain->SetBranchAddress("met_phi", &met_phi, &b_met_phi);
   fChain->SetBranchAddress("met_mass", &met_mass, &b_met_mass);
   fChain->SetBranchAddress("met_sumEt", &met_sumEt, &b_met_sumEt);
   fChain->SetBranchAddress("met_genPt", &met_genPt, &b_met_genPt);
   fChain->SetBranchAddress("met_genPhi", &met_genPhi, &b_met_genPhi);
   fChain->SetBranchAddress("met_genEta", &met_genEta, &b_met_genEta);
   fChain->SetBranchAddress("metraw_pt", &metraw_pt, &b_metraw_pt);
   fChain->SetBranchAddress("metraw_eta", &metraw_eta, &b_metraw_eta);
   fChain->SetBranchAddress("metraw_phi", &metraw_phi, &b_metraw_phi);
   fChain->SetBranchAddress("metraw_mass", &metraw_mass, &b_metraw_mass);
   fChain->SetBranchAddress("metraw_sumEt", &metraw_sumEt, &b_metraw_sumEt);
   fChain->SetBranchAddress("metraw_genPt", &metraw_genPt, &b_metraw_genPt);
   fChain->SetBranchAddress("metraw_genPhi", &metraw_genPhi, &b_metraw_genPhi);
   fChain->SetBranchAddress("metraw_genEta", &metraw_genEta, &b_metraw_genEta);
   fChain->SetBranchAddress("pumet_pt", &pumet_pt, &b_pumet_pt);
   fChain->SetBranchAddress("pumet_eta", &pumet_eta, &b_pumet_eta);
   fChain->SetBranchAddress("pumet_phi", &pumet_phi, &b_pumet_phi);
   fChain->SetBranchAddress("pumet_mass", &pumet_mass, &b_pumet_mass);
   fChain->SetBranchAddress("metNoPU_pt", &metNoPU_pt, &b_metNoPU_pt);
   fChain->SetBranchAddress("metNoPU_eta", &metNoPU_eta, &b_metNoPU_eta);
   fChain->SetBranchAddress("metNoPU_phi", &metNoPU_phi, &b_metNoPU_phi);
   fChain->SetBranchAddress("metNoPU_mass", &metNoPU_mass, &b_metNoPU_mass);
   fChain->SetBranchAddress("nLepOther", &nLepOther, &b_nLepOther);
   fChain->SetBranchAddress("LepOther_pt", LepOther_pt, &b_LepOther_pt);
   fChain->SetBranchAddress("LepOther_eta", LepOther_eta, &b_LepOther_eta);
   fChain->SetBranchAddress("LepOther_phi", LepOther_phi, &b_LepOther_phi);
   fChain->SetBranchAddress("LepOther_mass", LepOther_mass, &b_LepOther_mass);
   fChain->SetBranchAddress("LepOther_pdgId", LepOther_pdgId, &b_LepOther_pdgId);
   fChain->SetBranchAddress("LepOther_charge", LepOther_charge, &b_LepOther_charge);
   fChain->SetBranchAddress("LepOther_dxy", LepOther_dxy, &b_LepOther_dxy);
   fChain->SetBranchAddress("LepOther_dz", LepOther_dz, &b_LepOther_dz);
   fChain->SetBranchAddress("LepOther_edxy", LepOther_edxy, &b_LepOther_edxy);
   fChain->SetBranchAddress("LepOther_edz", LepOther_edz, &b_LepOther_edz);
   fChain->SetBranchAddress("LepOther_ip3d", LepOther_ip3d, &b_LepOther_ip3d);
   fChain->SetBranchAddress("LepOther_sip3d", LepOther_sip3d, &b_LepOther_sip3d);
   fChain->SetBranchAddress("LepOther_tightId", LepOther_tightId, &b_LepOther_tightId);
   fChain->SetBranchAddress("LepOther_convVeto", LepOther_convVeto, &b_LepOther_convVeto);
   fChain->SetBranchAddress("LepOther_lostHits", LepOther_lostHits, &b_LepOther_lostHits);
   fChain->SetBranchAddress("LepOther_looseIdSusy", LepOther_looseIdSusy, &b_LepOther_looseIdSusy);
   fChain->SetBranchAddress("LepOther_relIso03", LepOther_relIso03, &b_LepOther_relIso03);
   fChain->SetBranchAddress("LepOther_relIso04", LepOther_relIso04, &b_LepOther_relIso04);
   fChain->SetBranchAddress("LepOther_chargedHadRelIso03", LepOther_chargedHadRelIso03, &b_LepOther_chargedHadRelIso03);
   fChain->SetBranchAddress("LepOther_chargedHadRelIso04", LepOther_chargedHadRelIso04, &b_LepOther_chargedHadRelIso04);
   fChain->SetBranchAddress("LepOther_convVetoFull", LepOther_convVetoFull, &b_LepOther_convVetoFull);
   fChain->SetBranchAddress("LepOther_eleCutId", LepOther_eleCutId, &b_LepOther_eleCutId);
   fChain->SetBranchAddress("LepOther_eleMVAId", LepOther_eleMVAId, &b_LepOther_eleMVAId);
   fChain->SetBranchAddress("LepOther_tightCharge", LepOther_tightCharge, &b_LepOther_tightCharge);
   fChain->SetBranchAddress("LepOther_mvaId", LepOther_mvaId, &b_LepOther_mvaId);
   fChain->SetBranchAddress("LepOther_mvaIdTrig", LepOther_mvaIdTrig, &b_LepOther_mvaIdTrig);
   fChain->SetBranchAddress("LepOther_nStations", LepOther_nStations, &b_LepOther_nStations);
   fChain->SetBranchAddress("LepOther_trkKink", LepOther_trkKink, &b_LepOther_trkKink);
   fChain->SetBranchAddress("LepOther_caloCompatibility", LepOther_caloCompatibility, &b_LepOther_caloCompatibility);
   fChain->SetBranchAddress("LepOther_globalTrackChi2", LepOther_globalTrackChi2, &b_LepOther_globalTrackChi2);
   fChain->SetBranchAddress("LepOther_trackerLayers", LepOther_trackerLayers, &b_LepOther_trackerLayers);
   fChain->SetBranchAddress("LepOther_pixelLayers", LepOther_pixelLayers, &b_LepOther_pixelLayers);
   fChain->SetBranchAddress("LepOther_mvaTTH", LepOther_mvaTTH, &b_LepOther_mvaTTH);
   fChain->SetBranchAddress("LepOther_jetPtRatio", LepOther_jetPtRatio, &b_LepOther_jetPtRatio);
   fChain->SetBranchAddress("LepOther_jetBTagCSV", LepOther_jetBTagCSV, &b_LepOther_jetBTagCSV);
   fChain->SetBranchAddress("LepOther_jetDR", LepOther_jetDR, &b_LepOther_jetDR);
   fChain->SetBranchAddress("LepOther_mcMatchId", LepOther_mcMatchId, &b_LepOther_mcMatchId);
   fChain->SetBranchAddress("LepOther_mcMatchAny", LepOther_mcMatchAny, &b_LepOther_mcMatchAny);
   fChain->SetBranchAddress("LepOther_mcMatchAny2", LepOther_mcMatchAny2, &b_LepOther_mcMatchAny2);
   fChain->SetBranchAddress("LepOther_mcMatchTau", LepOther_mcMatchTau, &b_LepOther_mcMatchTau);
   fChain->SetBranchAddress("LepOther_softMuID", LepOther_softMuID, &b_LepOther_softMuID);
   fChain->SetBranchAddress("LepOther_trgMatch", LepOther_trgMatch, &b_LepOther_trgMatch);
   fChain->SetBranchAddress("LepOther_scEta", LepOther_scEta, &b_LepOther_scEta);
   fChain->SetBranchAddress("LepOther_r9", LepOther_r9, &b_LepOther_r9);
   fChain->SetBranchAddress("LepOther_classification", LepOther_classification, &b_LepOther_classification);
   fChain->SetBranchAddress("LepOther_detaIn", LepOther_detaIn, &b_LepOther_detaIn);
   fChain->SetBranchAddress("LepOther_dphiIn", LepOther_dphiIn, &b_LepOther_dphiIn);
   fChain->SetBranchAddress("LepOther_sigmaIetaIeta", LepOther_sigmaIetaIeta, &b_LepOther_sigmaIetaIeta);
   fChain->SetBranchAddress("LepOther_sigmaIphiIphi", LepOther_sigmaIphiIphi, &b_LepOther_sigmaIphiIphi);
   fChain->SetBranchAddress("LepOther_hcalOverEcal", LepOther_hcalOverEcal, &b_LepOther_hcalOverEcal);
   fChain->SetBranchAddress("LepOther_correctedEcalEnergyError", LepOther_correctedEcalEnergyError, &b_LepOther_correctedEcalEnergyError);
   fChain->SetBranchAddress("LepOther_regressionEnergyError", LepOther_regressionEnergyError, &b_LepOther_regressionEnergyError);
   fChain->SetBranchAddress("LepOther_ecalEnergy", LepOther_ecalEnergy, &b_LepOther_ecalEnergy);
   fChain->SetBranchAddress("LepOther_superCluster_rawEnergy", LepOther_superCluster_rawEnergy, &b_LepOther_superCluster_rawEnergy);
   fChain->SetBranchAddress("LepOther_superCluster_preshowerEnergy", LepOther_superCluster_preshowerEnergy, &b_LepOther_superCluster_preshowerEnergy);
   fChain->SetBranchAddress("LepOther_superCluster_energy", LepOther_superCluster_energy, &b_LepOther_superCluster_energy);
   fChain->SetBranchAddress("LepOther_eSuperClusterOverP", LepOther_eSuperClusterOverP, &b_LepOther_eSuperClusterOverP);
   fChain->SetBranchAddress("nLepGood", &nLepGood, &b_nLepGood);
   fChain->SetBranchAddress("LepGood_pt", LepGood_pt, &b_LepGood_pt);
   fChain->SetBranchAddress("LepGood_eta", LepGood_eta, &b_LepGood_eta);
   fChain->SetBranchAddress("LepGood_phi", LepGood_phi, &b_LepGood_phi);
   fChain->SetBranchAddress("LepGood_mass", LepGood_mass, &b_LepGood_mass);
   fChain->SetBranchAddress("LepGood_pdgId", LepGood_pdgId, &b_LepGood_pdgId);
   fChain->SetBranchAddress("LepGood_charge", LepGood_charge, &b_LepGood_charge);
   fChain->SetBranchAddress("LepGood_dxy", LepGood_dxy, &b_LepGood_dxy);
   fChain->SetBranchAddress("LepGood_dz", LepGood_dz, &b_LepGood_dz);
   fChain->SetBranchAddress("LepGood_edxy", LepGood_edxy, &b_LepGood_edxy);
   fChain->SetBranchAddress("LepGood_edz", LepGood_edz, &b_LepGood_edz);
   fChain->SetBranchAddress("LepGood_ip3d", LepGood_ip3d, &b_LepGood_ip3d);
   fChain->SetBranchAddress("LepGood_sip3d", LepGood_sip3d, &b_LepGood_sip3d);
   fChain->SetBranchAddress("LepGood_tightId", LepGood_tightId, &b_LepGood_tightId);
   fChain->SetBranchAddress("LepGood_convVeto", LepGood_convVeto, &b_LepGood_convVeto);
   fChain->SetBranchAddress("LepGood_lostHits", LepGood_lostHits, &b_LepGood_lostHits);
   fChain->SetBranchAddress("LepGood_looseIdSusy", LepGood_looseIdSusy, &b_LepGood_looseIdSusy);
   fChain->SetBranchAddress("LepGood_relIso03", LepGood_relIso03, &b_LepGood_relIso03);
   fChain->SetBranchAddress("LepGood_relIso04", LepGood_relIso04, &b_LepGood_relIso04);
   fChain->SetBranchAddress("LepGood_chargedHadRelIso03", LepGood_chargedHadRelIso03, &b_LepGood_chargedHadRelIso03);
   fChain->SetBranchAddress("LepGood_chargedHadRelIso04", LepGood_chargedHadRelIso04, &b_LepGood_chargedHadRelIso04);
   fChain->SetBranchAddress("LepGood_convVetoFull", LepGood_convVetoFull, &b_LepGood_convVetoFull);
   fChain->SetBranchAddress("LepGood_eleCutId", LepGood_eleCutId, &b_LepGood_eleCutId);
   fChain->SetBranchAddress("LepGood_eleMVAId", LepGood_eleMVAId, &b_LepGood_eleMVAId);
   fChain->SetBranchAddress("LepGood_tightCharge", LepGood_tightCharge, &b_LepGood_tightCharge);
   fChain->SetBranchAddress("LepGood_mvaId", LepGood_mvaId, &b_LepGood_mvaId);
   fChain->SetBranchAddress("LepGood_mvaIdTrig", LepGood_mvaIdTrig, &b_LepGood_mvaIdTrig);
   fChain->SetBranchAddress("LepGood_nStations", LepGood_nStations, &b_LepGood_nStations);
   fChain->SetBranchAddress("LepGood_trkKink", LepGood_trkKink, &b_LepGood_trkKink);
   fChain->SetBranchAddress("LepGood_caloCompatibility", LepGood_caloCompatibility, &b_LepGood_caloCompatibility);
   fChain->SetBranchAddress("LepGood_globalTrackChi2", LepGood_globalTrackChi2, &b_LepGood_globalTrackChi2);
   fChain->SetBranchAddress("LepGood_trackerLayers", LepGood_trackerLayers, &b_LepGood_trackerLayers);
   fChain->SetBranchAddress("LepGood_pixelLayers", LepGood_pixelLayers, &b_LepGood_pixelLayers);
   fChain->SetBranchAddress("LepGood_mvaTTH", LepGood_mvaTTH, &b_LepGood_mvaTTH);
   fChain->SetBranchAddress("LepGood_jetPtRatio", LepGood_jetPtRatio, &b_LepGood_jetPtRatio);
   fChain->SetBranchAddress("LepGood_jetBTagCSV", LepGood_jetBTagCSV, &b_LepGood_jetBTagCSV);
   fChain->SetBranchAddress("LepGood_jetDR", LepGood_jetDR, &b_LepGood_jetDR);
   fChain->SetBranchAddress("LepGood_mcMatchId", LepGood_mcMatchId, &b_LepGood_mcMatchId);
   fChain->SetBranchAddress("LepGood_mcMatchAny", LepGood_mcMatchAny, &b_LepGood_mcMatchAny);
   fChain->SetBranchAddress("LepGood_mcMatchAny2", LepGood_mcMatchAny2, &b_LepGood_mcMatchAny2);
   fChain->SetBranchAddress("LepGood_mcMatchTau", LepGood_mcMatchTau, &b_LepGood_mcMatchTau);
   fChain->SetBranchAddress("LepGood_softMuID", LepGood_softMuID, &b_LepGood_softMuID);
   fChain->SetBranchAddress("LepGood_trgMatch", LepGood_trgMatch, &b_LepGood_trgMatch);
   fChain->SetBranchAddress("LepGood_scEta", LepGood_scEta, &b_LepGood_scEta);
   fChain->SetBranchAddress("LepGood_r9", LepGood_r9, &b_LepGood_r9);
   fChain->SetBranchAddress("LepGood_classification", LepGood_classification, &b_LepGood_classification);
   fChain->SetBranchAddress("LepGood_detaIn", LepGood_detaIn, &b_LepGood_detaIn);
   fChain->SetBranchAddress("LepGood_dphiIn", LepGood_dphiIn, &b_LepGood_dphiIn);
   fChain->SetBranchAddress("LepGood_sigmaIetaIeta", LepGood_sigmaIetaIeta, &b_LepGood_sigmaIetaIeta);
   fChain->SetBranchAddress("LepGood_sigmaIphiIphi", LepGood_sigmaIphiIphi, &b_LepGood_sigmaIphiIphi);
   fChain->SetBranchAddress("LepGood_hcalOverEcal", LepGood_hcalOverEcal, &b_LepGood_hcalOverEcal);
   fChain->SetBranchAddress("LepGood_correctedEcalEnergyError", LepGood_correctedEcalEnergyError, &b_LepGood_correctedEcalEnergyError);
   fChain->SetBranchAddress("LepGood_regressionEnergyError", LepGood_regressionEnergyError, &b_LepGood_regressionEnergyError);
   fChain->SetBranchAddress("LepGood_ecalEnergy", LepGood_ecalEnergy, &b_LepGood_ecalEnergy);
   fChain->SetBranchAddress("LepGood_superCluster_rawEnergy", LepGood_superCluster_rawEnergy, &b_LepGood_superCluster_rawEnergy);
   fChain->SetBranchAddress("LepGood_superCluster_preshowerEnergy", LepGood_superCluster_preshowerEnergy, &b_LepGood_superCluster_preshowerEnergy);
   fChain->SetBranchAddress("LepGood_superCluster_energy", LepGood_superCluster_energy, &b_LepGood_superCluster_energy);
   fChain->SetBranchAddress("LepGood_eSuperClusterOverP", LepGood_eSuperClusterOverP, &b_LepGood_eSuperClusterOverP);
   fChain->SetBranchAddress("nJetFwd", &nJetFwd, &b_nJetFwd);
   fChain->SetBranchAddress("JetFwd_pt", JetFwd_pt, &b_JetFwd_pt);
   fChain->SetBranchAddress("JetFwd_eta", JetFwd_eta, &b_JetFwd_eta);
   fChain->SetBranchAddress("JetFwd_phi", JetFwd_phi, &b_JetFwd_phi);
   fChain->SetBranchAddress("JetFwd_mass", JetFwd_mass, &b_JetFwd_mass);
   fChain->SetBranchAddress("JetFwd_btagCSV", JetFwd_btagCSV, &b_JetFwd_btagCSV);
   fChain->SetBranchAddress("JetFwd_rawPt", JetFwd_rawPt, &b_JetFwd_rawPt);
   fChain->SetBranchAddress("JetFwd_mcPt", JetFwd_mcPt, &b_JetFwd_mcPt);
   fChain->SetBranchAddress("JetFwd_mcFlavour", JetFwd_mcFlavour, &b_JetFwd_mcFlavour);
   fChain->SetBranchAddress("JetFwd_quarkGluonID", JetFwd_quarkGluonID, &b_JetFwd_quarkGluonID);
   fChain->SetBranchAddress("JetFwd_mcMatchId", JetFwd_mcMatchId, &b_JetFwd_mcMatchId);
   fChain->SetBranchAddress("JetFwd_mcMatchFlav", JetFwd_mcMatchFlav, &b_JetFwd_mcMatchFlav);
   fChain->SetBranchAddress("JetFwd_puId", JetFwd_puId, &b_JetFwd_puId);
   fChain->SetBranchAddress("JetFwd_area", JetFwd_area, &b_JetFwd_area);
   fChain->SetBranchAddress("JetFwd_id", JetFwd_id, &b_JetFwd_id);
   fChain->SetBranchAddress("JetFwd_CHEF", JetFwd_CHEF, &b_JetFwd_CHEF);
   fChain->SetBranchAddress("JetFwd_NHEF", JetFwd_NHEF, &b_JetFwd_NHEF);
   fChain->SetBranchAddress("JetFwd_PHEF", JetFwd_PHEF, &b_JetFwd_PHEF);
   fChain->SetBranchAddress("JetFwd_MUEF", JetFwd_MUEF, &b_JetFwd_MUEF);
   fChain->SetBranchAddress("JetFwd_ELEF", JetFwd_ELEF, &b_JetFwd_ELEF);
   fChain->SetBranchAddress("nJet", &nJet, &b_nJet);
   fChain->SetBranchAddress("Jet_pt", Jet_pt, &b_Jet_pt);
   fChain->SetBranchAddress("Jet_eta", Jet_eta, &b_Jet_eta);
   fChain->SetBranchAddress("Jet_phi", Jet_phi, &b_Jet_phi);
   fChain->SetBranchAddress("Jet_mass", Jet_mass, &b_Jet_mass);
   fChain->SetBranchAddress("Jet_btagCSV", Jet_btagCSV, &b_Jet_btagCSV);
   fChain->SetBranchAddress("Jet_rawPt", Jet_rawPt, &b_Jet_rawPt);
   fChain->SetBranchAddress("Jet_mcPt", Jet_mcPt, &b_Jet_mcPt);
   fChain->SetBranchAddress("Jet_mcFlavour", Jet_mcFlavour, &b_Jet_mcFlavour);
   fChain->SetBranchAddress("Jet_quarkGluonID", Jet_quarkGluonID, &b_Jet_quarkGluonID);
   fChain->SetBranchAddress("Jet_mcMatchId", Jet_mcMatchId, &b_Jet_mcMatchId);
   fChain->SetBranchAddress("Jet_mcMatchFlav", Jet_mcMatchFlav, &b_Jet_mcMatchFlav);
   fChain->SetBranchAddress("Jet_puId", Jet_puId, &b_Jet_puId);
   fChain->SetBranchAddress("Jet_area", Jet_area, &b_Jet_area);
   fChain->SetBranchAddress("Jet_id", Jet_id, &b_Jet_id);
   fChain->SetBranchAddress("Jet_CHEF", Jet_CHEF, &b_Jet_CHEF);
   fChain->SetBranchAddress("Jet_NHEF", Jet_NHEF, &b_Jet_NHEF);
   fChain->SetBranchAddress("Jet_PHEF", Jet_PHEF, &b_Jet_PHEF);
   fChain->SetBranchAddress("Jet_MUEF", Jet_MUEF, &b_Jet_MUEF);
   fChain->SetBranchAddress("Jet_ELEF", Jet_ELEF, &b_Jet_ELEF);
   fChain->SetBranchAddress("nGenLep", &nGenLep, &b_nGenLep);
   fChain->SetBranchAddress("GenLep_pt", GenLep_pt, &b_GenLep_pt);
   fChain->SetBranchAddress("GenLep_eta", GenLep_eta, &b_GenLep_eta);
   fChain->SetBranchAddress("GenLep_phi", GenLep_phi, &b_GenLep_phi);
   fChain->SetBranchAddress("GenLep_mass", GenLep_mass, &b_GenLep_mass);
   fChain->SetBranchAddress("GenLep_pdgId", GenLep_pdgId, &b_GenLep_pdgId);
   fChain->SetBranchAddress("GenLep_charge", GenLep_charge, &b_GenLep_charge);
   fChain->SetBranchAddress("GenLep_sourceId", GenLep_sourceId, &b_GenLep_sourceId);
   fChain->SetBranchAddress("nGenP6StatusThree", &nGenP6StatusThree, &b_nGenP6StatusThree);
   fChain->SetBranchAddress("GenP6StatusThree_pt", GenP6StatusThree_pt, &b_GenP6StatusThree_pt);
   fChain->SetBranchAddress("GenP6StatusThree_eta", GenP6StatusThree_eta, &b_GenP6StatusThree_eta);
   fChain->SetBranchAddress("GenP6StatusThree_phi", GenP6StatusThree_phi, &b_GenP6StatusThree_phi);
   fChain->SetBranchAddress("GenP6StatusThree_mass", GenP6StatusThree_mass, &b_GenP6StatusThree_mass);
   fChain->SetBranchAddress("GenP6StatusThree_pdgId", GenP6StatusThree_pdgId, &b_GenP6StatusThree_pdgId);
   fChain->SetBranchAddress("GenP6StatusThree_charge", GenP6StatusThree_charge, &b_GenP6StatusThree_charge);
   fChain->SetBranchAddress("GenP6StatusThree_motherId", GenP6StatusThree_motherId, &b_GenP6StatusThree_motherId);
   fChain->SetBranchAddress("GenP6StatusThree_grandmaId", GenP6StatusThree_grandmaId, &b_GenP6StatusThree_grandmaId);
   Notify();
}

Bool_t treeProducerWMassEle::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void treeProducerWMassEle::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t treeProducerWMassEle::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef treeProducerWMassEle_cxx
