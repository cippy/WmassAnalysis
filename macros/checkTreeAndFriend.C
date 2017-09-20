#include <stdio.h>
#include <stdlib.h>
#include <cstdlib> //as stdlib.h
#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TList.h>

using namespace std;


bool fileExists(const string& fname) {
  
  std::ifstream infile(fname);
  return infile.good();

}
     
string getNewString(const string& strToModify, const string& match, const string& newMatch) {
  
  string prematch = "";
  string postmatch = "";
  
  size_t pos = strToModify.find(match);  // look for match in strToModify
  if (pos != string::npos) {
    
    prematch.assign(strToModify, 0, pos);   // copy strToModify up to last character before match
    postmatch.assign(strToModify,pos+match.size(),string::npos);  // copy strToModify from first character after match
    // cout << "prematch = " << prematch << endl;
    // cout << "postmatch = " << postmatch << endl;

  } else {

    //cout << match << " not found in " << strToModify << endl; 
    return strToModify;   // if no match return original strToModify

  }
  
  string finalstrToModify = prematch + newMatch + postmatch;
  //cout << "finalstrToModify = " << finalstrToModify << endl;
  return finalstrToModify;

}


Long64_t readTree(const string& fname, const string& sampleName, const string& whichTree = "base", const string& baseTreeName = "") {

  // which tree can be "base", "evVarFriend", "sfFriend"
  
  TFile* infile = new TFile(fname.c_str());
  if (!infile || infile->IsZombie()) {
    cout << "Missing " << whichTree << " file ==> " << sampleName << endl;
    return -1;
  } else {
    string treePath = "tree"; // dafault, overwritten below
    if (whichTree == "base" && baseTreeName != "") treePath = baseTreeName;
    if (whichTree == "evVarFriend") treePath = "mjvars/t";
    if (whichTree == "sfFriend") treePath = "sf/t";
    TTree* intree = (TTree*) infile->Get(treePath.c_str());
    if (!intree) {
      cout << "Missing " << whichTree << " tree ==> " << sampleName << endl;
      return -2;
    } else {
      return intree->GetEntries();
    }
  }

}



void doCheck(const string& path             = "/u2/emanuele/TREES_1LEP_53X_V2/", 
	     const string& sampleName       = "", 
	     const string& dirPattern       = "/treeProducerWMassEle/", 
	     const string& treeFileName     = "treeProducerWMassEle_tree.root", 
	     const string& treeName         = "treeProducerWMassEle", 
	     const string& evVarFriend_path = "friends/",
	     const Int_t   has_sfFriend     = 0,
	     const string& sfFriend_path    = "friends/") 
{

  // dirPatternTag says if base tree are to be taken inside /treeProducerDarkMatterMonoJet/ folder or if their name has just 'treeProducerDarkMatterMonoJet'
  // e.g. , can have:
  // <path>/QCD/treeProducerDarkMatterMonoJet/tree.root
  // or
  // <path>/QCD_treeProducerDarkMatterMonoJet_tree.root

  // evVarFriend_path and sfFriend_path tell say if friends are inside a folder (you pass the name) or at the same level of the base trees (default)

 
  // try to manage the fact that sometimes the tree is named "tree" and is inside tree.root file despite the fact that "treePatternTag = treeProducerWMassEle"
  // usually "treePatternTag = treeProducerWMassEle" --> "infileName = /.../treeProducerWMassEle_tree.root" and "treeName = tree"  
  // this happens with the skimmed sample
  // if (path.find("WSKIM_V") != string::npos || path.find("QCDSKIM_V") != string::npos || path.find("ZEESKIM_V") != string::npos) {
  //   dirPattern = dirPatternInit;
  //   treeName = "tree";
  // }


  string infileName         = path +                                     sampleName + dirPattern + treeFileName;  // cout << "infileName --> " << infileName << endl;
  string infileFriendName   = path + evVarFriend_path + "evVarFriend_" + sampleName +                  ".root";
  string infileSfFriendName = path + sfFriend_path    + "sfFriend_"    + sampleName +                  ".root";
  
  Long64_t nentries = -1;
  Long64_t nentriesFriend = -1;
  Long64_t nentriesSfFriend = -1;

  cout << "======================" << endl;
  cout << "======================" << endl;
  cout << sampleName << endl;
  cout << "----------------------" << endl;
  
  nentries                           = readTree(infileName,        sampleName,"base", treeName);
  // try to manage the fact that sometimes the tree is named "tree" and is inside tree.root file despite the fact that "treePatternTag = treeProducerWMassEle"
  // usually "treePatternTag = treeProducerWMassEle" --> "infileName = /.../treeProducerWMassEle_tree.root" and "treeName = tree"  
  // if (nentries == -1 && treeName == "treeProducerWMassEle") {
  //   cout << "\n#######\n";
  //   cout << "### WARNING: attempt to open file 'treeProducerWMassEle_tree.root' looking for TTree named 'treeProducerWMassEle' was unsuccessful.\n";
  //   cout << "### Trying again with 'infileName = tree.root' and 'treeName = tree' ...\n";
  //   // setting again parameters
  //   dirPattern = dirPatternInit;
  //   infileName = path + sampleName + dirPattern + "tree.root";
  //   treeName = "tree";
  //   // check again with new parameters
  //   cout << "#######\n" << endl;
  //   nentries = readTree(infileName, sampleName, "base", treeName);
  //   // if (nentries >= 0) cout << "File and tree read successfully!" << endl;
  // }

  nentriesFriend                     = readTree(infileFriendName,  sampleName,"evVarFriend");
  if (has_sfFriend) nentriesSfFriend = readTree(infileSfFriendName,sampleName,"sfFriend");

  if (nentries >= 0 && nentriesFriend >= 0 && nentries != nentriesFriend) {
    cout << sampleName << " ==> different nentries with evVarFriend" << endl;
  }
  if (has_sfFriend && nentries >= 0 && nentriesSfFriend >= 0 && nentries != nentriesSfFriend) {
    cout << sampleName << " ==> different nentries with sfFriend" << endl;
  }

  // sometimes the file exists but the name has HTXtY instead of HTXtoY ("to" is "t" in the name). 
  // if readTree returned -1, check for the presence of that file and use that for the following.
  // if not found, it is assumed the file name was correct and it just doesn't exist

  // passing name of sample, pattern to match and patten to substitutes the match with
  // function returns sampleName if the file with modified name doesn't exist
  // from now on in this function we will use the following string as sampleName
  string newSampleName = getNewString(sampleName, "to", "t");
 
  if (newSampleName != sampleName) {

    string newInfileName         = path +                                     newSampleName + dirPattern + treeFileName;
    string newInfileFriendName   = path + evVarFriend_path + "evVarFriend_" + newSampleName +                  ".root";
    string newInfileSfFriendName = path + sfFriend_path    + "sfFriend_"    + newSampleName +                  ".root";

    cout << endl;
    
    if ( nentries == -1 && fileExists(newInfileName) ) {

      cout << "Found base file with similar name ==> " << newSampleName << endl;  
      nentries = readTree(newInfileName, newSampleName,"base", treeName);

    }
    
    if ( nentriesFriend == -1 && fileExists(newInfileFriendName) ) {
      cout << "Found evVarFriend file with similar name ==> " << newSampleName << endl;  
      nentriesFriend = readTree(newInfileFriendName,newSampleName,"evVarFriend");
      if (nentries >= 0 && nentriesFriend >= 0 && nentries != nentriesFriend) {
	cout << newSampleName << " ==> different nentries with evVarFriend" << endl;
      }

    }
    
    if ( nentriesSfFriend == -1 && fileExists(newInfileSfFriendName) ) {

      cout << "Found sfFriend file with similar name ==> " << newSampleName << endl;  
      if (has_sfFriend) nentriesSfFriend = readTree(newInfileSfFriendName,newSampleName,"sfFriend");
      if (has_sfFriend && nentries >= 0 && nentriesSfFriend >= 0 && nentries != nentriesSfFriend) {
	cout << newSampleName << " ==> different nentries with sfFriend" << endl;
      }

    }
    
  }

  cout << endl;

}

//====================================================

// use as:
// root -l -b -q checkTreeAndFriend.C++
// you can store the output inside a file. 
// to use from terminal and change inputs:
// root -l -b -q  'checkTreeAndFriend.C++("/u2/emanuele/","TREES_MET_80X_V4/", "/", "treeProducerDarkMatterMonoJet", "friends_VM/", 1, "friends/")'

void checkTreeAndFriend(const string& path             = "/u2/emanuele/TREES_1LEP_53X_V2/",
			const string& dirPattern       = "/treeProducerWMassEle/",			
			const string& treeFileName     = "treeProducerWMassEle_tree.root",
			const string& treeName         = "treeProducerWMassEle",
			const string& evVarFriend_path = "friends/", 
			const Int_t   has_sfFriend     = 1,             // pass 0 if no sfFriend trees are used for MC (data will automatically use 0 in doCheck())
			const string& sfFriend_path    = "friends/")
{

  // dirPatternTag says if base tree are to be taken inside /treeProducerDarkMatterMonoJet/ folder or if their name has just 'treeProducerDarkMatterMonoJet'
  // e.g. , can have:
  // <path>/QCD/treeProducerDarkMatterMonoJet/tree.root
  // or
  // <path>/QCD_treeProducerDarkMatterMonoJet_tree.root

  // evVarFriend_path and sfFriend_path tell say if friends are inside a folder (you pass the name) or at the same level of the base trees (default)
  // evVarFriend_path can be "", "friends", "friends_SR", "friends_VM" or "friends_VE"
  // sfFriend_path can be "" or "friends"

  cout << endl;
  cout << "Tree location: " << path << endl;
  cout << "evVarFriend location: " << path + evVarFriend_path << endl;
  cout << "sfFriend location: " << path + sfFriend_path << endl;
  cout << "Tree pattern: " << dirPattern + treeFileName << endl;
  cout << endl;

  vector<string> sampleNameVector;

  /////////////
  // MC samples

  // // Zvv QCD
  // sampleNameVector.push_back("ZJetsToNuNu_HT100to200");
  // sampleNameVector.push_back("ZJetsToNuNu_HT100to200_ext");
  // sampleNameVector.push_back("ZJetsToNuNu_HT200to400");
  // sampleNameVector.push_back("ZJetsToNuNu_HT200to400_ext");
  // sampleNameVector.push_back("ZJetsToNuNu_HT400to600");
  // sampleNameVector.push_back("ZJetsToNuNu_HT400to600_ext");
  // sampleNameVector.push_back("ZJetsToNuNu_HT600to800");
  // sampleNameVector.push_back("ZJetsToNuNu_HT800to1200");
  // sampleNameVector.push_back("ZJetsToNuNu_HT800to1200_ext");
  // sampleNameVector.push_back("ZJetsToNuNu_HT1200to2500");
  // sampleNameVector.push_back("ZJetsToNuNu_HT1200to2500_ext");
  // sampleNameVector.push_back("ZJetsToNuNu_HT2500toInf");
  // sampleNameVector.push_back("ZJetsToNuNu_HT2500toInf_ext");
  // // Wlv QCD
  // sampleNameVector.push_back("WJetsToLNu_HT100to200");
  // sampleNameVector.push_back("WJetsToLNu_HT100to200_ext");
  // sampleNameVector.push_back("WJetsToLNu_HT200to400");
  // sampleNameVector.push_back("WJetsToLNu_HT200to400_ext");
  // sampleNameVector.push_back("WJetsToLNu_HT400to600");
  // sampleNameVector.push_back("WJetsToLNu_HT400to600_ext");
  // sampleNameVector.push_back("WJetsToLNu_HT600to800");
  // sampleNameVector.push_back("WJetsToLNu_HT800to1200");
  // sampleNameVector.push_back("WJetsToLNu_HT800to1200_ext");
  // sampleNameVector.push_back("WJetsToLNu_HT1200to2500");
  // sampleNameVector.push_back("WJetsToLNu_HT1200to2500_ext");
  // sampleNameVector.push_back("WJetsToLNu_HT2500toInf");
  // sampleNameVector.push_back("WJetsToLNu_HT2500toInf_ext");
  // // Top
  // sampleNameVector.push_back("TBar_tWch");
  // sampleNameVector.push_back("TTJets");
  // sampleNameVector.push_back("TToLeptons_sch");
  // sampleNameVector.push_back("T_tWch");
  // // Diboson
  // sampleNameVector.push_back("ZZ");
  // sampleNameVector.push_back("WW");
  // sampleNameVector.push_back("WZ");
  // // QCD
  // sampleNameVector.push_back("QCD_HT200to300_ext");
  // sampleNameVector.push_back("QCD_HT300to500");
  // sampleNameVector.push_back("QCD_HT300to500_ext");
  // sampleNameVector.push_back("QCD_HT500to700");
  // sampleNameVector.push_back("QCD_HT500to700_ext");
  // sampleNameVector.push_back("QCD_HT700to1000");
  // sampleNameVector.push_back("QCD_HT700to1000_ext");
  // sampleNameVector.push_back("QCD_HT1000to1500");
  // sampleNameVector.push_back("QCD_HT1000to1500_ext");
  // sampleNameVector.push_back("QCD_HT1500to2000");
  // sampleNameVector.push_back("QCD_HT1500to2000_ext");
  // sampleNameVector.push_back("QCD_HT2000toInf");
  // sampleNameVector.push_back("QCD_HT2000toInf_ext");
  // // Zll QCD
  // sampleNameVector.push_back("DYJetsToLL_M50_HT100to200");
  // sampleNameVector.push_back("DYJetsToLL_M50_HT100to200_ext");
  // sampleNameVector.push_back("DYJetsToLL_M50_HT200to400");
  // sampleNameVector.push_back("DYJetsToLL_M50_HT200to400_ext");
  // sampleNameVector.push_back("DYJetsToLL_M50_HT400to600");
  // sampleNameVector.push_back("DYJetsToLL_M50_HT400to600_ext");
  // sampleNameVector.push_back("DYJetsToLL_M50_HT600toInf");
  // sampleNameVector.push_back("DYJetsToLL_M50_HT600toInf_ext");
  // // Gamma
  // sampleNameVector.push_back("GJets_HT40to100");
  // sampleNameVector.push_back("GJets_HT100to200");
  // sampleNameVector.push_back("GJets_HT400to600");
  // sampleNameVector.push_back("GJets_HT600toInf");
  // // Zvv EWK
  // sampleNameVector.push_back("EWKZToNuNu2Jets");
  // // Wlv EWK
  // sampleNameVector.push_back("EWKWMinus2Jets");
  // sampleNameVector.push_back("EWKWPlus2Jets");
  // // Zll EWK
  // sampleNameVector.push_back("EWKZToLL2Jets");
  // // signals only for SR or VM (not VE)
  // if (evVarFriend_path.find("SR") != string::npos || evVarFriend_path.find("VM") != string::npos) {
  //   // GGH signal
  //   sampleNameVector.push_back("GluGlu_HToInvisible_M110");
  //   sampleNameVector.push_back("GluGlu_HToInvisible_M125");
  //   sampleNameVector.push_back("GluGlu_HToInvisible_M150");
  //   sampleNameVector.push_back("GluGlu_HToInvisible_M200");
  //   sampleNameVector.push_back("GluGlu_HToInvisible_M300");
  //   sampleNameVector.push_back("GluGlu_HToInvisible_M400");
  //   sampleNameVector.push_back("GluGlu_HToInvisible_M500");
  //   sampleNameVector.push_back("GluGlu_HToInvisible_M600");
  //   sampleNameVector.push_back("GluGlu_HToInvisible_M800");
  //   // VBF H signal
  //   sampleNameVector.push_back("VBF_HToInvisible_M110");
  //   sampleNameVector.push_back("VBF_HToInvisible_M125");
  //   sampleNameVector.push_back("VBF_HToInvisible_M150");
  //   sampleNameVector.push_back("VBF_HToInvisible_M200");
  //   sampleNameVector.push_back("VBF_HToInvisible_M300");
  //   sampleNameVector.push_back("VBF_HToInvisible_M400");
  //   sampleNameVector.push_back("VBF_HToInvisible_M500");
  //   sampleNameVector.push_back("VBF_HToInvisible_M600");
  //   sampleNameVector.push_back("VBF_HToInvisible_M800");
  // }

  // sampleNameVector.push_back("DYJetsToLL_M50_reHLT");
  // sampleNameVector.push_back("WJetsToLNu_reHLT");
  // sampleNameVector.push_back("QCD_Pt1000toInf_Mu5");
  // sampleNameVector.push_back("QCD_Pt120to170_EMEnriched");
  // sampleNameVector.push_back("QCD_Pt120to170_Mu5");
  // sampleNameVector.push_back("QCD_Pt15to20_Mu5");
  // sampleNameVector.push_back("QCD_Pt170to300_EMEnriched");
  // sampleNameVector.push_back("QCD_Pt170to300_Mu5");
  // sampleNameVector.push_back("QCD_Pt20to30_EMEnriched");
  // sampleNameVector.push_back("QCD_Pt20to30_Mu5");
  // sampleNameVector.push_back("QCD_Pt300to470_Mu5");
  // sampleNameVector.push_back("QCD_Pt300to470_Mu5_ext");
  // sampleNameVector.push_back("QCD_Pt300toInf_EMEnriched");
  // sampleNameVector.push_back("QCD_Pt30to50_EMEnriched");
  // sampleNameVector.push_back("QCD_Pt30to50_Mu5");
  // sampleNameVector.push_back("QCD_Pt470to600_Mu5");
  // sampleNameVector.push_back("QCD_Pt470to600_Mu5_ext");
  // sampleNameVector.push_back("QCD_Pt50to80_EMEnriched");
  // sampleNameVector.push_back("QCD_Pt50to80_Mu5");
  // sampleNameVector.push_back("QCD_Pt600to800_Mu5");
  // sampleNameVector.push_back("QCD_Pt600to800_Mu5_ext");
  // sampleNameVector.push_back("QCD_Pt800to1000_Mu5");
  // sampleNameVector.push_back("QCD_Pt800to1000_Mu5_ext");
  // sampleNameVector.push_back("QCD_Pt80to120_EMEnriched");
  // sampleNameVector.push_back("QCD_Pt80to120_Mu5");
  // sampleNameVector.push_back("QCD_Pt_170to250_bcToE");
  // sampleNameVector.push_back("QCD_Pt_250toInf_bcToE");
  // sampleNameVector.push_back("QCD_Pt_30to80_bcToE");

  // /u2/emanuele/TREES_1LEP_53X_V2/
  sampleNameVector.push_back("DYJetsM50");
  sampleNameVector.push_back("QCDMuPt15");
  sampleNameVector.push_back("TTJets");
  sampleNameVector.push_back("Tbarsch");
  sampleNameVector.push_back("TbartW");
  sampleNameVector.push_back("Tbartch");
  sampleNameVector.push_back("Tsch");
  sampleNameVector.push_back("TtW");
  sampleNameVector.push_back("Ttch");
  sampleNameVector.push_back("WJets");
  sampleNameVector.push_back("WWJets");
  sampleNameVector.push_back("WZJets");


  for (UInt_t i = 0; i < sampleNameVector.size(); i++) {
    doCheck(path, sampleNameVector[i], dirPattern, treeFileName, treeName, evVarFriend_path, has_sfFriend, sfFriend_path);
  }

  ////////
  // Data samples

  vector<string> sampleNameDataVector;

  // if (evVarFriend_path.find("SR") != string::npos || evVarFriend_path.find("VM") != string::npos) {
  //   sampleNameDataVector.push_back("MET_Run2016B_PromptReco_v1_runs_272023_273146");
  //   sampleNameDataVector.push_back("MET_Run2016B_PromptReco_v2_runs_273150_275376");
  //   sampleNameDataVector.push_back("MET_Run2016C_PromptReco_v2_runs_275420_276283");
  //   sampleNameDataVector.push_back("MET_Run2016D_PromptReco_v2_runs_276315_276811");
  //   sampleNameDataVector.push_back("MET_Run2016E_PromptReco_v2_runs_276830_277420");
  //   sampleNameDataVector.push_back("MET_Run2016G_PromptReco_v1_runs_278817_279931");
  // } else if (evVarFriend_path.find("VE") != string::npos) {
  //   sampleNameDataVector.push_back("SingleElectron_Run2016B_PromptReco_v1_runs_272023_273146");
  //   sampleNameDataVector.push_back("SingleElectron_Run2016B_PromptReco_v2_runs_273150_275376");
  //   sampleNameDataVector.push_back("SingleElectron_Run2016C_PromptReco_v2_runs_275420_276283");
  //   sampleNameDataVector.push_back("SingleElectron_Run2016D_PromptReco_v2_runs_276315_276811");
  //   sampleNameDataVector.push_back("SingleElectron_Run2016E_PromptReco_v2_runs_276830_277420");
  //   sampleNameDataVector.push_back("SingleElectron_Run2016F_PromptReco_v1_runs_277820_278808");
  //   sampleNameDataVector.push_back("SingleElectron_Run2016G_PromptReco_v1_runs_278817_279931");
  // }

   // sampleNameDataVector.push_back("DoubleEG_Run2016B_23Sep2016_v3_runs_273150_275376");
   // sampleNameDataVector.push_back("DoubleEG_Run2016C_23Sep2016_v1_runs_271036_284044");
   // sampleNameDataVector.push_back("DoubleEG_Run2016D_23Sep2016_v1_runs_271036_284044");
   // sampleNameDataVector.push_back("DoubleEG_Run2016E_23Sep2016_v1_runs_271036_284044");
   // sampleNameDataVector.push_back("DoubleEG_Run2016F_23Sep2016_v1_runs_271036_284044");
   // sampleNameDataVector.push_back("DoubleEG_Run2016G_23Sep2016_v1_runs_271036_284044");
   // sampleNameDataVector.push_back("DoubleEG_Run2016H_PromptReco_v1_runs_281085_281201");
   // sampleNameDataVector.push_back("DoubleEG_Run2016H_PromptReco_v2_runs_281207_284035");
   // sampleNameDataVector.push_back("DoubleEG_Run2016H_PromptReco_v3_runs_284036_284044");

   // /u2/emanuele/TREES_1LEP_53X_V2/
   sampleNameDataVector.push_back("DoubleElectronAB");
   sampleNameDataVector.push_back("DoubleElectronC");
   sampleNameDataVector.push_back("DoubleElectronD");
   sampleNameDataVector.push_back("DoubleMuAB");
   sampleNameDataVector.push_back("DoubleMuC");
   sampleNameDataVector.push_back("DoubleMuD");
   sampleNameDataVector.push_back("MuEGAB");
   sampleNameDataVector.push_back("MuEGC");
   sampleNameDataVector.push_back("MuEGD");
   sampleNameDataVector.push_back("SingleElectronAB");
   sampleNameDataVector.push_back("SingleElectronC");
   sampleNameDataVector.push_back("SingleElectronD");
   sampleNameDataVector.push_back("SingleMuAB");
   sampleNameDataVector.push_back("SingleMuC");
   sampleNameDataVector.push_back("SingleMuD");


  for (UInt_t i = 0; i < sampleNameDataVector.size(); i++) {
    doCheck(path, sampleNameDataVector[i], dirPattern, treeFileName, treeName, evVarFriend_path, 0, sfFriend_path);
  }



}
