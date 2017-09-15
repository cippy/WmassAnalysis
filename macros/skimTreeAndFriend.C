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

// this macro is used to copy trees with skim, but can be used also to copy a tree (e.g., from a version produced in 53X to another in 80X)

bool fileExists(const string& fname) {
  
  std::ifstream infile(fname);
  return infile.good();

}


void doSkim(const string& path = "./",
	    const string& outpath = "./",
	    const string& sampleName = "",
	    const string& dirPattern = "/treeProducerWMassEle/",
	    const string& treeFileName = "treeProducerWMassEle_tree.root",
	    const string& treeName = "treeProducerWMassEle",
	    const string& evVarFriend_path = "",
	    const Int_t has_sfFriend = 0,
	    const string& sfFriend_path = "")
{


  string infileName         = path +                                     sampleName + dirPattern + treeFileName;
  string infileFriendName   = path + evVarFriend_path + "evVarFriend_" + sampleName +                  ".root";
  string infileSfFriendName = path + sfFriend_path    + "sfFriend_"    + sampleName +                  ".root";


  // create directories for output
  string createDirCommand = "mkdir -p " + outpath + sampleName + dirPattern; 
  system(createDirCommand.c_str());
  
  string outfileName         = outpath +                                     sampleName + dirPattern + treeFileName;
  string outfileFriendName   = outpath + evVarFriend_path + "evVarFriend_" + sampleName +                  ".root";
  string outfileSfFriendName = outpath + sfFriend_path    + "sfFriend_"    + sampleName +                  ".root";


  cout << "======================" << endl;
  cout << "======================" << endl;
  cout << sampleName << endl;
  cout << "----------------------" << endl;

  ////////////////////////////
  // base tree

  TFile* infile = new TFile(infileName.c_str());
  if (!infile || infile->IsZombie()) {
    std::cout << "Cannot open file " << infileName << std::endl;
    exit(EXIT_FAILURE);
  }
  TTree* intree = (TTree*)infile->Get(treeName.c_str());
  if (!intree || intree == NULL) {
    std::cout << "Error: tree not found in file " << infileName << std::endl;
    exit(EXIT_FAILURE);
  }
  Long64_t nentries = intree->GetEntries();

  //////////////////////////////
  // evVarFriend

  TFile* infileFriend = new TFile(infileFriendName.c_str());
  if (!infileFriend || infileFriend->IsZombie()) {
    std::cout << "Cannot open file " << infileFriendName << std::endl;
    exit(EXIT_FAILURE);
  }
  TTree* infriend = (TTree*)infileFriend->Get("mjvars/t");
  if (!infriend || infriend == NULL) {
    std::cout << "Error: tree not found in file " << infileFriend << std::endl;
    exit(EXIT_FAILURE);
  }
  Long64_t nentriesFriend = infriend->GetEntries();


  //////////////////////////////
  // sfFriend if any

  TFile* infileSfFriend = NULL;
  TTree* insffriend = NULL;
  Long64_t nentriesSfFriend = -1;

  if (has_sfFriend) {

    infileSfFriend = new TFile(infileSfFriendName.c_str());
    if (!infileSfFriend || infileSfFriend->IsZombie()) {
      std::cout << "Cannot open file " << infileSfFriendName << std::endl;
      exit(EXIT_FAILURE);
    }
    insffriend = (TTree*)infileSfFriend->Get("sf/t");
    if (!insffriend || insffriend == NULL) {
      std::cout << "Error: tree not found in file " << infileSfFriend << std::endl;
      exit(EXIT_FAILURE);
    }
    nentriesSfFriend = insffriend->GetEntries();

  }


  if (nentries != nentriesFriend) {
    
    std::cout << "Warning: input tree and evVar friend have different number of entries" << std::endl;
    std::cout << "End of programme" << std::endl;
    exit(EXIT_FAILURE);

  }

  if (has_sfFriend && nentries != nentriesSfFriend) {
    
    std::cout << "Warning: input tree and sf friend have different number of entries" << std::endl;
    std::cout << "End of programme" << std::endl;
    exit(EXIT_FAILURE);

  }


  intree->SetName(treeName.c_str());
  infriend->SetName("t");
  if (has_sfFriend) insffriend->SetName("t");

  TFile* outfile = new TFile(outfileName.c_str(), "RECREATE");
  if (!outfile || outfile->IsZombie()) {
    std::cout << "Cannot open file " << outfileName << std::endl;
    exit(EXIT_FAILURE);
  }
  TTree* outtree = intree->CopyTree(""); // dummy selection

  cout << "Entries before skim = " << nentries << endl;
  cout << "Entries after skim = " << outtree->GetEntries() << "\t(" << (Double_t) 100. * outtree->GetEntries()/nentries << " % efficiency)" << endl;

  outfile->Write();
  outfile->Close();

  //std::cout << "==== check ====" << std::endl;

  ////////////////////////////
  // evVarFriend

  cout << "Now skimming evVarFriend" << endl;
  TFile* outfileFriend = new TFile(outfileFriendName.c_str(), "RECREATE");
  if (!outfileFriend || outfileFriend->IsZombie()) {
    std::cout << "Cannot open file " << outfileFriendName << std::endl;
    exit(EXIT_FAILURE);
  }
  TDirectory *dir = outfileFriend->mkdir("mjvars");
  dir->cd();
  TTree* outfriend = infriend->CloneTree(0);

  for (Long64_t i=0; i<nentries; i++) {
    intree->GetEntry(i);
    if (1) {  // dummy selection
	infriend->GetEntry(i);
	outfriend->Fill();
    }
  }

  outfileFriend->Write();
  outfileFriend->Close();

  //////////////////////////////////
  // sfFriend

  TFile* outfileSfFriend = NULL;
  TDirectory *sfdir = NULL;
  TTree* outsffriend = NULL;

  if (has_sfFriend) {

    cout << "Now skimming sfFriend" << endl;
    outfileSfFriend = new TFile(outfileSfFriendName.c_str(), "RECREATE");
    if (!outfileSfFriend || outfileSfFriend->IsZombie()) {
      std::cout << "Cannot open file " << outfileSfFriendName << std::endl;
      exit(EXIT_FAILURE);
    }
    sfdir = outfileSfFriend->mkdir("sf");
    sfdir->cd();
    outsffriend = insffriend->CloneTree(0);
    
    for (Long64_t i=0; i<nentries; i++) {
      intree->GetEntry(i);
      if (1) { // dummy selection
	insffriend->GetEntry(i);
	outsffriend->Fill();
      }
    }
    
    outfileSfFriend->Write();
    outfileSfFriend->Close();
    
  }


  infile->Close();
  infileFriend->Close();
  if (has_sfFriend) infileSfFriend->Close();



}



//====================================================

// use as:
// root -l -b -q checkTreeAndFriend.C++
// you can store the output inside a file. 
// to use from terminal and change inputs:
// root -l -b -q  'skimTreeAndFriend.C++("/u2/emanuele/","TREES_MET_80X_V4/", "/", "friends_SR/", "friends/","/u2/mciprian/tests/skimmedTrees/)'

void skimTreeAndFriend(const string& inpath           = "/u2/emanuele/TREES_1LEP_53X_V2/",
		       const string& outpath          = "./",
		       const string& dirPattern       = "/treeProducerWMassEle/",
		       const string& treeFileName     = "treeProducerWMassEle_tree.root",
		       const string& treeName         = "treeProducerWMassEle",
		       const string& evVarFriend_path = "friends/",  
		       const Int_t   has_sfFriend     = 1,        // pass 0 if no sfFriend trees are used for MC (data will automatically use 0)		 
		       const string& sfFriend_path    = "friends/")
{

  // when using MET trees, directory for evVar friends can be 'friends_SR' or 'friends_VM', so we skim both together (base tree is the same)

  // dirPatternTag says if base tree are to be taken inside /treeProducerDarkMatterMonoJet/ folder or if their name has just 'treeProducerDarkMatterMonoJet'
  // e.g. , can have:
  // <path>/QCD/treeProducerDarkMatterMonoJet/tree.root
  // or
  // <path>/QCD_treeProducerDarkMatterMonoJet_tree.root

  cout << endl;
  cout << "Starting skim ..." << endl;
  cout << endl;

  // create directories
  string createDirCommand = "";
  if (outpath != "./") {
    createDirCommand = "mkdir -p " + outpath;
    system(createDirCommand.c_str());
  }
  createDirCommand = "mkdir -p " + outpath + evVarFriend_path;
  system(createDirCommand.c_str());
  createDirCommand = "mkdir -p " + outpath + sfFriend_path;
  system(createDirCommand.c_str());

  cout << endl;
  cout << "Tree location: " << inpath << endl;
  cout << "evVarFriend location: " << inpath + evVarFriend_path << endl;
  cout << "sfFriend location: " << inpath + sfFriend_path << endl;
  cout << endl;

  vector<string> sampleNameVector;

  /////////////
  // MC samples

  // sampleNameVector.push_back("DYJetsM50");
  // sampleNameVector.push_back("QCDMuPt15");
  // sampleNameVector.push_back("TTJets");
  // sampleNameVector.push_back("Tbarsch");
  // sampleNameVector.push_back("TbartW");
  // sampleNameVector.push_back("Tbartch");
  // sampleNameVector.push_back("Tsch");
  // sampleNameVector.push_back("TtW");
  // sampleNameVector.push_back("Ttch");
  // sampleNameVector.push_back("WJets");
  // sampleNameVector.push_back("WWJets");
  // sampleNameVector.push_back("WZJets");
   

  for (UInt_t i = 0; i < sampleNameVector.size(); i++) {

    // check if tree exists (assume that if it exists then the friends are ok, so I recommend using checkTreeAndFriend.C before using this code)
    string fname =  inpath + sampleNameVector[i] + dirPattern + treeFileName;
    if (fileExists(fname)) {
      doSkim(inpath, outpath, sampleNameVector[i], dirPattern, treeFileName, treeName, evVarFriend_path, has_sfFriend, sfFriend_path);
    } else {
      cout << "Following file does not exist. Skipping it" << endl;
      cout << fname << endl;
    }

  }

  ////////
  // Data samples

  vector<string> sampleNameDataVector;

  // /u2/emanuele/TREES_1LEP_53X_V2/                                                                                                                                        
  sampleNameDataVector.push_back("DoubleElectronAB");
  sampleNameDataVector.push_back("DoubleElectronC");
  sampleNameDataVector.push_back("DoubleElectronD");
  sampleNameDataVector.push_back("DoubleMuAB");
  sampleNameDataVector.push_back("DoubleMuC");
  sampleNameDataVector.push_back("DoubleMuD");
  // sampleNameDataVector.push_back("MuEGAB");                                              
  // sampleNameDataVector.push_back("MuEGC");                                                          
  // sampleNameDataVector.push_back("MuEGD");                                                                              
  sampleNameDataVector.push_back("SingleElectronAB");
  sampleNameDataVector.push_back("SingleElectronC");
  sampleNameDataVector.push_back("SingleElectronD");
  sampleNameDataVector.push_back("SingleMuAB");
  sampleNameDataVector.push_back("SingleMuC");
  sampleNameDataVector.push_back("SingleMuD");
  

  for (UInt_t i = 0; i < sampleNameDataVector.size(); i++) {

    string fname =  inpath + sampleNameDataVector[i] + dirPattern + treeFileName;
    if (fileExists(fname)) {
      doSkim(inpath, outpath, sampleNameDataVector[i], dirPattern, treeFileName, treeName, evVarFriend_path, 0, sfFriend_path);
    } else {
      cout << "Following file does not exist. Skipping it" << endl;
      cout << fname << endl;
    }

  }



}
