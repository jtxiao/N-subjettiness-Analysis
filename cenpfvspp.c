#include "HiForestAnalysis/hiForest.h"
#include <TH1D.h>
#include <TFile.h>
#include <TCanvas.h>
#include "fastjet/JetDefinition.hh"
#include <TPad.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TProfile.h>
#include <iostream>
#include <math.h>
// #include "fastjet/contrib/SoftDrop.hh"
#include "softkiller.h"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/Njettiness.hh"
#include "fastjet/contrib/NjettinessPlugin.hh"
#include "fastjet/contrib/WinnerTakeAllRecombiner.hh"
#include <fstream>
#include "TLorentzVector.h"
#include <vector>
#include <string>
#ifdef __MAKECINT__
#pragma link C++ class vector<TLorentzVector>;
#endif

using namespace std;
using namespace fastjet;
using namespace fastjet::contrib;

class MyUserInfo: public PseudoJet::UserInfoBase{
public:
  MyUserInfo(const int & pdg_id_in) :
  _pdg_id(pdg_id_in){}

  int pdg_id() const { return _pdg_id;}
  
protected:
  int _pdg_id;        
};

void cenpfvspp()
{
  TH1::SetDefaultSumw2();
  // HiForest * c = new HiForest("/mnt/hadoop/cms/store/user/dgulhan/Private_PYTHIA_2p76TeV_Track9Jet29_merged/HiForest_Private_PYTHIA_pthat80_track8_jet29_merged_forest_0.root","forest",cPP, 0);
  HiForest * c = new HiForest("/mnt/hadoop/cms/store/user/dgulhan/pp2013/data/prod24/hiForest_merged/HiForest_pp_Jet80_v8_PP2013_HiForest_PromptReco_JsonPP_Jet80_PPReco_merged_forest_0.root","forest",cPP, 0);

  // c->LoadNoTrees();
  // c->hasTrackTree = true;
  // c->hasSkimTree = true;
  // c->hasEvtTree = true;
  // c->hasPFTree = true;
  // c->hasTowerTree = true;
  
  // TFile * outputfile = new TFile("taupythia8.root", "recreate");
  // TFile * outputfile = new TFile("taupp8.root", "recreate");
  // TFile * outputfile = new TFile("ptpythia2.root", "recreate");
  TFile * outputfile = new TFile("tautest.root", "recreate");
  // TFile * outputfile = new TFile("ptpp2.root", "recreate");
  
  TH1D * canpt = new TH1D("canpt","Pt spectra of PF candidates; Pt; N", 35, 0.0001, 249.99999);
  TH1D * trkcanpt = new TH1D("trkcanpt","Pt spectra of track PF candidates; Pt; N", 35, 0.0001, 249.99999);
  TH1D * ECALcanpt = new TH1D("ECALcanpt","Pt spectra of ECAL PF candidates; Pt; N", 35, 0.0001, 249.99999);
  TH1D * HCALcanpt = new TH1D("HCALcanpt","Pt spectra of HCAL PF candidates; Pt; N", 35, 0.0001, 249.99999);
  
  TH1D * jetpt = new TH1D("jetpt","Pt spectra of jets; Pt; N", 35, 0.0001, 299.99999);
  
  TH1D * jetFF = new TH1D("jetFF","Fragmentation Functions of jets; Pt; N", 20, 0.0001, 4.99999);
  TH1D * jetFFc1 = new TH1D("jetFFc1","Fragmentation Functions of jets; Pt; N", 20, 0.0001, 4.99999);
  TH1D * jetFFc2 = new TH1D("jetFFc2","Fragmentation Functions of jets; Pt; N", 20, 0.0001, 4.99999);
  TH1D * jetFFc3 = new TH1D("jetFFc3","Fragmentation Functions of jets; Pt; N", 20, 0.0001, 4.99999);
  TH1D * jetFFc4 = new TH1D("jetFFc4","Fragmentation Functions of jets; Pt; N", 20, 0.0001, 4.99999);
  
  TH1D * tauonec1 = new TH1D("tauonec1","N-subjettiness Distribution;;N", 25, 0.0001, 64.9999);
  TH1D * tauonec2 = new TH1D("tauonec2",";;", 25, 0.0001, 64.9999);
  TH1D * tauonec3 = new TH1D("tauonec3",";", 25, 0.0001, 64.9999);
  TH1D * tauonec4 = new TH1D("tauonec4",";Tau One;", 25, 0.0001, 64.9999);
  TH1D * tautwoc1 = new TH1D("tautwoc1","N-subjettiness Distribution;;N", 25, 0.0001, 64.9999);
  TH1D * tautwoc2 = new TH1D("tautwoc2",";;N", 25, 0.0001, 64.9999);
  TH1D * tautwoc3 = new TH1D("tautwoc3",";;N", 25, 0.0001, 64.9999);
  TH1D * tautwoc4 = new TH1D("tautwoc4",";Tau Two;N", 25, 0.0001, 64.9999);
  TH1D * tauthreec1 = new TH1D("tauthreec1","N-subjettiness Distribution;;N", 25, 0.0001, 64.9999);
  TH1D * tauthreec2 = new TH1D("tauthreec2",";;N", 25, 0.0001, 64.9999);
  TH1D * tauthreec3 = new TH1D("tauthreec3",";;N", 25, 0.0001, 64.9999);
  TH1D * tauthreec4 = new TH1D("tauthreec4",";Tau Three;N", 25, 0.0001, 64.9999);
  TH1D * tautworc1 = new TH1D("tautworc1","N-subjettiness Distribution;;N", 25, 0.0001, 1.4999);
  TH1D * tautworc2 = new TH1D("tautworc2",";;N", 25, 0.0001, 1.4999);
  TH1D * tautworc3 = new TH1D("tautworc3",";;N", 25, 0.0001, 1.4999);
  TH1D * tautworc4 = new TH1D("tautworc4",";Tau Two/Tau One;N", 25, 0.0001, 1.4999);
  TH1D * tauthreerc1 = new TH1D("tauthreerc1","N-subjettiness Distribution;;N", 25, 0.0001, 1.4999);
  TH1D * tauthreerc2 = new TH1D("tauthreerc2",";;N", 25, 0.0001, 1.4999);
  TH1D * tauthreerc3 = new TH1D("tauthreerc3",";;N", 25, 0.0001, 1.4999);
  TH1D * tauthreerc4 = new TH1D("tauthreerc4",";Tau Three/Tau Two;N", 25, 0.0001, 1.4999);
  
  TH1D * tauoneb2c1 = new TH1D("tauoneb2c1","N-subjettiness Distribution;;N", 25, 0.0001, 64.9999);
  TH1D * tauoneb2c2 = new TH1D("tauoneb2c2",";;N", 25, 0.0001, 64.9999);
  TH1D * tauoneb2c3 = new TH1D("tauoneb2c3",";;N", 25, 0.0001, 64.9999);
  TH1D * tauoneb2c4 = new TH1D("tauoneb2c4",";Tau One;N", 25, 0.0001, 64.9999);
  TH1D * tautwob2c1 = new TH1D("tautwob2c1","N-subjettiness Distribution;;N", 25, 0.0001, 64.9999);
  TH1D * tautwob2c2 = new TH1D("tautwob2c2",";;N", 25, 0.0001, 64.9999);
  TH1D * tautwob2c3 = new TH1D("tautwob2c3",";;N", 25, 0.0001, 64.9999);
  TH1D * tautwob2c4 = new TH1D("tautwob2c4",";Tau Two;N", 25, 0.0001, 64.9999);
  TH1D * tauthreeb2c1 = new TH1D("tauthreeb2c1","N-subjettiness Distribution;;N", 25, 0.0001, 64.9999);
  TH1D * tauthreeb2c2 = new TH1D("tauthreeb2c2",";;N", 25, 0.0001, 64.9999);
  TH1D * tauthreeb2c3 = new TH1D("tauthreeb2c3",";;N", 25, 0.0001, 64.9999);
  TH1D * tauthreeb2c4 = new TH1D("tauthreeb2c4",";Tau Three;N", 25, 0.0001, 64.9999);
  TH1D * tautworb2c1 = new TH1D("tautworb2c1","N-subjettiness Distribution;;N", 25, 0.0001, 1.4999);
  TH1D * tautworb2c2 = new TH1D("tautworb2c2",";;N", 25, 0.0001, 1.4999);
  TH1D * tautworb2c3 = new TH1D("tautworb2c3",";;N", 25, 0.0001, 1.4999);
  TH1D * tautworb2c4 = new TH1D("tautworb2c4",";Tau Two/Tau One;N", 25, 0.0001, 1.4999);
  TH1D * tauthreerb2c1 = new TH1D("tauthreerb2c1","N-subjettiness Distribution;;N", 25, 0.0001, 1.4999);
  TH1D * tauthreerb2c2 = new TH1D("tauthreerb2c2",";/;N", 25, 0.0001, 1.4999);
  TH1D * tauthreerb2c3 = new TH1D("tauthreerb2c3",";/;N", 25, 0.0001, 1.4999);
  TH1D * tauthreerb2c4 = new TH1D("tauthreerbc42",";Tau Three/Tau Two;N", 25, 0.0001, 1.4999);

  TH1D * tauoneb3c1 = new TH1D("tauoneb3c1","N-subjettiness Distribution;;N", 25, 0.0001, 64.9999);
  TH1D * tauoneb3c2 = new TH1D("tauoneb3c2",";;N", 25, 0.0001, 64.9999);
  TH1D * tauoneb3c3 = new TH1D("tauoneb3c3",";;N", 25, 0.0001, 64.9999);
  TH1D * tauoneb3c4 = new TH1D("tauoneb3c4",";Tau One;N", 25, 0.0001, 64.9999);
  TH1D * tautwob3c1 = new TH1D("tautwob3c1","N-subjettiness Distribution;;N", 25, 0.0001, 64.9999);
  TH1D * tautwob3c2 = new TH1D("tautwob3c2",";;N", 25, 0.0001, 64.9999);
  TH1D * tautwob3c3 = new TH1D("tautwob3c3",";;N", 25, 0.0001, 64.9999);
  TH1D * tautwob3c4 = new TH1D("tautwob3c4",";Tau Two;N", 25, 0.0001, 64.9999);
  TH1D * tauthreeb3c1 = new TH1D("tauthreeb3c1","N-subjettiness Distribution;;N", 25, 0.0001, 64.9999);
  TH1D * tauthreeb3c2 = new TH1D("tauthreeb3c2",";;N", 25, 0.0001, 64.9999);
  TH1D * tauthreeb3c3 = new TH1D("tauthreeb3c3",";;N", 25, 0.0001, 64.9999);
  TH1D * tauthreeb3c4 = new TH1D("tauthreeb3c4",";Tau Three;N", 25, 0.0001, 64.9999);
  TH1D * tautworb3c1 = new TH1D("tautworb3c1","N-subjettiness Distribution;;N", 25, 0.0001, 1.4999);
  TH1D * tautworb3c2 = new TH1D("tautworb3c2",";;N", 25, 0.0001, 1.4999);
  TH1D * tautworb3c3 = new TH1D("tautworb3c3",";;N", 25, 0.0001, 1.4999);
  TH1D * tautworb3c4 = new TH1D("tautworb3c4",";Tau Two/Tau one;N", 25, 0.0001, 1.4999);
  TH1D * tauthreerb3c1 = new TH1D("tauthreerb3c1","N-subjettiness Distribution;;N", 25, 0.0001, 1.4999);
  TH1D * tauthreerb3c2 = new TH1D("tauthreerb3c2",";;N", 25, 0.0001, 1.4999);
  TH1D * tauthreerb3c3 = new TH1D("tauthreerb3c3",";;N", 25, 0.0001, 1.4999);
  TH1D * tauthreerb3c4 = new TH1D("tauthreerb3c4",";Tau Three/Tau Two;N", 25, 0.0001, 1.4999);
  
  TH1D * tauoneb4c1 = new TH1D("tauoneb4c1","N-subjettiness Distribution;;N", 25, 0.0001, 64.9999);
  TH1D * tauoneb4c2 = new TH1D("tauoneb4c2",";;N", 25, 0.0001, 64.9999);
  TH1D * tauoneb4c3 = new TH1D("tauoneb4c3",";;N", 25, 0.0001, 64.9999);
  TH1D * tauoneb4c4 = new TH1D("tauoneb4c4",";Tau One;N", 25, 0.0001, 64.9999);
  TH1D * tautwob4c1 = new TH1D("tautwob4c1","N-subjettiness Distribution;;N", 25, 0.0001, 64.9999);
  TH1D * tautwob4c2 = new TH1D("tautwob4c2",";;N", 25, 0.0001, 64.9999);
  TH1D * tautwob4c3 = new TH1D("tautwob4c3",";;N", 25, 0.0001, 64.9999);
  TH1D * tautwob4c4 = new TH1D("tautwob4c4",";Tau Two;N", 25, 0.0001, 64.9999);
  TH1D * tauthreeb4c1 = new TH1D("tauthreeb4c1","N-subjettiness Distribution;;N", 25, 0.0001, 64.9999);
  TH1D * tauthreeb4c2 = new TH1D("tauthreeb4c2",";;N", 25, 0.0001, 64.9999);
  TH1D * tauthreeb4c3 = new TH1D("tauthreeb4c3",";;N", 25, 0.0001, 64.9999);
  TH1D * tauthreeb4c4 = new TH1D("tauthreeb4c4",";Tau Three;N", 25, 0.0001, 64.9999);
  TH1D * tautworb4c1 = new TH1D("tautworb4c1","N-subjettiness Distribution;;N", 25, 0.0001, 1.4999);
  TH1D * tautworb4c2 = new TH1D("tautworb4c2",";;N", 25, 0.0001, 1.4999);
  TH1D * tautworb4c3 = new TH1D("tautworb4c3",";;N", 25, 0.0001, 1.4999);
  TH1D * tautworb4c4 = new TH1D("tautworb4c4",";Tau Two/Tau One;N", 25, 0.0001, 1.4999);
  TH1D * tauthreerb4c1 = new TH1D("tauthreerb4c1","N-subjettiness Distribution;;N", 25, 0.0001, 1.4999);
  TH1D * tauthreerb4c2 = new TH1D("tauthreerb4c2",";;N", 25, 0.0001, 1.4999);
  TH1D * tauthreerb4c3 = new TH1D("tauthreerb4c3",";;N", 25, 0.0001, 1.4999);
  TH1D * tauthreerb4c4 = new TH1D("tauthreerb4c4",";Tau Three/Tau Two;N", 25, 0.0001, 1.4999);
  
  TH1D * tauoneb5c1 = new TH1D("tauoneb5c1","N-subjettiness Distribution;;N", 25, 0.0001, 64.9999);
  TH1D * tauoneb5c2 = new TH1D("tauoneb5c2",";;N", 25, 0.0001, 64.9999);
  TH1D * tauoneb5c3 = new TH1D("tauoneb5c3",";;N", 25, 0.0001, 64.9999);
  TH1D * tauoneb5c4 = new TH1D("tauoneb5c4",";Tau One;N", 25, 0.0001, 64.9999);
  TH1D * tautwob5c1 = new TH1D("tautwob5c1","N-subjettiness Distribution;;N", 25, 0.0001, 64.9999);
  TH1D * tautwob5c2 = new TH1D("tautwob5c2",";;N", 25, 0.0001, 64.9999);
  TH1D * tautwob5c3 = new TH1D("tautwob5c3",";;N", 25, 0.0001, 64.9999);
  TH1D * tautwob5c4 = new TH1D("tautwob5c4",";Tau Two;N", 25, 0.0001, 64.9999);
  TH1D * tauthreeb5c1 = new TH1D("tauthreeb5c1","N-subjettiness Distribution;;N", 25, 0.0001, 64.9999);
  TH1D * tauthreeb5c2 = new TH1D("tauthreeb5c2",";;N", 25, 0.0001, 64.9999);
  TH1D * tauthreeb5c3 = new TH1D("tauthreeb5c3",";;N", 25, 0.0001, 64.9999);
  TH1D * tauthreeb5c4 = new TH1D("tauthreeb5c4",";Tau Three;N", 25, 0.0001, 64.9999);
  TH1D * tautworb5c1 = new TH1D("tautworb5c1","N-subjettiness Distribution;;N", 25, 0.0001, 1.4999);
  TH1D * tautworb5c2 = new TH1D("tautworb5c2",";;N", 25, 0.0001, 1.4999);
  TH1D * tautworb5c3 = new TH1D("tautworb5c3",";;N", 25, 0.0001, 1.4999);
  TH1D * tautworb5c4 = new TH1D("tautworb5c4",";Tau Two/Tau One;N", 25, 0.0001, 1.4999);
  TH1D * tauthreerb5c1 = new TH1D("tauthreerb5c1","N-subjettiness Distribution;;N", 25, 0.0001, 1.4999);
  TH1D * tauthreerb5c2 = new TH1D("tauthreerb5c2",";;N", 25, 0.0001, 1.4999);
  TH1D * tauthreerb5c3 = new TH1D("tauthreerb5c3",";;N", 25, 0.0001, 1.4999);
  TH1D * tauthreerb5c4 = new TH1D("tauthreerb5c4",";Tau Three/Tau Two;N", 25, 0.0001, 1.4999);

  TH1D * tauoneb6c1 = new TH1D("tauoneb6c1","N-subjettiness Distribution;;N", 25, 0.0001, 64.9999);
  TH1D * tauoneb6c2 = new TH1D("tauoneb6c2",";;N", 25, 0.0001, 64.9999);
  TH1D * tauoneb6c3 = new TH1D("tauoneb6c3",";;N", 25, 0.0001, 64.9999);
  TH1D * tauoneb6c4 = new TH1D("tauoneb6c4",";Tau One;N", 25, 0.0001, 64.9999);
  TH1D * tautwob6c1 = new TH1D("tautwob6c1","N-subjettiness Distribution;;N", 25, 0.0001, 64.9999);
  TH1D * tautwob6c2 = new TH1D("tautwob6c2",";;N", 25, 0.0001, 64.9999);
  TH1D * tautwob6c3 = new TH1D("tautwob6c3",";;N", 25, 0.0001, 64.9999);
  TH1D * tautwob6c4 = new TH1D("tautwob6c4",";Tau Two;N", 25, 0.0001, 64.9999);
  TH1D * tauthreeb6c1 = new TH1D("tauthreeb6c1","N-subjettiness Distribution;;N", 25, 0.0001, 64.9999);
  TH1D * tauthreeb6c2 = new TH1D("tauthreeb6c2",";;N", 25, 0.0001, 64.9999);
  TH1D * tauthreeb6c3 = new TH1D("tauthreeb6c3",";;N", 25, 0.0001, 64.9999);
  TH1D * tauthreeb6c4 = new TH1D("tauthreeb6c4",";Tau Three;N", 25, 0.0001, 64.9999);
  TH1D * tautworb6c1 = new TH1D("tautworb6c1","N-subjettiness Distribution;;N", 25, 0.0001, 1.4999);
  TH1D * tautworb6c2 = new TH1D("tautworb6c2",";;N", 25, 0.0001, 1.4999);
  TH1D * tautworb6c3 = new TH1D("tautworb6c3",";;N", 25, 0.0001, 1.4999);
  TH1D * tautworb6c4 = new TH1D("tautworb6c4",";Tau Two/Tau One;N", 25, 0.0001, 1.4999);
  TH1D * tauthreerb6c1 = new TH1D("tauthreerb6c1","N-subjettiness Distribution;;N", 25, 0.0001, 1.4999);
  TH1D * tauthreerb6c2 = new TH1D("tauthreerb6c2",";;N", 25, 0.0001, 1.4999);
  TH1D * tauthreerb6c3 = new TH1D("tauthreerb6c3",";;N", 25, 0.0001, 1.4999);
  TH1D * tauthreerb6c4 = new TH1D("tauthreerb6c4",";Tau Three/Tau Two;N", 25, 0.0001, 1.4999);
  
  
  double beta = 1.0;
  double beta2 = 2.0;
  double beta3 = 0.5;
  double beta4 = 0.2;
  double beta5 = 1.5;
  double beta6 = 3.0;
  
  double alpha = 1.0;
  
  WinnerTakeAllRecombiner wta_alpha(alpha);
  WinnerTakeAllRecombiner *wta;
  wta = &wta_alpha;
  
  UnnormalizedMeasure measureSpec1(beta);
  UnnormalizedMeasure measureSpec2(beta2);
  UnnormalizedMeasure measureSpec3(beta3);
  UnnormalizedMeasure measureSpec4(beta4);
  UnnormalizedMeasure measureSpec5(beta5);
  UnnormalizedMeasure measureSpec6(beta6);
  
  WTA_KT_Axes axisMode;

  int bin = tauonec1->GetSize();
  int bins = bin - 2;
  Long64_t nevents = c->GetEntries();
  
  Double_t R = 0.4;
  
  vector<PseudoJet>* jets_unsort = new vector<PseudoJet>;
  vector<PseudoJet>* trkLvect = new vector<PseudoJet>;
  vector<PseudoJet>* jets = new vector<PseudoJet>;
  for(Long64_t jentry = 0; jentry<nevents; jentry++){
    
    c->GetEntry(jentry);
    
    // if( !c->selectEvent() && !c->ismc) continue;
    
    
    if(jentry%100 == 0){
      cout << jentry << "/" << nevents << endl;
    }
    TLorentzVector tempVect;
    
    for(int trkIter = 0; trkIter < c->pf.nPFpart; trkIter++)
    {     
      if(c->pf.pfEta[trkIter]<2.0)
      {
        if(c->pf.pfEta[trkIter]>-2.0)
        {
          float VsPt = c->pf.pfPt[trkIter];
          float Eta = c->pf.pfEta[trkIter];
          float Phi = c->pf.pfPhi[trkIter];
          
          
          if(VsPt > 2){
            tempVect.SetPtEtaPhiM(VsPt, Eta, Phi, 0);
            trkLvect->push_back(tempVect);
            PseudoJet p = trkLvect->back();
            p.set_user_info(new MyUserInfo(c->pf.pfId[trkIter]));
            jets_unsort->push_back(p);
          }
        }  
      }    
    }  
    
    JetDefinition jet_def(cambridge_algorithm, R, wta, Best);

    ClusterSequence cs(*jets_unsort, jet_def);
    
    vector<PseudoJet> soft_jets = cs.inclusive_jets();
    // for(unsigned int j = 0; j<soft_jets.size();j++)
    // {
    // PseudoJet dropped = soft_drop(soft_jets[j]);
    // soft_jets[j] = dropped;
    // }
    vector<PseudoJet> jets = sorted_by_pt(soft_jets);
    int size = jets.size();
    

    
    for(unsigned int i =0; i<min(2,size); i++)
    {
      if(jets[i].perp() >100.0)
      { 
        vector<PseudoJet> constituents = jets[i].constituents();
        
        for(int constIter = 0; constIter < constituents.size(); constIter++)
        {
          const int & constID = constituents[constIter].user_info<MyUserInfo>().pdg_id();
          double constPt = constituents[constIter].perp();
          canpt->Fill(constPt);
          double FF = log(jets[i].perp()/constPt);
          jetFF->Fill(FF);
          jetFFc1->Fill(FF);
          jetFFc2->Fill(FF);   
          jetFFc3->Fill(FF);          
          jetFFc4->Fill(FF);
          
          if(abs(constID) <=3 && abs(constID) > 0)
          {
            trkcanpt->Fill(constPt);
          }
          if(abs(constID) == 4)
          {
            ECALcanpt->Fill(constPt);
          }
          if(abs(constID) <=7 && abs(constID) > 4)
          {
            HCALcanpt->Fill(constPt);
          }
        }
        
        jetpt->Fill(jets[i].perp());
        
        Nsubjettiness nSub1_beta1(1, axisMode, measureSpec1);
        double tau1 = nSub1_beta1(jets[i]);
        Nsubjettiness nSub2_beta1(2, axisMode, measureSpec1);
        double tau2 = nSub2_beta1(jets[i]);
        Nsubjettiness nSub3_beta1(3, axisMode, measureSpec1);
        double tau3 = nSub3_beta1(jets[i]);
        NsubjettinessRatio nSub21_beta1(2,1,axisMode, measureSpec1);
        double tau2r = nSub21_beta1(jets[i]);
        NsubjettinessRatio nSub32_beta1(3,2,axisMode, measureSpec1);
        double tau3r = nSub32_beta1(jets[i]);
        
        
        Nsubjettiness nSub1_beta2(1, axisMode, measureSpec2);
        double tau1_beta2 = nSub1_beta2(jets[i]);
        Nsubjettiness nSub2_beta2(2, axisMode, measureSpec2);
        double tau2_beta2 = nSub2_beta2(jets[i]);
        Nsubjettiness nSub3_beta2(3, axisMode, measureSpec2);
        double tau3_beta2 = nSub3_beta2(jets[i]);
        NsubjettinessRatio nSub21_beta2(2,1,axisMode, measureSpec2);
        double tau2r_beta2 = nSub21_beta2(jets[i]);
        NsubjettinessRatio nSub32_beta2(3,2,axisMode, measureSpec2);
        double tau3r_beta2 = nSub32_beta2(jets[i]);
        
        
        Nsubjettiness nSub1_beta3(1, axisMode, measureSpec3);
        double tau1_beta3 = nSub1_beta3(jets[i]);
        Nsubjettiness nSub2_beta3(2, axisMode, measureSpec3);
        double tau2_beta3 = nSub2_beta3(jets[i]);
        Nsubjettiness nSub3_beta3(3, axisMode,measureSpec3);
        double tau3_beta3 = nSub3_beta3(jets[i]);
        NsubjettinessRatio nSub21_beta3(2,1,axisMode, measureSpec3);
        double tau2r_beta3 = nSub21_beta3(jets[i]);
        NsubjettinessRatio nSub32_beta3(3,2,axisMode, measureSpec3);
        double tau3r_beta3 = nSub32_beta3(jets[i]);
        
        Nsubjettiness nSub1_beta4(1, axisMode, measureSpec4);
        double tau1_beta4 = nSub1_beta4(jets[i]);
        Nsubjettiness nSub2_beta4(2, axisMode, measureSpec4);
        double tau2_beta4 = nSub2_beta4(jets[i]);
        Nsubjettiness nSub3_beta4(3, axisMode, measureSpec4);
        double tau3_beta4 = nSub3_beta4(jets[i]);
        NsubjettinessRatio nSub21_beta4(2,1,axisMode, measureSpec4);
        double tau2r_beta4 = nSub21_beta1(jets[i]);
        NsubjettinessRatio nSub32_beta4(3,2,axisMode, measureSpec4);
        double tau3r_beta4 = nSub32_beta4(jets[i]);
        
        
        Nsubjettiness nSub1_beta5(1, axisMode, measureSpec5);
        double tau1_beta5 = nSub1_beta5(jets[i]);
        Nsubjettiness nSub2_beta5(2, axisMode, measureSpec5);
        double tau2_beta5 = nSub2_beta5(jets[i]);
        Nsubjettiness nSub3_beta5(3, axisMode, measureSpec5);
        double tau3_beta5 = nSub3_beta5(jets[i]);
        NsubjettinessRatio nSub21_beta5(2,1,axisMode, measureSpec5);
        double tau2r_beta5 = nSub21_beta5(jets[i]);
        NsubjettinessRatio nSub32_beta5(3,2,axisMode, measureSpec5);
        double tau3r_beta5 = nSub32_beta5(jets[i]);
        
        
        Nsubjettiness nSub1_beta6(1, axisMode, measureSpec6);
        double tau1_beta6 = nSub1_beta6(jets[i]);
        Nsubjettiness nSub2_beta6(2, axisMode, measureSpec6);
        double tau2_beta6 = nSub2_beta6(jets[i]);
        Nsubjettiness nSub3_beta6(3, axisMode,measureSpec6);
        double tau3_beta6 = nSub3_beta6(jets[i]);
        NsubjettinessRatio nSub21_beta6(2,1,axisMode, measureSpec6);
        double tau2r_beta6 = nSub21_beta6(jets[i]);
        NsubjettinessRatio nSub32_beta6(3,2,axisMode, measureSpec6);
        double tau3r_beta6 = nSub32_beta6(jets[i]);
        
        
        tauonec1->Fill(tau1);
        tauonec2->Fill(tau1);       
        tauonec3->Fill(tau1);     
        tauonec4->Fill(tau1);
        
        tautwoc1->Fill(tau2);       
        tautwoc2->Fill(tau2);       
        tautwoc3->Fill(tau2);       
        tautwoc4->Fill(tau2);
        
        tauthreec1->Fill(tau3);      
        tauthreec2->Fill(tau3);       
        tauthreec3->Fill(tau3);
        tauthreec4->Fill(tau3);
        
        tautworc1->Fill(tau2r);      
        tautworc2->Fill(tau2r);        
        tautworc3->Fill(tau2r);              
        tautworc4->Fill(tau2r);
        
        tauthreerc1->Fill(tau3r);
        tauthreerc2->Fill(tau3r);      
        tauthreerc3->Fill(tau3r);      
        tauthreerc4->Fill(tau3r);
        
        tauoneb2c1->Fill(tau1_beta2);
        tauoneb2c2->Fill(tau1_beta2);  
        tauoneb2c3->Fill(tau1_beta2);    
        tauoneb2c4->Fill(tau1_beta2);

        tautwob2c1->Fill(tau2_beta2);
        tautwob2c2->Fill(tau2_beta2);       
        tautwob2c3->Fill(tau2_beta2);      
        tautwob2c4->Fill(tau2_beta2);
        
        tauthreeb2c1->Fill(tau3_beta2);        
        tauthreeb2c2->Fill(tau3_beta2);       
        tauthreeb2c3->Fill(tau3_beta2);        
        tauthreeb2c4->Fill(tau3_beta2);

        tautworb2c1->Fill(tau2r_beta2);      
        tautworb2c2->Fill(tau2r_beta2);       
        tautworb2c3->Fill(tau2r_beta2);      
        tautworb2c4->Fill(tau2r_beta2);
        
        tauthreerb2c1->Fill(tau3r_beta2);    
        tauthreerb2c2->Fill(tau3r_beta2);       
        tauthreerb2c3->Fill(tau3r_beta2);
        tauthreerb2c4->Fill(tau3r_beta2);

        tauoneb3c1->Fill(tau1_beta3);      
        tauoneb3c2->Fill(tau1_beta3);    
        tauoneb3c3->Fill(tau1_beta3);       
        tauoneb3c4->Fill(tau1_beta3);
        
        tautwob3c1->Fill(tau2_beta3);       
        tautwob3c2->Fill(tau2_beta3);      
        tautwob3c3->Fill(tau2_beta3);       
        tautwob3c4->Fill(tau2_beta3);
        
        tauthreeb3c1->Fill(tau3_beta3);       
        tauthreeb3c2->Fill(tau3_beta3);       
        tauthreeb3c3->Fill(tau3_beta3);      
        tauthreeb3c4->Fill(tau3_beta3);
        
        tautworb3c1->Fill(tau2r_beta3);       
        tautworb3c2->Fill(tau2r_beta3);       
        tautworb3c3->Fill(tau2r_beta3);      
        tautworb3c4->Fill(tau2r_beta3);
        
        tauthreerb3c1->Fill(tau3r_beta3);       
        tauthreerb3c2->Fill(tau3r_beta3);       
        tauthreerb3c3->Fill(tau3r_beta3);    
        tauthreerb3c4->Fill(tau3r_beta3);
        
        tauoneb4c1->Fill(tau1_beta4);       
        tauoneb4c2->Fill(tau1_beta4);     
        tauoneb4c3->Fill(tau1_beta4);      
        tauoneb4c4->Fill(tau1_beta4);
        
        tautwob4c1->Fill(tau2_beta4);       
        tautwob4c2->Fill(tau2_beta4);      
        tautwob4c3->Fill(tau2_beta4);       
        tautwob4c4->Fill(tau2_beta4);
        
        tauthreeb4c1->Fill(tau3_beta4);       
        tauthreeb4c2->Fill(tau3_beta4);       
        tauthreeb4c3->Fill(tau3_beta4);      
        tauthreeb4c4->Fill(tau3_beta4);
        
        tautworb4c1->Fill(tau2r_beta4);       
        tautworb4c2->Fill(tau2r_beta4);       
        tautworb4c3->Fill(tau2r_beta4);       
        tautworb4c4->Fill(tau2r_beta4);
        
        tauthreerb4c1->Fill(tau3r_beta4);   
        tauthreerb4c2->Fill(tau3r_beta4);      
        tauthreerb4c3->Fill(tau3r_beta4);       
        tauthreerb4c4->Fill(tau3r_beta4);
        
        tauoneb5c1->Fill(tau1_beta5);       
        tauoneb5c2->Fill(tau1_beta5);     
        tauoneb5c3->Fill(tau1_beta5);     
        tauoneb5c4->Fill(tau1_beta5);
        
        tautwob5c1->Fill(tau2_beta5);      
        tautwob5c2->Fill(tau2_beta5);       
        tautwob5c3->Fill(tau2_beta5);       
        tautwob5c4->Fill(tau2_beta5);
        
        tauthreeb5c1->Fill(tau3_beta5);      
        tauthreeb5c2->Fill(tau3_beta5);       
        tauthreeb5c3->Fill(tau3_beta5);       
        tauthreeb5c4->Fill(tau3_beta5);
        
        tautworb5c1->Fill(tau2r_beta5);
        tautworb5c2->Fill(tau2r_beta5); 
        tautworb5c3->Fill(tau2r_beta5);
        tautworb5c4->Fill(tau2r_beta5);
        
        tauthreerb5c1->Fill(tau3r_beta5);       
        tauthreerb5c2->Fill(tau3r_beta5);      
        tauthreerb5c3->Fill(tau3r_beta5);     
        tauthreerb5c4->Fill(tau3r_beta5);
        
        tauoneb6c1->Fill(tau1_beta6);    
        tauoneb6c2->Fill(tau1_beta6);      
        tauoneb6c3->Fill(tau1_beta6);       
        tauoneb6c4->Fill(tau1_beta6);
        
        tautwob6c1->Fill(tau2_beta6);       
        tautwob6c2->Fill(tau2_beta6);       
        tautwob6c3->Fill(tau2_beta6);      
        tautwob6c4->Fill(tau2_beta6);
        
        tauthreeb6c1->Fill(tau3_beta6);        
        tauthreeb6c2->Fill(tau3_beta6);  
        tauthreeb6c3->Fill(tau3_beta6);      
        tauthreeb6c4->Fill(tau3_beta6);
        
        tautworb6c1->Fill(tau2r_beta6);       
        tautworb6c2->Fill(tau2r_beta6);       
        tautworb6c3->Fill(tau2r_beta6);       
        tautworb6c4->Fill(tau2r_beta6);
        
        tauthreerb6c1->Fill(tau3r_beta6);        
        tauthreerb6c2->Fill(tau3r_beta6);       
        tauthreerb6c3->Fill(tau3r_beta6);             
        tauthreerb6c4->Fill(tau3r_beta6);   
      }     
    }
    jets_unsort->clear();
  }
  
  
  double norm1canpt = 1.0/jetpt->GetEntries();
  canpt->Scale(norm1canpt);
  double norm1trkcanpt = 1.0/jetpt->GetEntries();
  trkcanpt->Scale(norm1canpt);
  double norm1ECALcanpt = 1.0/jetpt->GetEntries();
  ECALcanpt->Scale(norm1canpt);
  double norm1HCALcanpt = 1.0/jetpt->GetEntries();
  HCALcanpt->Scale(norm1canpt);
  
  double norm1jetpt = 1.0/jetpt->GetEntries();
  jetpt->Scale(norm1jetpt);
  
  double norm1FF = 1.0/jetpt->GetEntries();
  jetFF->Scale(norm1FF);
  jetFFc1->Scale(norm1FF);
  jetFFc2->Scale(norm1FF);
  jetFFc3->Scale(norm1FF);
  jetFFc4->Scale(norm1FF);

  double norm1c1 = 1.0/tauonec1->Integral();
  tauonec1->Scale(norm1c1);
  double norm1c2 = 1.0/tauonec2->Integral();
  tauonec2->Scale(norm1c2);
  double norm1c3 = 1.0/tauonec3->Integral();
  tauonec3->Scale(norm1c3);
  double norm1c4 = 1.0/tauonec4->Integral();
  tauonec4->Scale(norm1c4);
  
  double norm2c1 = 1.0/tautwoc1->Integral();
  tautwoc1->Scale(norm2c1);  
  double norm2c2 = 1.0/tautwoc2->Integral();
  tautwoc2->Scale(norm2c2);  
  double norm2c3 = 1.0/tautwoc3->Integral();
  tautwoc3->Scale(norm2c3);  
  double norm2c4 = 1.0/tautwoc4->Integral();
  tautwoc4->Scale(norm2c4);
  
  double norm3c1 = 1.0/tauthreec1->Integral();
  tauthreec1->Scale(norm3c1);  
  double norm3c2 = 1.0/tauthreec2->Integral();
  tauthreec2->Scale(norm3c2);  
  double norm3c3 = 1.0/tauthreec3->Integral();
  tauthreec3->Scale(norm3c3);  
  double norm3c4 = 1.0/tauthreec4->Integral();
  tauthreec4->Scale(norm3c4);
  
  double norm4c1 = 1.0/tautworc1->Integral();
  tautworc1->Scale(norm4c1);  
  double norm4c2 = 1.0/tautworc2->Integral();
  tautworc2->Scale(norm4c2);  
  double norm4c3 = 1.0/tautworc3->Integral();
  tautworc3->Scale(norm4c3);  
  double norm4c4 = 1.0/tautworc4->Integral();
  tautworc4->Scale(norm4c4);
  
  double norm5c1 = 1.0/tauthreerc1->Integral();
  tauthreerc1->Scale(norm5c1);  
  double norm5c2 = 1.0/tauthreerc2->Integral();
  tauthreerc2->Scale(norm5c2);  
  double norm5c3 = 1.0/tauthreerc3->Integral();
  tauthreerc3->Scale(norm5c3);  
  double norm5c4 = 1.0/tauthreerc4->Integral();
  tauthreerc4->Scale(norm5c4);

  
  double norm1b2c1 = 1.0/tauoneb2c1->Integral();
  tauoneb2c1->Scale(norm1b2c1);  
  double norm1b2c2 = 1.0/tauoneb2c2->Integral();
  tauoneb2c2->Scale(norm1b2c2);  
  double norm1b2c3 = 1.0/tauoneb2c3->Integral();
  tauoneb2c3->Scale(norm1b2c3);  
  double norm1b2c4 = 1.0/tauoneb2c4->Integral();
  tauoneb2c4->Scale(norm1b2c4);
  
  double norm2b2c1 = 1.0/tautwob2c1->Integral();
  tautwob2c1->Scale(norm2b2c1);  
  double norm2b2c2 = 1.0/tautwob2c2->Integral();
  tautwob2c2->Scale(norm2b2c2);  
  double norm2b2c3 = 1.0/tautwob2c3->Integral();
  tautwob2c3->Scale(norm2b2c3);  
  double norm2b2c4 = 1.0/tautwob2c4->Integral();
  tautwob2c4->Scale(norm2b2c4);
  
  double norm3b2c1 = 1.0/tauthreeb2c1->Integral();
  tauthreeb2c1->Scale(norm3b2c1);  
  double norm3b2c2 = 1.0/tauthreeb2c2->Integral();
  tauthreeb2c2->Scale(norm3b2c2);  
  double norm3b2c3 = 1.0/tauthreeb2c3->Integral();
  tauthreeb2c3->Scale(norm3b2c3);  
  double norm3b2c4 = 1.0/tauthreeb2c4->Integral();
  tauthreeb2c4->Scale(norm3b2c4);
  
  double norm4b2c1 = 1.0/tautworb2c1->Integral();
  tautworb2c1->Scale(norm4b2c1);  
  double norm4b2c2 = 1.0/tautworb2c2->Integral();
  tautworb2c2->Scale(norm4b2c2);  
  double norm4b2c3 = 1.0/tautworb2c3->Integral();
  tautworb2c3->Scale(norm4b2c3);  
  double norm4b2c4 = 1.0/tautworb2c4->Integral();
  tautworb2c4->Scale(norm4b2c4);
  
  double norm5b2c1 = 1.0/tauthreerb2c1->Integral();
  tauthreerb2c1->Scale(norm5b2c1);  
  double norm5b2c2 = 1.0/tauthreerb2c2->Integral();
  tauthreerb2c2->Scale(norm5b2c2);  
  double norm5b2c3 = 1.0/tauthreerb2c3->Integral();
  tauthreerb2c3->Scale(norm5b2c3);  
  double norm5b2c4 = 1.0/tauthreerb2c4->Integral();
  tauthreerb2c4->Scale(norm5b2c4);
  
  
  double norm1b3c1 = 1.0/tauoneb3c1->Integral();
  tauoneb3c1->Scale(norm1b3c1);  
  double norm1b3c2 = 1.0/tauoneb3c2->Integral();
  tauoneb3c2->Scale(norm1b3c2);  
  double norm1b3c3 = 1.0/tauoneb3c3->Integral();
  tauoneb3c3->Scale(norm1b3c3);  
  double norm1b3c4 = 1.0/tauoneb3c4->Integral();
  tauoneb3c4->Scale(norm1b3c4);
  
  double norm2b3c1 = 1.0/tautwob3c1->Integral();
  tautwob3c1->Scale(norm2b3c1);  
  double norm2b3c2 = 1.0/tautwob3c2->Integral();
  tautwob3c2->Scale(norm2b3c2);  
  double norm2b3c3 = 1.0/tautwob3c3->Integral();
  tautwob3c3->Scale(norm2b3c3);  
  double norm2b3c4 = 1.0/tautwob3c4->Integral();
  tautwob3c4->Scale(norm2b3c4);
  
  double norm3b3c1 = 1.0/tauthreeb3c1->Integral();
  tauthreeb3c1->Scale(norm3b3c1);  
  double norm3b3c2 = 1.0/tauthreeb3c2->Integral();
  tauthreeb3c2->Scale(norm3b3c2);  
  double norm3b3c3 = 1.0/tauthreeb3c3->Integral();
  tauthreeb3c3->Scale(norm3b3c3);  
  double norm3b3c4 = 1.0/tauthreeb3c4->Integral();
  tauthreeb3c4->Scale(norm3b3c4);
  
  double norm4b3c1 = 1.0/tautworb3c1->Integral();
  tautworb3c1->Scale(norm4b3c1);  
  double norm4b3c2 = 1.0/tautworb3c2->Integral();
  tautworb3c2->Scale(norm4b3c2);  
  double norm4b3c3 = 1.0/tautworb3c3->Integral();
  tautworb3c3->Scale(norm4b3c3);  
  double norm4b3c4 = 1.0/tautworb3c4->Integral();
  tautworb3c4->Scale(norm4b3c4);
  
  double norm5b3c1 = 1.0/tauthreerb3c1->Integral();
  tauthreerb3c1->Scale(norm5b3c1);  
  double norm5b3c2 = 1.0/tauthreerb3c2->Integral();
  tauthreerb3c2->Scale(norm5b3c2);  
  double norm5b3c3 = 1.0/tauthreerb3c3->Integral();
  tauthreerb3c3->Scale(norm5b3c3);  
  double norm5b3c4 = 1.0/tauthreerb3c4->Integral();
  tauthreerb3c4->Scale(norm5b3c4);
  
  
  double norm1b4c1 = 1.0/tauoneb4c1->Integral();
  tauoneb4c1->Scale(norm1b4c1);  
  double norm1b4c2 = 1.0/tauoneb4c2->Integral();
  tauoneb4c2->Scale(norm1b4c2);  
  double norm1b4c3 = 1.0/tauoneb4c3->Integral();
  tauoneb4c3->Scale(norm1b4c3);  
  double norm1b4c4 = 1.0/tauoneb4c4->Integral();
  tauoneb4c4->Scale(norm1b4c4);
  
  double norm2b4c1 = 1.0/tautwob4c1->Integral();
  tautwob4c1->Scale(norm2b4c1);  
  double norm2b4c2 = 1.0/tautwob4c2->Integral();
  tautwob4c2->Scale(norm2b4c2);  
  double norm2b4c3 = 1.0/tautwob4c3->Integral();
  tautwob4c3->Scale(norm2b4c3);  
  double norm2b4c4 = 1.0/tautwob4c4->Integral();
  tautwob4c4->Scale(norm2b4c4);
  
  double norm3b4c1 = 1.0/tauthreeb4c1->Integral();
  tauthreeb4c1->Scale(norm3b4c1);  
  double norm3b4c2 = 1.0/tauthreeb4c2->Integral();
  tauthreeb4c2->Scale(norm3b4c2);  
  double norm3b4c3 = 1.0/tauthreeb4c3->Integral();
  tauthreeb4c3->Scale(norm3b4c3);  
  double norm3b4c4 = 1.0/tauthreeb4c4->Integral();
  tauthreeb4c4->Scale(norm3b4c4);
  
  double norm4b4c1 = 1.0/tautworb4c1->Integral();
  tautworb4c1->Scale(norm4b4c1);  
  double norm4b4c2 = 1.0/tautworb4c2->Integral();
  tautworb4c2->Scale(norm4b4c2);  
  double norm4b4c3 = 1.0/tautworb4c3->Integral();
  tautworb4c3->Scale(norm4b4c3);  
  double norm4b4c4 = 1.0/tautworb4c4->Integral();
  tautworb4c4->Scale(norm4b4c4);
  
  double norm5b4c1 = 1.0/tauthreerb4c1->Integral();
  tauthreerb4c1->Scale(norm5b4c1);  
  double norm5b4c2 = 1.0/tauthreerb4c2->Integral();
  tauthreerb4c2->Scale(norm5b4c2);  
  double norm5b4c3 = 1.0/tauthreerb4c3->Integral();
  tauthreerb4c3->Scale(norm5b4c3);  
  double norm5b4c4 = 1.0/tauthreerb4c4->Integral();
  tauthreerb4c4->Scale(norm5b4c4);

  
  double norm1b5c1 = 1.0/tauoneb5c1->Integral();
  tauoneb5c1->Scale(norm1b5c1);  
  double norm1b5c2 = 1.0/tauoneb5c2->Integral();
  tauoneb5c2->Scale(norm1b5c2);  
  double norm1b5c3 = 1.0/tauoneb5c3->Integral();
  tauoneb5c3->Scale(norm1b5c3);  
  double norm1b5c4 = 1.0/tauoneb5c4->Integral();
  tauoneb5c4->Scale(norm1b5c4);
  
  double norm2b5c1 = 1.0/tautwob5c1->Integral();
  tautwob5c1->Scale(norm2b5c1);  
  double norm2b5c2 = 1.0/tautwob5c2->Integral();
  tautwob5c2->Scale(norm2b5c2);  
  double norm2b5c3 = 1.0/tautwob5c3->Integral();
  tautwob5c3->Scale(norm2b5c3);  
  double norm2b5c4 = 1.0/tautwob5c4->Integral();
  tautwob5c4->Scale(norm2b5c4);
  
  double norm3b5c1 = 1.0/tauthreeb5c1->Integral();
  tauthreeb5c1->Scale(norm3b5c1);  
  double norm3b5c2 = 1.0/tauthreeb5c2->Integral();
  tauthreeb5c2->Scale(norm3b5c2);  
  double norm3b5c3 = 1.0/tauthreeb5c3->Integral();
  tauthreeb5c3->Scale(norm3b5c3);  
  double norm3b5c4 = 1.0/tauthreeb5c4->Integral();
  tauthreeb5c4->Scale(norm3b5c4);
  
  double norm4b5c1 = 1.0/tautworb5c1->Integral();
  tautworb5c1->Scale(norm4b5c1);  
  double norm4b5c2 = 1.0/tautworb5c2->Integral();
  tautworb5c2->Scale(norm4b5c2);  
  double norm4b5c3= 1.0/tautworb5c3->Integral();
  tautworb5c3->Scale(norm4b5c3);  
  double norm4b5c4 = 1.0/tautworb5c4->Integral();
  tautworb5c4->Scale(norm4b5c4);
  
  double norm5b5c1 = 1.0/tauthreerb5c1->Integral();
  tauthreerb5c1->Scale(norm5b5c1);  
  double norm5b5c2 = 1.0/tauthreerb5c2->Integral();
  tauthreerb5c2->Scale(norm5b5c2);  
  double norm5b5c3 = 1.0/tauthreerb5c3->Integral();
  tauthreerb5c3->Scale(norm5b5c3);  
  double norm5b5c4 = 1.0/tauthreerb5c4->Integral();
  tauthreerb5c4->Scale(norm5b5c4);
  
  
  double norm1b6c1 = 1.0/tauoneb6c1->Integral();
  tauoneb6c1->Scale(norm1b6c1);  
  double norm1b6c2 = 1.0/tauoneb6c2->Integral();
  tauoneb6c2->Scale(norm1b6c2);  
  double norm1b6c3 = 1.0/tauoneb6c3->Integral();
  tauoneb6c3->Scale(norm1b6c3);  
  double norm1b6c4 = 1.0/tauoneb6c4->Integral();
  tauoneb6c4->Scale(norm1b6c4);
  
  double norm2b6c1 = 1.0/tautwob6c1->Integral();
  tautwob6c1->Scale(norm2b6c1);  
  double norm2b6c2 = 1.0/tautwob6c2->Integral();
  tautwob6c2->Scale(norm2b6c2);  
  double norm2b6c3 = 1.0/tautwob6c3->Integral();
  tautwob6c3->Scale(norm2b6c3);  
  double norm2b6c4 = 1.0/tautwob6c4->Integral();
  tautwob6c4->Scale(norm2b6c4);
  
  double norm3b6c1 = 1.0/tauthreeb6c1->Integral();
  tauthreeb6c1->Scale(norm3b6c1);  
  double norm3b6c2 = 1.0/tauthreeb6c2->Integral();
  tauthreeb6c2->Scale(norm3b6c2);  
  double norm3b6c3 = 1.0/tauthreeb6c3->Integral();
  tauthreeb6c3->Scale(norm3b6c3);  
  double norm3b6c4 = 1.0/tauthreeb6c4->Integral();
  tauthreeb6c4->Scale(norm3b6c4);
  
  double norm4b6c1 = 1.0/tautworb6c1->Integral();
  tautworb6c1->Scale(norm4b6c1);  
  double norm4b6c2 = 1.0/tautworb6c2->Integral();
  tautworb6c2->Scale(norm4b6c2);  
  double norm4b6c3 = 1.0/tautworb6c3->Integral();
  tautworb6c3->Scale(norm4b6c3);  
  double norm4b6c4 = 1.0/tautworb6c4->Integral();
  tautworb6c4->Scale(norm4b6c4);
  
  double norm5b6c1 = 1.0/tauthreerb6c1->Integral();
  tauthreerb6c1->Scale(norm5b6c1);  
  double norm5b6c2 = 1.0/tauthreerb6c2->Integral();
  tauthreerb6c2->Scale(norm5b6c2);  
  double norm5b6c3 = 1.0/tauthreerb6c3->Integral();
  tauthreerb6c3->Scale(norm5b6c3);  
  double norm5b6c4 = 1.0/tauthreerb6c4->Integral();
  tauthreerb6c4->Scale(norm5b6c4);

  outputfile->Write();
  outputfile->Close();
}

int main()
{
  cenpfvspp();
  return 0;
}