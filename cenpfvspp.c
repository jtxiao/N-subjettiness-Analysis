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
#include "fastjet/contrib/SoftDrop.hh"
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


void cenpfvspp()
{
  TH1::SetDefaultSumw2();
  using namespace std;
  using namespace fastjet;
  using namespace fastjet::contrib;

  HiForest * c = new HiForest("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/v85/HiForest_pt80_merged/HiForest_pt80_PYTHIA_ppReco_JECv85_merged_forest_0.root","forest",cPP, 0);
  // HiForest * c = new HiForest("/mnt/hadoop/cms/store/user/dgulhan/pp2013/data/prod24/hiForest_merged/HiForest_pp_Jet80_v8_PP2013_HiForest_PromptReco_JsonPP_Jet80_PPReco_merged_forest_0.root","forest",cPP, 0);

  // c->LoadNoTrees();
  // c->hasTrackTree = true;
  // c->hasSkimTree = true;
  // c->hasEvtTree = true;
  // c->hasPFTree = true;
  // c->hasTowerTree = true;
  
  TFile * outputfile = new TFile("taupythia.root", "recreate");
  // TFile * outputfile = new TFile("taupp.root", "recreate");
  // TFile * outputfile = new TFile("wom.root", "recreate");
  
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
  
  
  double beta_sd = 2.0;
  double zcut = 0.1;
  contrib::SoftDrop soft_drop(beta_sd, zcut);

  
  
  double beta = 1.0;
  double beta2 = 2.0;
  double beta3 = 0.5;
  double beta4 = 0.2;
  double beta5 = 1.5;
  double beta6 = 3.0;
  
  // double alpha = 1.0;
  
  // WinnerTakeAllRecombiner wta_alpha(alpha);
  // WinnerTakeAllRecombiner *wta;
  // wta = &wta_alpha;
  
  UnnormalizedMeasure measureSpec1(beta);
  UnnormalizedMeasure measureSpec2(beta2);
  UnnormalizedMeasure measureSpec3(beta3);
  UnnormalizedMeasure measureSpec4(beta4);
  UnnormalizedMeasure measureSpec5(beta5);
  UnnormalizedMeasure measureSpec6(beta6);
  
  WTA_KT_Axes axisMode;
  
  InitEtaSKGrid();
  InitPhiSKGrid();

  int bin = tauonec1->GetSize();
  int bins = bin - 2;
  Long64_t nevents = c->GetEntries();
  
  Double_t R = 0.4;
  
  vector<PseudoJet>* jets_unsort = new vector<PseudoJet>;
  vector<PseudoJet>* jets = new vector<PseudoJet>;

  
  for(Long64_t jentry = 0; jentry<nevents; jentry++){
    
    c->GetEntry(jentry);
    
    if(!c->selectEvent())
    continue;
    
    int cen_bin = c->evt.hiBin;
    
    double centrality = cen_bin/2;
    
    
    InitSKGrid();
    
    
    if(jentry%100 == 0){
      cout << jentry << "/" << nevents << endl;
    }
    TLorentzVector tempVect;
    
    for(int trkIter = 0; trkIter < c->pf.nPFpart; trkIter++)
    {     
      if(c->pf.pfEta[trkIter]<2.3)
      {
        if(c->pf.pfEta[trkIter]>-2.3)
        {
          float VsPt = c->pf.pfPt[trkIter];
          float Eta = c->pf.pfEta[trkIter];
          float Phi = c->pf.pfPhi[trkIter];
          
         
          float eventSKPtCut = getSKPtCut(c->pf.nPFpart, &VsPt, &Phi, &Eta);
          // cout<<eventSKPtCut<<endl;
          if(c->pf.pfPt[trkIter] > eventSKPtCut)
          {
            
            if(VsPt > 0.01){
              tempVect.SetPtEtaPhiM(VsPt, Eta, Phi, 0);
              jets_unsort->push_back(tempVect);
            } 
          }
        }  
      }    
    }  
    
    
    // for(int trkIter = 0; trkIter < c->track.nTrk; trkIter++)
    // {
    // if(c->track.pPt[trkIter]>0.5)
    // {
    // if(c->track.pEta[trkIter]<2.3)
    // {
    // if(c->track.pEta[trkIter]>-2.3)
    // {
    
    // if(c->track.pPt[trkIter] > 0)
    // {
    // tempVect.SetPtEtaPhiM(c->track.pPt[trkIter], c->track.pEta[trkIter], c->track.pPhi[trkIter], 0);
    // jets_unsort->push_back(tempVect);
    
    // }
    // }  
    // }  
    // }
    // }  
    
    JetDefinition jet_def(cambridge_algorithm, R);

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
      // cout<<jets[i].perp();
      if(jets[i].perp() >100.0)
      { 
        
        
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
  

  double norm1c1 = 1.0/tauonec1->Integral(0,bins);
  tauonec1->Scale(norm1c1);
  double norm1c2 = 1.0/tauonec2->Integral(0,bins);
  tauonec2->Scale(norm1c2);
  double norm1c3 = 1.0/tauonec3->Integral(0,bins);
  tauonec3->Scale(norm1c3);
  double norm1c4 = 1.0/tauonec4->Integral(0,bins);
  tauonec4->Scale(norm1c4);
  
  double norm2c1 = 1.0/tautwoc1->Integral(0,bins);
  tautwoc1->Scale(norm2c1);  
  double norm2c2 = 1.0/tautwoc2->Integral(0,bins);
  tautwoc2->Scale(norm2c2);  
  double norm2c3 = 1.0/tautwoc3->Integral(0,bins);
  tautwoc3->Scale(norm2c3);  
  double norm2c4 = 1.0/tautwoc4->Integral(0,bins);
  tautwoc4->Scale(norm2c4);
  
  double norm3c1 = 1.0/tauthreec1->Integral(0,bins);
  tauthreec1->Scale(norm3c1);  
  double norm3c2 = 1.0/tauthreec2->Integral(0,bins);
  tauthreec2->Scale(norm3c2);  
  double norm3c3 = 1.0/tauthreec3->Integral(0,bins);
  tauthreec3->Scale(norm3c3);  
  double norm3c4 = 1.0/tauthreec4->Integral(0,bins);
  tauthreec4->Scale(norm3c4);
  
  double norm4c1 = 1.0/tautworc1->Integral(0,bins);
  tautworc1->Scale(norm4c1);  
  double norm4c2 = 1.0/tautworc2->Integral(0,bins);
  tautworc2->Scale(norm4c2);  
  double norm4c3 = 1.0/tautworc3->Integral(0,bins);
  tautworc3->Scale(norm4c3);  
  double norm4c4 = 1.0/tautworc4->Integral(0,bins);
  tautworc4->Scale(norm4c4);
  
  double norm5c1 = 1.0/tauthreerc1->Integral(0,bins);
  tauthreerc1->Scale(norm5c1);  
  double norm5c2 = 1.0/tauthreerc2->Integral(0,bins);
  tauthreerc2->Scale(norm5c2);  
  double norm5c3 = 1.0/tauthreerc3->Integral(0,bins);
  tauthreerc3->Scale(norm5c3);  
  double norm5c4 = 1.0/tauthreerc4->Integral(0,bins);
  tauthreerc4->Scale(norm5c4);

  
  double norm1b2c1 = 1.0/tauoneb2c1->Integral(0,bins);
  tauoneb2c1->Scale(norm1b2c1);  
  double norm1b2c2 = 1.0/tauoneb2c2->Integral(0,bins);
  tauoneb2c2->Scale(norm1b2c2);  
  double norm1b2c3 = 1.0/tauoneb2c3->Integral(0,bins);
  tauoneb2c3->Scale(norm1b2c3);  
  double norm1b2c4 = 1.0/tauoneb2c4->Integral(0,bins);
  tauoneb2c4->Scale(norm1b2c4);
  
  double norm2b2c1 = 1.0/tautwob2c1->Integral(0,bins);
  tautwob2c1->Scale(norm2b2c1);  
  double norm2b2c2 = 1.0/tautwob2c2->Integral(0,bins);
  tautwob2c2->Scale(norm2b2c2);  
  double norm2b2c3 = 1.0/tautwob2c3->Integral(0,bins);
  tautwob2c3->Scale(norm2b2c3);  
  double norm2b2c4 = 1.0/tautwob2c4->Integral(0,bins);
  tautwob2c4->Scale(norm2b2c4);
  
  double norm3b2c1 = 1.0/tauthreeb2c1->Integral(0,bins);
  tauthreeb2c1->Scale(norm3b2c1);  
  double norm3b2c2 = 1.0/tauthreeb2c2->Integral(0,bins);
  tauthreeb2c2->Scale(norm3b2c2);  
  double norm3b2c3 = 1.0/tauthreeb2c3->Integral(0,bins);
  tauthreeb2c3->Scale(norm3b2c3);  
  double norm3b2c4 = 1.0/tauthreeb2c4->Integral(0,bins);
  tauthreeb2c4->Scale(norm3b2c4);
  
  double norm4b2c1 = 1.0/tautworb2c1->Integral(0,bins);
  tautworb2c1->Scale(norm4b2c1);  
  double norm4b2c2 = 1.0/tautworb2c2->Integral(0,bins);
  tautworb2c2->Scale(norm4b2c2);  
  double norm4b2c3 = 1.0/tautworb2c3->Integral(0,bins);
  tautworb2c3->Scale(norm4b2c3);  
  double norm4b2c4 = 1.0/tautworb2c4->Integral(0,bins);
  tautworb2c4->Scale(norm4b2c4);
  
  double norm5b2c1 = 1.0/tauthreerb2c1->Integral(0,bins);
  tauthreerb2c1->Scale(norm5b2c1);  
  double norm5b2c2 = 1.0/tauthreerb2c2->Integral(0,bins);
  tauthreerb2c2->Scale(norm5b2c2);  
  double norm5b2c3 = 1.0/tauthreerb2c3->Integral(0,bins);
  tauthreerb2c3->Scale(norm5b2c3);  
  double norm5b2c4 = 1.0/tauthreerb2c4->Integral(0,bins);
  tauthreerb2c4->Scale(norm5b2c4);
  
  
  double norm1b3c1 = 1.0/tauoneb3c1->Integral(0,bins);
  tauoneb3c1->Scale(norm1b3c1);  
  double norm1b3c2 = 1.0/tauoneb3c2->Integral(0,bins);
  tauoneb3c2->Scale(norm1b3c2);  
  double norm1b3c3 = 1.0/tauoneb3c3->Integral(0,bins);
  tauoneb3c3->Scale(norm1b3c3);  
  double norm1b3c4 = 1.0/tauoneb3c4->Integral(0,bins);
  tauoneb3c4->Scale(norm1b3c4);
  
  double norm2b3c1 = 1.0/tautwob3c1->Integral(0,bins);
  tautwob3c1->Scale(norm2b3c1);  
  double norm2b3c2 = 1.0/tautwob3c2->Integral(0,bins);
  tautwob3c2->Scale(norm2b3c2);  
  double norm2b3c3 = 1.0/tautwob3c3->Integral(0,bins);
  tautwob3c3->Scale(norm2b3c3);  
  double norm2b3c4 = 1.0/tautwob3c4->Integral(0,bins);
  tautwob3c4->Scale(norm2b3c4);
  
  double norm3b3c1 = 1.0/tauthreeb3c1->Integral(0,bins);
  tauthreeb3c1->Scale(norm3b3c1);  
  double norm3b3c2 = 1.0/tauthreeb3c2->Integral(0,bins);
  tauthreeb3c2->Scale(norm3b3c2);  
  double norm3b3c3 = 1.0/tauthreeb3c3->Integral(0,bins);
  tauthreeb3c3->Scale(norm3b3c3);  
  double norm3b3c4 = 1.0/tauthreeb3c4->Integral(0,bins);
  tauthreeb3c4->Scale(norm3b3c4);
  
  double norm4b3c1 = 1.0/tautworb3c1->Integral(0,bins);
  tautworb3c1->Scale(norm4b3c1);  
  double norm4b3c2 = 1.0/tautworb3c2->Integral(0,bins);
  tautworb3c2->Scale(norm4b3c2);  
  double norm4b3c3 = 1.0/tautworb3c3->Integral(0,bins);
  tautworb3c3->Scale(norm4b3c3);  
  double norm4b3c4 = 1.0/tautworb3c4->Integral(0,bins);
  tautworb3c4->Scale(norm4b3c4);
  
  double norm5b3c1 = 1.0/tauthreerb3c1->Integral(0,bins);
  tauthreerb3c1->Scale(norm5b3c1);  
  double norm5b3c2 = 1.0/tauthreerb3c2->Integral(0,bins);
  tauthreerb3c2->Scale(norm5b3c2);  
  double norm5b3c3 = 1.0/tauthreerb3c3->Integral(0,bins);
  tauthreerb3c3->Scale(norm5b3c3);  
  double norm5b3c4 = 1.0/tauthreerb3c4->Integral(0,bins);
  tauthreerb3c4->Scale(norm5b3c4);
  
  
  double norm1b4c1 = 1.0/tauoneb4c1->Integral(0,bins);
  tauoneb4c1->Scale(norm1b4c1);  
  double norm1b4c2 = 1.0/tauoneb4c2->Integral(0,bins);
  tauoneb4c2->Scale(norm1b4c2);  
  double norm1b4c3 = 1.0/tauoneb4c3->Integral(0,bins);
  tauoneb4c3->Scale(norm1b4c3);  
  double norm1b4c4 = 1.0/tauoneb4c4->Integral(0,bins);
  tauoneb4c4->Scale(norm1b4c4);
  
  double norm2b4c1 = 1.0/tautwob4c1->Integral(0,bins);
  tautwob4c1->Scale(norm2b4c1);  
  double norm2b4c2 = 1.0/tautwob4c2->Integral(0,bins);
  tautwob4c2->Scale(norm2b4c2);  
  double norm2b4c3 = 1.0/tautwob4c3->Integral(0,bins);
  tautwob4c3->Scale(norm2b4c3);  
  double norm2b4c4 = 1.0/tautwob4c4->Integral(0,bins);
  tautwob4c4->Scale(norm2b4c4);
  
  double norm3b4c1 = 1.0/tauthreeb4c1->Integral(0,bins);
  tauthreeb4c1->Scale(norm3b4c1);  
  double norm3b4c2 = 1.0/tauthreeb4c2->Integral(0,bins);
  tauthreeb4c2->Scale(norm3b4c2);  
  double norm3b4c3 = 1.0/tauthreeb4c3->Integral(0,bins);
  tauthreeb4c3->Scale(norm3b4c3);  
  double norm3b4c4 = 1.0/tauthreeb4c4->Integral(0,bins);
  tauthreeb4c4->Scale(norm3b4c4);
  
  double norm4b4c1 = 1.0/tautworb4c1->Integral(0,bins);
  tautworb4c1->Scale(norm4b4c1);  
  double norm4b4c2 = 1.0/tautworb4c2->Integral(0,bins);
  tautworb4c2->Scale(norm4b4c2);  
  double norm4b4c3 = 1.0/tautworb4c3->Integral(0,bins);
  tautworb4c3->Scale(norm4b4c3);  
  double norm4b4c4 = 1.0/tautworb4c4->Integral(0,bins);
  tautworb4c4->Scale(norm4b4c4);
  
  double norm5b4c1 = 1.0/tauthreerb4c1->Integral(0,bins);
  tauthreerb4c1->Scale(norm5b4c1);  
  double norm5b4c2 = 1.0/tauthreerb4c2->Integral(0,bins);
  tauthreerb4c2->Scale(norm5b4c2);  
  double norm5b4c3 = 1.0/tauthreerb4c3->Integral(0,bins);
  tauthreerb4c3->Scale(norm5b4c3);  
  double norm5b4c4 = 1.0/tauthreerb4c4->Integral(0,bins);
  tauthreerb4c4->Scale(norm5b4c4);

  
  double norm1b5c1 = 1.0/tauoneb5c1->Integral(0,bins);
  tauoneb5c1->Scale(norm1b5c1);  
  double norm1b5c2 = 1.0/tauoneb5c2->Integral(0,bins);
  tauoneb5c2->Scale(norm1b5c2);  
  double norm1b5c3 = 1.0/tauoneb5c3->Integral(0,bins);
  tauoneb5c3->Scale(norm1b5c3);  
  double norm1b5c4 = 1.0/tauoneb5c4->Integral(0,bins);
  tauoneb5c4->Scale(norm1b5c4);
  
  double norm2b5c1 = 1.0/tautwob5c1->Integral(0,bins);
  tautwob5c1->Scale(norm2b5c1);  
  double norm2b5c2 = 1.0/tautwob5c2->Integral(0,bins);
  tautwob5c2->Scale(norm2b5c2);  
  double norm2b5c3 = 1.0/tautwob5c3->Integral(0,bins);
  tautwob5c3->Scale(norm2b5c3);  
  double norm2b5c4 = 1.0/tautwob5c4->Integral(0,bins);
  tautwob5c4->Scale(norm2b5c4);
  
  double norm3b5c1 = 1.0/tauthreeb5c1->Integral(0,bins);
  tauthreeb5c1->Scale(norm3b5c1);  
  double norm3b5c2 = 1.0/tauthreeb5c2->Integral(0,bins);
  tauthreeb5c2->Scale(norm3b5c2);  
  double norm3b5c3 = 1.0/tauthreeb5c3->Integral(0,bins);
  tauthreeb5c3->Scale(norm3b5c3);  
  double norm3b5c4 = 1.0/tauthreeb5c4->Integral(0,bins);
  tauthreeb5c4->Scale(norm3b5c4);
  
  double norm4b5c1 = 1.0/tautworb5c1->Integral(0,bins);
  tautworb5c1->Scale(norm4b5c1);  
  double norm4b5c2 = 1.0/tautworb5c2->Integral(0,bins);
  tautworb5c2->Scale(norm4b5c2);  
  double norm4b5c3= 1.0/tautworb5c3->Integral(0,bins);
  tautworb5c3->Scale(norm4b5c3);  
  double norm4b5c4 = 1.0/tautworb5c4->Integral(0,bins);
  tautworb5c4->Scale(norm4b5c4);
  
  double norm5b5c1 = 1.0/tauthreerb5c1->Integral(0,bins);
  tauthreerb5c1->Scale(norm5b5c1);  
  double norm5b5c2 = 1.0/tauthreerb5c2->Integral(0,bins);
  tauthreerb5c2->Scale(norm5b5c2);  
  double norm5b5c3 = 1.0/tauthreerb5c3->Integral(0,bins);
  tauthreerb5c3->Scale(norm5b5c3);  
  double norm5b5c4 = 1.0/tauthreerb5c4->Integral(0,bins);
  tauthreerb5c4->Scale(norm5b5c4);
  
  
  double norm1b6c1 = 1.0/tauoneb6c1->Integral(0,bins);
  tauoneb6c1->Scale(norm1b6c1);  
  double norm1b6c2 = 1.0/tauoneb6c2->Integral(0,bins);
  tauoneb6c2->Scale(norm1b6c2);  
  double norm1b6c3 = 1.0/tauoneb6c3->Integral(0,bins);
  tauoneb6c3->Scale(norm1b6c3);  
  double norm1b6c4 = 1.0/tauoneb6c4->Integral(0,bins);
  tauoneb6c4->Scale(norm1b6c4);
  
  double norm2b6c1 = 1.0/tautwob6c1->Integral(0,bins);
  tautwob6c1->Scale(norm2b6c1);  
  double norm2b6c2 = 1.0/tautwob6c2->Integral(0,bins);
  tautwob6c2->Scale(norm2b6c2);  
  double norm2b6c3 = 1.0/tautwob6c3->Integral(0,bins);
  tautwob6c3->Scale(norm2b6c3);  
  double norm2b6c4 = 1.0/tautwob6c4->Integral(0,bins);
  tautwob6c4->Scale(norm2b6c4);
  
  double norm3b6c1 = 1.0/tauthreeb6c1->Integral(0,bins);
  tauthreeb6c1->Scale(norm3b6c1);  
  double norm3b6c2 = 1.0/tauthreeb6c2->Integral(0,bins);
  tauthreeb6c2->Scale(norm3b6c2);  
  double norm3b6c3 = 1.0/tauthreeb6c3->Integral(0,bins);
  tauthreeb6c3->Scale(norm3b6c3);  
  double norm3b6c4 = 1.0/tauthreeb6c4->Integral(0,bins);
  tauthreeb6c4->Scale(norm3b6c4);
  
  double norm4b6c1 = 1.0/tautworb6c1->Integral(0,bins);
  tautworb6c1->Scale(norm4b6c1);  
  double norm4b6c2 = 1.0/tautworb6c2->Integral(0,bins);
  tautworb6c2->Scale(norm4b6c2);  
  double norm4b6c3 = 1.0/tautworb6c3->Integral(0,bins);
  tautworb6c3->Scale(norm4b6c3);  
  double norm4b6c4 = 1.0/tautworb6c4->Integral(0,bins);
  tautworb6c4->Scale(norm4b6c4);
  
  double norm5b6c1 = 1.0/tauthreerb6c1->Integral(0,bins);
  tauthreerb6c1->Scale(norm5b6c1);  
  double norm5b6c2 = 1.0/tauthreerb6c2->Integral(0,bins);
  tauthreerb6c2->Scale(norm5b6c2);  
  double norm5b6c3 = 1.0/tauthreerb6c3->Integral(0,bins);
  tauthreerb6c3->Scale(norm5b6c3);  
  double norm5b6c4 = 1.0/tauthreerb6c4->Integral(0,bins);
  tauthreerb6c4->Scale(norm5b6c4);

  TCanvas * c1 = new TCanvas();
  c1->Divide(4,1);
  c1->cd(1);
  tauonec1->Draw();  
  c1->cd(2);
  tauonec2->Draw();  
  c1->cd(3);
  tauonec3->Draw();
  c1->cd(4);  
  tauonec4->Draw();
  
  TCanvas * c2 = new TCanvas();
  c2->Divide(4,1);
  c2->cd(1);
  tautwoc1->Draw();  
  c2->cd(2);
  tautwoc2->Draw();  
  c2->cd(3);
  tautwoc3->Draw();  
  c2->cd(4);
  tautwoc4->Draw();
  
  TCanvas * c3 = new TCanvas();
  c3->Divide(4,1);
  c3->cd(1);
  tauthreec1->Draw();
  c3->cd(2);
  tauthreec2->Draw();
  c3->cd(3);
  tauthreec3->Draw();
  c3->cd(4);
  tauthreec4->Draw();
  
  TCanvas * c4 = new TCanvas();
  c3->Divide(4,1);
  c4->cd(1);
  tautworc1->Draw();  
  c4->cd(2);
  tautworc2->Draw();  
  c4->cd(3);
  tautworc3->Draw();  
  c4->cd(4);
  tautworc4->Draw();
  
  TCanvas * c5 = new TCanvas();
  c5->Divide(4,1);
  c5->cd(1);
  tauthreerc1->Draw();  
  c5->cd(2);
  tauthreerc2->Draw();  
  c5->cd(3);
  tauthreerc3->Draw();  
  c5->cd(4);
  tauthreerc4->Draw();
  
  
  TCanvas * c6 = new TCanvas();
  c6->Divide(4,1);
  c6->cd(1);
  tauoneb2c1->Draw();  
  c6->cd(2);
  tauoneb2c2->Draw();  
  c6->cd(3);
  tauoneb2c3->Draw();  
  c6->cd(4);
  tauoneb2c4->Draw();
  
  TCanvas * c7 = new TCanvas();
  c7->Divide(4,1);
  c7->cd(1);
  tautwob2c1->Draw();  
  c7->cd(2);
  tautwob2c2->Draw();  
  c7->cd(3);
  tautwob2c3->Draw();  
  c7->cd(4);
  tautwob2c4->Draw();
  
  TCanvas * c8 = new TCanvas();
  c8->Divide(4,1);  
  c8->cd(1);
  tauthreeb2c1->Draw();  
  c8->cd(2);
  tauthreeb2c2->Draw();  
  c8->cd(3);
  tauthreeb2c3->Draw();  
  c8->cd(4);
  tauthreeb2c4->Draw();
  
  TCanvas * c9 = new TCanvas();
  c9->Divide(4,1);
  c9->cd(1);
  tautworb2c1->Draw();  
  c9->cd(2);
  tautworb2c2->Draw();  
  c9->cd(3);
  tautworb2c3->Draw();  
  c9->cd(4);
  tautworb2c4->Draw();
  
  TCanvas * c10 = new TCanvas();
  c10->Divide(4,1);
  c10->cd(1);
  tauthreerb2c1->Draw();  
  c10->cd(2);
  tauthreerb2c2->Draw();  
  c10->cd(3);
  tauthreerb2c3->Draw();  
  c10->cd(4);
  tauthreerb2c4->Draw();
  
  TCanvas * c11 = new TCanvas();
  c11->Divide(4,1);
  c11->cd(1);
  tauoneb3c1->Draw();  
  c11->cd(2);
  tauoneb3c2->Draw();  
  c11->cd(3);
  tauoneb3c3->Draw();  
  c11->cd(4);
  tauoneb3c4->Draw();
  
  TCanvas * c12 = new TCanvas();
  c12->Divide(4,1);
  c12->cd(1);
  tautwob3c1->Draw();  
  c12->cd(2);
  tautwob3c2->Draw();  
  c12->cd(3);
  tautwob3c3->Draw();  
  c12->cd(4);
  tautwob3c4->Draw();
  
  
  TCanvas * c13 = new TCanvas();
  c13->Divide(4,1);
  c13->cd(1);
  tauthreeb3c1->Draw();  
  c13->cd(2);
  tauthreeb3c2->Draw();  
  c13->cd(3);
  tauthreeb3c3->Draw();  
  c13->cd(4);
  tauthreeb3c4->Draw();
  
  
  
  TCanvas * c15 = new TCanvas();
  c15->Divide(4,1);
  c15->cd(1);
  tautworb3c1->Draw();  
  c15->cd(2);
  tautworb3c2->Draw();  
  c15->cd(3);
  tautworb3c3->Draw();  
  c15->cd(4);
  tautworb3c4->Draw();
  
  TCanvas * c16 = new TCanvas();
  c16->Divide(4,1);
  c16->cd(1);
  tauthreerb3c1->Draw();  
  c16->cd(2);
  tauthreerb3c2->Draw();  
  c16->cd(3);
  tauthreerb3c3->Draw();  
  c16->cd(4);
  tauthreerb3c4->Draw();
  
  
  TCanvas * c17 = new TCanvas();
  c17->Divide(4,1);  
  c17->cd(1);
  tauoneb4c1->Draw();  
  c17->cd(2);
  tauoneb4c2->Draw();  
  c17->cd(3);
  tauoneb4c3->Draw();  
  c17->cd(4);
  tauoneb4c4->Draw();
  
  
  TCanvas * c18 = new TCanvas();
  c18->Divide(4,1);
  c18->cd(1);
  tautwob4c1->Draw();  
  c18->cd(2);
  tautwob4c2->Draw();  
  c18->cd(3);
  tautwob4c3->Draw();  
  c18->cd(4);
  tautwob4c4->Draw();
  
  
  TCanvas * c19 = new TCanvas();
  c19->Divide(4,1);
  c19->cd(1);
  tauthreeb4c1->Draw();  
  c19->cd(2);
  tauthreeb4c2->Draw();  
  c19->cd(3);
  tauthreeb4c3->Draw();  
  c19->cd(4);
  tauthreeb4c4->Draw();
  
  
  TCanvas * c20 = new TCanvas();
  c20->Divide(4,1);
  c20->cd(1);
  tautworb4c1->Draw();  
  c20->cd(2);
  tautworb4c2->Draw();  
  c20->cd(3);
  tautworb4c3->Draw();  
  c20->cd(4);
  tautworb4c4->Draw();
  
  
  TCanvas * c21 = new TCanvas();
  c21->Divide(4,1);
  c21->cd(1);
  tauthreerb4c1->Draw();  
  c21->cd(2);
  tauthreerb4c2->Draw();  
  c21->cd(3);
  tauthreerb4c3->Draw();  
  c21->cd(4);
  tauthreerb4c4->Draw();
  
  
  TCanvas * c22 = new TCanvas();
  c22->Divide(4,1);
  c22->cd(1);
  tauoneb5c1->Draw();  
  c22->cd(2);
  tauoneb5c2->Draw();  
  c22->cd(3);
  tauoneb5c3->Draw();  
  c22->cd(4);
  tauoneb5c4->Draw();
  
  
  TCanvas * c23 = new TCanvas();
  c23->Divide(4,1);  
  c23->cd(1);
  tautwob5c1->Draw();  
  c23->cd(2);
  tautwob5c2->Draw();  
  c23->cd(3);
  tautwob5c3->Draw();  
  c23->cd(4);
  tautwob5c4->Draw();
  
  
  
  TCanvas * c24 = new TCanvas();
  c24->Divide(4,1);  
  c24->cd(1);
  tauthreeb5c1->Draw();  
  c24->cd(2);
  tauthreeb5c2->Draw();  
  c24->cd(3);
  tauthreeb5c3->Draw();  
  c24->cd(4);
  tauthreeb5c4->Draw();
  
  TCanvas * c25 = new TCanvas();
  c25->Divide(4,1);
  c25->cd(1);
  tautworb5c1->Draw();  
  c25->cd(2);
  tautworb5c2->Draw();  
  c25->cd(3);
  tautworb5c3->Draw();  
  c25->cd(4);
  tautworb5c4->Draw();
  
  TCanvas * c26 = new TCanvas();
  c26->Divide(4,1);
  c26->cd(1);
  tauthreerb5c1->Draw();
  c26->cd(2);
  tauthreerb5c2->Draw();
  c26->cd(3);
  tauthreerb5c3->Draw();
  c26->cd(4);
  tauthreerb5c4->Draw();
  
  TCanvas * c27 = new TCanvas();
  c27->Divide(4,1);
  c27->cd(1);
  tauoneb6c1->Draw();  
  c27->cd(2);
  tauoneb6c2->Draw();  
  c27->cd(3);
  tauoneb6c3->Draw();  
  c27->cd(4);
  tauoneb6c4->Draw();
  
  TCanvas * c28 = new TCanvas();
  c28->Divide(4,1);
  c28->cd(1);
  tautwob6c1->Draw();
  c28->cd(2);
  tautwob6c2->Draw();
  c28->cd(3);
  tautwob6c3->Draw();
  c28->cd(4);
  tautwob6c4->Draw();
  
  TCanvas * c29 = new TCanvas();
  c29->Divide(4,1);
  c29->cd(1);
  tauthreeb6c1->Draw();  
  c29->cd(2);
  tauthreeb6c2->Draw();  
  c29->cd(3);
  tauthreeb6c3->Draw();  
  c29->cd(4);
  tauthreeb6c4->Draw();
  
  TCanvas * c30 = new TCanvas();
  c30->Divide(4,1);
  c30->cd(1);
  tautworb6c1->Draw();
  c30->cd(2);
  tautworb6c2->Draw();
  c30->cd(3);
  tautworb6c3->Draw();
  c30->cd(4);
  tautworb6c4->Draw();
  
  TCanvas * c31 = new TCanvas();
  c31->Divide(4,1);
  c31->cd(1);
  tauthreerb6c1->Draw();
  c31->cd(2);
  tauthreerb6c2->Draw();
  c31->cd(3);
  tauthreerb6c3->Draw();
  c31->cd(4);
  tauthreerb6c4->Draw();


  outputfile->Write();
  outputfile->Close();
}




int main(int argc, char *argv[])
{
  cenpfvspp();
  return 0;
}