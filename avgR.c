#include "HiForestAnalysis-master/hiForest.h"
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
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/Njettiness.hh"
#include "fastjet/contrib/NjettinessPlugin.hh"
#include "fastjet/contrib/WinnerTakeAllRecombiner.hh"
#include "HiForestAnalysis-master/commonUtility.h"
#include <fstream>
#include "TLorentzVector.h"
#include <vector>
#include <string>
#ifdef __MAKECINT__
#pragma link C++ class vector<TLorentzVector>;
#endif


void avgR()
{
  TH1::SetDefaultSumw2();
  using namespace std;
  using namespace fastjet;
  using namespace fastjet::contrib;

  // HiForest * c = new HiForest("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/v85/HiForest_pt80_merged/HiForest_pt80_PYTHIA_ppReco_JECv85_merged_forest_0.root","forest",cPP, 0);
  // HiForest * c = new HiForest("/mnt/hadoop/cms/store/user/velicanu/HydjetDrum_Pyquen_Dijet80_FOREST_Track8_Jet24_FixedPtHatJES_v0/0.root","forest",cPbPb, 0);
  HiForest * c = new HiForest("/mnt/hadoop/cms/store/user/velicanu/HIRun2011-14Mar2014-v2-6lumi-jet80-forest-v4ANDv9-merged/0.root","forest",cPbPb, 0);
  c->LoadNoTrees();
  c->hasTrackTree = true;
  c->hasSkimTree = true;
  c->hasEvtTree = true;
  c->hasPFTree = true;
  c->hasTowerTree = true;
  
  // TFile * outputfile = new TFile("Rpythia.root", "recreate");
  // TFile * outputfile = new TFile("Rpp.root", "recreate");
  // TFile * outputfile = new TFile("RPH.root", "recreate");
  TFile * outputfile = new TFile("RPbPb.root", "recreate");
  
  
  TH1D * Rc1 = new TH1D("Rc1", "R distribution between 0-10% centralities; ;N", 50, 0.0001, 0.99999);
  TH1D * Rc2 = new TH1D("Rc2", "R distribution between 10-30% centralities; ;N", 50, 0.0001, 0.99999);
  TH1D * Rc3 = new TH1D("Rc3", "R distribution between 30-50% centralities; ;N", 50, 0.0001, 0.99999);
  TH1D * Rc4 = new TH1D("Rc4", "R distribution between 50-100% centralities; R ;N", 50, 0.0001, 0.99999);
  TH1D * Rm = new TH1D("Rm", "; R ;N", 50, 0.0001, 0.0999);
 
 
  TProfile *delR= new TProfile("delR","Profile of delta R versus centrality",100,0,100,0,5);
  const double pi = M_PI;
  
  int size = 0;
  
  double alpha = 1.0;
  
  WinnerTakeAllRecombiner wta_alpha(alpha);
  WinnerTakeAllRecombiner *wta;
  wta = &wta_alpha;
  
  

  int bin = Rc1->GetSize();
  int bins = bin - 2;
  Long64_t nevents = c->GetEntries();
  
  Double_t R = 0.4;
  
  vector<PseudoJet>* jets_unsort = new vector<PseudoJet>;
  vector<PseudoJet>* jets = new vector<PseudoJet>;

  
  for(Long64_t jentry = 0; jentry<30000; jentry++){
    
    c->GetEntry(jentry);
    
    if(!c->selectEvent())
    continue;
    
    int cen_bin = c->evt.hiBin;
    
    double centrality = cen_bin/2;
    
    
    if(jentry%100 == 0){
      cout << jentry << "/" << 30000 << endl;
    }
    TLorentzVector tempVect;
    
    for(int trkIter = 0; trkIter < c->pf.nPFpart; trkIter++)
    {     
      if(c->pf.pfEta[trkIter]<2.3)
      {
        if(c->pf.pfEta[trkIter]>-2.3)
        {
          
          float VsPt = c->pf.pfVsPt[trkIter];
          float Eta = c->pf.pfEta[trkIter];
          float Phi = c->pf.pfPhi[trkIter];
          
          if(VsPt > 0.01){
          tempVect.SetPtEtaPhiM(VsPt, Eta, Phi, 0);
          jets_unsort->push_back(tempVect);
          } 
        }  
      }    
    }  
    
    
    
    
    JetDefinition jet_def(antikt_algorithm, R);

    ClusterSequence cs(*jets_unsort, jet_def);
        
    vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());
    
    int size = jets.size();
    
    for(unsigned int i =0; i<min(2,size); i++)
    {
      if(jets[i].perp() >30.0)
      { 
      
        vector<PseudoJet> constituents = jets[i].constituents();
        double phi;
        vector<double> R_vect;
        vector<double> pt_vect;
        double jet_phi = jets[i].phi();
        double jet_eta=jets[i].eta();
        for(int j= 0; j < constituents.size(); j++)
        {
          double part_phi = constituents[j].phi();
          double part_eta = constituents[j].eta();
         
          double distance = getDR(jet_eta, jet_phi, part_eta, part_phi);
          R_vect.push_back(distance*constituents[j].perp2());
          pt_vect.push_back(constituents[j].perp2());
        }
        
        double R_sum = 0.0;
        double pt_sum = 0.0;
        for(int q = 0; q < R_vect.size(); q++)
        {
          R_sum += R_vect[i];
          pt_sum += pt_vect[i];
        }
        
        double avgR = R_sum/pt_sum;
        
        cout<<"avgR: "<<avgR<<endl;
        
       
        Rm->Fill(avgR);
        
        if(centrality<=10.0 && centrality >= 0.0){
          Rc1->Fill(avgR);
        }
        else if(centrality<=30.0 && centrality > 10.0){
          Rc2->Fill(avgR);
        }
        else if(centrality<=50.0 && centrality > 30.0){
          Rc3->Fill(avgR);
        }
        else if(centrality<=100.0 && centrality > 50.0){
          Rc4->Fill(avgR);
        }
        
        
        
      }
      
      
    }
    jets_unsort->clear();

  }
  
  Rm->Draw();
  
  double norm1c1 = 1.0/Rc1->Integral(0,bins);
  Rc1->Scale(norm1c1);
  double norm1c2 = 1.0/Rc2->Integral(0,bins);
  Rc2->Scale(norm1c2);
  double norm1c3 = 1.0/Rc3->Integral(0,bins);
  Rc3->Scale(norm1c3);
  double norm1c4 = 1.0/Rc4->Integral(0,bins);
  Rc4->Scale(norm1c4);
    
  TCanvas * c32 = new TCanvas();
  c32->Divide(4,1);
  c32->cd(4);
  Rc1->Draw();
  c32->cd(3);
  Rc2->Draw();
  c32->cd(2);
  Rc3->Draw();
  c32->cd(1);
  Rc4->Draw();
  
  delR->Draw();

  outputfile->Write();
  outputfile->Close();
}




int main(int argc, char *argv[])
{
  avgR();
  return 0;
}