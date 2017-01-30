#include <iostream>
#include <iomanip>
#include <atomic>
#include <math.h>

#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TH1D.h"

#include "TLorentzVector.h"
#include "TVector2.h"
#include "TMatrixD.h"

#include "HHKinFit2/HHKinFit2Scenarios/interface/HHKinFitMasterSingleHiggsSoftLimits.h"
#include "HHKinFit2/HHKinFit2Scenarios/interface/HHKinFitMasterSingleHiggsSoftLimitsPS.h"

#include "TGraph.h" 
#include "TFile.h" 

#include <chrono>

using HHKinFit2::HHKinFitMasterSingleHiggsSoftLimits;
using HHKinFit2::HHKinFitMasterSingleHiggsSoftLimitsPS;

int main(int argc, char* argv[])
{

  TFile* file = new TFile("/nfs/dust/cms/user/tlenz/13TeV/2016/CMSSW_8_0_20/src/DesyTauAnalyses/NTupleMaker/test/NTuple_Dec_ReReco/etau/final/GluGluHToTauTau_M125.root","READ");
  TTree* tree = (TTree*) file->Get("TauCheck");
  float pt_1, phi_1, eta_1, m_1, pt_2, phi_2, eta_2, m_2, met, metphi, metcov00, metcov10, metcov01, metcov11;

  //TTreeReader reader("TauCheck", file);
  tree->SetBranchAddress("pt_1",&pt_1);
  tree->SetBranchAddress("phi_1",&phi_1);
  tree->SetBranchAddress("eta_1",&eta_1);
  tree->SetBranchAddress("m_1",&m_1);

  tree->SetBranchAddress("pt_2",&pt_2);
  tree->SetBranchAddress("phi_2",&phi_2);
  tree->SetBranchAddress("eta_2",&eta_2);
  tree->SetBranchAddress("m_2",&m_2);

  tree->SetBranchAddress("met",&met);
  tree->SetBranchAddress("metphi",&metphi);
  tree->SetBranchAddress("metcov00",&metcov00);
  tree->SetBranchAddress("metcov10",&metcov10);
  tree->SetBranchAddress("metcov01",&metcov01);
  tree->SetBranchAddress("metcov11",&metcov11);

  Int_t events = tree->GetEntriesFast();
  double timer=0;
  double timerPS=0;
  Int_t skipcounter = 0;
  Int_t skipcounterPS = 0;
  Int_t errorcounter = 0;
  Int_t errorcounterPS = 0;

  TH1D h_minchi2PS("minchi2PS","",100,0,10);
  TH1D h_probPS("probPS","",20,0,1);
  TH1D h_etauPS("etauPS","",30,0,150);

  TH1D h_minchi2("minchi2","",100,0,10);
  TH1D h_prob("prob","",20,0,1);
  TH1D h_etau("etau","",30,0,150);

  TH1D h_Dminchi2("Dminchi2","",200,-1,1);
  TH1D h_Dprob("Dprob","",200,-1,1);
  TH1D h_Detau("Detau","",200,-1,1);

  for (Int_t i=0;i<events;i++) {
    tree->GetEvent(i);

    TLorentzVector tauvis1(0,0,0,0);
    TLorentzVector tauvis2(0,0,0,0);
    TVector2 met_vec(0,0);
    TMatrixD met_cov(4,4);

    tauvis1.SetPtEtaPhiM(pt_1,eta_1,phi_1,m_1);
    tauvis2.SetPtEtaPhiM(pt_2,eta_2,phi_2,m_2);

    met_vec.Set(met*cos(metphi), met*sin(metphi));
    met_cov[0][0]= metcov00;
    met_cov[0][1]= metcov01;
    met_cov[1][0]= metcov10;
    met_cov[1][1]= metcov11;

    double minchi2 = -1;
    double prob = -1;
    double etau = -1;

    double minchi2PS = -1;
    double probPS = -1;
    double etauPS = -1;
    
    try{
      HHKinFitMasterSingleHiggsSoftLimitsPS singlehiggsfitPS(tauvis1,tauvis2,met_vec,met_cov);
      singlehiggsfitPS.addHypo(125);
      auto start = std::chrono::high_resolution_clock::now();
      singlehiggsfitPS.fit();
      timerPS+= std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()-start).count();
     
      minchi2PS = singlehiggsfitPS.getChi2(125);
      probPS = singlehiggsfitPS.getFitProb(125);
      etauPS = singlehiggsfitPS.getFittedTau1(125).E();
    }
    catch(...){
      errorcounterPS++;
    }

    try{
      HHKinFitMasterSingleHiggsSoftLimits singlehiggsfit(tauvis1,tauvis2,met_vec,met_cov);
      singlehiggsfit.addHypo(125);
      auto start = std::chrono::high_resolution_clock::now();
      singlehiggsfit.fit();
      timer+= std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()-start).count();

      minchi2 = singlehiggsfit.getChi2(125);
      prob = singlehiggsfit.getFitProb(125);
      etau = singlehiggsfit.getFittedTau1(125).E();
    }
    catch(...){
      errorcounter++;
    }
   
    if(minchi2<0 || minchi2PS<0) {
      if (minchi2<0)   skipcounter++;
      if (minchi2PS<0) skipcounterPS++;
      continue;
    }
     
    h_minchi2PS.Fill(minchi2PS);
    h_probPS.Fill(probPS);
    h_etauPS.Fill(etauPS);
    
    h_minchi2.Fill(minchi2);
    h_prob.Fill(prob);
    h_etau.Fill(etau);
    
    h_Dminchi2.Fill((minchi2PS-minchi2)/minchi2);
    h_Dprob.Fill((probPS-prob)/prob);
    h_Detau.Fill((etauPS-etau)/etau);  
  }

  std::cout << "processed events: " << events << std::endl;

  std::cout << "unhandled exceptions (Minuit): " << errorcounter << std::endl;
  std::cout << "skipped events (Minuit): " << skipcounter << std::endl;
  std::cout << "time (Minuit): " << timer/events << "us/event" << std::endl;

  std::cout << "unhandled exceptions (PS): " << errorcounterPS << std::endl;
  std::cout << "skipped events (PS): " << skipcounterPS << std::endl;
  std::cout << "time (PS): " << timerPS/events << "us/event" << std::endl;



  TFile fout("comparison.root","RECREATE");
  h_minchi2PS.Write();
  h_probPS.Write();
  h_etauPS.Write();
    
  h_minchi2.Write();
  h_prob.Write();
  h_etau.Write();
  
  h_Dminchi2.Write();
  h_Dprob.Write();
  h_Detau.Write();
    
  fout.Close();
  
  return (0);
}
