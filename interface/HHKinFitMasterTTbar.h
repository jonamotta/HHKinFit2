#ifndef HHKinFitMasterTTbar_H
#define HHKinFitMasterTTbar_H

#include <TLorentzVector.h>
#include <TMatrixD.h>
#include <TVector2.h>
#include <TCanvas.h>

#ifdef HHKINFIT2
#include "HHLorentzVector.h"
#include "HHKinFit.h"
#else
#include "HHKinFit2/HHKinFit2/interface/HHLorentzVector.h"
#include "HHKinFit2/HHKinFit2/interface/HHKinFit.h"
#endif

#include <stdio.h>
#include <utility>
#include <vector>
#include <sstream>

namespace HHKinFit2{

typedef std::pair< int, int >  HHFitHypothesis2M;

class HHKinFitMasterTTbar{
  public:
 
  HHKinFitMasterTTbar(TLorentzVector const& bjet1,
                      TLorentzVector const& bjet2,
                      TLorentzVector const& wvis1,
                      TLorentzVector const& wvis2,
                      TVector2 const& met, TMatrixD const& met_cov,
                      double sigmaEbjet1 = -1.0, double sigmaEbjet2 = -1.0);
  
  //the main action, runs over all hypotheses and performs the fit
  void fit();

  void useAdvancedBJetChi2(bool useAdvanceBJetChi2 = true);
  
  //Getters for fit results
  double getChi2();
  double getChi2BJet1();
  double getChi2BJet2();
  double getChi2Balance();
  double getFitProb();
  int getConvergence();
  TLorentzVector getFittedW1();
  TLorentzVector getFittedW2();
  TLorentzVector getFittedBJet1();
  TLorentzVector getFittedBJet2();
    
  //For Gen Studies to check smearing
  TLorentzVector getUnfittedBJet1(){return(m_bjet1);};
  TLorentzVector getUnfittedBJet2(){return(m_bjet2);};
  
  double getBJet1Resolution(){return m_sigma_bjet1;};
  double getBJet2Resolution(){return m_sigma_bjet2;};

  int m_loopsNeeded;
  int m_bestMethodFlag;
private:
  double GetPFBJetRes(double eta, double et);

  //hypotheses
  HHFitHypothesis2M m_hypo;

  //input vectors
  HHLorentzVector m_bjet1;
  HHLorentzVector m_bjet2;
  HHLorentzVector m_wvis1;
  HHLorentzVector m_wvis2;

  TVector2 m_MET;
  TMatrixD m_MET_COV;

  double m_sigma_bjet1;
  double m_sigma_bjet2;

  bool m_useAdveancedBJetChi2;
  //fit result
  int m_convergence;
  
  double m_chi2;
  double m_chi2BJet1;
  double m_chi2BJet2;
  double m_chi2Balance;
  double m_prob;

  TLorentzVector m_fittedW1;
  TLorentzVector m_fittedW2;
  TLorentzVector m_fittedB1;
  TLorentzVector m_fittedB2;

  //For ToyMC Studies 
  TMatrixD m_bjet1_COV;
  TMatrixD m_bjet2_COV;

};
}
#endif
