#ifndef HHKinFitMasterHeavyHiggs_H
#define HHKinFitMasterHeavyHiggs_H

#include <TLorentzVector.h>
#include <TMatrixD.h>
#include <TVector2.h>

#include "HHLorentzVector.h"
#include "HHKinFit.h"

#include <stdio.h>
#include <map>
#include <utility>
#include <vector>
#include <sstream>

namespace HHKinFit2{

typedef std::pair< int, int >                                 HHFitHypothesisHeavyHiggs;
typedef std::map< HHFitHypothesisHeavyHiggs, double >         HHFitResultD;
typedef std::map< HHFitHypothesisHeavyHiggs, int >            HHFitResultI;
typedef std::map< HHFitHypothesisHeavyHiggs, bool >           HHFitResultB;
typedef std::map< HHFitHypothesisHeavyHiggs, TLorentzVector > HHFitResultTLor;

class HHKinFitMasterHeavyHiggs{
  public:
  HHKinFitMasterHeavyHiggs(TLorentzVector bjet1, double sigmaEbjet1,
                           TLorentzVector bjet2, double sigmaEbjet2,
                           TLorentzVector tauvis1,
                           TLorentzVector tauvis2,
                           TVector2 met, TMatrixD met_cov, 
			   bool istruth=false, TLorentzVector* heavyhiggsgen=0);

  //the main action, runs over all hypotheses and performs the fit
  void fit();

  void addHypo(HHFitHypothesisHeavyHiggs hypo);
  void addHypo(int mh1, int mh2);

  //Getters
  HHFitHypothesisHeavyHiggs getBestHypothesis();
  HHFitHypothesisHeavyHiggs getHypothesis(int mh1, int mh2);

  //Hypotheses
  void addMh1Mh2Hypothesis(HHFitHypothesisHeavyHiggs hypo);
  void addMh1Mh2Hypothesis(int mh1, int mh2);

  //Getters for fit results
  double getChi2(HHFitHypothesisHeavyHiggs hypo);
  double getChi2(int mh1, int mh2);

  double getChi2BJet1(HHFitHypothesisHeavyHiggs hypo);
  double getChi2BJet1(int mh1, int mh2);

  double getChi2BJet2(HHFitHypothesisHeavyHiggs hypo);
  double getChi2BJet2(int mh1, int mh2);
    
  double getChi2Balance(HHFitHypothesisHeavyHiggs hypo);
  double getChi2Balance(int mh1, int mh2);
  
  double getFitProb(HHFitHypothesisHeavyHiggs hypo);
  double getFitProb(int mh1, int mh2);
  
  double getMH(HHFitHypothesisHeavyHiggs hypo);
  double getMH(int mh1, int mh2);

  int getConvergence(HHFitHypothesisHeavyHiggs hypo);
  int getConvergence(int mh1, int mh2);

  TLorentzVector getFittedTau1(HHFitHypothesisHeavyHiggs hypo);
  TLorentzVector getFittedTau1(int mh1, int mh2);

  TLorentzVector getFittedTau2(HHFitHypothesisHeavyHiggs hypo);
  TLorentzVector getFittedTau2(int mh1, int mh2);

  TLorentzVector getFittedBJet1(HHFitHypothesisHeavyHiggs hypo);
  TLorentzVector getFittedBJet1(int mh1, int mh2);
  
  TLorentzVector getFittedBJet2(HHFitHypothesisHeavyHiggs hypo);
  TLorentzVector getFittedBJet2(int mh1, int mh2);
    
  //For Gen Studies to check smearing
  TLorentzVector getUnfittedBjet1(){return m_bjet1;};
  TLorentzVector getUnfittedBjet2(){return m_bjet2;};
  
private:  
  //hypotheses
  std::vector< HHFitHypothesisHeavyHiggs > m_hypos;

  //input vectors
  HHLorentzVector m_bjet1;
  HHLorentzVector m_bjet2;
  HHLorentzVector m_tauvis1;
  HHLorentzVector m_tauvis2;

  TVector2 m_MET;
  TMatrixD m_MET_COV;

  double m_sigma_bjet1;
  double m_sigma_bjet2;

  //fit result
  HHFitResultI m_map_convergence;
  
  HHFitResultD m_map_chi2;
  HHFitResultD m_map_chi2BJet1;
  HHFitResultD m_map_chi2BJet2;
  HHFitResultD m_map_chi2Balance;
  HHFitResultD m_map_prob;
  HHFitResultD m_map_mH;

  HHFitResultTLor m_map_fittedTau1;
  HHFitResultTLor m_map_fittedTau2;
  HHFitResultTLor m_map_fittedB1;
  HHFitResultTLor m_map_fittedB2;

  HHFitHypothesisHeavyHiggs m_bestHypo;
  double m_chi2_best;
  double m_mH_best;
};
}
#endif
