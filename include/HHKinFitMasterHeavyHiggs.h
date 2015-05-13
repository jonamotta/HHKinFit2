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

typedef std::pair< int, int > HHFitHypothesisHeavyHiggs;
typedef std::map< HHFitHypothesisHeavyHiggs, double > HHFitResultD;
typedef std::map< HHFitHypothesisHeavyHiggs, int > HHFitResultI;
typedef std::map< HHFitHypothesisHeavyHiggs, bool > HHFitResultB;
typedef std::map< HHFitHypothesisHeavyHiggs, HHKinFit* > HHFitterMap;

class HHKinFitMasterHeavyHiggs{
  public:
  HHKinFitMasterHeavyHiggs(TLorentzVector &bjet1, double sigmaEbjet1,
                           TLorentzVector &bjet2, double sigmaEbjet2,
                           TLorentzVector &tauvis1,
                           TLorentzVector &tauvis2,
                           TVector2 &met, TMatrixD &met_cov, bool istruth=false);

  //the main action, runs over all hypotheses and performs the fit
  void fit();
  
  //Getters
  HHFitHypothesisHeavyHiggs getBestHypothesis();
  HHFitHypothesisHeavyHiggs getHypothesis(int mh1, int mh2);
  HHKinFit* getKinFit(HHFitHypothesisHeavyHiggs hypo);

  //Hypotheses
  void addMh1Mh2Hypothesis(int mh1, int mh2);

  //Getters for fit results
  double getChi2(HHFitHypothesisHeavyHiggs hypo);
  double getChi2(int mh1, int mh2);

  double getChi2BJet1(HHFitHypothesisHeavyHiggs hypo);
  double getChi2BJet1(int mh1, int mh2);

  double getChi2BJet2(HHFitHypothesisHeavyHiggs hypo);
  double getChi2Balance(HHFitHypothesisHeavyHiggs hypo);
  double getFitProb(HHFitHypothesisHeavyHiggs hypo);
  double getMH(HHFitHypothesisHeavyHiggs hypo);
  int getConvergence(HHFitHypothesisHeavyHiggs hypo);
  double getChi2BJet2(int mh1, int mh2);
  double getChi2Balance(int mh1, int mh2);
  double getFitProb(int mh1, int mh2);
  double getMH(int mh1, int mh2);
  int getConvergence(int mh1, int mh2);

private:
  //flag and variables needed for truth input
  bool m_truthInput;
  HHLorentzVector m_bjet1Smear;
  HHLorentzVector m_bjet2Smear;

  //input vectors
  HHLorentzVector m_bjet1;
  HHLorentzVector m_bjet2;
  HHLorentzVector m_tauvis1;
  HHLorentzVector m_tauvis2;

  TVector2 m_MET;
  TMatrixD m_MET_COV;

  //fitters
  HHFitterMap m_map_kinfits;

  //hypotheses
  std::vector< HHFitHypothesisHeavyHiggs > m_hypos;

  //fit result
  HHFitResultD m_map_chi2;
  HHFitResultD m_map_chi2BJet1;
  HHFitResultD m_map_chi2BJet2;
  HHFitResultD m_map_chi2Balance;
  HHFitResultD m_map_prob;
  HHFitResultD m_map_mH;

  HHFitHypothesisHeavyHiggs m_bestHypo;
};
}
#endif
