#ifndef HHKinFitMasterHeavyHiggs_H
#define HHKinFitMasterHeavyHiggs_H

#include <TLorentzVector.h>
#include <TMatrixD.h>
#include <TVector2.h>

#include <stdio.h>
#include <map>
#include <utility>
#include <vector>
#include <sstream>


namespace HHKinFit2{

  typedef std::pair< int, int > HHFitHypothesis;
  typedef std::map< HHFitHypothesis, double > HHFitResultD;
  typedef std::map< HHFitHypothesis, int > HHFitResultI;
  typedef std::map< HHFitHypothesis, bool > HHFitResultB;

  class HHKinFitMasterHeavyHiggs
  {
  public:
  HHKinFitMasterHeavyHiggs();

  void fit();
  
  //Setters
  void setBjet1(TLorentzVector &bjet1, double sigmaE);
  void setBjet2(TLorentzVector &bjet2, double sigmaE);
  void setTauvis1(TLorentzVector &tauvis1);
  void setTauvis2(TLorentzVector &tauvis2);
  void setMET(TVector2 &met, TMatrixD &met_cov);
  
  //Getters for fit results
  HHFitHypothesis getBestHypothesis();
  double getBestFitProb();
  double getBestChi2();
  double getBestChi2BJet1();
  double getBestChi2BJet2();
  double getBestChi2Balance();
  double getBestMH();
  double getConvergence();

  HHFitResultD getChi2();
  HHFitResultD getChi2BJet1();
  HHFitResultD getChi2BJet2();
  HHFitResultD getChi2Balance();
  HHFitResultD getFitProb();
  HHFitResultD getMH();
  HHFitResultI getConvergence();

  //Hypotheses
  void addMh1Hypothesis(std::vector<int> v);
  void addMh1Hypothesis(double m1, double m2=0, double m3=0, double m4=0, double m5=0, double m6=0, double m7=0, double m8=0, double m9=0, double m10=0);
  void addMh2Hypothesis(std::vector<int> v);
  void addMh2Hypothesis(double m1, double m2=0, double m3=0, double m4=0, double m5=0, double m6=0, double m7=0, double m8=0, double m9=0, double m10=0);


private:
  //hypotheses
  std::vector< int > m_mh1;
  std::vector< int > m_mh2;

  //input vectors
  TLorentzVector* m_bjet1;
  TLorentzVector* m_bjet2;
  TLorentzVector* m_tauvis1;
  TLorentzVector* m_tauvis2;

  TLorentzVector* m_MET;
  TMatrixD m_MET_COV;

  //full event fit
  bool m_truthInput;
  bool m_advancedBalance;
  double m_simpleBalancePt;
  double m_simpleBalanceUncert;
  HHFitResultD m_map_chi2;
  HHFitResultD m_map_chi2BJet1;
  HHFitResultD m_map_chi2BJet2;
  HHFitResultD m_map_chi2Balance;
  HHFitResultD m_map_prob;
  HHFitResultD m_map_mH;

  HHFitHypothesis m_bestHypo;
  double m_bestChi2;
  double m_bestChi2BJet1;
  double m_bestChi2BJet2;
  double m_bestChi2Balance;
  double m_bestProb;
  double m_bestMH;
};
}
#endif
