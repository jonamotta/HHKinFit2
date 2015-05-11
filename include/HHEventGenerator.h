#ifndef HHEventGenerator_
#define HHEventGenerator_

#include "HHLorentzVector.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TVectorD.h"
 
namespace HHKinFit2{
class HHEventGenerator {

public:
  HHEventGenerator(const double mH, const double mh1, const double mh2, TF1* const a,TF1* const b, int seed=4357);
  void generateEvent();
  double getMH();
  double getMh1();
  double getMh2();
  TMatrixD getJetCov(HHLorentzVector const& j) const;
  double getJetSigma(HHLorentzVector const& j) const;
  HHLorentzVector simulateJet(HHLorentzVector &j);

  int m_seed;
private:
 int m_eventnumber;
 TRandom3 m_randomnumber;
 TF1 m_PDF1;
 TF1 m_PDF2;

 double m_mH;
 double m_mh1;
 double m_mh2;

 double m_sigmab1;
 double m_sigmab2;
 double m_sigmaisr;
 TMatrixD m_covrecoil;
};
}
#endif /* HHEventGenerator_ */
