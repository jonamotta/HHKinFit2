/*
 * class for fit objects
 */

#ifndef HHFitObjectE_
#define HHFitObjectE_

#include "Rtypes.h"
#include "TLorentzVector.h"
#include "TMatrixD.h"
#include "HHFitObject.h"


class HHFitObjectE : public HHFitObject {
 public:
  HHFitObjectE(TLorentzVector initial4vector);
  
  virtual TLorentzVector changeE(Double_t E)=0;
  void changeEandSave(Double_t E);

  virtual TLorentzVector scaleE(Double_t scale)=0;
  void scaleEandSave(Double_t scale);

  virtual TLorentzVector constrainEtoMinv(Double_t minv, TLorentzVector other4vector)=0;
  void constrainEtoMinvandSave(Double_t minv, TLorentzVector other4vector);

  Double_t getUpperFitLimitE();
  Double_t getLowerFitLimitE();

  void setUpperFitLimitE(Double_t upperlimit);
  void setUpperFitLimitE(Double_t minv, TLorentzVector other4vectorMin);
  void setLowerFitLimitE(Double_t lowerlimit);
  void setLowerFitLimitE(TLorentzVector other4vectorMin);


  virtual void print();

 private:
  Double_t m_upperLimitE;
  Double_t m_lowerLimitE;
};

#endif /* HHFitObjectE_ */
