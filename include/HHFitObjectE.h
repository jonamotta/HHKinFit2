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
  HHFitObjectE(TLorentzVector const& initial4vector);
  
  Double_t getE() const;
  virtual TLorentzVector changeE(Double_t E) const = 0;
  void changeEandSave(Double_t E);
  virtual TLorentzVector scaleE(Double_t scale) const = 0;
  void scaleEandSave(Double_t scale);
  virtual TLorentzVector constrainEtoMinv(Double_t minv, TLorentzVector const& other4vector) const =0;
  void constrainEtoMinvandSave(Double_t minv, TLorentzVector const& other4vector);

  Double_t getUpperFitLimitE() const;
  Double_t getLowerFitLimitE() const;

  void setUpperFitLimitE(Double_t const upperlimit);
  void setUpperFitLimitE(Double_t const minv, TLorentzVector const& other4vectorMin);
  void setLowerFitLimitE(Double_t const lowerlimit);
  void setLowerFitLimitE(TLorentzVector const& other4vectorMin);

  virtual void print() const;

 private:
  Double_t m_upperLimitE;
  Double_t m_lowerLimitE;
};

#endif /* HHFitObjectE_ */
