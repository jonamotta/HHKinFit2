/*
 * class for fit objects
 */

#ifndef HHFitObjectE_
#define HHFitObjectE_

#include "Rtypes.h"
#include "HHLorentzVector.h"
#include "TMatrixD.h"
#include "HHFitObject.h"


class HHFitObjectE : public HHFitObject {
 public:
  HHFitObjectE(HHLorentzVector const& initial4vector);
  
  Double_t getE() const;
  virtual HHLorentzVector changeE(Double_t E) const = 0;
  void changeEandSave(Double_t E);
  virtual HHLorentzVector scaleE(Double_t scale) const = 0;
  void scaleEandSave(Double_t scale);
  virtual HHLorentzVector constrainEtoMinv(Double_t minv, HHLorentzVector const& other4vector) const =0;
  void constrainEtoMinvandSave(Double_t minv, HHLorentzVector const& other4vector);

  Double_t getUpperFitLimitE() const;
  Double_t getLowerFitLimitE() const;

  void setUpperFitLimitE(Double_t const upperlimit);
  void setUpperFitLimitE(Double_t const minv, HHLorentzVector const& other4vectorMin);
  void setLowerFitLimitE(Double_t const lowerlimit);
  void setLowerFitLimitE(HHLorentzVector const& other4vectorMin);

  virtual void print() const;
  void printLimits() const;

 private:
  Double_t m_upperLimitE;
  Double_t m_lowerLimitE;
};

#endif /* HHFitObjectE_ */
