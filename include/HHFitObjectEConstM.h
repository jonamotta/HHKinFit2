/*
 * class for fit objects with constant mass
 */

#ifndef HHFitObjectConstM_
#define HHFitObjectConstM_

#include "Rtypes.h"
#include "TLorentzVector.h"
#include "TMatrixD.h"

#include "HHFitObjectE.h"

class HHFitObjectEConstM : public HHFitObjectE {
 public:
  HHFitObjectEConstM(TLorentzVector const& initial4vector);
  TLorentzVector constrainEtoMinv(Double_t m, TLorentzVector const& other4vector) const;
  TLorentzVector changeE(Double_t E) const;
  TLorentzVector scaleE(Double_t scale) const;

  void print() const;

};

#endif /* HHFitObjectConstM_ */
