/*
 * class for fit objects with constant mass
 */

#ifndef HHFitObjectConstM_
#define HHFitObjectConstM_

#include "Rtypes.h"
#include "HHLorentzVector.h"
#include "TMatrixD.h"

#include "HHFitObjectE.h"

class HHFitObjectEConstM : public HHFitObjectE {
 public:
  HHFitObjectEConstM(HHLorentzVector const& initial4vector);
  HHLorentzVector constrainEtoMinv(Double_t m, HHLorentzVector const& other4vector) const;
  HHLorentzVector changeE(Double_t E) const;
  HHLorentzVector scaleE(Double_t scale) const;

  void print() const;

};

#endif /* HHFitObjectConstM_ */
