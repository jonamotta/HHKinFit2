/*
 * class for fit objects with constant beta (=p/E)
 */

#ifndef HHFitObjectConstBeta_
#define HHFitObjectConstBeta_

#include "Rtypes.h"
#include "TLorentzVector.h"
#include "TMatrixD.h"

#include "HHFitObjectE.h"

class HHFitObjectEConstBeta : public HHFitObjectE {
 public:
  HHFitObjectEConstBeta(TLorentzVector const& initial4vector);
  TLorentzVector constrainEtoMinv(Double_t m, TLorentzVector const& other4vector) const;

  void print() const;

};

#endif /* HHFitObjectConstBeta_ */
