/*
 * class for fit objects with constant beta (=p/E)
 */

#ifndef HHFitObjectConstBeta_
#define HHFitObjectConstBeta_

#include "Rtypes.h"
#include "HHLorentzVector.h"
#include "TMatrixD.h"

#include "HHFitObjectE.h"

class HHFitObjectEConstBeta : public HHFitObjectE {
 public:
  HHFitObjectEConstBeta(HHLorentzVector const& initial4vector);
  HHLorentzVector constrainEtoMinv(Double_t m, HHLorentzVector const& other4vector) const;

  void print() const;

};

#endif /* HHFitObjectConstBeta_ */
