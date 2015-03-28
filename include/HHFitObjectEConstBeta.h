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
  HHFitObjectEConstBeta(TLorentzVector initial4vector);
  TLorentzVector constrainEtoMinv(Double_t m, TLorentzVector other4vector);

  void print();

};

#endif /* HHFitObjectConstBeta_ */
