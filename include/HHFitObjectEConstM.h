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
  HHFitObjectEConstM(TLorentzVector initial4vector);
  TLorentzVector constrainEtoMinv(Double_t m, TLorentzVector other4vector);
  TLorentzVector changeE(Double_t E);
  TLorentzVector scaleE(Double_t scale);

  void print();

};

#endif /* HHFitObjectConstM_ */
