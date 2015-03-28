/*
 * class for additional likelihood-based constraints
 */

#ifndef HHFitConstraintLikelihood_
#define HHFitConstraintLikelihood_

#include "HHFitConstraint.h"
#include "TF1.h"

class HHFitConstraintLikelihood : public HHFitConstraint {
 public:
  HHFitConstraintLikelihood(HHFitObject* object, TF1* likelihood);

  Double_t getChi2();
  void setLikelihood(TF1* likelihood);

 private:
  TF1* m_likelihood;
};

#endif /* HHFitConstraintLikelihood_ */
