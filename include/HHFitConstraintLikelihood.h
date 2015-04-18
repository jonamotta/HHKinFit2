/*
 * class for additional likelihood-based constraints
 */

#ifndef HHFitConstraintLikelihood_
#define HHFitConstraintLikelihood_

#include "HHFitConstraint.h"
#include "TF1.h"
#include "HHFitObject.h"


namespace HHKinFit2{
class HHFitConstraintLikelihood : public HHFitConstraint {
 public:
  HHFitConstraintLikelihood(HHFitObject* object, TF1* likelihood);

  Double_t getChi2() const;

 private:
  TF1* const m_likelihood;
  HHFitObject* m_object;

};
}
#endif /* HHFitConstraintLikelihood_ */
