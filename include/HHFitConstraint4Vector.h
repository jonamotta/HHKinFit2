/*
 * class for constraint on 4 momentum conservation
 */

#ifndef HHFitConstraint4Vector_
#define HHFitConstraint4Vector_

#include "Rtypes.h"
#include "HHFitConstraint.h"
#include "HHFitObject.h"
#include <vector>

class HHFitConstraint4Vector : public HHFitConstraint {
 public:
  HHFitConstraint4Vector(HHFitObject* object, Bool_t px, Bool_t py, Bool_t pz, Bool_t E);

  Double_t getChi2() const;

 private:
  Bool_t m_components[4];
};

#endif /* HHFitConstraint4Vector_ */
