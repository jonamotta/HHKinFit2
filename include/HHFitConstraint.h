/*
 * class for general fit constraint
 */

#ifndef HHFitConstraint_
#define HHFitConstraint_

#include "Rtypes.h"
#include "HHFitObject.h"

class HHFitConstraint {
 public:
  HHFitConstraint(HHFitObject* fitobject);
  virtual ~HHFitConstraint() {};

  virtual Double_t getChi2() const = 0;
  
 protected:
  HHFitObject* const m_fitobject;

};

#endif /* HHFitConstraint_ */
