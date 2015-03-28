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

  virtual Double_t getChi2()=0;
  void setFitObject(HHFitObject* fitobject);
  
 protected:
  HHFitObject* m_fitobject;

};

#endif /* HHFitConstraint_ */
