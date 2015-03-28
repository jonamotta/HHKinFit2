/*
 * class for general fit constraint
 */

#ifndef HHKinFit_
#define HHKinFit_

#include "Rtypes.h"
#include <vector>
#include "HHFitObject.h"
#include "HHFitConstraint.h"

class HHKinFit {
 public:
  HHKinFit();

  void Fit();
  Double_t getChi2();
  std::vector<HHFitObject*> getListOfFitObjects();
  std::vector<HHFitConstraint*> getListOfConstraints();
  
  void addFitObject(HHFitObject* fitobject);
  void addConstraint(HHFitConstraint* constraint);
  
 private:
  std::vector<HHFitObject*> m_fitobjects;
  std::vector<HHFitConstraint*> m_constraints;

};

#endif /* HHKinFit_ */
