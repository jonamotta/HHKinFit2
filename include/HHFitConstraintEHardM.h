/*
 * class for constraint on 4 momentum conservation
 */

#ifndef HHFitConstraintEHardM_
#define HHFitConstraintEHardM_

#include "HHFitConstraint.h"
#include "HHFitObject.h"
#include <vector>

class HHFitConstraintEHardM : public HHFitConstraint {
 public:
  HHFitConstraintEHardM(HHFitObject* fitobject, HHFitObject* constrainedobject, Double_t mass);

  Double_t getChi2();
  void setMass(Double_t mass);

 private:
  HHFitObject* m_constrainedobject;
  Double_t m_mass;
  
};

#endif /* HHFitConstraintEHardM_ */
