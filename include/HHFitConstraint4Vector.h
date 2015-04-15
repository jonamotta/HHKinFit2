/*
 * class for constraint on 4 momentum conservation
 */

#ifndef HHFitConstraint4Vector_
#define HHFitConstraint4Vector_

#include "HHFitConstraint.h"
#include "HHFitObject.h"
#include <vector>

namespace HHKinFit2{
class HHFitConstraint4Vector : public HHFitConstraint {
 public:
  HHFitConstraint4Vector(HHFitObject* object, bool px, bool py, bool pz, bool E);

  double getChi2() const;

 private:
  bool m_components[4];
};
}
#endif /* HHFitConstraint4Vector_ */
