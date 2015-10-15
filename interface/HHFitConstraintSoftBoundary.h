/*
 * class for soft energy limit
 */

#ifndef HHFitConstraintSoftBoundary_
#define HHFitConstraintSoftBoundary_

#include "HHFitConstraint.h"
#include "HHFitObject.h"
#include <vector>

namespace HHKinFit2{
class HHFitConstraintSoftBoundary : public HHFitConstraint {
 public:
  HHFitConstraintSoftBoundary(HHFitObject* object, double percent);
  double getChi2() const;
  double getLikelihood() const;

 private:
  double m_percent;
};
}
#endif /* HHFitConstraintSoftBoundary_ */
