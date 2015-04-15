#include "HHFitConstraintLikelihood.h"

HHKinFit2::HHFitConstraintLikelihood::HHFitConstraintLikelihood(HHFitObject* object, TF1* likelihood)
  : HHFitConstraint(object),
    m_likelihood(likelihood){

}

double
HHKinFit2::HHFitConstraintLikelihood::getChi2() const{
  return(0);
}
