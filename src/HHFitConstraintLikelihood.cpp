#include "HHFitConstraintLikelihood.h"

HHFitConstraintLikelihood::HHFitConstraintLikelihood(HHFitObject* object, TF1* likelihood)
  : HHFitConstraint(object),
    m_likelihood(likelihood){

}

Double_t 
HHFitConstraintLikelihood::getChi2() const{
  return(0);
}
