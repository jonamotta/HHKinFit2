#include "HHFitConstraintLikelihood.h"

HHFitConstraintLikelihood::HHFitConstraintLikelihood(HHFitObject* object, TF1* likelihood)
  : HHFitConstraint(object),
    m_likelihood(likelihood){

}

Double_t 
HHFitConstraintLikelihood::getChi2(){
  return(0);
}

void 
HHFitConstraintLikelihood::setLikelihood(TF1* likelihood){
  m_likelihood = likelihood;
}
