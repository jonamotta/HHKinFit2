#include "HHFitConstraintLikelihood.h"
#include <cmath>


HHKinFit2::HHFitConstraintLikelihood::HHFitConstraintLikelihood(HHFitObject* object, TF1* likelihood)
  : HHFitConstraint(object),
    m_likelihood(likelihood),
    m_object(object){

}

Double_t 
HHKinFit2::HHFitConstraintLikelihood::getChi2() const{
	double frac=m_object->getInitial4Vector().E()/m_object->getFit4Vector().E();
	TF1 pdfunction=*m_likelihood;
	double pdf=pdfunction(frac);
	double chi2=-2*log(pdf);
  return(chi2);
}
