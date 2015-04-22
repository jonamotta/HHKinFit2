#include "HHFitConstraintLikelihood.h"
#include <cmath>
#include <iostream>


HHKinFit2::HHFitConstraintLikelihood::HHFitConstraintLikelihood(HHFitObject* object1,HHFitObject* object2, TF1* likelihood1,TF1* likelihood2 )
  : HHFitConstraint(object1),
    m_likelihood1(likelihood1),
    m_likelihood2(likelihood2),
    m_object2(object2){

}

Double_t 
HHKinFit2::HHFitConstraintLikelihood::getChi2() const{
	double frac1=m_fitobject->getInitial4Vector().E()/m_fitobject->getFit4Vector().E();
	TF1 pdfunction1=*m_likelihood1;
	double pdf1=pdfunction1(frac1);
	double frac2=m_object2->getInitial4Vector().E()/m_object2->getFit4Vector().E();
	TF1 pdfunction2=*m_likelihood2;
	double pdf2=pdfunction2(frac2);
	double chi2=-2*log(pdf1)-2*log(pdf2);
  return(chi2);
}
