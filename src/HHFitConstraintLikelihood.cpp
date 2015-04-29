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
	double chi2=-2*log(m_likelihood1->Eval(m_fitobject->getInitial4Vector().E()/m_fitobject->getFit4Vector().E()))-2*log(m_likelihood2->Eval(m_object2->getInitial4Vector().E()/m_object2->getFit4Vector().E()));
  return(chi2);
}

double
HHKinFit2::HHFitConstraintLikelihood::getLikelihood() const{
  return(m_likelihood1->Eval(m_fitobject->getInitial4Vector().E()/m_fitobject->getFit4Vector().E())*m_likelihood2->Eval(m_object2->getInitial4Vector().E()/m_object2->getFit4Vector().E()));
}

