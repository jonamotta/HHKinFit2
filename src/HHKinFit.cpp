#include "HHKinFit.h"
#include "HHFitConstraint4Vector.h"
#include "HHFitConstraintEHardM.h"
#include "HHFitConstraint.h"
#include "HHFitConstraintLikelihood.h"
#include "HHFitObjectEConstM.h"
#include "HHFitObjectE.h"
#include "HHFitObject.h"
#include "HHFitObjectComposite.h"


HHKinFit::HHKinFit()
: m_fitobjects(std::vector<HHFitObjectE*>()),
  m_constraints(std::vector<HHFitConstraint*>()){
}

void 
HHKinFit::Fit(){

}

Double_t 
HHKinFit::getChi2() const{
  Double_t chi2=0;
  for(std::vector<HHFitConstraint*>::const_iterator it = m_constraints.begin();it != m_constraints.end(); ++it)
    chi2 += (*it)->getChi2();
  return(chi2);
}

std::vector<HHFitObjectE*>
HHKinFit::getListOfFitObjects() const{
  return(m_fitobjects);
}

std::vector<HHFitConstraint*> 
HHKinFit::getListOfConstraints() const{
  return(m_constraints);
}

void
HHKinFit::addFitObjectE(HHFitObjectE* fitobject){
  m_fitobjects.push_back(fitobject);
}

void
HHKinFit::addConstraint(HHFitConstraint* constraint){
  m_constraints.push_back(constraint);
}
