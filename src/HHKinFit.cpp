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
: m_fitobjects(std::vector<HHFitObject*>()),
  m_constraints(std::vector<HHFitConstraint*>()){
}

void 
HHKinFit::Fit(){

}

Double_t 
HHKinFit::getChi2(){
  Double_t chi2=0;
  for(std::vector<HHFitConstraint*>::iterator it = m_constraints.begin();it != m_constraints.end(); ++it)
    chi2 += (*it)->getChi2();
  return(chi2);
}

std::vector<HHFitObject*> 
HHKinFit::getListOfFitObjects(){
  return(m_fitobjects);
}

std::vector<HHFitConstraint*> 
HHKinFit::getListOfConstraints(){
  return(m_constraints);
}

void
HHKinFit::addFitObject(HHFitObject* fitobject){
  m_fitobjects.push_back(fitobject);
}

void
HHKinFit::addConstraint(HHFitConstraint* constraint){
  m_constraints.push_back(constraint);
}
