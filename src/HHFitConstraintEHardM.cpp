#include "HHFitConstraintEHardM.h"
#include "HHFitObjectE.h"

HHFitConstraintEHardM::HHFitConstraintEHardM(HHFitObject* fitobject, HHFitObject* constrainedobject, Double_t mass)
  : HHFitConstraint(fitobject),
    m_constrainedobject(constrainedobject),
    m_mass(mass){

}

Double_t
HHFitConstraintEHardM::getChi2(){
  HHFitObjectE* new4momentum2 = static_cast<HHFitObjectE*>(m_constrainedobject);
  new4momentum2->constrainEtoMinvandSave(m_mass, m_fitobject->getFit4Vector());
  return(0);
}

void 
HHFitConstraintEHardM::setMass(Double_t mass){
  m_mass=mass;
}


