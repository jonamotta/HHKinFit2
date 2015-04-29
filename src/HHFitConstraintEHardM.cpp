#include "HHFitConstraintEHardM.h"
#include "HHFitObjectE.h"

HHKinFit2::HHFitConstraintEHardM::HHFitConstraintEHardM(HHFitObject* fitobject, HHFitObject* constrainedobject, double mass)
  : HHFitConstraint(fitobject),
    m_constrainedobject(constrainedobject),
    m_mass(mass){

}

double
HHKinFit2::HHFitConstraintEHardM::getChi2() const{
  HHFitObjectE* new4momentum2 = static_cast<HHFitObjectE*>(m_constrainedobject);
  new4momentum2->constrainEtoMinvandSave(m_mass, m_fitobject->getFit4Vector());
  return(0);
}

double
HHKinFit2::HHFitConstraintEHardM::getLikelihood() const{
  return(0);
}
