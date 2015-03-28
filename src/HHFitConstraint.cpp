#include "HHFitConstraint.h"

HHFitConstraint::HHFitConstraint(HHFitObject* fitobject)
  :m_fitobject(fitobject){

}

void HHFitConstraint::setFitObject(HHFitObject* fitobject){
  m_fitobject=fitobject;
}
