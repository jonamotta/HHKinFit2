#include "HHFitObjectComposite.h"
#include <iostream>

HHFitObjectComposite::HHFitObjectComposite(std::vector<HHFitObject*> subobjects)
  : m_subobjects(subobjects){

}

TLorentzVector
HHFitObjectComposite::getFit4Vector(){
  TLorentzVector p(0,0,0,0);
  for (std::vector<HHFitObject*>::iterator it = m_subobjects.begin() ; it != m_subobjects.end(); ++it)
    p += (*it)->getFit4Vector();
  m_fit4vector = p;
  return(p);
}


TLorentzVector
HHFitObjectComposite::getInitial4Vector(){
  TLorentzVector p(0,0,0,0);
  for (std::vector<HHFitObject*>::iterator it = m_subobjects.begin() ; it != m_subobjects.end(); ++it)
    p += (*it)->getInitial4Vector();
  m_initial4vector = p;
  return(p);
}

TMatrixD
HHFitObjectComposite::getCovMatrix(){
  TMatrixD cov(4,4);
  for (std::vector<HHFitObject*>::iterator it = m_subobjects.begin() ; it != m_subobjects.end(); ++it)
    cov += (*it)->getCovMatrix();
  m_covmatrix=cov;
  return(cov);
}


void
HHFitObjectComposite::setSubobjects(std::vector<HHFitObject*> subobjects){
  m_subobjects=subobjects;
}


void
HHFitObjectComposite::addSubobject(HHFitObject* subobject){
  m_subobjects.push_back(subobject);
}

void
HHFitObjectComposite::print(){
  std::cout << "---" << std::endl;
  std::cout <<  "composite fit object:" << std::endl;
  this->printInitial4Vector();
  this->printFit4Vector();
  this->printCovMatrix();
}
