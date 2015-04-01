#include "HHFitObjectComposite.h"
#include <iostream>

HHFitObjectComposite::HHFitObjectComposite(std::vector<HHFitObject*> const& subobjects)
  : m_subobjects(subobjects){
}

HHFitObjectComposite::HHFitObjectComposite(HHFitObject* subobject1, HHFitObject* subobject2, HHFitObject* subobject3)
  : m_subobjects(){
  this->addSubobject(subobject1);
  this->addSubobject(subobject2);
  this->addSubobject(subobject3);
}


TLorentzVector
HHFitObjectComposite::getFit4Vector() const{
  TLorentzVector p(0,0,0,0);
  for (std::vector<HHFitObject*>::const_iterator it = m_subobjects.begin() ; it != m_subobjects.end(); ++it)
    p += (*it)->getFit4Vector();
//  m_fit4vector = p;
  return(p);
}


TLorentzVector
HHFitObjectComposite::getInitial4Vector() const {
  TLorentzVector p(0,0,0,0);
  for (std::vector<HHFitObject*>::const_iterator it = m_subobjects.begin() ; it != m_subobjects.end(); ++it)
    p += (*it)->getInitial4Vector();
//  m_initial4vector = p;
  return(p);
}

TMatrixD
HHFitObjectComposite::getCovMatrix() const{
  TMatrixD cov(4,4);
  for (std::vector<HHFitObject*>::const_iterator it = m_subobjects.begin() ; it != m_subobjects.end(); ++it)
    cov += (*it)->getCovMatrix();
//  m_covmatrix=cov;
  return(cov);
}


void
HHFitObjectComposite::setSubobjects(std::vector<HHFitObject*> const& subobjects){
  m_subobjects=subobjects;
}


void
HHFitObjectComposite::addSubobject(HHFitObject* subobject){
  m_subobjects.push_back(subobject);
}

void
HHFitObjectComposite::print() const{
  std::cout << "---" << std::endl;
  std::cout <<  "composite fit object:" << std::endl;
  this->printInitial4Vector();
  this->printFit4Vector();
  this->printCovMatrix();
}
