#include "HHFitObjectComposite.h"
#include <iostream>

HHKinFit2::HHFitObjectComposite::HHFitObjectComposite(std::vector<HHFitObject*> const& subobjects)
  : m_subobjects(subobjects){
}

HHKinFit2::HHFitObjectComposite::HHFitObjectComposite(HHFitObject* subobject1, HHFitObject* subobject2, HHFitObject* subobject3)
  : m_subobjects(){
  this->addSubobject(subobject1);
  this->addSubobject(subobject2);
  this->addSubobject(subobject3);
}


HHKinFit2::HHLorentzVector
HHKinFit2::HHFitObjectComposite::getFit4Vector() const{
  HHLorentzVector p(0,0,0,0);
  for (std::vector<HHFitObject*>::const_iterator it = m_subobjects.begin() ; it != m_subobjects.end(); ++it)
    p += (*it)->getFit4Vector();
//  m_fit4vector = p;
  return(p);
}


HHKinFit2::HHLorentzVector
HHKinFit2::HHFitObjectComposite::getInitial4Vector() const {
  HHLorentzVector p(0,0,0,0);
  for (std::vector<HHFitObject*>::const_iterator it = m_subobjects.begin() ; it != m_subobjects.end(); ++it)
    p += (*it)->getInitial4Vector();
//  m_initial4vector = p;
  return(p);
}

TMatrixD
HHKinFit2::HHFitObjectComposite::getCovMatrix() const{
  TMatrixD cov(4,4);
  for (std::vector<HHFitObject*>::const_iterator it = m_subobjects.begin() ; it != m_subobjects.end(); ++it)
    cov += (*it)->getCovMatrix();
//  m_covmatrix=cov;
  return(cov);
}


void
HHKinFit2::HHFitObjectComposite::setSubobjects(std::vector<HHFitObject*> const& subobjects){
  m_subobjects=subobjects;
}


void
HHKinFit2::HHFitObjectComposite::addSubobject(HHFitObject* subobject){
  m_subobjects.push_back(subobject);
}

void
HHKinFit2::HHFitObjectComposite::print() const{
  std::cout << "---" << std::endl;
  std::cout <<  "composite fit object:" << std::endl;
  this->printInitial4Vector();
  this->printFit4Vector();
  std::cout <<  "-----------------" << std::endl;
  for (std::vector<HHFitObject*>::const_iterator it = m_subobjects.begin() ; it != m_subobjects.end(); ++it){
    (*it)->printFit4Vector();
  }
  std::cout <<  "-----------------" << std::endl;
  this->printCovMatrix();
}
