#include "HHFitConstraint4Vector.h"
#include "TMatrixD.h"
#include "TMatrixDEigen.h"
#include <iostream>
#include <sstream>
#include "exceptions/HHCovarianceMatrixException.h"


HHFitConstraint4Vector::HHFitConstraint4Vector(HHFitObject* object, Bool_t px, Bool_t py, Bool_t pz, Bool_t E)
  : HHFitConstraint(object){
  m_components[0]=px;
  m_components[1]=py;
  m_components[2]=pz;
  m_components[3]=E;
}

Double_t 
HHFitConstraint4Vector::getChi2(){
  int ncomp = 0;
  for (int i=0; i<4; i++) if (m_components[i])  ncomp++;
  TMatrixD cov(ncomp,ncomp);
  for (int i=0;i<4;i++){
    for (int j=0; j<4; j++){
      if (m_components[i]&&m_components[j]) cov(i,j)=m_fitobject->getCovMatrix()(i,j);
    }
  }

  TMatrixDEigen eigenmatrix(cov);
  for (int i=0; i<ncomp; i++)
  if (eigenmatrix.GetEigenValues()(i,i)<0){
    std::stringstream msg;
    msg << "covariance matrix is not positive-definite, but has at least one negative eigenvalue";
    throw(HHCovarianceMatrixException(msg.str()));
  }

  TMatrixD invcov = cov.Invert();
  TLorentzVector res = m_fitobject->getFit4Vector()-m_fitobject->getInitial4Vector();

  Double_t chi2sum=0;
  for(int i=0; i<ncomp; i++){
    for(int j=0;j<ncomp;j++){
      chi2sum+=res(i)*res(j)*invcov(i,j);
      //std::cout << res(i) << "*" <<res(j) << "*" << invcov(i,j) << "=" << res(i)*res(j)*invcov(i,j) << std::endl;
    }
  }
  return(chi2sum);
}

void 
HHFitConstraint4Vector::setUsedComponents(Bool_t px, Bool_t py, Bool_t pz, Bool_t E){
  m_components[0]=px;
  m_components[1]=py;
  m_components[2]=pz;
  m_components[3]=E;
}

