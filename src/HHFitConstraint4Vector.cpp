#include "HHFitConstraint4Vector.h"
#include "TMatrixD.h"
#include "TMatrixDEigen.h"
#include <iostream>
#include <sstream>
#include "exceptions/HHCovarianceMatrixException.h"
#include <cmath>
#include "TMath.h"


HHKinFit2::HHFitConstraint4Vector::HHFitConstraint4Vector(HHFitObject* object, bool px, bool py, bool pz, bool E)
  : HHFitConstraint(object),
    m_ncomp(0),
    m_cov(),
    m_invcov()
{
  m_components[0]=px;
  m_components[1]=py;
  m_components[2]=pz;
  m_components[3]=E;
  for (int i=0; i<4; i++) if (m_components[i]) m_ncomp++;
  m_cov.ResizeTo(m_ncomp,m_ncomp);

  for (int i=0;i<4;i++){
    for (int j=0; j<4; j++){
      if (m_components[i]&&m_components[j]) m_cov(i,j)=m_fitobject->getCovMatrix()(i,j);
    }
  }

  TMatrixDEigen eigenmatrix(m_cov);
  for (int i=0; i<m_ncomp; i++){
    if (eigenmatrix.GetEigenValues()(i,i)<0){
      std::stringstream msg;
      msg << "covariance matrix is not positive-definite, but has at least one negative eigenvalue";
      throw(HHCovarianceMatrixException(msg.str()));
    }
  }

  m_invcov.ResizeTo(m_ncomp,m_ncomp);
  m_invcov = m_cov;
  m_invcov.Invert();
}

double
HHKinFit2::HHFitConstraint4Vector::getChi2() const{
	 HHLorentzVector res = m_fitobject->getFit4Vector()-m_fitobject->getInitial4Vector();
	      double chi2sum=0;
	      for(int i=0; i<m_ncomp; i++){
	        for(int j=0;j<m_ncomp;j++){
	          chi2sum+=res(i)*res(j)*m_invcov(i,j);
	          //std::cout << res(i) << "*" <<res(j) << "*" << invcov(i,j) << "=" << res(i)*res(j)*invcov(i,j) << std::endl;
	        }
	      }
  return(chi2sum);
}

double
HHKinFit2::HHFitConstraint4Vector::getLikelihood() const{
    double likelihood=(1.0/pow(2*TMath::Pi(),m_ncomp*1.0/2.0))*(1.0/sqrt(m_cov.Determinant()))*exp(-0.5*this->getChi2());
    return(likelihood);



}
