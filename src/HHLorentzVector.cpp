#include "HHLorentzVector.h"
#include <cmath>
#include "exceptions/HHEnergyRangeException.h"
#include <sstream>

HHKinFit2::HHLorentzVector::HHLorentzVector(double x, double y, double z, double t)
  : TLorentzVector(x,y,z,t){
}

HHKinFit2::HHLorentzVector::HHLorentzVector()
  : TLorentzVector(){
}

void
HHKinFit2::HHLorentzVector::SetEEtaPhiM(double E, double eta, double phi, double m){
  double p = sqrt(E*E-m*m);
  double sinth = sin(2*atan(exp(-eta)));
  double pt = sinth*p;
  this->SetPtEtaPhiE(pt,eta,phi,E);
}

HHKinFit2::HHLorentzVector
HHKinFit2::HHLorentzVector::operator+(HHLorentzVector const& rhs) const{
  return(HHLorentzVector(rhs.Px()+this->Px(),rhs.Py()+this->Py(),rhs.Pz()+this->Pz(),rhs.E()+this->E()));
}

HHKinFit2::HHLorentzVector
HHKinFit2::HHLorentzVector::operator-(HHLorentzVector const& rhs) const{
  return(HHLorentzVector(rhs.Px()-this->Px(),rhs.Py()-this->Py(),rhs.Pz()-this->Pz(),rhs.E()-this->E()));
}


void
HHKinFit2::HHLorentzVector::SetEkeepM(double E){
  if(E<M()){
    std::stringstream msg;
    msg << "target energy is smaller than the particle mass: "<<"E(set)="<<E<<" "<<"m="<<M();
    throw(HHEnergyRangeException(msg.str()));
  }
  
  double pnew = sqrt(pow(E,2)-pow(M(),2));
  double ptnew = pnew * sin(2.*atan(exp(-Eta())));
  SetPtEtaPhiE(ptnew,Eta(),Phi(),E);

}

void
HHKinFit2::HHLorentzVector::SetEkeepBeta(double E){
  double pnew = Beta()*E;
  double ptnew = pnew * sin(2.*atan(exp(-Eta())));
  SetPtEtaPhiE(ptnew,Eta(),Phi(),E);
}
