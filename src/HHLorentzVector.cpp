#include "HHLorentzVector.h"
#include <cmath>
#include "exceptions/HHEnergyRangeException.h"
#include <sstream>

HHLorentzVector::HHLorentzVector(double x, double y, double z, double t)
  : TLorentzVector(x,y,z,t){
}

HHLorentzVector::HHLorentzVector()
  : TLorentzVector(){
}

void HHLorentzVector::SetEEtaPhiM(double E, double eta, double phi, double m){
  double p = sqrt(E*E-m*m);
  double sinth = sin(2*atan(exp(-eta)));
  double pt = sinth*p;
  this->SetPtEtaPhiE(pt,eta,phi,E);
}

HHLorentzVector HHLorentzVector::operator+(HHLorentzVector const& rhs) const{
  return(HHLorentzVector(rhs.Px()+this->Px(),rhs.Py()+this->Py(),rhs.Pz()+this->Pz(),rhs.E()+this->E()));
}

HHLorentzVector HHLorentzVector::operator-(HHLorentzVector const& rhs) const{
  return(HHLorentzVector(rhs.Px()-this->Px(),rhs.Py()-this->Py(),rhs.Pz()-this->Pz(),rhs.E()-this->E()));
}


void HHLorentzVector::SetEkeepM(double E){
  if(E<M()){
    std::stringstream msg;
    msg << "target energy is smaller than the particle mass: "<<"E(set)="<<E<<" "<<"m="<<M();
    throw(HHEnergyRangeException(msg.str()));
  }
  
  Double_t pnew = sqrt(pow(E,2)-pow(M(),2));
  Double_t ptnew = pnew * sin(2.*atan(exp(-Eta())));
  SetPtEtaPhiE(ptnew,Eta(),Phi(),E);

}
