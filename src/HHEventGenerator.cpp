#include "HHEventGenerator.h"
#include "HHLorentzVector.h"
#include <cmath>
#include <iostream>
#include "TVector3.h"

HHKinFit2::HHEventGenerator::HHEventGenerator(const double mH, const double mh1, const double mh2, TF1* const a,TF1* const b, int seed):
  m_seed(seed),
  m_PDF1(*a),
  m_PDF2(*b),
  m_mH(mH),
  m_mh1(mh1),
  m_mh2(mh2),
  m_randomnumber(m_seed),
  m_sigmab1(0),
  m_sigmab2(0),
  m_sigmaisr(0),
  m_covrecoil(2,2),
  m_eventnumber(0)
  {
}
void HHKinFit2::HHEventGenerator::generateEvent() {
  m_seed=m_randomnumber.GetSeed();
  m_eventnumber++;

  /////////////////////////////////////////////////////////
  //recoil and heavy higgs
  double recoil_pt=fabs(m_randomnumber.Gaus(0,30)); // absolut Value of Pt
  double recoil_phi=m_randomnumber.Uniform(0,2*TMath::Pi());
  double recoil_eta=m_randomnumber.Uniform(-5,5);

  double higgs_cth =m_randomnumber.Uniform(-1,1);
  double higgs_eta =-log(tan(acos(higgs_cth)/2));
  double higgs_pt  =recoil_pt;
  double higgs_phi =recoil_phi+TMath::Pi();

  HHLorentzVector recoil;
  HHLorentzVector higgs;
  recoil.SetPtEtaPhiM(recoil_pt,recoil_eta,recoil_phi,0);
  higgs.SetPtEtaPhiM(higgs_pt,higgs_eta,higgs_phi,m_mH);


  /////////////////////////////////////////////////////////
  //higgs1 and higgs2
  HHLorentzVector higgs1;
  HHLorentzVector higgs2;

  double higgs1_cth=m_randomnumber.Uniform(-1,1);
  double higgs1_eta=log(tan(acos(higgs1_cth)/2));
  double higgs2_eta=-higgs1_eta;
  double higgs1_phi=m_randomnumber.Uniform(0,2*TMath::Pi());
  double higgs2_phi=higgs1_phi+TMath::Pi();
  double higgs1_E=(m_mH*m_mH-m_mh2*m_mh2 +m_mh1*m_mh1)/(2*m_mH);
  double higgs2_E=(m_mH*m_mH-m_mh1*m_mh1 +m_mh2*m_mh2)/(2*m_mH);

  higgs1.SetEEtaPhiM(higgs1_E, higgs1_eta, higgs1_phi, m_mh1);
  higgs2.SetEEtaPhiM(higgs2_E, higgs2_eta, higgs2_phi, m_mh2);

  double bx=higgs.Px()/higgs.E();
  double by=higgs.Py()/higgs.E();
  double bz=higgs.Pz()/higgs.E();

  higgs1.Boost(bx,by,bz);
  higgs2.Boost(bx,by,bz);

  /////////////////////////////////////////////////////////
  //tau1/2 of higgs1
  double mtau=1.77682; //GeV
  HHLorentzVector tau1;
  HHLorentzVector tau2;

  double tau1_cth=m_randomnumber.Uniform(-1,1);
  double tau1_eta=log(tan(acos(tau1_cth)/2));
  double tau2_eta=-tau1_eta;
  double tau1_phi=m_randomnumber.Uniform(0,2*TMath::Pi());
  double tau2_phi=tau1_phi+TMath::Pi();
  double tau1_E=m_mh1/2.0;
  double tau2_E=m_mh1/2.0;

  tau1.SetEEtaPhiM(tau1_E, tau1_eta, tau1_phi, mtau);
  tau2.SetEEtaPhiM(tau2_E, tau2_eta, tau2_phi, mtau);

  double bh1x=higgs1.Px()/higgs1.E();
  double bh1y=higgs1.Py()/higgs1.E();
  double bh1z=higgs1.Pz()/higgs1.E();
  tau1.Boost(bh1x,bh1y,bh1z);
  tau2.Boost(bh1x,bh1y,bh1z);

  double frac1=m_PDF1.GetRandom(0,1);
  double frac2=m_PDF2.GetRandom(0,1);
  tau1.SetEkeepM(frac1*tau1.E());
  tau2.SetEkeepM(frac2*tau2.E());

  /////////////////////////////////////////////////////////
  //bjet1/2 of higgs2
  double mb=4.18; //GeV
  HHLorentzVector b1;
  HHLorentzVector b2;

  double b1_cth=m_randomnumber.Uniform(-1,1);
  double b1_eta=log(tan(acos(b1_cth)/2));
  double b2_eta=-b1_eta;
  double b1_phi=m_randomnumber.Uniform(0,2*TMath::Pi());
  double b2_phi=b1_phi+TMath::Pi();
  double b1_E=m_mh2/2.0;
  double b2_E=m_mh2/2.0;

  b1.SetEEtaPhiM(b1_E, b1_eta, b1_phi, mb);
  b2.SetEEtaPhiM(b2_E, b2_eta, b2_phi, mb);

  double bh2x=higgs2.Px()/higgs2.E();
  double bh2y=higgs2.Py()/higgs2.E();
  double bh2z=higgs2.Pz()/higgs2.E();
  b1.Boost(bh2x,bh2y,bh2z);
  b2.Boost(bh2x,bh2y,bh2z);

  /////////////////////////////////////////////////////////
  //find MET
  TVectorD met(2);
  met(0) = -(tau1.Px()+tau2.Px()+b1.Px()+b2.Px()+recoil.Px());
  met(1) = -(tau1.Py()+tau2.Py()+b1.Py()+b2.Py()+recoil.Py());


  HHLorentzVector b1sim = simulateJet(b1);
  HHLorentzVector b2sim = simulateJet(b2);
  HHLorentzVector recoilsim = simulateJet(recoil);

  TVectorD metsim(2);
  metsim(0) = -(tau1.Px()+tau2.Px()+b1sim.Px()+b2sim.Px()+recoilsim.Px());
  metsim(1) = -(tau1.Py()+tau2.Py()+b1sim.Py()+b2sim.Py()+recoilsim.Py());

  TMatrixD cov(2,2);
  cov = getJetCov(b1) + getJetCov(b2) +getJetCov(recoil);
  cov.Print();

//  TMatrixD m_L;
//  m_L[0][0]=sqrt(c[0][0]);
//  m_L[0][1]=(c[0][1]-c[1][0])/sqrt(c[1][1]-c[1][0]/sqrt(c[0][0]));
//  m_L[1][0]=c[1][0]/sqrt(c[0][0]);
//  m_L[1][1]=sqrt(c[1][1]-c[1][0]/sqrt(c[0][0]));
//
//
//  // modify MET with a covariance-matrix
//  TVectorD Gausvector(2);
//  Gausvector[0]=m_randomnumber.Gaus(0,1);
//  Gausvector[1]=m_randomnumber.Gaus(0,1);
//  m_METwithsigma=m_MET+m_L*Gausvector;
//  HHLorentzVector sum=m_higgs+m_isr;
//  if(sum.M()>13000){
//	  std::cout << "Ecm in Event "<<m_eventnumber<<"  over 13 TeV " << sum.M() <<std::endl;
//      this->generateEvent();
//  }
}

double
HHKinFit2::HHEventGenerator::getJetSigma(HHLorentzVector const&j) const{
  return(j.E() * ( 0.05 + 0.05 * sqrt(50./j.E())));
}

TMatrixD
HHKinFit2::HHEventGenerator::getJetCov(HHLorentzVector const& j) const{
  double dE = getJetSigma(j);
  double dp = j.E()/j.P()*dE; // error propagation p=sqrt(e^2-m^2)
  double dpt = sin(j.Theta())*dp;

  TMatrixD cov(2,2);
  cov(0,0)=pow(dp*cos(j.Phi()),2);
  cov(0,1)=pow(dp,2)*cos(j.Phi())*sin(j.Phi());
  cov(1,0)=cov(0,1);
  cov(1,1)=pow(dp*sin(j.Phi()),2);

  return(cov);
}

HHKinFit2::HHLorentzVector
HHKinFit2::HHEventGenerator::simulateJet(HHLorentzVector &j) {
  HHLorentzVector newj=j;
  newj.SetEkeepBeta(m_randomnumber.Gaus(j.E(),getJetSigma(j)));
  return(newj);
}
