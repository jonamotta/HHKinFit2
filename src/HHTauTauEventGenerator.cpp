#include "HHTauTauEventGenerator.h"
#include <cmath>
#include <iostream>
#include "TVector3.h"

HHKinFit2::HHTauTauEventGenerator::HHTauTauEventGenerator(TF1* a,TF1* b, TMatrixD c, int seed):
  m_seed(seed),
  m_PDF1(*a),
  m_PDF2(*b),
  m_randomnumber(m_seed),
  m_tau1(),
  m_tau2(),
  m_isr(),
  m_higgs(),
  m_tau1boosted(),
  m_tau2boosted(),
  m_visfrac1(),
  m_visfrac2(),
  m_tau1vis(),
  m_tau2vis(),
  m_MET(2),
  m_mhiggs(125.7),
  m_mtau1(1.77682),
  m_mtau2(1.77682),
  m_pi(3.14159265359),
  m_misr(0),
  m_covarmatrix(c),
  m_METwithsigma(2),
  m_L(2,2),
  m_isrwithsigma(),
  m_eventnumber(0)
  {
	m_L[0][0]=sqrt(c[0][0]);
	m_L[0][1]=(c[0][1]-c[1][0])/sqrt(c[1][1]-c[1][0]/sqrt(c[0][0]));
	m_L[1][0]=c[1][0]/sqrt(c[0][0]);
	m_L[1][1]=sqrt(c[1][1]-c[1][0]/sqrt(c[0][0]));
}
void HHKinFit2::HHTauTauEventGenerator::generateEvent() {
  m_seed=m_randomnumber.GetSeed();
  // Generate lorentzvectors for the two tau in the rest frame of the higgs
  m_eventnumber++;
  double cthtau1=m_randomnumber.Uniform(-1,1);
  double etatau1=log(tan(acos(cthtau1)/2));
  double phitau1=m_randomnumber.Uniform(0,2*m_pi);
  //double etau1=(m_mhiggs*m_mhiggs-m_mtau2*m_mtau2 +m_mtau1*m_mtau1)/(2*m_mhiggs);
  double etau1=m_mhiggs/2;
  m_tau1.SetEEtaPhiM(etau1, etatau1, phitau1, m_mtau1);
  //Tau 2

  //double etau2=(m_mhiggs*m_mhiggs-m_mtau1*m_mtau1 +m_mtau2*m_mtau2)/(2*m_mhiggs);
  double etau2=m_mhiggs/2;
  double etatau2=-etatau1;
  double phitau2=phitau1+m_pi;
  m_tau2.SetEEtaPhiM(etau2, etatau2, phitau2, m_mtau2);
  //generate lorentzvectors for isr-jet

  double ptisr=fabs(m_randomnumber.Gaus(0,30)); // absolut Value of Pt
  double phiisr=m_randomnumber.Uniform(0,2*m_pi);
  double etaisr=m_randomnumber.Uniform(-5,5);
  m_isr.SetPtEtaPhiM(ptisr,etaisr,phiisr,m_misr);
  //generate lorentzvector for higgs with lorentzvector from isr-jet

  double cthhiggs=m_randomnumber.Uniform(-1,1);
  double etahiggs=-log(tan(acos(cthhiggs)/2));
  double pthiggs=ptisr;
  double phihiggs=phiisr+m_pi;
  m_higgs.SetPtEtaPhiM(pthiggs,etahiggs,phihiggs,m_mhiggs);
  //boosting the two tau vectors

  double bx=m_higgs.Px()/m_higgs.E();
  double by=m_higgs.Py()/m_higgs.E();
  double bz=m_higgs.Pz()/m_higgs.E();
  m_tau1boosted=m_tau1;
  m_tau1boosted.Boost(bx,by,bz);
  m_tau2boosted=m_tau2;
  m_tau2boosted.Boost(bx,by,bz);
  //generate lorenzvector for tauvis
  //generate visfrac

  m_visfrac1=m_PDF1.GetRandom(0,1);
  //m_visfrac1=m_PDF1(m_randomnumber.Uniform(0,1));
  m_visfrac2=m_PDF2.GetRandom(0,1);
  //m_visfrac2=m_PDF2(m_randomnumber.Uniform(0,1));
  //tau1vis

  m_tau1vis=m_tau1boosted;
  m_tau1vis.SetEkeepM(m_visfrac1*m_tau1boosted.E());
  //tau2vis

  m_tau2vis=m_tau2boosted;
  m_tau2vis.SetEkeepM(m_visfrac2*m_tau2boosted.E());
  //find MET
  m_MET[0]= -(m_tau1vis.Px()+m_tau2vis.Px()+m_isr.Px());
  m_MET[1]= -(m_tau1vis.Py()+m_tau2vis.Py()+m_isr.Py());
  // modify MET with a covariance-matrix
  TVectorD Gausvector(2);
  Gausvector[0]=m_randomnumber.Gaus(0,1);
  Gausvector[1]=m_randomnumber.Gaus(0,1);
  m_METwithsigma=m_MET+m_L*Gausvector;



  HHLorentzVector sum=m_higgs+m_isr;
  if(sum.M()>13000){
	  std::cout << "Ecm in Event "<<m_eventnumber<<"  over 13 TeV " << sum.M() <<std::endl;
      this->generateEvent();
  }
}

HHKinFit2::HHLorentzVector HHKinFit2::HHTauTauEventGenerator::getTau1(){
  return(m_tau1);
}

HHKinFit2::HHLorentzVector HHKinFit2::HHTauTauEventGenerator::getTau2(){
  return(m_tau2);
}

HHKinFit2::HHLorentzVector  HHKinFit2::HHTauTauEventGenerator::getTau1boosted(){
  return(m_tau1boosted);
}

HHKinFit2::HHLorentzVector  HHKinFit2::HHTauTauEventGenerator::getTau2boosted(){
  return(m_tau2boosted);
}

HHKinFit2::HHLorentzVector  HHKinFit2::HHTauTauEventGenerator::getTau1Vis(){
  return(m_tau1vis);
}

HHKinFit2::HHLorentzVector  HHKinFit2::HHTauTauEventGenerator::getTau2Vis(){
  return(m_tau2vis);
}
HHKinFit2::HHLorentzVector HHKinFit2::HHTauTauEventGenerator::getISR(){
  return(m_isr);
}

HHKinFit2::HHLorentzVector HHKinFit2::HHTauTauEventGenerator::getHiggs(){
	return(m_higgs);
}

TVectorD  HHKinFit2::HHTauTauEventGenerator::getMET(){
  return(m_MET);
}

TVectorD  HHKinFit2::HHTauTauEventGenerator::getMETwithsigma(){
  return(m_METwithsigma);
}

double HHKinFit2::HHTauTauEventGenerator::getvisfrac1(){
  return(m_visfrac1);
}

double HHKinFit2::HHTauTauEventGenerator::getvisfrac2(){
  return(m_visfrac2);
}

void HHKinFit2::HHTauTauEventGenerator::setMhiggs(double M){
  m_mhiggs=M;
}

void HHKinFit2::HHTauTauEventGenerator::setMtau1(double M){
  m_mtau1=M;
}

void HHKinFit2::HHTauTauEventGenerator::setMtau2(double M){
  m_mtau2=M;
}

double HHKinFit2::HHTauTauEventGenerator::getMhiggs(){
  return(m_mhiggs);
}

double HHKinFit2::HHTauTauEventGenerator::getMtau1(){
  return(m_mtau1);
}

double HHKinFit2::HHTauTauEventGenerator::getMtau2(){
  return(m_mtau2);
}

void HHKinFit2::HHTauTauEventGenerator::PrintCovarmatrix(){
	m_covarmatrix.Print();
}

void HHKinFit2::HHTauTauEventGenerator::PrintLmatrix(){
	m_L.Print();
}

double HHKinFit2::HHTauTauEventGenerator::getInveriantMass(){
    HHLorentzVector temp=m_tau1+m_tau2;
	return(temp.M());
}

double HHKinFit2::HHTauTauEventGenerator::getAbsPtMET(){
	double temp=sqrt(m_MET[0]*m_MET[0]+m_MET[1]*m_MET[1]);
	return(temp);
}

double HHKinFit2::HHTauTauEventGenerator::getPhiMET(){
	TVector2 temp(m_MET[0],m_MET[1]);
	return(temp.Phi());
}

 double HHKinFit2::HHTauTauEventGenerator::getAbsPtMETwithsigma(){
	 double temp=sqrt(m_METwithsigma[0]*m_METwithsigma[0]+m_METwithsigma[1]*m_METwithsigma[1]);
	 return(temp);
}

 double HHKinFit2::HHTauTauEventGenerator::getPhiMETwithsigma(){
 	TVector2 temp(m_METwithsigma[0],m_METwithsigma[1]);
 	return(temp.Phi());
 }

 TMatrixD HHKinFit2::HHTauTauEventGenerator::getCovarmatrix(){
	 return(m_covarmatrix);
 }





