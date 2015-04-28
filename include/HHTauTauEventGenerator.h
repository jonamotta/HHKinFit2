/*
class for generating lorentzvectors for Tau in higgsdecay: Higgs->TauTau

*/

#ifndef HHTauTauEventGenerator_
#define HHTauTauEventGenerator_

#include "HHLorentzVector.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TVectorD.h"
 
namespace HHKinFit2{
class HHTauTauEventGenerator {

public:

HHTauTauEventGenerator(TF1 *a,TF1 *b, TMatrixD c, int seed=4357);
 HHLorentzVector getTau1boosted();
 HHLorentzVector getTau2boosted();
 HHLorentzVector getISR();
 HHLorentzVector getHiggs();
 HHLorentzVector getTau1();
 HHLorentzVector getTau2();
 HHLorentzVector getTau1Vis();
 HHLorentzVector getTau2Vis();
 void setMhiggs(double M);
 void setMtau1(double M);
 void setMtau2(double M);
 double getMhiggs();
 double getMtau1();
 double getMtau2();
 double getInveriantMass();
 double getAbsPtMET();
 double getPhiMET();
 double getAbsPtMETwithsigma();
 double getPhiMETwithsigma();
 TVectorD getMET();
 TVectorD getMETwithsigma();
 void generateEvent();
 double getvisfrac1();
 double getvisfrac2();
 void PrintLmatrix();
 void PrintCovarmatrix();
 TMatrixD getCovarmatrix();
 int m_seed;


private:
 int m_eventnumber;
 TRandom3 m_randomnumber;
 TF1 m_PDF1;
 TF1 m_PDF2;
 HHLorentzVector m_tau1;
 HHLorentzVector m_tau2;
 HHLorentzVector m_isr;
 HHLorentzVector m_higgs;
 HHLorentzVector m_tau1boosted;
 HHLorentzVector m_tau2boosted;
 HHLorentzVector m_isrwithsigma;
 double m_visfrac1;
 double m_visfrac2;
 HHLorentzVector m_tau1vis;
 HHLorentzVector m_tau2vis;
 TVectorD m_MET;
 TVectorD m_METwithsigma;
 double m_mhiggs;
 double m_mtau1;
 double m_mtau2;
 double m_pi;
 double m_misr;
 TMatrixD m_covarmatrix;
 TMatrixD m_L;


};
}
#endif /* HHTauTauEventGenerator_ */
