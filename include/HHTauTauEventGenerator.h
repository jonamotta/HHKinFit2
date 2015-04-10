/*
class for generating lorentzvectors for Tau in higgsdecay: Higgs->TauTau

*/

#include "HHLorentzVector.h"
#include "TF1.h"
#include "TRandom.h"
#include "TVector2.h"
 
class HHTauTauEventGenerator {

public:

HHTauTauEventGenerator(TF1* a);
 HHLorentzVector getTau1boosted();
 HHLorentzVector getTau2boosted();
 HHLorentzVector getTau1();
 HHLorentzVector getTau2();
 HHLorentzVector getTau1Vis();
 HHLorentzVector getTau2Vis();
 TVector2 getMET();
 void generateEvent();


private:
 TRandom m_randomnumber;
 TF1* m_PDF;
 HHLorentzVector m_tau1;
 HHLorentzVector m_tau2;
 HHLorentzVector m_isr;
 HHLorentzVector m_higgs;
 HHLorentzVector m_tau1boosted;
 HHLorentzVector m_tau2boosted;
 double m_visfrac1;
 double m_visfrac2;
 HHLorentzVector m_tau1vis;
 HHLorentzVector m_tau2vis;
 TVector2 m_MET;




};
