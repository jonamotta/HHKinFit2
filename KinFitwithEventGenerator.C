#include "HHTauTauEventGenerator.h"
#include "TFile.h"
#include "TH1D.h"
#include "exceptions/HHEnergyRangeException.h"
#include "TDirectory.h"
#include <iostream>
#include "HHKinFit.h"
#include "HHFitConstraint4Vector.h"
#include "HHFitConstraintEHardM.h"
#include "HHFitConstraint.h"
#include "HHFitConstraintLikelihood.h"
#include "HHFitObjectEConstM.h"
#include "HHFitObjectE.h"
#include "HHFitObjectMET.h"
#include "HHFitObject.h"
#include "HHFitObjectComposite.h"
#include "HHLorentzVector.h"
#include "exceptions/HHCovarianceMatrixException.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TFile.h"
#include "HHFitConstraintLikelihood.h"

using namespace HHKinFit2;

int main(int argc, char* argv[])
{
  TF1 PDF1("PDF1","2*x",0,2);
  TF1 PDF2("PDF2","2*x",0,2);
  TF1* pdf1=&PDF1;
  TMatrixD covarmatrix(2,2);
  covarmatrix[0][0]=130;
  covarmatrix[0][1]=0;
  covarmatrix[1][0]=0;
  covarmatrix[1][1]=130;
  HHTauTauEventGenerator testgenerator(PDF1,PDF2,covarmatrix);
  TH1D h_invariantmass("h_invariantmass","invariant mass of the two tau-vectors in GeV",10,120,130);
  TH1D h_EtaTau1("h_EtaTau1","Eta from Tau1",100,-5,5);
  TH1D h_CthTau1("h_CthTau1","Cth from Tau1",100,-1,1);
  TH1D h_PhiTau1("h_PhiTau1","Phi from Tau1",100,-3.2,3.2);
  TH1D h_PtTau1("h_PtTau1","Pt from Tau1",100,0,65);
  TH1D h_ETau1("h_ETau1","Energy from Tau1",10,60,65);
  TH1D h_EtaTau2("h_EtaTau2","Eta from Tau2",100,-5,5);
  TH1D h_CthTau2("h_CthTau2","Cth from Tau2",100,-1,1);
  TH1D h_PhiTau2("h_PhiTau2","Phi from Tau2",100,-3.2,3.2);
  TH1D h_PtTau2("h_PtTau2","Pt from Tau2",100,0,65);
  TH1D h_ETau2("h_ETau2","Energy from Tau2",10,60,65);
  TH1D h_EtaISR("h_EtaISR","Eta from the ISR-jet",100,-6,6);
  TH1D h_CthISR("h_CthISR","Cth from ISR-jet",100,-1,1);
  TH1D h_PhiISR("h_PhiISR","Phi from ISR-jet",100,-3.2,3.2);
  TH1D h_PtISR("h_PtISR","Pt from ISR-jet",100,0,100);
  TH1D h_EtaHiggs("h_EtaHiggs","Eta from the Higgs",100,-5,5);
  TH1D h_CthHiggs("h_CthHiggs","Cth from Higgs",100,-1,1);
  TH1D h_PhiHiggs("h_PhiHiggs","Phi from Higgs",100,-3.2,3.2);
  TH1D h_PtHiggs("h_PtHiggs","Pt from Higgs",100,0,100);
  TH1D h_EtaTau1boosted("h_EtaTau1boosted","Eta from the boosted Tau1",100,-5,5);
  TH1D h_CthTau1boosted("h_CthTau1boosted","Cth from the boosted Tau1",100,-1,1);
  TH1D h_PhiTau1boosted("h_PhiTau1boosted","Phi from the boosted Tau1",100,-3.2,3.2);
  TH1D h_PtTau1boosted("h_PtTau1boosted","Pt from the boosted Tau1",100,0,100);
  TH1D h_EtaTau2boosted("h_EtaTau2boosted","Eta from the boosted Tau2",100,-5,5);
  TH1D h_CthTau2boosted("h_CthTau2boosted","Cth from the boosted Tau2",100,-1,1);
  TH1D h_PhiTau2boosted("h_PhiTau2boosted","Phi from the boosted Tau2",100,-3.2,3.2);
  TH1D h_PtTau2boosted("h_PtTau2boosted","Pt from the boosted Tau2",100,0,100);
  TH1D h_VisFracTau1("h_VisFracTau1","The Energy fraction from the visible Tau-Component form Tau1",100,0,1);
  TH1D h_VisFracTau2("h_VisFracTau2","The Energy fraction from the visible Tau-Component form Tau2",100,0,1);
  TH1D h_EtaTau1Vis("h_EtaTau1Vis","Eta from the visible component of Tau1",100,-5,5);
  TH1D h_CthTau1Vis("h_CthTau1Vis","Cth from the visible component of Tau1",100,-1,1);
  TH1D h_PhiTau1Vis("h_PhiTau1Vis","Phi from the visible component of Tau1",100,-3.2,3.2);
  TH1D h_PtTau1Vis("h_PtTau1Vis","Pt from the the visible component of Tau1",100,0,100);
  TH1D h_ETau1Vis("h_ETau1Vis","Energy from the the visible component of Tau1",100,0,100);
  TH1D h_EtaTau2Vis("h_EtaTau2Vis","Eta from the visible component of Tau2",100,-5,5);
  TH1D h_CthTau2Vis("h_CthTau2Vis","Cth from the visible component of Tau2",100,-1,1);
  TH1D h_PhiTau2Vis("h_PhiTau2Vis","Phi from the visible component of Tau2",100,-3.2,3.2);
  TH1D h_PtTau2Vis("h_PtTau2Vis","Pt from the the visible component of Tau2",100,0,100);
  TH1D h_ETau2Vis("h_ETau2Vis","Energy from the the visible component of Tau2",100,0,100);
  TH1D h_AbsPtmiss("h_AbsPtmiss","Absolute Value of the missing transverse momentum",100,0,100);
  TH1D h_PhiMET("h_PhiMET","Angle Phi of the missing transversal momentum vector",100,-3.2,3.2);
  TH1D h_AbsPtmisswithsigma("h_AbsPtmisswithsigma","Absolute Value of the missing transverse momentum with measurmental errors",100,0,100);
  TH1D h_PhiMETwithsigma("h_PhiMETwithsigma","Angle Phi of the missing transversal momentum vector with measurmental errors",100,-3.2,3.2);
  TH1D h_FitFinalChi2("h_FitFinalChi2","The Final chi2 from the KinFit",100,0,0.0025);
  TH1D h_FitlikelihoodFinalChi2("h_FitFinalChi2likelihood","The Final chi2 from the KinFit",100,0,0.0025);

  for(unsigned int i=0; i<100000; i++){
    try{
      testgenerator.generateEvent();
    }
    catch(const HHEnergyRangeException& e){
      //std::cout << e.what() << std::endl;
      i--;
      continue;
    }
    //EventGenerator control Plots
    h_invariantmass.Fill(testgenerator.getInveriantMass());
    h_EtaTau1.Fill(testgenerator.getTau1().Eta());
    h_CthTau1.Fill(testgenerator.getTau1().CosTheta());
    h_PhiTau1.Fill(testgenerator.getTau1().Phi());
    h_PtTau1.Fill(testgenerator.getTau1().Pt());
    h_ETau1.Fill(testgenerator.getTau1().E());
    h_EtaTau2.Fill(testgenerator.getTau2().Eta());
    h_CthTau2.Fill(testgenerator.getTau2().CosTheta());
    h_PhiTau2.Fill(testgenerator.getTau2().Phi());
    h_PtTau2.Fill(testgenerator.getTau2().Pt());
    h_ETau2.Fill(testgenerator.getTau2().E());
    h_EtaISR.Fill(testgenerator.getISR().Eta());
    h_CthISR.Fill(testgenerator.getISR().CosTheta());
    h_PhiISR.Fill(testgenerator.getISR().Phi());
    h_PtISR.Fill(testgenerator.getISR().Pt());
    h_EtaHiggs.Fill(testgenerator.getHiggs().Eta());
    h_CthHiggs.Fill(testgenerator.getHiggs().CosTheta());
    h_PhiHiggs.Fill(testgenerator.getHiggs().Phi());
    h_PtHiggs.Fill(testgenerator.getHiggs().Pt());
    h_EtaTau1boosted.Fill(testgenerator.getTau1boosted().Eta());
    h_CthTau1boosted.Fill(testgenerator.getTau1boosted().CosTheta());
    h_PhiTau1boosted.Fill(testgenerator.getTau1boosted().Phi());
    h_PtTau1boosted.Fill(testgenerator.getTau1boosted().Pt());
    h_EtaTau2boosted.Fill(testgenerator.getTau2boosted().Eta());
    h_CthTau2boosted.Fill(testgenerator.getTau2boosted().CosTheta());
    h_PhiTau2boosted.Fill(testgenerator.getTau2boosted().Phi());
    h_PtTau2boosted.Fill(testgenerator.getTau2boosted().Pt());
    h_VisFracTau1.Fill(testgenerator.getvisfrac1());
    h_VisFracTau2.Fill(testgenerator.getvisfrac2());
    h_EtaTau1Vis.Fill(testgenerator.getTau1Vis().Eta());
    h_CthTau1Vis.Fill(testgenerator.getTau1Vis().CosTheta());
    h_PhiTau1Vis.Fill(testgenerator.getTau1Vis().Phi());
    h_PtTau1Vis.Fill(testgenerator.getTau1Vis().Pt());
    h_ETau1Vis.Fill(testgenerator.getTau1Vis().E());
    h_EtaTau2Vis.Fill(testgenerator.getTau2Vis().Eta());
    h_CthTau2Vis.Fill(testgenerator.getTau2Vis().CosTheta());
    h_PhiTau2Vis.Fill(testgenerator.getTau2Vis().Phi());
    h_PtTau2Vis.Fill(testgenerator.getTau2Vis().Pt());
    h_ETau2Vis.Fill(testgenerator.getTau2Vis().E());
    h_AbsPtmiss.Fill(testgenerator.getAbsPtMET());
    h_PhiMET.Fill(testgenerator.getPhiMET());
    h_AbsPtmisswithsigma.Fill(testgenerator.getAbsPtMETwithsigma());
    h_PhiMETwithsigma.Fill(testgenerator.getPhiMETwithsigma());
    
    //KinFit:
    
    double mass = testgenerator.getMhiggs();

    //prepare tau objects
    HHFitObjectE* tau1 = new HHFitObjectEConstM(testgenerator.getTau1Vis());//  visible Tau1-component from HHTauTauEventGenerator
    HHFitObjectE* tau2 = new HHFitObjectEConstM(testgenerator.getTau2Vis());//  visible Tau2-component from HHTauTauEventGenerator

    //prepare MET object
    HHFitObjectMET* met = new HHFitObjectMET(TVector2(testgenerator.getMETwithsigma()[0],testgenerator.getMETwithsigma()[1]));//Use Met components from HHTauTauEventGenerator
    //met->setCovMatrix(100,-100,50);// set Covarmatrix with Matrix in HHTauTauEventGenerator
    met->setCovMatrix(testgenerator.getCovarmatrix()[0][0],testgenerator.getCovarmatrix()[1][1],testgenerator.getCovarmatrix()[1][0]); // set Covarmatrix with Matrix inserted in HHTauTauEventGenerator
    
    //prepare composite object: Higgs
    HHFitObject* higgs  = new HHFitObjectComposite(tau1, tau2, met);

   //std::cout << "#####################" << i << "#####################" << std::endl;
    //tau1->print();
    //tau2->print();
    //met->print();
    //higgs->print();
    tau1->setLowerFitLimitE(tau1->getInitial4Vector());
    tau1->setUpperFitLimitE(mass,tau2->getInitial4Vector());
    tau2->setLowerFitLimitE(tau2->getInitial4Vector());
    tau2->setUpperFitLimitE(mass,tau1->getInitial4Vector());

    
    //prepare constraints
    HHFitConstraint* invm = new HHFitConstraintEHardM(tau1, tau2, mass);
    HHFitConstraint* balance = new HHFitConstraint4Vector(higgs, true, true, false, false);
    HHFitConstraint* Likelihood = new HHFitConstraintLikelihood(tau1,pdf1);
    
    //fit
    HHKinFit* singlefit = new HHKinFit();
    singlefit->addFitObjectE(tau1);
    singlefit->addConstraint(invm);
    singlefit->addConstraint(balance);
    singlefit->fit();
    h_FitFinalChi2.Fill(singlefit->getChi2());
    
    //fit with likelihood
    HHKinFit* singlefitliekelihood = new HHKinFit();
    singlefitliekelihood->addFitObjectE(tau1);
    singlefitliekelihood->addConstraint(invm);
    singlefitliekelihood->addConstraint(balance);
    singlefitliekelihood->addConstraint(Likelihood);
    singlefitliekelihood->fit();
    h_FitlikelihoodFinalChi2.Fill(singlefit->getChi2());

    
  }
  TFile controlplots("controlplots.root","RECREATE");
  
  h_invariantmass.Write();
  h_EtaTau1.Write();
  h_CthTau1.Write();
  h_PhiTau1.Write();
  h_PtTau1.Write();
  h_ETau1.Write();
  h_EtaTau2.Write();
  h_CthTau2.Write();
  h_PhiTau2.Write();
  h_PtTau2.Write();
  h_ETau2.Write();
  h_EtaISR.Write();
  h_CthISR.Write();
  h_PhiISR.Write();
  h_PtISR.Write();
  h_EtaHiggs.Write();
  h_CthHiggs.Write();
  h_PhiHiggs.Write();
  h_PtHiggs.Write();
  h_EtaTau1boosted.Write();
  h_CthTau1boosted.Write();
  h_PhiTau1boosted.Write();
  h_PtTau1boosted.Write();
  h_EtaTau2boosted.Write();
  h_CthTau2boosted.Write();
  h_PhiTau2boosted.Write();
  h_PtTau2boosted.Write();
  h_VisFracTau1.Write();
  h_VisFracTau2.Write();
  h_EtaTau1Vis.Write();
  h_CthTau1Vis.Write();
  h_PhiTau1Vis.Write();
  h_PtTau1Vis.Write();
  h_ETau1Vis.Write();
  h_EtaTau2Vis.Write();
  h_CthTau2Vis.Write();
  h_PhiTau2Vis.Write();
  h_PtTau2Vis.Write();
  h_ETau2Vis.Write();
  h_AbsPtmiss.Write();
  h_PhiMET.Write();
  h_AbsPtmisswithsigma.Write();
  h_PhiMETwithsigma.Write();
  h_FitFinalChi2.Write();
  h_FitlikelihoodFinalChi2.Write();
  controlplots.Close();
  
  return(0);
}
