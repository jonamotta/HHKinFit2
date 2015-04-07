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
#include "exceptions/HHEnergyRangeException.h"
#include "exceptions/HHCovarianceMatrixException.h"
#include "TGraph.h"
#include "TCanvas.h"

int main(int argc, char* argv[])
{
  double mass = 125.0;

  //prepare tau objects
  HHFitObjectE* tau1 = new HHFitObjectEConstM(TLorentzVector(71,0,0,80));
  HHFitObjectE* tau2 = new HHFitObjectEConstM(TLorentzVector(0,58,0,60));
  tau1->setLowerFitLimitE(tau1->getInitial4Vector());
  tau1->setUpperFitLimitE(mass,tau2->getInitial4Vector());
  tau2->setLowerFitLimitE(tau2->getInitial4Vector());
  tau2->setUpperFitLimitE(mass,tau1->getInitial4Vector());

  //prepare MET object
  HHFitObjectMET* met = new HHFitObjectMET(TVector2(10,20));
  //met->setCovMatrix(100,-100,50);
  met->setCovMatrix(100,100,0);

  //prepare composite object: Higgs
  HHFitObject* higgs  = new HHFitObjectComposite(tau1, tau2, met);

  //prepare constraints
  HHFitConstraint* invm = new HHFitConstraintEHardM(tau1, tau2, mass);
  HHFitConstraint* balance = new HHFitConstraint4Vector(higgs, true, true, false, false);

  //fit
  HHKinFit* singlefit = new HHKinFit();
  singlefit->addFitObjectE(tau1);
  singlefit->addConstraint(invm);
  singlefit->addConstraint(balance);
  singlefit->fit();
  std::cout << "final chi2: " << singlefit->getChi2() << std::endl;

  TCanvas c("c","c",700,700);
  TGraph gr(singlefit->getChi2Function(100));
  gr.Draw();
  c.Print("mychi2function.pdf");
  c.Print("mychi2function.png");

//
//  try{
//    std::cout << singlefit->getChi2() << std::endl;
//  }
//  catch(const HHCovarianceMatrixException& e){
//    std::cout << e.what() << std::endl;
//    std::cout << "will fix it manually" << std::endl;
//    met->setCovMatrix(100,100,0);
//    std::cout << fit->getChi2() << std::endl;
//  }
//  catch(...){
//    std::cout << "caught an unexpected exception" << std::endl;
//  }

  return(0);
}
