#include "HHTauTauEventGenerator.h"
#include "TFile.h"
#include "TH1D.h"
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
#include "exceptions/HHEnergyRangeException.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TFile.h"
#include "HHFitConstraintLikelihood.h"
#include "TMath.h"
#include "TH2D.h"
#include "TGraph.h"

using namespace HHKinFit2;


int main(){
  
  TF1 PDF1("PDF1","2*x",0,1);
  TF1 PDF2("PDF2","2-2*x",0,1);
  TF1* pdf1=&PDF1;
  TF1* pdf2=&PDF2;
  TMatrixD covarmatrix(2,2);
  covarmatrix[0][0]=130;
  covarmatrix[0][1]=0;
  covarmatrix[1][0]=0;
  covarmatrix[1][1]=130;
  HHTauTauEventGenerator higgsgenerator(PDF1,PDF2,covarmatrix);
  higgsgenerator.setMhiggs(125.7);
  HHTauTauEventGenerator BackroundGenerator(PDF1,PDF1,covarmatrix);
  BackroundGenerator.setMhiggs(91.1876);
  //Histogramms for h-fit
  TH1D h_FitFinalChi2h("h_FitFinalChi2h","The Final chi2 from the KinFit for higgsdecay",50,-5,20);
  TH1D h_FitFinalChi2probh("h_FitFinalChi2probh","The Final chi2 from the KinFit for higgsdecay",20,0,1);
  TH1D h_fracresulutiontau1h("h_fracresulutiontau1h","resulution of the energyfraction from tau1vis from the KinFit for higgsdecay",100,-1,1);
  TH1D h_fracresulutiontau2h("h_fracresulutiontau2h","resulution of the energyfraction from tau2vis from the KinFit for higgsdecay",100,-1,1);
  TH1D h_fracresulutiontau1weightedprobh("h_fracresulutiontau1weightedprobh","weighted(chi2prob) resulution of the energyfraction from tau1vis from the KinFit for higgsdecay",100,-1,1);
  TH1D h_fracresulutiontau2weightedprobh("h_fracresulutiontau2weightedprobh","weighted(chi2prob) resulution of the energyfraction from tau2vis from the KinFit for higgsdecay",100,-1,1);
  TH1D h_fracresulutiontau1weightedpdfh("h_fracresulutiontau1weightedpdfh","weighted(pdf) resulution of the energyfraction from tau1vis from the KinFit for higgsdecay",100,-1,1);
  TH1D h_fracresulutiontau2weightedpdfh("h_fracresulutiontau2weightedpdfh","weighted(pdf) resulution of the energyfraction from tau2vis from the KinFit for higgsdecay",100,-1,1);
  TH1D h_energyresulution1h("h_energyresulution1h","the energyresulution of tau 1 in the KinFit for higgsdecay",100,-1.3,1.3);
  TH1D h_energyresulution2h("h_energyresulution2h","the energyresulution of tau 2 in the KinFit for higgsdecay",100,-1.3,1.3);
  TH1D h_energyresulution1weightedprobh("h_energyresulution1weightedprobh","the weighted(chi2prob) energyresulution of tau 1 in the KinFit for higgsdecay",100,-1.3,1.3);
  TH1D h_energyresulution2weightedprobh("h_energyresulution2weightedprobh","the weighted(chi2prob) energyresulution of tau 2 in the KinFit for higgsdecay",100,-1.3,1.3);
  TH1D h_energyresulution1weightedpdfh("h_energyresulution1weightedpdfh","the weighted(pdf) energyresulution of tau 1 in the KinFit for higgsdecay",100,-1.3,1.3);
  TH1D h_energyresulution2weightedpdfh("h_energyresulution2weightedpdfh","the weighted(pdf) energyresulution of tau 2 in the KinFit for higgsdecay",100,-1.3,1.3);
  //Histogramms for h-fitwithlikelihood
  
  TH1D h_FitFinalChi2hlikelihood("h_FitFinalChi2hlikelihood","The Final chi2 from the likelihood-KinFit for higgsdecay",50,-5,20);
  TH1D h_fracresulutiontau1hlikelihood("h_fracresulutiontau1hlikelihood","resulution of the energyfraction from tau1vis from the likelihood-KinFit for higgsdecay",100,-1,1);
  TH1D h_fracresulutiontau2hlikelihood("h_fracresulutiontau2hlikelihood","resulution of the energyfraction from tau2vis from the likelihood-KinFit for higgsdecay",100,-1,1);
  TH1D h_fracresulutiontau1weightedpdfhlikelihood("h_fracresulutiontau1weightedpdfhlikelihood","weighted(pdf) resulution of the energyfraction from tau1vis from the likelihood-KinFit for higgsdecay",100,-1,1);
  TH1D h_fracresulutiontau2weightedpdfhlikelihood("h_fracresulutiontau2weightedpdfhlikelihood","weighted(pdf) resulution of the energyfraction from tau2vis from the likelihood-KinFit for higgsdecay",100,-1,1);
  TH1D h_energyresulution1hlikelihood("h_energyresulution1hlikelihood","the energyresulution of tau 1 in the likelihood-KinFit for higgsdecay",100,-1.3,1.3);
  TH1D h_energyresulution2hlikelihood("h_energyresulution2hlikelihood","the energyresulution of tau 2 in the likelihood-KinFit for higgsdecay",100,-1.3,1.3);
  TH1D h_energyresulution1weightedpdfhlikelihood("h_energyresulution1weightedpdfhlikelihood","the weighted(pdf) energyresulution of tau 1 in the likelihood-KinFit for higgsdecay",100,-1.3,1.3);
  TH1D h_energyresulution2weightedpdfhlikelihood("h_energyresulution2weightedpdfhlikelihood","the weighted(pdf) energyresulution of tau 2 in the likelihood-KinFit for higgsdecay",100,-1.3,1.3);
  
  //Histogramms for background-fit
  TH1D h_FitFinalChi2b("h_FitFinalChi2b","The Final chi2 from the KinFit for background-events",50,-5,20);
  TH1D h_FitFinalChi2probb("h_FitFinalChi2probb","The Final chi2 from the KinFit for background-events",20,0,1);
  TH1D h_fracresulutiontau1b("h_fracresulutiontau1b","resulution of the energyfraction from tau1vis from the KinFit for background-events",100,-1,1);
  TH1D h_fracresulutiontau2b("h_fracresulutiontau2b","resulution of the energyfraction from tau2vis from the KinFit for background-events",100,-1,1);
  TH1D h_fracresulutiontau1weightedprobb("h_fracresulutiontau1weightedprobb","weighted(chi2prob) resulution of the energyfraction from tau1vis from the KinFit for background-events",100,-1,1);
  TH1D h_fracresulutiontau2weightedprobb("h_fracresulutiontau2weightedprobb","weighted(chi2prob) resulution of the energyfraction from tau2vis from the KinFit for background-events",100,-1,1);
  TH1D h_fracresulutiontau1weightedpdfb("h_fracresulutiontau1weightedpdfb","weighted(pdf) resulution of the energyfraction from tau1vis from the KinFit for background-events",100,-1,1);
  TH1D h_fracresulutiontau2weightedpdfb("h_fracresulutiontau2weightedpdfb","weighted(pdf) resulution of the energyfraction from tau2vis from the KinFit for background-events",100,-1,1);
  TH1D h_energyresulution1b("h_energyresulution1b","the energyresulution of tau 1 in the KinFit for background-events",100,-1.3,1.3);
  TH1D h_energyresulution2b("h_energyresulution2b","the energyresulution of tau 2 in the KinFit for background-events",100,-1.3,1.3);
  TH1D h_energyresulution1weightedprobb("h_energyresulution1weightedprobb","the weighted(chi2prob) energyresulution of tau 1 in the KinFit for background-events",100,-1.3,1.3);
  TH1D h_energyresulution2weightedprobb("h_energyresulution2weightedprobb","the weighted(chi2prob) energyresulution of tau 2 in the KinFit for background-events",100,-1.3,1.3);
  TH1D h_energyresulution1weightedpdfb("h_energyresulution1weightedpdfb","the weighted(pdf) energyresulution of tau 1 in the KinFit for background-events",100,-1.3,1.3);
  TH1D h_energyresulution2weightedpdfb("h_energyresulution2weightedpdfb","the weighted(pdf) energyresulution of tau 2 in the KinFit for background-events",100,-1.3,1.3);
  //Histogramms for background-fitwithlikelihood
  
  TH1D h_FitFinalChi2blikelihood("h_FitFinalChi2blikelihood","The Final chi2 from the likelihood-KinFit for background-events",50,-5,20);
  TH1D h_fracresulutiontau1blikelihood("h_fracresulutiontau1blikelihood","resulution of the energyfraction from tau1vis from the likelihood-KinFit for background-events",100,-1,1);
  TH1D h_fracresulutiontau2blikelihood("h_fracresulutiontau2blikelihood","resulution of the energyfraction from tau2vis from the likelihood-KinFit for background-events",100,-1,1);
  TH1D h_fracresulutiontau1weightedpdfblikelihood("h_fracresulutiontau1weightedpdfblikelihood","weighted(pdf) resulution of the energyfraction from tau1vis from the likelihood-KinFit for background-events",100,-1,1);
  TH1D h_fracresulutiontau2weightedpdfblikelihood("h_fracresulutiontau2weightedpdfblikelihood","weighted(pdf) resulution of the energyfraction from tau2vis from the likelihood-KinFit for background-events",100,-1,1);
  TH1D h_energyresulution1blikelihood("h_energyresulution1blikelihood","the energyresulution of tau 1 in the likelihood-KinFit for background-events",100,-1.3,1.3);
  TH1D h_energyresulution2blikelihood("h_energyresulution2blikelihood","the energyresulution of tau 2 in the likelihood-KinFit for background-events",100,-1.3,1.3);
  TH1D h_energyresulution1weightedpdfblikelihood("h_energyresulution1weightedpdfblikelihood","the weighted(pdf) energyresulution of tau 1 in the likelihood-KinFit for background-events",100,-1.3,1.3);
  TH1D h_energyresulution2weightedpdfblikelihood("h_energyresulution2weightedpdfblikelihood","the weighted(pdf) energyresulution of tau 2 in the likelihood-KinFit for background-events",100,-1.3,1.3);
  
  
  
  for(unsigned int i=0; i<50000; i++){
    try{
      higgsgenerator.generateEvent();
    }
    catch(const HHEnergyRangeException& e){
      //std::cout << e.what() << std::endl;
      i--;
      continue;
    }
    try{
      BackroundGenerator.generateEvent();
    }
    catch(const HHEnergyRangeException& e){
      //std::cout << e.what() << std::endl;
      i--;
      continue;
    }
    //Fitting higgsevents:
    
    //without likelihood---------------------------------------------------------------------------------------------------------------------------
    //KinFit:
    double mass = higgsgenerator.getMhiggs();
    
    //prepare tau objects
    HHFitObjectE* tau1h = new HHFitObjectEConstM(higgsgenerator.getTau1Vis());//  visible Tau1-component from HHTauTauEventGenerator
    HHFitObjectE* tau2h = new HHFitObjectEConstM(higgsgenerator.getTau2Vis());//  visible Tau2-component from HHTauTauEventGenerator
    
    //prepare MET object
    HHFitObjectMET* meth = new HHFitObjectMET(TVector2(higgsgenerator.getMETwithsigma()[0],higgsgenerator.getMETwithsigma()[1]));//Use Met components from HHTauTauEventGenerator
    meth->setCovMatrix(higgsgenerator.getCovarmatrix()[0][0],higgsgenerator.getCovarmatrix()[1][1],higgsgenerator.getCovarmatrix()[1][0]); // set Covarmatrix with Matrix inserted in HHTauTauEventGenerator
    
    //prepare composite object: Higgs
    HHFitObject* higgsh  = new HHFitObjectComposite(tau1h, tau2h, meth);
    
    tau1h->setLowerFitLimitE(tau1h->getInitial4Vector());
    tau1h->setUpperFitLimitE(mass,tau2h->getInitial4Vector());
    tau2h->setLowerFitLimitE(tau2h->getInitial4Vector());
    tau2h->setUpperFitLimitE(mass,tau1h->getInitial4Vector());
    
    
    //prepare constraints
    HHFitConstraint* invmh = new HHFitConstraintEHardM(tau1h, tau2h, mass);
    HHFitConstraint* balanceh = new HHFitConstraint4Vector(higgsh, true, true, false, false);
    
    
    //fit
    HHKinFit* singlefith = new HHKinFit();
    singlefith->addFitObjectE(tau1h);
    singlefith->addConstraint(invmh);
    singlefith->addConstraint(balanceh);
    try {
      singlefith->fit();
    }
    catch(HHEnergyRangeException const& e){
      std::cout << i << std::endl;
      tau1h->print();
      tau2h->print();
      meth->print();
      higgsh->print();
      std::cout << e.what() << std::endl;
      std::cout << higgsgenerator.m_seed << std::endl;
      //throw(e);
      std::cout << "-----------------------------------------------" << std::endl;
      continue;
    }
    // histogramms
    h_FitFinalChi2h.Fill(singlefith->getChi2());
    h_FitFinalChi2probh.Fill(TMath::Prob(singlefith->getChi2(),1));
    
    double  fitfractau1h=tau1h->getInitial4Vector().E()/tau1h->getFit4Vector().E();
    double genfractau1h=higgsgenerator.getvisfrac1();
    double comparefrac1h=(genfractau1h-fitfractau1h)/genfractau1h;
    h_fracresulutiontau1h.Fill(comparefrac1h);
    double  fitfractau2h=tau2h->getInitial4Vector().E()/tau2h->getFit4Vector().E();
    double genfractau2h=higgsgenerator.getvisfrac2();
    double comparefrac2h=(genfractau2h-fitfractau2h)/genfractau2h;
    h_fracresulutiontau2h.Fill(comparefrac2h);
    h_energyresulution1h.Fill((higgsgenerator.getTau1boosted().E()-tau1h->getFit4Vector().E())/higgsgenerator.getTau1boosted().E());
    h_energyresulution2h.Fill((higgsgenerator.getTau2boosted().E()-tau2h->getFit4Vector().E())/higgsgenerator.getTau2boosted().E());
    h_fracresulutiontau1weightedprobh.Fill(comparefrac1h,TMath::Prob(singlefith->getChi2(),1));
    h_fracresulutiontau2weightedprobh.Fill(comparefrac2h,TMath::Prob(singlefith->getChi2(),1));
    h_fracresulutiontau1weightedpdfh.Fill(comparefrac1h,PDF1(comparefrac1h));
    h_fracresulutiontau2weightedpdfh.Fill(comparefrac2h,PDF2(comparefrac2h));
    h_energyresulution1weightedprobh.Fill((higgsgenerator.getTau1boosted().E()-tau1h->getFit4Vector().E())/higgsgenerator.getTau1boosted().E(),TMath::Prob(singlefith->getChi2(),1));
    h_energyresulution2weightedprobh.Fill((higgsgenerator.getTau2boosted().E()-tau2h->getFit4Vector().E())/higgsgenerator.getTau2boosted().E(),TMath::Prob(singlefith->getChi2(),1));
    h_energyresulution1weightedpdfh.Fill((higgsgenerator.getTau1boosted().E()-tau1h->getFit4Vector().E())/higgsgenerator.getTau1boosted().E(),PDF1(comparefrac1h));
    h_energyresulution2weightedpdfh.Fill((higgsgenerator.getTau2boosted().E()-tau2h->getFit4Vector().E())/higgsgenerator.getTau2boosted().E(),PDF2(comparefrac2h));
    

    
    
    //with likelihood-----------------------------------------------------------------------------------------------------------------------------------
    //prepare tau objects
    HHFitObjectE* tau1hlikelihood = new HHFitObjectEConstM(higgsgenerator.getTau1Vis());//  visible Tau1-component from HHTauTauEventGenerator
    HHFitObjectE* tau2hlikelihood = new HHFitObjectEConstM(higgsgenerator.getTau2Vis());//  visible Tau2-component from HHTauTauEventGenerator
    
    //prepare MET object
    HHFitObjectMET* methlikelihood = new HHFitObjectMET(TVector2(higgsgenerator.getMETwithsigma()[0],higgsgenerator.getMETwithsigma()[1]));//Use Met components from HHTauTauEventGenerator
    meth->setCovMatrix(higgsgenerator.getCovarmatrix()[0][0],higgsgenerator.getCovarmatrix()[1][1],higgsgenerator.getCovarmatrix()[1][0]); // set Covarmatrix with Matrix inserted in HHTauTauEventGenerator
    
    //prepare composite object: Higgs
    HHFitObject* higgshlikelihood  = new HHFitObjectComposite(tau1h, tau2h, meth);
    
    tau1hlikelihood->setLowerFitLimitE(tau1hlikelihood->getInitial4Vector());
    tau1hlikelihood->setUpperFitLimitE(mass,tau2hlikelihood->getInitial4Vector());
    tau2hlikelihood->setLowerFitLimitE(tau2hlikelihood->getInitial4Vector());
    tau2hlikelihood->setUpperFitLimitE(mass,tau1hlikelihood->getInitial4Vector());
    
    //prepare constraints
    HHFitConstraint* invmhlikelihood = new HHFitConstraintEHardM(tau1hlikelihood, tau2hlikelihood, mass);
    HHFitConstraint* balancehlikelihood = new HHFitConstraint4Vector(higgshlikelihood, true, true, false, false);
    HHFitConstraint* Likelihoodh = new HHFitConstraintLikelihood(tau1hlikelihood,tau2hlikelihood,pdf1,pdf2);
    
    //fit
    HHKinFit* singlefithlikelihood = new HHKinFit();
    singlefithlikelihood->addFitObjectE(tau1hlikelihood);
    singlefithlikelihood->addConstraint(invmhlikelihood);
    singlefithlikelihood->addConstraint(balancehlikelihood);
    singlefithlikelihood->addConstraint(Likelihoodh);
    try {
      singlefithlikelihood->fit();
    }
    catch(HHEnergyRangeException const& e){
      std::cout << i << std::endl;
      tau1hlikelihood->print();
      tau2hlikelihood->print();
      methlikelihood->print();
      higgshlikelihood->print();
      std::cout << e.what() << std::endl;
      std::cout << higgsgenerator.m_seed << std::endl;
      //throw(e);
      std::cout << "-----------------------------------------------" << std::endl;
      continue;
    }
    // histogramms
    h_FitFinalChi2hlikelihood.Fill(singlefithlikelihood->getChi2());
    
    double  fitfractau1hlikelihood=tau1hlikelihood->getInitial4Vector().E()/tau1hlikelihood->getFit4Vector().E();
    double genfractau1hlikelihood=higgsgenerator.getvisfrac1();
    double comparefrac1hlikelihood=(genfractau1hlikelihood-fitfractau1hlikelihood)/genfractau1hlikelihood;
    h_fracresulutiontau1hlikelihood.Fill(comparefrac1hlikelihood);
    double  fitfractau2hlikelihood=tau2hlikelihood->getInitial4Vector().E()/tau2hlikelihood->getFit4Vector().E();
    double genfractau2hlikelihood=higgsgenerator.getvisfrac2();
    double comparefrac2hlikelihood=(genfractau2hlikelihood-fitfractau2hlikelihood)/genfractau2hlikelihood;
    h_fracresulutiontau2hlikelihood.Fill(comparefrac2hlikelihood);
    h_energyresulution1hlikelihood.Fill((higgsgenerator.getTau1boosted().E()-tau1hlikelihood->getFit4Vector().E())/higgsgenerator.getTau1boosted().E());
    h_energyresulution2hlikelihood.Fill((higgsgenerator.getTau2boosted().E()-tau2hlikelihood->getFit4Vector().E())/higgsgenerator.getTau2boosted().E());
    h_fracresulutiontau1weightedpdfhlikelihood.Fill(comparefrac1hlikelihood,PDF1(comparefrac1hlikelihood));
    h_fracresulutiontau2weightedpdfhlikelihood.Fill(comparefrac2hlikelihood,PDF2(comparefrac2hlikelihood));
    h_energyresulution1weightedpdfhlikelihood.Fill((higgsgenerator.getTau1boosted().E()-tau1hlikelihood->getFit4Vector().E())/higgsgenerator.getTau1boosted().E(),PDF1(comparefrac1hlikelihood));
    h_energyresulution2weightedpdfhlikelihood.Fill((higgsgenerator.getTau2boosted().E()-tau2hlikelihood->getFit4Vector().E())/higgsgenerator.getTau2boosted().E(),PDF2(comparefrac2hlikelihood));
    //--------------------------------------------------------------------------------------------------------------------------------------
    //Fitting background-events:
    
    //without likelihood---------------------------------------------------------------------------------------------------------------------------
    //KinFit:
    
    
    //prepare tau objects
    HHFitObjectE* tau1b = new HHFitObjectEConstM(BackroundGenerator.getTau1Vis());//  visible Tau1-component from HHTauTauEventGenerator
    HHFitObjectE* tau2b = new HHFitObjectEConstM(BackroundGenerator.getTau2Vis());//  visible Tau2-component from HHTauTauEventGenerator
    
    //prepare MET object
    HHFitObjectMET* metb = new HHFitObjectMET(TVector2(BackroundGenerator.getMETwithsigma()[0],BackroundGenerator.getMETwithsigma()[1]));//Use Met components from HHTauTauEventGenerator
    metb->setCovMatrix(BackroundGenerator.getCovarmatrix()[0][0],BackroundGenerator.getCovarmatrix()[1][1],BackroundGenerator.getCovarmatrix()[1][0]); // set Covarmatrix with Matrix inserted in HHTauTauEventGenerator
    
    //prepare composite object: Higgs
    HHFitObject* higgsb  = new HHFitObjectComposite(tau1b, tau2b, metb);
    
    tau1b->setLowerFitLimitE(tau1b->getInitial4Vector());
    tau1b->setUpperFitLimitE(mass,tau2b->getInitial4Vector());
    tau2b->setLowerFitLimitE(tau2b->getInitial4Vector());
    tau2b->setUpperFitLimitE(mass,tau1b->getInitial4Vector());

    
    //prepare constraints
    HHFitConstraint* invmb = new HHFitConstraintEHardM(tau1b, tau2b, mass);
    HHFitConstraint* balanceb = new HHFitConstraint4Vector(higgsb, true, true, false, false);
    
    
    //fit
    HHKinFit* singlefitb = new HHKinFit();
    singlefitb->addFitObjectE(tau1b);
    singlefitb->addConstraint(invmb);
    singlefitb->addConstraint(balanceb);
    try {
      singlefitb->fit();
    }
    catch(HHEnergyRangeException const& e){
      std::cout << i << std::endl;
      tau1b->print();
      tau2b->print();
      metb->print();
      higgsb->print();
      std::cout << e.what() << std::endl;
      std::cout << BackroundGenerator.m_seed << std::endl;
      //throw(e);
      std::cout << "-----------------------------------------------" << std::endl;
      continue;
    }
    // histogramms
    h_FitFinalChi2b.Fill(singlefitb->getChi2());
    h_FitFinalChi2probb.Fill(TMath::Prob(singlefitb->getChi2(),1));
    
    double  fitfractau1b=tau1b->getInitial4Vector().E()/tau1b->getFit4Vector().E();
    double genfractau1b=BackroundGenerator.getvisfrac1();
    double comparefrac1b=(genfractau1b-fitfractau1b)/genfractau1b;
    h_fracresulutiontau1b.Fill(comparefrac1b);
    double  fitfractau2b=tau2b->getInitial4Vector().E()/tau2b->getFit4Vector().E();
    double genfractau2b=BackroundGenerator.getvisfrac2();
    double comparefrac2b=(genfractau2b-fitfractau2b)/genfractau2b;
    h_fracresulutiontau2b.Fill(comparefrac2b);
    h_energyresulution1b.Fill((BackroundGenerator.getTau1boosted().E()-tau1b->getFit4Vector().E())/BackroundGenerator.getTau1boosted().E());
    h_energyresulution2b.Fill((BackroundGenerator.getTau2boosted().E()-tau2b->getFit4Vector().E())/BackroundGenerator.getTau2boosted().E());
    h_fracresulutiontau1weightedprobb.Fill(comparefrac1b,TMath::Prob(singlefitb->getChi2(),1));
    h_fracresulutiontau2weightedprobb.Fill(comparefrac2b,TMath::Prob(singlefitb->getChi2(),1));
    h_fracresulutiontau1weightedpdfb.Fill(comparefrac1b,PDF1(comparefrac1b));
    h_fracresulutiontau2weightedpdfb.Fill(comparefrac2b,PDF2(comparefrac2b));
    h_energyresulution1weightedprobb.Fill((BackroundGenerator.getTau1boosted().E()-tau1b->getFit4Vector().E())/BackroundGenerator.getTau1boosted().E(),TMath::Prob(singlefitb->getChi2(),1));
    h_energyresulution2weightedprobb.Fill((BackroundGenerator.getTau2boosted().E()-tau2b->getFit4Vector().E())/BackroundGenerator.getTau2boosted().E(),TMath::Prob(singlefitb->getChi2(),1));
    h_energyresulution1weightedpdfb.Fill((BackroundGenerator.getTau1boosted().E()-tau1b->getFit4Vector().E())/BackroundGenerator.getTau1boosted().E(),PDF1(comparefrac1b));
    h_energyresulution2weightedpdfb.Fill((BackroundGenerator.getTau2boosted().E()-tau2b->getFit4Vector().E())/BackroundGenerator.getTau2boosted().E(),PDF2(comparefrac2b));
    
    
    
    
    //with likelihood-----------------------------------------------------------------------------------------------------------------------------------
    //prepare tau objects
    HHFitObjectE* tau1blikelihood = new HHFitObjectEConstM(BackroundGenerator.getTau1Vis());//  visible Tau1-component from HHTauTauEventGenerator
    HHFitObjectE* tau2blikelihood = new HHFitObjectEConstM(BackroundGenerator.getTau2Vis());//  visible Tau2-component from HHTauTauEventGenerator
    
    //prepare MET object
    HHFitObjectMET* metblikelihood = new HHFitObjectMET(TVector2(BackroundGenerator.getMETwithsigma()[0],BackroundGenerator.getMETwithsigma()[1]));//Use Met components from HHTauTauEventGenerator
    metb->setCovMatrix(BackroundGenerator.getCovarmatrix()[0][0],BackroundGenerator.getCovarmatrix()[1][1],BackroundGenerator.getCovarmatrix()[1][0]); // set Covarmatrix with Matrix inserted in HHTauTauEventGenerator
    
    //prepare composite object: Higgs
    HHFitObject* higgsblikelihood  = new HHFitObjectComposite(tau1b, tau2b, metb);
    
    tau1blikelihood->setLowerFitLimitE(tau1blikelihood->getInitial4Vector());
    tau1blikelihood->setUpperFitLimitE(mass,tau2blikelihood->getInitial4Vector());
    tau2blikelihood->setLowerFitLimitE(tau2blikelihood->getInitial4Vector());
    tau2blikelihood->setUpperFitLimitE(mass,tau1blikelihood->getInitial4Vector());
    
    //prepare constraints
    HHFitConstraint* invmblikelihood = new HHFitConstraintEHardM(tau1blikelihood, tau2blikelihood, mass);
    HHFitConstraint* balanceblikelihood = new HHFitConstraint4Vector(higgsblikelihood, true, true, false, false);
    HHFitConstraint* Likelihoodb = new HHFitConstraintLikelihood(tau1blikelihood,tau2blikelihood,pdf1,pdf2);
    
    //fit
    HHKinFit* singlefitblikelihood = new HHKinFit();
    singlefitblikelihood->addFitObjectE(tau1blikelihood);
    singlefitblikelihood->addConstraint(invmblikelihood);
    singlefitblikelihood->addConstraint(balanceblikelihood);
    singlefitblikelihood->addConstraint(Likelihoodb);
    try {
      singlefitblikelihood->fit();
    }
    catch(HHEnergyRangeException const& e){
      std::cout << i << std::endl;
      tau1blikelihood->print();
      tau2blikelihood->print();
      metblikelihood->print();
      higgsblikelihood->print();
      std::cout << e.what() << std::endl;
      std::cout << BackroundGenerator.m_seed << std::endl;
      //throw(e);
      std::cout << "-----------------------------------------------" << std::endl;
      continue;
    }
    // histogramms
    h_FitFinalChi2blikelihood.Fill(singlefitblikelihood->getChi2());
    
    double  fitfractau1blikelihood=tau1blikelihood->getInitial4Vector().E()/tau1blikelihood->getFit4Vector().E();
    double genfractau1blikelihood=BackroundGenerator.getvisfrac1();
    double comparefrac1blikelihood=(genfractau1blikelihood-fitfractau1blikelihood)/genfractau1blikelihood;
    h_fracresulutiontau1blikelihood.Fill(comparefrac1blikelihood);
    double  fitfractau2blikelihood=tau2blikelihood->getInitial4Vector().E()/tau2blikelihood->getFit4Vector().E();
    double genfractau2blikelihood=BackroundGenerator.getvisfrac2();
    double comparefrac2blikelihood=(genfractau2blikelihood-fitfractau2blikelihood)/genfractau2blikelihood;
    h_fracresulutiontau2blikelihood.Fill(comparefrac2blikelihood);
    h_energyresulution1blikelihood.Fill((BackroundGenerator.getTau1boosted().E()-tau1blikelihood->getFit4Vector().E())/BackroundGenerator.getTau1boosted().E());
    h_energyresulution2blikelihood.Fill((BackroundGenerator.getTau2boosted().E()-tau2blikelihood->getFit4Vector().E())/BackroundGenerator.getTau2boosted().E());
    h_fracresulutiontau1weightedpdfblikelihood.Fill(comparefrac1blikelihood,PDF1(comparefrac1blikelihood));
    h_fracresulutiontau2weightedpdfblikelihood.Fill(comparefrac2blikelihood,PDF2(comparefrac2blikelihood));
    h_energyresulution1weightedpdfblikelihood.Fill((BackroundGenerator.getTau1boosted().E()-tau1blikelihood->getFit4Vector().E())/BackroundGenerator.getTau1boosted().E(),PDF1(comparefrac1blikelihood));
    h_energyresulution2weightedpdfblikelihood.Fill((BackroundGenerator.getTau2boosted().E()-tau2blikelihood->getFit4Vector().E())/BackroundGenerator.getTau2boosted().E(),PDF2(comparefrac2blikelihood));
    //--------------------------------------------------------------------------------------------------------------------------------------
  }
  
  
  TFile backgroundtest("backgroundtest.root","RECREATE");
  
  h_FitFinalChi2h.Write();
  h_FitFinalChi2probh.Write();
  h_fracresulutiontau1h.Write();
  h_fracresulutiontau2h.Write();
  h_fracresulutiontau1weightedprobh.Write();
  h_fracresulutiontau2weightedprobh.Write();
  h_fracresulutiontau1weightedpdfh.Write();
  h_fracresulutiontau2weightedpdfh.Write();
  h_energyresulution1h.Write();
  h_energyresulution2h.Write();
  h_energyresulution1weightedprobh.Write();
  h_energyresulution2weightedprobh.Write();
  h_energyresulution1weightedpdfh.Write();
  h_energyresulution2weightedpdfh.Write();
  //Histogramms for h-fitwithlikelihood
  
  h_FitFinalChi2hlikelihood.Write();
  h_fracresulutiontau1hlikelihood.Write();
  h_fracresulutiontau2hlikelihood.Write();
  h_fracresulutiontau1weightedpdfhlikelihood.Write();
  h_fracresulutiontau2weightedpdfhlikelihood.Write();
  h_energyresulution1hlikelihood.Write();
  h_energyresulution2hlikelihood.Write();
  h_energyresulution1weightedpdfhlikelihood.Write();
  h_energyresulution2weightedpdfhlikelihood.Write();
  
  //Histogramms for background-fit
  h_FitFinalChi2b.Write();
  h_FitFinalChi2probb.Write();
  h_fracresulutiontau1b.Write();
  h_fracresulutiontau2b.Write();
  h_fracresulutiontau1weightedprobb.Write();
  h_fracresulutiontau2weightedprobb.Write();
  h_fracresulutiontau1weightedpdfb.Write();
  h_fracresulutiontau2weightedpdfb.Write();
  h_energyresulution1b.Write();
  h_energyresulution2b.Write();
  h_energyresulution1weightedprobb.Write();
  h_energyresulution2weightedprobb.Write();
  h_energyresulution1weightedpdfb.Write();
  h_energyresulution2weightedpdfb.Write();
  //Histogramms for background-fitwithlikelihood
  
  h_FitFinalChi2blikelihood.Write();
  h_fracresulutiontau1blikelihood.Write();
  h_fracresulutiontau2blikelihood.Write();
  h_fracresulutiontau1weightedpdfblikelihood.Write();
  h_fracresulutiontau2weightedpdfblikelihood.Write();
  h_energyresulution2blikelihood.Write();
  h_energyresulution1weightedpdfblikelihood.Write();
  h_energyresulution2weightedpdfblikelihood.Write();
  
  backgroundtest.Close();
  
  return(0);
}
