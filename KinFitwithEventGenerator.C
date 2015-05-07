/*
 * This is a little program to test the HHTauTauEventGenerator
 * and combine the data from the generated Events with the HHKinFit.
 * Also testing the implementation of an Likelihood-Constraint
 * in the KinFit.
 */



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
#include "TRandom3.h"
#include <cmath>


using namespace HHKinFit2;

int main(int argc, char* argv[])
{

  //preparing objects for the HHTauTauEventGenerator

  TF1 *PDF1 = new TF1("PDF1","2*x",0,1);
  //TF1 *PDF1 = new TF1("PDF1","sqrt(x)",0,2);
  PDF1->SetNpx(100000);
  TF1 *PDF2= new TF1("PDF2","2-2*x",0,1);
  //TF1 *PDF2=new TF1("PDF2","1-sqrt(1-x)",0,2);
  PDF2->SetNpx(100000);
  TMatrixD covarmatrix(2,2);
  covarmatrix[0][0]=130;
  covarmatrix[0][1]=0;
  covarmatrix[1][0]=0;
  covarmatrix[1][1]=130;
  HHTauTauEventGenerator testgenerator(PDF1,PDF2,covarmatrix);
  //-------------------------------------------------------------------------
  //-------------------------------------------------------------------------
  // Define Histograms for observing the eventgenerator and the kinfit
  // Histos for the Generator data:
  TH1D h_invariantmass("h_invariantmass","invariant mass of the two tau-vectors in GeV",100,125.5,126);
  TH1D h_EtaTau1("h_EtaTau1","Eta from Tau1",100,-5,5);
  TH1D h_CthTau1("h_CthTau1","Cth from Tau1",100,-1,1);
  TH1D h_PhiTau1("h_PhiTau1","Phi from Tau1",100,-3.2,3.2);
  TH1D h_PtTau1("h_PtTau1","Pt from Tau1",100,0,100);
  TH1D h_ETau1("h_ETau1","Energy from Tau1",100,0,100);
  TH1D h_EtaTau2("h_EtaTau2","Eta from Tau2",100,-5,5);
  TH1D h_CthTau2("h_CthTau2","Cth from Tau2",100,-1,1);
  TH1D h_PhiTau2("h_PhiTau2","Phi from Tau2",100,-3.2,3.2);
  TH1D h_PtTau2("h_PtTau2","Pt from Tau2",100,0,100);
  TH1D h_invmassboosted("h_invmassboosted","invariant mass of the two boosted tau-vectors in GeV",100,125.5,126);
  TH1D h_ETau2("h_ETau2","Energy from Tau2",100,0,100);
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
  TH1D h_ETau1boosted("h_ETau1boosted","Energy from Tau1boosted",100,0,100);
  TH1D h_EtaTau2boosted("h_EtaTau2boosted","Eta from the boosted Tau2",100,-5,5);
  TH1D h_CthTau2boosted("h_CthTau2boosted","Cth from the boosted Tau2",100,-1,1);
  TH1D h_PhiTau2boosted("h_PhiTau2boosted","Phi from the boosted Tau2",100,-3.2,3.2);
  TH1D h_PtTau2boosted("h_PtTau2boosted","Pt from the boosted Tau2",100,0,100);
  TH1D h_ETau2boosted("h_ETau2boosted","Energy from Tau2boosted",100,0,100);
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
  TH1D h_Pxmiss("h_Pxmiss","Missing momentum in x-direction",200,-100,100);
  TH1D h_Pymiss("h_Pymiss","Missing momentum in y-direction",200,-100,100);
  TH1D h_PhiMET("h_PhiMET","Angle Phi of the missing transversal momentum vector",100,0,6.4);
  TH1D h_AbsPtmisswithsigma("h_AbsPtmisswithsigma","Absolute Value of the missing transverse momentum with measurmental errors",100,0,100);
  TH1D h_PhiMETwithsigma("h_PhiMETwithsigma","Angle Phi of the missing transversal momentum vector with measurmental errors",100,0,6.4);
  TH1D h_Pxmisswithsigma("h_Pxmisswithsigma","Missing momentum in x-direction with measurmental errors ",200,-100,100);
  TH1D h_Pymisswithsigma("h_Pymisswithsigma","Missing momentum in y-direction with measurmental errors",200,-100,100);
  TH1D h_EFracTau1Fit("h_EFracTau1Fit","h_EFracTau1Fit",100,0,1);
  TH1D h_EFracTau2Fit("h_EFracTau2Fit","h_EFracTau2Fit",100,0,1);
  TH1D h_EFracTau1Gen("h_EFracTau1Gen","h_EFracTau1Gen",100,0,1);
  TH1D h_EFracTau2Gen("h_EFracTau2Gen","h_EFracTau2Gen",100,0,1);
  //---------------------------------------------------------------------------------------------------------------
  //Histos for the kinfit:
  TH1D h_FitFinalChi2("h_FitFinalChi2","The Final chi2 from the KinFit",100,-6,6);
  TH1D h_FitFinalChi2prob("h_FitFinalChi2prob","The Final chi2 from the KinFit",20,0,1);
  TH1D h_FitlikelihoodFinalChi2("h_FitFinalChi2likelihood","The Final chi2 from the KinFit",50,-5,20);
  TH1D h_FitlikelihoodFinalChi2prob("h_FitFinalChi2likelihoodprob","The Final probability of the chi2 from the KinFit",20,0,1);
  TH1D h_fracresulutiontau1("h_fracresulutiontau1","resulution of the energyfraction from tau1vis",100,-1,1);
  TH1D h_fracresulutiontau2("h_fracresulutiontau2","resulution of the energyfraction from tau2vis",100,-1,1);
  TH1D h_fracresulutiontau1likelihood("h_fracresulutiontau1likelihood","resulution of the energyfraction from tau1vis in likelihood-fit",100,-1,1);
  TH1D h_fracresulutiontau2likelihood("h_fracresulutiontau2likelihood","resulution of the energyfraction from tau2vis in likelihood-fit",100,-1,1);
  TH1D h_energyresulution1("h_energyresulution1","the energyresulution of tau 1 in the KinFit",100,-1.3,1.3);
  TH1D h_energyresulution1likelihood("h_energyresulution1likelihood","the energyresulution of tau 1 in the Likelihood-KinFit",100,-1.3,1.3);
  TH1D h_energyresulution2("h_energyresulution2","the energyresulution of tau 2 in the KinFit",100,-1.3,1.3);
  TH1D h_energyresulution2likelihood("h_energyresulution2likelihood","the energyresulution of tau 2 in the Likelihood-KinFit",100,-1.3,1.3);
  TH1D h_fracresulutiontau1weighted("h_fracresulutiontau1weighted","weighted resulution of the energyfraction from tau1vis",100,-1,1);
  TH1D h_fracresulutiontau2weighted("h_fracresulutiontau2weighted","weighted resulution of the energyfraction from tau2vis",100,-1,1);
  TH1D h_energyresulution1weighted("h_energyresulutiontau1weighted","the weighted energyresulution of tau 1 in the KinFit",100,-1.3,1.3);
  TH1D h_energyresulution2weighted("h_energyresulutiontau2weighted","the weighted energyresulution of tau 2 in the KinFit",100,-1.3,1.3);
  TH1D h_testingfraclikelihood("h_testingfraclikelihood","visfracs in likelihood-fit",100,0,1.5);
  TH1D h_balanceconstraintwithlikelihood("h_balanceconstraintwithlikelihood","the chi2term for the balanceconstraint",50,-5,20);
  TH1D h_likelihoodconstraint("h_likelihoodconstraint","the chi2term for the likelihoodconstraint",100,-2,2);

  TH2D h_comparemine1("h_comparemine1","blub",100,0,1500,100,0,1500);
  TH2D h_L2D("h_L2D","h_L2D",100,0,0.005,100,0,0.005);





  //------------------------------------------------------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------------------------------------------------------
  // Generator and KinFit-Loop:
  for(unsigned int i=0; i<1000000; i++){
    try{
      testgenerator.generateEvent();
    }
    catch(const HHEnergyRangeException& e){
      i--;
      continue;
    }


    //Filling the Histograms for the eventgenerator
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
    h_ETau1boosted.Fill(testgenerator.getTau1boosted().E());
    h_EtaTau2boosted.Fill(testgenerator.getTau2boosted().Eta());
    h_CthTau2boosted.Fill(testgenerator.getTau2boosted().CosTheta());
    h_PhiTau2boosted.Fill(testgenerator.getTau2boosted().Phi());
    h_PtTau2boosted.Fill(testgenerator.getTau2boosted().Pt());
    h_ETau2boosted.Fill(testgenerator.getTau2boosted().E());
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
    h_Pxmiss.Fill(testgenerator.getMET()[0]);
    h_Pymiss.Fill(testgenerator.getMET()[1]);
    h_PhiMET.Fill(testgenerator.getPhiMET());
    h_AbsPtmisswithsigma.Fill(testgenerator.getAbsPtMETwithsigma());
    h_Pxmisswithsigma.Fill(testgenerator.getMETwithsigma()[0]);
    h_Pymisswithsigma.Fill(testgenerator.getMETwithsigma()[1]);
    h_PhiMETwithsigma.Fill(testgenerator.getPhiMETwithsigma());
    HHLorentzVector boostedsum= (testgenerator.getTau1boosted()+testgenerator.getTau2boosted());
    h_invmassboosted.Fill(boostedsum.M());


 //-----------------------------------------------------------------------------------------------------------------------------------
 //normal KinFit:
    
    //double mass = testgenerator.getMhiggs();
    double mass = 91;

    //prepare tau objects
    HHFitObjectE* tau1 = new HHFitObjectEConstM(testgenerator.getTau1Vis());//  visible Tau1-component from HHTauTauEventGenerator
    HHFitObjectE* tau2 = new HHFitObjectEConstM(testgenerator.getTau2Vis());//  visible Tau2-component from HHTauTauEventGenerator

    //prepare MET object
    HHFitObjectMET* met = new HHFitObjectMET(TVector2(testgenerator.getMETwithsigma()[0],testgenerator.getMETwithsigma()[1]));//Use Met components from HHTauTauEventGenerator
    met->setCovMatrix(testgenerator.getCovarmatrix()[0][0],testgenerator.getCovarmatrix()[1][1],testgenerator.getCovarmatrix()[1][0]); // set Covarmatrix with Matrix inserted in HHTauTauEventGenerator

    //prepare composite object: Higgs
    HHFitObject* higgs  = new HHFitObjectComposite(tau1, tau2, met);

   tau1->setLowerFitLimitE(tau1->getInitial4Vector());
   tau1->setUpperFitLimitE(mass,tau2->getInitial4Vector());
   tau2->setLowerFitLimitE(tau2->getInitial4Vector());
   tau2->setUpperFitLimitE(mass,tau1->getInitial4Vector());

    //prepare constraints
    HHFitConstraint* invm = new HHFitConstraintEHardM(tau1, tau2, mass);
    HHFitConstraint* balance = new HHFitConstraint4Vector(higgs, true, true, false, false);
    
    //fit
    HHKinFit* singlefit = new HHKinFit();
    singlefit->addFitObjectE(tau1);
    singlefit->addConstraint(invm);
    singlefit->addConstraint(balance);

    try {
    singlefit->fit();
    if (!((singlefit->getConvergence()==1)||(singlefit->getConvergence()==2))) continue;
    }
    catch(HHEnergyRangeException const& e){
    	std::cout << i << std::endl;
    	tau1->print();
    	tau2->print();
    	met->print();
    	higgs->print();
    	std::cout << e.what() << std::endl;
    	std::cout << testgenerator.m_seed << std::endl;
    	//throw(e);
    	std::cout << "-----------------------------------------------" << std::endl;
    	//continue;
    }

    //filling Fit-Histos
    h_FitFinalChi2.Fill(singlefit->getChi2());
    h_FitFinalChi2prob.Fill(TMath::Prob(singlefit->getChi2(),1));
	double  fitfractau1=tau1->getInitial4Vector().E()/tau1->getFit4Vector().E();
    double genfractau1=testgenerator.getvisfrac1();
    double comparefrac1=(genfractau1-fitfractau1)/genfractau1;
    h_fracresulutiontau1.Fill(comparefrac1);
    double  fitfractau2=tau2->getInitial4Vector().E()/tau2->getFit4Vector().E();
    double genfractau2=testgenerator.getvisfrac2();
    double comparefrac2=(genfractau2-fitfractau2)/genfractau2;
    h_fracresulutiontau2.Fill(comparefrac2);
    h_energyresulution1.Fill((testgenerator.getTau1boosted().E()-tau1->getFit4Vector().E())/testgenerator.getTau1boosted().E());
    h_energyresulution2.Fill((testgenerator.getTau2boosted().E()-tau2->getFit4Vector().E())/testgenerator.getTau2boosted().E());
    h_fracresulutiontau1weighted.Fill(comparefrac1,TMath::Prob(singlefit->getChi2(),1));
    h_fracresulutiontau2weighted.Fill(comparefrac2,TMath::Prob(singlefit->getChi2(),1));
    h_energyresulution1weighted.Fill((testgenerator.getTau1boosted().E()-tau1->getFit4Vector().E())/testgenerator.getTau1boosted().E(),TMath::Prob(singlefit->getChi2(),1));
    h_energyresulution2weighted.Fill((testgenerator.getTau2boosted().E()-tau2->getFit4Vector().E())/testgenerator.getTau2boosted().E(),TMath::Prob(singlefit->getChi2(),1));
    h_EFracTau1Fit.Fill(tau1->getInitial4Vector().E()/tau1->getFit4Vector().E());
    h_EFracTau2Fit.Fill(tau2->getInitial4Vector().E()/tau2->getFit4Vector().E());
    h_EFracTau1Gen.Fill(testgenerator.getTau1Vis().E()/testgenerator.getTau1boosted().E());
    h_EFracTau2Gen.Fill(testgenerator.getTau2Vis().E()/testgenerator.getTau2boosted().E());

    if(tau2->getInitial4Vector().E()/tau2->getFit4Vector().E()==1)
    {
    TCanvas *chi2canvas=new TCanvas("chi2canvas","chi2canvas",1000,1000);
    TGraph* g = singlefit->getChi2Function(100);
    g->Draw("APL");
    chi2canvas->Print("chi2.pdf");
    chi2canvas->Close();
    singlefit->getChi2();
    break;
    };

    //---------------------------------------------------------------------------------------------------------------------
	//#######################testing#likelihood########################################

    //prepare tau objects
       HHFitObjectE* tau1likelihood = new HHFitObjectEConstM(testgenerator.getTau1Vis());//  visible Tau1-component from HHTauTauEventGenerator
       HHFitObjectE* tau2likelihood = new HHFitObjectEConstM(testgenerator.getTau2Vis());//  visible Tau2-component from HHTauTauEventGenerator

       //prepare MET object
       HHFitObjectMET* metlikelihood = new HHFitObjectMET(TVector2(testgenerator.getMETwithsigma()[0],testgenerator.getMETwithsigma()[1]));//Use Met components from HHTauTauEventGenerator
       //met->setCovMatrix(100,-100,50);// set Covarmatrix with Matrix in HHTauTauEventGenerator
       metlikelihood->setCovMatrix(testgenerator.getCovarmatrix()[0][0],testgenerator.getCovarmatrix()[1][1],testgenerator.getCovarmatrix()[1][0]); // set Covarmatrix with Matrix inserted in HHTauTauEventGenerator

       //prepare composite object: Higgs
       HHFitObject* higgslikelihood  = new HHFitObjectComposite(tau1likelihood, tau2likelihood, metlikelihood);

       tau1likelihood->setLowerFitLimitE(tau1likelihood->getInitial4Vector());
       tau1likelihood->setUpperFitLimitE(mass,tau2likelihood->getInitial4Vector());
       tau2likelihood->setLowerFitLimitE(tau2likelihood->getInitial4Vector());
       tau2likelihood->setUpperFitLimitE(mass,tau1likelihood->getInitial4Vector());


       //prepare constraints
       HHFitConstraint* invmlikelihood = new HHFitConstraintEHardM(tau1likelihood, tau2likelihood, mass);
       HHFitConstraint* balancelikelihood = new HHFitConstraint4Vector(higgslikelihood, true, true, false, false);
       //HHFitConstraint* Likelihood = new HHFitConstraintLikelihood(tau1likelihood,tau2likelihood,PDF1,PDF2);
       HHFitConstraint* Likelihood = new HHFitConstraintLikelihood(tau1likelihood,tau2likelihood,PDF1,PDF1);


    //fit with likelihood
    HHKinFit* singlefitliekelihood = new HHKinFit();
    singlefitliekelihood->addFitObjectE(tau1likelihood);
    singlefitliekelihood->addConstraint(invmlikelihood);
    singlefitliekelihood->addConstraint(balancelikelihood);
    singlefitliekelihood->addConstraint(Likelihood);
    try {
    	singlefitliekelihood->fit();
    	 if (!((singlefitliekelihood->getConvergence()==1)||(singlefitliekelihood->getConvergence()==2))) continue;
        }
        catch(HHEnergyRangeException const& e){
        	std::cout << i << std::endl;
        	tau1->print();
        	tau2->print();
        	met->print();
        	higgs->print();
        	std::cout << e.what() << std::endl;
        	std::cout << testgenerator.m_seed << std::endl;
        	throw(e);
        	std::cout << "-----------------------------------------------" << std::endl;
        	continue;
        }
    h_FitlikelihoodFinalChi2.Fill(singlefitliekelihood->getChi2());
    h_FitlikelihoodFinalChi2prob.Fill(TMath::Prob(singlefitliekelihood->getChi2(),1));
    h_balanceconstraintwithlikelihood.Fill(balancelikelihood->getChi2());
	h_likelihoodconstraint.Fill(Likelihood->getChi2());


    double  fitfractau1likelihood=tau1likelihood->getInitial4Vector().E()/tau1likelihood->getFit4Vector().E();
        double genfractau1likelihood=testgenerator.getvisfrac1();
        double comparefrac1likelihood=(genfractau1likelihood-fitfractau1likelihood)/genfractau1likelihood;
        h_fracresulutiontau1likelihood.Fill(comparefrac1likelihood);
        double  fitfractau2likelihood=tau2likelihood->getInitial4Vector().E()/tau2likelihood->getFit4Vector().E();
               double genfractau2likelihood=testgenerator.getvisfrac2();
               double comparefrac2likelihood=(genfractau2likelihood-fitfractau2likelihood)/genfractau2likelihood;
               h_fracresulutiontau2likelihood.Fill(comparefrac2likelihood);
        h_testingfraclikelihood.Fill(fitfractau1likelihood);
        h_energyresulution1likelihood.Fill((testgenerator.getTau1boosted().E()-tau1likelihood->getFit4Vector().E())/testgenerator.getTau1boosted().E());
        h_energyresulution2likelihood.Fill((testgenerator.getTau2boosted().E()-tau2likelihood->getFit4Vector().E())/testgenerator.getTau2boosted().E());


        h_comparemine1.Fill(tau1->getE(),tau1likelihood->getE());
        h_L2D.Fill(singlefit->getL(),singlefitliekelihood->getL());



  }

  //------------------------------------------------------------------------------------------------------------------------------------
  //writing the Histos into a pdf-file
  TFile controlplots("controlplots.root","RECREATE");
  
  h_L2D.Write();
  h_comparemine1.Write();
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
    h_ETau1boosted.Write();
    h_EtaTau2boosted.Write();
    h_CthTau2boosted.Write();
    h_PhiTau2boosted.Write();
    h_PtTau2boosted.Write();
    h_ETau2boosted.Write();
    h_invmassboosted.Write();
    h_VisFracTau1.Write();
    h_VisFracTau2.Write();
    h_EFracTau1Fit.Write();
    h_EFracTau2Fit.Write();
    h_EFracTau1Gen.Write();
    h_EFracTau2Gen.Write();
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
    h_Pxmiss.Write();
    h_Pymiss.Write();
    h_PhiMET.Write();
    h_AbsPtmisswithsigma.Write();
    h_Pxmisswithsigma.Write();
    h_Pymisswithsigma.Write();
    h_PhiMETwithsigma.Write();
    h_FitFinalChi2.Write();
    h_FitFinalChi2prob.Write();
    h_fracresulutiontau1.Write();
    h_fracresulutiontau1weighted.Write();
    h_fracresulutiontau2.Write();
    h_fracresulutiontau2weighted.Write();
    h_energyresulution1.Write();
    h_energyresulution1weighted.Write();
    h_energyresulution2.Write();
    h_energyresulution2weighted.Write();
    h_FitlikelihoodFinalChi2.Write();
    h_FitlikelihoodFinalChi2prob.Write();
    h_fracresulutiontau1likelihood.Write();
    h_fracresulutiontau2likelihood.Write();
    h_energyresulution1likelihood.Write();
    h_energyresulution2likelihood.Write();
    h_testingfraclikelihood.Write();
    h_balanceconstraintwithlikelihood.Write();
    h_likelihoodconstraint.Write();


  controlplots.Close();



  //---------------------------------------------------------------------------------------------------------------------------------
  return(0);
}
