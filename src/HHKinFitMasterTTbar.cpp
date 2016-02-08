#ifdef HHKINFIT2
#include "HHKinFitMasterTTbar.h"
#include "HHKinFit.h"
#include "HHFitObjectEConstM.h"
#include "HHFitObjectEConstBeta.h"
#include "HHFitObjectMET.h"
#include "HHFitObjectComposite.h"
#include "HHFitConstraintLikelihood.h"
#include "HHFitConstraint.h"
#include "HHFitConstraintEHardM.h"
#include "HHFitObjectEConstM.h"
#include "HHFitObjectEConstBeta.h"
#include "HHFitConstraint4Vector.h"
#include "HHFitConstraint4VectorBJet.h"
#include "exceptions/HHEnergyRangeException.h"
#include "exceptions/HHLimitSettingException.h"
#include "exceptions/HHCovarianceMatrixException.h"
#else
#include "HHKinFit2/HHKinFit2/interface/HHKinFitMasterTTbar.h"
#include "HHKinFit2/HHKinFit2/interface/HHKinFit.h"
#include "HHKinFit2/HHKinFit2/interface/HHFitObjectEConstM.h"
#include "HHKinFit2/HHKinFit2/interface/HHFitObjectEConstBeta.h"
#include "HHKinFit2/HHKinFit2/interface/HHFitObjectMET.h"
#include "HHKinFit2/HHKinFit2/interface/HHFitObjectComposite.h"
#include "HHKinFit2/HHKinFit2/interface/HHFitConstraintLikelihood.h"
#include "HHKinFit2/HHKinFit2/interface/HHFitConstraint.h"
#include "HHKinFit2/HHKinFit2/interface/HHFitConstraintEHardM.h"
#include "HHKinFit2/HHKinFit2/interface/HHFitObjectEConstM.h"
#include "HHKinFit2/HHKinFit2/interface/HHFitObjectEConstBeta.h"
#include "HHKinFit2/HHKinFit2/interface/HHFitConstraint4Vector.h"
#include "HHKinFit2/HHKinFit2/interface/HHFitConstraint4VectorBJet.h"
#include "HHKinFit2/HHKinFit2/interface/exceptions/HHEnergyRangeException.h"
#include "HHKinFit2/HHKinFit2/interface/exceptions/HHLimitSettingException.h"
#include "HHKinFit2/HHKinFit2/interface/exceptions/HHCovarianceMatrixException.h"
#endif

#include "TMatrixD.h"
#include "TRandom3.h"

#include <TMath.h>
#include <cmath>
#include <cstdlib>
#include <iterator>
#include <sstream>

HHKinFit2::HHKinFitMasterTTbar::HHKinFitMasterTTbar(TLorentzVector const& bjet1, 
                                                    TLorentzVector const& bjet2, 
                                                    TLorentzVector const& wvis1,
                                                    TLorentzVector const& wvis2,
                                                    TVector2 const& met, 
                                                    TMatrixD const& met_cov, 
                                                    double sigmaEbjet1,
                                                    double sigmaEbjet2)
  :m_MET_COV(TMatrixD(4,4)), m_bjet1_COV(TMatrixD(4,4)), m_bjet2_COV(TMatrixD(4,4))
{
  m_bjet1 = HHLorentzVector(bjet1.Px(), bjet1.Py(), bjet1.Pz(), bjet1.E());
  m_bjet2 = HHLorentzVector(bjet2.Px(), bjet2.Py(), bjet2.Pz(), bjet2.E());
  m_wvis1 = HHLorentzVector(wvis1.Px(), wvis1.Py(), wvis1.Pz(), wvis1.E());  
  m_wvis2 = HHLorentzVector(wvis2.Px(), wvis2.Py(), wvis2.Pz(), wvis2.E());
  
  m_wvis1.SetMkeepE(80.385);
  m_wvis2.SetMkeepE(80.385);

  m_bestMethodFlag = 0;
  m_loopsNeeded = 0;
  m_useAdveancedBJetChi2 = false;

  m_MET = TVector2(met.Px(), met.Py());
  
  m_MET_COV(0,0) = met_cov(0,0);
  m_MET_COV(1,0) = met_cov(1,0);
  m_MET_COV(0,1) = met_cov(0,1);
  m_MET_COV(1,1) = met_cov(1,1);

  if(sigmaEbjet1 >= 0.0) 
    m_sigma_bjet1 = sigmaEbjet1;
  else
    m_sigma_bjet1 = GetPFBJetRes(m_bjet1.Eta(), m_bjet1.Et());
  
  if(sigmaEbjet2 >= 0.0)
    m_sigma_bjet2 = sigmaEbjet2;
  else
    m_sigma_bjet2 = GetPFBJetRes(m_bjet2.Eta(), m_bjet2.Et());

  m_chi2 = 99999.0;

  m_hypo = HHFitHypothesis2M(173, 173);

}

void HHKinFit2::HHKinFitMasterTTbar::fit()
{
  //TCanvas* c1 = new TCanvas();

    HHFitObjectE* w1Fit = new HHFitObjectEConstM(m_wvis1);
    HHFitObjectE* w2Fit = new HHFitObjectEConstM(m_wvis2);

    HHLorentzVector w1min = m_wvis1;
    w1min.SetEkeepM(1.01*m_wvis1.M());      
    HHLorentzVector w2min = m_wvis2;
    w2min.SetEkeepM(1.01*m_wvis2.M());

    HHFitObjectE* b1Fit = new HHFitObjectEConstBeta(m_bjet1);
    HHFitObjectE* b2Fit = new HHFitObjectEConstBeta(m_bjet2);

    HHLorentzVector b1min = m_bjet1;
    b1min.SetEkeepM(5);      
    HHLorentzVector b2min = m_bjet2;
    b2min.SetEkeepM(5);
  

    //prepare MET object
    HHFitObjectMET* metFit = new HHFitObjectMET(m_MET);

    //prepare composite object: ttbarsystem
    HHFitObject* ttbarsystem  = new HHFitObjectComposite(w1Fit, w2Fit, b1Fit, b2Fit, metFit);
    HHFitObject* t1  = new HHFitObjectComposite(b1Fit, w1Fit);  
    HHFitObject* t2  = new HHFitObjectComposite(b2Fit, w2Fit);

    int mt1 = m_hypo.first;
    int mt2 = m_hypo.second;
  
    try
    {
      w1Fit->setFitLimitsE(w1min,mt1,b1min);
      w2Fit->setFitLimitsE(w2min,mt2,b2min);
    }
    catch(HHLimitSettingException const& e)
    {
      std::cout << "Exception while setting tau limits:" << std::endl;
      std::cout << e.what() << std::endl;
      std::cout << "Tau energies are not compatible with invariant mass constraint." << std::endl;

      m_prob = -pow(10,10);
      m_chi2 = -pow(10,10);
      m_convergence = -1;
    }

    if(!m_useAdveancedBJetChi2)
    {
      if(fabs(m_bjet1_COV(0,0)) < 0.001)
      {
	b1Fit->setCovMatrix(m_sigma_bjet1);
	b2Fit->setCovMatrix(m_sigma_bjet2);
      }
      else
      {
	std::cout << "Cov MATRIX set manually for ToyMC study." << std::endl;      
	b1Fit->setCovMatrix(m_bjet1_COV);
	b2Fit->setCovMatrix(m_bjet2_COV);
      }
    }
    

    metFit->setCovMatrix(m_MET_COV);
    ttbarsystem->setCovMatrix(m_MET_COV);// - m_bjet1_COV - m_bjet2_COV);

    //prepare constraints
    HHFitConstraint* c_invmt1 = new HHFitConstraintEHardM(w1Fit, b1Fit, mt1);
    HHFitConstraint* c_invmt2 = new HHFitConstraintEHardM(w2Fit, b2Fit, mt2);    

    HHFitConstraint* c_b1;
    HHFitConstraint* c_b2;
    if(m_useAdveancedBJetChi2)
    {
       c_b1 = new HHFitConstraint4VectorBJet(b1Fit);
       c_b2 = new HHFitConstraint4VectorBJet(b2Fit);
    }
    else
    {
      c_b1 = new HHFitConstraint4Vector(b1Fit, false, false, false, true);
      c_b2 = new HHFitConstraint4Vector(b2Fit, false, false, false, true);
    }

    HHFitConstraint* c_balance = new HHFitConstraint4Vector(ttbarsystem, true, true, false, false);

    //fit
    HHKinFit2::HHKinFit* fitObject = new HHKinFit2::HHKinFit();

    w1Fit->setInitDirection(1.0);
    w1Fit->setInitPrecision(0.1);
    w1Fit->setInitStepWidth(0.1*(w1Fit->getUpperFitLimitE() -  w1Fit->getLowerFitLimitE()));

    w2Fit->setInitDirection(1.0);
    w2Fit->setInitPrecision(0.1);
    w2Fit->setInitStepWidth(0.1*(w2Fit->getUpperFitLimitE() -  w2Fit->getLowerFitLimitE()));
 
    fitObject->addFitObjectE(w1Fit);
    fitObject->addFitObjectE(w2Fit);

    fitObject->addConstraint(c_invmt1);
    fitObject->addConstraint(c_invmt2);
    fitObject->addConstraint(c_b1);
    fitObject->addConstraint(c_b2);
    fitObject->addConstraint(c_balance);

//    double b1NegTauUpFit,b1PosTauDownFit, b1NegTauDownFit;
//    double tau1NegTauUpFit, tau1PosTauDownFit, tau1NegTauDownFit;
    double chi2NegTauUp = 99999;
    double chi2PosTauDown = 99999;
    double chi2NegTauDown = 99999;
    double chi2Min = 99999;

    //--NegTauDown
    w1Fit->setInitStart(w1Fit->getUpperFitLimitE());
    try
    {
      fitObject->fit();
    }
    catch(HHLimitSettingException const& e)
    {
      std::cout << e.what() << std::endl;
      m_chi2 = -pow(10,10);
      m_prob = -pow(10,10);
      m_chi2 = -pow(10,10);
      if((b1Fit->getUpperFitLimitE() - b1Fit->getLowerFitLimitE()) <
	 (w1Fit->getUpperFitLimitE() - w1Fit->getLowerFitLimitE()))
      {
	m_convergence = -2;
      }
      else
      {
	m_convergence = -1;
      }
    }
    catch(HHKinFit2::HHEnergyRangeException const& e){
      std::cout << "Energy Range Exception" << std::endl;
      m_prob = -pow(10,10);
      m_chi2 = -pow(10,10);
      m_convergence = 0;
    }

    if(fitObject->getChi2() > 0)
    {
      chi2NegTauDown = fitObject->getChi2();
      if(chi2NegTauDown < chi2Min)
      {
	chi2Min = chi2NegTauDown;
      }
      //b1NegTauDownFit = b1Fit->getFit4Vector().E();
      //tau1NegTauDownFit = w1Fit->getFit4Vector().E();
    }


    //--NegTauUp
//      b1Fit->setInitStart( 0.9*b1Fit->getInitial4Vector().E());
    w1Fit->setInitStart((w1Fit->getUpperFitLimitE() - 
                         w1Fit->getLowerFitLimitE())/2.0);
    try
    {
      fitObject->fit();
    }
    catch(HHLimitSettingException const& e)
    {
      std::cout << e.what() << std::endl;
      m_prob = -pow(10,10);
      m_chi2 = -pow(10,10);
      if((b1Fit->getUpperFitLimitE() - b1Fit->getLowerFitLimitE()) <
	 (w1Fit->getUpperFitLimitE() - w1Fit->getLowerFitLimitE()))
      {
	m_convergence = -2;
      }
      else
      {
	m_convergence = -1;
      }
    }
    catch(HHKinFit2::HHEnergyRangeException const& e){
      std::cout << "Energy Range Exception" << std::endl;
      m_prob = -pow(10,10);
      m_chi2 = -pow(10,10);
      m_convergence = 0;
    }

    if(fitObject->getChi2() > 0)
    {
      chi2NegTauUp = fitObject->getChi2();
      if(chi2NegTauUp < chi2Min)
      {
	chi2Min = chi2NegTauUp;
      }
      //b1NegTauUpFit = b1Fit->getFit4Vector().E();
      //tau1NegTauUpFit = w1Fit->getFit4Vector().E();
    }

    //--PosTauDown
    w1Fit->setInitStart(w1Fit->getLowerFitLimitE());
    try
    {
      fitObject->fit();
    }
    catch(HHLimitSettingException const& e)
    {
      std::cout << e.what() << std::endl;
      m_prob = -pow(10,10);
      m_chi2 = -pow(10,10);
      if((b1Fit->getUpperFitLimitE() - b1Fit->getLowerFitLimitE()) <
	 (w1Fit->getUpperFitLimitE() - w1Fit->getLowerFitLimitE()))
      {
	m_convergence = -2;
      }
      else
      {
	m_convergence = -1;
      }
    }
    catch(HHKinFit2::HHEnergyRangeException const& e){
      std::cout << "Energy Range Exception" << std::endl;
      m_prob = -pow(10,10);
      m_chi2 = -pow(10,10);
      m_convergence = 0;
    }

    if(fitObject->getChi2() > 0)
    {
      chi2PosTauDown = fitObject->getChi2();
      if(chi2PosTauDown < chi2Min)
      {
	chi2Min = chi2PosTauDown;
      }
      //b1PosTauDownFit = b1Fit->getFit4Vector().E();
      //tau1PosTauDownFit = w1Fit->getFit4Vector().E();
    }
       
    bool improvementDetected = false;
    if(fabs(chi2NegTauDown - chi2Min) > 0.1)
    {
      improvementDetected = true;
    }
    if(fabs(chi2PosTauDown - chi2Min) > 0.1)
    {
      improvementDetected = true;
    }
    if(fabs(chi2NegTauUp - chi2Min) > 0.1)
    {
      improvementDetected = true;
    }

    if(improvementDetected) 
    {
      /*
	std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" 
	<< std::endl;
	std::cout << "Difference in B-Fit directions detected: " << std::endl;
	std::cout << "Chi2 PosTauUp: " << chi2PosTauUp << std::endl
	<< "Chi2 NegTauUp: " << chi2NegTauUp << std::endl
	<< "Chi2 PosTauDown: " << chi2PosTauDown << std::endl
	<< "Chi2 negTauDown: " << chi2NegTauDown << std::endl;
	mins[0] = b1PosTauDownFit;
	maxs[0] = b1PosTauDownFit;
      if(b1NegTauDownFit < mins[0])
      {
	mins[0] = b1NegTauDownFit;
      }
      if(b1NegTauDownFit > maxs[0])
      {
	maxs[0] = b1NegTauDownFit;
      }
      
      if(doTauUpStart)
      {
	if(b1NegTauUpFit < mins[0])
	{
	  mins[0] = b1NegTauUpFit;
	}
	if(b1NegTauUpFit > maxs[0])
	{
	  maxs[0] = b1NegTauUpFit;
	}
      
	if(b1PosTauUpFit < mins[0])
	{
	  mins[0] = b1PosTauUpFit;
	}
	if(b1PosTauUpFit > maxs[0])
	{
	  maxs[0] = b1PosTauUpFit;
	}
      }
      mins[1] = tau1PosTauDownFit;
      maxs[1] = tau1PosTauDownFit;
      if(tau1NegTauDownFit < mins[1])
      {
	mins[1] = tau1NegTauDownFit;
      }
      if(tau1NegTauDownFit > maxs[1])
      {
	maxs[1] = tau1NegTauDownFit;
      }
      
      if(doTauUpStart)
      {
	if(tau1NegTauUpFit < mins[1])
	{
	  mins[1] = tau1NegTauUpFit;
	}
	if(tau1NegTauUpFit > maxs[1])
	{
	  maxs[1] = tau1NegTauUpFit;
	}

	if(tau1PosTauUpFit < mins[1])
	{
	  mins[1] = tau1PosTauUpFit;
	}
	if(tau1PosTauUpFit > maxs[1])
	{
	  maxs[1] = tau1PosTauUpFit;
	}
      }
      mins[0] = mins[0]*0.9;
      mins[1] = mins[1]*0.9;
      maxs[0] = maxs[0]*1.1;
      maxs[1] = maxs[1]*1.1;

      TGraph2D* graph = fitObject->getChi2Function(steps, mins, maxs);

      graph->GetXaxis()->SetTitle("X");
      //graph->GetXaxis()->SetTitleOffset(0.1);
      graph->GetYaxis()->SetTitle("Y");
      //graph->GetYaxis()->SetTitleOffset(0.1);
      graph->GetZaxis()->SetTitle("Z");
      //graph->GetZaxis()->SetTitleOffset(1.0);
      graph->Draw("Cont3COLZ");

      TString b1PosTauUp;
      TString b1NegTauUp;
      TString b1PosTauDown;
      TString b1NegTauDown;
      TString tau1PosTauUp;
      TString tau1NegTauUp;
      TString tau1PosTauDown;
      TString tau1NegTauDown;

      std::stringstream stream;

      if(doTauUpStart)
      {
	stream << b1PosTauUpFit;
	stream >> b1PosTauUp;
	stream.clear();

	stream << b1NegTauUpFit;
	stream >> b1NegTauUp;
	stream.clear();
      }
      stream << b1PosTauDownFit;
      stream >> b1PosTauDown;
      stream.clear();

      stream << b1NegTauDownFit;
      stream >> b1NegTauDown;
      stream.clear();
      
      if(doTauUpStart)
      {
	stream << tau1PosTauUpFit;
	stream >> tau1PosTauUp;
	stream.clear();
	
	stream << tau1NegTauUpFit;
	stream >> tau1NegTauUp;
	stream.clear();
      }
      stream << tau1PosTauDownFit;
      stream >> tau1PosTauDown;
      stream.clear();

      stream << tau1NegTauDownFit;
      stream >> tau1NegTauDown;
      stream.clear();

      if(doTauUpStart)
      {
	c1->SaveAs("PosUp" + b1PosTauUp + "-" + tau1PosTauUp +
		   "NegUp" + b1NegTauUp + "-" + tau1NegTauUp +
		   "PosDown" + b1PosTauDown + "-" + tau1PosTauDown + 
		   "NegDown" + b1NegTauDown + "-" + tau1NegTauDown +
		   ".pdf");
      }
      else
      {
	c1->SaveAs("PosDown" + b1PosTauDown + "-" + tau1PosTauDown + 
		   "NegDown" + b1NegTauDown + "-" + tau1NegTauDown +
		   ".pdf");
      }
      */
      bool hasRerun = false;
      m_bestMethodFlag = 0;
      //Mark all methods that got the best chi2
      //1 = chi2NegTauDown
      //2 = chi2PosTauDown
      //3 = chi2NegTauDown + chi2PosTauDown
      //4 = chi2NegTauUp
      //5 = chi2NegTauUp + chi2NegTauDown
      //6 = chi2NegTauUp + chi2PosTauDown
      //7 = All, should not happen
      
      if(fabs(chi2NegTauDown - chi2Min) <= 0.1)
      {
	m_bestMethodFlag = m_bestMethodFlag | (1 << 0); 
	w1Fit->setInitStart(w1Fit->getUpperFitLimitE());
	fitObject->fit();
	hasRerun = true;
      }
      if(fabs(chi2PosTauDown - chi2Min) <= 0.1)
      {
	m_bestMethodFlag = m_bestMethodFlag | (1 << 1); 
	hasRerun = true; //No need to rerun as this was the last config to run
      }
   
      if(fabs(chi2NegTauUp - chi2Min) <= 0.1)
      {
	m_bestMethodFlag = m_bestMethodFlag | (1 << 2); 
	if(!hasRerun)
	{
	  w1Fit->setInitStart((w1Fit->getUpperFitLimitE() - 
				 w1Fit->getLowerFitLimitE())/2.0);
	  fitObject->fit();
	  hasRerun = true;
	}
      }
      
      if(!hasRerun)
      {
	std::cout << "--------ALARM! ALARM! SOMETHING HAS GONE HORRIBLY WRONG!-------" 
		  << std::endl;
      }
    }

    /* No Convergence MAP
    if(fitObject->getConvergence() == 0)
    {
      std::cout << "####################################################" << std::endl;
      std::cout << "No Convergence! Printing how we got here: " << std::endl;
      fitObject->setPrintLevel(3);
      fitObject->fit();
      fitObject->setPrintLevel(0);
      std::cout << "####################################################" << std::endl;

      mins[0] = 0.9*b1Fit->getFit4Vector().E();
      mins[1] = 0.9*w1Fit->getFit4Vector().E();
      maxs[0] = 1.1*b1Fit->getFit4Vector().E();
      maxs[1] = 1.1*w1Fit->getFit4Vector().E();
 
      TString b1Point;
      TString tau1Point;
      std::stringstream stream;

      stream << b1Fit->getFit4Vector().E();
      stream >> b1Point;
      stream.clear();

      stream << w1Fit->getFit4Vector().E();
      stream >> tau1Point;
      stream.clear();

      TGraph2D* graph = (TGraph2D*)fitObject->getChi2Function(steps, mins, maxs);
      graph->Draw("Cont3COLZ");

      c1->SaveAs("NoConvergence" + b1Point + "-" + tau1Point + ".pdf");
    }
    */

    {//Filling Results
      m_convergence = fitObject->getConvergence();    
    
      double chi2 = fitObject->getChi2();
      m_chi2 = chi2;
      m_prob = TMath::Prob(chi2, 2);
      m_loopsNeeded = fitObject->m_loopsNeeded;

      if(chi2 < m_chi2)
      {
	m_chi2 = chi2;
      }
    
      m_chi2BJet1 = c_b1->getChi2();
      m_chi2BJet2 = c_b2->getChi2();
      m_chi2Balance = c_balance->getChi2();
    
      TLorentzVector fittedW1 =  ( (TLorentzVector)w1Fit->getFit4Vector()  );
      m_fittedW1 = fittedW1;
      TLorentzVector fittedW2 =  ( (TLorentzVector)w2Fit->getFit4Vector()  );
      m_fittedW2 = fittedW2;
      TLorentzVector fittedB1 =  ( (TLorentzVector)b1Fit->getFit4Vector() ) ;
      m_fittedB1 = fittedB1;
      TLorentzVector fittedB2 =  ( (TLorentzVector)b2Fit->getFit4Vector() ) ;
      m_fittedB2 = fittedB2;
    }

    delete c_invmt1;
    delete c_invmt2;
    delete c_b1;
    delete c_b2;
    delete c_balance;
    
    delete fitObject;
   
    delete w1Fit;
    delete w2Fit;
    delete b1Fit;
    delete b2Fit;
    delete metFit;
    delete ttbarsystem;
    delete t1;
    delete t2;
}

void HHKinFit2::HHKinFitMasterTTbar::useAdvancedBJetChi2(bool useAdvancedBJetChi2)
{
  m_useAdveancedBJetChi2 = useAdvancedBJetChi2;
}

//Getters for fit results

double HHKinFit2::HHKinFitMasterTTbar::getChi2(){
  return(m_chi2);
}

  
double HHKinFit2::HHKinFitMasterTTbar::getChi2BJet1(){
  return(m_chi2BJet1);
}


double HHKinFit2::HHKinFitMasterTTbar::getChi2BJet2(){
  return(m_chi2BJet2);
}
    

double HHKinFit2::HHKinFitMasterTTbar::getChi2Balance(){
  return(m_chi2Balance);
}
  
double HHKinFit2::HHKinFitMasterTTbar::getFitProb(){
  return(m_prob);
}
  
int HHKinFit2::HHKinFitMasterTTbar::getConvergence(){
  return(m_convergence);
}

TLorentzVector HHKinFit2::HHKinFitMasterTTbar::getFittedW1(){
  return(m_fittedW1);
}


TLorentzVector HHKinFit2::HHKinFitMasterTTbar::getFittedW2(){
  return(m_fittedW2);
}

TLorentzVector HHKinFit2::HHKinFitMasterTTbar::getFittedBJet1(){
  return(m_fittedB1);
}
  

TLorentzVector HHKinFit2::HHKinFitMasterTTbar::getFittedBJet2(){
  return(m_fittedB2);
}

double HHKinFit2::HHKinFitMasterTTbar::GetPFBJetRes(double eta, double et){
  double det=0;
  double de=10;

  if(0.000<=fabs(eta) && fabs(eta)<0.087){
    det = et * (sqrt(0.0686*0.0686 + (1.03/sqrt(et))*(1.03/sqrt(et)) + (1.68/et)*(1.68/et)));
    de = 1.0/sin(2 * atan(exp(-(0.000+0.087)/2))) * det;
  }

  if(0.087<=fabs(eta) && fabs(eta)<0.174){
    det = et * (sqrt(0.0737*0.0737 + (1.01/sqrt(et))*(1.01/sqrt(et)) + (1.74/et)*(1.74/et)));
    de = 1.0/sin(2 * atan(exp(-(0.087+0.174)/2))) * det;
  }

  if(0.174<=fabs(eta) && fabs(eta)<0.261){
    det = et * (sqrt(0.0657*0.0657 + (1.07/sqrt(et))*(1.07/sqrt(et)) + (5.16e-06/et)*(5.16e-06/et)));
    de = 1.0/sin(2 * atan(exp(-(0.174+0.261)/2))) * det;
  }

  if(0.261<=fabs(eta) && fabs(eta)<0.348){
    det = et * (sqrt(0.062*0.062 + (1.07/sqrt(et))*(1.07/sqrt(et)) + (0.000134/et)*(0.000134/et)));
    de = 1.0/sin(2 * atan(exp(-(0.261+0.348)/2))) * det;
  }

  if(0.348<=fabs(eta) && fabs(eta)<0.435){
    det = et * (sqrt(0.0605*0.0605 + (1.07/sqrt(et))*(1.07/sqrt(et)) + (1.84e-07/et)*(1.84e-07/et)));
    de = 1.0/sin(2 * atan(exp(-(0.348+0.435)/2))) * det;
  }

  if(0.435<=fabs(eta) && fabs(eta)<0.522){
    det = et * (sqrt(0.059*0.059 + (1.08/sqrt(et))*(1.08/sqrt(et)) + (9.06e-09/et)*(9.06e-09/et)));
    de = 1.0/sin(2 * atan(exp(-(0.435+0.522)/2))) * det;
  }

  if(0.522<=fabs(eta) && fabs(eta)<0.609){
    det = et * (sqrt(0.0577*0.0577 + (1.08/sqrt(et))*(1.08/sqrt(et)) + (5.46e-06/et)*(5.46e-06/et)));
    de = 1.0/sin(2 * atan(exp(-(0.522+0.609)/2))) * det;
  }

  if(0.609<=fabs(eta) && fabs(eta)<0.696){
    det = et * (sqrt(0.0525*0.0525 + (1.09/sqrt(et))*(1.09/sqrt(et)) + (4.05e-05/et)*(4.05e-05/et)));
    de = 1.0/sin(2 * atan(exp(-(0.609+0.696)/2))) * det;
  }

  if(0.696<=fabs(eta) && fabs(eta)<0.783){
    det = et * (sqrt(0.0582*0.0582 + (1.09/sqrt(et))*(1.09/sqrt(et)) + (1.17e-05/et)*(1.17e-05/et)));
    de = 1.0/sin(2 * atan(exp(-(0.696+0.783)/2))) * det;
  }

  if(0.783<=fabs(eta) && fabs(eta)<0.870){
    det = et * (sqrt(0.0649*0.0649 + (1.08/sqrt(et))*(1.08/sqrt(et)) + (7.85e-06/et)*(7.85e-06/et)));
    de = 1.0/sin(2 * atan(exp(-(0.783+0.870)/2))) * det;
  }

  if(0.870<=fabs(eta) && fabs(eta)<0.957){
    det = et * (sqrt(0.0654*0.0654 + (1.1/sqrt(et))*(1.1/sqrt(et)) + (1.09e-07/et)*(1.09e-07/et)));
    de = 1.0/sin(2 * atan(exp(-(0.870+0.957)/2))) * det;
  }

  if(0.957<=fabs(eta) && fabs(eta)<1.044){
    det = et * (sqrt(0.0669*0.0669 + (1.11/sqrt(et))*(1.11/sqrt(et)) + (1.87e-06/et)*(1.87e-06/et)));
    de = 1.0/sin(2 * atan(exp(-(0.957+1.044)/2))) * det;
  }

  if(1.044<=fabs(eta) && fabs(eta)<1.131){
    det = et * (sqrt(0.0643*0.0643 + (1.15/sqrt(et))*(1.15/sqrt(et)) + (2.76e-05/et)*(2.76e-05/et)));
    de = 1.0/sin(2 * atan(exp(-(1.044+1.131)/2))) * det;
  }

  if(1.131<=fabs(eta) && fabs(eta)<1.218){
    det = et * (sqrt(0.0645*0.0645 + (1.16/sqrt(et))*(1.16/sqrt(et)) + (1.04e-06/et)*(1.04e-06/et)));
    de = 1.0/sin(2 * atan(exp(-(1.131+1.218)/2))) * det;
  }

  if(1.218<=fabs(eta) && fabs(eta)<1.305){
    det = et * (sqrt(0.0637*0.0637 + (1.19/sqrt(et))*(1.19/sqrt(et)) + (1.08e-07/et)*(1.08e-07/et)));
    de = 1.0/sin(2 * atan(exp(-(1.218+1.305)/2))) * det;
  }

  if(1.305<=fabs(eta) && fabs(eta)<1.392){
    det = et * (sqrt(0.0695*0.0695 + (1.21/sqrt(et))*(1.21/sqrt(et)) + (5.75e-06/et)*(5.75e-06/et)));
    de = 1.0/sin(2 * atan(exp(-(1.305+1.392)/2))) * det;
  }

  if(1.392<=fabs(eta) && fabs(eta)<1.479){
    det = et * (sqrt(0.0748*0.0748 + (1.2/sqrt(et))*(1.2/sqrt(et)) + (5.15e-08/et)*(5.15e-08/et)));
    de = 1.0/sin(2 * atan(exp(-(1.392+1.479)/2))) * det;
  }

  if(1.479<=fabs(eta) && fabs(eta)<1.566){
    det = et * (sqrt(0.0624*0.0624 + (1.23/sqrt(et))*(1.23/sqrt(et)) + (2.28e-05/et)*(2.28e-05/et)));
    de = 1.0/sin(2 * atan(exp(-(1.479+1.566)/2))) * det;
  }

  if(1.566<=fabs(eta) && fabs(eta)<1.653){
    det = et * (sqrt(0.0283*0.0283 + (1.25/sqrt(et))*(1.25/sqrt(et)) + (4.79e-07/et)*(4.79e-07/et)));
    de = 1.0/sin(2 * atan(exp(-(1.566+1.653)/2))) * det;
  }

  if(1.653<=fabs(eta) && fabs(eta)<1.740){
    det = et * (sqrt(0.0316*0.0316 + (1.21/sqrt(et))*(1.21/sqrt(et)) + (5e-05/et)*(5e-05/et)));
    de = 1.0/sin(2 * atan(exp(-(1.653+1.740)/2))) * det;
  }

  if(1.740<=fabs(eta) && fabs(eta)<1.830){
    det = et * (sqrt(2.29e-07*2.29e-07 + (1.2/sqrt(et))*(1.2/sqrt(et)) + (1.71e-05/et)*(1.71e-05/et)));
    de = 1.0/sin(2 * atan(exp(-(1.740+1.830)/2))) * det;
  }

  if(1.830<=fabs(eta) && fabs(eta)<1.930){
    det = et * (sqrt(5.18e-09*5.18e-09 + (1.14/sqrt(et))*(1.14/sqrt(et)) + (1.7/et)*(1.7/et)));
    de = 1.0/sin(2 * atan(exp(-(1.830+1.930)/2))) * det;
  }

  if(1.930<=fabs(eta) && fabs(eta)<2.043){
    det = et * (sqrt(2.17e-07*2.17e-07 + (1.09/sqrt(et))*(1.09/sqrt(et)) + (2.08/et)*(2.08/et)));
    de = 1.0/sin(2 * atan(exp(-(1.930+2.043)/2))) * det;
  }

  if(2.043<=fabs(eta) && fabs(eta)<2.172){
    det = et * (sqrt(3.65e-07*3.65e-07 + (1.09/sqrt(et))*(1.09/sqrt(et)) + (1.63/et)*(1.63/et)));
    de = 1.0/sin(2 * atan(exp(-(2.043+2.172)/2))) * det;
  }

  if(2.172<=fabs(eta) && fabs(eta)<2.322){
    det = et * (sqrt(2.02e-07*2.02e-07 + (1.09/sqrt(et))*(1.09/sqrt(et)) + (1.68/et)*(1.68/et)));
    de = 1.0/sin(2 * atan(exp(-(2.172+2.322)/2))) * det;
  }

  if(2.322<=fabs(eta) && fabs(eta)<2.500){
    det = et * (sqrt(5.27e-07*5.27e-07 + (1.12/sqrt(et))*(1.12/sqrt(et)) + (1.78/et)*(1.78/et)));
    de = 1.0/sin(2 * atan(exp(-(2.322+2.500)/2))) * det;
  }

  return de;

}    

/*double HHKinFit2::crystBallLikePDF(double x, double alpha, double n, double sigma, 
				   double mean, double beta, double normalization) 
{
  if (sigma < 0.)     return 0.;
  double fitVal;
  double z = (x - mean)/sigma; 
  //if (alpha < 0) z = -z; 
  double abs_alpha = std::abs(alpha);
  double abs_beta = std::abs(beta);
  if (x < 0)
  {
    fitVal = 0;
  }
  else if (z  > - abs_alpha && z < abs_beta)
    fitVal = std::exp(- 0.5 * z * z);
  else if(z  <= - abs_alpha)
  {
    double nDivAlpha = n/abs_alpha;
    double AA =  std::exp(-0.5*abs_alpha*abs_alpha);
    double B = nDivAlpha -abs_alpha;
    double arg = nDivAlpha/(B-z);
    fitVal = AA * std::pow(arg,n);
  }
  else
  {
    double nDivBeta = n/abs_beta;
    double AA =  std::exp(-0.5*abs_beta*abs_beta);
    double B = nDivBeta -abs_beta;
    double arg = nDivBeta/(B+z);
    fitVal = AA * std::pow(arg,n);
  }

  return normalization * fitVal;
}

double HHKinFit2::crystalBallLikePDFROOT(double* x, double *par)
{
  double fitVal = crystBallLikePDF(x[0], par[0], par[1], par[2], par[3],
				   par[4], par[5]);
  return fitVal;
}
*/
