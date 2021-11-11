#ifdef HHKINFIT2
#include "HHKinFitMasterHeavyHiggs.h"
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
#include "HHKinFit2/HHKinFit2/interface/HHKinFitMasterHeavyHiggs.h"
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

//namespace {
//    double roundToNdigits(double x, int n = 3)
//    {
//        if(x == 0.)
//            return 0.;
//        const int p = std::min<int>(n, n - std::ceil(std::log10(std::abs(x))));
//        const double scale = std::pow(10., p);
//        return std::floor(x * scale + 0.5) / scale;
//    }
//}

HHKinFit2::HHKinFitMasterHeavyHiggs::HHKinFitMasterHeavyHiggs(TLorentzVector const& bjet1, 
                                                              TLorentzVector const& bjet2, 
                                                              TLorentzVector const& tauvis1,
                                                              TLorentzVector const& tauvis2,
                                                              TVector2 const& met, 
                                                              TMatrixD const& met_cov, 
                                                              double sigmaEbjet1,
                                                              double sigmaEbjet2,
                                                              bool istruth, 
                                                              TLorentzVector const&  heavyhiggsgen)
  :m_MET_COV(TMatrixD(4,4)), m_bjet1_COV(TMatrixD(4,4)), m_bjet2_COV(TMatrixD(4,4))
{
  m_bjet1 = HHLorentzVector(bjet1.Px(), bjet1.Py(), bjet1.Pz(), bjet1.E());
  m_bjet2 = HHLorentzVector(bjet2.Px(), bjet2.Py(), bjet2.Pz(), bjet2.E());
  m_tauvis1 = HHLorentzVector(tauvis1.Px(), tauvis1.Py(), tauvis1.Pz(), tauvis1.E());
  m_tauvis2 = HHLorentzVector(tauvis2.Px(), tauvis2.Py(), tauvis2.Pz(), tauvis2.E());

  //m_bjet1 = HHLorentzVector(roundToNdigits(bjet1.Px()), roundToNdigits(bjet1.Py()), roundToNdigits(bjet1.Pz()), roundToNdigits(bjet1.E()));
  //m_bjet2 = HHLorentzVector(roundToNdigits(bjet2.Px()), roundToNdigits(bjet2.Py()), roundToNdigits(bjet2.Pz()), roundToNdigits(bjet2.E()));
  //m_tauvis1 = HHLorentzVector(roundToNdigits(tauvis1.Px()), roundToNdigits(tauvis1.Py()), roundToNdigits(tauvis1.Pz()), roundToNdigits(tauvis1.E()));
  //m_tauvis2 = HHLorentzVector(roundToNdigits(tauvis2.Px()), roundToNdigits(tauvis2.Py()), roundToNdigits(tauvis2.Pz()), roundToNdigits(tauvis2.E()));

  m_tauvis1.SetMkeepE(1.77682);
  m_tauvis2.SetMkeepE(1.77682);

  m_bestMethodFlag = 0;
  m_loopsNeeded = 0;
  m_useAdveancedBJetChi2 = false;

  //m_MET = TVector2(roundToNdigits(met.Px()), roundToNdigits(met.Py()));
  m_MET = TVector2(met.Px(), met.Py());
  
  //m_MET_COV(0,0) = roundToNdigits(met_cov(0,0));
  //m_MET_COV(1,0) = roundToNdigits(met_cov(1,0));
  //m_MET_COV(0,1) = roundToNdigits(met_cov(0,1));
  //m_MET_COV(1,1) = roundToNdigits(met_cov(1,1));
  m_MET_COV(0,0) = met_cov(0,0);
  m_MET_COV(1,0) = met_cov(1,0);
  m_MET_COV(0,1) = met_cov(0,1);
  m_MET_COV(1,1) = met_cov(1,1);

  if(sigmaEbjet1 >= 0.0)
  {
    //m_sigma_bjet1 = roundToNdigits(sigmaEbjet1);
    m_sigma_bjet1 = sigmaEbjet1;
  }
  else
  {
    m_sigma_bjet1 = GetPFBJetRes(m_bjet1.Eta(), m_bjet1.Et());
  }
  
  if(sigmaEbjet2 >= 0.0)
  {
    //m_sigma_bjet2 = roundToNdigits(sigmaEbjet2);
    m_sigma_bjet2 = sigmaEbjet2;
  }
  else
  {
    m_sigma_bjet2 = GetPFBJetRes(m_bjet2.Eta(), m_bjet2.Et());
  }

  m_bestChi2 = 99999.0;
  m_bestHypo = HHFitHypothesisHeavyHiggs(0,0);

  if (istruth)
  {
    /*Old Systamtics Test
    Double_t bjet1_E  = r.Gaus(m_bjet1.E(), m_sigma_bjet1);
    Double_t bjet1_P  = sqrt(pow(bjet1_E,2) - pow(m_bjet1.M(),2));
    Double_t bjet1_Pt = sin(m_bjet1.Theta())*bjet1_P;
       
    double b1Px_beforeSmear = m_bjet1.Px();

    // NOT NEEDED ANYMORE. COV Matr. is calculated in HHFitObjectE
    TMatrixD bjet1Cov(4,4);
    // error propagation p=sqrt(e^2-m^2)
    Double_t bjet1_dpt = sin(m_bjet1.Theta())*m_bjet1.E()/m_bjet1.P()*m_sigma_bjet1;
    Double_t bjet1_dpx = cos(m_bjet1.Phi())*bjet1_dpt;
    Double_t bjet1_dpy = sin(m_bjet1.Phi())*bjet1_dpt;

    bjet1Cov(0,0) = pow(bjet1_dpx,2);                           
    bjet1Cov(0,1) = bjet1_dpx * bjet1_dpy;
    bjet1Cov(1,0) = bjet1_dpx * bjet1_dpy; 
    bjet1Cov(1,1) = pow(bjet1_dpy,2);
    bjet1Cov(3,3) = m_sigma_bjet1*m_sigma_bjet1;
    
    m_bjet1.SetPtEtaPhiE(bjet1_Pt, m_bjet1.Eta(), m_bjet1.Phi(), bjet1_E);
    
    //m_bjet1_COV = bjet1Cov;
    
    double b1Px_afterSmear = m_bjet1.Px();
    b1Px_standardDevs = (b1Px_beforeSmear - b1Px_afterSmear)/bjet1_dpx;

    Double_t bjet2_E  = r.Gaus(m_bjet2.E(), m_sigma_bjet2);
    Double_t bjet2_P  = sqrt(pow(bjet2_E,2) - pow(m_bjet2.M(),2));
    Double_t bjet2_Pt = sin(m_bjet2.Theta())*bjet2_P;

    
    TMatrixD bjet2Cov(4,4);

    Double_t bjet2_dpt = sin(m_bjet2.Theta())*m_bjet2.E()/m_bjet2.P()*m_sigma_bjet2;  
    Double_t bjet2_dpx = cos(m_bjet2.Phi())*bjet2_dpt;
    Double_t bjet2_dpy = sin(m_bjet2.Phi())*bjet2_dpt;

    bjet2Cov(0,0) = pow(bjet2_dpx,2); 
    bjet2Cov(0,1) = bjet2_dpx * bjet2_dpy;
    bjet2Cov(1,0) = bjet2_dpx * bjet2_dpy;
    bjet2Cov(1,1) = pow(bjet2_dpy,2);
    bjet2Cov(3,3) = m_sigma_bjet2*m_sigma_bjet2;
    
    
    m_bjet2.SetPtEtaPhiE(bjet2_Pt, m_bjet2.Eta(), m_bjet2.Phi(), bjet2_E);
    //m_bjet2_COV = bjet2Cov;
    */
    TF1* crystalBall = new TF1("crystalBall", crystalBallLikePDFROOT, 0.0, 10.0, 6);
    crystalBall->SetParameter(0,2.06421);
    crystalBall->SetParameter(1,99.92);
    crystalBall->SetParameter(2,0.17);
    crystalBall->SetParameter(3,1.15);
    crystalBall->SetParameter(4,0.466);
    crystalBall->SetParameter(5,1.608);

    /*
    TCanvas* c1 = new TCanvas;
    crystalBall.Draw();
    c1.SaveAs("crystBall.pdf");
    delete c1;
    */

    TRandom3 r(0);
    gRandom = new TRandom3(0);
    Double_t randomSmear = crystalBall->GetRandom();
    if(randomSmear < 0.1)
    {
      randomSmear = 0.1;
    }
    Double_t eTemp = m_bjet1.E() * 1/randomSmear;
    if(pow(eTemp,2) < pow(m_bjet1.M(),2))
    {
      eTemp = m_bjet1.M()*1.1;
    }

    Double_t bjet1_E  = eTemp;
    Double_t bjet1_P  = sqrt(pow(bjet1_E,2) - pow(m_bjet1.M(),2));
    Double_t bjet1_Pt = sin(m_bjet1.Theta())*bjet1_P;    
    m_bjet1.SetPtEtaPhiE(bjet1_Pt, m_bjet1.Eta(), m_bjet1.Phi(), bjet1_E);
    randomSmear = crystalBall->GetRandom();
    if(randomSmear < 0.1)
    {
      randomSmear = 0.1;
    }
    eTemp = m_bjet2.E() * 1/randomSmear;
    if(pow(eTemp,2) < pow(m_bjet2.M(),2))
    {
      eTemp = m_bjet2.M()*1.1;
    }

    Double_t bjet2_E  = eTemp;
    Double_t bjet2_P  = sqrt(pow(bjet2_E,2) - pow(m_bjet2.M(),2));
    Double_t bjet2_Pt = sin(m_bjet2.Theta())*bjet2_P;
    m_bjet2.SetPtEtaPhiE(bjet2_Pt, m_bjet2.Eta(), m_bjet2.Phi(), bjet2_E);

    HHLorentzVector heavyH;
    heavyH = HHLorentzVector(heavyhiggsgen.Px(), heavyhiggsgen.Py(),
			     heavyhiggsgen.Pz(), heavyhiggsgen.E());

    double recoilXSmeared = r.Gaus(heavyhiggsgen.Px(), 20.0);
    double recoilYSmeared = r.Gaus(heavyhiggsgen.Py(), 20.0);

    TMatrixD recoil_COV = TMatrixD(4,4);
    recoil_COV(0,0)=400;    recoil_COV(0,1)=0;
    recoil_COV(1,0)=0;      recoil_COV(1,1)=400;

    TVector2 recoilVec2 = TVector2(recoilXSmeared, recoilYSmeared);
    TVector2 b1Vec2 = TVector2(m_bjet1.Px(), m_bjet1.Py());
    TVector2 b2Vec2 = TVector2(m_bjet2.Px(), m_bjet2.Py());
    TVector2 tauVis1Vec2 = TVector2(m_tauvis1.Px(), m_tauvis1.Py());
    TVector2 tauVis2Vec2 = TVector2(m_tauvis2.Px(), m_tauvis2.Py());

    m_MET = recoilVec2 - b1Vec2 - b2Vec2 - tauVis1Vec2 - tauVis2Vec2; 
    smearedMET = m_MET;

    //m_MET_COV = recoil_COV + bjet1Cov + bjet2Cov;
    m_MET_COV = recoil_COV;
  }  
}


void HHKinFit2::HHKinFitMasterHeavyHiggs::doFit()
{
  std::cout << "DEPRECATED! Please use fit()." << std::endl;
  this->fit();
}

void HHKinFit2::HHKinFitMasterHeavyHiggs::fit()
{
  if(m_hypos.size() == 0)
  {
    addHypo(125, 125);
  }
  
  //TCanvas* c1 = new TCanvas();

  for(unsigned int i = 0; i < m_hypos.size(); ++i)
  { 
    HHFitObjectE* tau1Fit = new HHFitObjectEConstM(m_tauvis1);
    HHFitObjectE* tau2Fit = new HHFitObjectEConstM(m_tauvis2);

    double bJet1Emin;
    HHLorentzVector bJet2min = m_bjet2;
    bJet1Emin = 5.0;
    bJet2min.SetEkeepBeta(5);

    HHLorentzVector tau1min = m_tauvis1;
    tau1min.SetEkeepM(0.9*m_tauvis1.E());
      
    HHLorentzVector tau2min = m_tauvis2;
    tau2min.SetEkeepM(0.9*m_tauvis2.E());

    HHFitObjectE* b1Fit = new HHFitObjectEConstBeta(m_bjet1);
    HHFitObjectE* b2Fit = new HHFitObjectEConstBeta(m_bjet2);
  
    //prepare MET object
    HHFitObjectMET* metFit = new HHFitObjectMET(m_MET);

    //prepare composite object: Higgs
    HHFitObject* heavyHiggs  = new HHFitObjectComposite(tau1Fit, tau2Fit, 
							b1Fit, b2Fit,
							metFit);
    HHFitObject* higgs1  = new HHFitObjectComposite(tau1Fit, tau2Fit);  
    HHFitObject* higgs2  = new HHFitObjectComposite(b1Fit, b2Fit);

    int mh1 = m_hypos[i].first;
    int mh2 = m_hypos[i].second;
  
    try
    {
      tau1Fit->setFitLimitsE(tau1min,mh1,tau2min);
      tau2Fit->setFitLimitsE(tau2min,mh1,tau1min);
    }
    catch(HHLimitSettingException const& e)
    {
      std::cout << "Exception while setting tau limits:" << std::endl;
      std::cout << e.what() << std::endl;
      std::cout << "Tau energies are not compatible with invariant mass constraint." << std::endl;

      m_map_chi2[m_hypos[i]] = -pow(10,10);
      m_map_prob[m_hypos[i]] = -pow(10,10);
      m_bestHypo = HHFitHypothesisHeavyHiggs(-1,-1);
      m_bestChi2 = -pow(10,10);
      m_map_convergence[m_hypos[i]] = -1;
      throw;
      continue;
    }

    try
    {
      b1Fit->setLowerFitLimitE(bJet1Emin);
      b1Fit->setUpperFitLimitE(mh2, bJet2min);
    }
    catch(HHLimitSettingException const& e)
    {
      std::cout << "Exception while setting b-jet limits" << std::endl;
      std::cout << e.what() << std::endl;
      std::cout << "Bjet energies are not compatible with invariant mass constraint." << std::endl;
      std::cout << "B1 Momentum: " << m_bjet1.P() << std::endl;
      std::cout << "B2 Momentum: " << m_bjet2.P() << std::endl;

      m_map_chi2[m_hypos[i]] = -pow(10,10);
      m_map_prob[m_hypos[i]] = -pow(10,10);
      m_bestHypo = HHFitHypothesisHeavyHiggs(-1,-1);
      m_bestChi2 = -pow(10,10);
      m_map_convergence[m_hypos[i]] = -2;
      throw;
      continue;
    }
    
    
    if(!m_useAdveancedBJetChi2)
    {
      if(fabs(m_bjet1_COV(0,0)) < 0.001)
      {
	b1Fit->setCovMatrix(m_sigma_bjet1);
	b2Fit->setCovMatrix(m_sigma_bjet2);
	//b1Fit->printCovMatrix();
	//b2Fit->printCovMatrix();
      }
      else
      {
	std::cout << "Cov MATRIX set manually for ToyMC study." << std::endl;      
	b1Fit->setCovMatrix(m_bjet1_COV);
	b2Fit->setCovMatrix(m_bjet2_COV);
      }
    }
    

    metFit->setCovMatrix(m_MET_COV);
    heavyHiggs->setCovMatrix(m_MET_COV);// - m_bjet1_COV - m_bjet2_COV);

    //prepare constraints
    HHFitConstraint* c_invmh1 = new HHFitConstraintEHardM(tau1Fit, tau2Fit, mh1);
    HHFitConstraint* c_invmh2 = new HHFitConstraintEHardM(b1Fit, b2Fit, mh2);
    

    HHFitConstraint* c_b1;
    HHFitConstraint* c_b2;
    if(m_useAdveancedBJetChi2)
    {
       c_b1 = new HHFitConstraint4VectorBJet(b1Fit);
       c_b2 = new HHFitConstraint4VectorBJet(b2Fit);
    }
    else
    {
      c_b1 = new HHFitConstraint4Vector(b1Fit, false, false, 
					false, true);
      c_b2 = new HHFitConstraint4Vector(b2Fit, false, false, 
					false, true);
    }
//   HHFitConstraint* c_b1 = new HHFitConstraint4VectorBJet(b1Fit, 2.06421, 99.92,
//							   0.17, 1.15, 0.466, 1.608);
//    HHFitConstraint* c_b2 = new HHFitConstraint4VectorBJet(b2Fit, 2.06421, 99.92,
//							   0.17, 1.15, 0.466, 1.608);
    HHFitConstraint* c_balance = new HHFitConstraint4Vector(heavyHiggs, true, true, 
							    false, false);

    //fit
    HHKinFit2::HHKinFit* fitObject = new HHKinFit2::HHKinFit();


    tau1Fit->setInitDirection(1.0);
    tau1Fit->setInitPrecision(0.1);
    tau1Fit->setInitStepWidth(0.1*(tau1Fit->getUpperFitLimitE() - 
				   tau1Fit->getLowerFitLimitE()));

    b1Fit->setInitPrecision(0.1);
    b1Fit->setInitStepWidth(0.02*(b1Fit->getUpperFitLimitE() - 
				  b1Fit->getLowerFitLimitE()));
    b1Fit->setInitStart( 0.9*b1Fit->getInitial4Vector().E());
    b1Fit->setInitDirection(-1.0);

    fitObject->addFitObjectE(b1Fit);
    fitObject->addFitObjectE(tau1Fit);

    fitObject->addConstraint(c_invmh1);
    fitObject->addConstraint(c_invmh2);
    fitObject->addConstraint(c_b1);
    fitObject->addConstraint(c_b2);
    fitObject->addConstraint(c_balance);

    /*
    //For Chi2Map
    int steps[2];
    steps[0] = 30;
    steps[1] = 30;
    
    double mins[2];
    double maxs[2];
    //END For Chi2Map
    */

//    double b1NegTauUpFit,b1PosTauDownFit, b1NegTauDownFit;
//    double tau1NegTauUpFit, tau1PosTauDownFit, tau1NegTauDownFit;
    double chi2NegTauUp = 99999;
    double chi2PosTauDown = 99999;
    double chi2NegTauDown = 99999;
    double chi2Min = 99999;

    //--NegTauDown
    tau1Fit->setInitStart(tau1Fit->getUpperFitLimitE());
    try
    {
      fitObject->fit();
    }
    catch(HHLimitSettingException const& e)
    {
      std::cout << e.what() << std::endl;
      m_map_chi2[m_hypos[i]] = -pow(10,10);
      m_map_prob[m_hypos[i]] = -pow(10,10);
      m_bestHypo = HHFitHypothesisHeavyHiggs(-1,-1);
      m_bestChi2 = -pow(10,10);
      if((b1Fit->getUpperFitLimitE() - b1Fit->getLowerFitLimitE()) <
	 (tau1Fit->getUpperFitLimitE() - tau1Fit->getLowerFitLimitE()))
      {
	m_map_convergence[m_hypos[i]] = -2;
      }
      else
      {
	m_map_convergence[m_hypos[i]] = -1;
      }
      throw;
      continue;
    }
    catch(HHKinFit2::HHEnergyRangeException const& e){
      std::cout << "Energy Range Exception" << std::endl;
      m_map_chi2[m_hypos[i]] = -pow(10,10);
      m_map_prob[m_hypos[i]] = -pow(10,10);
      m_bestHypo = HHFitHypothesisHeavyHiggs(-1,-1);
      m_bestChi2 = -pow(10,10);
      m_map_convergence[m_hypos[i]] = 0;
      throw;
      continue;
    }

    if(fitObject->getChi2() > 0)
    {
      chi2NegTauDown = fitObject->getChi2();
      if(chi2NegTauDown < chi2Min)
      {
	chi2Min = chi2NegTauDown;
      }
      //b1NegTauDownFit = b1Fit->getFit4Vector().E();
      //tau1NegTauDownFit = tau1Fit->getFit4Vector().E();
    }


    //--NegTauUp
//      b1Fit->setInitStart( 0.9*b1Fit->getInitial4Vector().E());
    tau1Fit->setInitStart((tau1Fit->getUpperFitLimitE() - 
			   tau1Fit->getLowerFitLimitE())/2.0);
    try
    {
      fitObject->fit();
    }
    catch(HHLimitSettingException const& e)
    {
      std::cout << e.what() << std::endl;
      m_map_chi2[m_hypos[i]] = -pow(10,10);
      m_map_prob[m_hypos[i]] = -pow(10,10);
      m_bestHypo = HHFitHypothesisHeavyHiggs(-1,-1);
      m_bestChi2 = -pow(10,10);
      if((b1Fit->getUpperFitLimitE() - b1Fit->getLowerFitLimitE()) <
	 (tau1Fit->getUpperFitLimitE() - tau1Fit->getLowerFitLimitE()))
      {
	m_map_convergence[m_hypos[i]] = -2;
      }
      else
      {
	m_map_convergence[m_hypos[i]] = -1;
      }
      throw;
      continue;
    }
    catch(HHKinFit2::HHEnergyRangeException const& e){
      std::cout << "Energy Range Exception" << std::endl;
      m_map_chi2[m_hypos[i]] = -pow(10,10);
      m_map_prob[m_hypos[i]] = -pow(10,10);
      m_bestHypo = HHFitHypothesisHeavyHiggs(-1,-1);
      m_bestChi2 = -pow(10,10);
      m_map_convergence[m_hypos[i]] = 0;
      throw;
      continue;
    }

    if(fitObject->getChi2() > 0)
    {
      chi2NegTauUp = fitObject->getChi2();
      if(chi2NegTauUp < chi2Min)
      {
	chi2Min = chi2NegTauUp;
      }
      //b1NegTauUpFit = b1Fit->getFit4Vector().E();
      //tau1NegTauUpFit = tau1Fit->getFit4Vector().E();
    }

    //--PosTauDown
    tau1Fit->setInitStart(tau1Fit->getLowerFitLimitE());
    try
    {
      fitObject->fit();
    }
    catch(HHLimitSettingException const& e)
    {
      std::cout << e.what() << std::endl;
      m_map_chi2[m_hypos[i]] = -pow(10,10);
      m_map_prob[m_hypos[i]] = -pow(10,10);
      m_bestHypo = HHFitHypothesisHeavyHiggs(-1,-1);
      m_bestChi2 = -pow(10,10);
      if((b1Fit->getUpperFitLimitE() - b1Fit->getLowerFitLimitE()) <
	 (tau1Fit->getUpperFitLimitE() - tau1Fit->getLowerFitLimitE()))
      {
	m_map_convergence[m_hypos[i]] = -2;
      }
      else
      {
	m_map_convergence[m_hypos[i]] = -1;
      }
      throw;
      continue;
    }
    catch(HHKinFit2::HHEnergyRangeException const& e){
      std::cout << "Energy Range Exception" << std::endl;
      m_map_chi2[m_hypos[i]] = -pow(10,10);
      m_map_prob[m_hypos[i]] = -pow(10,10);
      m_bestHypo = HHFitHypothesisHeavyHiggs(-1,-1);
      m_bestChi2 = -pow(10,10);
      m_map_convergence[m_hypos[i]] = 0;
      throw;
      continue;
    }

    if(fitObject->getChi2() > 0)
    {
      chi2PosTauDown = fitObject->getChi2();
      if(chi2PosTauDown < chi2Min)
      {
	chi2Min = chi2PosTauDown;
      }
      //b1PosTauDownFit = b1Fit->getFit4Vector().E();
      //tau1PosTauDownFit = tau1Fit->getFit4Vector().E();
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
	tau1Fit->setInitStart(tau1Fit->getUpperFitLimitE());
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
	  tau1Fit->setInitStart((tau1Fit->getUpperFitLimitE() - 
				 tau1Fit->getLowerFitLimitE())/2.0);
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
      mins[1] = 0.9*tau1Fit->getFit4Vector().E();
      maxs[0] = 1.1*b1Fit->getFit4Vector().E();
      maxs[1] = 1.1*tau1Fit->getFit4Vector().E();
 
      TString b1Point;
      TString tau1Point;
      std::stringstream stream;

      stream << b1Fit->getFit4Vector().E();
      stream >> b1Point;
      stream.clear();

      stream << tau1Fit->getFit4Vector().E();
      stream >> tau1Point;
      stream.clear();

      TGraph2D* graph = (TGraph2D*)fitObject->getChi2Function(steps, mins, maxs);
      graph->Draw("Cont3COLZ");

      c1->SaveAs("NoConvergence" + b1Point + "-" + tau1Point + ".pdf");
    }
    */

    {//Filling Results
      initialHH = (TLorentzVector)heavyHiggs->getInitial4Vector();
      finalHH = (TLorentzVector)heavyHiggs->getFit4Vector();

      m_map_convergence[m_hypos[i]] = fitObject->getConvergence();    
    
      double chi2 = fitObject->getChi2();
      m_map_chi2[m_hypos[i]] = chi2;
      m_map_prob[m_hypos[i]] = TMath::Prob(chi2, 2);
      m_loopsNeeded = fitObject->m_loopsNeeded;

      if(chi2 < m_bestChi2)
      {
	m_bestHypo = m_hypos[i];
	m_bestChi2 = chi2;
      }
    
      m_map_mH[m_hypos[i]] = heavyHiggs->getFit4Vector().M();
      m_map_chi2BJet1[m_hypos[i]] = c_b1->getChi2();
      m_map_chi2BJet2[m_hypos[i]] = c_b2->getChi2();
      m_map_chi2Balance[m_hypos[i]] = c_balance->getChi2();
    
      TLorentzVector fittedTau1 =  ( (TLorentzVector)tau1Fit->getFit4Vector()  );
      m_map_fittedTau1[m_hypos[i]] = fittedTau1;
      TLorentzVector fittedTau2 =  ( (TLorentzVector)tau2Fit->getFit4Vector()  );
      m_map_fittedTau2[m_hypos[i]] = fittedTau2;
      TLorentzVector fittedB1 =  ( (TLorentzVector)b1Fit->getFit4Vector() ) ;
      m_map_fittedB1[m_hypos[i]] = fittedB1;
      TLorentzVector fittedB2 =  ( (TLorentzVector)b2Fit->getFit4Vector() ) ;
      m_map_fittedB2[m_hypos[i]]= fittedB2;
    }

    delete c_invmh1;
    delete c_invmh2;
    delete c_b1;
    delete c_b2;
    delete c_balance;
    
    delete fitObject;
   
    delete tau1Fit;
    delete tau2Fit;
    delete b1Fit;
    delete b2Fit;
    delete metFit;
    delete heavyHiggs;
    delete higgs1;
    delete higgs2;
  }
}

void HHKinFit2::HHKinFitMasterHeavyHiggs::addHypo(HHFitHypothesisHeavyHiggs hypo)
{
  m_hypos.push_back(hypo);
}

void HHKinFit2::HHKinFitMasterHeavyHiggs::addHypo(int mh1, int mh2)
{
  HHFitHypothesisHeavyHiggs hypo(mh1, mh2);
  m_hypos.push_back(hypo);
}

void HHKinFit2::HHKinFitMasterHeavyHiggs::useAdvancedBJetChi2(bool useAdvancedBJetChi2)
{
  m_useAdveancedBJetChi2 = useAdvancedBJetChi2;
}

//Getters
HHKinFit2::HHFitHypothesisHeavyHiggs HHKinFit2::HHKinFitMasterHeavyHiggs::getBestHypothesis(){
  return(m_bestHypo);
}

HHKinFit2::HHKinFitMasterHeavyHiggs::HHKinFitMasterHeavyHiggs(const TLorentzVector* bjet1, 
                                                              const TLorentzVector* bjet2, 
                                                              const TLorentzVector* tauvis1,
                                                              const TLorentzVector* tauvis2,
                                                              const TLorentzVector* met, 
                                                              TMatrixD met_cov, 
                                                              double sigmaEbjet1,
                                                              double sigmaEbjet2,
                                                              bool istruth, 
                                                              TLorentzVector* heavyhiggsgen)
  :m_MET_COV(TMatrixD(4,4)), m_bjet1_COV(TMatrixD(4,4)), m_bjet2_COV(TMatrixD(4,4))
					  
{
  std::cout << "DEPRECATED! Please use the new Constructor taking references instead of pointers." << std::endl;
  m_bjet1 = HHLorentzVector(bjet1->Px(), bjet1->Py(), 
			    bjet1->Pz(), bjet1->E());
  m_bjet2 = HHLorentzVector(bjet2->Px(), bjet2->Py(), 
			    bjet2->Pz(), bjet2->E());
  m_tauvis1 = HHLorentzVector(tauvis1->Px(), tauvis1->Py(), 
			      tauvis1->Pz(), tauvis1->E());  
  m_tauvis2 = HHLorentzVector(tauvis2->Px(), tauvis2->Py(), 
			      tauvis2->Pz(), tauvis2->E());
 
  m_tauvis1.SetMkeepE(1.77682);
  m_tauvis2.SetMkeepE(1.77682);

  m_bestMethodFlag = 0;
  m_loopsNeeded = 0;
  m_useAdveancedBJetChi2 = false;

  if(met != 0)
  {
    m_MET = TVector2(met->Px(), met->Py());
  }
  m_MET_COV(0,0) = met_cov(0,0);
  m_MET_COV(1,0) = met_cov(1,0);
  m_MET_COV(0,1) = met_cov(0,1);
  m_MET_COV(1,1) = met_cov(1,1);

  if(sigmaEbjet1 >= 0.0)
  {
    m_sigma_bjet1 = sigmaEbjet1;
  }
  else
  {
    m_sigma_bjet1 = GetPFBJetRes(m_bjet1.Eta(), m_bjet1.Et());
  }
  
  if(sigmaEbjet2 >= 0.0)
  {
    m_sigma_bjet2 = sigmaEbjet2;
  }
  else
  {
    m_sigma_bjet2 = GetPFBJetRes(m_bjet2.Eta(), m_bjet2.Et());
  }

  m_bestChi2 = 99999.0;
  m_bestHypo = HHFitHypothesisHeavyHiggs(0,0);

  if (istruth)
  {
    /*Old Systamtics Test
    Double_t bjet1_E  = r.Gaus(m_bjet1.E(), m_sigma_bjet1);
    Double_t bjet1_P  = sqrt(pow(bjet1_E,2) - pow(m_bjet1.M(),2));
    Double_t bjet1_Pt = sin(m_bjet1.Theta())*bjet1_P;
       
    double b1Px_beforeSmear = m_bjet1.Px();

    // NOT NEEDED ANYMORE. COV Matr. is calculated in HHFitObjectE
    TMatrixD bjet1Cov(4,4);
    // error propagation p=sqrt(e^2-m^2)
    Double_t bjet1_dpt = sin(m_bjet1.Theta())*m_bjet1.E()/m_bjet1.P()*m_sigma_bjet1;
    Double_t bjet1_dpx = cos(m_bjet1.Phi())*bjet1_dpt;
    Double_t bjet1_dpy = sin(m_bjet1.Phi())*bjet1_dpt;

    bjet1Cov(0,0) = pow(bjet1_dpx,2);                           
    bjet1Cov(0,1) = bjet1_dpx * bjet1_dpy;
    bjet1Cov(1,0) = bjet1_dpx * bjet1_dpy; 
    bjet1Cov(1,1) = pow(bjet1_dpy,2);
    bjet1Cov(3,3) = m_sigma_bjet1*m_sigma_bjet1;
    
    m_bjet1.SetPtEtaPhiE(bjet1_Pt, m_bjet1.Eta(), m_bjet1.Phi(), bjet1_E);
    
    //m_bjet1_COV = bjet1Cov;
    
    double b1Px_afterSmear = m_bjet1.Px();
    b1Px_standardDevs = (b1Px_beforeSmear - b1Px_afterSmear)/bjet1_dpx;

    Double_t bjet2_E  = r.Gaus(m_bjet2.E(), m_sigma_bjet2);
    Double_t bjet2_P  = sqrt(pow(bjet2_E,2) - pow(m_bjet2.M(),2));
    Double_t bjet2_Pt = sin(m_bjet2.Theta())*bjet2_P;

    
    TMatrixD bjet2Cov(4,4);

    Double_t bjet2_dpt = sin(m_bjet2.Theta())*m_bjet2.E()/m_bjet2.P()*m_sigma_bjet2;  
    Double_t bjet2_dpx = cos(m_bjet2.Phi())*bjet2_dpt;
    Double_t bjet2_dpy = sin(m_bjet2.Phi())*bjet2_dpt;

    bjet2Cov(0,0) = pow(bjet2_dpx,2); 
    bjet2Cov(0,1) = bjet2_dpx * bjet2_dpy;
    bjet2Cov(1,0) = bjet2_dpx * bjet2_dpy;
    bjet2Cov(1,1) = pow(bjet2_dpy,2);
    bjet2Cov(3,3) = m_sigma_bjet2*m_sigma_bjet2;
    
    
    m_bjet2.SetPtEtaPhiE(bjet2_Pt, m_bjet2.Eta(), m_bjet2.Phi(), bjet2_E);
    //m_bjet2_COV = bjet2Cov;
    */
    TF1* crystalBall = new TF1("crystalBall", crystalBallLikePDFROOT, 0.0, 10.0, 6);
    crystalBall->SetParameter(0,2.06421);
    crystalBall->SetParameter(1,99.92);
    crystalBall->SetParameter(2,0.17);
    crystalBall->SetParameter(3,1.15);
    crystalBall->SetParameter(4,0.466);
    crystalBall->SetParameter(5,1.608);

    /*
    TCanvas* c1 = new TCanvas;
    crystalBall->Draw();
    c1->SaveAs("crystBall.pdf");
    delete c1;
    */

    TRandom3 r(0);
    gRandom = new TRandom3(0);
    Double_t randomSmear = crystalBall->GetRandom();
    if(randomSmear < 0.1)
    {
      randomSmear = 0.1;
    }
    Double_t eTemp = m_bjet1.E() * 1/randomSmear;
    if(pow(eTemp,2) < pow(m_bjet1.M(),2))
    {
      eTemp = m_bjet1.M()*1.1;
    }

    Double_t bjet1_E  = eTemp;
    Double_t bjet1_P  = sqrt(pow(bjet1_E,2) - pow(m_bjet1.M(),2));
    Double_t bjet1_Pt = sin(m_bjet1.Theta())*bjet1_P;    
    m_bjet1.SetPtEtaPhiE(bjet1_Pt, m_bjet1.Eta(), m_bjet1.Phi(), bjet1_E);
    randomSmear = crystalBall->GetRandom();
    if(randomSmear < 0.1)
    {
      randomSmear = 0.1;
    }
    eTemp = m_bjet2.E() * 1/randomSmear;
    if(pow(eTemp,2) < pow(m_bjet2.M(),2))
    {
      eTemp = m_bjet2.M()*1.1;
    }

    Double_t bjet2_E  = eTemp;
    Double_t bjet2_P  = sqrt(pow(bjet2_E,2) - pow(m_bjet2.M(),2));
    Double_t bjet2_Pt = sin(m_bjet2.Theta())*bjet2_P;
    m_bjet2.SetPtEtaPhiE(bjet2_Pt, m_bjet2.Eta(), m_bjet2.Phi(), bjet2_E);

    HHLorentzVector heavyH;
    heavyH = HHLorentzVector(heavyhiggsgen->Px(), heavyhiggsgen->Py(),
			     heavyhiggsgen->Pz(), heavyhiggsgen->E());

    double recoilXSmeared = r.Gaus(heavyhiggsgen->Px(), 20.0);
    double recoilYSmeared = r.Gaus(heavyhiggsgen->Py(), 20.0);

    TMatrixD recoil_COV = TMatrixD(4,4);
    recoil_COV(0,0)=400;    recoil_COV(0,1)=0;
    recoil_COV(1,0)=0;      recoil_COV(1,1)=400;

    TVector2 recoilVec2 = TVector2(recoilXSmeared, recoilYSmeared);
    TVector2 b1Vec2 = TVector2(m_bjet1.Px(), m_bjet1.Py());
    TVector2 b2Vec2 = TVector2(m_bjet2.Px(), m_bjet2.Py());
    TVector2 tauVis1Vec2 = TVector2(m_tauvis1.Px(), m_tauvis1.Py());
    TVector2 tauVis2Vec2 = TVector2(m_tauvis2.Px(), m_tauvis2.Py());

    m_MET = recoilVec2 - b1Vec2 - b2Vec2 - tauVis1Vec2 - tauVis2Vec2; 
    smearedMET = m_MET;

    //m_MET_COV = recoil_COV + bjet1Cov + bjet2Cov;
    m_MET_COV = recoil_COV;
  }  
}

HHKinFit2::HHFitHypothesisHeavyHiggs HHKinFit2::HHKinFitMasterHeavyHiggs::getLowestChi2Hypothesis(){
  std::cout << "DEPRECATED! Please use getBestHypothesis()." << std::endl;
  return(getBestHypothesis());
}


double HHKinFit2::HHKinFitMasterHeavyHiggs::getBestChi2(){
  return(m_bestChi2);
}

//Getters for fit results

double HHKinFit2::HHKinFitMasterHeavyHiggs::getChi2(int mh1, int mh2){
  HHFitHypothesisHeavyHiggs hypo(mh1, mh2);
  return(m_map_chi2[hypo]);
}

  
double HHKinFit2::HHKinFitMasterHeavyHiggs::getChi2BJet1(int mh1, int mh2){
  HHFitHypothesisHeavyHiggs hypo(mh1, mh2);
  return(m_map_chi2BJet1[hypo]);
}



double HHKinFit2::HHKinFitMasterHeavyHiggs::getChi2BJet2(int mh1, int mh2){
  HHFitHypothesisHeavyHiggs hypo(mh1, mh2);
  return(m_map_chi2BJet2[hypo]);
}
    

double HHKinFit2::HHKinFitMasterHeavyHiggs::getChi2Balance(int mh1, int mh2){
  HHFitHypothesisHeavyHiggs hypo(mh1, mh2);
  return(m_map_chi2Balance[hypo]);
}
  
double HHKinFit2::HHKinFitMasterHeavyHiggs::getFitProb(int mh1, int mh2){
  HHFitHypothesisHeavyHiggs hypo(mh1, mh2);
  return( m_map_prob[hypo]);
}
  
double HHKinFit2::HHKinFitMasterHeavyHiggs::getMH(int mh1, int mh2){
  HHFitHypothesisHeavyHiggs hypo(mh1, mh2);
  return(m_map_mH[hypo]);
}


int HHKinFit2::HHKinFitMasterHeavyHiggs::getConvergence(int mh1, int mh2){
  HHFitHypothesisHeavyHiggs hypo(mh1, mh2);
  return(m_map_convergence[hypo]);
}

TLorentzVector HHKinFit2::HHKinFitMasterHeavyHiggs::getFittedTau1(int mh1, int mh2){
  HHFitHypothesisHeavyHiggs hypo(mh1, mh2);
  return(m_map_fittedTau1[hypo]);
}


TLorentzVector HHKinFit2::HHKinFitMasterHeavyHiggs::getFittedTau2(int mh1, int mh2){
  HHFitHypothesisHeavyHiggs hypo(mh1, mh2);
  return(m_map_fittedTau2[hypo]);
}

TLorentzVector HHKinFit2::HHKinFitMasterHeavyHiggs::getFittedBJet1(int mh1, int mh2){
  HHFitHypothesisHeavyHiggs hypo(mh1, mh2);
  return(m_map_fittedB1[hypo]);
}
  

TLorentzVector HHKinFit2::HHKinFitMasterHeavyHiggs::getFittedBJet2(int mh1, int mh2){
  HHFitHypothesisHeavyHiggs hypo(mh1, mh2);
  return(m_map_fittedB2[hypo]);
}

void HHKinFit2::HHKinFitMasterHeavyHiggs::setAdvancedBalance(const TLorentzVector* met, 
                                                             TMatrixD met_cov)
{
  std::cout << "DEPRECATED! Please hand over MET and the covariance matrix already in the constructor." << std::endl;
  if(met != 0)
  {
    m_MET = TVector2(met->Px(), met->Py());
  }
  m_MET_COV(0,0) = met_cov(0,0);
  m_MET_COV(1,0) = met_cov(1,0);
  m_MET_COV(0,1) = met_cov(0,1);
  m_MET_COV(1,1) = met_cov(1,1);
}

double HHKinFit2::HHKinFitMasterHeavyHiggs::GetPFBJetRes(double eta, double et){
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

double HHKinFit2::crystBallLikePDF(double x, double alpha, double n, double sigma, 
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
