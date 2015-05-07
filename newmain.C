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

  TF1 *PDF1 = new TF1("PDF1","2*x",0,1);
  PDF1->SetNpx(100000);
  TF1 *PDF2= new TF1("PDF2","2-2*x",0,1);
  PDF2->SetNpx(100000);
  TMatrixD covarmatrix(2,2);
  covarmatrix[0][0]=130;
  covarmatrix[0][1]=0;
  covarmatrix[1][0]=0;
  covarmatrix[1][1]=130;
  HHTauTauEventGenerator Higgsgenerator(PDF1,PDF2,covarmatrix);
  //HHTauTauEventGenerator Higgsgenerator(PDF1,PDF1,covarmatrix);
  //Higgsgenerator.setMhiggs(91.1876);
  unsigned int brake=0;

  TH1D h_hfFinalChi2("h_hfFinalChi2","The Final chi2 from the higgs-KinFit",100,-6,6);
  TH1D h_hflikelihoodFinalChi2("h_hflikelihoodFinalChi2","The Final chi2 from the higgslikelihood-KinFit",100,-6,6);
  TH1D h_zfFinalChi2("h_zfFinalChi2","The Final chi2 from the Z-KinFit",100,-6,6);
  TH1D h_zflikelihoodFinalChi2("h_zflikelihoodFinalChi2","The Final chi2 from the Zlikelihood-KinFit",100,-6,6);

  TH2D h_chi2comparehfzf("h_chi2comparehfzf","comparison of hfchi2 and zfchi2 ",100,-6,6,100,-6,6);
  TH2D h_chi2comparehflikelihoodzflikelihood("h_chi2comparehflikelihoodzflikelihood","comparison of hflikelihoodchi2 and zflikelihoodchi2 ",100,-6,6,100,-6,6);
  TH2D h_chi2comparehfhflikelihood("h_chi2comparehfhflikelihood","comparison of hfchi2 and hflikelihoodchi2 ",100,-6,6,100,-6,6);
  TH2D h_chi2comparehfzflikelihood("h_chi2comparehfzflikelihood","comparison of zfchi2 and zflikelihoodchi2 ",100,-6,6,100,-6,6);

  TH1D h_Lhf("h_Lhf","hf likelihood",100,0,0.005);
  TH1D h_Lhflikelihood("h_Lhflikelihood","hflikelihood likelihood",100,0,0.005);
  TH1D h_Lzf("h_Lzf","zf likelihood",100,0,0.005);
  TH1D h_Lzflikelihood("h_Lzflikelihood","zflikelihood likelihood",100,0,0.005);

  TH2D h_hfLhflikelihoodL("h_hfLhflikelihoodL","compare hfL and hflikelihoodL",100,0,0.005,100,0,0.005);
  TH2D h_zfLzflikelihoodL("h_zfLzflikelihoodL","compare zfL and zflikelihoodL",100,0,0.005,100,0,0.005);

  TH2D h_hfLkdkincompareLkddyn("h_hfLkdkincompareLkddyn","compare the kinematic part of L with the dynamic part of L",1000,0,0.00125,100,0,4);
  TH2D h_zfLkdkincompareLkddyn("h_zfLkdkincompareLkddyn","compare the kinematic part of L with the dynamic part of L",1000,0,0.00125,100,0,4);

  TH2D h_EkcompareEkdtau1h("h_EkcompareEkdtau1h","h_EkcompareEkdtau1h",1000,0,1000,1000,0,1000);

  TH1D h_Eresulutionkhtau1("h_Eresulutionkhtau1","h_Eresulutionkhtau1",1000,-1.5,1.5);
  TH1D h_Eresulutionkdhtau1("h_Eresulutionkdhtau1","h_Eresulutionkdhtau1",1000,-1.5,1.5);

  TH2D h_EkcompareEkdtau2h("h_EkcompareEkdtau2h","h_EkcompareEkdtau2h",1000,0,1000,1000,0,1000);

  TH1D h_Eresulutionkhtau2("h_Eresulutionkhtau2","h_Eresulutionkhtau2",1000,-1.5,1.5);
  TH1D h_Eresulutionkdhtau2("h_Eresulutionkdhtau2","h_Eresulutionkdhtau2",1000,-1.5,1.5);

  TH1D h_EkdminusEknormtau1h("h_EkdminusEknormtau1h","h_EkdminusEknormtau1h",1000,-1.5,1.5);
  TH1D h_EkdminusEknormtau2h("h_EkdminusEknormtau2h","h_EkdminusEknormtau2h",1000,-1.5,1.5);

  for(unsigned int i=0; i<1000000; i++){

	  //if (1000000%i==0) std::cout << 100*(i/1000000) << "%" << std::endl;
    HHKinFit* singlefithf = new HHKinFit();
    HHKinFit* singlefithflikelihood = new HHKinFit();
    HHKinFit* singlefitzf = new HHKinFit();
    HHKinFit* singlefitzflikelihood = new HHKinFit();
    HHFitConstraint* balancehflikelihood;
    HHFitConstraint* Likelihoodhf;
    HHFitConstraint* balancezflikelihood;
    HHFitConstraint* Likelihoodzf;
    HHFitObjectE* tau1hf;
    HHFitObjectE* tau1hflikelihood;
    HHFitObjectE* tau2hf;
    HHFitObjectE* tau2hflikelihood;

    try {
    try{
	Higgsgenerator.generateEvent();
      }
      catch(const HHEnergyRangeException& e){
	i--;
	continue;
      }
      //Higgs-fit:

      double masshf=125.7;

      tau1hf = new HHFitObjectEConstM(Higgsgenerator.getTau1Vis());
      tau2hf = new HHFitObjectEConstM(Higgsgenerator.getTau2Vis());

      tau1hf->setLowerFitLimitE(tau1hf->getInitial4Vector());
      tau1hf->setUpperFitLimitE(masshf,tau2hf->getInitial4Vector());
      tau2hf->setLowerFitLimitE(tau2hf->getInitial4Vector());
      tau2hf->setUpperFitLimitE(masshf,tau1hf->getInitial4Vector());


      HHFitObjectMET* methf = new HHFitObjectMET(TVector2(Higgsgenerator.getMETwithsigma()[0],Higgsgenerator.getMETwithsigma()[1]));
      methf->setCovMatrix(Higgsgenerator.getCovarmatrix()[0][0],Higgsgenerator.getCovarmatrix()[1][1],Higgsgenerator.getCovarmatrix()[1][0]);

      HHFitObject* higgshf  = new HHFitObjectComposite(tau1hf, tau2hf, methf);

      HHFitConstraint* invmhf = new HHFitConstraintEHardM(tau1hf, tau2hf, masshf);
      HHFitConstraint* balancehf = new HHFitConstraint4Vector(higgshf, true, true, false, false);


      singlefithf->addFitObjectE(tau1hf);
      singlefithf->addConstraint(invmhf);
      singlefithf->addConstraint(balancehf);

      try {
	singlefithf->fit();
	if (!((singlefithf->getConvergence()==1)||(singlefithf->getConvergence()==2))) {
	  i--;
	  continue;
	}
      }
      catch(HHEnergyRangeException const& e){
	//std::cout << i << std::endl;
	i--;
	continue;
      }


      //Higgs likelihood fit
      tau1hflikelihood = new HHFitObjectEConstM(Higgsgenerator.getTau1Vis());
      tau2hflikelihood = new HHFitObjectEConstM(Higgsgenerator.getTau2Vis());

      tau1hflikelihood->setLowerFitLimitE(tau1hflikelihood->getInitial4Vector());
      tau1hflikelihood->setUpperFitLimitE(masshf,tau2hflikelihood->getInitial4Vector());
      tau2hflikelihood->setLowerFitLimitE(tau2hflikelihood->getInitial4Vector());
      tau2hflikelihood->setUpperFitLimitE(masshf,tau1hflikelihood->getInitial4Vector());

      HHFitObjectMET* methflikelihood = new HHFitObjectMET(TVector2(Higgsgenerator.getMETwithsigma()[0],Higgsgenerator.getMETwithsigma()[1]));
      methflikelihood->setCovMatrix(Higgsgenerator.getCovarmatrix()[0][0],Higgsgenerator.getCovarmatrix()[1][1],Higgsgenerator.getCovarmatrix()[1][0]);

      HHFitObject* higgshflikelihood  = new HHFitObjectComposite(tau1hflikelihood, tau2hflikelihood, methflikelihood);

      HHFitConstraint* invmhflikelihood = new HHFitConstraintEHardM(tau1hflikelihood, tau2hflikelihood, masshf);
      balancehflikelihood = new HHFitConstraint4Vector(higgshflikelihood, true, true, false, false);
      Likelihoodhf = new HHFitConstraintLikelihood(tau1hflikelihood,tau2hflikelihood,PDF1,PDF2);


      singlefithflikelihood->addFitObjectE(tau1hflikelihood);
      singlefithflikelihood->addConstraint(invmhflikelihood);
      singlefithflikelihood->addConstraint(balancehflikelihood);
      singlefithflikelihood->addConstraint(Likelihoodhf);

      try {
	singlefithflikelihood->fit();
	if (!((singlefithflikelihood->getConvergence()==1)||(singlefithflikelihood->getConvergence()==2))){
	  i--;
	  continue;
	}
      }
      catch(HHEnergyRangeException const& e){
	//std::cout << i << std::endl;
	i--;
	continue;
      }

      //Z-Fit
      double masszf=91.1876;

      HHFitObjectE* tau1zf = new HHFitObjectEConstM(Higgsgenerator.getTau1Vis());
      HHFitObjectE* tau2zf = new HHFitObjectEConstM(Higgsgenerator.getTau2Vis());


      tau1zf->setLowerFitLimitE(tau1zf->getInitial4Vector());
      tau1zf->setUpperFitLimitE(masszf,tau2zf->getInitial4Vector());
      tau2zf->setLowerFitLimitE(tau2zf->getInitial4Vector());
      tau2zf->setUpperFitLimitE(masszf,tau1zf->getInitial4Vector());


      HHFitObjectMET* metzf = new HHFitObjectMET(TVector2(Higgsgenerator.getMETwithsigma()[0],Higgsgenerator.getMETwithsigma()[1]));
      metzf->setCovMatrix(Higgsgenerator.getCovarmatrix()[0][0],Higgsgenerator.getCovarmatrix()[1][1],Higgsgenerator.getCovarmatrix()[1][0]);

      HHFitObject* higgszf  = new HHFitObjectComposite(tau1zf, tau2zf, metzf);

      HHFitConstraint* invmzf = new HHFitConstraintEHardM(tau1zf, tau2zf, masszf);
      HHFitConstraint* balancezf = new HHFitConstraint4Vector(higgszf, true, true, false, false);


      singlefitzf->addFitObjectE(tau1zf);
      singlefitzf->addConstraint(invmzf);
      singlefitzf->addConstraint(balancezf);

      try {
	singlefitzf->fit();
	if (!((singlefitzf->getConvergence()==1)||(singlefitzf->getConvergence()==2))){
	  i--;
	  continue;
	}
      }
      catch(HHEnergyRangeException const& e){
	//std::cout << i << std::endl;
	i--;
	continue;
      }

      //Z likelihood fit
      HHFitObjectE* tau1zflikelihood = new HHFitObjectEConstM(Higgsgenerator.getTau1Vis());
      HHFitObjectE* tau2zflikelihood = new HHFitObjectEConstM(Higgsgenerator.getTau2Vis());

      tau1zflikelihood->setLowerFitLimitE(tau1zflikelihood->getInitial4Vector());
      tau1zflikelihood->setUpperFitLimitE(masszf,tau2zflikelihood->getInitial4Vector());
      tau2zflikelihood->setLowerFitLimitE(tau2zflikelihood->getInitial4Vector());
      tau2zflikelihood->setUpperFitLimitE(masszf,tau1zflikelihood->getInitial4Vector());

      HHFitObjectMET* metzflikelihood = new HHFitObjectMET(TVector2(Higgsgenerator.getMETwithsigma()[0],Higgsgenerator.getMETwithsigma()[1]));
      metzflikelihood->setCovMatrix(Higgsgenerator.getCovarmatrix()[0][0],Higgsgenerator.getCovarmatrix()[1][1],Higgsgenerator.getCovarmatrix()[1][0]);

      HHFitObject* higgszflikelihood  = new HHFitObjectComposite(tau1zflikelihood, tau2zflikelihood, metzflikelihood);

      HHFitConstraint* invmzflikelihood = new HHFitConstraintEHardM(tau1zflikelihood, tau2zflikelihood, masszf);
      balancezflikelihood = new HHFitConstraint4Vector(higgszflikelihood, true, true, false, false);
      Likelihoodzf = new HHFitConstraintLikelihood(tau1zflikelihood,tau2zflikelihood,PDF1,PDF1);


      singlefitzflikelihood->addFitObjectE(tau1zflikelihood);
      singlefitzflikelihood->addConstraint(invmzflikelihood);
      singlefitzflikelihood->addConstraint(balancezflikelihood);
      singlefitzflikelihood->addConstraint(Likelihoodzf);

      try {
	singlefitzflikelihood->fit();
	if (!((singlefitzflikelihood->getConvergence()==1)||(singlefitzflikelihood->getConvergence()==2))){
	  i--;
	  continue;
	}
      }
      catch(HHEnergyRangeException const& e){
	//std::cout << i << std::endl;
	i--;
	continue;
      }


    }
    catch(HHEnergyRangeException const& e){
      //std::cout << i << std::endl;
      i--;
      brake++;
      continue;

    }

    h_hfFinalChi2.Fill(singlefithf->getChi2());
    h_hflikelihoodFinalChi2.Fill(singlefithflikelihood->getChi2());
    h_zfFinalChi2.Fill(singlefitzf->getChi2());
    h_zflikelihoodFinalChi2.Fill(singlefitzflikelihood->getChi2());

    h_chi2comparehfzf.Fill(singlefithf->getChi2(),singlefitzf->getChi2());
    h_chi2comparehflikelihoodzflikelihood.Fill(singlefithflikelihood->getChi2(),singlefitzflikelihood->getChi2());
    h_chi2comparehfhflikelihood.Fill(singlefithf->getChi2(),singlefithflikelihood->getChi2());
    h_chi2comparehfzflikelihood.Fill(singlefitzf->getChi2(),singlefitzflikelihood->getChi2());

    h_Lhf.Fill(singlefithf->getL());
    h_Lhflikelihood.Fill(singlefithflikelihood->getL());
    h_Lzf.Fill(singlefitzf->getL());
    h_Lzflikelihood.Fill(singlefitzflikelihood->getL());

    h_hfLhflikelihoodL.Fill(singlefithf->getL(),singlefithflikelihood->getL());
    h_zfLzflikelihoodL.Fill(singlefitzf->getL(),singlefitzflikelihood->getL());



    h_hfLkdkincompareLkddyn.Fill(balancehflikelihood->getLikelihood(),Likelihoodhf->getLikelihood());
    h_zfLkdkincompareLkddyn.Fill(balancezflikelihood->getLikelihood(),Likelihoodzf->getLikelihood());

    h_EkcompareEkdtau1h.Fill(tau1hf->getFit4Vector().E(),tau1hflikelihood->getFit4Vector().E());

    h_Eresulutionkhtau1.Fill((Higgsgenerator.getTau1boosted().E()-tau1hf->getFit4Vector().E())/Higgsgenerator.getTau1boosted().E());
    h_Eresulutionkdhtau1.Fill((Higgsgenerator.getTau1boosted().E()-tau1hflikelihood->getFit4Vector().E())/Higgsgenerator.getTau1boosted().E());

    h_EkcompareEkdtau2h.Fill(tau2hf->getFit4Vector().E(),tau2hflikelihood->getFit4Vector().E());

    h_Eresulutionkhtau2.Fill((Higgsgenerator.getTau2boosted().E()-tau2hf->getFit4Vector().E())/Higgsgenerator.getTau2boosted().E());
    h_Eresulutionkdhtau2.Fill((Higgsgenerator.getTau2boosted().E()-tau2hflikelihood->getFit4Vector().E())/Higgsgenerator.getTau2boosted().E());

    h_EkdminusEknormtau1h.Fill((tau1hflikelihood->getFit4Vector().E()-tau1hf->getFit4Vector().E())/Higgsgenerator.getTau1boosted().E());
    h_EkdminusEknormtau2h.Fill((tau2hflikelihood->getFit4Vector().E()-tau2hf->getFit4Vector().E())/Higgsgenerator.getTau2boosted().E());



  }

  std::cout << brake << std::endl;
  TFile controlplots("multiplehiggssamplefit.root","RECREATE");

  h_hfFinalChi2.Write();
  h_hflikelihoodFinalChi2.Write();
  h_zfFinalChi2.Write();
  h_zflikelihoodFinalChi2.Write();

  h_chi2comparehfzf.Write();
  h_chi2comparehflikelihoodzflikelihood.Write();
  h_chi2comparehfhflikelihood.Write();
  h_chi2comparehfzflikelihood.Write();

  h_Lhf.Write();
  h_Lhflikelihood.Write();
  h_Lzf.Write();
  h_Lzflikelihood.Write();

  h_hfLhflikelihoodL.Write();
  h_zfLzflikelihoodL.Write();

  h_hfLkdkincompareLkddyn.Write();
  h_zfLkdkincompareLkddyn.Write();

  h_EkcompareEkdtau1h.Write();

  h_Eresulutionkhtau1.Write();
  h_Eresulutionkdhtau1.Write();

  h_EkcompareEkdtau2h.Write();

  h_Eresulutionkhtau2.Write();
  h_Eresulutionkdhtau2.Write();

  h_EkdminusEknormtau1h.Write();
  h_EkdminusEknormtau2h.Write();



  //---------------------------------------------------------------------------------------------
  // L(E)

  for(unsigned int i=0; i<3; i++){
    try{
      Higgsgenerator.generateEvent();

    }
    catch(const HHEnergyRangeException& e){
      i--;
      continue;
    }
    double massh=125.7;

    HHFitObjectE*  tau1k = new HHFitObjectEConstM(Higgsgenerator.getTau1Vis());
    HHFitObjectE* tau2k = new HHFitObjectEConstM(Higgsgenerator.getTau2Vis());

    tau1k->setLowerFitLimitE(tau1k->getInitial4Vector());
    tau1k->setUpperFitLimitE(massh,tau2k->getInitial4Vector());
    tau2k->setLowerFitLimitE(tau2k->getInitial4Vector());
    tau2k->setUpperFitLimitE(massh,tau1k->getInitial4Vector());

    HHFitObjectMET* metk = new HHFitObjectMET(TVector2(Higgsgenerator.getMETwithsigma()[0],Higgsgenerator.getMETwithsigma()[1]));
    metk->setCovMatrix(Higgsgenerator.getCovarmatrix()[0][0],Higgsgenerator.getCovarmatrix()[1][1],Higgsgenerator.getCovarmatrix()[1][0]);

    HHFitObject* higgsk  = new HHFitObjectComposite(tau1k, tau2k, metk);

    HHFitConstraint* invmk = new HHFitConstraintEHardM(tau1k, tau2k, massh);
    HHFitConstraint* balancek = new HHFitConstraint4Vector(higgsk, true, true, false, false);


    HHKinFit* singlefitk = new HHKinFit();

    singlefitk->addFitObjectE(tau1k);
    singlefitk->addConstraint(invmk);
    singlefitk->addConstraint(balancek);


    try {
      singlefitk->fit();
      if (!((singlefitk->getConvergence()==1)||(singlefitk->getConvergence()==2))){
	i--;
	continue;
      }		}
    catch(HHEnergyRangeException const& e){
      std::cout << i << std::endl;
      i--;
      continue;
    }
    HHFitObjectE*  tau1kd = new HHFitObjectEConstM(Higgsgenerator.getTau1Vis());
    HHFitObjectE* tau2kd = new HHFitObjectEConstM(Higgsgenerator.getTau2Vis());

    tau1kd->setLowerFitLimitE(tau1kd->getInitial4Vector());
    tau1kd->setUpperFitLimitE(massh,tau2kd->getInitial4Vector());
    tau2kd->setLowerFitLimitE(tau2kd->getInitial4Vector());
    tau2kd->setUpperFitLimitE(massh,tau1kd->getInitial4Vector());

    HHFitObjectMET* metkd = new HHFitObjectMET(TVector2(Higgsgenerator.getMETwithsigma()[0],Higgsgenerator.getMETwithsigma()[1]));
    metkd->setCovMatrix(Higgsgenerator.getCovarmatrix()[0][0],Higgsgenerator.getCovarmatrix()[1][1],Higgsgenerator.getCovarmatrix()[1][0]);

    HHFitObject* higgskd  = new HHFitObjectComposite(tau1kd, tau2kd, metkd);

    HHFitConstraint* invmkd = new HHFitConstraintEHardM(tau1kd, tau2kd, massh);
    HHFitConstraint* balancekd = new HHFitConstraint4Vector(higgskd, true, true, false, false);
    HHFitConstraint* Likelihoodkd = new HHFitConstraintLikelihood(tau1kd,tau2kd,PDF1,PDF2);

    HHKinFit* singlefitkd = new HHKinFit();

    singlefitkd->addFitObjectE(tau1kd);
    singlefitkd->addConstraint(invmkd);
    singlefitkd->addConstraint(balancekd);
    singlefitkd->addConstraint(Likelihoodkd);

    try {
      singlefitkd->fit();
      if (!((singlefitkd->getConvergence()==1)||(singlefitkd->getConvergence()==2))){
	i--;
	continue;
      }		}
    catch(HHEnergyRangeException const& e){
      std::cout << i << std::endl;
      i--;
      continue;
    }

    std::cout <<i << "tau1 low- and up-frac likelihood "<< tau1kd->getInitial4Vector().E()/tau1kd->getLowerFitLimitE() <<" " <<tau1kd->getInitial4Vector().E()/tau1kd->getUpperFitLimitE() << std::endl;
    std::cout <<i<<"tau2 low- and up-frac likelihood "<< tau2kd->getInitial4Vector().E()/tau2kd->getLowerFitLimitE() <<" " <<tau2kd->getInitial4Vector().E()/tau2kd->getUpperFitLimitE() << std::endl;
TGraph* h_LkfromE=singlefitk->getLFunction(1000);
TGraph* h_LkdfromE=singlefitkd->getLFunction(1000);

h_LkfromE->Write();
h_LkdfromE->Write();

  }
  controlplots.Close();

  return(0);

}
