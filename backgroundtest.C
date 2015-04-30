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
  
  unsigned int numberofevents=1000000;
  TF1 *PDF1=new TF1("PDF1","2*x",0,1);
  TF1 *PDF2=new TF1("PDF2","2-2*x",0,1);
  TF1* pdf1=PDF1;
  TF1* pdf2=PDF2;
  TMatrixD covarmatrix(2,2);
  covarmatrix[0][0]=130;
  covarmatrix[0][1]=0;
  covarmatrix[1][0]=0;
  covarmatrix[1][1]=130;
  HHTauTauEventGenerator higgsgenerator(PDF1,PDF2,covarmatrix);
  higgsgenerator.setMhiggs(125.7);
  HHTauTauEventGenerator BackroundGenerator(PDF1,PDF1,covarmatrix);
  BackroundGenerator.setMhiggs(91.1876);
  HHTauTauEventGenerator higgsgenerator2(PDF1,PDF1,covarmatrix);
  higgsgenerator2.setMhiggs(125.7);
  HHTauTauEventGenerator BackroundGenerator2(PDF1,PDF2,covarmatrix);
  BackroundGenerator2.setMhiggs(91.1876);
  HHTauTauEventGenerator higgsgenerator3(PDF2,PDF1,covarmatrix);
  higgsgenerator3.setMhiggs(125.7);
  HHTauTauEventGenerator BackroundGenerator3(PDF2,PDF1,covarmatrix);
  BackroundGenerator3.setMhiggs(91.1876);
  HHTauTauEventGenerator higgsgenerator4(PDF2,PDF2,covarmatrix);
  higgsgenerator4.setMhiggs(125.7);
  HHTauTauEventGenerator BackroundGenerator4(PDF2,PDF2,covarmatrix);
  BackroundGenerator4.setMhiggs(91.1876);


  //Histogramms for h-fit
  TH1D h_FitFinalChi2h("h_FitFinalChi2h","The Final chi2 from the KinFit for higgsdecay",200,-5,20);
  TH1D h_FitFinalChi2probh("h_FitFinalChi2probh","The Final chi2 from the KinFit for higgsdecay",100,0,1);
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
  TH1D h_Likelihoodh("h_Likelihoodh","The likelihood of the fit",100,0,0.005);
  //Histogramms for h-fitwithlikelihood
  
  TH1D h_FitFinalChi2hlikelihood("h_FitFinalChi2hlikelihood","The Final chi2 from the likelihood-KinFit for higgsdecay",200,-5,20);
  TH1D h_fracresulutiontau1hlikelihood("h_fracresulutiontau1hlikelihood","resulution of the energyfraction from tau1vis from the likelihood-KinFit for higgsdecay",100,-1,1);
  TH1D h_fracresulutiontau2hlikelihood("h_fracresulutiontau2hlikelihood","resulution of the energyfraction from tau2vis from the likelihood-KinFit for higgsdecay",100,-1,1);
  TH1D h_fracresulutiontau1weightedpdfhlikelihood("h_fracresulutiontau1weightedpdfhlikelihood","weighted(pdf) resulution of the energyfraction from tau1vis from the likelihood-KinFit for higgsdecay",100,-1,1);
  TH1D h_fracresulutiontau2weightedpdfhlikelihood("h_fracresulutiontau2weightedpdfhlikelihood","weighted(pdf) resulution of the energyfraction from tau2vis from the likelihood-KinFit for higgsdecay",100,-1,1);
  TH1D h_energyresulution1hlikelihood("h_energyresulution1hlikelihood","the energyresulution of tau 1 in the likelihood-KinFit for higgsdecay",100,-1.3,1.3);
  TH1D h_energyresulution2hlikelihood("h_energyresulution2hlikelihood","the energyresulution of tau 2 in the likelihood-KinFit for higgsdecay",100,-1.3,1.3);
  TH1D h_energyresulution1weightedpdfhlikelihood("h_energyresulution1weightedpdfhlikelihood","the weighted(pdf) energyresulution of tau 1 in the likelihood-KinFit for higgsdecay",100,-1.3,1.3);
  TH1D h_energyresulution2weightedpdfhlikelihood("h_energyresulution2weightedpdfhlikelihood","the weighted(pdf) energyresulution of tau 2 in the likelihood-KinFit for higgsdecay",100,-1.3,1.3);
  TH1D h_Likelihoodhlikelihood("h_Likelihoodhlikelihood","The likelihood of the fit",100,0,0.005);
  
  //Histogramms for background-fit
  TH1D h_FitFinalChi2b("h_FitFinalChi2b","The Final chi2 from the KinFit for background-events",200,-5,20);
  TH1D h_FitFinalChi2probb("h_FitFinalChi2probb","The Final chi2 from the KinFit for background-events",100,0,1);
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
  TH1D h_Likelihoodb("h_Likelihoodb","The likelihood of the fit",100,0,0.005);
  //Histogramms for background-fitwithlikelihood
  
  TH1D h_FitFinalChi2blikelihood("h_FitFinalChi2blikelihood","The Final chi2 from the likelihood-KinFit for background-events",200,-5,20);
  TH1D h_fracresulutiontau1blikelihood("h_fracresulutiontau1blikelihood","resulution of the energyfraction from tau1vis from the likelihood-KinFit for background-events",100,-1,1);
  TH1D h_fracresulutiontau2blikelihood("h_fracresulutiontau2blikelihood","resulution of the energyfraction from tau2vis from the likelihood-KinFit for background-events",100,-1,1);
  TH1D h_fracresulutiontau1weightedpdfblikelihood("h_fracresulutiontau1weightedpdfblikelihood","weighted(pdf) resulution of the energyfraction from tau1vis from the likelihood-KinFit for background-events",100,-1,1);
  TH1D h_fracresulutiontau2weightedpdfblikelihood("h_fracresulutiontau2weightedpdfblikelihood","weighted(pdf) resulution of the energyfraction from tau2vis from the likelihood-KinFit for background-events",100,-1,1);
  TH1D h_energyresulution1blikelihood("h_energyresulution1blikelihood","the energyresulution of tau 1 in the likelihood-KinFit for background-events",100,-1.3,1.3);
  TH1D h_energyresulution2blikelihood("h_energyresulution2blikelihood","the energyresulution of tau 2 in the likelihood-KinFit for background-events",100,-1.3,1.3);
  TH1D h_energyresulution1weightedpdfblikelihood("h_energyresulution1weightedpdfblikelihood","the weighted(pdf) energyresulution of tau 1 in the likelihood-KinFit for background-events",100,-1.3,1.3);
  TH1D h_energyresulution2weightedpdfblikelihood("h_energyresulution2weightedpdfblikelihood","the weighted(pdf) energyresulution of tau 2 in the likelihood-KinFit for background-events",100,-1.3,1.3);
  TH1D h_Likelihoodblikelihood("h_Likelihoodblikelihood","The likelihood of the fit",100,0,0.005);
  //---------------------------------------------------------------------------------------------------------------------------------------------------
  // alternative events:
  //Histogramms for h-fit with same pdf for both tau
  TH1D h_FitFinalChi2halt("h_FitFinalChi2halt","The Final chi2 from the KinFit for higgsdecay",200,-5,20);
  TH1D h_FitFinalChi2probhalt("h_FitFinalChi2probhalt","The Final chi2 from the KinFit for higgsdecay",100,0,1);
  TH1D h_fracresulutiontau1halt("h_fracresulutiontau1halt","resulution of the energyfraction from tau1vis from the KinFit for higgsdecay",100,-1,1);
  TH1D h_fracresulutiontau2halt("h_fracresulutiontau2halt","resulution of the energyfraction from tau2vis from the KinFit for higgsdecay",100,-1,1);
  TH1D h_fracresulutiontau1weightedprobhalt("h_fracresulutiontau1weightedprobhalt","weighted(chi2prob) resulution of the energyfraction from tau1vis from the KinFit for higgsdecay",100,-1,1);
  TH1D h_fracresulutiontau2weightedprobhalt("h_fracresulutiontau2weightedprobhalt","weighted(chi2prob) resulution of the energyfraction from tau2vis from the KinFit for higgsdecay",100,-1,1);
  TH1D h_fracresulutiontau1weightedpdfhalt("h_fracresulutiontau1weightedpdfhalt","weighted(pdf) resulution of the energyfraction from tau1vis from the KinFit for higgsdecay",100,-1,1);
  TH1D h_fracresulutiontau2weightedpdfhalt("h_fracresulutiontau2weightedpdfhalt","weighted(pdf) resulution of the energyfraction from tau2vis from the KinFit for higgsdecay",100,-1,1);
  TH1D h_energyresulution1halt("h_energyresulution1halt","the energyresulution of tau 1 in the KinFit for higgsdecay",100,-1.3,1.3);
  TH1D h_energyresulution2halt("h_energyresulution2halt","the energyresulution of tau 2 in the KinFit for higgsdecay",100,-1.3,1.3);
  TH1D h_energyresulution1weightedprobhalt("h_energyresulution1weightedprobhalt","the weighted(chi2prob) energyresulution of tau 1 in the KinFit for higgsdecay",100,-1.3,1.3);
  TH1D h_energyresulution2weightedprobhalt("h_energyresulution2weightedprobhalt","the weighted(chi2prob) energyresulution of tau 2 in the KinFit for higgsdecay",100,-1.3,1.3);
  TH1D h_energyresulution1weightedpdfhalt("h_energyresulution1weightedpdfhalt","the weighted(pdf) energyresulution of tau 1 in the KinFit for higgsdecay",100,-1.3,1.3);
  TH1D h_energyresulution2weightedpdfhalt("h_energyresulution2weightedpdfhalt","the weighted(pdf) energyresulution of tau 2 in the KinFit for higgsdecay",100,-1.3,1.3);
  TH1D h_Likelihoodhalt("h_Likelihoodhalt","The likelihood of the fit",100,0,0.005);
  //Histogramms for h-fitwithlikelihood with same pdf for both tau

  TH1D h_FitFinalChi2haltlikelihood("h_FitFinalChi2haltlikelihood","The Final chi2 from the likelihood-KinFit for higgsdecay",200,-5,20);
  TH1D h_fracresulutiontau1haltlikelihood("h_fracresulutiontau1haltlikelihood","resulution of the energyfraction from tau1vis from the likelihood-KinFit for higgsdecay",100,-1,1);
  TH1D h_fracresulutiontau2haltlikelihood("h_fracresulutiontau2haltlikelihood","resulution of the energyfraction from tau2vis from the likelihood-KinFit for higgsdecay",100,-1,1);
  TH1D h_fracresulutiontau1weightedpdfhaltlikelihood("h_fracresulutiontau1weightedpdfhaltlikelihood","weighted(pdf) resulution of the energyfraction from tau1vis from the likelihood-KinFit for higgsdecay",100,-1,1);
  TH1D h_fracresulutiontau2weightedpdfhaltlikelihood("h_fracresulutiontau2weightedpdfhaltlikelihood","weighted(pdf) resulution of the energyfraction from tau2vis from the likelihood-KinFit for higgsdecay",100,-1,1);
  TH1D h_energyresulution1haltlikelihood("h_energyresulution1haltlikelihood","the energyresulution of tau 1 in the likelihood-KinFit for higgsdecay",100,-1.3,1.3);
  TH1D h_energyresulution2haltlikelihood("h_energyresulution2haltlikelihood","the energyresulution of tau 2 in the likelihood-KinFit for higgsdecay",100,-1.3,1.3);
  TH1D h_energyresulution1weightedpdfhaltlikelihood("h_energyresulution1weightedpdfhaltlikelihood","the weighted(pdf) energyresulution of tau 1 in the likelihood-KinFit for higgsdecay",100,-1.3,1.3);
  TH1D h_energyresulution2weightedpdfhaltlikelihood("h_energyresulution2weightedpdfhaltlikelihood","the weighted(pdf) energyresulution of tau 2 in the likelihood-KinFit for higgsdecay",100,-1.3,1.3);
  TH1D h_Likelihoodhaltlikelihood("h_Likelihoodhaltlikelihood","The likelihood of the fit",100,0,0.005);

  //Histogramms for background-fit with higgs pdfs
  TH1D h_FitFinalChi2balt("h_FitFinalChi2balt","The Final chi2 from the KinFit for background-events",200,-5,20);
  TH1D h_FitFinalChi2probbalt("h_FitFinalChi2probbalt","The Final chi2 from the KinFit for background-events",100,0,1);
  TH1D h_fracresulutiontau1balt("h_fracresulutiontau1balt","resulution of the energyfraction from tau1vis from the KinFit for background-events",100,-1,1);
  TH1D h_fracresulutiontau2balt("h_fracresulutiontau2balt","resulution of the energyfraction from tau2vis from the KinFit for background-events",100,-1,1);
  TH1D h_fracresulutiontau1weightedprobbalt("h_fracresulutiontau1weightedprobbalt","weighted(chi2prob) resulution of the energyfraction from tau1vis from the KinFit for background-events",100,-1,1);
  TH1D h_fracresulutiontau2weightedprobbalt("h_fracresulutiontau2weightedprobbalt","weighted(chi2prob) resulution of the energyfraction from tau2vis from the KinFit for background-events",100,-1,1);
  TH1D h_fracresulutiontau1weightedpdfbalt("h_fracresulutiontau1weightedpdfbalt","weighted(pdf) resulution of the energyfraction from tau1vis from the KinFit for background-events",100,-1,1);
  TH1D h_fracresulutiontau2weightedpdfbalt("h_fracresulutiontau2weightedpdfbalt","weighted(pdf) resulution of the energyfraction from tau2vis from the KinFit for background-events",100,-1,1);
  TH1D h_energyresulution1balt("h_energyresulution1balt","the energyresulution of tau 1 in the KinFit for background-events",100,-1.3,1.3);
  TH1D h_energyresulution2balt("h_energyresulution2balt","the energyresulution of tau 2 in the KinFit for background-events",100,-1.3,1.3);
  TH1D h_energyresulution1weightedprobbalt("h_energyresulution1weightedprobbalt","the weighted(chi2prob) energyresulution of tau 1 in the KinFit for background-events",100,-1.3,1.3);
  TH1D h_energyresulution2weightedprobbalt("h_energyresulution2weightedprobbalt","the weighted(chi2prob) energyresulution of tau 2 in the KinFit for background-events",100,-1.3,1.3);
  TH1D h_energyresulution1weightedpdfbalt("h_energyresulution1weightedpdfbalt","the weighted(pdf) energyresulution of tau 1 in the KinFit for background-events",100,-1.3,1.3);
  TH1D h_energyresulution2weightedpdfbalt("h_energyresulution2weightedpdfbalt","the weighted(pdf) energyresulution of tau 2 in the KinFit for background-events",100,-1.3,1.3);
  TH1D h_Likelihoodbalt("h_Likelihoodbalt","The likelihood of the fit",100,0,0.005);
  //Histogramms for background-fitwithlikelihood with higgs pdfs

  TH1D h_FitFinalChi2baltlikelihood("h_FitFinalChi2baltlikelihood","The Final chi2 from the likelihood-KinFit for background-events",200,-5,20);
  TH1D h_fracresulutiontau1baltlikelihood("h_fracresulutiontau1baltlikelihood","resulution of the energyfraction from tau1vis from the likelihood-KinFit for background-events",100,-1,1);
  TH1D h_fracresulutiontau2baltlikelihood("h_fracresulutiontau2baltlikelihood","resulution of the energyfraction from tau2vis from the likelihood-KinFit for background-events",100,-1,1);
  TH1D h_fracresulutiontau1weightedpdfbaltlikelihood("h_fracresulutiontau1weightedpdfbaltlikelihood","weighted(pdf) resulution of the energyfraction from tau1vis from the likelihood-KinFit for background-events",100,-1,1);
  TH1D h_fracresulutiontau2weightedpdfbaltlikelihood("h_fracresulutiontau2weightedpdfbaltlikelihood","weighted(pdf) resulution of the energyfraction from tau2vis from the likelihood-KinFit for background-events",100,-1,1);
  TH1D h_energyresulution1baltlikelihood("h_energyresulution1baltlikelihood","the energyresulution of tau 1 in the likelihood-KinFit for background-events",100,-1.3,1.3);
  TH1D h_energyresulution2baltlikelihood("h_energyresulution2baltlikelihood","the energyresulution of tau 2 in the likelihood-KinFit for background-events",100,-1.3,1.3);
  TH1D h_energyresulution1weightedpdfbaltlikelihood("h_energyresulution1weightedpdfbaltlikelihood","the weighted(pdf) energyresulution of tau 1 in the likelihood-KinFit for background-events",100,-1.3,1.3);
  TH1D h_energyresulution2weightedpdfbaltlikelihood("h_energyresulution2weightedpdfbaltlikelihood","the weighted(pdf) energyresulution of tau 2 in the likelihood-KinFit for background-events",100,-1.3,1.3);
  TH1D h_Likelihoodbaltlikelihood("h_Likelihoodbaltlikelihood","The likelihood of the fit",100,0,0.005);
  //-------------------------------------------------------------------------------------------------------------------------------------------------
  // alt2ernative events2:
  //Histogramms for h-fit with same pdf for both tau
  TH1D h_FitFinalChi2halt2("h_FitFinalChi2halt2","The Final chi2 from the KinFit for higgsdecay",200,-5,20);
  TH1D h_FitFinalChi2probhalt2("h_FitFinalChi2probhalt2","The Final chi2 from the KinFit for higgsdecay",100,0,1);
  TH1D h_fracresulutiontau1halt2("h_fracresulutiontau1halt2","resulution of the energyfraction from tau1vis from the KinFit for higgsdecay",100,-1,1);
  TH1D h_fracresulutiontau2halt2("h_fracresulutiontau2halt2","resulution of the energyfraction from tau2vis from the KinFit for higgsdecay",100,-1,1);
  TH1D h_fracresulutiontau1weightedprobhalt2("h_fracresulutiontau1weightedprobhalt2","weighted(chi2prob) resulution of the energyfraction from tau1vis from the KinFit for higgsdecay",100,-1,1);
  TH1D h_fracresulutiontau2weightedprobhalt2("h_fracresulutiontau2weightedprobhalt2","weighted(chi2prob) resulution of the energyfraction from tau2vis from the KinFit for higgsdecay",100,-1,1);
  TH1D h_fracresulutiontau1weightedpdfhalt2("h_fracresulutiontau1weightedpdfhalt2","weighted(pdf) resulution of the energyfraction from tau1vis from the KinFit for higgsdecay",100,-1,1);
  TH1D h_fracresulutiontau2weightedpdfhalt2("h_fracresulutiontau2weightedpdfhalt2","weighted(pdf) resulution of the energyfraction from tau2vis from the KinFit for higgsdecay",100,-1,1);
  TH1D h_energyresulution1halt2("h_energyresulution1halt2","the energyresulution of tau 1 in the KinFit for higgsdecay",100,-1.3,1.3);
  TH1D h_energyresulution2halt2("h_energyresulution2halt2","the energyresulution of tau 2 in the KinFit for higgsdecay",100,-1.3,1.3);
  TH1D h_energyresulution1weightedprobhalt2("h_energyresulution1weightedprobhalt2","the weighted(chi2prob) energyresulution of tau 1 in the KinFit for higgsdecay",100,-1.3,1.3);
  TH1D h_energyresulution2weightedprobhalt2("h_energyresulution2weightedprobhalt2","the weighted(chi2prob) energyresulution of tau 2 in the KinFit for higgsdecay",100,-1.3,1.3);
  TH1D h_energyresulution1weightedpdfhalt2("h_energyresulution1weightedpdfhalt2","the weighted(pdf) energyresulution of tau 1 in the KinFit for higgsdecay",100,-1.3,1.3);
  TH1D h_energyresulution2weightedpdfhalt2("h_energyresulution2weightedpdfhalt2","the weighted(pdf) energyresulution of tau 2 in the KinFit for higgsdecay",100,-1.3,1.3);
  TH1D h_Likelihoodhalt2("h_Likelihoodhalt2","The likelihood of the fit",100,0,0.005);
  //Histogramms for h-fitwithlikelihood with same pdf for both tau

  TH1D h_FitFinalChi2halt2likelihood("h_FitFinalChi2halt2likelihood","The Final chi2 from the likelihood-KinFit for higgsdecay",200,-5,20);
  TH1D h_fracresulutiontau1halt2likelihood("h_fracresulutiontau1halt2likelihood","resulution of the energyfraction from tau1vis from the likelihood-KinFit for higgsdecay",100,-1,1);
  TH1D h_fracresulutiontau2halt2likelihood("h_fracresulutiontau2halt2likelihood","resulution of the energyfraction from tau2vis from the likelihood-KinFit for higgsdecay",100,-1,1);
  TH1D h_fracresulutiontau1weightedpdfhalt2likelihood("h_fracresulutiontau1weightedpdfhalt2likelihood","weighted(pdf) resulution of the energyfraction from tau1vis from the likelihood-KinFit for higgsdecay",100,-1,1);
  TH1D h_fracresulutiontau2weightedpdfhalt2likelihood("h_fracresulutiontau2weightedpdfhalt2likelihood","weighted(pdf) resulution of the energyfraction from tau2vis from the likelihood-KinFit for higgsdecay",100,-1,1);
  TH1D h_energyresulution1halt2likelihood("h_energyresulution1halt2likelihood","the energyresulution of tau 1 in the likelihood-KinFit for higgsdecay",100,-1.3,1.3);
  TH1D h_energyresulution2halt2likelihood("h_energyresulution2halt2likelihood","the energyresulution of tau 2 in the likelihood-KinFit for higgsdecay",100,-1.3,1.3);
  TH1D h_energyresulution1weightedpdfhalt2likelihood("h_energyresulution1weightedpdfhalt2likelihood","the weighted(pdf) energyresulution of tau 1 in the likelihood-KinFit for higgsdecay",100,-1.3,1.3);
  TH1D h_energyresulution2weightedpdfhalt2likelihood("h_energyresulution2weightedpdfhalt2likelihood","the weighted(pdf) energyresulution of tau 2 in the likelihood-KinFit for higgsdecay",100,-1.3,1.3);
  TH1D h_Likelihoodhalt2likelihood("h_Likelihoodhalt2likelihood","The likelihood of the fit",100,0,0.005);

  //Histogramms for background-fit with higgs pdfs
  TH1D h_FitFinalChi2balt2("h_FitFinalChi2balt2","The Final chi2 from the KinFit for background-events",200,-5,20);
  TH1D h_FitFinalChi2probbalt2("h_FitFinalChi2probbalt2","The Final chi2 from the KinFit for background-events",100,0,1);
  TH1D h_fracresulutiontau1balt2("h_fracresulutiontau1balt2","resulution of the energyfraction from tau1vis from the KinFit for background-events",100,-1,1);
  TH1D h_fracresulutiontau2balt2("h_fracresulutiontau2balt2","resulution of the energyfraction from tau2vis from the KinFit for background-events",100,-1,1);
  TH1D h_fracresulutiontau1weightedprobbalt2("h_fracresulutiontau1weightedprobbalt2","weighted(chi2prob) resulution of the energyfraction from tau1vis from the KinFit for background-events",100,-1,1);
  TH1D h_fracresulutiontau2weightedprobbalt2("h_fracresulutiontau2weightedprobbalt2","weighted(chi2prob) resulution of the energyfraction from tau2vis from the KinFit for background-events",100,-1,1);
  TH1D h_fracresulutiontau1weightedpdfbalt2("h_fracresulutiontau1weightedpdfbalt2","weighted(pdf) resulution of the energyfraction from tau1vis from the KinFit for background-events",100,-1,1);
  TH1D h_fracresulutiontau2weightedpdfbalt2("h_fracresulutiontau2weightedpdfbalt2","weighted(pdf) resulution of the energyfraction from tau2vis from the KinFit for background-events",100,-1,1);
  TH1D h_energyresulution1balt2("h_energyresulution1balt2","the energyresulution of tau 1 in the KinFit for background-events",100,-1.3,1.3);
  TH1D h_energyresulution2balt2("h_energyresulution2balt2","the energyresulution of tau 2 in the KinFit for background-events",100,-1.3,1.3);
  TH1D h_energyresulution1weightedprobbalt2("h_energyresulution1weightedprobbalt2","the weighted(chi2prob) energyresulution of tau 1 in the KinFit for background-events",100,-1.3,1.3);
  TH1D h_energyresulution2weightedprobbalt2("h_energyresulution2weightedprobbalt2","the weighted(chi2prob) energyresulution of tau 2 in the KinFit for background-events",100,-1.3,1.3);
  TH1D h_energyresulution1weightedpdfbalt2("h_energyresulution1weightedpdfbalt2","the weighted(pdf) energyresulution of tau 1 in the KinFit for background-events",100,-1.3,1.3);
  TH1D h_energyresulution2weightedpdfbalt2("h_energyresulution2weightedpdfbalt2","the weighted(pdf) energyresulution of tau 2 in the KinFit for background-events",100,-1.3,1.3);
  TH1D h_Likelihoodbalt2("h_Likelihoodbalt2","The likelihood of the fit",100,0,0.005);
  //Histogramms for background-fitwithlikelihood with higgs pdfs

  TH1D h_FitFinalChi2balt2likelihood("h_FitFinalChi2balt2likelihood","The Final chi2 from the likelihood-KinFit for background-events",200,-5,20);
  TH1D h_fracresulutiontau1balt2likelihood("h_fracresulutiontau1balt2likelihood","resulution of the energyfraction from tau1vis from the likelihood-KinFit for background-events",100,-1,1);
  TH1D h_fracresulutiontau2balt2likelihood("h_fracresulutiontau2balt2likelihood","resulution of the energyfraction from tau2vis from the likelihood-KinFit for background-events",100,-1,1);
  TH1D h_fracresulutiontau1weightedpdfbalt2likelihood("h_fracresulutiontau1weightedpdfbalt2likelihood","weighted(pdf) resulution of the energyfraction from tau1vis from the likelihood-KinFit for background-events",100,-1,1);
  TH1D h_fracresulutiontau2weightedpdfbalt2likelihood("h_fracresulutiontau2weightedpdfbalt2likelihood","weighted(pdf) resulution of the energyfraction from tau2vis from the likelihood-KinFit for background-events",100,-1,1);
  TH1D h_energyresulution1balt2likelihood("h_energyresulution1balt2likelihood","the energyresulution of tau 1 in the likelihood-KinFit for background-events",100,-1.3,1.3);
  TH1D h_energyresulution2balt2likelihood("h_energyresulution2balt2likelihood","the energyresulution of tau 2 in the likelihood-KinFit for background-events",100,-1.3,1.3);
  TH1D h_energyresulution1weightedpdfbalt2likelihood("h_energyresulution1weightedpdfbalt2likelihood","the weighted(pdf) energyresulution of tau 1 in the likelihood-KinFit for background-events",100,-1.3,1.3);
  TH1D h_energyresulution2weightedpdfbalt2likelihood("h_energyresulution2weightedpdfbalt2likelihood","the weighted(pdf) energyresulution of tau 2 in the likelihood-KinFit for background-events",100,-1.3,1.3);
  TH1D h_Likelihoodbalt2likelihood("h_Likelihoodbalt2likelihood","The likelihood of the fit",100,0,0.005);
  //-------------------------------------------------------------------------------------------------------------------------------------------------
  // alt3ernative events2:
  //Histogramms for h-fit with same pdf for both tau
  TH1D h_FitFinalChi2halt3("h_FitFinalChi2halt3","The Final chi2 from the KinFit for higgsdecay",200,-5,20);
  TH1D h_FitFinalChi2probhalt3("h_FitFinalChi2probhalt3","The Final chi2 from the KinFit for higgsdecay",100,0,1);
  TH1D h_fracresulutiontau1halt3("h_fracresulutiontau1halt3","resulution of the energyfraction from tau1vis from the KinFit for higgsdecay",100,-1,1);
  TH1D h_fracresulutiontau2halt3("h_fracresulutiontau2halt3","resulution of the energyfraction from tau2vis from the KinFit for higgsdecay",100,-1,1);
  TH1D h_fracresulutiontau1weightedprobhalt3("h_fracresulutiontau1weightedprobhalt3","weighted(chi2prob) resulution of the energyfraction from tau1vis from the KinFit for higgsdecay",100,-1,1);
  TH1D h_fracresulutiontau2weightedprobhalt3("h_fracresulutiontau2weightedprobhalt3","weighted(chi2prob) resulution of the energyfraction from tau2vis from the KinFit for higgsdecay",100,-1,1);
  TH1D h_fracresulutiontau1weightedpdfhalt3("h_fracresulutiontau1weightedpdfhalt3","weighted(pdf) resulution of the energyfraction from tau1vis from the KinFit for higgsdecay",100,-1,1);
  TH1D h_fracresulutiontau2weightedpdfhalt3("h_fracresulutiontau2weightedpdfhalt3","weighted(pdf) resulution of the energyfraction from tau2vis from the KinFit for higgsdecay",100,-1,1);
  TH1D h_energyresulution1halt3("h_energyresulution1halt3","the energyresulution of tau 1 in the KinFit for higgsdecay",100,-1.3,1.3);
  TH1D h_energyresulution2halt3("h_energyresulution2halt3","the energyresulution of tau 2 in the KinFit for higgsdecay",100,-1.3,1.3);
  TH1D h_energyresulution1weightedprobhalt3("h_energyresulution1weightedprobhalt3","the weighted(chi2prob) energyresulution of tau 1 in the KinFit for higgsdecay",100,-1.3,1.3);
  TH1D h_energyresulution2weightedprobhalt3("h_energyresulution2weightedprobhalt3","the weighted(chi2prob) energyresulution of tau 2 in the KinFit for higgsdecay",100,-1.3,1.3);
  TH1D h_energyresulution1weightedpdfhalt3("h_energyresulution1weightedpdfhalt3","the weighted(pdf) energyresulution of tau 1 in the KinFit for higgsdecay",100,-1.3,1.3);
  TH1D h_energyresulution2weightedpdfhalt3("h_energyresulution2weightedpdfhalt3","the weighted(pdf) energyresulution of tau 2 in the KinFit for higgsdecay",100,-1.3,1.3);
  TH1D h_Likelihoodhalt3("h_Likelihoodhalt3","The likelihood of the fit",100,0,0.005);
  //Histogramms for h-fitwithlikelihood with same pdf for both tau

  TH1D h_FitFinalChi2halt3likelihood("h_FitFinalChi2halt3likelihood","The Final chi2 from the likelihood-KinFit for higgsdecay",200,-5,20);
  TH1D h_fracresulutiontau1halt3likelihood("h_fracresulutiontau1halt3likelihood","resulution of the energyfraction from tau1vis from the likelihood-KinFit for higgsdecay",100,-1,1);
  TH1D h_fracresulutiontau2halt3likelihood("h_fracresulutiontau2halt3likelihood","resulution of the energyfraction from tau2vis from the likelihood-KinFit for higgsdecay",100,-1,1);
  TH1D h_fracresulutiontau1weightedpdfhalt3likelihood("h_fracresulutiontau1weightedpdfhalt3likelihood","weighted(pdf) resulution of the energyfraction from tau1vis from the likelihood-KinFit for higgsdecay",100,-1,1);
  TH1D h_fracresulutiontau2weightedpdfhalt3likelihood("h_fracresulutiontau2weightedpdfhalt3likelihood","weighted(pdf) resulution of the energyfraction from tau2vis from the likelihood-KinFit for higgsdecay",100,-1,1);
  TH1D h_energyresulution1halt3likelihood("h_energyresulution1halt3likelihood","the energyresulution of tau 1 in the likelihood-KinFit for higgsdecay",100,-1.3,1.3);
  TH1D h_energyresulution2halt3likelihood("h_energyresulution2halt3likelihood","the energyresulution of tau 2 in the likelihood-KinFit for higgsdecay",100,-1.3,1.3);
  TH1D h_energyresulution1weightedpdfhalt3likelihood("h_energyresulution1weightedpdfhalt3likelihood","the weighted(pdf) energyresulution of tau 1 in the likelihood-KinFit for higgsdecay",100,-1.3,1.3);
  TH1D h_energyresulution2weightedpdfhalt3likelihood("h_energyresulution2weightedpdfhalt3likelihood","the weighted(pdf) energyresulution of tau 2 in the likelihood-KinFit for higgsdecay",100,-1.3,1.3);
  TH1D h_Likelihoodhalt3likelihood("h_Likelihoodhalt3likelihood","The likelihood of the fit",100,0,0.005);

  //Histogramms for background-fit with higgs pdfs
  TH1D h_FitFinalChi2balt3("h_FitFinalChi2balt3","The Final chi2 from the KinFit for background-events",200,-5,20);
  TH1D h_FitFinalChi2probbalt3("h_FitFinalChi2probbalt3","The Final chi2 from the KinFit for background-events",100,0,1);
  TH1D h_fracresulutiontau1balt3("h_fracresulutiontau1balt3","resulution of the energyfraction from tau1vis from the KinFit for background-events",100,-1,1);
  TH1D h_fracresulutiontau2balt3("h_fracresulutiontau2balt3","resulution of the energyfraction from tau2vis from the KinFit for background-events",100,-1,1);
  TH1D h_fracresulutiontau1weightedprobbalt3("h_fracresulutiontau1weightedprobbalt3","weighted(chi2prob) resulution of the energyfraction from tau1vis from the KinFit for background-events",100,-1,1);
  TH1D h_fracresulutiontau2weightedprobbalt3("h_fracresulutiontau2weightedprobbalt3","weighted(chi2prob) resulution of the energyfraction from tau2vis from the KinFit for background-events",100,-1,1);
  TH1D h_fracresulutiontau1weightedpdfbalt3("h_fracresulutiontau1weightedpdfbalt3","weighted(pdf) resulution of the energyfraction from tau1vis from the KinFit for background-events",100,-1,1);
  TH1D h_fracresulutiontau2weightedpdfbalt3("h_fracresulutiontau2weightedpdfbalt3","weighted(pdf) resulution of the energyfraction from tau2vis from the KinFit for background-events",100,-1,1);
  TH1D h_energyresulution1balt3("h_energyresulution1balt3","the energyresulution of tau 1 in the KinFit for background-events",100,-1.3,1.3);
  TH1D h_energyresulution2balt3("h_energyresulution2balt3","the energyresulution of tau 2 in the KinFit for background-events",100,-1.3,1.3);
  TH1D h_energyresulution1weightedprobbalt3("h_energyresulution1weightedprobbalt3","the weighted(chi2prob) energyresulution of tau 1 in the KinFit for background-events",100,-1.3,1.3);
  TH1D h_energyresulution2weightedprobbalt3("h_energyresulution2weightedprobbalt3","the weighted(chi2prob) energyresulution of tau 2 in the KinFit for background-events",100,-1.3,1.3);
  TH1D h_energyresulution1weightedpdfbalt3("h_energyresulution1weightedpdfbalt3","the weighted(pdf) energyresulution of tau 1 in the KinFit for background-events",100,-1.3,1.3);
  TH1D h_energyresulution2weightedpdfbalt3("h_energyresulution2weightedpdfbalt3","the weighted(pdf) energyresulution of tau 2 in the KinFit for background-events",100,-1.3,1.3);
  TH1D h_Likelihoodbalt3("h_Likelihoodbalt3","The likelihood of the fit",100,0,0.005);
  //Histogramms for background-fitwithlikelihood with higgs pdfs

  TH1D h_FitFinalChi2balt3likelihood("h_FitFinalChi2balt3likelihood","The Final chi2 from the likelihood-KinFit for background-events",200,-5,20);
  TH1D h_fracresulutiontau1balt3likelihood("h_fracresulutiontau1balt3likelihood","resulution of the energyfraction from tau1vis from the likelihood-KinFit for background-events",100,-1,1);
  TH1D h_fracresulutiontau2balt3likelihood("h_fracresulutiontau2balt3likelihood","resulution of the energyfraction from tau2vis from the likelihood-KinFit for background-events",100,-1,1);
  TH1D h_fracresulutiontau1weightedpdfbalt3likelihood("h_fracresulutiontau1weightedpdfbalt3likelihood","weighted(pdf) resulution of the energyfraction from tau1vis from the likelihood-KinFit for background-events",100,-1,1);
  TH1D h_fracresulutiontau2weightedpdfbalt3likelihood("h_fracresulutiontau2weightedpdfbalt3likelihood","weighted(pdf) resulution of the energyfraction from tau2vis from the likelihood-KinFit for background-events",100,-1,1);
  TH1D h_energyresulution1balt3likelihood("h_energyresulution1balt3likelihood","the energyresulution of tau 1 in the likelihood-KinFit for background-events",100,-1.3,1.3);
  TH1D h_energyresulution2balt3likelihood("h_energyresulution2balt3likelihood","the energyresulution of tau 2 in the likelihood-KinFit for background-events",100,-1.3,1.3);
  TH1D h_energyresulution1weightedpdfbalt3likelihood("h_energyresulution1weightedpdfbalt3likelihood","the weighted(pdf) energyresulution of tau 1 in the likelihood-KinFit for background-events",100,-1.3,1.3);
  TH1D h_energyresulution2weightedpdfbalt3likelihood("h_energyresulution2weightedpdfbalt3likelihood","the weighted(pdf) energyresulution of tau 2 in the likelihood-KinFit for background-events",100,-1.3,1.3);
  TH1D h_Likelihoodbalt3likelihood("h_Likelihoodbalt3likelihood","The likelihood of the fit",100,0,0.005);

  
  TH1D h_ztest("h_ztest","h_ztest",100,80,100);

  for(unsigned int i=0; i<numberofevents; i++){
    try{
      higgsgenerator.generateEvent();
    }
    catch(const HHEnergyRangeException& e){
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
      if (!((singlefith->getConvergence()==1)||(singlefith->getConvergence()==2))){
	i--;
	continue;
      }
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

    h_fracresulutiontau1weightedpdfh.Fill(comparefrac1h,PDF1->Eval(fitfractau1h));
    h_fracresulutiontau2weightedpdfh.Fill(comparefrac2h,PDF2->Eval(fitfractau2h));
    h_energyresulution1weightedprobh.Fill((higgsgenerator.getTau1boosted().E()-tau1h->getFit4Vector().E())/higgsgenerator.getTau1boosted().E(),TMath::Prob(singlefith->getChi2(),1));
    h_energyresulution2weightedprobh.Fill((higgsgenerator.getTau2boosted().E()-tau2h->getFit4Vector().E())/higgsgenerator.getTau2boosted().E(),TMath::Prob(singlefith->getChi2(),1));
    h_energyresulution1weightedpdfh.Fill((higgsgenerator.getTau1boosted().E()-tau1h->getFit4Vector().E())/higgsgenerator.getTau1boosted().E(),PDF1->Eval(fitfractau1h));
    h_energyresulution2weightedpdfh.Fill((higgsgenerator.getTau2boosted().E()-tau2h->getFit4Vector().E())/higgsgenerator.getTau2boosted().E(),PDF2->Eval(fitfractau2h));
    h_Likelihoodh.Fill(balanceh->getLikelihood()*PDF1->Eval(fitfractau1h)*PDF2->Eval(fitfractau2h));

  }
  for(unsigned int i=0; i<numberofevents; i++){
    try{
      higgsgenerator.generateEvent();
    }
    catch(const HHEnergyRangeException& e){
      i--;
      continue;
    }


    //with likelihood-----------------------------------------------------------------------------------------------------------------------------------
    double mass = higgsgenerator.getMhiggs();
    //prepare tau objects
    HHFitObjectE* tau1hlikelihood = new HHFitObjectEConstM(higgsgenerator.getTau1Vis());//  visible Tau1-component from HHTauTauEventGenerator
    HHFitObjectE* tau2hlikelihood = new HHFitObjectEConstM(higgsgenerator.getTau2Vis());//  visible Tau2-component from HHTauTauEventGenerator
    
    //prepare MET object
    HHFitObjectMET* methlikelihood = new HHFitObjectMET(TVector2(higgsgenerator.getMETwithsigma()[0],higgsgenerator.getMETwithsigma()[1]));//Use Met components from HHTauTauEventGenerator
    methlikelihood->setCovMatrix(higgsgenerator.getCovarmatrix()[0][0],higgsgenerator.getCovarmatrix()[1][1],higgsgenerator.getCovarmatrix()[1][0]); // set Covarmatrix with Matrix inserted in HHTauTauEventGenerator
    
    //prepare composite object: Higgs
    HHFitObject* higgshlikelihood  = new HHFitObjectComposite(tau1hlikelihood, tau2hlikelihood, methlikelihood);
    
    tau1hlikelihood->setLowerFitLimitE(tau1hlikelihood->getInitial4Vector());
    tau1hlikelihood->setUpperFitLimitE(mass,tau2hlikelihood->getInitial4Vector());
    tau2hlikelihood->setLowerFitLimitE(tau2hlikelihood->getInitial4Vector());
    tau2hlikelihood->setUpperFitLimitE(mass,tau1hlikelihood->getInitial4Vector());
    
    //prepare constraints
    HHFitConstraint* invmhlikelihood = new HHFitConstraintEHardM(tau1hlikelihood, tau2hlikelihood, mass);
    HHFitConstraint* balancehlikelihood = new HHFitConstraint4Vector(higgshlikelihood, true, true, false, false);
    HHFitConstraint* Likelihoodh = new HHFitConstraintLikelihood(tau1hlikelihood,tau2hlikelihood,PDF1,PDF2);
    
    //fit
    HHKinFit* singlefithlikelihood = new HHKinFit();
    singlefithlikelihood->addFitObjectE(tau1hlikelihood);
    singlefithlikelihood->addConstraint(invmhlikelihood);
    singlefithlikelihood->addConstraint(balancehlikelihood);
    singlefithlikelihood->addConstraint(Likelihoodh);
    try {
      singlefithlikelihood->fit();
      if (!((singlefithlikelihood->getConvergence()==1)||(singlefithlikelihood->getConvergence()==2))) {
	i--;
	continue;
      }
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
    h_fracresulutiontau1weightedpdfhlikelihood.Fill(comparefrac1hlikelihood,PDF1->Eval(fitfractau1hlikelihood));
    h_fracresulutiontau2weightedpdfhlikelihood.Fill(comparefrac2hlikelihood,PDF2->Eval(fitfractau2hlikelihood));
    h_energyresulution1weightedpdfhlikelihood.Fill((higgsgenerator.getTau1boosted().E()-tau1hlikelihood->getFit4Vector().E())/higgsgenerator.getTau1boosted().E(),PDF1->Eval(fitfractau1hlikelihood));
    h_energyresulution2weightedpdfhlikelihood.Fill((higgsgenerator.getTau2boosted().E()-tau2hlikelihood->getFit4Vector().E())/higgsgenerator.getTau2boosted().E(),PDF2->Eval(fitfractau2hlikelihood));
    h_Likelihoodhlikelihood.Fill(Likelihoodh->getLikelihood()*balancehlikelihood->getLikelihood());

  }
  //--------------------------------------------------------------------------------------------------------------------------------------
  //Fitting background-events:
  for(unsigned int i=0; i<numberofevents; i++){
    try{
      BackroundGenerator.generateEvent();
      HHLorentzVector test=BackroundGenerator.getTau1boosted()+BackroundGenerator.getTau2boosted();
      h_ztest.Fill(test.M());
    }
    catch(const HHEnergyRangeException& e){
      i--;
      continue;
    }
    
    //without likelihood---------------------------------------------------------------------------------------------------------------------------
    double mass = higgsgenerator.getMhiggs();
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
      if (!((singlefitb->getConvergence()==1)||(singlefitb->getConvergence()==2))) {
	i--;
	continue;
      }
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
    h_fracresulutiontau1weightedpdfb.Fill(comparefrac1b,PDF1->Eval(fitfractau1b));
    h_fracresulutiontau2weightedpdfb.Fill(comparefrac2b,PDF2->Eval(fitfractau2b));
    h_energyresulution1weightedprobb.Fill((BackroundGenerator.getTau1boosted().E()-tau1b->getFit4Vector().E())/BackroundGenerator.getTau1boosted().E(),TMath::Prob(singlefitb->getChi2(),1));
    h_energyresulution2weightedprobb.Fill((BackroundGenerator.getTau2boosted().E()-tau2b->getFit4Vector().E())/BackroundGenerator.getTau2boosted().E(),TMath::Prob(singlefitb->getChi2(),1));
    h_energyresulution1weightedpdfb.Fill((BackroundGenerator.getTau1boosted().E()-tau1b->getFit4Vector().E())/BackroundGenerator.getTau1boosted().E(),PDF1->Eval(fitfractau1b));
    h_energyresulution2weightedpdfb.Fill((BackroundGenerator.getTau2boosted().E()-tau2b->getFit4Vector().E())/BackroundGenerator.getTau2boosted().E(),PDF2->Eval(fitfractau2b));
    h_Likelihoodb.Fill(balanceb->getLikelihood()*PDF1->Eval(fitfractau1b)*PDF2->Eval(fitfractau2b));
  }

  //with likelihood-----------------------------------------------------------------------------------------------------------------------------------
  for(unsigned int i=0; i<numberofevents; i++){
    try{
      BackroundGenerator.generateEvent();
      HHLorentzVector test=BackroundGenerator.getTau1boosted()+BackroundGenerator.getTau2boosted();
      h_ztest.Fill(test.M());
    }
    catch(const HHEnergyRangeException& e){
      i--;
      continue;
    }
    double mass = higgsgenerator.getMhiggs();
    //prepare tau objects
    HHFitObjectE* tau1blikelihood = new HHFitObjectEConstM(BackroundGenerator.getTau1Vis());//  visible Tau1-component from HHTauTauEventGenerator
    HHFitObjectE* tau2blikelihood = new HHFitObjectEConstM(BackroundGenerator.getTau2Vis());//  visible Tau2-component from HHTauTauEventGenerator
    
    //prepare MET object
    HHFitObjectMET* metblikelihood = new HHFitObjectMET(TVector2(BackroundGenerator.getMETwithsigma()[0],BackroundGenerator.getMETwithsigma()[1]));//Use Met components from HHTauTauEventGenerator
    metblikelihood->setCovMatrix(BackroundGenerator.getCovarmatrix()[0][0],BackroundGenerator.getCovarmatrix()[1][1],BackroundGenerator.getCovarmatrix()[1][0]); // set Covarmatrix with Matrix inserted in HHTauTauEventGenerator
    
    //prepare composite object: Higgs
    HHFitObject* higgsblikelihood  = new HHFitObjectComposite(tau1blikelihood, tau2blikelihood, metblikelihood);
    
    tau1blikelihood->setLowerFitLimitE(tau1blikelihood->getInitial4Vector());
    tau1blikelihood->setUpperFitLimitE(mass,tau2blikelihood->getInitial4Vector());
    tau2blikelihood->setLowerFitLimitE(tau2blikelihood->getInitial4Vector());
    tau2blikelihood->setUpperFitLimitE(mass,tau1blikelihood->getInitial4Vector());
    
    //prepare constraints
    HHFitConstraint* invmblikelihood = new HHFitConstraintEHardM(tau1blikelihood, tau2blikelihood, mass);
    HHFitConstraint* balanceblikelihood = new HHFitConstraint4Vector(higgsblikelihood, true, true, false, false);
    HHFitConstraint* Likelihoodb = new HHFitConstraintLikelihood(tau1blikelihood,tau2blikelihood,PDF1,PDF2);
    
    //fit
    HHKinFit* singlefitblikelihood = new HHKinFit();
    singlefitblikelihood->addFitObjectE(tau1blikelihood);
    singlefitblikelihood->addConstraint(invmblikelihood);
    singlefitblikelihood->addConstraint(balanceblikelihood);
    singlefitblikelihood->addConstraint(Likelihoodb);
    try {
      singlefitblikelihood->fit();
      if (!((singlefitblikelihood->getConvergence()==1)||(singlefitblikelihood->getConvergence()==2))) {
	i--;
	continue;
      }
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
    h_fracresulutiontau1weightedpdfblikelihood.Fill(comparefrac1blikelihood,PDF1->Eval(fitfractau1blikelihood));
    h_fracresulutiontau2weightedpdfblikelihood.Fill(comparefrac2blikelihood,PDF2->Eval(fitfractau2blikelihood));
    h_energyresulution1weightedpdfblikelihood.Fill((BackroundGenerator.getTau1boosted().E()-tau1blikelihood->getFit4Vector().E())/BackroundGenerator.getTau1boosted().E(),PDF1->Eval(fitfractau1blikelihood));
    h_energyresulution2weightedpdfblikelihood.Fill((BackroundGenerator.getTau2boosted().E()-tau2blikelihood->getFit4Vector().E())/BackroundGenerator.getTau2boosted().E(),PDF2->Eval(fitfractau2blikelihood));
    h_Likelihoodblikelihood.Fill(Likelihoodb->getLikelihood()*balanceblikelihood->getLikelihood());
    //--------------------------------------------------------------------------------------------------------------------------------------
  }
  
  //alternative events
  //----------------------------------------------------------------------------------
  for(unsigned int i=0; i<numberofevents; i++){
    try{
      higgsgenerator2.generateEvent();
    }
    catch(const HHEnergyRangeException& e){
      i--;
      continue;
    }



    //Fitting higgsevents:

    //without likelihood---------------------------------------------------------------------------------------------------------------------------
    //KinFit:
    double mass = higgsgenerator.getMhiggs();

    //prepare tau objects
    HHFitObjectE* tau1halt = new HHFitObjectEConstM(higgsgenerator2.getTau1Vis());//  visible Tau1-component from HHTauTauEventGenerator
    HHFitObjectE* tau2halt = new HHFitObjectEConstM(higgsgenerator2.getTau2Vis());//  visible Tau2-component from HHTauTauEventGenerator

    //prepare MET object
    HHFitObjectMET* methalt = new HHFitObjectMET(TVector2(higgsgenerator2.getMETwithsigma()[0],higgsgenerator2.getMETwithsigma()[1]));//Use Met components from HHTauTauEventGenerator
    methalt->setCovMatrix(higgsgenerator2.getCovarmatrix()[0][0],higgsgenerator2.getCovarmatrix()[1][1],higgsgenerator2.getCovarmatrix()[1][0]); // set Covarmatrix with Matrix inserted in HHTauTauEventGenerator

    //prepare composite object: Higgs
    HHFitObject* higgshalt  = new HHFitObjectComposite(tau1halt, tau2halt, methalt);

    tau1halt->setLowerFitLimitE(tau1halt->getInitial4Vector());
    tau1halt->setUpperFitLimitE(mass,tau2halt->getInitial4Vector());
    tau2halt->setLowerFitLimitE(tau2halt->getInitial4Vector());
    tau2halt->setUpperFitLimitE(mass,tau1halt->getInitial4Vector());


    //prepare constraints
    HHFitConstraint* invmhalt = new HHFitConstraintEHardM(tau1halt, tau2halt, mass);
    HHFitConstraint* balancehalt = new HHFitConstraint4Vector(higgshalt, true, true, false, false);


    //fit
    HHKinFit* singlefithalt = new HHKinFit();
    singlefithalt->addFitObjectE(tau1halt);
    singlefithalt->addConstraint(invmhalt);
    singlefithalt->addConstraint(balancehalt);
    try {
      singlefithalt->fit();
      if (!((singlefithalt->getConvergence()==1)||(singlefithalt->getConvergence()==2))) {
	i--;
	continue;
      }
    }
    catch(HHEnergyRangeException const& e){
      std::cout << i << std::endl;
      tau1halt->print();
      tau2halt->print();
      methalt->print();
      higgshalt->print();
      std::cout << e.what() << std::endl;
      std::cout << higgsgenerator2.m_seed << std::endl;
      //throw(e);
      std::cout << "-----------------------------------------------" << std::endl;
      continue;
    }
    // histogramms
    h_FitFinalChi2halt.Fill(singlefithalt->getChi2());
    h_FitFinalChi2probhalt.Fill(TMath::Prob(singlefithalt->getChi2(),1));

    double  fitfractau1halt=tau1halt->getInitial4Vector().E()/tau1halt->getFit4Vector().E();
    double genfractau1halt=higgsgenerator2.getvisfrac1();
    double comparefrac1halt=(genfractau1halt-fitfractau1halt)/genfractau1halt;
    h_fracresulutiontau1halt.Fill(comparefrac1halt);
    double  fitfractau2halt=tau2halt->getInitial4Vector().E()/tau2halt->getFit4Vector().E();
    double genfractau2halt=higgsgenerator2.getvisfrac2();
    double comparefrac2halt=(genfractau2halt-fitfractau2halt)/genfractau2halt;
    h_fracresulutiontau2halt.Fill(comparefrac2halt);
    h_energyresulution1halt.Fill((higgsgenerator2.getTau1boosted().E()-tau1halt->getFit4Vector().E())/higgsgenerator2.getTau1boosted().E());
    h_energyresulution2halt.Fill((higgsgenerator2.getTau2boosted().E()-tau2halt->getFit4Vector().E())/higgsgenerator2.getTau2boosted().E());
    h_fracresulutiontau1weightedprobhalt.Fill(comparefrac1halt,TMath::Prob(singlefithalt->getChi2(),1));
    h_fracresulutiontau2weightedprobhalt.Fill(comparefrac2halt,TMath::Prob(singlefithalt->getChi2(),1));

    h_fracresulutiontau1weightedpdfhalt.Fill(comparefrac1halt,PDF1->Eval(fitfractau1halt));
    h_fracresulutiontau2weightedpdfhalt.Fill(comparefrac2halt,PDF2->Eval(fitfractau2halt));
    h_energyresulution1weightedprobhalt.Fill((higgsgenerator2.getTau1boosted().E()-tau1halt->getFit4Vector().E())/higgsgenerator2.getTau1boosted().E(),TMath::Prob(singlefithalt->getChi2(),1));
    h_energyresulution2weightedprobhalt.Fill((higgsgenerator2.getTau2boosted().E()-tau2halt->getFit4Vector().E())/higgsgenerator2.getTau2boosted().E(),TMath::Prob(singlefithalt->getChi2(),1));
    h_energyresulution1weightedpdfhalt.Fill((higgsgenerator2.getTau1boosted().E()-tau1halt->getFit4Vector().E())/higgsgenerator2.getTau1boosted().E(),PDF1->Eval(fitfractau1halt));
    h_energyresulution2weightedpdfhalt.Fill((higgsgenerator2.getTau2boosted().E()-tau2halt->getFit4Vector().E())/higgsgenerator2.getTau2boosted().E(),PDF2->Eval(fitfractau2halt));
    h_Likelihoodhalt.Fill(balancehalt->getLikelihood()*PDF1->Eval(fitfractau1halt)*PDF2->Eval(fitfractau2halt));

  }
  for(unsigned int i=0; i<numberofevents; i++){
    try{
      higgsgenerator2.generateEvent();
    }
    catch(const HHEnergyRangeException& e){
      i--;
      continue;
    }

    //with likelihood-----------------------------------------------------------------------------------------------------------------------------------
    double mass = higgsgenerator.getMhiggs();
    //prepare tau objects
    HHFitObjectE* tau1haltlikelihood = new HHFitObjectEConstM(higgsgenerator2.getTau1Vis());//  visible Tau1-component from HHTauTauEventGenerator
    HHFitObjectE* tau2haltlikelihood = new HHFitObjectEConstM(higgsgenerator2.getTau2Vis());//  visible Tau2-component from HHTauTauEventGenerator

    //prepare MET object
    HHFitObjectMET* methaltlikelihood = new HHFitObjectMET(TVector2(higgsgenerator2.getMETwithsigma()[0],higgsgenerator2.getMETwithsigma()[1]));//Use Met components from HHTauTauEventGenerator
    methaltlikelihood->setCovMatrix(higgsgenerator2.getCovarmatrix()[0][0],higgsgenerator2.getCovarmatrix()[1][1],higgsgenerator2.getCovarmatrix()[1][0]); // set Covarmatrix with Matrix inserted in HHTauTauEventGenerator

    //prepare composite object: Higgs
    HHFitObject* higgshaltlikelihood  = new HHFitObjectComposite(tau1haltlikelihood, tau2haltlikelihood, methaltlikelihood);

    tau1haltlikelihood->setLowerFitLimitE(tau1haltlikelihood->getInitial4Vector());
    tau1haltlikelihood->setUpperFitLimitE(mass,tau2haltlikelihood->getInitial4Vector());
    tau2haltlikelihood->setLowerFitLimitE(tau2haltlikelihood->getInitial4Vector());
    tau2haltlikelihood->setUpperFitLimitE(mass,tau1haltlikelihood->getInitial4Vector());

    //prepare constraints
    HHFitConstraint* invmhaltlikelihood = new HHFitConstraintEHardM(tau1haltlikelihood, tau2haltlikelihood, mass);
    HHFitConstraint* balancehaltlikelihood = new HHFitConstraint4Vector(higgshaltlikelihood, true, true, false, false);
    HHFitConstraint* Likelihoodhalt = new HHFitConstraintLikelihood(tau1haltlikelihood,tau2haltlikelihood,PDF1,PDF2);

    //fit
    HHKinFit* singlefithaltlikelihood = new HHKinFit();
    singlefithaltlikelihood->addFitObjectE(tau1haltlikelihood);
    singlefithaltlikelihood->addConstraint(invmhaltlikelihood);
    singlefithaltlikelihood->addConstraint(balancehaltlikelihood);
    singlefithaltlikelihood->addConstraint(Likelihoodhalt);
    try {
      singlefithaltlikelihood->fit();
      if (!((singlefithaltlikelihood->getConvergence()==1)||(singlefithaltlikelihood->getConvergence()==2))) {
	i--;
	continue;
      }
    }
    catch(HHEnergyRangeException const& e){
      std::cout << i << std::endl;
      tau1haltlikelihood->print();
      tau2haltlikelihood->print();
      methaltlikelihood->print();
      higgshaltlikelihood->print();
      std::cout << e.what() << std::endl;
      std::cout << higgsgenerator2.m_seed << std::endl;
      //throw(e);
      std::cout << "-----------------------------------------------" << std::endl;
      continue;
    }
    // histogramms
    h_FitFinalChi2haltlikelihood.Fill(singlefithaltlikelihood->getChi2());

    double  fitfractau1haltlikelihood=tau1haltlikelihood->getInitial4Vector().E()/tau1haltlikelihood->getFit4Vector().E();
    double genfractau1haltlikelihood=higgsgenerator2.getvisfrac1();
    double comparefrac1haltlikelihood=(genfractau1haltlikelihood-fitfractau1haltlikelihood)/genfractau1haltlikelihood;
    h_fracresulutiontau1haltlikelihood.Fill(comparefrac1haltlikelihood);
    double  fitfractau2haltlikelihood=tau2haltlikelihood->getInitial4Vector().E()/tau2haltlikelihood->getFit4Vector().E();
    double genfractau2haltlikelihood=higgsgenerator2.getvisfrac2();
    double comparefrac2haltlikelihood=(genfractau2haltlikelihood-fitfractau2haltlikelihood)/genfractau2haltlikelihood;
    h_fracresulutiontau2haltlikelihood.Fill(comparefrac2haltlikelihood);
    h_energyresulution1haltlikelihood.Fill((higgsgenerator2.getTau1boosted().E()-tau1haltlikelihood->getFit4Vector().E())/higgsgenerator2.getTau1boosted().E());
    h_energyresulution2haltlikelihood.Fill((higgsgenerator2.getTau2boosted().E()-tau2haltlikelihood->getFit4Vector().E())/higgsgenerator2.getTau2boosted().E());
    h_fracresulutiontau1weightedpdfhaltlikelihood.Fill(comparefrac1haltlikelihood,PDF1->Eval(fitfractau1haltlikelihood));
    h_fracresulutiontau2weightedpdfhaltlikelihood.Fill(comparefrac2haltlikelihood,PDF2->Eval(fitfractau2haltlikelihood));
    h_energyresulution1weightedpdfhaltlikelihood.Fill((higgsgenerator2.getTau1boosted().E()-tau1haltlikelihood->getFit4Vector().E())/higgsgenerator2.getTau1boosted().E(),PDF1->Eval(fitfractau1haltlikelihood));
    h_energyresulution2weightedpdfhaltlikelihood.Fill((higgsgenerator2.getTau2boosted().E()-tau2haltlikelihood->getFit4Vector().E())/higgsgenerator2.getTau2boosted().E(),PDF2->Eval(fitfractau2haltlikelihood));
    h_Likelihoodhaltlikelihood.Fill(Likelihoodhalt->getLikelihood()*balancehaltlikelihood->getLikelihood());
    //--------------------------------------------------------------------------------------------------------------------------------------
    //Fitting background-events:
  }
  for(unsigned int i=0; i<numberofevents; i++){
    try{
      BackroundGenerator2.generateEvent();
      HHLorentzVector test=BackroundGenerator2.getTau1boosted()+BackroundGenerator2.getTau2boosted();
      h_ztest.Fill(test.M());
    }
    catch(const HHEnergyRangeException& e){
      i--;
      continue;
    }
    //without likelihood---------------------------------------------------------------------------------------------------------------------------
    //KinFit:
    double mass = higgsgenerator.getMhiggs();


    //prepare tau objects
    HHFitObjectE* tau1balt = new HHFitObjectEConstM(BackroundGenerator2.getTau1Vis());//  visible Tau1-component from HHTauTauEventGenerator
    HHFitObjectE* tau2balt = new HHFitObjectEConstM(BackroundGenerator2.getTau2Vis());//  visible Tau2-component from HHTauTauEventGenerator

    //prepare MET object
    HHFitObjectMET* metbalt = new HHFitObjectMET(TVector2(BackroundGenerator2.getMETwithsigma()[0],BackroundGenerator2.getMETwithsigma()[1]));//Use Met components from HHTauTauEventGenerator
    metbalt->setCovMatrix(BackroundGenerator2.getCovarmatrix()[0][0],BackroundGenerator2.getCovarmatrix()[1][1],BackroundGenerator2.getCovarmatrix()[1][0]); // set Covarmatrix with Matrix inserted in HHTauTauEventGenerator

    //prepare composite object: Higgs
    HHFitObject* higgsbalt  = new HHFitObjectComposite(tau1balt, tau2balt, metbalt);

    tau1balt->setLowerFitLimitE(tau1balt->getInitial4Vector());
    tau1balt->setUpperFitLimitE(mass,tau2balt->getInitial4Vector());
    tau2balt->setLowerFitLimitE(tau2balt->getInitial4Vector());
    tau2balt->setUpperFitLimitE(mass,tau1balt->getInitial4Vector());


    //prepare constraints
    HHFitConstraint* invmbalt = new HHFitConstraintEHardM(tau1balt, tau2balt, mass);
    HHFitConstraint* balancebalt = new HHFitConstraint4Vector(higgsbalt, true, true, false, false);


    //fit
    HHKinFit* singlefitbalt = new HHKinFit();
    singlefitbalt->addFitObjectE(tau1balt);
    singlefitbalt->addConstraint(invmbalt);
    singlefitbalt->addConstraint(balancebalt);
    try {
      singlefitbalt->fit();
      if (!((singlefitbalt->getConvergence()==1)||(singlefitbalt->getConvergence()==2))) {
	i--;
	continue;
      }
    }
    catch(HHEnergyRangeException const& e){
      std::cout << i << std::endl;
      tau1balt->print();
      tau2balt->print();
      metbalt->print();
      higgsbalt->print();
      std::cout << e.what() << std::endl;
      std::cout << BackroundGenerator2.m_seed << std::endl;
      //throw(e);
      std::cout << "-----------------------------------------------" << std::endl;
      continue;
    }
    // histogramms
    h_FitFinalChi2balt.Fill(singlefitbalt->getChi2());
    h_FitFinalChi2probbalt.Fill(TMath::Prob(singlefitbalt->getChi2(),1));

    double  fitfractau1balt=tau1balt->getInitial4Vector().E()/tau1balt->getFit4Vector().E();
    double genfractau1balt=BackroundGenerator2.getvisfrac1();
    double comparefrac1balt=(genfractau1balt-fitfractau1balt)/genfractau1balt;
    h_fracresulutiontau1balt.Fill(comparefrac1balt);
    double  fitfractau2balt=tau2balt->getInitial4Vector().E()/tau2balt->getFit4Vector().E();
    double genfractau2balt=BackroundGenerator2.getvisfrac2();
    double comparefrac2balt=(genfractau2balt-fitfractau2balt)/genfractau2balt;
    h_fracresulutiontau2balt.Fill(comparefrac2balt);
    h_energyresulution1balt.Fill((BackroundGenerator2.getTau1boosted().E()-tau1balt->getFit4Vector().E())/BackroundGenerator2.getTau1boosted().E());
    h_energyresulution2balt.Fill((BackroundGenerator2.getTau2boosted().E()-tau2balt->getFit4Vector().E())/BackroundGenerator2.getTau2boosted().E());
    h_fracresulutiontau1weightedprobbalt.Fill(comparefrac1balt,TMath::Prob(singlefitbalt->getChi2(),1));
    h_fracresulutiontau2weightedprobbalt.Fill(comparefrac2balt,TMath::Prob(singlefitbalt->getChi2(),1));
    h_fracresulutiontau1weightedpdfbalt.Fill(comparefrac1balt,PDF1->Eval(fitfractau1balt));
    h_fracresulutiontau2weightedpdfbalt.Fill(comparefrac2balt,PDF2->Eval(fitfractau2balt));
    h_energyresulution1weightedprobbalt.Fill((BackroundGenerator2.getTau1boosted().E()-tau1balt->getFit4Vector().E())/BackroundGenerator2.getTau1boosted().E(),TMath::Prob(singlefitbalt->getChi2(),1));
    h_energyresulution2weightedprobbalt.Fill((BackroundGenerator2.getTau2boosted().E()-tau2balt->getFit4Vector().E())/BackroundGenerator2.getTau2boosted().E(),TMath::Prob(singlefitbalt->getChi2(),1));
    h_energyresulution1weightedpdfbalt.Fill((BackroundGenerator2.getTau1boosted().E()-tau1balt->getFit4Vector().E())/BackroundGenerator2.getTau1boosted().E(),PDF1->Eval(fitfractau1balt));
    h_energyresulution2weightedpdfbalt.Fill((BackroundGenerator2.getTau2boosted().E()-tau2balt->getFit4Vector().E())/BackroundGenerator2.getTau2boosted().E(),PDF2->Eval(fitfractau2balt));
    h_Likelihoodbalt.Fill(balancebalt->getLikelihood()*PDF1->Eval(fitfractau1balt)*PDF2->Eval(fitfractau2balt));

  }
  for(unsigned int i=0; i<numberofevents; i++){
    try{
      BackroundGenerator2.generateEvent();
      HHLorentzVector test=BackroundGenerator2.getTau1boosted()+BackroundGenerator2.getTau2boosted();
      h_ztest.Fill(test.M());
    }
    catch(const HHEnergyRangeException& e){
      i--;
      continue;
    }

    //with likelihood-----------------------------------------------------------------------------------------------------------------------------------
    double mass = higgsgenerator.getMhiggs();
    //prepare tau objects
    HHFitObjectE* tau1baltlikelihood = new HHFitObjectEConstM(BackroundGenerator2.getTau1Vis());//  visible Tau1-component from HHTauTauEventGenerator
    HHFitObjectE* tau2baltlikelihood = new HHFitObjectEConstM(BackroundGenerator2.getTau2Vis());//  visible Tau2-component from HHTauTauEventGenerator

    //prepare MET object
    HHFitObjectMET* metbaltlikelihood = new HHFitObjectMET(TVector2(BackroundGenerator2.getMETwithsigma()[0],BackroundGenerator2.getMETwithsigma()[1]));//Use Met components from HHTauTauEventGenerator
    metbaltlikelihood->setCovMatrix(BackroundGenerator2.getCovarmatrix()[0][0],BackroundGenerator2.getCovarmatrix()[1][1],BackroundGenerator2.getCovarmatrix()[1][0]); // set Covarmatrix with Matrix inserted in HHTauTauEventGenerator

    //prepare composite object: Higgs
    HHFitObject* higgsbaltlikelihood  = new HHFitObjectComposite(tau1baltlikelihood, tau2baltlikelihood, metbaltlikelihood);

    tau1baltlikelihood->setLowerFitLimitE(tau1baltlikelihood->getInitial4Vector());
    tau1baltlikelihood->setUpperFitLimitE(mass,tau2baltlikelihood->getInitial4Vector());
    tau2baltlikelihood->setLowerFitLimitE(tau2baltlikelihood->getInitial4Vector());
    tau2baltlikelihood->setUpperFitLimitE(mass,tau1baltlikelihood->getInitial4Vector());

    //prepare constraints
    HHFitConstraint* invmbaltlikelihood = new HHFitConstraintEHardM(tau1baltlikelihood, tau2baltlikelihood, mass);
    HHFitConstraint* balancebaltlikelihood = new HHFitConstraint4Vector(higgsbaltlikelihood, true, true, false, false);
    HHFitConstraint* Likelihoodbalt = new HHFitConstraintLikelihood(tau1baltlikelihood,tau2baltlikelihood,PDF1,PDF2);

    //fit
    HHKinFit* singlefitbaltlikelihood = new HHKinFit();
    singlefitbaltlikelihood->addFitObjectE(tau1baltlikelihood);
    singlefitbaltlikelihood->addConstraint(invmbaltlikelihood);
    singlefitbaltlikelihood->addConstraint(balancebaltlikelihood);
    singlefitbaltlikelihood->addConstraint(Likelihoodbalt);
    try {
      singlefitbaltlikelihood->fit();
      if (!((singlefitbaltlikelihood->getConvergence()==1)||(singlefitbaltlikelihood->getConvergence()==2))) {
	i--;
	continue;
      }
    }
    catch(HHEnergyRangeException const& e){
      std::cout << i << std::endl;
      tau1baltlikelihood->print();
      tau2baltlikelihood->print();
      metbaltlikelihood->print();
      higgsbaltlikelihood->print();
      std::cout << e.what() << std::endl;
      std::cout << BackroundGenerator2.m_seed << std::endl;
      //throw(e);
      std::cout << "-----------------------------------------------" << std::endl;
      continue;
    }
    // histogramms
    h_FitFinalChi2baltlikelihood.Fill(singlefitbaltlikelihood->getChi2());

    double  fitfractau1baltlikelihood=tau1baltlikelihood->getInitial4Vector().E()/tau1baltlikelihood->getFit4Vector().E();
    double genfractau1baltlikelihood=BackroundGenerator2.getvisfrac1();
    double comparefrac1baltlikelihood=(genfractau1baltlikelihood-fitfractau1baltlikelihood)/genfractau1baltlikelihood;
    h_fracresulutiontau1baltlikelihood.Fill(comparefrac1baltlikelihood);
    double  fitfractau2baltlikelihood=tau2baltlikelihood->getInitial4Vector().E()/tau2baltlikelihood->getFit4Vector().E();
    double genfractau2baltlikelihood=BackroundGenerator2.getvisfrac2();
    double comparefrac2baltlikelihood=(genfractau2baltlikelihood-fitfractau2baltlikelihood)/genfractau2baltlikelihood;
    h_fracresulutiontau2baltlikelihood.Fill(comparefrac2baltlikelihood);
    h_energyresulution1baltlikelihood.Fill((BackroundGenerator2.getTau1boosted().E()-tau1baltlikelihood->getFit4Vector().E())/BackroundGenerator2.getTau1boosted().E());
    h_energyresulution2baltlikelihood.Fill((BackroundGenerator2.getTau2boosted().E()-tau2baltlikelihood->getFit4Vector().E())/BackroundGenerator2.getTau2boosted().E());
    h_fracresulutiontau1weightedpdfbaltlikelihood.Fill(comparefrac1baltlikelihood,PDF1->Eval(fitfractau1baltlikelihood));
    h_fracresulutiontau2weightedpdfbaltlikelihood.Fill(comparefrac2baltlikelihood,PDF2->Eval(fitfractau2baltlikelihood));
    h_energyresulution1weightedpdfbaltlikelihood.Fill((BackroundGenerator2.getTau1boosted().E()-tau1baltlikelihood->getFit4Vector().E())/BackroundGenerator2.getTau1boosted().E(),PDF1->Eval(fitfractau1baltlikelihood));
    h_energyresulution2weightedpdfbaltlikelihood.Fill((BackroundGenerator2.getTau2boosted().E()-tau2baltlikelihood->getFit4Vector().E())/BackroundGenerator2.getTau2boosted().E(),PDF2->Eval(fitfractau2baltlikelihood));
    h_Likelihoodbaltlikelihood.Fill(Likelihoodbalt->getLikelihood()*balancebaltlikelihood->getLikelihood());
    //--------------------------------------------------------------------------------------------------------------------------------------


  }
  //-------------------------------------------------------------------------------------------------------------------------------------------
  //alt2ernative events2
  //----------------------------------------------------------------------------------
  for(unsigned int i=0; i<numberofevents; i++){
    try{
      higgsgenerator3.generateEvent();
    }
    catch(const HHEnergyRangeException& e){
      i--;
      continue;
    }



    //Fitting higgsevents:

    //without likelihood---------------------------------------------------------------------------------------------------------------------------
    //KinFit:
    double mass = higgsgenerator.getMhiggs();

    //prepare tau objects
    HHFitObjectE* tau1halt2 = new HHFitObjectEConstM(higgsgenerator3.getTau1Vis());//  visible Tau1-component from HHTauTauEventGenerator
    HHFitObjectE* tau2halt2 = new HHFitObjectEConstM(higgsgenerator3.getTau2Vis());//  visible Tau2-component from HHTauTauEventGenerator

    //prepare MET object
    HHFitObjectMET* methalt2 = new HHFitObjectMET(TVector2(higgsgenerator3.getMETwithsigma()[0],higgsgenerator3.getMETwithsigma()[1]));//Use Met components from HHTauTauEventGenerator
    methalt2->setCovMatrix(higgsgenerator3.getCovarmatrix()[0][0],higgsgenerator3.getCovarmatrix()[1][1],higgsgenerator3.getCovarmatrix()[1][0]); // set Covarmatrix with Matrix inserted in HHTauTauEventGenerator

    //prepare composite object: Higgs
    HHFitObject* higgshalt2  = new HHFitObjectComposite(tau1halt2, tau2halt2, methalt2);

    tau1halt2->setLowerFitLimitE(tau1halt2->getInitial4Vector());
    tau1halt2->setUpperFitLimitE(mass,tau2halt2->getInitial4Vector());
    tau2halt2->setLowerFitLimitE(tau2halt2->getInitial4Vector());
    tau2halt2->setUpperFitLimitE(mass,tau1halt2->getInitial4Vector());


    //prepare constraints
    HHFitConstraint* invmhalt2 = new HHFitConstraintEHardM(tau1halt2, tau2halt2, mass);
    HHFitConstraint* balancehalt2 = new HHFitConstraint4Vector(higgshalt2, true, true, false, false);


    //fit
    HHKinFit* singlefithalt2 = new HHKinFit();
    singlefithalt2->addFitObjectE(tau1halt2);
    singlefithalt2->addConstraint(invmhalt2);
    singlefithalt2->addConstraint(balancehalt2);
    try {
      singlefithalt2->fit();
      if (!((singlefithalt2->getConvergence()==1)||(singlefithalt2->getConvergence()==2))) {
	i--;
	continue;
      }
    }
    catch(HHEnergyRangeException const& e){
      std::cout << i << std::endl;
      tau1halt2->print();
      tau2halt2->print();
      methalt2->print();
      higgshalt2->print();
      std::cout << e.what() << std::endl;
      std::cout << higgsgenerator3.m_seed << std::endl;
      //throw(e);
      std::cout << "-----------------------------------------------" << std::endl;
      continue;
    }
    // histogramms
    h_FitFinalChi2halt2.Fill(singlefithalt2->getChi2());
    h_FitFinalChi2probhalt2.Fill(TMath::Prob(singlefithalt2->getChi2(),1));

    double  fitfractau1halt2=tau1halt2->getInitial4Vector().E()/tau1halt2->getFit4Vector().E();
    double genfractau1halt2=higgsgenerator3.getvisfrac1();
    double comparefrac1halt2=(genfractau1halt2-fitfractau1halt2)/genfractau1halt2;
    h_fracresulutiontau1halt2.Fill(comparefrac1halt2);
    double  fitfractau2halt2=tau2halt2->getInitial4Vector().E()/tau2halt2->getFit4Vector().E();
    double genfractau2halt2=higgsgenerator3.getvisfrac2();
    double comparefrac2halt2=(genfractau2halt2-fitfractau2halt2)/genfractau2halt2;
    h_fracresulutiontau2halt2.Fill(comparefrac2halt2);
    h_energyresulution1halt2.Fill((higgsgenerator3.getTau1boosted().E()-tau1halt2->getFit4Vector().E())/higgsgenerator3.getTau1boosted().E());
    h_energyresulution2halt2.Fill((higgsgenerator3.getTau2boosted().E()-tau2halt2->getFit4Vector().E())/higgsgenerator3.getTau2boosted().E());
    h_fracresulutiontau1weightedprobhalt2.Fill(comparefrac1halt2,TMath::Prob(singlefithalt2->getChi2(),1));
    h_fracresulutiontau2weightedprobhalt2.Fill(comparefrac2halt2,TMath::Prob(singlefithalt2->getChi2(),1));

    h_fracresulutiontau1weightedpdfhalt2.Fill(comparefrac1halt2,PDF1->Eval(fitfractau1halt2));
    h_fracresulutiontau2weightedpdfhalt2.Fill(comparefrac2halt2,PDF2->Eval(fitfractau2halt2));
    h_energyresulution1weightedprobhalt2.Fill((higgsgenerator3.getTau1boosted().E()-tau1halt2->getFit4Vector().E())/higgsgenerator3.getTau1boosted().E(),TMath::Prob(singlefithalt2->getChi2(),1));
    h_energyresulution2weightedprobhalt2.Fill((higgsgenerator3.getTau2boosted().E()-tau2halt2->getFit4Vector().E())/higgsgenerator3.getTau2boosted().E(),TMath::Prob(singlefithalt2->getChi2(),1));
    h_energyresulution1weightedpdfhalt2.Fill((higgsgenerator3.getTau1boosted().E()-tau1halt2->getFit4Vector().E())/higgsgenerator3.getTau1boosted().E(),PDF1->Eval(fitfractau1halt2));
    h_energyresulution2weightedpdfhalt2.Fill((higgsgenerator3.getTau2boosted().E()-tau2halt2->getFit4Vector().E())/higgsgenerator3.getTau2boosted().E(),PDF2->Eval(fitfractau2halt2));
    h_Likelihoodhalt2.Fill(balancehalt2->getLikelihood()*PDF1->Eval(fitfractau1halt2)*PDF2->Eval(fitfractau2halt2));

  }
  for(unsigned int i=0; i<numberofevents; i++){
    try{
      higgsgenerator3.generateEvent();
    }
    catch(const HHEnergyRangeException& e){
      i--;
      continue;
    }

    //with likelihood-----------------------------------------------------------------------------------------------------------------------------------
    double mass = higgsgenerator.getMhiggs();
    //prepare tau objects
    HHFitObjectE* tau1halt2likelihood = new HHFitObjectEConstM(higgsgenerator3.getTau1Vis());//  visible Tau1-component from HHTauTauEventGenerator
    HHFitObjectE* tau2halt2likelihood = new HHFitObjectEConstM(higgsgenerator3.getTau2Vis());//  visible Tau2-component from HHTauTauEventGenerator

    //prepare MET object
    HHFitObjectMET* methalt2likelihood = new HHFitObjectMET(TVector2(higgsgenerator3.getMETwithsigma()[0],higgsgenerator3.getMETwithsigma()[1]));//Use Met components from HHTauTauEventGenerator
    methalt2likelihood->setCovMatrix(higgsgenerator3.getCovarmatrix()[0][0],higgsgenerator3.getCovarmatrix()[1][1],higgsgenerator3.getCovarmatrix()[1][0]); // set Covarmatrix with Matrix inserted in HHTauTauEventGenerator

    //prepare composite object: Higgs
    HHFitObject* higgshalt2likelihood  = new HHFitObjectComposite(tau1halt2likelihood, tau2halt2likelihood, methalt2likelihood);

    tau1halt2likelihood->setLowerFitLimitE(tau1halt2likelihood->getInitial4Vector());
    tau1halt2likelihood->setUpperFitLimitE(mass,tau2halt2likelihood->getInitial4Vector());
    tau2halt2likelihood->setLowerFitLimitE(tau2halt2likelihood->getInitial4Vector());
    tau2halt2likelihood->setUpperFitLimitE(mass,tau1halt2likelihood->getInitial4Vector());

    //prepare constraints
    HHFitConstraint* invmhalt2likelihood = new HHFitConstraintEHardM(tau1halt2likelihood, tau2halt2likelihood, mass);
    HHFitConstraint* balancehalt2likelihood = new HHFitConstraint4Vector(higgshalt2likelihood, true, true, false, false);
    HHFitConstraint* Likelihoodhalt2 = new HHFitConstraintLikelihood(tau1halt2likelihood,tau2halt2likelihood,PDF1,PDF2);

    //fit
    HHKinFit* singlefithalt2likelihood = new HHKinFit();
    singlefithalt2likelihood->addFitObjectE(tau1halt2likelihood);
    singlefithalt2likelihood->addConstraint(invmhalt2likelihood);
    singlefithalt2likelihood->addConstraint(balancehalt2likelihood);
    singlefithalt2likelihood->addConstraint(Likelihoodhalt2);
    try {
      singlefithalt2likelihood->fit();
      if (!((singlefithalt2likelihood->getConvergence()==1)||(singlefithalt2likelihood->getConvergence()==2))) {
	i--;
	continue;
      }
    }
    catch(HHEnergyRangeException const& e){
      std::cout << i << std::endl;
      tau1halt2likelihood->print();
      tau2halt2likelihood->print();
      methalt2likelihood->print();
      higgshalt2likelihood->print();
      std::cout << e.what() << std::endl;
      std::cout << higgsgenerator3.m_seed << std::endl;
      //throw(e);
      std::cout << "-----------------------------------------------" << std::endl;
      continue;
    }
    // histogramms
    h_FitFinalChi2halt2likelihood.Fill(singlefithalt2likelihood->getChi2());

    double  fitfractau1halt2likelihood=tau1halt2likelihood->getInitial4Vector().E()/tau1halt2likelihood->getFit4Vector().E();
    double genfractau1halt2likelihood=higgsgenerator3.getvisfrac1();
    double comparefrac1halt2likelihood=(genfractau1halt2likelihood-fitfractau1halt2likelihood)/genfractau1halt2likelihood;
    h_fracresulutiontau1halt2likelihood.Fill(comparefrac1halt2likelihood);
    double  fitfractau2halt2likelihood=tau2halt2likelihood->getInitial4Vector().E()/tau2halt2likelihood->getFit4Vector().E();
    double genfractau2halt2likelihood=higgsgenerator3.getvisfrac2();
    double comparefrac2halt2likelihood=(genfractau2halt2likelihood-fitfractau2halt2likelihood)/genfractau2halt2likelihood;
    h_fracresulutiontau2halt2likelihood.Fill(comparefrac2halt2likelihood);
    h_energyresulution1halt2likelihood.Fill((higgsgenerator3.getTau1boosted().E()-tau1halt2likelihood->getFit4Vector().E())/higgsgenerator3.getTau1boosted().E());
    h_energyresulution2halt2likelihood.Fill((higgsgenerator3.getTau2boosted().E()-tau2halt2likelihood->getFit4Vector().E())/higgsgenerator3.getTau2boosted().E());
    h_fracresulutiontau1weightedpdfhalt2likelihood.Fill(comparefrac1halt2likelihood,PDF1->Eval(fitfractau1halt2likelihood));
    h_fracresulutiontau2weightedpdfhalt2likelihood.Fill(comparefrac2halt2likelihood,PDF2->Eval(fitfractau2halt2likelihood));
    h_energyresulution1weightedpdfhalt2likelihood.Fill((higgsgenerator3.getTau1boosted().E()-tau1halt2likelihood->getFit4Vector().E())/higgsgenerator3.getTau1boosted().E(),PDF1->Eval(fitfractau1halt2likelihood));
    h_energyresulution2weightedpdfhalt2likelihood.Fill((higgsgenerator3.getTau2boosted().E()-tau2halt2likelihood->getFit4Vector().E())/higgsgenerator3.getTau2boosted().E(),PDF2->Eval(fitfractau2halt2likelihood));
    h_Likelihoodhalt2likelihood.Fill(Likelihoodhalt2->getLikelihood()*balancehalt2likelihood->getLikelihood());
    //--------------------------------------------------------------------------------------------------------------------------------------
    //Fitting background-events:
  }
  for(unsigned int i=0; i<numberofevents; i++){
    try{
      BackroundGenerator3.generateEvent();
      HHLorentzVector test=BackroundGenerator3.getTau1boosted()+BackroundGenerator3.getTau2boosted();
      h_ztest.Fill(test.M());
    }
    catch(const HHEnergyRangeException& e){
      i--;
      continue;
    }
    //without likelihood---------------------------------------------------------------------------------------------------------------------------
    //KinFit:
    double mass = higgsgenerator.getMhiggs();


    //prepare tau objects
    HHFitObjectE* tau1balt2 = new HHFitObjectEConstM(BackroundGenerator3.getTau1Vis());//  visible Tau1-component from HHTauTauEventGenerator
    HHFitObjectE* tau2balt2 = new HHFitObjectEConstM(BackroundGenerator3.getTau2Vis());//  visible Tau2-component from HHTauTauEventGenerator

    //prepare MET object
    HHFitObjectMET* metbalt2 = new HHFitObjectMET(TVector2(BackroundGenerator3.getMETwithsigma()[0],BackroundGenerator3.getMETwithsigma()[1]));//Use Met components from HHTauTauEventGenerator
    metbalt2->setCovMatrix(BackroundGenerator3.getCovarmatrix()[0][0],BackroundGenerator3.getCovarmatrix()[1][1],BackroundGenerator3.getCovarmatrix()[1][0]); // set Covarmatrix with Matrix inserted in HHTauTauEventGenerator

    //prepare composite object: Higgs
    HHFitObject* higgsbalt2  = new HHFitObjectComposite(tau1balt2, tau2balt2, metbalt2);

    tau1balt2->setLowerFitLimitE(tau1balt2->getInitial4Vector());
    tau1balt2->setUpperFitLimitE(mass,tau2balt2->getInitial4Vector());
    tau2balt2->setLowerFitLimitE(tau2balt2->getInitial4Vector());
    tau2balt2->setUpperFitLimitE(mass,tau1balt2->getInitial4Vector());


    //prepare constraints
    HHFitConstraint* invmbalt2 = new HHFitConstraintEHardM(tau1balt2, tau2balt2, mass);
    HHFitConstraint* balancebalt2 = new HHFitConstraint4Vector(higgsbalt2, true, true, false, false);


    //fit
    HHKinFit* singlefitbalt2 = new HHKinFit();
    singlefitbalt2->addFitObjectE(tau1balt2);
    singlefitbalt2->addConstraint(invmbalt2);
    singlefitbalt2->addConstraint(balancebalt2);
    try {
      singlefitbalt2->fit();
      if (!((singlefitbalt2->getConvergence()==1)||(singlefitbalt2->getConvergence()==2))) {
	i--;
	continue;
      }
    }
    catch(HHEnergyRangeException const& e){
      std::cout << i << std::endl;
      tau1balt2->print();
      tau2balt2->print();
      metbalt2->print();
      higgsbalt2->print();
      std::cout << e.what() << std::endl;
      std::cout << BackroundGenerator3.m_seed << std::endl;
      //throw(e);
      std::cout << "-----------------------------------------------" << std::endl;
      continue;
    }
    // histogramms
    h_FitFinalChi2balt2.Fill(singlefitbalt2->getChi2());
    h_FitFinalChi2probbalt2.Fill(TMath::Prob(singlefitbalt2->getChi2(),1));

    double  fitfractau1balt2=tau1balt2->getInitial4Vector().E()/tau1balt2->getFit4Vector().E();
    double genfractau1balt2=BackroundGenerator3.getvisfrac1();
    double comparefrac1balt2=(genfractau1balt2-fitfractau1balt2)/genfractau1balt2;
    h_fracresulutiontau1balt2.Fill(comparefrac1balt2);
    double  fitfractau2balt2=tau2balt2->getInitial4Vector().E()/tau2balt2->getFit4Vector().E();
    double genfractau2balt2=BackroundGenerator3.getvisfrac2();
    double comparefrac2balt2=(genfractau2balt2-fitfractau2balt2)/genfractau2balt2;
    h_fracresulutiontau2balt2.Fill(comparefrac2balt2);
    h_energyresulution1balt2.Fill((BackroundGenerator3.getTau1boosted().E()-tau1balt2->getFit4Vector().E())/BackroundGenerator3.getTau1boosted().E());
    h_energyresulution2balt2.Fill((BackroundGenerator3.getTau2boosted().E()-tau2balt2->getFit4Vector().E())/BackroundGenerator3.getTau2boosted().E());
    h_fracresulutiontau1weightedprobbalt2.Fill(comparefrac1balt2,TMath::Prob(singlefitbalt2->getChi2(),1));
    h_fracresulutiontau2weightedprobbalt2.Fill(comparefrac2balt2,TMath::Prob(singlefitbalt2->getChi2(),1));
    h_fracresulutiontau1weightedpdfbalt2.Fill(comparefrac1balt2,PDF1->Eval(fitfractau1balt2));
    h_fracresulutiontau2weightedpdfbalt2.Fill(comparefrac2balt2,PDF2->Eval(fitfractau2balt2));
    h_energyresulution1weightedprobbalt2.Fill((BackroundGenerator3.getTau1boosted().E()-tau1balt2->getFit4Vector().E())/BackroundGenerator3.getTau1boosted().E(),TMath::Prob(singlefitbalt2->getChi2(),1));
    h_energyresulution2weightedprobbalt2.Fill((BackroundGenerator3.getTau2boosted().E()-tau2balt2->getFit4Vector().E())/BackroundGenerator3.getTau2boosted().E(),TMath::Prob(singlefitbalt2->getChi2(),1));
    h_energyresulution1weightedpdfbalt2.Fill((BackroundGenerator3.getTau1boosted().E()-tau1balt2->getFit4Vector().E())/BackroundGenerator3.getTau1boosted().E(),PDF1->Eval(fitfractau1balt2));
    h_energyresulution2weightedpdfbalt2.Fill((BackroundGenerator3.getTau2boosted().E()-tau2balt2->getFit4Vector().E())/BackroundGenerator3.getTau2boosted().E(),PDF2->Eval(fitfractau2balt2));
    h_Likelihoodbalt2.Fill(balancebalt2->getLikelihood()*PDF1->Eval(fitfractau1balt2)*PDF2->Eval(fitfractau2balt2));

  }
  for(unsigned int i=0; i<numberofevents; i++){
    try{
      BackroundGenerator.generateEvent();
      HHLorentzVector test=BackroundGenerator3.getTau1boosted()+BackroundGenerator3.getTau2boosted();
      h_ztest.Fill(test.M());
    }
    catch(const HHEnergyRangeException& e){
      i--;
      continue;
    }

    //with likelihood-----------------------------------------------------------------------------------------------------------------------------------
    double mass = higgsgenerator.getMhiggs();
    //prepare tau objects
    HHFitObjectE* tau1balt2likelihood = new HHFitObjectEConstM(BackroundGenerator3.getTau1Vis());//  visible Tau1-component from HHTauTauEventGenerator
    HHFitObjectE* tau2balt2likelihood = new HHFitObjectEConstM(BackroundGenerator3.getTau2Vis());//  visible Tau2-component from HHTauTauEventGenerator

    //prepare MET object
    HHFitObjectMET* metbalt2likelihood = new HHFitObjectMET(TVector2(BackroundGenerator3.getMETwithsigma()[0],BackroundGenerator3.getMETwithsigma()[1]));//Use Met components from HHTauTauEventGenerator
    metbalt2likelihood->setCovMatrix(BackroundGenerator3.getCovarmatrix()[0][0],BackroundGenerator3.getCovarmatrix()[1][1],BackroundGenerator3.getCovarmatrix()[1][0]); // set Covarmatrix with Matrix inserted in HHTauTauEventGenerator

    //prepare composite object: Higgs
    HHFitObject* higgsbalt2likelihood  = new HHFitObjectComposite(tau1balt2likelihood, tau2balt2likelihood, metbalt2likelihood);

    tau1balt2likelihood->setLowerFitLimitE(tau1balt2likelihood->getInitial4Vector());
    tau1balt2likelihood->setUpperFitLimitE(mass,tau2balt2likelihood->getInitial4Vector());
    tau2balt2likelihood->setLowerFitLimitE(tau2balt2likelihood->getInitial4Vector());
    tau2balt2likelihood->setUpperFitLimitE(mass,tau1balt2likelihood->getInitial4Vector());

    //prepare constraints
    HHFitConstraint* invmbalt2likelihood = new HHFitConstraintEHardM(tau1balt2likelihood, tau2balt2likelihood, mass);
    HHFitConstraint* balancebalt2likelihood = new HHFitConstraint4Vector(higgsbalt2likelihood, true, true, false, false);
    HHFitConstraint* Likelihoodbalt2 = new HHFitConstraintLikelihood(tau1balt2likelihood,tau2balt2likelihood,PDF1,PDF2);

    //fit
    HHKinFit* singlefitbalt2likelihood = new HHKinFit();
    singlefitbalt2likelihood->addFitObjectE(tau1balt2likelihood);
    singlefitbalt2likelihood->addConstraint(invmbalt2likelihood);
    singlefitbalt2likelihood->addConstraint(balancebalt2likelihood);
    singlefitbalt2likelihood->addConstraint(Likelihoodbalt2);
    try {
      singlefitbalt2likelihood->fit();
      if (!((singlefitbalt2likelihood->getConvergence()==1)||(singlefitbalt2likelihood->getConvergence()==2))) {
	i--;
	continue;
      }
    }
    catch(HHEnergyRangeException const& e){
      std::cout << i << std::endl;
      tau1balt2likelihood->print();
      tau2balt2likelihood->print();
      metbalt2likelihood->print();
      higgsbalt2likelihood->print();
      std::cout << e.what() << std::endl;
      std::cout << BackroundGenerator3.m_seed << std::endl;
      //throw(e);
      std::cout << "-----------------------------------------------" << std::endl;
      continue;
    }
    // histogramms
    h_FitFinalChi2balt2likelihood.Fill(singlefitbalt2likelihood->getChi2());

    double  fitfractau1balt2likelihood=tau1balt2likelihood->getInitial4Vector().E()/tau1balt2likelihood->getFit4Vector().E();
    double genfractau1balt2likelihood=BackroundGenerator3.getvisfrac1();
    double comparefrac1balt2likelihood=(genfractau1balt2likelihood-fitfractau1balt2likelihood)/genfractau1balt2likelihood;
    h_fracresulutiontau1balt2likelihood.Fill(comparefrac1balt2likelihood);
    double  fitfractau2balt2likelihood=tau2balt2likelihood->getInitial4Vector().E()/tau2balt2likelihood->getFit4Vector().E();
    double genfractau2balt2likelihood=BackroundGenerator3.getvisfrac2();
    double comparefrac2balt2likelihood=(genfractau2balt2likelihood-fitfractau2balt2likelihood)/genfractau2balt2likelihood;
    h_fracresulutiontau2balt2likelihood.Fill(comparefrac2balt2likelihood);
    h_energyresulution1balt2likelihood.Fill((BackroundGenerator3.getTau1boosted().E()-tau1balt2likelihood->getFit4Vector().E())/BackroundGenerator3.getTau1boosted().E());
    h_energyresulution2balt2likelihood.Fill((BackroundGenerator3.getTau2boosted().E()-tau2balt2likelihood->getFit4Vector().E())/BackroundGenerator3.getTau2boosted().E());
    h_fracresulutiontau1weightedpdfbalt2likelihood.Fill(comparefrac1balt2likelihood,PDF1->Eval(fitfractau1balt2likelihood));
    h_fracresulutiontau2weightedpdfbalt2likelihood.Fill(comparefrac2balt2likelihood,PDF2->Eval(fitfractau2balt2likelihood));
    h_energyresulution1weightedpdfbalt2likelihood.Fill((BackroundGenerator3.getTau1boosted().E()-tau1balt2likelihood->getFit4Vector().E())/BackroundGenerator3.getTau1boosted().E(),PDF1->Eval(fitfractau1balt2likelihood));
    h_energyresulution2weightedpdfbalt2likelihood.Fill((BackroundGenerator3.getTau2boosted().E()-tau2balt2likelihood->getFit4Vector().E())/BackroundGenerator3.getTau2boosted().E(),PDF2->Eval(fitfractau2balt2likelihood));
    h_Likelihoodbalt2likelihood.Fill(Likelihoodbalt2->getLikelihood()*balancebalt2likelihood->getLikelihood());
    //--------------------------------------------------------------------------------------------------------------------------------------


  }
  //--------------------------------------------------------------------------------------------------------------------------------------------
  //alt3ernative events
  //----------------------------------------------------------------------------------
  for(unsigned int i=0; i<numberofevents; i++){
    try{
      higgsgenerator4.generateEvent();
    }
    catch(const HHEnergyRangeException& e){
      i--;
      continue;
    }



    //Fitting higgsevents:

    //without likelihood---------------------------------------------------------------------------------------------------------------------------
    //KinFit:
    double mass = higgsgenerator.getMhiggs();

    //prepare tau objects
    HHFitObjectE* tau1halt3 = new HHFitObjectEConstM(higgsgenerator4.getTau1Vis());//  visible Tau1-component from HHTauTauEventGenerator
    HHFitObjectE* tau2halt3 = new HHFitObjectEConstM(higgsgenerator4.getTau2Vis());//  visible Tau2-component from HHTauTauEventGenerator

    //prepare MET object
    HHFitObjectMET* methalt3 = new HHFitObjectMET(TVector2(higgsgenerator4.getMETwithsigma()[0],higgsgenerator4.getMETwithsigma()[1]));//Use Met components from HHTauTauEventGenerator
    methalt3->setCovMatrix(higgsgenerator4.getCovarmatrix()[0][0],higgsgenerator4.getCovarmatrix()[1][1],higgsgenerator4.getCovarmatrix()[1][0]); // set Covarmatrix with Matrix inserted in HHTauTauEventGenerator

    //prepare composite object: Higgs
    HHFitObject* higgshalt3  = new HHFitObjectComposite(tau1halt3, tau2halt3, methalt3);

    tau1halt3->setLowerFitLimitE(tau1halt3->getInitial4Vector());
    tau1halt3->setUpperFitLimitE(mass,tau2halt3->getInitial4Vector());
    tau2halt3->setLowerFitLimitE(tau2halt3->getInitial4Vector());
    tau2halt3->setUpperFitLimitE(mass,tau1halt3->getInitial4Vector());


    //prepare constraints
    HHFitConstraint* invmhalt3 = new HHFitConstraintEHardM(tau1halt3, tau2halt3, mass);
    HHFitConstraint* balancehalt3 = new HHFitConstraint4Vector(higgshalt3, true, true, false, false);


    //fit
    HHKinFit* singlefithalt3 = new HHKinFit();
    singlefithalt3->addFitObjectE(tau1halt3);
    singlefithalt3->addConstraint(invmhalt3);
    singlefithalt3->addConstraint(balancehalt3);
    try {
      singlefithalt3->fit();
      if (!((singlefithalt3->getConvergence()==1)||(singlefithalt3->getConvergence()==2))) {
	i--;
	continue;
      }
    }
    catch(HHEnergyRangeException const& e){
      std::cout << i << std::endl;
      tau1halt3->print();
      tau2halt3->print();
      methalt3->print();
      higgshalt3->print();
      std::cout << e.what() << std::endl;
      std::cout << higgsgenerator4.m_seed << std::endl;
      //throw(e);
      std::cout << "-----------------------------------------------" << std::endl;
      continue;
    }
    // histogramms
    h_FitFinalChi2halt3.Fill(singlefithalt3->getChi2());
    h_FitFinalChi2probhalt3.Fill(TMath::Prob(singlefithalt3->getChi2(),1));

    double  fitfractau1halt3=tau1halt3->getInitial4Vector().E()/tau1halt3->getFit4Vector().E();
    double genfractau1halt3=higgsgenerator4.getvisfrac1();
    double comparefrac1halt3=(genfractau1halt3-fitfractau1halt3)/genfractau1halt3;
    h_fracresulutiontau1halt3.Fill(comparefrac1halt3);
    double  fitfractau2halt3=tau2halt3->getInitial4Vector().E()/tau2halt3->getFit4Vector().E();
    double genfractau2halt3=higgsgenerator4.getvisfrac2();
    double comparefrac2halt3=(genfractau2halt3-fitfractau2halt3)/genfractau2halt3;
    h_fracresulutiontau2halt3.Fill(comparefrac2halt3);
    h_energyresulution1halt3.Fill((higgsgenerator4.getTau1boosted().E()-tau1halt3->getFit4Vector().E())/higgsgenerator4.getTau1boosted().E());
    h_energyresulution2halt3.Fill((higgsgenerator4.getTau2boosted().E()-tau2halt3->getFit4Vector().E())/higgsgenerator4.getTau2boosted().E());
    h_fracresulutiontau1weightedprobhalt3.Fill(comparefrac1halt3,TMath::Prob(singlefithalt3->getChi2(),1));
    h_fracresulutiontau2weightedprobhalt3.Fill(comparefrac2halt3,TMath::Prob(singlefithalt3->getChi2(),1));

    h_fracresulutiontau1weightedpdfhalt3.Fill(comparefrac1halt3,PDF1->Eval(fitfractau1halt3));
    h_fracresulutiontau2weightedpdfhalt3.Fill(comparefrac2halt3,PDF2->Eval(fitfractau2halt3));
    h_energyresulution1weightedprobhalt3.Fill((higgsgenerator4.getTau1boosted().E()-tau1halt3->getFit4Vector().E())/higgsgenerator4.getTau1boosted().E(),TMath::Prob(singlefithalt3->getChi2(),1));
    h_energyresulution2weightedprobhalt3.Fill((higgsgenerator4.getTau2boosted().E()-tau2halt3->getFit4Vector().E())/higgsgenerator4.getTau2boosted().E(),TMath::Prob(singlefithalt3->getChi2(),1));
    h_energyresulution1weightedpdfhalt3.Fill((higgsgenerator4.getTau1boosted().E()-tau1halt3->getFit4Vector().E())/higgsgenerator4.getTau1boosted().E(),PDF1->Eval(fitfractau1halt3));
    h_energyresulution2weightedpdfhalt3.Fill((higgsgenerator4.getTau2boosted().E()-tau2halt3->getFit4Vector().E())/higgsgenerator4.getTau2boosted().E(),PDF2->Eval(fitfractau2halt3));
    h_Likelihoodhalt3.Fill(balancehalt3->getLikelihood()*PDF1->Eval(fitfractau1halt3)*PDF2->Eval(fitfractau2halt3));

  }
  for(unsigned int i=0; i<numberofevents; i++){
    try{
      higgsgenerator4.generateEvent();
    }
    catch(const HHEnergyRangeException& e){
      i--;
      continue;
    }

    //with likelihood-----------------------------------------------------------------------------------------------------------------------------------
    double mass = higgsgenerator.getMhiggs();
    //prepare tau objects
    HHFitObjectE* tau1halt3likelihood = new HHFitObjectEConstM(higgsgenerator4.getTau1Vis());//  visible Tau1-component from HHTauTauEventGenerator
    HHFitObjectE* tau2halt3likelihood = new HHFitObjectEConstM(higgsgenerator4.getTau2Vis());//  visible Tau2-component from HHTauTauEventGenerator

    //prepare MET object
    HHFitObjectMET* methalt3likelihood = new HHFitObjectMET(TVector2(higgsgenerator4.getMETwithsigma()[0],higgsgenerator4.getMETwithsigma()[1]));//Use Met components from HHTauTauEventGenerator
    methalt3likelihood->setCovMatrix(higgsgenerator4.getCovarmatrix()[0][0],higgsgenerator4.getCovarmatrix()[1][1],higgsgenerator4.getCovarmatrix()[1][0]); // set Covarmatrix with Matrix inserted in HHTauTauEventGenerator

    //prepare composite object: Higgs
    HHFitObject* higgshalt3likelihood  = new HHFitObjectComposite(tau1halt3likelihood, tau2halt3likelihood, methalt3likelihood);

    tau1halt3likelihood->setLowerFitLimitE(tau1halt3likelihood->getInitial4Vector());
    tau1halt3likelihood->setUpperFitLimitE(mass,tau2halt3likelihood->getInitial4Vector());
    tau2halt3likelihood->setLowerFitLimitE(tau2halt3likelihood->getInitial4Vector());
    tau2halt3likelihood->setUpperFitLimitE(mass,tau1halt3likelihood->getInitial4Vector());

    //prepare constraints
    HHFitConstraint* invmhalt3likelihood = new HHFitConstraintEHardM(tau1halt3likelihood, tau2halt3likelihood, mass);
    HHFitConstraint* balancehalt3likelihood = new HHFitConstraint4Vector(higgshalt3likelihood, true, true, false, false);
    HHFitConstraint* Likelihoodhalt3 = new HHFitConstraintLikelihood(tau1halt3likelihood,tau2halt3likelihood,PDF1,PDF2);

    //fit
    HHKinFit* singlefithalt3likelihood = new HHKinFit();
    singlefithalt3likelihood->addFitObjectE(tau1halt3likelihood);
    singlefithalt3likelihood->addConstraint(invmhalt3likelihood);
    singlefithalt3likelihood->addConstraint(balancehalt3likelihood);
    singlefithalt3likelihood->addConstraint(Likelihoodhalt3);
    try {
      singlefithalt3likelihood->fit();
      if (!((singlefithalt3likelihood->getConvergence()==1)||(singlefithalt3likelihood->getConvergence()==2))) {
	i--;
	continue;
      }
    }
    catch(HHEnergyRangeException const& e){
      std::cout << i << std::endl;
      tau1halt3likelihood->print();
      tau2halt3likelihood->print();
      methalt3likelihood->print();
      higgshalt3likelihood->print();
      std::cout << e.what() << std::endl;
      std::cout << higgsgenerator4.m_seed << std::endl;
      //throw(e);
      std::cout << "-----------------------------------------------" << std::endl;
      continue;
    }
    // histogramms
    h_FitFinalChi2halt3likelihood.Fill(singlefithalt3likelihood->getChi2());

    double  fitfractau1halt3likelihood=tau1halt3likelihood->getInitial4Vector().E()/tau1halt3likelihood->getFit4Vector().E();
    double genfractau1halt3likelihood=higgsgenerator4.getvisfrac1();
    double comparefrac1halt3likelihood=(genfractau1halt3likelihood-fitfractau1halt3likelihood)/genfractau1halt3likelihood;
    h_fracresulutiontau1halt3likelihood.Fill(comparefrac1halt3likelihood);
    double  fitfractau2halt3likelihood=tau2halt3likelihood->getInitial4Vector().E()/tau2halt3likelihood->getFit4Vector().E();
    double genfractau2halt3likelihood=higgsgenerator4.getvisfrac2();
    double comparefrac2halt3likelihood=(genfractau2halt3likelihood-fitfractau2halt3likelihood)/genfractau2halt3likelihood;
    h_fracresulutiontau2halt3likelihood.Fill(comparefrac2halt3likelihood);
    h_energyresulution1halt3likelihood.Fill((higgsgenerator4.getTau1boosted().E()-tau1halt3likelihood->getFit4Vector().E())/higgsgenerator4.getTau1boosted().E());
    h_energyresulution2halt3likelihood.Fill((higgsgenerator4.getTau2boosted().E()-tau2halt3likelihood->getFit4Vector().E())/higgsgenerator4.getTau2boosted().E());
    h_fracresulutiontau1weightedpdfhalt3likelihood.Fill(comparefrac1halt3likelihood,PDF1->Eval(fitfractau1halt3likelihood));
    h_fracresulutiontau2weightedpdfhalt3likelihood.Fill(comparefrac2halt3likelihood,PDF2->Eval(fitfractau2halt3likelihood));
    h_energyresulution1weightedpdfhalt3likelihood.Fill((higgsgenerator4.getTau1boosted().E()-tau1halt3likelihood->getFit4Vector().E())/higgsgenerator4.getTau1boosted().E(),PDF1->Eval(fitfractau1halt3likelihood));
    h_energyresulution2weightedpdfhalt3likelihood.Fill((higgsgenerator4.getTau2boosted().E()-tau2halt3likelihood->getFit4Vector().E())/higgsgenerator4.getTau2boosted().E(),PDF2->Eval(fitfractau2halt3likelihood));
    h_Likelihoodhalt3likelihood.Fill(Likelihoodhalt3->getLikelihood()*balancehalt3likelihood->getLikelihood());
    //--------------------------------------------------------------------------------------------------------------------------------------
    //Fitting background-events:
  }
  for(unsigned int i=0; i<numberofevents; i++){
    try{
      BackroundGenerator4.generateEvent();
      HHLorentzVector test=BackroundGenerator4.getTau1boosted()+BackroundGenerator4.getTau2boosted();
      h_ztest.Fill(test.M());
    }
    catch(const HHEnergyRangeException& e){
      i--;
      continue;
    }
    //without likelihood---------------------------------------------------------------------------------------------------------------------------
    //KinFit:
    double mass = higgsgenerator.getMhiggs();


    //prepare tau objects
    HHFitObjectE* tau1balt3 = new HHFitObjectEConstM(BackroundGenerator4.getTau1Vis());//  visible Tau1-component from HHTauTauEventGenerator
    HHFitObjectE* tau2balt3 = new HHFitObjectEConstM(BackroundGenerator4.getTau2Vis());//  visible Tau2-component from HHTauTauEventGenerator

    //prepare MET object
    HHFitObjectMET* metbalt3 = new HHFitObjectMET(TVector2(BackroundGenerator4.getMETwithsigma()[0],BackroundGenerator4.getMETwithsigma()[1]));//Use Met components from HHTauTauEventGenerator
    metbalt3->setCovMatrix(BackroundGenerator4.getCovarmatrix()[0][0],BackroundGenerator4.getCovarmatrix()[1][1],BackroundGenerator4.getCovarmatrix()[1][0]); // set Covarmatrix with Matrix inserted in HHTauTauEventGenerator

    //prepare composite object: Higgs
    HHFitObject* higgsbalt3  = new HHFitObjectComposite(tau1balt3, tau2balt3, metbalt3);

    tau1balt3->setLowerFitLimitE(tau1balt3->getInitial4Vector());
    tau1balt3->setUpperFitLimitE(mass,tau2balt3->getInitial4Vector());
    tau2balt3->setLowerFitLimitE(tau2balt3->getInitial4Vector());
    tau2balt3->setUpperFitLimitE(mass,tau1balt3->getInitial4Vector());


    //prepare constraints
    HHFitConstraint* invmbalt3 = new HHFitConstraintEHardM(tau1balt3, tau2balt3, mass);
    HHFitConstraint* balancebalt3 = new HHFitConstraint4Vector(higgsbalt3, true, true, false, false);


    //fit
    HHKinFit* singlefitbalt3 = new HHKinFit();
    singlefitbalt3->addFitObjectE(tau1balt3);
    singlefitbalt3->addConstraint(invmbalt3);
    singlefitbalt3->addConstraint(balancebalt3);
    try {
      singlefitbalt3->fit();
      if (!((singlefitbalt3->getConvergence()==1)||(singlefitbalt3->getConvergence()==2))) {
	i--;
	continue;
      }
    }
    catch(HHEnergyRangeException const& e){
      std::cout << i << std::endl;
      tau1balt3->print();
      tau2balt3->print();
      metbalt3->print();
      higgsbalt3->print();
      std::cout << e.what() << std::endl;
      std::cout << BackroundGenerator4.m_seed << std::endl;
      //throw(e);
      std::cout << "-----------------------------------------------" << std::endl;
      continue;
    }
    // histogramms
    h_FitFinalChi2balt3.Fill(singlefitbalt3->getChi2());
    h_FitFinalChi2probbalt3.Fill(TMath::Prob(singlefitbalt3->getChi2(),1));

    double  fitfractau1balt3=tau1balt3->getInitial4Vector().E()/tau1balt3->getFit4Vector().E();
    double genfractau1balt3=BackroundGenerator4.getvisfrac1();
    double comparefrac1balt3=(genfractau1balt3-fitfractau1balt3)/genfractau1balt3;
    h_fracresulutiontau1balt3.Fill(comparefrac1balt3);
    double  fitfractau2balt3=tau2balt3->getInitial4Vector().E()/tau2balt3->getFit4Vector().E();
    double genfractau2balt3=BackroundGenerator4.getvisfrac2();
    double comparefrac2balt3=(genfractau2balt3-fitfractau2balt3)/genfractau2balt3;
    h_fracresulutiontau2balt3.Fill(comparefrac2balt3);
    h_energyresulution1balt3.Fill((BackroundGenerator4.getTau1boosted().E()-tau1balt3->getFit4Vector().E())/BackroundGenerator4.getTau1boosted().E());
    h_energyresulution2balt3.Fill((BackroundGenerator4.getTau2boosted().E()-tau2balt3->getFit4Vector().E())/BackroundGenerator4.getTau2boosted().E());
    h_fracresulutiontau1weightedprobbalt3.Fill(comparefrac1balt3,TMath::Prob(singlefitbalt3->getChi2(),1));
    h_fracresulutiontau2weightedprobbalt3.Fill(comparefrac2balt3,TMath::Prob(singlefitbalt3->getChi2(),1));
    h_fracresulutiontau1weightedpdfbalt3.Fill(comparefrac1balt3,PDF1->Eval(fitfractau1balt3));
    h_fracresulutiontau2weightedpdfbalt3.Fill(comparefrac2balt3,PDF2->Eval(fitfractau2balt3));
    h_energyresulution1weightedprobbalt3.Fill((BackroundGenerator4.getTau1boosted().E()-tau1balt3->getFit4Vector().E())/BackroundGenerator4.getTau1boosted().E(),TMath::Prob(singlefitbalt3->getChi2(),1));
    h_energyresulution2weightedprobbalt3.Fill((BackroundGenerator4.getTau2boosted().E()-tau2balt3->getFit4Vector().E())/BackroundGenerator4.getTau2boosted().E(),TMath::Prob(singlefitbalt3->getChi2(),1));
    h_energyresulution1weightedpdfbalt3.Fill((BackroundGenerator4.getTau1boosted().E()-tau1balt3->getFit4Vector().E())/BackroundGenerator4.getTau1boosted().E(),PDF1->Eval(fitfractau1balt3));
    h_energyresulution2weightedpdfbalt3.Fill((BackroundGenerator4.getTau2boosted().E()-tau2balt3->getFit4Vector().E())/BackroundGenerator4.getTau2boosted().E(),PDF2->Eval(fitfractau2balt3));
    h_Likelihoodbalt3.Fill(balancebalt3->getLikelihood()*PDF1->Eval(fitfractau1balt3)*PDF2->Eval(fitfractau2balt3));

  }
  for(unsigned int i=0; i<numberofevents; i++){
    try{
      BackroundGenerator4.generateEvent();
      HHLorentzVector test=BackroundGenerator4.getTau1boosted()+BackroundGenerator4.getTau2boosted();
      h_ztest.Fill(test.M());
    }
    catch(const HHEnergyRangeException& e){
      i--;
      continue;
    }

    //with likelihood-----------------------------------------------------------------------------------------------------------------------------------
    double mass = higgsgenerator.getMhiggs();
    //prepare tau objects
    HHFitObjectE* tau1balt3likelihood = new HHFitObjectEConstM(BackroundGenerator4.getTau1Vis());//  visible Tau1-component from HHTauTauEventGenerator
    HHFitObjectE* tau2balt3likelihood = new HHFitObjectEConstM(BackroundGenerator4.getTau2Vis());//  visible Tau2-component from HHTauTauEventGenerator

    //prepare MET object
    HHFitObjectMET* metbalt3likelihood = new HHFitObjectMET(TVector2(BackroundGenerator4.getMETwithsigma()[0],BackroundGenerator4.getMETwithsigma()[1]));//Use Met components from HHTauTauEventGenerator
    metbalt3likelihood->setCovMatrix(BackroundGenerator4.getCovarmatrix()[0][0],BackroundGenerator4.getCovarmatrix()[1][1],BackroundGenerator4.getCovarmatrix()[1][0]); // set Covarmatrix with Matrix inserted in HHTauTauEventGenerator

    //prepare composite object: Higgs
    HHFitObject* higgsbalt3likelihood  = new HHFitObjectComposite(tau1balt3likelihood, tau2balt3likelihood, metbalt3likelihood);

    tau1balt3likelihood->setLowerFitLimitE(tau1balt3likelihood->getInitial4Vector());
    tau1balt3likelihood->setUpperFitLimitE(mass,tau2balt3likelihood->getInitial4Vector());
    tau2balt3likelihood->setLowerFitLimitE(tau2balt3likelihood->getInitial4Vector());
    tau2balt3likelihood->setUpperFitLimitE(mass,tau1balt3likelihood->getInitial4Vector());

    //prepare constraints
    HHFitConstraint* invmbalt3likelihood = new HHFitConstraintEHardM(tau1balt3likelihood, tau2balt3likelihood, mass);
    HHFitConstraint* balancebalt3likelihood = new HHFitConstraint4Vector(higgsbalt3likelihood, true, true, false, false);
    HHFitConstraint* Likelihoodbalt3 = new HHFitConstraintLikelihood(tau1balt3likelihood,tau2balt3likelihood,PDF1,PDF2);

    //fit
    HHKinFit* singlefitbalt3likelihood = new HHKinFit();
    singlefitbalt3likelihood->addFitObjectE(tau1balt3likelihood);
    singlefitbalt3likelihood->addConstraint(invmbalt3likelihood);
    singlefitbalt3likelihood->addConstraint(balancebalt3likelihood);
    singlefitbalt3likelihood->addConstraint(Likelihoodbalt3);
    try {
      singlefitbalt3likelihood->fit();
      if (!((singlefitbalt3likelihood->getConvergence()==1)||(singlefitbalt3likelihood->getConvergence()==2))) {
	i--;
	continue;
      }
    }
    catch(HHEnergyRangeException const& e){
      std::cout << i << std::endl;
      tau1balt3likelihood->print();
      tau2balt3likelihood->print();
      metbalt3likelihood->print();
      higgsbalt3likelihood->print();
      std::cout << e.what() << std::endl;
      std::cout << BackroundGenerator4.m_seed << std::endl;
      //throw(e);
      std::cout << "-----------------------------------------------" << std::endl;
      continue;
    }
    // histogramms
    h_FitFinalChi2balt3likelihood.Fill(singlefitbalt3likelihood->getChi2());

    double  fitfractau1balt3likelihood=tau1balt3likelihood->getInitial4Vector().E()/tau1balt3likelihood->getFit4Vector().E();
    double genfractau1balt3likelihood=BackroundGenerator4.getvisfrac1();
    double comparefrac1balt3likelihood=(genfractau1balt3likelihood-fitfractau1balt3likelihood)/genfractau1balt3likelihood;
    h_fracresulutiontau1balt3likelihood.Fill(comparefrac1balt3likelihood);
    double  fitfractau2balt3likelihood=tau2balt3likelihood->getInitial4Vector().E()/tau2balt3likelihood->getFit4Vector().E();
    double genfractau2balt3likelihood=BackroundGenerator4.getvisfrac2();
    double comparefrac2balt3likelihood=(genfractau2balt3likelihood-fitfractau2balt3likelihood)/genfractau2balt3likelihood;
    h_fracresulutiontau2balt3likelihood.Fill(comparefrac2balt3likelihood);
    h_energyresulution1balt3likelihood.Fill((BackroundGenerator4.getTau1boosted().E()-tau1balt3likelihood->getFit4Vector().E())/BackroundGenerator4.getTau1boosted().E());
    h_energyresulution2balt3likelihood.Fill((BackroundGenerator4.getTau2boosted().E()-tau2balt3likelihood->getFit4Vector().E())/BackroundGenerator4.getTau2boosted().E());
    h_fracresulutiontau1weightedpdfbalt3likelihood.Fill(comparefrac1balt3likelihood,PDF1->Eval(fitfractau1balt3likelihood));
    h_fracresulutiontau2weightedpdfbalt3likelihood.Fill(comparefrac2balt3likelihood,PDF2->Eval(fitfractau2balt3likelihood));
    h_energyresulution1weightedpdfbalt3likelihood.Fill((BackroundGenerator4.getTau1boosted().E()-tau1balt3likelihood->getFit4Vector().E())/BackroundGenerator4.getTau1boosted().E(),PDF1->Eval(fitfractau1balt3likelihood));
    h_energyresulution2weightedpdfbalt3likelihood.Fill((BackroundGenerator4.getTau2boosted().E()-tau2balt3likelihood->getFit4Vector().E())/BackroundGenerator4.getTau2boosted().E(),PDF2->Eval(fitfractau2balt3likelihood));
    h_Likelihoodbalt3likelihood.Fill(Likelihoodbalt3->getLikelihood()*balancebalt3likelihood->getLikelihood());
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
  h_Likelihoodh.Write();
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
  h_Likelihoodhlikelihood.Write();
  
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
  h_Likelihoodb.Write();
  //Histogramms for background-fitwithlikelihood
  
  h_FitFinalChi2blikelihood.Write();
  h_fracresulutiontau1blikelihood.Write();
  h_fracresulutiontau2blikelihood.Write();
  h_fracresulutiontau1weightedpdfblikelihood.Write();
  h_fracresulutiontau2weightedpdfblikelihood.Write();
  h_energyresulution2blikelihood.Write();
  h_energyresulution1weightedpdfblikelihood.Write();
  h_energyresulution2weightedpdfblikelihood.Write();
  h_Likelihoodblikelihood.Write();
  //---------------------------------------------------------------------------------------------------------------------
  h_FitFinalChi2halt.Write();
  h_FitFinalChi2probhalt.Write();
  h_fracresulutiontau1halt.Write();
  h_fracresulutiontau2halt.Write();
  h_fracresulutiontau1weightedprobhalt.Write();
  h_fracresulutiontau2weightedprobhalt.Write();
  h_fracresulutiontau1weightedpdfhalt.Write();
  h_fracresulutiontau2weightedpdfhalt.Write();
  h_energyresulution1halt.Write();
  h_energyresulution2halt.Write();
  h_energyresulution1weightedprobhalt.Write();
  h_energyresulution2weightedprobhalt.Write();
  h_energyresulution1weightedpdfhalt.Write();
  h_energyresulution2weightedpdfhalt.Write();
  h_Likelihoodhalt.Write();
  //Histogramms for h-fitwithlikelihood

  h_FitFinalChi2haltlikelihood.Write();
  h_fracresulutiontau1haltlikelihood.Write();
  h_fracresulutiontau2haltlikelihood.Write();
  h_fracresulutiontau1weightedpdfhaltlikelihood.Write();
  h_fracresulutiontau2weightedpdfhaltlikelihood.Write();
  h_energyresulution1haltlikelihood.Write();
  h_energyresulution2haltlikelihood.Write();
  h_energyresulution1weightedpdfhaltlikelihood.Write();
  h_energyresulution2weightedpdfhaltlikelihood.Write();
  h_Likelihoodhaltlikelihood.Write();

  //Histogramms for background-fit
  h_FitFinalChi2balt.Write();
  h_FitFinalChi2probbalt.Write();
  h_fracresulutiontau1balt.Write();
  h_fracresulutiontau2balt.Write();
  h_fracresulutiontau1weightedprobbalt.Write();
  h_fracresulutiontau2weightedprobbalt.Write();
  h_fracresulutiontau1weightedpdfbalt.Write();
  h_fracresulutiontau2weightedpdfbalt.Write();
  h_energyresulution1balt.Write();
  h_energyresulution2balt.Write();
  h_energyresulution1weightedprobbalt.Write();
  h_energyresulution2weightedprobbalt.Write();
  h_energyresulution1weightedpdfbalt.Write();
  h_energyresulution2weightedpdfbalt.Write();
  h_Likelihoodbalt.Write();
  //Histogramms for background-fitwithlikelihood

  h_FitFinalChi2baltlikelihood.Write();
  h_fracresulutiontau1baltlikelihood.Write();
  h_fracresulutiontau2baltlikelihood.Write();
  h_fracresulutiontau1weightedpdfbaltlikelihood.Write();
  h_fracresulutiontau2weightedpdfbaltlikelihood.Write();
  h_energyresulution2baltlikelihood.Write();
  h_energyresulution1weightedpdfbaltlikelihood.Write();
  h_energyresulution2weightedpdfbaltlikelihood.Write();
  h_Likelihoodbaltlikelihood.Write();
  //---------------------------------------------------------------------------------------------------------------------
  h_FitFinalChi2halt2.Write();
  h_FitFinalChi2probhalt2.Write();
  h_fracresulutiontau1halt2.Write();
  h_fracresulutiontau2halt2.Write();
  h_fracresulutiontau1weightedprobhalt2.Write();
  h_fracresulutiontau2weightedprobhalt2.Write();
  h_fracresulutiontau1weightedpdfhalt2.Write();
  h_fracresulutiontau2weightedpdfhalt2.Write();
  h_energyresulution1halt2.Write();
  h_energyresulution2halt2.Write();
  h_energyresulution1weightedprobhalt2.Write();
  h_energyresulution2weightedprobhalt2.Write();
  h_energyresulution1weightedpdfhalt2.Write();
  h_energyresulution2weightedpdfhalt2.Write();
  h_Likelihoodhalt2.Write();
  //Histogramms for h-fitwithlikelihood

  h_FitFinalChi2halt2likelihood.Write();
  h_fracresulutiontau1halt2likelihood.Write();
  h_fracresulutiontau2halt2likelihood.Write();
  h_fracresulutiontau1weightedpdfhalt2likelihood.Write();
  h_fracresulutiontau2weightedpdfhalt2likelihood.Write();
  h_energyresulution1halt2likelihood.Write();
  h_energyresulution2halt2likelihood.Write();
  h_energyresulution1weightedpdfhalt2likelihood.Write();
  h_energyresulution2weightedpdfhalt2likelihood.Write();
  h_Likelihoodhalt2likelihood.Write();

  //Histogramms for background-fit
  h_FitFinalChi2balt2.Write();
  h_FitFinalChi2probbalt2.Write();
  h_fracresulutiontau1balt2.Write();
  h_fracresulutiontau2balt2.Write();
  h_fracresulutiontau1weightedprobbalt2.Write();
  h_fracresulutiontau2weightedprobbalt2.Write();
  h_fracresulutiontau1weightedpdfbalt2.Write();
  h_fracresulutiontau2weightedpdfbalt2.Write();
  h_energyresulution1balt2.Write();
  h_energyresulution2balt2.Write();
  h_energyresulution1weightedprobbalt2.Write();
  h_energyresulution2weightedprobbalt2.Write();
  h_energyresulution1weightedpdfbalt2.Write();
  h_energyresulution2weightedpdfbalt2.Write();
  h_Likelihoodbalt2.Write();
  //Histogramms for background-fitwithlikelihood

  h_FitFinalChi2balt2likelihood.Write();
  h_fracresulutiontau1balt2likelihood.Write();
  h_fracresulutiontau2balt2likelihood.Write();
  h_fracresulutiontau1weightedpdfbalt2likelihood.Write();
  h_fracresulutiontau2weightedpdfbalt2likelihood.Write();
  h_energyresulution2balt2likelihood.Write();
  h_energyresulution1weightedpdfbalt2likelihood.Write();
  h_energyresulution2weightedpdfbalt2likelihood.Write();
  h_Likelihoodbalt2likelihood.Write();
  //---------------------------------------------------------------------------------------------------------------------
  h_FitFinalChi2halt3.Write();
  h_FitFinalChi2probhalt3.Write();
  h_fracresulutiontau1halt3.Write();
  h_fracresulutiontau2halt3.Write();
  h_fracresulutiontau1weightedprobhalt3.Write();
  h_fracresulutiontau2weightedprobhalt3.Write();
  h_fracresulutiontau1weightedpdfhalt3.Write();
  h_fracresulutiontau2weightedpdfhalt3.Write();
  h_energyresulution1halt3.Write();
  h_energyresulution2halt3.Write();
  h_energyresulution1weightedprobhalt3.Write();
  h_energyresulution2weightedprobhalt3.Write();
  h_energyresulution1weightedpdfhalt3.Write();
  h_energyresulution2weightedpdfhalt3.Write();
  h_Likelihoodhalt3.Write();
  //Histogramms for h-fitwithlikelihood

  h_FitFinalChi2halt3likelihood.Write();
  h_fracresulutiontau1halt3likelihood.Write();
  h_fracresulutiontau2halt3likelihood.Write();
  h_fracresulutiontau1weightedpdfhalt3likelihood.Write();
  h_fracresulutiontau2weightedpdfhalt3likelihood.Write();
  h_energyresulution1halt3likelihood.Write();
  h_energyresulution2halt3likelihood.Write();
  h_energyresulution1weightedpdfhalt3likelihood.Write();
  h_energyresulution2weightedpdfhalt3likelihood.Write();
  h_Likelihoodhalt3likelihood.Write();

  //Histogramms for background-fit
  h_FitFinalChi2balt3.Write();
  h_FitFinalChi2probbalt3.Write();
  h_fracresulutiontau1balt3.Write();
  h_fracresulutiontau2balt3.Write();
  h_fracresulutiontau1weightedprobbalt3.Write();
  h_fracresulutiontau2weightedprobbalt3.Write();
  h_fracresulutiontau1weightedpdfbalt3.Write();
  h_fracresulutiontau2weightedpdfbalt3.Write();
  h_energyresulution1balt3.Write();
  h_energyresulution2balt3.Write();
  h_energyresulution1weightedprobbalt3.Write();
  h_energyresulution2weightedprobbalt3.Write();
  h_energyresulution1weightedpdfbalt3.Write();
  h_energyresulution2weightedpdfbalt3.Write();
  h_Likelihoodbalt3.Write();
  //Histogramms for background-fitwithlikelihood

  h_FitFinalChi2balt3likelihood.Write();
  h_fracresulutiontau1balt3likelihood.Write();
  h_fracresulutiontau2balt3likelihood.Write();
  h_fracresulutiontau1weightedpdfbalt3likelihood.Write();
  h_fracresulutiontau2weightedpdfbalt3likelihood.Write();
  h_energyresulution2balt3likelihood.Write();
  h_energyresulution1weightedpdfbalt3likelihood.Write();
  h_energyresulution2weightedpdfbalt3likelihood.Write();
  h_Likelihoodbalt3likelihood.Write();

  h_ztest.Write();

  backgroundtest.Close();
  
  return(0);
}
