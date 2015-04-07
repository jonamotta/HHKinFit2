#include "HHKinFit.h"
#include "HHFitConstraint4Vector.h"
#include "HHFitConstraintEHardM.h"
#include "HHFitConstraint.h"
#include "HHFitConstraintLikelihood.h"
#include "HHFitObjectEConstM.h"
#include "HHFitObjectE.h"
#include "HHFitObject.h"
#include "HHFitObjectComposite.h"
#include "PSMath.h"
#include <iomanip>
#include "TAxis.h"


HHKinFit::HHKinFit()
: m_fitobjects(std::vector<HHFitObjectE*>()),
  m_constraints(std::vector<HHFitConstraint*>()){
}

///TODO: make it more general!
void 
HHKinFit::fit(){

	  //  ----------  for PSfit -----
	  const Int_t np = 1;
	  Double_t a[np];
	  Double_t astart[np];
	  Double_t alimit[np][2];
	  Double_t aprec[np];
	  Double_t daN[np];
	  Double_t h[np];
	  Double_t chi2iter[1], aMemory[np][5], g[np], H[np * np], Hinv[np * np];
	  Bool_t noNewtonShifts = false;

	  Int_t iter = 0;             //  number of iterations
	  Int_t method = 1;           //  initial fit method, see PSfit()
	  Int_t mode = 1;             //  mode =1 for start of a new fit by PSfit()
	//   Int_t icallNewton = -1;     //  init start of Newton Method
	//   Int_t iloop = 0;            // counter for falls to fit function




	  // fill initial tau fit parameters
	  astart[0] = m_fitobjects[0]->getE();    // energy of first tau
	  aprec[0]  = 0.1;                        //0.1                 // precision for fit
	  // fill initial step width
	  h[0] = 0.1*m_fitobjects[0]->getLowerFitLimitE();   //
	  daN[0] = 1.0;   //0.0                 // initial search direction



	  alimit[0][0] = m_fitobjects[0]->getLowerFitLimitE();              // tau: minimum is visible tau1 energy
	  alimit[0][1] = m_fitobjects[0]->getUpperFitLimitE();              //      maximum as computed above

	  // tau: check initial values against fit range
	  if (astart[0] - h[0] < alimit[0][0]) {
	    astart[0] = alimit[0][0] + h[0];
	  }
	  else if (astart[0] + h[0] > alimit[0][1]) {
	    astart[0] = alimit[0][1] - h[0];
	  }

	  for (Int_t ip = 0; ip < np; ip++) {
	    a[ip] = astart[ip];
	  }
	  for (Int_t ip = 0; ip < np; ip++) {
	    aMemory[ip][0] = -999.0;
	    aMemory[ip][1] = -995.0;
	    aMemory[ip][2] = -990.0;
	    aMemory[ip][3] = -985.0;
	    aMemory[ip][3] = -980.0;
	  }

	  Double_t chi2(99999);
	  Int_t convergence(0);
	  Int_t printlevel(0);
	  Int_t m_maxloops(10000);
	  for (Int_t iloop = 0; iloop < m_maxloops * 10 && iter < m_maxloops; iloop++) { // FIT loop
	    m_fitobjects[0]->changeEandSave(a[0]);
	    chi2=this->getChi2();
//	    std::cout << iloop << " a[0]: " << a[0] << " chi2: " << std::fixed << std::setprecision(8) << chi2 << std::endl;
//	    m_fitobjects[0]->print();



	    if (convergence != 0) break;
	    convergence = PSMath::PSfit(iloop, iter, method, mode, noNewtonShifts, printlevel,
	                                  np, a, astart, alimit, aprec,
	                                  daN, h, aMemory, chi2, chi2iter, g, H,
	                                  Hinv);
	  }
	  // ------ end of FIT loop

	  if(convergence != 0 && convergence != 5){
	    if(a[0] < (alimit[0][0] + 2*aprec[0]) ){
	      if(convergence == 3)
	        convergence = 5;
	      else{
//	        if (logLevel>1) std::cout << "Convergence at lower tau limit!" << std::endl;
	        convergence = 4;
	      }
	    }
	    if(a[0] > (alimit[0][1] - 2*aprec[0]) ){
	      if(convergence == 3)
		convergence = 5;
	      else{
//		if (logLevel>1)
//		  std::cout << "Convergence at upper tau limit!" << std::endl;
		convergence = 4;
	      }
	    }
	  }
//	  if (m_logLevel>1)
//	    std::cout << "Convergence is " << m_convergence << std::endl;

}

Double_t 
HHKinFit::getChi2() const{
  Double_t chi2=0;
  for(std::vector<HHFitConstraint*>::const_iterator it = m_constraints.begin();it != m_constraints.end(); ++it)
    chi2 += (*it)->getChi2();
  return(chi2);
}

std::vector<HHFitObjectE*>
HHKinFit::getListOfFitObjects() const{
  return(m_fitobjects);
}

std::vector<HHFitConstraint*> 
HHKinFit::getListOfConstraints() const{
  return(m_constraints);
}

void
HHKinFit::addFitObjectE(HHFitObjectE* fitobject){
  m_fitobjects.push_back(fitobject);
}

void
HHKinFit::addConstraint(HHFitConstraint* constraint){
  m_constraints.push_back(constraint);
}

TGraph
HHKinFit::getChi2Function(int steps){
  Int_t npoints(m_fitobjects.size()*steps);
  TGraph gr(npoints);
  gr.SetName("chi2function");
  Double_t stepsize((m_fitobjects[0]->getUpperFitLimitE() - m_fitobjects[0]->getLowerFitLimitE())/steps);
  for (unsigned int i=0; i<npoints; i++){
	Double_t e = m_fitobjects[0]->getLowerFitLimitE()+ i*stepsize;
    m_fitobjects[0]->changeEandSave(e);
    Double_t chi2(this->getChi2());
//    std::cout << i << " " << e << " " << chi2 << std::endl;
    gr.SetPoint(i,e,chi2);
  }
  gr.SetMinimum(0);
  gr.GetXaxis()->SetTitle("E_{1}[GeV]");
  gr.GetYaxis()->SetTitle("#chi^{2}");
  return(gr);
}
