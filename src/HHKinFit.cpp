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

HHKinFit2::HHKinFit::HHKinFit()
: m_fitobjects(std::vector<HHFitObjectE*>()),
  m_constraints(std::vector<HHFitConstraint*>()){
}

///TODO: make it more general!
void 
HHKinFit2::HHKinFit::fit(){

	  //  ----------  for PSfit -----
	  const int np = 1;
	  double a[np];
	  double astart[np];
	  double alimit[np][2];
	  double aprec[np];
	  double daN[np];
	  double h[np];
	  double chi2iter[1], aMemory[np][5], g[np], H[np * np], Hinv[np * np];
	  bool noNewtonShifts = false;

	  int iter = 0;             //  number of iterations
	  int method = 1;           //  initial fit method, see PSfit()
	  int mode = 1;             //  mode =1 for start of a new fit by PSfit()
	//   int icallNewton = -1;     //  init start of Newton Method
	//   int iloop = 0;            // counter for falls to fit function




	  // fill initial tau fit parameters
	  astart[0] = m_fitobjects[0]->getE();    // energy of first tau
	  aprec[0]  = 0.1;                        //0.1                 // precision for fit
	  // fill initial step width
	  h[0] = 0.2*(m_fitobjects[0]->getUpperFitLimitE()-m_fitobjects[0]->getLowerFitLimitE());   //
	  daN[0] = 1.0;   //0.0                 // initial search direction



	  alimit[0][0] = 1.00001* m_fitobjects[0]->getLowerFitLimitE();              // tau: minimum is visible tau1 energy
	  alimit[0][1] = 0.99999* m_fitobjects[0]->getUpperFitLimitE();              //      maximum as computed above

	  // tau: check initial values against fit range
	  if (astart[0] - h[0] < alimit[0][0]) {
	    astart[0] = alimit[0][0] + h[0];
	  }
	  else if (astart[0] + h[0] > alimit[0][1]) {
	    astart[0] = alimit[0][1] - h[0];
	  }

	  for (int ip = 0; ip < np; ip++) {
	    a[ip] = astart[ip];
	  }
	  for (int ip = 0; ip < np; ip++) {
	    aMemory[ip][0] = -999.0;
	    aMemory[ip][1] = -995.0;
	    aMemory[ip][2] = -990.0;
	    aMemory[ip][3] = -985.0;
	    aMemory[ip][3] = -980.0;
	  }

	  double chi2(99999);
	  int convergence(0);
	  int printlevel(0);
	  int m_maxloops(500);
	  for (int iloop = 0; iloop < m_maxloops * 10 && iter < m_maxloops; iloop++) { // FIT loop
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

double
HHKinFit2::HHKinFit::getChi2() const{
  double chi2=0;
  for(std::vector<HHFitConstraint*>::const_iterator it = m_constraints.begin();it != m_constraints.end(); ++it)
    chi2 += (*it)->getChi2();
  return(chi2);
}

std::vector<HHKinFit2::HHFitObjectE*>
HHKinFit2::HHKinFit::getListOfFitObjects() const{
  return(m_fitobjects);
}

std::vector<HHKinFit2::HHFitConstraint*>
HHKinFit2::HHKinFit::getListOfConstraints() const{
  return(m_constraints);
}

void
HHKinFit2::HHKinFit::addFitObjectE(HHFitObjectE* fitobject){
  m_fitobjects.push_back(fitobject);
}

void
HHKinFit2::HHKinFit::addConstraint(HHFitConstraint* constraint){
  m_constraints.push_back(constraint);
}

TGraph
HHKinFit2::HHKinFit::getChi2Function(int steps){
  int npoints(m_fitobjects.size()*steps);
  TGraph gr(npoints);
  gr.SetName("chi2function");
  double stepsize((m_fitobjects[0]->getUpperFitLimitE() - m_fitobjects[0]->getLowerFitLimitE())/steps);
  for (unsigned int i=0; i<npoints; i++){
	double e = m_fitobjects[0]->getLowerFitLimitE()+ i*stepsize;
    m_fitobjects[0]->changeEandSave(e);
    double chi2(this->getChi2());
//    std::cout << i << " " << e << " " << chi2 << std::endl;
    gr.SetPoint(i,e,chi2);
  }
  gr.SetMinimum(0);
  gr.GetXaxis()->SetTitle("E_{1}[GeV]");
  gr.GetYaxis()->SetTitle("#chi^{2}");
  return(gr);
}
