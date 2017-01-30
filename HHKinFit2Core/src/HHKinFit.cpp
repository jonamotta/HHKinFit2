#include "HHKinFit2/HHKinFit2Core/interface/HHKinFit.h"
#include "HHKinFit2/HHKinFit2Core/interface/HHFitConstraint4Vector.h"
#include "HHKinFit2/HHKinFit2Core/interface/HHFitConstraintEHardM.h"
#include "HHKinFit2/HHKinFit2Core/interface/HHFitConstraint.h"
#include "HHKinFit2/HHKinFit2Core/interface/HHFitConstraintLikelihood.h"
#include "HHKinFit2/HHKinFit2Core/interface/HHFitObjectEConstM.h"
#include "HHKinFit2/HHKinFit2Core/interface/HHFitObjectE.h"
#include "HHKinFit2/HHKinFit2Core/interface/HHFitObject.h"
#include "HHKinFit2/HHKinFit2Core/interface/HHFitObjectComposite.h"
#include "HHKinFit2/HHKinFit2Core/interface/PSMath.h"
#include "HHKinFit2/HHKinFit2Core/interface/exceptions/HHEnergyRangeException.h"
#include "HHKinFit2/HHKinFit2Core/interface/exceptions/HHInvMConstraintException.h"
#include "HHKinFit2/HHKinFit2Core/interface/exceptions/HHLimitSettingException.h"

#include "TAxis.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include <iomanip>
#include <iostream>
#include <sstream>

HHKinFit2::HHKinFit::HHKinFit()
: m_fitobjects(std::vector<HHFitObjectE*>()),
  m_constraints(std::vector<HHFitConstraint*>()),
  m_chi2(99999),
  m_convergence(0),
  m_printlevel(0),
  m_maxloops(10000){
}

void
HHKinFit2::HHKinFit::setPrintLevel(int printlevel){
  m_printlevel=printlevel;
}

void 
HHKinFit2::HHKinFit::fit(){
  m_chi2 = 99999;
  m_convergence = -10;
  const int np = m_fitobjects.size();

  ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Fumili2");
  //ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  min->SetMaxFunctionCalls(1000000);
  min->SetMaxIterations(100000);
  min->SetTolerance(0.0001);
 
  ROOT::Math::Functor f(this,&HHKinFit2::HHKinFit::minfunction,np); 
  min->SetFunction(f);
 
  // Set the free variables to be minimized!
  for (int i=0; i<np; i++){
    //min->SetVariable(0,Form("E%d",i),m_fitobjects[i]->getInitStart(), m_fitobjects[i]->getInitStepWidth());
    min->SetLimitedVariable(i,Form("E%d",i),m_fitobjects[i]->getInitStart(), m_fitobjects[i]->getInitStepWidth(), m_fitobjects[i]->getLowerFitLimitE(), m_fitobjects[i]->getUpperFitLimitE());
  }
 
  m_convergence = min->Minimize(); 
  bool respectLimits=false;
  m_chi2=this->getChi2(respectLimits);
}

double
HHKinFit2::HHKinFit::getChi2(bool respectLimits) const{
  double chi2=0;

  for(std::vector<HHFitConstraint*>::const_iterator it = m_constraints.begin();
      it != m_constraints.end(); 
      ++it)
  {
    (*it)->prepare(respectLimits);
  }
  for(std::vector<HHFitConstraint*>::const_iterator it = m_constraints.begin();it != m_constraints.end(); ++it)
    chi2 += (*it)->getChi2();

  return(chi2);
}

double HHKinFit2::HHKinFit::minfunction(const double *xx){
  double chi2=-1;
  bool respectLimits=false;

  for (unsigned int i=0; i<m_fitobjects.size();i++){
    try{
      m_fitobjects[i]->changeEandSave(xx[i],respectLimits);
      chi2=this->getChi2(respectLimits);
    }
    catch(HHKinFit2::HHEnergyRangeException const& e){
      std::cout << e.what() << std::endl;
      throw(e);
    }
  }

  return(chi2);
}

void
HHKinFit2::HHKinFit::printChi2() const{
  double chi2=0;

  for(std::vector<HHFitConstraint*>::const_iterator it = m_constraints.begin();it != m_constraints.end(); ++it)
    (*it)->prepare();

  for(std::vector<HHFitConstraint*>::const_iterator it = m_constraints.begin();it != m_constraints.end(); ++it){
    std::cout << (*it)->getChi2() << std::endl;
    chi2 += (*it)->getChi2();
    (*it)->printChi2();
  }
  std::cout << chi2 << std::endl;
  std::cout << "----------------------------------------------------------------------------------------------"<<std::endl;
}

double
HHKinFit2::HHKinFit::getL(bool respectLimits) const{
  double L=1;

  for(std::vector<HHFitConstraint*>::const_iterator it = m_constraints.begin();it != m_constraints.end(); ++it)
    (*it)->prepare(respectLimits);

  for(std::vector<HHFitConstraint*>::const_iterator it = m_constraints.begin();it != m_constraints.end(); ++it)
    L *= (*it)->getLikelihood();

  return(L);
}


int
HHKinFit2::HHKinFit::getConvergence() const{
  return(m_convergence);
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

TGraph*
HHKinFit2::HHKinFit::getChi2Function(int steps){
  int npoints(m_fitobjects.size()*steps);
  TGraph* gr = new TGraph(npoints);
  gr->SetName("chi2function");
  double stepsize((m_fitobjects[0]->getUpperFitLimitE() - m_fitobjects[0]->getLowerFitLimitE())/steps);
  for (int i=0; i<npoints; i++){
    double e = 1.00001*m_fitobjects[0]->getLowerFitLimitE()+ i*stepsize;
    m_fitobjects[0]->changeEandSave(e);
    double chi2(this->getChi2());
    //    std::cout << i << " " << e << " " << chi2 << std::endl;
    gr->SetPoint(i,e,chi2);
  }
  gr->SetMinimum(0);
  gr->GetXaxis()->SetTitle("E_{1}[GeV]");
  gr->GetYaxis()->SetTitle("#chi^{2}");
  return(gr);
}

TGraph2D*
HHKinFit2::HHKinFit::getChi2Function(int* steps, double* mins, double* maxs)
{
  TGraph2D* gr=NULL;
//  TObject* obj;
  /*
  if(m_fitobjects.size() == 1)
  {
    int npoints = steps[0];
    obj = new TGraph(npoints);
    TGraph* gr = (TGraph*)obj;
    gr->SetName("chi2function");
    double stepsize((m_fitobjects[0]->getUpperFitLimitE() - m_fitobjects[0]->getLowerFitLimitE())/npoints);
    for (int i=0; i<npoints; i++){
      double e = 1.00001*m_fitobjects[0]->getLowerFitLimitE()+ i*stepsize;
      m_fitobjects[0]->changeEandSave(e);
      double chi2(this->getChi2());
      //    std::cout << i << " " << e << " " << chi2 << std::endl;
      gr->SetPoint(i,e,chi2);
    }
    gr->SetMinimum(0);
    gr->GetXaxis()->SetTitle("E_{1}[GeV]");
    gr->GetYaxis()->SetTitle("#chi^{2}");
  }
  */
  if(m_fitobjects.size()==2)
  {
    if(mins[0] < m_fitobjects[0]->getLowerFitLimitE())
    {
      mins[0] = m_fitobjects[0]->getLowerFitLimitE();
    }
    if(mins[1] < m_fitobjects[1]->getLowerFitLimitE())
    {
      mins[1] = m_fitobjects[1]->getLowerFitLimitE();
    }
    if(maxs[0] > m_fitobjects[0]->getUpperFitLimitE())
    {
      maxs[0] = m_fitobjects[0]->getUpperFitLimitE();
    }
    if(maxs[1] > m_fitobjects[1]->getUpperFitLimitE())
    {
      maxs[1] = m_fitobjects[1]->getUpperFitLimitE();
    }
    
    int stepsPar1 = steps[0];
    int stepsPar2 = steps[1];

    gr = new TGraph2D(stepsPar1 * stepsPar2);
//    TGraph2D* gr = (TGraph2D*)obj;
    gr->SetName("chi2function");
    double stepsize1 = (maxs[0] - mins[0])/stepsPar1;
    double stepsize2 = (maxs[1] - mins[1])/stepsPar2;
    for(int i = 0;
	i < stepsPar1;
	i++)
    {
      double e1 = 1.00001*mins[0] + i*stepsize1;
      m_fitobjects[0]->changeEandSave(e1);
      for(int j = 0;
	  j < stepsPar2;
	  j++)
      {
	double e2 = 1.00001*mins[1] + j*stepsize2;
	m_fitobjects[1]->changeEandSave(e2);
	double chi2 = this->getChi2();
	gr->SetPoint(i*stepsPar1+j, e1, e2, chi2);
      }
    }  
    gr->GetZaxis()->SetRangeUser(0.0, 100.0);
  }
  return(gr);
}


TGraph*
HHKinFit2::HHKinFit::getLFunction(int steps){
  int npoints(m_fitobjects.size()*steps);
  TGraph* gr = new TGraph(npoints);
  gr->SetName("Lfunction");
  double stepsize((m_fitobjects[0]->getUpperFitLimitE() - m_fitobjects[0]->getLowerFitLimitE())/steps);
  for (int i=0; i<npoints; i++){
    double e = 1.00001*m_fitobjects[0]->getLowerFitLimitE()+ i*stepsize;
    m_fitobjects[0]->changeEandSave(e);
    double L(this->getL());
    //    std::cout << i << " " << e << " " << chi2 << std::endl;
    gr->SetPoint(i,e,L);
  }
  gr->SetMinimum(0);
  gr->GetXaxis()->SetTitle("E_{1}[GeV]");
  gr->GetYaxis()->SetTitle("L");
  return(gr);
}
