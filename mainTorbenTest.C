#include <iostream>
#include "HHTauTauEventGenerator.h"


int main(int argc, char* argv[])
{
	HHTauTauEventGenerator testgenerator;
    HHLorentzVector tau1test(testgenerator.getTau1());
	tau1test.Print();

  return(0);
}
