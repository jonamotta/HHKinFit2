#!/bin/bash
echo "removing old files"
rm -f libHHKinFit.so run

echo "creating shared library"
g++ -fPIC -shared src/*.cpp `root-config --cflags --glibs` -I ./include -o libHHKinFit.so

echo "creating executable"
g++ main.C `root-config --cflags --glibs` -I ./include -L . -lHHKinFit  -o runHHKinFit
g++ mainTorbenTest.C `root-config --cflags --glibs` -I ./include -L . -lHHKinFit  -o runHHTauTauEventGenerator
g++ controlplots.C `root-config --cflags --glibs` -I ./include -L . -lHHKinFit  -o createcontrolplots
g++ KinFitwithEventGenerator.C `root-config --cflags --glibs` -I ./include -L . -lHHKinFit  -o KinFitwithEventGenerator
#g++ -std=c++11 compareKinFits.C `root-config --cflags --glibs` -I ./include -I ../HHKinFit/interface -L . -L ../HHKinFit -lHHKinFit2 -lHHKinFit -o compareKinFits
