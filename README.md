[![Build Status](https://travis-ci.org/bvormwald/HHKinFit2.svg?branch=master)](https://travis-ci.org/bvormwald/HHKinFit2)
# HHKinFit2

HHKinFit2 is a tool for kinematically fitting events of heavy particles via a resonance. It is especially optimized for the case where tau leptons are involved in the final state.

HHKinFit2 consists of several subpackages:
*  HHKinFit2Core
*  HHKinFit2Scenarios
*  HHKinFit2Testing
*  HHKinFit2CMSSWPlugins

## HHKinFit2Core
This package contains all the basic classes used to build a full fit.

## HHKinFit2Scenarios
This package provides fit master classes which setup a fit for a given event topology. You can get inspired by these classes if you need to setup a new topology.

## HHKinFit2Testing
This package provides some test programs.

## HHKinFit2CMSSWPlugins
This package contains CMSSW plugins which can be used to run HHKinFit2 directy within CMSSW.

# Installation
HHKinFit2 can be used in standalone mode or within CMSSW. If you want to compile the code in standalone mode, just run setup.sh as well as compile.sh to compile the code. Compiling within CMSSW works just using the scram build system.