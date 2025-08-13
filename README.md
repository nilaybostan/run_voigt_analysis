For CMS Dimuon Analysis, these codes can be applied to data and Drell-Yan (MC) files:

1) Muon Efficiency Measurement for Higgs to Dimuon Analysis in CMS (with Tag-and-Probe)

Execute the code by using the command: 

g++ -o run_voigt_analysis run_voigt_analysis.cpp `root-config --cflags --libs` -lRooFit -lRooFitCore

Then run the code by using: 

./run_voigt_analysis 

2) Momentum corrections, applying Rochester corrections

you need to be in cmsssw to run this code, so you should do first `cmsenv` and then, run the code by using the command below: 

root -l

.L testmacro_Rocco_Zregion.C

testmacro_Rocco_Zregion("root://xrootd-cms.infn.it//store/data/Run2022F/Muon/NANOAOD/*.root", false); 

(this is for data file)

For Monte Carlo, make the second parameter `true` as following: 

testmacro_Rocco_Zregion("root://xrootd-cms.infn.it//store/data/Run2022F/Muon/NANOAOD/DY_.root", true); 
