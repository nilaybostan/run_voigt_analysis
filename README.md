For CMS Dimuon Analysis, these codes can be applied to data and Drell-Yan (MC) files:

1) Muon Efficiency Measurement for Higgs to Dimuon Analysis in CMS (with Tag-and-Probe):

Execute the code by using the command: 

g++ -o run_voigt_analysis run_voigt_analysis.cpp `root-config --cflags --libs` -lRooFit -lRooFitCore

Then run the code by using: 

./run_voigt_analysis 

2) Momentum corrections, applying Rochester corrections:

first, compile (use the `Makefile` in the folder)

make
export LD_LIBRARY_PATH=$PWD:$LD_LIBRARY_PATH
export ROOT_INCLUDE_PATH=$PWD:$ROOT_INCLUDE_PATH

you need to be in cmssw to run this code, so you should do `cmsenv` and then, run the code by using the command below: 

root -l

.L testmacro_Rocco_Zregion.C

testmacro_Rocco_Zregion("root://xrootd-cms.infn.it//store/data/Run2022F/Muon/NANOAOD/*.root", false); 

(this is for data file)

For Monte Carlo, make the second parameter `true` as following: 

testmacro_Rocco_Zregion("root://xrootd-cms.infn.it//store/data/Run2022F/Muon/NANOAOD/DY_.root", true); 

3) All Muon Corrections applied code: testmacro_FullCorrections_NoJets_2.C, just run using below:
   
root -l

root [0] gSystem->Load("$CMSSW_BASE/src/RoccoR/libRoccoR.so");

root [1] .L testmacro_FullCorrections_NoJets_2.C+

root [2] std::vector<std::string> files = {"root://xrootd-cms.infn.it//store/data/Run2022F/Muon/NANOAOD/22Sep2023-v2/2550000/187620f5-a4ae-467a-a56f-4569e7b32801.root"};

root [3] testmacro_FullCorrections_NoJets_2(files, true);

Note: testmacro_FullCorrections_NoJets_2(files, false);  -> For MC files (for instance Drell-Yan)

(now the correct macro file: testmacro_FullCorrections_NoJets_corrected.C)


4) Condor Job Submission Files:

runMacro.sh

runMacro.sub

runTestMacro.C

Condor job submission files can be run by using the command: condor_submit runMacro.sub

then, you can check your job summary by using the command in bash: condor_q 

