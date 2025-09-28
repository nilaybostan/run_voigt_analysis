#!/bin/bash
# runMacro.sh
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/user/n/nbostan/new_CMS/CMSSW_14_0_18/src
eval `scramv1 runtime -sh`
cd RoccoR

root -l -b -q runTestMacro.C
