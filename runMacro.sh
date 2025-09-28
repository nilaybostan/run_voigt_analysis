#!/bin/bash
# runMacro.sh

# Set up CMSSW environment
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/user/n/nbostan/new_CMS/CMSSW_14_0_18/src
eval `scramv1 runtime -sh`
cd RoccoR

# Make sure output folder exists
mkdir -p output

# Run your ROOT macro
root -l -b -q 'testmacro_FullCorrections_NoJets_2.C+'

# Move all PDFs into output/
mv *.pdf output/ 2>/dev/null || true
