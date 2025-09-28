#!/bin/bash
# runMacro.sh
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/user/n/nbostan/new_CMS/CMSSW_14_0_18/src
eval `scramv1 runtime -sh`
cd RoccoR

#!/bin/bash
mkdir -p output   # make sure the folder exists

# Run your ROOT macro
root -l -b -q runTestMacro.C

# Move any PDFs into output/
mv *.pdf output/ 2>/dev/null || true
