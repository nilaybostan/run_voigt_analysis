#include "testmacro_FullCorrections_NoJets.C"  

void runTestmacro_MC() {
    std::vector<std::string> files = {
        "root://xrootd-cms.infn.it//store/mc/Run3Summer22NanoAODv12/DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/NANOAODSIM/130X_mcRun3_2022_realistic_v5-v2/2520000/181836d0-879c-4c05-93c1-69956def2efb.root"
    };
    testmacro_FullCorrections_NoJets(files, false);
}
