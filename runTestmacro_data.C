#include "testmacro_FullCorrections_NoJets.C"  

void runTestmacro() {
    std::vector<std::string> files = {
        "root://xrootd-cms.infn.it//store/data/Run2022F/Muon/NANOAOD/22Sep2023-v2/2550000/187620f5-a4ae-467a-a56f-4569e7b32801.root"
    };
    testmacro_FullCorrections_NoJets(files, true);
}
