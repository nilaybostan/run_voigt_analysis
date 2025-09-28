// runTestMacro.C
#include "testmacro_FullCorrections_NoJets_2.C"  // include your main macro
void runTestMacro() {
    gSystem->Load("libRoccoR.so");

    std::vector<std::string> files = {
        "root://xrootd-cms.infn.it//store/data/Run2022F/Muon/NANOAOD/22Sep2023-v2/2550000/187620f5-a4ae-467a-a56f-4569e7b32801.root",
        "root://xrootd-cms.infn.it//store/data/Run2022F/Muon/NANOAOD/22Sep2023-v2/2550000/65234aec-f4a2-4fb8-bfb8-33c63b033f9f.root",
        "root://xrootd-cms.infn.it//store/data/Run2022F/Muon/NANOAOD/22Sep2023-v2/2550000/70ef9ab5-5109-423a-9953-d11e695fb123.root"
    };

    testmacro_FullCorrections_NoJets_2(files, true);
}
