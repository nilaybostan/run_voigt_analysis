#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TCut.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"

// how to run: 
// root -l
// root [0] .L train_BDT_HTo2Mu.C+
// root [1] train_BDT_HTo2Mu()

void train_BDT_HTo2Mu() {

    // --------------------------------------------------------------
    // TMVA başlat
    // --------------------------------------------------------------
    TMVA::Tools::Instance();
    TString outfileName("TMVA_HTo2Mu_BDT.root");
    TFile* outputFile = TFile::Open(outfileName, "RECREATE");

    TMVA::Factory factory("TMVAClassification", outputFile,
        "!V:!Silent:Color:DrawProgressBar:AnalysisType=Classification");

    TMVA::DataLoader* dataloader = new TMVA::DataLoader("dataset");

    // --------------------------------------------------------------
    // Sinyal (H→μμ) ve background (Z→μμ, DY, ttbar, vb.) dosyalarını aç
    // --------------------------------------------------------------
    TFile* inputSig = TFile::Open("root://xrootd-cms.infn.it//store/mc/Run3Summer22NanoAODv12/VBFHto2Mu_M-125_TuneCP5_withDipoleRecoil_13p6TeV_powheg-pythia8/NANOAODSIM/130X_mcRun3_2022_realistic_v5-v3/50000/07efcc80-1811-4648-90bb-9a52a423a65f.root");
    TFile* inputBkg = TFile::Open("root://xrootd-cms.infn.it//store/mc/Run3Summer22NanoAODv12/DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/NANOAODSIM/130X_mcRun3_2022_realistic_v5-v2/2520000/aec6684e-3416-4b7c-a96f-78f261659521.root");

    TTree* signalTree = (TTree*)inputSig->Get("Events");
    TTree* backgroundTree = (TTree*)inputBkg->Get("Events");

    dataloader->AddSignalTree(signalTree, 1.0);
    dataloader->AddBackgroundTree(backgroundTree, 1.0);

    // --------------------------------------------------------------
    // Kullanılacak değişkenleri tanımla (NanoAOD uyumlu)
    // --------------------------------------------------------------

    // m(μμ)
    dataloader->AddVariable(
        "sqrt(2*Muon_pt[0]*Muon_pt[1]*(cosh(Muon_eta[0]-Muon_eta[1]) - cos(Muon_phi[0]-Muon_phi[1])))",
        "dimuon_mass", "GeV", 'F');

    // pT(μμ)
    dataloader->AddVariable(
        "sqrt( pow(Muon_pt[0]*cos(Muon_phi[0]) + Muon_pt[1]*cos(Muon_phi[1]), 2) + "
        "pow(Muon_pt[0]*sin(Muon_phi[0]) + Muon_pt[1]*sin(Muon_phi[1]), 2) )",
        "dimuon_pt", "GeV", 'F');

    // Tekil muon değişkenleri
    dataloader->AddVariable("Muon_pt[0]", "muon1_pt", "GeV", 'F');
    dataloader->AddVariable("Muon_pt[1]", "muon2_pt", "GeV", 'F');
    dataloader->AddVariable("Muon_eta[0]", "muon1_eta", "", 'F');
    dataloader->AddVariable("Muon_eta[1]", "muon2_eta", "", 'F');

    // ΔR(μ, μ)
    dataloader->AddVariable(
        "sqrt( pow(Muon_eta[0]-Muon_eta[1],2) + pow(acos(cos(Muon_phi[0]-Muon_phi[1])),2) )",
        "deltaR_mumu", "", 'F');

    // Jet sayısı
    dataloader->AddVariable("nJet", "nJets", "", 'I');

    // Missing ET (isteğe bağlı)
    dataloader->AddVariable("MET_pt", "met", "GeV", 'F');

    // --------------------------------------------------------------
    // Basit olay seçimi (preselection)
    // --------------------------------------------------------------
    //
    // TMVA ROOT’un TCut sınıfını kullanır. TString yerine TCut olmalı.

    TCut preselection =
        "(nMuon >= 2)"
        " && (Muon_charge[0] != Muon_charge[1])"
        " && (Muon_pt[0] > 20) && (Muon_pt[1] > 15)"
        " && (abs(Muon_eta[0]) < 2.4) && (abs(Muon_eta[1]) < 2.4)"
        " && (sqrt(2*Muon_pt[0]*Muon_pt[1]*(cosh(Muon_eta[0]-Muon_eta[1]) - cos(Muon_phi[0]-Muon_phi[1]))) > 70)"
        " && (sqrt(2*Muon_pt[0]*Muon_pt[1]*(cosh(Muon_eta[0]-Muon_eta[1]) - cos(Muon_phi[0]-Muon_phi[1]))) < 115)";

    // --------------------------------------------------------------
    // Eğitim / test ayrımı
    // --------------------------------------------------------------
    dataloader->PrepareTrainingAndTestTree(
        preselection,
        "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V"
    );

    // --------------------------------------------------------------
    // BDT yöntemi (Boosted Decision Trees)
    // --------------------------------------------------------------
    factory.BookMethod(dataloader, TMVA::Types::kBDT, "BDT",
        "!H:!V:NTrees=500:MaxDepth=3:BoostType=Grad:"
        "Shrinkage=0.1:SeparationType=GiniIndex:nCuts=20");

    // --------------------------------------------------------------
    // Eğit, test et, değerlendir
    // --------------------------------------------------------------
    factory.TrainAllMethods();
    factory.TestAllMethods();
    factory.EvaluateAllMethods();

    // --------------------------------------------------------------
    // Çıktı ve ROC GUI
    // --------------------------------------------------------------
    outputFile->Close();
    std::cout << "BDT eğitimi tamamlandı! Sonuç dosyası: " << outfileName << std::endl;
    TMVA::TMVAGui(outfileName);
}
