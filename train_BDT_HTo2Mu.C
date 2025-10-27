#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TCut.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TCanvas.h"
#include "TGraph.h"

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
    // Sinyal ve background dosyalarını aç
    // --------------------------------------------------------------
    TFile* inputSig = TFile::Open("root://xrootd-cms.infn.it//store/mc/Run3Summer22NanoAODv12/VBFHto2Mu_M-125_TuneCP5_withDipoleRecoil_13p6TeV_powheg-pythia8/NANOAODSIM/130X_mcRun3_2022_realistic_v5-v3/50000/0ea52cdc-3d79-497f-9275-b6fd2041e9a0.root");
    TFile* inputBkg = TFile::Open("root://xrootd-cms.infn.it//store/mc/Run3Summer22NanoAODv12/DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/NANOAODSIM/130X_mcRun3_2022_realistic_v5-v2/2520000/aec6684e-3416-4b7c-a96f-78f261659521.root");

    TTree* signalTree = (TTree*)inputSig->Get("Events");
    TTree* backgroundTree = (TTree*)inputBkg->Get("Events");

    dataloader->AddSignalTree(signalTree, 1.0);
    dataloader->AddBackgroundTree(backgroundTree, 1.0);

    // --------------------------------------------------------------
    // Değişkenler (Alt$() ile güvenli)
    // --------------------------------------------------------------

    dataloader->AddVariable(
        "sqrt(2*Alt$(Muon_pt[0],0)*Alt$(Muon_pt[1],0)*(cosh(Alt$(Muon_eta[0],0)-Alt$(Muon_eta[1],0)) - cos(Alt$(Muon_phi[0],0)-Alt$(Muon_phi[1],0))))",
        "dimuon_mass", "GeV", 'F');

    dataloader->AddVariable(
        "sqrt(pow(Alt$(Muon_pt[0],0)*cos(Alt$(Muon_phi[0],0)) + Alt$(Muon_pt[1],0)*cos(Alt$(Muon_phi[1],0)),2) + "
        "pow(Alt$(Muon_pt[0],0)*sin(Alt$(Muon_phi[0],0)) + Alt$(Muon_pt[1],0)*sin(Alt$(Muon_phi[1],0)),2))",
        "dimuon_pt", "GeV", 'F');

    dataloader->AddVariable("Alt$(Muon_pt[0],0)", "muon1_pt", "GeV", 'F');
    dataloader->AddVariable("Alt$(Muon_pt[1],0)", "muon2_pt", "GeV", 'F');
    dataloader->AddVariable("Alt$(Muon_eta[0],0)", "muon1_eta", "", 'F');
    dataloader->AddVariable("Alt$(Muon_eta[1],0)", "muon2_eta", "", 'F');

    dataloader->AddVariable(
        "sqrt(pow(Alt$(Muon_eta[0],0)-Alt$(Muon_eta[1],0),2)+pow(acos(cos(Alt$(Muon_phi[0],0)-Alt$(Muon_phi[1],0))),2))",
        "deltaR_mumu", "", 'F');

    dataloader->AddVariable("Alt$(nJet,0)", "nJets", "", 'I');
    dataloader->AddVariable("Alt$(MET_pt,0)", "met", "GeV", 'F');

    // --------------------------------------------------------------
    // Preselection
    // --------------------------------------------------------------
    TCut preselection =
        "(Alt$(nMuon,0) >= 2)"
        " && (Alt$(Muon_charge[0],0) != Alt$(Muon_charge[1],0))"
        " && (Alt$(Muon_pt[0],0) > 20) && (Alt$(Muon_pt[1],0) > 15)"
        " && (abs(Alt$(Muon_eta[0],0)) < 2.4) && (abs(Alt$(Muon_eta[1],0)) < 2.4)"
        " && (sqrt(2*Alt$(Muon_pt[0],0)*Alt$(Muon_pt[1],0)*(cosh(Alt$(Muon_eta[0],0)-Alt$(Muon_eta[1],0)) - cos(Alt$(Muon_phi[0],0)-Alt$(Muon_phi[1],0)))) > 70)"
        " && (sqrt(2*Alt$(Muon_pt[0],0)*Alt$(Muon_pt[1],0)*(cosh(Alt$(Muon_eta[0],0)-Alt$(Muon_eta[1],0)) - cos(Alt$(Muon_phi[0],0)-Alt$(Muon_phi[1],0)))) < 115)";

    // --------------------------------------------------------------
    // Eğitim / test ayrımı
    // --------------------------------------------------------------
    dataloader->PrepareTrainingAndTestTree(
        preselection,
        "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V"
    );

    // --------------------------------------------------------------
    // BDT yöntemi
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
    // ROC eğrisi çizimi
    // --------------------------------------------------------------
    TFile* resultsFile = TFile::Open(outfileName);
    TKey* key = (TKey*)resultsFile->GetListOfKeys()->FindObject("dataset/Method_BDT/MVA_BDT");
    if (key) {
        TH1* hROC = (TH1*)resultsFile->Get("dataset/Method_BDT/MVA_BDT/MVA_BDT_TrainingROC");
        if (hROC) {
            TCanvas* cROC = new TCanvas("cROC", "ROC Curve", 800, 600);
            hROC->Draw("AL");
            cROC->SaveAs("ROC_curve.png");
        }
    }

    // --------------------------------------------------------------
    // Çıktıyı kapat
    // --------------------------------------------------------------
    outputFile->Close();
    std::cout << "BDT eğitimi tamamlandı! Sonuç dosyası: " << outfileName << std::endl;
    TMVA::TMVAGui(outfileName);
}
