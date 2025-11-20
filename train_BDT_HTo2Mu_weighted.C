#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TH1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include <iostream>
#include <vector>
#include <cmath>

void train_BDT_HTo2Mu_weighted() {

    TMVA::Tools::Instance();
    TString outfileName("TMVA_HTo2Mu_BDT.root");
    TFile* outputFile = TFile::Open(outfileName, "RECREATE");

    TMVA::Factory factory("TMVAClassification", outputFile,
        "!V:Color:DrawProgressBar:AnalysisType=Classification");

    TMVA::DataLoader* dataloader = new TMVA::DataLoader("dataset");

    // --------------------------------------------------------------
    // SIGNAL FILES
    // --------------------------------------------------------------
    std::vector<TString> sigFiles = {
        "root://xrootd-cms.infn.it//store/mc/Run3Summer22NanoAODv12/VBFHto2Mu_M-125_TuneCP5_withDipoleRecoil_13p6TeV_powheg-pythia8/NANOAODSIM/130X_mcRun3_2022_realistic_v5-v3/50000/07efcc80-1811-4648-90bb-9a52a423a65f.root",
        "root://xrootd-cms.infn.it//store/mc/Run3Summer22NanoAODv12/GluGluHto2Mu_M-125_TuneCP5_13p6TeV_powheg-pythia8/NANOAODSIM/130X_mcRun3_2022_realistic_v5-v3/50000/154d3c87-08bf-4277-ba6f-d798a21fa339.root"
    };

    TChain* sigTree = new TChain("Events");
    for (auto &file : sigFiles) sigTree->Add(file);
    dataloader->AddSignalTree(sigTree, 1.0);

    // --------------------------------------------------------------
    // BACKGROUND FILES
    // --------------------------------------------------------------
    std::vector<TString> bkgFiles = {
        "root://xrootd-cms.infn.it//store/mc/Run3Summer22NanoAODv12/DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/NANOAODSIM/130X_mcRun3_2022_realistic_v5-v2/30000/a642010b-525a-4da9-bb72-ad5e204c12c2.root",
   "root://xrootd-cms.infn.it//store/mc/Run3Summer22NanoAODv12/TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8/NANOAODSIM/130X_mcRun3_2022_realistic_v5-v2/2520000/66b834d6-61f7-4109-b5ae-54a150d4814b.root",
   "root://xrootd-cms.infn.it//store/mc/Run3Summer22NanoAODv12/EWK_2L2J_TuneCH3_13p6TeV_madgraph-herwig7/NANOAODSIM/130X_mcRun3_2022_realistic_v5-v2/130000/f7ad2334-785b-4377-bbf9-136725a1892f.root",
     "root://xrootd-cms.infn.it//store/mc/Run3Summer22NanoAODv12/WZto2L2Q_TuneCP5_13p6TeV_powheg-pythia8/NANOAODSIM/130X_mcRun3_2022_realistic_v5-v2/2540000/43df5e8d-7e75-4f0f-a119-c7a12fde30bc.root",
    "root://xrootd-cms.infn.it//store/mc/Run3Summer22NanoAODv12/ZZto2L2Q_TuneCP5_13p6TeV_powheg-pythia8/NANOAODSIM/130X_mcRun3_2022_realistic_v5-v2/2520000/46fe5fe8-c2ec-4df3-a7a6-32ce6cf32ca5.root",
    "root://xrootd-cms.infn.it//store/mc/Run3Summer22NanoAODv12/ZZto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8/NANOAODSIM/130X_mcRun3_2022_realistic_v5-v2/2520000/6aeb6e62-3d36-453b-bdb6-0f503de49661.root"
    };

    TChain* bkgTree = new TChain("Events");
    for (auto &file : bkgFiles) bkgTree->Add(file);
    dataloader->AddBackgroundTree(bkgTree, 1.0);

    // --------------------------------------------------------------
    // EVENT WEIGHTS (absolute value to remove negative weights)
    // --------------------------------------------------------------
    dataloader->SetSignalWeightExpression("abs(genWeight)");
    dataloader->SetBackgroundWeightExpression("abs(genWeight)");

    // --------------------------------------------------------------
    // VARIABLES (Alt$ for safe access)
    // --------------------------------------------------------------
    dataloader->AddVariable(
        "sqrt(2*Alt$(Muon_pt[0],0)*Alt$(Muon_pt[1],0)*(cosh(Alt$(Muon_eta[0],0)-Alt$(Muon_eta[1],0)) - cos(Alt$(Muon_phi[0],0)-Alt$(Muon_phi[1],0))))",
        "dimuon_mass", "GeV", 'F'
    );
    dataloader->AddVariable("Alt$(Muon_pt[0],0)", "muon1_pt", "GeV", 'F');
    dataloader->AddVariable("Alt$(Muon_pt[1],0)", "muon2_pt", "GeV", 'F');
    dataloader->AddVariable("Alt$(Muon_eta[0],0)", "muon1_eta", "", 'F');
    dataloader->AddVariable("Alt$(Muon_eta[1],0)", "muon2_eta", "", 'F');
    dataloader->AddVariable(
        "sqrt(pow(Alt$(Muon_eta[0],0)-Alt$(Muon_eta[1],0),2)+pow(acos(cos(Alt$(Muon_phi[0],0)-Alt$(Muon_phi[1],0))),2))",
        "deltaR_mumu", "", 'F'
    );
    dataloader->AddVariable("Alt$(nJet,0)", "nJets", "", 'I');
    dataloader->AddVariable("Alt$(MET_pt,0)", "met", "GeV", 'F');

    // --------------------------------------------------------------
    // PRESELECTION
    // --------------------------------------------------------------
    TCut preselection =
        "(Alt$(nMuon,0) >= 2)"
        " && (Alt$(Muon_charge[0],0) != Alt$(Muon_charge[1],0))"
        " && (Alt$(Muon_pt[0],0) > 20) && (Alt$(Muon_pt[1],0) > 15)"
        " && (abs(Alt$(Muon_eta[0],0)) < 2.4) && (abs(Alt$(Muon_eta[1],0)) < 2.4)"
        " && (sqrt(2*Alt$(Muon_pt[0],0)*Alt$(Muon_pt[1],0)*(cosh(Alt$(Muon_eta[0],0)-Alt$(Muon_eta[1],0)) - cos(Alt$(Muon_phi[0],0)-Alt$(Muon_phi[1],0)))) > 70)"
        " && (sqrt(2*Alt$(Muon_pt[0],0)*Alt$(Muon_pt[1],0)*(cosh(Alt$(Muon_eta[0],0)-Alt$(Muon_eta[1],0)) - cos(Alt$(Muon_phi[0],0)-Alt$(Muon_phi[1],0)))) < 115)";

    dataloader->PrepareTrainingAndTestTree(
        preselection,
        "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V"
    );

    // --------------------------------------------------------------
    // BOOK BDT
    // --------------------------------------------------------------
    factory.BookMethod(
        dataloader,
        TMVA::Types::kBDT,
        "BDT",
        "!H:!V:NTrees=600:MaxDepth=4:MinNodeSize=2.5%:"
        "BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:"
        "BaggedSampleFraction=0.5:SeparationType=GiniIndex"
    );

    // --------------------------------------------------------------
    // TRAIN / TEST / EVALUATE
    // --------------------------------------------------------------
    factory.TrainAllMethods();
    factory.TestAllMethods();
    factory.EvaluateAllMethods();

    outputFile->Close();
    std::cout << "\n=== TMVA training complete ===\n";

    // --------------------------------------------------------------
    // DRAW BDT SCORE HISTOGRAMS
    // --------------------------------------------------------------
    TFile* f = TFile::Open("TMVA_HTo2Mu_BDT.root");
    if (!f) { std::cerr << "Cannot open TMVA output.\n"; return; }

    TH1* hTrainSig = (TH1*)f->Get("dataset/Method_BDT/BDT/MVA_BDT_Train_S");
    TH1* hTrainBkg = (TH1*)f->Get("dataset/Method_BDT/BDT/MVA_BDT_Train_B");
    TH1* hTestSig  = (TH1*)f->Get("dataset/Method_BDT/BDT/MVA_BDT_S");
    TH1* hTestBkg  = (TH1*)f->Get("dataset/Method_BDT/BDT/MVA_BDT_B");

    TCanvas* c1 = new TCanvas("c1", "BDT Score", 900, 700);
    hTrainSig->SetLineColor(kBlue);
    hTrainBkg->SetLineColor(kRed);
    hTestSig->SetLineColor(kBlue+1);
    hTestBkg->SetLineColor(kRed+1);

    hTrainSig->SetLineWidth(2);
    hTrainBkg->SetLineWidth(2);
    hTestSig->SetLineWidth(2);
    hTestBkg->SetLineWidth(2);

    hTrainSig->Draw("HIST");
    hTrainBkg->Draw("HIST SAME");
    hTestSig->Draw("HIST SAME");
    hTestBkg->Draw("HIST SAME");

    c1->SaveAs("BDT_score.png");
    std::cout << "BDT score histogram saved: BDT_score.png\n";

    // --------------------------------------------------------------
    // DRAW ROC CURVE AND CALCULATE AUC
    // --------------------------------------------------------------
    TObject* objROC = f->Get("dataset/Method_BDT/BDT/MVA_BDT_TrainingROC");

    if(objROC){
        TCanvas* c2 = new TCanvas("c2", "ROC Curve", 800, 600);

        if(objROC->InheritsFrom("TGraph")){
            TGraph* rocGraph = (TGraph*)objROC;
            rocGraph->SetLineWidth(3);
            rocGraph->SetLineColor(kBlack);
            rocGraph->Draw("AL");
            c2->SaveAs("ROC_curve.png");
            std::cout << "ROC curve saved: ROC_curve.png\n";

            double auc = 0.0;
            for(int i=1; i<rocGraph->GetN(); i++){
                double x1,y1,x2,y2;
                rocGraph->GetPoint(i-1,x1,y1);
                rocGraph->GetPoint(i,x2,y2);
                auc += 0.5*(y1+y2)*(x2-x1);
            }
            std::cout << "AUC = " << auc << std::endl;

        } else if(objROC->InheritsFrom("TH1")){
            TH1* hROC = (TH1*)objROC;
            hROC->SetLineWidth(3);
            hROC->SetLineColor(kBlack);
            hROC->Draw("HIST");
            c2->SaveAs("ROC_curve.png");
            std::cout << "ROC curve saved: ROC_curve.png\n";

            int nBins = hROC->GetNbinsX();
            double auc = 0.0;
            for(int i=1; i<nBins; i++){
                double x1 = hROC->GetBinLowEdge(i);
                double x2 = hROC->GetBinLowEdge(i+1);
                double y1 = hROC->GetBinContent(i);
                double y2 = hROC->GetBinContent(i+1);
                auc += 0.5*(y1+y2)*(x2-x1);
            }
            std::cout << "AUC (approx.) = " << auc << std::endl;
        } else {
            std::cout << "ROC object is not TGraph or TH1, cannot compute AUC.\n";
        }
    } else {
        std::cout << "WARNING: ROC curve object not found.\n";
    }

    TMVA::TMVAGui("TMVA_HTo2Mu_BDT.root");
}
