#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TCut.h"
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

void train_DNN_HTo2Mu_weighted_2()
{
    TMVA::Tools::Instance();

    // ------------------------------------------------------
    // OUTPUT FILE
    // ------------------------------------------------------
    TMVA::Tools::Instance();
    TString outfileName("TMVA_HTo2Mu_DNN.root");
    TFile* outputFile = TFile::Open(outfileName, "RECREATE");

    TMVA::Factory factory("TMVAClassification", outputFile,
        "!V:Color:DrawProgressBar:AnalysisType=Classification");

    TMVA::DataLoader* dataloader = new TMVA::DataLoader("dataset");


    // ------------------------------------------------------
    // SIGNAL
    // ------------------------------------------------------
     
    std::vector<TString> sigFiles = {"*.root"};

    TChain* sigTree = new TChain("Events");
    for (auto &file : sigFiles) sigTree->Add(file);
    dataloader->AddSignalTree(sigTree, 1.0);

    // --------------------------------------------------------------
    // BACKGROUND FILES
    // --------------------------------------------------------------
    std::vector<TString> bkgFiles = {"*.root"};

    TChain* bkgTree = new TChain("Events");
    for (auto &file : bkgFiles) bkgTree->Add(file);

    dataloader->AddBackgroundTree(bkgTree, 1.0);

    // ------------------------------------------------------
    // WEIGHTS
    // ------------------------------------------------------
    dataloader->SetSignalWeightExpression("abs(genWeight)");
    dataloader->SetBackgroundWeightExpression("abs(genWeight)");

    // ------------------------------------------------------
    // VARIABLES
    // ------------------------------------------------------
    TString dimuonMassExpr =
"sqrt(2*Alt$(Muon_pt[0],0)*Alt$(Muon_pt[1],0)*"
"(cosh(Alt$(Muon_eta[0],0)-Alt$(Muon_eta[1],0))"
"-cos(Alt$(Muon_phi[0],0)-Alt$(Muon_phi[1],0))))";

dataloader->AddVariable(dimuonMassExpr, "dimuon_mass", "GeV", 'F');
    
    dataloader->AddVariable("Alt$(Muon_pt[0],0)", "muon1_pt", "GeV", 'F');
    dataloader->AddVariable("Alt$(Muon_pt[1],0)", "muon2_pt", "GeV", 'F');
    dataloader->AddVariable("Alt$(Muon_eta[0],0)", "muon1_eta", "", 'F');
    dataloader->AddVariable("Alt$(Muon_eta[1],0)", "muon2_eta", "", 'F');
    dataloader->AddVariable(
        "sqrt(pow(Alt$(Muon_eta[0],0)-Alt$(Muon_eta[1],0),2)+pow(acos(cos(Alt$(Muon_phi[0],0)-Alt$(Muon_phi[1],0))),2))",
        "deltaR_mumu", "", 'F'
    );
    
    dataloader->AddVariable("Alt$(PuppiMET_pt,0)",  "met_pt",  "GeV", 'F');
    dataloader->AddVariable("Alt$(PuppiMET_phi,0)", "met_phi", "", 'F');
    
    // Leading jets
dataloader->AddVariable("Alt$(Jet_pt[0],0)", "jet1_pt", "GeV", 'F');
dataloader->AddVariable("Alt$(Jet_pt[1],0)", "jet2_pt", "GeV", 'F');

dataloader->AddVariable(
    "abs(Alt$(Jet_eta[0],0)-Alt$(Jet_eta[1],0))",
    "detajj", "", 'F'
);

TString mjjExpr =
"sqrt(2*Alt$(Jet_pt[0],0)*Alt$(Jet_pt[1],0)*"
"(cosh(Alt$(Jet_eta[0],0)-Alt$(Jet_eta[1],0))"
"-cos(Alt$(Jet_phi[0],0)-Alt$(Jet_phi[1],0))))";

dataloader->AddVariable(mjjExpr, "mjj", "GeV", 'F');


dataloader->AddVariable("nJet", "nJets", "", 'I');


    // ------------------------------------------------------
    // PRESELECTION
    // ------------------------------------------------------
TCut preselection =
TCut("(nMuon>=2)"
     "&&(nElectron==0)"
     "&&(Alt$(Muon_charge[0],0)!=Alt$(Muon_charge[1],0))"
     "&&(Alt$(Muon_pt[0],0)>26)"
     "&&(Alt$(Muon_pt[1],0)>15)"
     "&&(abs(Alt$(Muon_eta[0],0))<2.4)"
     "&&(abs(Alt$(Muon_eta[1],0))<2.4)"
     "&&(Sum$(Jet_pt > 25 && abs(Jet_eta) < 4.7) >= 2)"
     "&&(Alt$(Jet_pt[0],0)>35)");

preselection = preselection &&
               TCut(Form("(%s>110)&&(%s<150)",
                         dimuonMassExpr.Data(),
                         dimuonMassExpr.Data()));

    
   dataloader->PrepareTrainingAndTestTree(
    preselection,
    "nTrain_Signal=0:"
    "nTrain_Background=0:"
    "SplitMode=Random:"
    "NormMode=NumEvents:"
    "!V"
);

    // ------------------------------------------------------
    // DNN METHOD (TMVA DL)
    // ------------------------------------------------------
    // ------------------------------------------------------
// DNN METHOD (TMVA DL)
// ------------------------------------------------------
factory.BookMethod(
    dataloader,
    TMVA::Types::kDL,
    "DL",
    "!H:V:"
    "VarTransform=N:" //burasi G`de olabilir, Gaussian Transformations icin
    "ErrorStrategy=CROSSENTROPY:"
    "TrainingStrategy="
        "LearningRate=0.001,"
        "Momentum=0.9,"
        "MaxEpochs=100,"
        "BatchSize=1000:"
    "Layout=RELU|225,RELU|100,RELU|64,SIGMOID"
);


    factory.TrainAllMethods();
    factory.TestAllMethods();
    factory.EvaluateAllMethods();

    outputFile->Close();

    std::cout << "\n=== Training finished ===\n";

    // ======================================================
    // 🔍 ANALYSIS PART (SAFE + CORRECT PATHS)
    // ======================================================
   TFile* f = TFile::Open("TMVA_HTo2Mu_DNN.root");
if (!f || f->IsZombie()) {
    std::cerr << "ERROR: cannot open output file\n";
    return;
}
    f->ls();

    // ROOT directory check
    TDirectory* dir = (TDirectory*)f->Get("dataset");
    if (!dir) {
        std::cerr << "ERROR: dataset directory not found\n";
        return;
    }

    dir->ls();

    // ------------------------------------------------------
    // SAFE GET HELPERS
    // ------------------------------------------------------
    auto getH = [&](const char* name) -> TH1* {
        TH1* h = (TH1*)f->Get(name);
        if (!h) std::cerr << "WARNING: missing histogram " << name << std::endl;
        return h;
    };

    // ⚠️ TMVA DL correct naming (important fix)
    TH1* hTrainSig = getH("dataset/Method_DL/DL/MVA_DL_Train_S");
    TH1* hTrainBkg = getH("dataset/Method_DL/DL/MVA_DL_Train_B");
    TH1* hTestSig  = getH("dataset/Method_DL/DL/MVA_DL_S");
    TH1* hTestBkg  = getH("dataset/Method_DL/DL/MVA_DL_B");

    if (!hTrainSig || !hTrainBkg || !hTestSig || !hTestBkg) {
        std::cerr << "ERROR: missing TMVA histograms -> check file structure\n";
        return;
    }

    // ------------------------------------------------------
    // DRAW SCORE DISTRIBUTIONS
    // ------------------------------------------------------
    // ------------------------------------------------------
// DRAW SCORE DISTRIBUTIONS
// ------------------------------------------------------
TCanvas* c1 = new TCanvas("c1", "DNN Score", 900, 700);

hTrainSig->SetLineColor(kBlue);
hTrainBkg->SetLineColor(kRed);
hTestSig->SetLineColor(kBlue+1);
hTestBkg->SetLineColor(kRed+1);

hTrainSig->SetLineWidth(2);
hTrainBkg->SetLineWidth(2);
hTestSig->SetLineWidth(2);
hTestBkg->SetLineWidth(2);

// force visible range
hTrainSig->GetXaxis()->SetRangeUser(0.0,1.0);

hTrainSig->SetTitle("DNN Response");
hTrainSig->GetXaxis()->SetTitle("DNN Score");
hTrainSig->GetYaxis()->SetTitle("Events");

hTrainSig->Draw("HIST");
hTrainBkg->Draw("HIST SAME");
hTestSig->Draw("E SAME");
hTestBkg->Draw("E SAME");

TLegend *leg = new TLegend(0.60,0.70,0.88,0.88);
leg->AddEntry(hTrainSig,"Signal (train)","l");
leg->AddEntry(hTrainBkg,"Background (train)","l");
leg->AddEntry(hTestSig ,"Signal (test)","lep");
leg->AddEntry(hTestBkg ,"Background (test)","lep");
leg->Draw();

c1->SaveAs("DNN_score.png");

    // ------------------------------------------------------
    // ROC CURVE (SAFE)
    // ------------------------------------------------------
    TObject* objROC = f->Get("dataset/Method_DL/DL/MVA_DL_trainingEffBvsRejBf");

    if(objROC){
        TCanvas* c2 = new TCanvas("c2", "ROC Curve", 800, 600);

        if(objROC->InheritsFrom("TGraph")){
            TGraph* rocGraph = (TGraph*)objROC;
            rocGraph->SetLineWidth(3);
            rocGraph->SetLineColor(kBlack);
            rocGraph->Draw("AL");
            c2->SaveAs("ROC_curve_DNN.png");
            std::cout << "ROC curve saved: ROC_curve.png\n";

            double auc = 0.0;
            for(int i=1; i<rocGraph->GetN(); i++){
                double x1,y1,x2,y2;
                rocGraph->GetPoint(i-1,x1,y1);
                rocGraph->GetPoint(i,x2,y2);
                auc += 0.5*(y1+y2)*fabs(x2-x1);
            }
            std::cout << "AUC = " << auc << std::endl;

        } else if(objROC->InheritsFrom("TH1")){
            TH1* hROC = (TH1*)objROC;
            hROC->SetLineWidth(3);
            hROC->SetLineColor(kBlack);
            hROC->Draw("HIST");
            c2->SaveAs("ROC_curve_DNN.png");
            std::cout << "ROC curve saved: ROC_curve.png\n";

            int nBins = hROC->GetNbinsX();
            double auc = 0.0;
            for(int i=1; i<nBins; i++){
                double x1 = hROC->GetBinLowEdge(i);
                double x2 = hROC->GetBinLowEdge(i+1);
                double y1 = hROC->GetBinContent(i);
                double y2 = hROC->GetBinContent(i+1);
                auc += 0.5*(y1+y2)*fabs(x2-x1);
            }
            std::cout << "AUC (approx.) = " << auc << std::endl;
        } else {
            std::cout << "ROC object is not TGraph or TH1, cannot compute AUC.\n";
        }
    } else {
        std::cout << "WARNING: ROC curve object not found.\n";
    }


    TMVA::TMVAGui("TMVA_HTo2Mu_DNN.root");
}
