#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TCut.h"

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"

#include "TH1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"

#include <iostream>
#include <vector>
#include <string>


void train_DNN_ggH_HTo2Mu_2025()
{
    // ========================================================================
    // INITIALIZATION
    // ========================================================================

    TMVA::Tools::Instance();

    TString outputName =
        "TMVA_ggH_HTo2Mu_DNN_2025.root";


    TFile *outputFile =
        TFile::Open(
            outputName,
            "RECREATE"
        );


    if(!outputFile || outputFile->IsZombie())
    {
        std::cerr
            << "ERROR: Cannot create TMVA output file."
            << std::endl;

        return;
    }


    // ========================================================================
    // TMVA FACTORY
    // ========================================================================

    TMVA::Factory factory(
        "TMVAClassification",
        outputFile,
        "!V:"
        "Color:"
        "DrawProgressBar:"
        "AnalysisType=Classification"
    );


    TMVA::DataLoader *dataloader =
        new TMVA::DataLoader(
            "dataset"
        );


    // ========================================================================
    // SIGNAL
    //
    // ggH -> H -> mumu
    //
    // ========================================================================

    TString signalFileName =
        "TMVA_2025_ggH_signal.root";


    TFile *signalFile =
        TFile::Open(
            signalFileName
        );


    if(!signalFile || signalFile->IsZombie())
    {
        std::cerr
            << "ERROR: Cannot open signal file: "
            << signalFileName
            << std::endl;

        outputFile->Close();

        return;
    }


    TTree *signalTree =
        dynamic_cast<TTree*>(
            signalFile->Get(
                "TMVATree"
            )
        );


    if(!signalTree)
    {
        std::cerr
            << "ERROR: TMVATree not found in signal file."
            << std::endl;

        outputFile->Close();

        return;
    }


    std::cout
        << "\n============================================"
        << "\n SIGNAL"
        << "\n============================================"
        << "\n Process  : ggH -> H -> mu mu"
        << "\n File     : "
        << signalFileName
        << "\n Entries  : "
        << signalTree->GetEntries()
        << "\n"
        << std::endl;


    dataloader->AddSignalTree(
        signalTree,
        1.0
    );


    // ========================================================================
    // BACKGROUND FILES
    //
    // DY
    // EWK
    // VV
    // TTbar
    //
    // Each process is stored separately.
    // ========================================================================

    std::vector<TString> backgroundFiles =
    {
        "TMVA_2025_DY_background.root",

        "TMVA_2025_EWK_background.root",

        "TMVA_2025_VV_background.root",

        "TMVA_2025_TTbar_background.root"
    };


    std::vector<TFile*> backgroundROOTFiles;


    std::cout
        << "\n============================================"
        << "\n BACKGROUND"
        << "\n============================================"
        << std::endl;


    for(const auto &fileName : backgroundFiles)
    {
        TFile *file =
            TFile::Open(
                fileName
            );


        if(!file || file->IsZombie())
        {
            std::cerr
                << "ERROR: Cannot open background file: "
                << fileName
                << std::endl;

            outputFile->Close();

            return;
        }


        TTree *tree =
            dynamic_cast<TTree*>(
                file->Get(
                    "TMVATree"
                )
            );


        if(!tree)
        {
            std::cerr
                << "ERROR: TMVATree not found in: "
                << fileName
                << std::endl;

            outputFile->Close();

            return;
        }


        std::cout
            << "  "
            << fileName
            << " : "
            << tree->GetEntries()
            << " entries"
            << std::endl;


        dataloader->AddBackgroundTree(
            tree,
            1.0
        );


        backgroundROOTFiles.push_back(
            file
        );
    }


    // ========================================================================
    // EVENT WEIGHT
    //
    // eventWeight is assumed to already contain:
    //
    // genWeight
    // xsec * lumi / sumGenWeight
    // PU weight
    // muon SF
    // b-tag SF
    //
    // Therefore we use eventWeight directly.
    // ========================================================================

    dataloader->SetSignalWeightExpression(
        "eventWeight"
    );


    dataloader->SetBackgroundWeightExpression(
        "eventWeight"
    );


    // ========================================================================
    // INPUT VARIABLES
    //
    // IMPORTANT:
    //
    // m_mumu is NOT included.
    //
    // This prevents the network from learning the Higgs mass peak directly.
    //
    // 13 variables:
    //
    //  1  mu1_pt
    //  2  mu2_pt
    //  3  mu1_eta
    //  4  mu2_eta
    //  5  dR_mumu
    //  6  dimuon_pt
    //  7  dimuon_eta
    //  8  met_pt
    //  9  nJet
    // 10  jet1_pt
    // 11  jet2_pt
    // 12  dEta_jj
    // 13  mjj
    //
    // ========================================================================


    // ------------------------------------------------------------------------
    // MUONS
    // ------------------------------------------------------------------------

    dataloader->AddVariable(
        "mu1_pt",
        "p_{T}^{#mu1}",
        "GeV",
        'F'
    );


    dataloader->AddVariable(
        "mu2_pt",
        "p_{T}^{#mu2}",
        "GeV",
        'F'
    );


    dataloader->AddVariable(
        "mu1_eta",
        "#eta_{#mu1}",
        "",
        'F'
    );


    dataloader->AddVariable(
        "mu2_eta",
        "#eta_{#mu2}",
        "",
        'F'
    );


    dataloader->AddVariable(
        "dR_mumu",
        "#Delta R_{#mu#mu}",
        "",
        'F'
    );


    dataloader->AddVariable(
        "dimuon_pt",
        "p_{T}^{#mu#mu}",
        "GeV",
        'F'
    );


    dataloader->AddVariable(
        "dimuon_eta",
        "#eta_{#mu#mu}",
        "",
        'F'
    );


    // ------------------------------------------------------------------------
    // MET
    // ------------------------------------------------------------------------

    dataloader->AddVariable(
        "met_pt",
        "p_{T}^{miss}",
        "GeV",
        'F'
    );


    // ------------------------------------------------------------------------
    // JETS
    // ------------------------------------------------------------------------

    dataloader->AddVariable(
        "nJet",
        "N_{jets}",
        "",
        'F'
    );


    dataloader->AddVariable(
        "jet1_pt",
        "p_{T}^{j1}",
        "GeV",
        'F'
    );


    dataloader->AddVariable(
        "jet2_pt",
        "p_{T}^{j2}",
        "GeV",
        'F'
    );


    dataloader->AddVariable(
        "dEta_jj",
        "|#Delta#eta_{jj}|",
        "",
        'F'
    );


    dataloader->AddVariable(
        "mjj",
        "m_{jj}",
        "GeV",
        'F'
    );


    // ========================================================================
    // PRESELECTION
    //
    // TMVATree is assumed to already contain events after:
    //
    // nMuon == 2
    // opposite sign
    // |eta| < 2.4
    // Medium ID
    // isolation < 0.25
    // leading pT > 26 GeV
    // subleading pT > 20 GeV
    // HLT_IsoMu24
    // trigger matching
    //
    // ========================================================================

    TCut preselection =
        "mu1_pt > 0"
        " && mu2_pt > 0";


    // ========================================================================
    // TRAIN / TEST SPLIT
    // ========================================================================

    dataloader->PrepareTrainingAndTestTree(
        preselection,
        "nTrain_Signal=0:"
        "nTrain_Background=0:"
        "SplitMode=Random:"
        "NormMode=NumEvents:"
        "!V"
    );


    // ========================================================================
    // DNN CONFIGURATION
    //
    // Architecture:
    //
    // 13
    //  |
    // 225
    //  |
    // 100
    //  |
    // 64
    //  |
    // 1
    //
    // Activation:
    //   ReLU hidden layers
    //
    // Output:
    //   Linear + cross entropy
    //
    // Training:
    //   Learning rate = 0.001
    //   Momentum      = 0.9
    //   Batch size    = 1000
    //   Epochs        = 100
    //
    // ========================================================================

    factory.BookMethod(
        dataloader,
        TMVA::Types::kDL,
        "DNN",
        "!H:"
        "!V:"
        "ErrorStrategy=CROSSENTROPY:"
        "VarTransform=N:"
        "WeightInitialization=XAVIERUNIFORM:"
        "Layout=RELU|225|RELU|100|RELU|64|LINEAR:"
        "TrainingStrategy="
        "LearningRate=1e-3,"
        "Momentum=0.9,"
        "Repetitions=1,"
        "ConvergenceSteps=10,"
        "BatchSize=1000,"
        "TestRepetitions=1,"
        "WeightDecay=1e-4,"
        "MaxEpochs=100"
    );


    // ========================================================================
    // TRAIN
    // ========================================================================

    std::cout
        << "\n"
        << "===================================================="
        << "\n"
        << " TRAINING ggH -> H -> mu mu DNN"
        << "\n"
        << "===================================================="
        << std::endl;


    factory.TrainAllMethods();


    // ========================================================================
    // TEST
    // ========================================================================

    std::cout
        << "\n"
        << "===================================================="
        << "\n"
        << " TESTING"
        << "\n"
        << "===================================================="
        << std::endl;


    factory.TestAllMethods();


    // ========================================================================
    // EVALUATE
    // ========================================================================

    std::cout
        << "\n"
        << "===================================================="
        << "\n"
        << " EVALUATING"
        << "\n"
        << "===================================================="
        << std::endl;


    factory.EvaluateAllMethods();


    // ========================================================================
    // SAVE
    // ========================================================================

    outputFile->Write();
    outputFile->Close();


    std::cout
        << "\n"
        << "===================================================="
        << "\n"
        << " TMVA DNN TRAINING COMPLETE"
        << "\n"
        << "===================================================="
        << "\n"
        << "Signal:"
        << "\n"
        << "  ggH -> H -> mu mu"
        << "\n"
        << "\n"
        << "Background:"
        << "\n"
        << "  DY"
        << "\n"
        << "  EWK"
        << "\n"
        << "  VV"
        << "\n"
        << "  TTbar"
        << "\n"
        << "\n"
        << "Architecture:"
        << "\n"
        << "  13 -> 225 -> 100 -> 64 -> 1"
        << "\n"
        << "\n"
        << "Training:"
        << "\n"
        << "  LR        = 0.001"
        << "\n"
        << "  Momentum  = 0.9"
        << "\n"
        << "  BatchSize = 1000"
        << "\n"
        << "  Epochs    = 100"
        << "\n"
        << "\n"
        << "Output:"
        << "\n"
        << "  "
        << outputName
        << "\n"
        << "===================================================="
        << std::endl;


    // ========================================================================
    // REOPEN TMVA OUTPUT
    // ========================================================================

    TFile *f =
        TFile::Open(
            outputName
        );


    if(!f || f->IsZombie())
    {
        std::cerr
            << "ERROR: Cannot reopen TMVA output."
            << std::endl;

        return;
    }


    // ========================================================================
    // DNN RESPONSE HISTOGRAMS
    // ========================================================================

    TH1 *hTrainSig =
        dynamic_cast<TH1*>(
            f->Get(
                "dataset/Method_DL/DNN/"
                "MVA_DNN_Train_S"
            )
        );


    TH1 *hTrainBkg =
        dynamic_cast<TH1*>(
            f->Get(
                "dataset/Method_DL/DNN/"
                "MVA_DNN_Train_B"
            )
        );


    TH1 *hTestSig =
        dynamic_cast<TH1*>(
            f->Get(
                "dataset/Method_DL/DNN/"
                "MVA_DNN_S"
            )
        );


    TH1 *hTestBkg =
        dynamic_cast<TH1*>(
            f->Get(
                "dataset/Method_DL/DNN/"
                "MVA_DNN_B"
            )
        );


    if(
        hTrainSig &&
        hTrainBkg &&
        hTestSig &&
        hTestBkg
    )
    {
        TCanvas *c1 =
            new TCanvas(
                "c1",
                "DNN response",
                900,
                700
            );


        hTrainSig->SetLineColor(kBlue);
        hTrainBkg->SetLineColor(kRed);

        hTestSig->SetLineColor(kBlue+1);
        hTestBkg->SetLineColor(kRed+1);


        hTrainSig->SetLineWidth(2);
        hTrainBkg->SetLineWidth(2);

        hTestSig->SetLineWidth(2);
        hTestBkg->SetLineWidth(2);


        hTrainSig->SetTitle(
            "ggH #rightarrow #mu#mu DNN response"
        );


        hTrainSig->GetXaxis()->SetTitle(
            "DNN response"
        );


        hTrainSig->GetYaxis()->SetTitle(
            "Events"
        );


        hTrainSig->Draw("HIST");

        hTrainBkg->Draw(
            "HIST SAME"
        );

        hTestSig->Draw(
            "HIST SAME"
        );

        hTestBkg->Draw(
            "HIST SAME"
        );


        TLegend *legend =
            new TLegend(
                0.55,
                0.60,
                0.88,
                0.88
            );


        legend->AddEntry(
            hTrainSig,
            "ggH signal - Train",
            "l"
        );


        legend->AddEntry(
            hTrainBkg,
            "DY+EWK+VV+TTbar - Train",
            "l"
        );


        legend->AddEntry(
            hTestSig,
            "ggH signal - Test",
            "l"
        );


        legend->AddEntry(
            hTestBkg,
            "DY+EWK+VV+TTbar - Test",
            "l"
        );


        legend->Draw();


        c1->SaveAs(
            "DNN_score_ggH_2025.png"
        );


        std::cout
            << "DNN score plot saved: "
            << "DNN_score_ggH_2025.png"
            << std::endl;
    }
    else
    {
        std::cerr
            << "WARNING: Could not find all DNN response histograms."
            << std::endl;
    }


    // ========================================================================
    // ROC CURVE
    // ========================================================================

    TObject *objROC =
        f->Get(
            "dataset/Method_DL/DNN/"
            "MVA_DNN_TrainingROC"
        );


    if(
        objROC &&
        objROC->InheritsFrom("TGraph")
    )
    {
        TCanvas *c2 =
            new TCanvas(
                "c2",
                "DNN ROC curve",
                800,
                600
            );


        TGraph *rocGraph =
            dynamic_cast<TGraph*>(
                objROC
            );


        rocGraph->SetLineWidth(3);


        rocGraph->SetTitle(
            "ggH #rightarrow #mu#mu DNN ROC"
        );


        rocGraph->GetXaxis()->SetTitle(
            "Background efficiency"
        );


        rocGraph->GetYaxis()->SetTitle(
            "Signal efficiency"
        );


        rocGraph->Draw(
            "AL"
        );


        c2->SaveAs(
            "ROC_curve_ggH_DNN_2025.png"
        );


        // --------------------------------------------------------------------
        // AUC
        // --------------------------------------------------------------------

        double auc = 0.0;


        for(
            int i=1;
            i<rocGraph->GetN();
            ++i
        )
        {
            double x1;
            double y1;

            double x2;
            double y2;


            rocGraph->GetPoint(
                i-1,
                x1,
                y1
            );


            rocGraph->GetPoint(
                i,
                x2,
                y2
            );


            auc +=
                0.5 *
                (y1+y2) *
                (x2-x1);
        }


        std::cout
            << "\n"
            << "===================================================="
            << "\n"
            << " DNN ROC / AUC"
            << "\n"
            << "===================================================="
            << "\n"
            << "AUC = "
            << auc
            << "\n"
            << std::endl;
    }
    else
    {
        std::cerr
            << "WARNING: DNN ROC curve not found."
            << std::endl;
    }


    // ========================================================================
    // CLOSE INPUT FILES
    // ========================================================================

    signalFile->Close();


    for(auto file : backgroundROOTFiles)
    {
        if(file)
            file->Close();
    }


    f->Close();


    // ========================================================================
    // TMVA GUI
    // ========================================================================

    TMVA::TMVAGui(
        outputName
    );
}
