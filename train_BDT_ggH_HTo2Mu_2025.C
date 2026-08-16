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


void train_BDT_ggH_HTo2Mu_2025()
{
    // ========================================================================
    // INITIALIZATION
    // ========================================================================

    TMVA::Tools::Instance();

    TString outputName =
        "TMVA_ggH_HTo2Mu_BDT_2025.root";

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
    // One separate signal file.
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
            << "ERROR: TMVATree not found in signal file: "
            << signalFileName
            << std::endl;

        outputFile->Close();

        return;
    }


    std::cout
        << "\nSignal:"
        << "\n  ggH"
        << "\n  File = "
        << signalFileName
        << "\n  Entries = "
        << signalTree->GetEntries()
        << std::endl;


    dataloader->AddSignalTree(
        signalTree,
        1.0
    );


    // ========================================================================
    // BACKGROUND FILES
    //
    // Each process is kept in a separate ROOT file.
    //
    // DY
    // EWK
    // VV = WZ + ZZ
    // TTbar
    //
    // TMVA can directly use several background trees.
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
        << "\nBackground samples:"
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


        // ------------------------------------------------------------
        // Add this process as a background tree.
        //
        // eventWeight is already stored in the tree.
        // ------------------------------------------------------------

        dataloader->AddBackgroundTree(
            tree,
            1.0
        );


        backgroundROOTFiles.push_back(
            file
        );
    }


    // ========================================================================
    // EVENT WEIGHTS
    //
    // eventWeight should already contain:
    //
    //   genWeight
    //   xsec * lumi / sumGenWeight
    //   PU weight
    //   muon SF
    //   b-tag SF
    //
    // Therefore:
    //
    //     weight = eventWeight
    //
    // DO NOT use:
    //
    //     abs(genWeight)
    //
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
    // m_mumu is deliberately NOT included.
    //
    // The BDT should learn kinematic differences between:
    //
    //     ggH signal
    //
    // and
    //
    //     DY + EWK + VV + TTbar
    //
    // without learning the Higgs mass peak itself.
    // ========================================================================


    // ------------------------------------------------------------------------
    // MUON VARIABLES
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
    // TMVATree is assumed to have already been produced after:
    //
    //   nMuon == 2
    //   opposite sign
    //   |eta_mu| < 2.4
    //   Medium ID
    //   isolation < 0.25
    //   leading muon pT > 26 GeV
    //   subleading muon pT > 20 GeV
    //   HLT_IsoMu24
    //   trigger matching
    //
    // We only require physically valid variables here.
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
    // BDT CONFIGURATION
    //
    // Your established configuration:
    //
    // NTrees      = 500
    // MaxDepth    = 5
    // MinNodeSize = 3%
    // BoostType   = Grad
    // Shrinkage   = 0.10
    //
    // ========================================================================

    factory.BookMethod(
        dataloader,
        TMVA::Types::kBDT,
        "BDT",
        "!H:"
        "!V:"
        "NTrees=500:"
        "MinNodeSize=3%:"
        "MaxDepth=5:"
        "BoostType=Grad:"
        "Shrinkage=0.10:"
        "SeparationType=GiniIndex"
    );


    // ========================================================================
    // TRAIN
    // ========================================================================

    std::cout
        << "\n"
        << "===================================================="
        << "\n"
        << " TRAINING ggH -> H -> mu mu BDT"
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
        << " TMVA TRAINING COMPLETE"
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
    // BDT SCORE HISTOGRAMS
    // ========================================================================

    TH1 *hTrainSig =
        dynamic_cast<TH1*>(
            f->Get(
                "dataset/Method_BDT/BDT/"
                "MVA_BDT_Train_S"
            )
        );


    TH1 *hTrainBkg =
        dynamic_cast<TH1*>(
            f->Get(
                "dataset/Method_BDT/BDT/"
                "MVA_BDT_Train_B"
            )
        );


    TH1 *hTestSig =
        dynamic_cast<TH1*>(
            f->Get(
                "dataset/Method_BDT/BDT/"
                "MVA_BDT_S"
            )
        );


    TH1 *hTestBkg =
        dynamic_cast<TH1*>(
            f->Get(
                "dataset/Method_BDT/BDT/"
                "MVA_BDT_B"
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
                "BDT response",
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
            "ggH #rightarrow #mu#mu BDT response"
        );


        hTrainSig->GetXaxis()->SetTitle(
            "BDT response"
        );


        hTrainSig->GetYaxis()->SetTitle(
            "Events"
        );


        hTrainSig->Draw("HIST");
        hTrainBkg->Draw("HIST SAME");

        hTestSig->Draw("HIST SAME");
        hTestBkg->Draw("HIST SAME");


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
            "BDT_score_ggH_2025.png"
        );


        std::cout
            << "BDT score plot saved: "
            << "BDT_score_ggH_2025.png"
            << std::endl;
    }
    else
    {
        std::cerr
            << "WARNING: Could not find all BDT score histograms."
            << std::endl;
    }


    // ========================================================================
    // ROC CURVE
    // ========================================================================

    TObject *objROC =
        f->Get(
            "dataset/Method_BDT/BDT/"
            "MVA_BDT_TrainingROC"
        );


    if(
        objROC &&
        objROC->InheritsFrom("TGraph")
    )
    {
        TCanvas *c2 =
            new TCanvas(
                "c2",
                "ROC curve",
                800,
                600
            );


        TGraph *rocGraph =
            dynamic_cast<TGraph*>(
                objROC
            );


        rocGraph->SetLineWidth(3);


        rocGraph->SetTitle(
            "ggH #rightarrow #mu#mu BDT ROC"
        );


        rocGraph->GetXaxis()->SetTitle(
            "Background efficiency"
        );


        rocGraph->GetYaxis()->SetTitle(
            "Signal efficiency"
        );


        rocGraph->Draw("AL");


        c2->SaveAs(
            "ROC_curve_ggH_2025.png"
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
            << " ROC / AUC"
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
            << "WARNING: ROC curve not found."
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
