#include "RoccoR.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include "TH1D.h"
#include "TFile.h"
#include "TChain.h"

//run the code: 
//.L testRoccoR_NanoAOD.C+
//testmacro_Rocco_NanoAOD("root://xrootd-cms.infn.it//store/data/Run2022F/Muon/NANOAOD/*.root", false);
//(MC için ikinci parametreyi true yap)

struct Muon {
    int charge;
    double pt, eta, phi;
    double genPt;
};

double invariantMass(double pt1, double eta1, double phi1,
                     double pt2, double eta2, double phi2) {
    return sqrt(2 * pt1 * pt2 * (cosh(eta1 - eta2) - cos(phi1 - phi2)));
}

void testmacro_Rocco_Zregion(const std::string &inputFileName, bool isMC = false) {
    RoccoR rc;
    rc.init(edm::FileInPath("RoccoR/data/RoccoR2022.txt").fullPath());

    TH1D* hMassUncorrected = new TH1D("hMassUncorrected", "Dimuon Mass Uncorrected (70-110 GeV);Mass (GeV);Events", 80, 70, 110);
    TH1D* hMassCorrected   = new TH1D("hMassCorrected",   "Dimuon Mass Rochester-Corrected (70-110 GeV);Mass (GeV);Events", 80, 70, 110);

    // Open NanoAOD file
    TFile *file = TFile::Open(inputFileName.c_str());
    if (!file || file->IsZombie()) {
        std::cerr << "Cannot open file: " << inputFileName << std::endl;
        return;
    }

    TTree *tree = (TTree*)file->Get("Events");
    if (!tree) {
        std::cerr << "Cannot find Events tree in file" << std::endl;
        file->Close();
        return;
    }

    // Variables for fixed-size arrays and count
    Int_t nMuon = 0;
    Float_t Muon_pt[100];
    Float_t Muon_eta[100];
    Float_t Muon_phi[100];
    Int_t Muon_charge[100];
    Float_t Muon_genPt[100]; // only for MC

    // Set branch addresses
    tree->SetBranchAddress("nMuon", &nMuon);
    tree->SetBranchAddress("Muon_pt", Muon_pt);
    tree->SetBranchAddress("Muon_eta", Muon_eta);
    tree->SetBranchAddress("Muon_phi", Muon_phi);
    tree->SetBranchAddress("Muon_charge", Muon_charge);
    if (isMC) {
        tree->SetBranchAddress("Muon_genPt", Muon_genPt);
    }

    Long64_t nEntries = tree->GetEntries();
    std::cout << "Processing " << nEntries << " events..." << std::endl;

    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);

        if (nMuon < 2) continue;  // skip events with less than 2 muons

        // Fill muons vector
        std::vector<Muon> muons;
        for (int j = 0; j < nMuon; ++j) {
            Muon mu;
            mu.pt = Muon_pt[j];
            mu.eta = Muon_eta[j];
            mu.phi = Muon_phi[j];
            mu.charge = Muon_charge[j];
            mu.genPt = (isMC) ? Muon_genPt[j] : 0.0;
            muons.push_back(mu);
        }

        // Loop over all opposite-charge muon pairs
        for (size_t j = 0; j < muons.size(); ++j) {
            for (size_t k = j + 1; k < muons.size(); ++k) {
                if (muons[j].charge * muons[k].charge >= 0) continue; // skip same sign

                double massUncorr = invariantMass(muons[j].pt, muons[j].eta, muons[j].phi,
                                                 muons[k].pt, muons[k].eta, muons[k].phi);

                if (massUncorr < 70 || massUncorr > 110) continue; // Z window cut

                // Rochester correction
                double sf_j = 1.0, sf_k = 1.0;
                if (isMC) {
                    if (muons[j].genPt > 0)
                        sf_j = rc.kSpreadMC(muons[j].charge, muons[j].pt, muons[j].eta, muons[j].phi, muons[j].genPt, 0, 0);
                    else
                        sf_j = rc.kScaleMC(muons[j].charge, muons[j].pt, muons[j].eta, muons[j].phi, 0, 0);

                    if (muons[k].genPt > 0)
                        sf_k = rc.kSpreadMC(muons[k].charge, muons[k].pt, muons[k].eta, muons[k].phi, muons[k].genPt, 0, 0);
                    else
                        sf_k = rc.kScaleMC(muons[k].charge, muons[k].pt, muons[k].eta, muons[k].phi, 0, 0);
                } else {
                    sf_j = rc.kScaleDT(muons[j].charge, muons[j].pt, muons[j].eta, muons[j].phi, 0, 0);
                    sf_k = rc.kScaleDT(muons[k].charge, muons[k].pt, muons[k].eta, muons[k].phi, 0, 0);
                }

                double corrPt_j = muons[j].pt * sf_j;
                double corrPt_k = muons[k].pt * sf_k;

                double massCorr = invariantMass(corrPt_j, muons[j].eta, muons[j].phi,
                                               corrPt_k, muons[k].eta, muons[k].phi);

                hMassUncorrected->Fill(massUncorr);
                hMassCorrected->Fill(massCorr);
            }
        }
    }

    TFile outFile(isMC ? "DimuonMass_MC_Rochester.root" : "DimuonMass_Data_Rochester.root", "RECREATE");
    hMassUncorrected->Write();
    hMassCorrected->Write();
    outFile.Close();

    std::cout << "Histograms saved to " << (isMC ? "DimuonMass_MC_Rochester.root" : "DimuonMass_Data_Rochester.root") << std::endl;
}
