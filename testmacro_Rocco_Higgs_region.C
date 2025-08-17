#include "RoccoR.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TLorentzVector.h"

struct Muon {
    int charge;
    double pt, eta, phi;
    double genPt;
    bool tightId;
    bool isGlobal;
    bool isTracker;
    bool isStandalone;
    double iso04;
    UChar_t nStations;
};

double invariantMass(double pt1, double eta1, double phi1,
                     double pt2, double eta2, double phi2) {
    return sqrt(2 * pt1 * pt2 * (cosh(eta1 - eta2) - cos(phi1 - phi2)));
}

void testmacro_Rocco_Higgs_region(const std::vector<std::string> &inputFiles, bool isMC = false) {

    RoccoR rc;
    rc.init(edm::FileInPath("RoccoR/data/RoccoR2022.txt").fullPath());

    // Histograms for Higgs region
    TH1D* hMassUncorrected = new TH1D("hMassUncorrected", "Dimuon Mass Uncorrected (110-150 GeV);Mass (GeV);Events", 80, 110, 150);
    TH1D* hMassCorrected   = new TH1D("hMassCorrected",   "Dimuon Mass Rochester-Corrected (110-150 GeV);Mass (GeV);Events", 80, 110, 150);

    // Variables
    Int_t   nMuon;
    Float_t Muon_pt[100], Muon_eta[100], Muon_phi[100], Muon_pfRelIso04_all[100];
    Int_t   Muon_charge[100];
    Bool_t  Muon_tightId[100], Muon_isGlobal[100], Muon_isTracker[100], Muon_isStandalone[100];
    UChar_t Muon_nStations[100];
    Float_t Muon_genPt[100]; // only for MC if exists

    Bool_t HLT_IsoMu24;

    Int_t   nTrigObj;
    Float_t TrigObj_pt[100], TrigObj_eta[100], TrigObj_phi[100];
    UShort_t TrigObj_id[100];

    // Build chain
    TChain chain("Events");
    for (const auto &fname : inputFiles) {
        chain.Add(fname.c_str());
    }

    // Check if Muon_genPt branch exists
    bool hasMuonGenPt = (chain.GetBranch("Muon_genPt") != nullptr);

    // Set branches
    chain.SetBranchAddress("nMuon", &nMuon);
    chain.SetBranchAddress("Muon_pt", Muon_pt);
    chain.SetBranchAddress("Muon_eta", Muon_eta);
    chain.SetBranchAddress("Muon_phi", Muon_phi);
    chain.SetBranchAddress("Muon_charge", Muon_charge);
    chain.SetBranchAddress("Muon_tightId", Muon_tightId);
    chain.SetBranchAddress("Muon_isGlobal", Muon_isGlobal);
    chain.SetBranchAddress("Muon_isTracker", Muon_isTracker);
    chain.SetBranchAddress("Muon_isStandalone", Muon_isStandalone);
    chain.SetBranchAddress("Muon_nStations", Muon_nStations);
    chain.SetBranchAddress("Muon_pfRelIso04_all", Muon_pfRelIso04_all);

    if (isMC && hasMuonGenPt) {
        chain.SetBranchAddress("Muon_genPt", Muon_genPt);
    } else if (isMC && !hasMuonGenPt) {
        std::cerr << "Warning: Muon_genPt branch missing in MC file. Using kScaleMC only." << std::endl;
    }

    chain.SetBranchAddress("HLT_IsoMu24", &HLT_IsoMu24);

    chain.SetBranchAddress("nTrigObj", &nTrigObj);
    chain.SetBranchAddress("TrigObj_pt", TrigObj_pt);
    chain.SetBranchAddress("TrigObj_eta", TrigObj_eta);
    chain.SetBranchAddress("TrigObj_phi", TrigObj_phi);
    chain.SetBranchAddress("TrigObj_id", TrigObj_id);

    Long64_t nEntries = chain.GetEntries();
    std::cout << "Processing " << nEntries << " events from " << inputFiles.size() << " files..." << std::endl;

    for (Long64_t i = 0; i < nEntries; ++i) {
        chain.GetEntry(i);

        // Trigger and at least 2 muons
        if (!HLT_IsoMu24 || nMuon < 2) continue;

        int tagIdx = -1;
        float maxTagPt = -1;

        // Tag muon selection
        for (int m = 0; m < nMuon; ++m) {
            if (!Muon_tightId[m]) continue;
            if (Muon_pt[m] < 29 || fabs(Muon_eta[m]) >= 2.4 || Muon_pfRelIso04_all[m] >= 0.15) continue;
            if (!Muon_isGlobal[m] || !Muon_isTracker[m]) continue;

            bool matched = false;
            for (int t = 0; t < nTrigObj; ++t) {
                if (TrigObj_id[t] != 13 || TrigObj_pt[t] < 24) continue;
                float deta = Muon_eta[m] - TrigObj_eta[t];
                float dphi = fabs(Muon_phi[m] - TrigObj_phi[t]);
                if (dphi > M_PI) dphi = 2 * M_PI - dphi;
                float dr2 = deta*deta + dphi*dphi;
                if (dr2 < 0.09) matched = true;
            }

            if (matched && Muon_pt[m] > maxTagPt) {
                maxTagPt = Muon_pt[m];
                tagIdx = m;
            }
        }

        if (tagIdx < 0) continue;

        // Probe selection
        for (int p = 0; p < nMuon; ++p) {
            if (p == tagIdx) continue;
            if (Muon_charge[tagIdx] * Muon_charge[p] >= 0) continue;
            if (!Muon_isTracker[p] || !Muon_isStandalone[p] || Muon_pt[p] < 20 || fabs(Muon_eta[p]) >= 2.4 || Muon_nStations[p] <= 1) continue;

            TLorentzVector tag, probe;
            tag.SetPtEtaPhiM(Muon_pt[tagIdx], Muon_eta[tagIdx], Muon_phi[tagIdx], 0.105);
            probe.SetPtEtaPhiM(Muon_pt[p], Muon_eta[p], Muon_phi[p], 0.105);
            float massUncorr = (tag + probe).M();

            // Higgs mass window
            if (massUncorr < 110 || massUncorr > 150) continue;

            // Rochester correction
            double sf_tag = 1.0, sf_probe = 1.0;
            if (isMC) {
                if (hasMuonGenPt && Muon_genPt[tagIdx] > 0)
                    sf_tag = rc.kSpreadMC(Muon_charge[tagIdx], Muon_pt[tagIdx], Muon_eta[tagIdx], Muon_phi[tagIdx], Muon_genPt[tagIdx], 0, 0);
                else
                    sf_tag = rc.kScaleMC(Muon_charge[tagIdx], Muon_pt[tagIdx], Muon_eta[tagIdx], Muon_phi[tagIdx], 0, 0);

                if (hasMuonGenPt && Muon_genPt[p] > 0)
                    sf_probe = rc.kSpreadMC(Muon_charge[p], Muon_pt[p], Muon_eta[p], Muon_phi[p], Muon_genPt[p], 0, 0);
                else
                    sf_probe = rc.kScaleMC(Muon_charge[p], Muon_pt[p], Muon_eta[p], Muon_phi[p], 0, 0);
            } else {
                sf_tag = rc.kScaleDT(Muon_charge[tagIdx], Muon_pt[tagIdx], Muon_eta[tagIdx], Muon_phi[tagIdx], 0, 0);
                sf_probe = rc.kScaleDT(Muon_charge[p], Muon_pt[p], Muon_eta[p], Muon_phi[p], 0, 0);
            }

            TLorentzVector tagCorr, probeCorr;
            tagCorr.SetPtEtaPhiM(Muon_pt[tagIdx]*sf_tag, Muon_eta[tagIdx], Muon_phi[tagIdx], 0.105);
            probeCorr.SetPtEtaPhiM(Muon_pt[p]*sf_probe, Muon_eta[p], Muon_phi[p], 0.105);
            float massCorr = (tagCorr + probeCorr).M();

            // Fill histograms
            hMassUncorrected->Fill(massUncorr);
            hMassCorrected->Fill(massCorr);
        }
    }

    TFile outFile(isMC ? "Dimuon_TagProbe_MC_Higgs.root" : "Dimuon_TagProbe_Data_Higgs.root", "RECREATE");
    hMassUncorrected->Write();
    hMassCorrected->Write();
    outFile.Close();

    std::cout << "Histograms saved to " << (isMC ? "Dimuon_TagProbe_MC_Higgs.root" : "Dimuon_TagProbe_Data_Higgs.root") << std::endl;
}
