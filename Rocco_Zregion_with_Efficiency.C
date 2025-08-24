#include "RoccoR.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TLorentzVector.h"
#include "TEfficiency.h"

double invariantMass(double pt1, double eta1, double phi1,
                     double pt2, double eta2, double phi2) {
    return sqrt(2 * pt1 * pt2 * (cosh(eta1 - eta2) - cos(phi1 - phi2)));
}

void Rocco_Zregion_with_Efficiency(const std::vector<std::string> &inputFiles, bool isMC = false) {

    RoccoR rc;
    rc.init(edm::FileInPath("RoccoR/data/RoccoR2022.txt").fullPath());

    // --- 1D dimuon mass histograms ---
    TH1D* hMassUncorrected = new TH1D("hMassUncorrected", "Dimuon Mass Uncorrected (70-110 GeV);Mass (GeV);Events", 80, 70, 110);
    TH1D* hMassCorrected   = new TH1D("hMassCorrected",   "Dimuon Mass Rochester-Corrected (70-110 GeV);Mass (GeV);Events", 80, 70, 110);

    // --- Efficiency binning ---
    int nPtBins = 10;
    double ptBins[11] = {20, 25, 30, 35, 40, 50, 60, 80, 100, 120, 200};
    int nEtaBins = 8;
    double etaBins[9] = {0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.2, 2.4, 2.5};

    // --- 1D probe histograms ---
    TH1D* hProbe_pt_total = new TH1D("hProbe_pt_total", "Probe pt total;pt (GeV);# of probes", nPtBins, ptBins);
    TH1D* hProbe_pt_pass  = new TH1D("hProbe_pt_pass",  "Probe pt passing tightId;pt (GeV);# of probes", nPtBins, ptBins);
    TH1D* hProbe_eta_total = new TH1D("hProbe_eta_total", "Probe eta total;#eta;# of probes", nEtaBins, etaBins);
    TH1D* hProbe_eta_pass  = new TH1D("hProbe_eta_pass",  "Probe eta passing tightId;#eta;# of probes", nEtaBins, etaBins);

    // --- 2D probe histograms ---
    TH2D* hProbe_pteta_total = new TH2D("hProbe_pteta_total", "Probe pT vs eta total;pt (GeV);|eta|", nPtBins, ptBins, nEtaBins, etaBins);
    TH2D* hProbe_pteta_pass  = new TH2D("hProbe_pteta_pass",  "Probe pT vs eta passing tightId;pt (GeV);|eta|", nPtBins, ptBins, nEtaBins, etaBins);

    // --- NanoAOD variables ---
    Int_t nMuon;
    Float_t Muon_pt[100], Muon_eta[100], Muon_phi[100], Muon_pfRelIso04_all[100];
    Int_t Muon_charge[100];
    Bool_t Muon_tightId[100], Muon_isGlobal[100], Muon_isTracker[100], Muon_isStandalone[100];
    UChar_t Muon_nStations[100];
    Float_t Muon_genPt[100];
    Bool_t HLT_IsoMu24;
    Int_t nTrigObj;
    Float_t TrigObj_pt[100], TrigObj_eta[100], TrigObj_phi[100];
    UShort_t TrigObj_id[100];

    // --- Build TChain ---
    TChain chain("Events");
    for (const auto &fname : inputFiles) chain.Add(fname.c_str());
    bool hasMuonGenPt = (chain.GetBranch("Muon_genPt") != nullptr);

    // --- Set branches ---
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
    if (isMC && hasMuonGenPt) chain.SetBranchAddress("Muon_genPt", Muon_genPt);
    chain.SetBranchAddress("HLT_IsoMu24", &HLT_IsoMu24);
    chain.SetBranchAddress("nTrigObj", &nTrigObj);
    chain.SetBranchAddress("TrigObj_pt", TrigObj_pt);
    chain.SetBranchAddress("TrigObj_eta", TrigObj_eta);
    chain.SetBranchAddress("TrigObj_phi", TrigObj_phi);
    chain.SetBranchAddress("TrigObj_id", &TrigObj_id);

    Long64_t nEntries = chain.GetEntries();
    std::cout << "Processing " << nEntries << " events from " << inputFiles.size() << " files..." << std::endl;

    // --- Event loop ---
    for (Long64_t i = 0; i < nEntries; ++i) {
        chain.GetEntry(i);
        if (!HLT_IsoMu24 || nMuon < 2) continue;

        // --- Tag muon selection ---
        int tagIdx = -1;
        float maxTagPt = -1;
        for (int m = 0; m < nMuon; ++m) {
            if (!Muon_tightId[m]) continue;
            if (Muon_pt[m] < 29 || fabs(Muon_eta[m]) >= 2.4 || Muon_pfRelIso04_all[m] >= 0.15) continue;
            if (!Muon_isGlobal[m] || !Muon_isTracker[m]) continue;
            bool matched = false;
            for (int t = 0; t < nTrigObj; ++t) {
                if (TrigObj_id[t] != 13 || TrigObj_pt[t] < 24) continue;
                float deta = Muon_eta[m] - TrigObj_eta[t];
                float dphi = fabs(Muon_phi[m] - TrigObj_phi[t]);
                if (dphi > M_PI) dphi = 2*M_PI - dphi;
                if (deta*deta + dphi*dphi < 0.01) matched = true;
            }
            if (matched && Muon_pt[m] > maxTagPt) { maxTagPt = Muon_pt[m]; tagIdx = m; }
        }
        if (tagIdx < 0) continue;

        // --- Probe muons ---
        for (int p = 0; p < nMuon; ++p) {
            if (p == tagIdx) continue;
            if (Muon_charge[tagIdx]*Muon_charge[p] >= 0) continue;
            if (!Muon_isTracker[p] || !Muon_isStandalone[p] || Muon_nStations[p] <= 1) continue;

            // --- Rochester corrections ---
            double sf_tag=1.0, sf_probe=1.0;
            if (isMC) {
                if (hasMuonGenPt && Muon_genPt[tagIdx]>0)
                    sf_tag = rc.kSpreadMC(Muon_charge[tagIdx], Muon_pt[tagIdx], Muon_eta[tagIdx], Muon_phi[tagIdx], Muon_genPt[tagIdx],0,0);
                else
                    sf_tag = rc.kScaleMC(Muon_charge[tagIdx], Muon_pt[tagIdx], Muon_eta[tagIdx], Muon_phi[tagIdx],0,0);
                if (hasMuonGenPt && Muon_genPt[p]>0)
                    sf_probe = rc.kSpreadMC(Muon_charge[p], Muon_pt[p], Muon_eta[p], Muon_phi[p], Muon_genPt[p],0,0);
                else
                    sf_probe = rc.kScaleMC(Muon_charge[p], Muon_pt[p], Muon_eta[p], Muon_phi[p],0,0);
            } else {
                sf_tag = rc.kScaleDT(Muon_charge[tagIdx], Muon_pt[tagIdx], Muon_eta[tagIdx], Muon_phi[tagIdx],0,0);
                sf_probe = rc.kScaleDT(Muon_charge[p], Muon_pt[p], Muon_eta[p], Muon_phi[p],0,0);
            }

            // --- Corrected probe pT for MC efficiency ---
            double pt_probe_corr = Muon_pt[p];
            if (isMC) pt_probe_corr *= sf_probe;

            // --- Probe selection using corrected pT ---
            if (pt_probe_corr < 20 || fabs(Muon_eta[p]) >= 2.4) continue;

            TLorentzVector tag, probe, tagCorr, probeCorr;
            tag.SetPtEtaPhiM(Muon_pt[tagIdx], Muon_eta[tagIdx], Muon_phi[tagIdx], 0.105);
            probe.SetPtEtaPhiM(Muon_pt[p], Muon_eta[p], Muon_phi[p], 0.105);
            tagCorr.SetPtEtaPhiM(Muon_pt[tagIdx]*sf_tag, Muon_eta[tagIdx], Muon_phi[tagIdx], 0.105);
            probeCorr.SetPtEtaPhiM(pt_probe_corr, Muon_eta[p], Muon_phi[p], 0.105);

            float massUncorr = (tag+probe).M();
            float massCorr = (tagCorr + probeCorr).M();
            if (massUncorr < 70 || massUncorr > 115) continue;

            // --- Fill 1D & 2D histograms ---
            hProbe_pt_total->Fill(pt_probe_corr);
            hProbe_eta_total->Fill(fabs(Muon_eta[p]));
            hProbe_pteta_total->Fill(pt_probe_corr, fabs(Muon_eta[p]));

            hMassUncorrected->Fill(massUncorr);
            hMassCorrected->Fill(massCorr);

            if (Muon_tightId[p]) {
                hProbe_pt_pass->Fill(pt_probe_corr);
                hProbe_eta_pass->Fill(fabs(Muon_eta[p]));
                hProbe_pteta_pass->Fill(pt_probe_corr, fabs(Muon_eta[p]));
            }
        }
    }

    // --- TEfficiency objects ---
    TEfficiency* eff_pt  = new TEfficiency(*hProbe_pt_pass, *hProbe_pt_total);
    TEfficiency* eff_eta = new TEfficiency(*hProbe_eta_pass, *hProbe_eta_total);
    TEfficiency* eff_pteta = new TEfficiency(*hProbe_pteta_pass, *hProbe_pteta_total);

    // --- Save ROOT file ---
    TFile outFile(isMC ? "Dimuon_TagProbe_MC.root" : "Dimuon_TagProbe_Data.root","RECREATE");
    hMassUncorrected->Write();
    hMassCorrected->Write();
    hProbe_pt_total->Write();
    hProbe_pt_pass->Write();
    hProbe_eta_total->Write();
    hProbe_eta_pass->Write();
    hProbe_pteta_total->Write();
    hProbe_pteta_pass->Write();
    eff_pt->Write("Efficiency_pt");
    eff_eta->Write("Efficiency_eta");
    eff_pteta->Write("Efficiency_pt_eta");
    outFile.Close();

    std::cout << "Histograms and efficiencies saved to " << (isMC ? "Dimuon_TagProbe_MC.root" : "Dimuon_TagProbe_Data.root") << std::endl;
}
