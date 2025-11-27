// ============================================================================
//  FULL ANALYSIS MACRO — CLEAN, FORMATTED, WITH B-TAG VETO + ggH FIX
// ============================================================================

#include "RoccoR.h"
#include "JetCorrections.h"

#include <TRandom.h>
#include <TChain.h>
#include <TError.h>
#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TLorentzVector.h>
#include <TH1F.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

// ============================================================================
//  STRUCTURES & HELPERS
// ============================================================================
struct MuonSF {
    double ptMin, ptMax;
    double etaMin, etaMax;
    double sf;
};

struct PUWeight {
    int nMin, nMax;
    double weight;
};

double getMuonSF(double pt, double eta, const std::vector<MuonSF>& table) {
    for (const auto& b : table)
        if (pt >= b.ptMin && pt < b.ptMax && eta >= b.etaMin && eta < b.etaMax)
            return b.sf;
    return 1.0;
}

double getPUWeight(float nTruePU, const std::vector<PUWeight>& table) {
    for (const auto& b : table)
        if (nTruePU >= b.nMin && nTruePU < b.nMax)
            return b.weight;
    return 1.0;
}

std::vector<MuonSF> readMuonSF(const std::string& filename) {
    std::vector<MuonSF> table;
    std::ifstream f(filename);
    if (!f.is_open()) {
        std::cerr << "ERROR: cannot open muon SF file " << filename << std::endl;
        return table;
    }
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        MuonSF b;
        std::stringstream ss(line);
        ss >> b.ptMin >> b.ptMax >> b.etaMin >> b.etaMax >> b.sf;
        table.push_back(b);
    }
    return table;
}

std::vector<PUWeight> readPUWeights(const std::string& filename) {
    std::vector<PUWeight> table;
    std::ifstream f(filename);
    if (!f.is_open()) {
        std::cerr << "ERROR: cannot open PU weight file " << filename << std::endl;
        return table;
    }
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        PUWeight b;
        std::stringstream ss(line);
        ss >> b.nMin >> b.nMax >> b.weight;
        table.push_back(b);
    }
    return table;
}

// ============================================================================
//  MAIN ANALYSIS
// ============================================================================
void testmacro_FullCorrections_FullJetsMuons_weightMC_last_update(
    std::vector<std::string> inputFiles,
    bool isData = true,
    const std::string& puFileName =
        "/eos/user/c/cmsdqm/www/CAF/certification/Collisions22/PileUp/EFG/pileup_JSON.txt",
    const std::string& muPOGDir =
        "/afs/cern.ch/user/n/nbostan/new_CMS/CMSSW_14_0_18/src/RoccoR/corrections/2025-08-14/")
{
    gErrorIgnoreLevel = kError;

    const bool isMC = !isData;

    // Lumi + xsec must be set for MC
    double lumi_fb = 26.6717; // 2022postEE
    //double lumi_fb = 7.9804; // 2022preEE
    //double lumi_fb = 17.794; // 2023
    //double lumi_fb = 9.451; // 2023BPix
    //double lumi_fb = 108.960; // 2024
    
    
    
    double xsec_pb = 0.00105742; // Run 3 VBF cross section
    
    //double xsec_pb = 0.0135096; // Run 3 ggH cross section
    
    //double xsec_pb = 359.05; // Run 3 dyTo2L_M-50_2jets cross section
    
    //double xsec_pb = 98.04; // Run 3 TTto2L2Nu cross section
    
    //double xsec_pb = 1.421; // Run 3 EWK_2l_2jets cross section
    
    //double xsec_pb = 7.568;  // Run 3 WZ_2l2q cross section
    
    //double xsec_pb = 6.788; //Run 3 ZZ_2l2q cross section
    
     //double xsec_pb = 4.924; //Run 3 WZ_3lnu cross section
     
     //double xsec_pb = 1.39; //Run 3 ZZ_4l  cross section
     
     //double xsec_pb = 1.031; //Run 3 ZZ_2l2nu  cross section
     
     //double xsec_pb = 0.2328; //Run 3 WWW  cross section
     
     //double xsec_pb = 405.69; //Run 3 TTtoLNu2Q  cross section
     
     //double xsec_pb = 0.1851; //Run 3 WWZ  cross section
     
     //double xsec_pb = 0.062; //Run 3 WZZ  cross section
     
     //double xsec_pb = 0.0159; //Run 3 ZZZ  cross section

    // ----------------------------------------------------------------------------
    // Initialize RoccoR + jet corrections
    // ----------------------------------------------------------------------------
    RoccoR rc;
    rc.init("/afs/cern.ch/user/n/nbostan/new_CMS/CMSSW_14_0_18/src/RoccoR/corrections/RoccoR2022EE.txt");

    JetCorrections jetCorr(
        "/afs/cern.ch/user/n/nbostan/jet_jerc.txt",
        "/afs/cern.ch/user/n/nbostan/jetid.txt",
        "/afs/cern.ch/user/n/nbostan/jetvetomaps.txt"
    );


    // ----------------------------------------------------------------------------
    // Load externally-provided muon & PU SF tables
    // ----------------------------------------------------------------------------
    std::vector<MuonSF> muonSF_JPsi       = readMuonSF(muPOGDir + "muon_JPsi.txt");
    std::vector<MuonSF> muonSF_Z          = readMuonSF(muPOGDir + "muon_Z.txt");
    std::vector<MuonSF> muonSF_HighPt     = readMuonSF(muPOGDir + "muon_HighPt.txt");
    std::vector<MuonSF> muonSF_ScaleSmear = readMuonSF(muPOGDir + "muon_scalesmearing.txt");

    std::vector<PUWeight> puWeights;
    if (isMC) puWeights = readPUWeights(puFileName);

    // ========================================================================
    // Histograms
    // ========================================================================
    TH1F *h_mass          = new TH1F("h_mass",          "Dimuon mass;Mass [GeV];Events", 100, 0, 200);
    TH1F *h_dimuonPt      = new TH1F("h_dimuonPt",      "Dimuon pT; pT [GeV];Events", 100, 0, 200);
    TH1F *h_dimuonEta     = new TH1F("h_dimuonEta",     "Dimuon eta; eta; Events", 50, -2.5, 2.5);
    TH1F *h_leadJetPt     = new TH1F("h_leadJetPt",     "Leading jet pT; pT [GeV]; Events", 100, 0, 200);
    TH1F *h_leadJetEta    = new TH1F("h_leadJetEta",    "Leading jet eta; eta; Events", 50, -5, 5);
    TH1F *h_dijetPt       = new TH1F("h_dijetPt",       "Dijet pT; pT [GeV]; Events", 100, 0, 1000);
    TH1F *h_dijetMass     = new TH1F("h_dijetMass",     "Dijet mass; Mjj [GeV]; Events", 100, 0, 2000);
    TH1F *h_jetPtCorr     = new TH1F("h_jetPtCorr",     "Corrected Jet pT; pT [GeV]; Jets", 100, 0, 1000);

    TH1F *h_mass_VBF      = new TH1F("h_mass_VBF",      "Dimuon mass (VBF);Mass [GeV];Events", 100, 0, 200);
    TH1F *h_mass_ggH      = new TH1F("h_mass_ggH",      "Dimuon mass (ggH);Mass [GeV];Events", 100, 0, 200);

    // ========================================================================
    // Build TChain
    // ========================================================================
    TChain chain("Events");
    for (auto& f : inputFiles) chain.Add(f.c_str());

    // ========================================================================
    // Compute MC normalization
    // ========================================================================
    double normFactor = 1.0;

    if (isMC) {
        double lumi_pb = lumi_fb * 1000.0;
        double nGenWeighted = 0;

        TTreeReader rcnt(&chain);
        TTreeReaderValue<float> genW(rcnt, "genWeight");

        while (rcnt.Next()) {
            nGenWeighted += (*genW >= 0 ? 1.0 : -1.0);
        }

        normFactor = (xsec_pb * lumi_pb) / nGenWeighted;

        std::cout << "Generated (signed): " << nGenWeighted << "\n"
                  << "Normalization factor: " << normFactor << std::endl;
    }

    // ========================================================================
    // Reader Setup
    // ========================================================================
    TTreeReader reader(&chain);

    TTreeReaderValue<float> genWeight(reader, "genWeight");
    TTreeReaderValue<bool>  HLT_IsoMu24(reader, "HLT_IsoMu24");

    TTreeReaderValue<int>   nMuon(reader, "nMuon");
    TTreeReaderArray<float> Muon_pt(reader, "Muon_pt");
    TTreeReaderArray<float> Muon_eta(reader, "Muon_eta");
    TTreeReaderArray<float> Muon_phi(reader, "Muon_phi");
    TTreeReaderArray<float> Muon_iso(reader, "Muon_pfRelIso04_all");
    TTreeReaderArray<int>   Muon_charge(reader, "Muon_charge");
    TTreeReaderArray<bool>  Muon_mediumID(reader, "Muon_mediumId");

    TTreeReaderValue<int>   nJet(reader, "nJet");
    TTreeReaderArray<float> Jet_pt(reader, "Jet_pt");
    TTreeReaderArray<float> Jet_eta(reader, "Jet_eta");
    TTreeReaderArray<float> Jet_phi(reader, "Jet_phi");
    TTreeReaderArray<float> Jet_mass(reader, "Jet_mass");
    TTreeReaderArray<float> Jet_raw(reader, "Jet_rawFactor");

    TTreeReaderArray<float>* Jet_btag = nullptr;
    if (chain.GetBranch("Jet_btagDeepFlavB"))
        Jet_btag = new TTreeReaderArray<float>(reader, "Jet_btagDeepFlavB");
    else if (chain.GetBranch("Jet_btagDeepB"))
        Jet_btag = new TTreeReaderArray<float>(reader, "Jet_btagDeepB");

    TTreeReaderValue<float>* nTruePU = nullptr;
    if (isMC && chain.GetBranch("Pileup_nTrueInt"))
        nTruePU = new TTreeReaderValue<float>(reader, "Pileup_nTrueInt");

    // ========================================================================
    // EVENT LOOP
    // ========================================================================
    int nPassed = 0;

    while (reader.Next()) {

        if (!(*HLT_IsoMu24)) continue;
        if (*nMuon != 2) continue;

        double w = 1.0;
        if (isMC) {
            double pu = getPUWeight(**nTruePU, puWeights);
            double gen = (*genWeight >= 0 ? 1.0 : -1.0);
            w = pu * gen * normFactor;
        }

        // --------------------------------------------------
        // TAG MUON
        // --------------------------------------------------
        int tagIdx = -1;
        float bestPt = -1;

        for (int i = 0; i < *nMuon; i++) {
            if (!Muon_mediumID[i]) continue;
            if (Muon_pt[i] < 26) continue;
            if (std::fabs(Muon_eta[i]) > 2.4) continue;
            if (Muon_iso[i] > 0.15) continue;

            if (Muon_pt[i] > bestPt) {
                bestPt = Muon_pt[i];
                tagIdx = i;
            }
        }
        if (tagIdx < 0) continue;

        // --------------------------------------------------
        // BUILD DIMUON
        // --------------------------------------------------
        TLorentzVector dimuon;
        bool foundPair = false;

        for (int j = 0; j < *nMuon; j++) {
            if (j == tagIdx) continue;
            if (Muon_charge[tagIdx] * Muon_charge[j] >= 0) continue;

            TLorentzVector t, p;
            t.SetPtEtaPhiM(Muon_pt[tagIdx], Muon_eta[tagIdx], Muon_phi[tagIdx], 0.105);
            p.SetPtEtaPhiM(Muon_pt[j], Muon_eta[j], Muon_phi[j], 0.105);
            dimuon = t + p;
            foundPair = true;
        }
        if (!foundPair) continue;

        nPassed++;

        h_mass->Fill(dimuon.M(), w);
        h_dimuonPt->Fill(dimuon.Pt(), w);
        h_dimuonEta->Fill(dimuon.Eta(), w);

        // ====================================================================
        //  JET COLLECTION (+ b-tag veto logic)
        // ====================================================================
        std::vector<TLorentzVector> jets;
        int nLooseB = 0;
        int nMedB   = 0;

        for (int j = 0; j < *nJet; j++) {
            double corrPt =
                jetCorr.getCorrectedJetPt(Jet_pt[j], Jet_eta[j], Jet_phi[j], Jet_raw[j], isMC);

            if (corrPt < 25) continue;
            if (std::fabs(Jet_eta[j]) > 4.7) continue;

            TLorentzVector jet;
            jet.SetPtEtaPhiM(corrPt, Jet_eta[j], Jet_phi[j], Jet_mass[j]);
            jets.push_back(jet);

            h_jetPtCorr->Fill(corrPt, w);

            if (Jet_btag) {
                float b = (*Jet_btag)[j];
                if      (b > 0.2783) nMedB++;
                else if (b > 0.0490) nLooseB++;
            }
        }

        std::sort(jets.begin(), jets.end(),
                  [](auto& a, auto& b){ return a.Pt() > b.Pt(); });

        if (!jets.empty()) {
            h_leadJetPt->Fill(jets[0].Pt(), w);
            h_leadJetEta->Fill(jets[0].Eta(), w);
        }

        // ====================================================================
        //  VBF CLASSIFICATION (WITH ggH FIX)
        // ====================================================================
        bool isVBF = false;
        bool isGGH = false;

        if (jets.size() >= 2) {
            TLorentzVector j1 = jets[0];
            TLorentzVector j2 = jets[1];

            double mjj = (j1 + j2).M();
            double deta = std::fabs(j1.Eta() - j2.Eta());
            double dipt = (j1 + j2).Pt();

            bool passVBFkin =
                (mjj > 400) &&
                (deta > 2.5) &&
                (j1.Pt() > 35) &&
                (j2.Pt() > 25);

            bool passBtagVeto =
                (nLooseB < 2 && nMedB < 1);

            bool exceptionRule = false;
            if (!passBtagVeto && passVBFkin) {
                // keep VBF if kinematics are strong
                exceptionRule = true;
            }

            // First decide VBF
            if (passVBFkin && (passBtagVeto || exceptionRule)) {
                isVBF = true;
            } else {
                // Candidate for ggH — now apply ggH-specific b-tag veto:
                // - reject if >=2 loose b-tags
                // - reject if >=1 medium AND Njets >= 2
                bool passGgHBtag = true;

                if (nLooseB >= 2) passGgHBtag = false;

                if (jets.size() >= 2 && nMedB >= 1) passGgHBtag = false;

                // Exception: if Njets == 1 and exactly 1 medium b-tag -> do NOT reject
                // (this branch is Njets>=2 so nothing extra to do here)

                if (passGgHBtag) isGGH = true;
            }

            // Fill dijet histos for 2+ jet events
            h_dijetPt->Fill(dipt, w);
            h_dijetMass->Fill(mjj, w);
        }
        else if (jets.size() == 1) {
            // 1-jet events: treat as ggH candidate unless vetoed by >=2 loose b-tags (impossible here)
            // Special rule: Njets == 1 and nMedB == 1 are NOT rejected → accept
            bool passGgHBtag = true;
            if (nLooseB >= 2) passGgHBtag = false; // unlikely with 1 jet
            // for 1-jet, a single medium b-tag is allowed
            if (passGgHBtag) isGGH = true;
        }

        if (isVBF) h_mass_VBF->Fill(dimuon.M(), w);
        if (isGGH) h_mass_ggH->Fill(dimuon.M(), w);
    }

    std::cout << "Events passing all selections: " << nPassed << std::endl;

    // ========================================================================
    //  SAVE OUTPUT
    // ========================================================================
    TFile out("output_histos.root", "RECREATE");

    h_mass->Write();
    h_dimuonPt->Write();
    h_dimuonEta->Write();
    h_leadJetPt->Write();
    h_leadJetEta->Write();
    h_dijetPt->Write();
    h_dijetMass->Write();
    h_jetPtCorr->Write();

    h_mass_VBF->Write();
    h_mass_ggH->Write();

    out.Close();
    std::cout << "Histograms saved to output_histos.root" << std::endl;
}
