#include "RoccoR.h"

#include <TRandom.h>
#include <TChain.h>
#include <TError.h>
#include "JetCorrections.h"
#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TLorentzVector.h>
#include <TH1F.h>
#include <TH2F.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

// -----------------------------
// Structs and helper functions
// -----------------------------
struct MuonSF {
    double ptMin, ptMax;
    double etaMin, etaMax;
    double sf;
};

struct PUWeight {
    int nMin, nMax;
    double weight;
};

double getMuonSF(double pt, double eta, const std::vector<MuonSF>& sfTable) {
    for (const auto &bin : sfTable)
        if (pt >= bin.ptMin && pt < bin.ptMax && eta >= bin.etaMin && eta < bin.etaMax)
            return bin.sf;
    return 1.0;
}

double getPUWeight(float nTruePU, const std::vector<PUWeight>& puTable) {
    for (const auto &bin : puTable)
        if (nTruePU >= bin.nMin && nTruePU < bin.nMax)
            return bin.weight;
    return 1.0;
}

std::vector<MuonSF> readMuonSF(const std::string& filename) {
    std::vector<MuonSF> table;
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Cannot open " << filename << std::endl;
        return table;
    }
    std::string line;
    while (std::getline(infile, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::stringstream ss(line);
        MuonSF bin;
        ss >> bin.ptMin >> bin.ptMax >> bin.etaMin >> bin.etaMax >> bin.sf;
        table.push_back(bin);
    }
    infile.close();
    return table;
}

std::vector<PUWeight> readPUWeights(const std::string& filename) {
    std::vector<PUWeight> table;
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Cannot open " << filename << std::endl;
        return table;
    }
    std::string line;
    while (std::getline(infile, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::stringstream ss(line);
        PUWeight bin;
        ss >> bin.nMin >> bin.nMax >> bin.weight;
        table.push_back(bin);
    }
    infile.close();
    return table;
}

// =====================================================
// Main analysis macro
// =====================================================
void testmacro_FullCorrections_FullJetsMuons_Data_last_update(std::vector<std::string> inputFiles,
                                             bool isData = true,
                                             const std::string &puFileName = "/eos/user/c/cmsdqm/www/CAF/certification/Collisions22/PileUp/EFG/pileup_JSON.txt",
                                             const std::string &muPOGDir   = "/afs/cern.ch/user/n/nbostan/new_CMS/CMSSW_14_0_18/src/RoccoR/corrections/2025-08-14/")
{
    bool isMC = !isData;

    gErrorIgnoreLevel = kError; // suppress harmless ROOT warnings

    // Initialize RoccoR (muon momentum corrections)
    RoccoR rc;
    rc.init("/afs/cern.ch/user/n/nbostan/new_CMS/CMSSW_14_0_18/src/RoccoR/corrections/RoccoR2022EE.txt");

    // Read muon scale/SF files (optional)
    std::vector<MuonSF> muonJPsi       = readMuonSF(muPOGDir + "muon_JPsi.txt");
    std::vector<MuonSF> muonZ          = readMuonSF(muPOGDir + "muon_Z.txt");
    std::vector<MuonSF> muonHighPt     = readMuonSF(muPOGDir + "muon_HighPt.txt");
    std::vector<MuonSF> muonScaleSmear = readMuonSF(muPOGDir + "muon_scalesmearing.txt");

    std::vector<PUWeight> puWeights;
    if (isMC) puWeights = readPUWeights(puFileName);

    // Jet corrections
    JetCorrections jetCorr(
        "/afs/cern.ch/user/n/nbostan/jet_jerc.txt",
        "/afs/cern.ch/user/n/nbostan/jetid.txt",
        "/afs/cern.ch/user/n/nbostan/jetvetomaps.txt"
    );

    // --------------------------------------------------
    // Histograms (inclusive + categories)
    // --------------------------------------------------
    // Inclusive
    TH1F *h_mass       = new TH1F("h_mass",       "Dimuon invariant mass;M_{#mu#mu} [GeV];Events", 100, 0, 200);
    TH1F *h_dimuonPt   = new TH1F("h_dimuonPt",   "Dimuon p_{T};p_{T}^{#mu#mu} [GeV];Events", 100, 0, 200);
    TH1F *h_dimuonEta  = new TH1F("h_dimuonEta",  "Dimuon #eta;#eta^{#mu#mu};Events", 50, -2.5, 2.5);
    TH1F *h_leadJetPt  = new TH1F("h_leadJetPt",  "Leading jet p_{T};p_{T}^{jet1} [GeV];Events", 100, 0, 200);
    TH1F *h_leadJetEta = new TH1F("h_leadJetEta", "Leading jet #eta;#eta^{jet1};Events", 50, -5, 5);
    TH1F *h_dijetPt    = new TH1F("h_dijetPt",    "Dijet p_{T};p_{T}^{jj} [GeV];Events", 100, 0, 1000);
    TH1F *h_dijetMass  = new TH1F("h_dijetMass",  "Dijet invariant mass;M_{jj} [GeV];Events", 100, 0, 2000);
    TH1F *h_jetPtCorr  = new TH1F("h_jetPtCorr",  "Corrected jet p_{T};p_{T}^{corr} [GeV];Jets", 100, 0, 1000);

    // VBF
    TH1F *h_mass_VBF       = new TH1F("h_mass_VBF",       "Dimuon invariant mass (VBF);M_{#mu#mu} [GeV];Events", 100, 0, 200);
    TH1F *h_dimuonPt_VBF   = new TH1F("h_dimuonPt_VBF",   "Dimuon p_{T} (VBF);p_{T}^{#mu#mu} [GeV];Events", 100, 0, 200);
    TH1F *h_dimuonEta_VBF  = new TH1F("h_dimuonEta_VBF",  "Dimuon #eta (VBF);#eta^{#mu#mu};Events", 50, -2.5, 2.5);
    TH1F *h_leadJetPt_VBF  = new TH1F("h_leadJetPt_VBF",  "Leading jet p_{T} (VBF);p_{T}^{jet1} [GeV];Events", 100, 0, 200);
    TH1F *h_leadJetEta_VBF = new TH1F("h_leadJetEta_VBF", "Leading jet #eta (VBF);#eta^{jet1};Events", 50, -5, 5);
    TH1F *h_dijetPt_VBF    = new TH1F("h_dijetPt_VBF",    "Dijet p_{T} (VBF);p_{T}^{jj} [GeV];Events", 100, 0, 1000);
    TH1F *h_dijetMass_VBF  = new TH1F("h_dijetMass_VBF",  "Dijet invariant mass (VBF);M_{jj} [GeV];Events", 100, 0, 2000);

    // ggH
    TH1F *h_mass_ggH       = new TH1F("h_mass_ggH",       "Dimuon invariant mass (ggH);M_{#mu#mu} [GeV];Events", 100, 0, 200);
    TH1F *h_dimuonPt_ggH   = new TH1F("h_dimuonPt_ggH",   "Dimuon p_{T} (ggH);p_{T}^{#mu#mu} [GeV];Events", 100, 0, 200);
    TH1F *h_dimuonEta_ggH  = new TH1F("h_dimuonEta_ggH",  "Dimuon #eta (ggH);#eta^{#mu#mu};Events", 50, -2.5, 2.5);
    TH1F *h_leadJetPt_ggH  = new TH1F("h_leadJetPt_ggH",  "Leading jet p_{T} (ggH);p_{T}^{jet1} [GeV];Events", 100, 0, 200);
    TH1F *h_leadJetEta_ggH = new TH1F("h_leadJetEta_ggH", "Leading jet #eta (ggH);#eta^{jet1};Events", 50, -5, 5);
    TH1F *h_dijetPt_ggH    = new TH1F("h_dijetPt_ggH",    "Dijet p_{T} (ggH);p_{T}^{jj} [GeV];Events", 100, 0, 1000);
    TH1F *h_dijetMass_ggH  = new TH1F("h_dijetMass_ggH",  "Dijet invariant mass (ggH);M_{jj} [GeV];Events", 100, 0, 2000);

    // --------------------------------------------------
    // Tree setup
    // --------------------------------------------------
    TChain chain("Events");
    for (auto &f : inputFiles) chain.Add(f.c_str());

    TTreeReader reader(&chain);

    TTreeReaderValue<Int_t> nMuon(reader, "nMuon");
    TTreeReaderArray<Float_t> Muon_pt(reader, "Muon_pt");
    TTreeReaderArray<Float_t> Muon_eta(reader, "Muon_eta");
    TTreeReaderArray<Float_t> Muon_phi(reader, "Muon_phi");
    TTreeReaderArray<Int_t> Muon_charge(reader, "Muon_charge");
    TTreeReaderArray<Bool_t> Muon_mediumId(reader, "Muon_mediumId");
    TTreeReaderArray<Float_t> Muon_pfRelIso04_all(reader, "Muon_pfRelIso04_all");
    TTreeReaderValue<Bool_t> HLT_IsoMu24(reader, "HLT_IsoMu24");

    TTreeReaderValue<Int_t> nJet(reader, "nJet");
    TTreeReaderArray<Float_t> Jet_pt(reader, "Jet_pt");
    TTreeReaderArray<Float_t> Jet_eta(reader, "Jet_eta");
    TTreeReaderArray<Float_t> Jet_phi(reader, "Jet_phi");
    TTreeReaderArray<Float_t> Jet_mass(reader, "Jet_mass");
    TTreeReaderArray<Float_t> Jet_rawFactor(reader, "Jet_rawFactor");

    // Handle multiple possible b-tag branches
    TTreeReaderArray<Float_t> *Jet_btag = nullptr;
    if (chain.GetBranch("Jet_btagDeepB"))
        Jet_btag = new TTreeReaderArray<Float_t>(reader, "Jet_btagDeepB"); // DeepCSV-like (DeepB)
    else if (chain.GetBranch("Jet_btagDeepFlavB"))
        Jet_btag = new TTreeReaderArray<Float_t>(reader, "Jet_btagDeepFlavB"); // DeepJet fallback
    else
        std::cerr << "⚠️ Warning: No known b-tag branch found (DeepB or DeepFlavB)!" << std::endl;

    TTreeReaderValue<Float_t> *Pileup_nTrueInt = nullptr;
    if (isMC && chain.GetBranch("Pileup_nTrueInt"))
        Pileup_nTrueInt = new TTreeReaderValue<Float_t>(reader, "Pileup_nTrueInt");

    int nEventPassed = 0;

    // b-tag thresholds (DeepCSV loose & medium working points assumed here)
    const float btagLooseWP  = 0.0490; // loose
    const float btagMediumWP = 0.2783; // medium

    // =====================================================
    // Event loop
    // =====================================================
    while (reader.Next()) {
        if (!(*HLT_IsoMu24) || *nMuon != 2) continue; 
        nEventPassed++;

        double puWeight = 1.0;
        if (isMC && Pileup_nTrueInt)
            puWeight = getPUWeight(**Pileup_nTrueInt, puWeights);

        // --- Muon selection: tag ---
        int tagIdx = -1; float maxPt = -1;
        for (int i = 0; i < *nMuon; i++) {
            if (!Muon_mediumId[i]) continue;
            if (Muon_pt[i] < 26 || fabs(Muon_eta[i]) > 2.4 || Muon_pfRelIso04_all[i] > 0.15) continue;
            if (Muon_pt[i] > maxPt) { maxPt = Muon_pt[i]; tagIdx = i; }
        }
        if (tagIdx < 0) continue;

        // --- Find dimuon pairs (tag-probe style) ---
        TLorentzVector dimuon;
        bool foundPair = false;
        for (int p = 0; p < *nMuon; p++) {
            if (p == tagIdx) continue;
            if (Muon_charge[tagIdx] * Muon_charge[p] >= 0) continue;

            TLorentzVector tag, probe;
            tag.SetPtEtaPhiM(Muon_pt[tagIdx], Muon_eta[tagIdx], Muon_phi[tagIdx], 0.105);
            probe.SetPtEtaPhiM(Muon_pt[p], Muon_eta[p], Muon_phi[p], 0.105);
            dimuon = tag + probe;
            foundPair = true;

            // Inclusive fills
            h_mass->Fill(dimuon.M(), puWeight);
            h_dimuonPt->Fill(dimuon.Pt(), puWeight);
            h_dimuonEta->Fill(dimuon.Eta(), puWeight);
        }
        if (!foundPair) continue;

        // --- Jet corrections and b-tagging ---
        std::vector<TLorentzVector> jets;
        int nLooseBtag = 0, nMediumBtag = 0;

        for (int j = 0; j < *nJet; j++) {
            double corrPt = jetCorr.getCorrectedJetPt(Jet_pt[j], Jet_eta[j], Jet_phi[j], Jet_rawFactor[j], isMC);
            if (corrPt < 25 || fabs(Jet_eta[j]) > 4.7) continue; // keep jets down to 25 GeV (we'll apply tighter requirements for categories)

            TLorentzVector jet;
            jet.SetPtEtaPhiM(corrPt, Jet_eta[j], Jet_phi[j], Jet_mass[j]);
            jets.push_back(jet);
            h_jetPtCorr->Fill(corrPt, puWeight);

            if (Jet_btag) {
                float bval = (*Jet_btag)[j];
                if (bval > btagLooseWP)  nLooseBtag++;
                if (bval > btagMediumWP) nMediumBtag++;
            }
        }

        // Sort jets by pT descending
        if (!jets.empty()) {
            std::sort(jets.begin(), jets.end(), [](const TLorentzVector &a, const TLorentzVector &b){ return a.Pt() > b.Pt(); });
            // Inclusive leading jet fill
            h_leadJetPt->Fill(jets[0].Pt(), puWeight);
            h_leadJetEta->Fill(jets[0].Eta(), puWeight);
        }

        // --- Categorization ---
        bool isVBF = false;
        bool isGGH = false;

        // VBF selection according to provided criteria:
        // • At least two jets with pT > 35 GeV and > 25 GeV (leading >35, subleading >25)
        // • mjj > 400 GeV and |Δη_jj| > 2.5
        // • Veto events with: >=2 loose b-tagged jets OR >=1 medium b-tagged jet
        if (jets.size() >= 2) {
            const TLorentzVector &j1 = jets[0];
            const TLorentzVector &j2 = jets[1];

            bool jetPtSel = (j1.Pt() > 35.0 && j2.Pt() > 25.0);
            double mjj  = (j1 + j2).M();
            double deta = fabs(j1.Eta() - j2.Eta());
            bool dijetSel = (mjj > 400.0 && deta > 2.5);

            bool bVeto = (nLooseBtag >= 2) || (nMediumBtag >= 1);

            if (jetPtSel && dijetSel && !bVeto) {
                isVBF = true;
            }
        }

        // ggH selection according to provided criteria:
        // • Accept events failing VBF category selection
        // • Veto events with: >=2 loose b-tagged jets OR >=1 medium b-tagged jet with at least two jets
        // • Events with Njets == 1 and 1 medium b-tagged jet are not rejected
        if (!isVBF) {
            bool bVeto = false;
            if (nLooseBtag >= 2) bVeto = true;
            if (nMediumBtag >= 1 && jets.size() >= 2) bVeto = true;
            // keep events where Njets == 1 and exactly 1 medium b-tag
            if (jets.size() == 1 && nMediumBtag == 1) bVeto = false;

            if (!bVeto) isGGH = true;
        }

        // --- Fill category histograms ---
        if (isVBF && jets.size() >= 2) {
            TLorentzVector dijet = jets[0] + jets[1];
            h_mass_VBF->Fill(dimuon.M(), puWeight);
            h_dimuonPt_VBF->Fill(dimuon.Pt(), puWeight);
            h_dimuonEta_VBF->Fill(dimuon.Eta(), puWeight);
            h_dijetPt_VBF->Fill(dijet.Pt(), puWeight);
            h_dijetMass_VBF->Fill(dijet.M(), puWeight);

            // Leading jet in VBF
            h_leadJetPt_VBF->Fill(jets[0].Pt(), puWeight);
            h_leadJetEta_VBF->Fill(jets[0].Eta(), puWeight);
        }

        if (isGGH) {
            h_mass_ggH->Fill(dimuon.M(), puWeight);
            h_dimuonPt_ggH->Fill(dimuon.Pt(), puWeight);
            h_dimuonEta_ggH->Fill(dimuon.Eta(), puWeight);

            if (jets.size() >= 1) {
                h_leadJetPt_ggH->Fill(jets[0].Pt(), puWeight);
                h_leadJetEta_ggH->Fill(jets[0].Eta(), puWeight);
            }
            if (jets.size() >= 2) {
                TLorentzVector dijet = jets[0] + jets[1];
                h_dijetPt_ggH->Fill(dijet.Pt(), puWeight);
                h_dijetMass_ggH->Fill(dijet.M(), puWeight);
            }
        }
    } // end event loop

    // =====================================================
    // Write output ROOT file
    // =====================================================
    std::string outFileName = isMC ? "output_FullCorrections_MC.root" : "output_FullCorrections_Data.root";
    TFile outFile(outFileName.c_str(), "RECREATE");

    // Inclusive
    h_mass->Write();
    h_dimuonPt->Write();
    h_dimuonEta->Write();
    h_leadJetPt->Write();
    h_leadJetEta->Write();
    h_dijetPt->Write();
    h_dijetMass->Write();
    h_jetPtCorr->Write();

    // VBF
    h_mass_VBF->Write();
    h_dimuonPt_VBF->Write();
    h_dimuonEta_VBF->Write();
    h_leadJetPt_VBF->Write();
    h_leadJetEta_VBF->Write();
    h_dijetPt_VBF->Write();
    h_dijetMass_VBF->Write();

    // ggH
    h_mass_ggH->Write();
    h_dimuonPt_ggH->Write();
    h_dimuonEta_ggH->Write();
    h_leadJetPt_ggH->Write();
    h_leadJetEta_ggH->Write();
    h_dijetPt_ggH->Write();
    h_dijetMass_ggH->Write();

    outFile.Close();

    std::cout << "\n✅ Analysis complete.\n";
    std::cout << "Output file written: " << outFileName << std::endl;
    std::cout << "Events passing trigger & tag selection: " << nEventPassed << std::endl;
    std::cout << "Histograms included:\n"
              << "  • Inclusive: h_mass, h_dimuonPt, h_dimuonEta, h_leadJetPt, h_leadJetEta, h_dijetPt, h_dijetMass, h_jetPtCorr\n"
              << "  • VBF: h_mass_VBF, h_dimuonPt_VBF, h_dimuonEta_VBF, h_leadJetPt_VBF, h_leadJetEta_VBF, h_dijetPt_VBF, h_dijetMass_VBF\n"
              << "  • ggH: h_mass_ggH, h_dimuonPt_ggH, h_dimuonEta_ggH, h_leadJetPt_ggH, h_leadJetEta_ggH, h_dijetPt_ggH, h_dijetMass_ggH\n"
              << std::endl;
}
