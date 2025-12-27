// ============================================================================
// FULL ANALYSIS MACRO — ROOT ONLY, KIT TEXT MUON CORRECTIONS
// ============================================================================

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
// STRUCTS
// ============================================================================
struct MuonSF { double ptMin, ptMax, etaMin, etaMax, sf; };
struct MuonScale { double ptMin, ptMax, etaMin, etaMax, scale; };
struct PUWeight { int nMin, nMax; double weight; };

// ============================================================================
// HELPERS
// ============================================================================
double getMuonSF(double pt, double eta, const std::vector<MuonSF>& table) {
    for (const auto& b : table)
        if (pt >= b.ptMin && pt < b.ptMax &&
            eta >= b.etaMin && eta < b.etaMax)
            return b.sf;
    return 1.0;
}

double getPUWeight(float nTruePU, const std::vector<PUWeight>& table) {
    for (const auto& b : table)
        if (nTruePU >= b.nMin && nTruePU < b.nMax)
            return b.weight;
    return 1.0;
}

double getCorrectedMuonPt(
    double pt, double eta, bool isMC,
    const std::vector<MuonScale>& scaleTable
) {
    if (!isMC) return pt;
    for (const auto& b : scaleTable)
        if (pt >= b.ptMin && pt < b.ptMax &&
            eta >= b.etaMin && eta < b.etaMax)
            return pt * b.scale;
    return pt;
}

// ============================================================================
// FILE READERS
// ============================================================================
std::vector<MuonSF> readMuonSF(const std::string& filename) {
    std::vector<MuonSF> table;
    std::ifstream f(filename);
    MuonSF b;
    while (f >> b.ptMin >> b.ptMax >> b.etaMin >> b.etaMax >> b.sf)
        table.push_back(b);
    return table;
}

std::vector<MuonScale> readMuonScale(const std::string& filename) {
    std::vector<MuonScale> table;
    std::ifstream f(filename);
    MuonScale b;
    while (f >> b.ptMin >> b.ptMax >> b.etaMin >> b.etaMax >> b.scale)
        table.push_back(b);
    return table;
}

std::vector<PUWeight> readPUWeights(const std::string& filename) {
    std::vector<PUWeight> table;
    std::ifstream f(filename);
    PUWeight b;
    while (f >> b.nMin >> b.nMax >> b.weight)
        table.push_back(b);
    return table;
}

// ============================================================================
// MAIN
// ============================================================================
void testmacro_FullCorrections_FullJetsMuons_weightMC_lastv_2024_kit(
    std::vector<std::string> inputFiles,
    bool isData = true
) {
    gErrorIgnoreLevel = kError;
    const bool isMC = !isData;

    double lumi_fb = 11.47;
    double xsec_pb = 0.06206;

    // =========================================================================
    // INPUT TABLES
    // =========================================================================
    auto muonScale = readMuonScale("/*/src/RoccoR/post2022E-update/2024_Summer24.txt");
    auto muonSF_Z  = readMuonSF("/*/src/RoccoR/post2022E-update/2024kit/muon_Z.txt");
    auto puWeights = readPUWeights("/eos/user/c/cmsdqm/www/CAF/certification/Collisions24/PileUp/pileup_JSON-2024BCDEFGHI.txt");

    JetCorrections jetCorr(
        "/*/src/RoccoR/post2022E-update/2024kit/jet_jerc.txt",
        "/*/src/RoccoR/post2022E-update/2024kit/jetid.txt",
        "/*/src/RoccoR/post2022E-update/2024kit/jetvetomaps.txt"
    );

    // =========================================================================
    // HISTOGRAMS — EXACTLY AS REQUESTED
    // =========================================================================
    TH1F *h_mass        = new TH1F("h_mass","Dimuon mass",100,0,200);
    TH1F *h_dimuonPt   = new TH1F("h_dimuonPt","Dimuon pT",100,0,200);
    TH1F *h_dimuonEta  = new TH1F("h_dimuonEta","Dimuon eta",50,-2.5,2.5);
    TH1F *h_leadJetPt  = new TH1F("h_leadJetPt","Leading jet pT",100,0,200);
    TH1F *h_leadJetEta = new TH1F("h_leadJetEta","Leading jet eta",50,-5,5);
    TH1F *h_dijetPt    = new TH1F("h_dijetPt","Dijet pT",100,0,1000);
    TH1F *h_dijetMass  = new TH1F("h_dijetMass","Dijet mass",100,0,2000);
    TH1F *h_jetPtCorr  = new TH1F("h_jetPtCorr","Corrected jet pT",100,0,1000);
    TH1F *h_mass_VBF   = new TH1F("h_mass_VBF","VBF dimuon mass",100,0,200);
    TH1F *h_mass_ggH   = new TH1F("h_mass_ggH","ggH dimuon mass",100,0,200);

    for (auto h : {h_mass,h_dimuonPt,h_dimuonEta,h_leadJetPt,h_leadJetEta,
                   h_dijetPt,h_dijetMass,h_jetPtCorr,h_mass_VBF,h_mass_ggH})
        h->Sumw2();

    // =========================================================================
    // CHAIN
    // =========================================================================
    TChain chain("Events");
    for (auto& f : inputFiles) chain.Add(f.c_str());

    // =========================================================================
    // MC NORMALIZATION
    // =========================================================================
    double normFactor = 1.0;
    if (isMC) {
        double lumi_pb = lumi_fb * 1000.0;
        double sumW = 0;
        TTreeReader r(&chain);
        TTreeReaderValue<float> genW(r,"genWeight");
        while (r.Next()) sumW += *genW;
        normFactor = (xsec_pb * lumi_pb) / sumW;
    }

    // =========================================================================
    // TREE READER
    // =========================================================================
    TTreeReader reader(&chain);

    TTreeReaderValue<float> genWeight(reader,"genWeight");
    TTreeReaderValue<bool>  HLT_IsoMu24(reader,"HLT_IsoMu24");

    TTreeReaderValue<int> nMuon(reader,"nMuon");
    TTreeReaderArray<float> Muon_pt(reader,"Muon_pt");
    TTreeReaderArray<float> Muon_eta(reader,"Muon_eta");
    TTreeReaderArray<float> Muon_phi(reader,"Muon_phi");
    TTreeReaderArray<float> Muon_iso(reader,"Muon_pfRelIso04_all");
    TTreeReaderArray<int>   Muon_charge(reader,"Muon_charge");
    TTreeReaderArray<bool>  Muon_mediumID(reader,"Muon_mediumId");

    TTreeReaderValue<int> nJet(reader,"nJet");
    TTreeReaderArray<float> Jet_pt(reader,"Jet_pt");
    TTreeReaderArray<float> Jet_eta(reader,"Jet_eta");
    TTreeReaderArray<float> Jet_phi(reader,"Jet_phi");
    TTreeReaderArray<float> Jet_mass(reader,"Jet_mass");
    TTreeReaderArray<float> Jet_raw(reader,"Jet_rawFactor");

    TTreeReaderArray<float>* Jet_btag = nullptr;
    if (chain.GetBranch("Jet_btagDeepFlavB"))
        Jet_btag = new TTreeReaderArray<float>(reader,"Jet_btagDeepFlavB");

    TTreeReaderValue<float>* nTruePU = nullptr;
    if (isMC && chain.GetBranch("Pileup_nTrueInt"))
        nTruePU = new TTreeReaderValue<float>(reader,"Pileup_nTrueInt");

    // =========================================================================
    // EVENT LOOP — UNCHANGED
    // =========================================================================
    while (reader.Next()) {

        if (!(*HLT_IsoMu24)) continue;
        if (*nMuon != 2) continue;

        double w = 1.0;
        if (isMC) {
            w = (*genWeight) * normFactor;
            if (nTruePU) w *= getPUWeight(**nTruePU, puWeights);
        }

        int tag = -1;
        double bestPt = -1;

        for (int i = 0; i < 2; i++) {
            if (!Muon_mediumID[i]) continue;
            if (Muon_iso[i] > 0.25) continue;
            if (fabs(Muon_eta[i]) > 2.4) continue;

            double pt = getCorrectedMuonPt(Muon_pt[i], Muon_eta[i], isMC, muonScale);
            if (pt < 26) continue;

            if (pt > bestPt) { bestPt = pt; tag = i; }
        }

        if (tag < 0) continue;
        int probe = 1 - tag;
        if (Muon_charge[tag] * Muon_charge[probe] >= 0) continue;

        double pt1 = getCorrectedMuonPt(Muon_pt[tag], Muon_eta[tag], isMC, muonScale);
        double pt2 = getCorrectedMuonPt(Muon_pt[probe], Muon_eta[probe], isMC, muonScale);

        if (isMC) {
            w *= getMuonSF(pt1, Muon_eta[tag], muonSF_Z);
            w *= getMuonSF(pt2, Muon_eta[probe], muonSF_Z);
        }

        TLorentzVector m1, m2;
        m1.SetPtEtaPhiM(pt1, Muon_eta[tag], Muon_phi[tag], 0.105);
        m2.SetPtEtaPhiM(pt2, Muon_eta[probe], Muon_phi[probe], 0.105);
        TLorentzVector dimuon = m1 + m2;

        h_mass->Fill(dimuon.M(), w);
        h_dimuonPt->Fill(dimuon.Pt(), w);
        h_dimuonEta->Fill(dimuon.Eta(), w);

        // ======================= JETS (UNCHANGED) =======================
        std::vector<TLorentzVector> jets;
        int nLooseB = 0, nMedB = 0;

        for (int j = 0; j < *nJet; j++) {
            double corrPt = jetCorr.getCorrectedJetPt(
                Jet_pt[j], Jet_eta[j], Jet_phi[j], Jet_raw[j], isMC
            );
            if (corrPt < 25 || fabs(Jet_eta[j]) > 4.7) continue;

            TLorentzVector jet;
            jet.SetPtEtaPhiM(corrPt, Jet_eta[j], Jet_phi[j], Jet_mass[j]);
            jets.push_back(jet);
            h_jetPtCorr->Fill(corrPt, w);

            if (Jet_btag) {
                float b = (*Jet_btag)[j];
                if (b > 0.2783) nMedB++;
                else if (b > 0.0490) nLooseB++;
            }
        }

        std::sort(jets.begin(), jets.end(),
                  [](auto& a, auto& b){ return a.Pt() > b.Pt(); });

        if (!jets.empty()) {
            h_leadJetPt->Fill(jets[0].Pt(), w);
            h_leadJetEta->Fill(jets[0].Eta(), w);
        }

        bool isVBF = false, isGGH = false;

        if (jets.size() >= 2) {
            TLorentzVector j1 = jets[0], j2 = jets[1];
            double mjj = (j1+j2).M();
            double deta = fabs(j1.Eta()-j2.Eta());

            bool passVBF =
                mjj > 400 && deta > 2.5 &&
                j1.Pt() > 35 && j2.Pt() > 25;

            bool passBtag = (nLooseB < 2 && nMedB < 1);

            if (passVBF && passBtag) isVBF = true;
            else if (passBtag) isGGH = true;

            h_dijetPt->Fill((j1+j2).Pt(), w);
            h_dijetMass->Fill(mjj, w);
        }
        else if (jets.size() == 1 && nLooseB < 2)
            isGGH = true;

        if (isVBF) h_mass_VBF->Fill(dimuon.M(), w);
        if (isGGH) h_mass_ggH->Fill(dimuon.M(), w);
    }

    // =========================================================================
    // OUTPUT
    // =========================================================================
    TFile out("output_histos_KIT_full_WZZ-4F_2024.root","RECREATE");
    for (auto h : {h_mass,h_dimuonPt,h_dimuonEta,h_leadJetPt,h_leadJetEta,
                   h_dijetPt,h_dijetMass,h_jetPtCorr,h_mass_VBF,h_mass_ggH})
        h->Write();
    out.Close();

    std::cout << "Finished MC Analysis with KIT" << std::endl;
}
