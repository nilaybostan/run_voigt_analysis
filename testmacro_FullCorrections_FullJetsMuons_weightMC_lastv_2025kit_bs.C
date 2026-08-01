// ============================================================================
// FULL ANALYSIS MACRO — ROOT ONLY, KIT TEXT MUON CORRECTIONS
// + BEAM-SPOT MUONS ADDED (NO PHYSICS CHANGES)
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
void testmacro_FullCorrections_FullJetsMuons_weightMC_lastv_2025kit_bs(
    std::vector<std::string> inputFiles,
    bool isData = true
) {
    gErrorIgnoreLevel = kError;
    const bool isMC = !isData;

    double lumi_fb = 26.50;
    double xsec_pb = 1.39;

    // =========================================================================
    // INPUT TABLES
    // =========================================================================
    auto muonScale = readMuonScale("/afs/cern.ch/user/n/nbostan/new_CMS/CMSSW_14_0_18/src/RoccoR/post2022E-update/schemaV2_2025.txt");
    auto muonSF_Z  = readMuonSF("/afs/cern.ch/user/n/nbostan/new_CMS/CMSSW_14_0_18/src/RoccoR/post2022E-update/muon_Z_2025.txt");
    auto puWeights = readPUWeights("/afs/cern.ch/user/n/nbostan/new_CMS/CMSSW_14_0_18/src/RoccoR/post2022E-update/puWeights_2025pp_Golden_Summer24_25ns_69200ub.txt");

    JetCorrections jetCorr(
        "/afs/cern.ch/user/n/nbostan/new_CMS/CMSSW_14_0_18/"
        "src/RoccoR/post2022E-update/jet_jerc.txt",
         "/afs/cern.ch/user/n/nbostan/new_CMS/CMSSW_14_0_18/"
        "src/RoccoR/post2022E-update/jetvetomaps.txt",
        "/afs/cern.ch/user/n/nbostan/new_CMS/CMSSW_14_0_18/"
        "src/RoccoR/post2022E-update/jetid.txt"
    );

    // =========================================================================
    // HISTOGRAMS (UNCHANGED)
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
    // MC NORMALIZATION (UNCHANGED)
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

    // =================== ADDED: BEAM-SPOT MUONS ===================
    TTreeReaderArray<float> Muon_bsConstrainedPt(reader,"Muon_bsConstrainedPt");
    TTreeReaderArray<float> Muon_bsConstrainedChi2(reader,"Muon_bsConstrainedChi2");
    // ===============================================================

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
    // EVENT LOOP
    // =========================================================================
    while (reader.Next()) {

        if (!(*HLT_IsoMu24)) continue;
        if (*nMuon != 2) continue;
        //if (*nMuon < 2) continue;

        // ----------- ADDED: beam-spot pT construction -----------
        std::vector<double> pt_bs(*nMuon);
        for (int i = 0; i < *nMuon; i++) {
            if (Muon_bsConstrainedChi2[i] < 30)
                pt_bs[i] = Muon_bsConstrainedPt[i];
            else
                pt_bs[i] = Muon_pt[i];
        }
        // ---------------------------------------------------------

        double w = 1.0;
        if (isMC) {
            w = (*genWeight) * normFactor;
            if (nTruePU) w *= getPUWeight(**nTruePU, puWeights);
        }

// =====================================================================
// MUON SELECTION
// =====================================================================

double ptCorr[2];

bool passBothMuons = true;

for (int i = 0; i < 2; i++) {

    // Corrected pT using beam-spot constrained pT when valid
    ptCorr[i] = getCorrectedMuonPt(
        pt_bs[i],
        Muon_eta[i],
        isMC,
        muonScale
    );

    // Both signal muons must satisfy the baseline selection
    if (ptCorr[i] < 20.0)
        passBothMuons = false;

    if (fabs(Muon_eta[i]) > 2.4)
        passBothMuons = false;

    if (!Muon_mediumID[i])
        passBothMuons = false;

    if (Muon_iso[i] > 0.25)
        passBothMuons = false;
}

// Reject event if either muon fails baseline selection
if (!passBothMuons)
    continue;


// =====================================================================
// TAG MUON: pT > 26 GeV
// =====================================================================

int tag = -1;
double bestPt = -1.0;

for (int i = 0; i < 2; i++) {

    if (ptCorr[i] < 26.0)
        continue;

    // If both muons pass pT > 26, take the highest-pT one as tag
    if (ptCorr[i] > bestPt) {
        bestPt = ptCorr[i];
        tag = i;
    }
}

if (tag < 0)
    continue;

int probe = 1 - tag;


// =====================================================================
// OPPOSITE SIGN
// =====================================================================

if (Muon_charge[tag] * Muon_charge[probe] >= 0)
    continue;


double pt1 = ptCorr[tag];
double pt2 = ptCorr[probe];

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
            if (corrPt < 20.0 || fabs(Jet_eta[j]) > 4.7) continue;

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
    // =======================================================================
    TFile out("/eos/user/n/nbostan/2025_Samples/output_histos_ZZto4L_2025kit_bs_test_upd.root","RECREATE");
    for (auto h : {h_mass,h_dimuonPt,h_dimuonEta,h_leadJetPt,h_leadJetEta,
                   h_dijetPt,h_dijetMass,h_jetPtCorr,h_mass_VBF,h_mass_ggH})
        h->Write();
    out.Close();

    std::cout << " Finished — beam-spot muons added, physics unchanged" << std::endl;
}
