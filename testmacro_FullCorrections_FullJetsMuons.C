#include "RoccoR.h"

#include <TRandom.h>    
#include <TChain.h>     

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
void testmacro_FullCorrections_FullJetsMuons(std::vector<std::string> inputFiles,
                                             bool isData = true,
                                             const std::string &puFileName = "/afs/cern.ch/user/n/nbostan/new_CMS/CMSSW_14_0_18/src/RoccoR/data/puWeights.txt",
                                             const std::string &muPOGDir   = "/afs/cern.ch/user/n/nbostan/new_CMS/CMSSW_14_0_18/src/RoccoR/data/POG/MUO/")
{
    bool isMC = !isData;

    // === Initialize Rochester correction ===
    RoccoR rc;
    rc.init("/afs/cern.ch/user/n/nbostan/new_CMS/CMSSW_14_0_18/src/RoccoR/data/RoccoR2022.txt");

    // === Load Muon POG SF tables ===
    std::vector<MuonSF> muonJPsi       = readMuonSF(muPOGDir + "muon_JPsi.txt");
    std::vector<MuonSF> muonZ          = readMuonSF(muPOGDir + "muon_Z.txt");
    std::vector<MuonSF> muonHighPt     = readMuonSF(muPOGDir + "muon_HighPt.txt");
    std::vector<MuonSF> muonScaleSmear = readMuonSF(muPOGDir + "muon_scalesmearing.txt");

    // === Load PU weights ===
    std::vector<PUWeight> puWeights;
    if (isMC) puWeights = readPUWeights(puFileName);

    // === Load jet corrections ===
    JetCorrections jetCorr(
        "/afs/cern.ch/user/n/nbostan/new_CMS/CMSSW_14_0_18/src/RoccoR/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/jet_jerc.txt",
        "/afs/cern.ch/user/n/nbostan/new_CMS/CMSSW_14_0_18/src/RoccoR/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/jetid.txt",
        "/afs/cern.ch/user/n/nbostan/new_CMS/CMSSW_14_0_18/src/RoccoR/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/jetvetomaps.txt"
    );

    // === Histograms ===
    TH1F *h_mass      = new TH1F("h_mass",      "Dimuon invariant mass;M_{#mu#mu} [GeV];Events", 80, 70, 110);
    TH1F *h_dimuonPt  = new TH1F("h_dimuonPt",  "Dimuon p_{T};p_{T}^{#mu#mu} [GeV];Events", 100, 0, 500);
    TH1F *h_dijetPt   = new TH1F("h_dijetPt",   "Dijet p_{T};p_{T}^{jj} [GeV];Events", 100, 0, 1000);
    TH1F *h_dijetMass = new TH1F("h_dijetMass", "Dijet invariant mass;M_{jj} [GeV];Events", 100, 0, 2000);
    TH1F *h_jetPtCorr = new TH1F("h_jetPtCorr", "Corrected jet p_{T};p_{T}^{corr} [GeV];Jets", 100, 0, 1000);

    // === 2D histograms ===
    TH2F *h_mass_vs_pt = new TH2F("h_mass_vs_pt",
                                  "Dimuon invariant mass vs p_{T}^{#mu#mu};M_{#mu#mu} [GeV];p_{T}^{#mu#mu} [GeV]",
                                  80, 70, 110, 100, 0, 500);

    TH2F *h_dijetMass_vs_pt = new TH2F("h_dijetMass_vs_pt",
                                       "Dijet invariant mass vs p_{T}^{jj};M_{jj} [GeV];p_{T}^{jj} [GeV]",
                                       100, 0, 2000, 100, 0, 1000);

    // === Setup tree reader ===
    TChain chain("Events");
    for (auto &f : inputFiles) chain.Add(f.c_str());
    TTreeReader reader(&chain);

    // === Define mandatory branches ===
    TTreeReaderValue<Int_t> nMuon(reader, "nMuon");
    TTreeReaderArray<Float_t> Muon_pt(reader, "Muon_pt");
    TTreeReaderArray<Float_t> Muon_eta(reader, "Muon_eta");
    TTreeReaderArray<Float_t> Muon_phi(reader, "Muon_phi");
    TTreeReaderArray<Int_t> Muon_charge(reader, "Muon_charge");
    TTreeReaderArray<Bool_t> Muon_mediumId(reader, "Muon_mediumId");
    TTreeReaderArray<Float_t> Muon_pfRelIso04_all(reader, "Muon_pfRelIso04_all");
    TTreeReaderValue<Bool_t> HLT_IsoMu24(reader, "HLT_IsoMu24");

    // Optional branches
    TTreeReaderArray<Float_t> *Muon_genPt = nullptr;
    TTreeReaderValue<Float_t> *Pileup_nTrueInt = nullptr;

    bool hasGenPt = chain.GetBranch("Muon_genPt") != nullptr;
    bool hasPileup = chain.GetBranch("Pileup_nTrueInt") != nullptr;

    if (hasGenPt) Muon_genPt = new TTreeReaderArray<Float_t>(reader, "Muon_genPt");
    else std::cerr << "Warning: Muon_genPt branch not found, skipping gen-based smearing.\n";

    if (hasPileup && isMC) Pileup_nTrueInt = new TTreeReaderValue<Float_t>(reader, "Pileup_nTrueInt");
    else if (isMC) std::cerr << "Warning: Pileup_nTrueInt branch not found, using weight=1.\n";

    // === Jet branches ===
    TTreeReaderValue<Int_t> nJet(reader, "nJet");
    TTreeReaderArray<Float_t> Jet_pt(reader, "Jet_pt");
    TTreeReaderArray<Float_t> Jet_eta(reader, "Jet_eta");
    TTreeReaderArray<Float_t> Jet_phi(reader, "Jet_phi");
    TTreeReaderArray<Float_t> Jet_mass(reader, "Jet_mass");
    TTreeReaderArray<Float_t> Jet_rawFactor(reader, "Jet_rawFactor");

    // ------------------------
    // Event loop
    // ------------------------
    while (reader.Next()) {
        if (!(*HLT_IsoMu24) || *nMuon < 2) continue;

        double puWeight = 1.0;
        if (isMC && hasPileup && Pileup_nTrueInt)
            puWeight = getPUWeight(**Pileup_nTrueInt, puWeights);

        // --- Tag muon ---
        int tagIdx = -1; float maxPt = -1;
        for (int i = 0; i < *nMuon; i++) {
            if (!Muon_mediumId[i]) continue; 
            if (Muon_pt[i] < 26 || fabs(Muon_eta[i]) > 2.4 || Muon_pfRelIso04_all[i] > 0.15) continue;
            if (Muon_pt[i] > maxPt) { maxPt = Muon_pt[i]; tagIdx = i; }
        }
        if (tagIdx < 0) continue;

        // --- Probe muon ---
        for (int p = 0; p < *nMuon; p++) {
            if (p == tagIdx) continue;
            if (Muon_charge[tagIdx] * Muon_charge[p] >= 0) continue;

            TLorentzVector tag, probe;
            tag.SetPtEtaPhiM(Muon_pt[tagIdx], Muon_eta[tagIdx], Muon_phi[tagIdx], 0.105);
            probe.SetPtEtaPhiM(Muon_pt[p], Muon_eta[p], Muon_phi[p], 0.105);

            // --- Apply RoccoR safely ---
            double sf_tag = 1.0, sf_probe = 1.0;
            if (isMC) {
                float genPt_tag   = (hasGenPt ? (*Muon_genPt)[tagIdx] : -1);
                float genPt_probe = (hasGenPt ? (*Muon_genPt)[p] : -1);

                sf_tag = (genPt_tag > 0)
                    ? rc.kSpreadMC(Muon_charge[tagIdx], Muon_pt[tagIdx], Muon_eta[tagIdx], Muon_phi[tagIdx], genPt_tag, 0, 0)
                    : rc.kScaleMC(Muon_charge[tagIdx], Muon_pt[tagIdx], Muon_eta[tagIdx], Muon_phi[tagIdx], 0, 0);

                sf_probe = (genPt_probe > 0)
                    ? rc.kSpreadMC(Muon_charge[p], Muon_pt[p], Muon_eta[p], Muon_phi[p], genPt_probe, 0, 0)
                    : rc.kScaleMC(Muon_charge[p], Muon_pt[p], Muon_eta[p], Muon_phi[p], 0, 0);
            } else {
                sf_tag   = rc.kScaleDT(Muon_charge[tagIdx], Muon_pt[tagIdx], Muon_eta[tagIdx], Muon_phi[tagIdx], 0, 0);
                sf_probe = rc.kScaleDT(Muon_charge[p], Muon_pt[p], Muon_eta[p], Muon_phi[p], 0, 0);
            }

            // --- Combined Muon SFs ---
            double sf_tag_pog   = getMuonSF(Muon_pt[tagIdx], Muon_eta[tagIdx], muonZ)
                                * getMuonSF(Muon_pt[tagIdx], Muon_eta[tagIdx], muonJPsi)
                                * getMuonSF(Muon_pt[tagIdx], Muon_eta[tagIdx], muonHighPt)
                                * getMuonSF(Muon_pt[tagIdx], Muon_eta[tagIdx], muonScaleSmear);

            double sf_probe_pog = getMuonSF(Muon_pt[p], Muon_eta[p], muonZ)
                                * getMuonSF(Muon_pt[p], Muon_eta[p], muonJPsi)
                                * getMuonSF(Muon_pt[p], Muon_eta[p], muonHighPt)
                                * getMuonSF(Muon_pt[p], Muon_eta[p], muonScaleSmear);

            double totalWeight = puWeight * sf_tag * sf_probe * sf_tag_pog * sf_probe_pog;

            tag.SetPtEtaPhiM(Muon_pt[tagIdx]*sf_tag, Muon_eta[tagIdx], Muon_phi[tagIdx], 0.105);
            probe.SetPtEtaPhiM(Muon_pt[p]*sf_probe, Muon_eta[p], Muon_phi[p], 0.105);

            TLorentzVector dimuon = tag + probe;
            h_mass->Fill(dimuon.M(), totalWeight);
            h_dimuonPt->Fill(dimuon.Pt(), totalWeight);
            h_mass_vs_pt->Fill(dimuon.M(), dimuon.Pt(), totalWeight);
        }

        // --- Jet loop ---
        std::vector<TLorentzVector> jets;
        for (int j = 0; j < *nJet; j++) {
            double corrPt = jetCorr.getCorrectedJetPt(Jet_pt[j], Jet_eta[j], Jet_phi[j], Jet_rawFactor[j], isMC);
            if (corrPt < 25 || fabs(Jet_eta[j]) > 4.7) continue;
            TLorentzVector jet; jet.SetPtEtaPhiM(corrPt, Jet_eta[j], Jet_phi[j], Jet_mass[j]);
            jets.push_back(jet);
            h_jetPtCorr->Fill(corrPt, puWeight);
        }

        if (jets.size() >= 2) {
            TLorentzVector dijet = jets[0] + jets[1];
            h_dijetPt->Fill(dijet.Pt(), puWeight);
            h_dijetMass->Fill(dijet.M(), puWeight);
            h_dijetMass_vs_pt->Fill(dijet.M(), dijet.Pt(), puWeight);
        }
    }

    // === Save ===
    std::string outFileName = isMC ? "output_FullCorrections_MC.root" : "output_FullCorrections_Data.root";
    TFile outFile(outFileName.c_str(),"RECREATE");
    h_mass->Write();
    h_dimuonPt->Write();
    h_dijetPt->Write();
    h_dijetMass->Write();
    h_jetPtCorr->Write();
    h_mass_vs_pt->Write();
    h_dijetMass_vs_pt->Write(); // NEW
    outFile.Close();
    std::cout << "Output written to " << outFileName << std::endl;
}
