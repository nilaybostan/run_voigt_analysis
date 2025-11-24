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
void testmacro_FullCorrections_FullJetsMuons_weightMC(std::vector<std::string> inputFiles,
                                             bool isData = true,
                                             const std::string &puFileName = "/afs/cern.ch/user/n/nbostan/new_CMS/CMSSW_14_0_18/src/RoccoR/data/puWeights.txt",
                                             const std::string &muPOGDir   = "/afs/cern.ch/user/n/nbostan/new_CMS/CMSSW_14_0_18/src/RoccoR/data/POG/MUO/")
{
    bool isMC = !isData;
    double lumi_fb = 26.6717;
    //double xsec_pb = 0.00105742; // Run 3 VBF cross section
    
    //double xsec_pb = 0.0135096; // Run 3 ggH cross section
    
    //double xsec_pb = 359.05; // Run 3 dyTo2L_M-50_2jets cross section
    
    //double xsec_pb = 98.04; // Run 3 TTto2L2Nu cross section
    
    //double xsec_pb = 1.421; // Run 3 EWK cross section
    
    //double xsec_pb = 7.568;  // Run 3 WZ_2l2q cross section
    
    //double xsec_pb = 6.788; //Run 3 ZZ_2l2q cross section
    
     //double xsec_pb = 4.924; //Run 3 WZ_3lnu cross section
     
     //double xsec_pb = 1.39; //Run 3 ZZ_4l  cross section
     
     double xsec_pb = 1.031; //Run 3 ZZ_2l2nu  cross section

    gErrorIgnoreLevel = kError; // suppress harmless ROOT warnings

    // Initialize RoccoR (muon momentum corrections)
    RoccoR rc;
    rc.init("/afs/cern.ch/user/n/nbostan/new_CMS/CMSSW_14_0_18/src/RoccoR/data/RoccoR2022.txt");

    // Read muon scale/SF files (optional)
    std::vector<MuonSF> muonJPsi       = readMuonSF(muPOGDir + "muon_JPsi.txt");
    std::vector<MuonSF> muonZ          = readMuonSF(muPOGDir + "muon_Z.txt");
    std::vector<MuonSF> muonHighPt     = readMuonSF(muPOGDir + "muon_HighPt.txt");
    std::vector<MuonSF> muonScaleSmear = readMuonSF(muPOGDir + "muon_scalesmearing.txt");

    std::vector<PUWeight> puWeights;
    if (isMC) puWeights = readPUWeights(puFileName);

    // Jet corrections
    JetCorrections jetCorr(
        "/afs/cern.ch/user/n/nbostan/new_CMS/CMSSW_14_0_18/src/RoccoR/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/jet_jerc.txt",
        "/afs/cern.ch/user/n/nbostan/new_CMS/CMSSW_14_0_18/src/RoccoR/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/jetid.txt",
        "/afs/cern.ch/user/n/nbostan/new_CMS/CMSSW_14_0_18/src/RoccoR/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/jetvetomaps.txt"
    );

    // --------------------------------------------------
    // Histograms (inclusive + categories)
    // --------------------------------------------------
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

    // --- Compute MC normalization factor ---
    double normFactor = 1.0;
    if (isMC) {
        double lumi_pb = lumi_fb * 1000.0; // convert fb^-1 to pb^-1
        double nGenEvents = 0;
        TTreeReader readerCount(&chain);
        TTreeReaderValue<Float_t> genW(readerCount, "genWeight");
        while (readerCount.Next()) {
            nGenEvents += (*genW >= 0) ? 1.0 : -1.0;
        }
        normFactor = xsec_pb * lumi_pb / nGenEvents;
        std::cout << "Total signed generated events: " << nGenEvents << std::endl;
        std::cout << "MC normalization factor applied: " << normFactor << std::endl;
    }

    TTreeReader reader(&chain);
    
    TTreeReaderValue<Float_t> genWeight(reader, "genWeight");
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

    TTreeReaderArray<Float_t> *Jet_btag = nullptr;
    if (chain.GetBranch("Jet_btagDeepB"))
        Jet_btag = new TTreeReaderArray<Float_t>(reader, "Jet_btagDeepB");
    else if (chain.GetBranch("Jet_btagDeepFlavB"))
        Jet_btag = new TTreeReaderArray<Float_t>(reader, "Jet_btagDeepFlavB");
    else
        std::cerr << "⚠️ Warning: No known b-tag branch found!" << std::endl;

    TTreeReaderValue<Float_t> *Pileup_nTrueInt = nullptr;
    if (isMC && chain.GetBranch("Pileup_nTrueInt"))
        Pileup_nTrueInt = new TTreeReaderValue<Float_t>(reader, "Pileup_nTrueInt");

    int nEventPassed = 0;

    // =====================================================
    // Event loop
    // =====================================================
    while (reader.Next()) {
        if (!(*HLT_IsoMu24) || *nMuon < 2) continue;
        nEventPassed++;

        double puWeight = 1.0;
        if (isMC && Pileup_nTrueInt)
            puWeight = getPUWeight(**Pileup_nTrueInt, puWeights);

        double evWeight = 1.0;
        if (isMC) evWeight = ((*genWeight >= 0) ? 1.0 : -1.0);

        double totalWeight = evWeight * puWeight;
        if (isMC) totalWeight *= normFactor;

        // --- Muon selection: tag ---
        int tagIdx = -1; float maxPt = -1;
        for (int i = 0; i < *nMuon; i++) {
            if (!Muon_mediumId[i]) continue;
            if (Muon_pt[i] < 26 || fabs(Muon_eta[i]) > 2.4 || Muon_pfRelIso04_all[i] > 0.15) continue;
            if (Muon_pt[i] > maxPt) { maxPt = Muon_pt[i]; tagIdx = i; }
        }
        if (tagIdx < 0) continue;

        // --- Find dimuon pairs ---
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

            h_mass->Fill(dimuon.M(), totalWeight);
            h_dimuonPt->Fill(dimuon.Pt(), totalWeight);
            h_dimuonEta->Fill(dimuon.Eta(), totalWeight);
        }
        if (!foundPair) continue;

        // --- Jets ---
        std::vector<TLorentzVector> jets;
        int nLooseBtag = 0, nMediumBtag = 0;

        for (int j = 0; j < *nJet; j++) {
            double corrPt = jetCorr.getCorrectedJetPt(Jet_pt[j], Jet_eta[j], Jet_phi[j], Jet_rawFactor[j], isMC);
            if (corrPt < 25 || fabs(Jet_eta[j]) > 4.7) continue;

            TLorentzVector jet; jet.SetPtEtaPhiM(corrPt, Jet_eta[j], Jet_phi[j], Jet_mass[j]);
            jets.push_back(jet);
            h_jetPtCorr->Fill(corrPt, totalWeight);

            if (Jet_btag) {
                float bval = (*Jet_btag)[j];
                if (bval > 0.0490) nLooseBtag++;
                if (bval > 0.2783) nMediumBtag++;
            }
        }

        if (!jets.empty()) {
            std::sort(jets.begin(), jets.end(), [](const TLorentzVector &a, const TLorentzVector &b){ return a.Pt() > b.Pt(); });
            h_leadJetPt->Fill(jets[0].Pt(), totalWeight);
            h_leadJetEta->Fill(jets[0].Eta(), totalWeight);
        }

        // --- Categorization ---
        bool isVBF = false, isGGH = false;

        if (jets.size() >= 2) {
            TLorentzVector j1 = jets[0], j2 = jets[1];
            double mjj = (j1 + j2).M();
            double deta = fabs(j1.Eta() - j2.Eta());
            if (mjj > 400 && deta > 2.5 && j1.Pt() > 35 && j2.Pt() > 25 && nLooseBtag < 2 && nMediumBtag < 1) isVBF = true;
            else isGGH = true;

            // Fill VBF histograms
            if (isVBF) {
                h_mass_VBF->Fill(dimuon.M(), totalWeight);
                h_dimuonPt_VBF->Fill(dimuon.Pt(), totalWeight);
                h_dimuonEta_VBF->Fill(dimuon.Eta(), totalWeight);
                h_leadJetPt_VBF->Fill(j1.Pt(), totalWeight);
                h_leadJetEta_VBF->Fill(j1.Eta(), totalWeight);
                h_dijetPt_VBF->Fill((j1+j2).Pt(), totalWeight);
                h_dijetMass_VBF->Fill(mjj, totalWeight);
            }

            // Fill ggH histograms
            if (isGGH) {
                h_mass_ggH->Fill(dimuon.M(), totalWeight);
                h_dimuonPt_ggH->Fill(dimuon.Pt(), totalWeight);
                h_dimuonEta_ggH->Fill(dimuon.Eta(), totalWeight);
                h_leadJetPt_ggH->Fill(j1.Pt(), totalWeight);
                h_leadJetEta_ggH->Fill(j1.Eta(), totalWeight);
                h_dijetPt_ggH->Fill((j1+j2).Pt(), totalWeight);
                h_dijetMass_ggH->Fill(mjj, totalWeight);
            }
        }
    } // end event loop

    std::cout << "Events passed selection: " << nEventPassed << std::endl;

    // --------------------------------------------------
    // Save output histograms
    // --------------------------------------------------
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
    h_dimuonPt_VBF->Write();
    h_dimuonEta_VBF->Write();
    h_leadJetPt_VBF->Write();
    h_leadJetEta_VBF->Write();
    h_dijetPt_VBF->Write();
    h_dijetMass_VBF->Write();

    h_mass_ggH->Write();
    h_dimuonPt_ggH->Write();
    h_dimuonEta_ggH->Write();
    h_leadJetPt_ggH->Write();
    h_leadJetEta_ggH->Write();
    h_dijetPt_ggH->Write();
    h_dijetMass_ggH->Write();

    out.Close();
    std::cout << "Histograms saved to output_histos.root" << std::endl;
}
