// ============================================================================
// FULL ANALYSIS MACRO — DATA ONLY
// NO correctionlib, SAME SELECTIONS
// + BEAM-SPOT MUON INFORMATION ADDED (NO OTHER CHANGE)
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
#include <TVector2.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

// ============================================================================
// STRUCTURES & HELPERS (UNCHANGED)
// ============================================================================
struct MuonScale {
    double ptMin, ptMax;
    double etaMin, etaMax;
    double scale;
};

std::vector<MuonScale> readMuonScale(const std::string& filename) {
    std::vector<MuonScale> table;
    std::ifstream f(filename);
    if (!f.is_open()) return table;

    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        MuonScale b;
        std::stringstream ss(line);
        ss >> b.ptMin >> b.ptMax >> b.etaMin >> b.etaMax >> b.scale;
        table.push_back(b);
    }
    return table;
}

// ============================================================================
// KIT MUON CORRECTION — DATA RETURNS RAW PT
// ============================================================================
double getCorrectedMuonPt(
    double pt, double /*eta*/,
    const std::vector<MuonScale>& /*scaleTable*/
) {
    return pt;
}


double deltaPhi(double phi1, double phi2)
{
    double dphi = phi1 - phi2;

    while (dphi > M_PI)
        dphi -= 2.0 * M_PI;

    while (dphi <= -M_PI)
        dphi += 2.0 * M_PI;

    return dphi;
}

double deltaR(
    double eta1, double phi1,
    double eta2, double phi2
)
{
    double deta = eta1 - eta2;
    double dphi = deltaPhi(phi1, phi2);

    return std::sqrt(deta * deta + dphi * dphi);
}

// ============================================================================
// MAIN ANALYSIS — DATA ONLY
// ============================================================================
void testmacro_FullCorrections_FullJetsMuons_2025kit_bs(
    std::vector<std::string> inputFiles,
    bool isData = true,
    const std::string& KITDir =
        "/afs/cern.ch/user/n/nbostan/new_CMS/CMSSW_14_0_18/"
        "src/RoccoR/post2022E-update/"
) {
    gErrorIgnoreLevel = kError;

    // =========================================================================
    // KIT SCALE FILE (READ BUT UNUSED FOR DATA)
    // =========================================================================
    std::vector<MuonScale> muonScale =
        readMuonScale(KITDir + "muon_Z_2025.txt");

    // =========================================================================
    // JET CORRECTIONS
    // =========================================================================
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
TH1F *h_dimuonPt    = new TH1F("h_dimuonPt","Dimuon pT",100,0,200);
TH1F *h_dimuonEta   = new TH1F("h_dimuonEta","Dimuon eta",50,-2.5,2.5);
TH1F *h_leadJetPt   = new TH1F("h_leadJetPt","Leading jet pT",100,0,200);
TH1F *h_leadJetEta  = new TH1F("h_leadJetEta","Leading jet eta",50,-5,5);
TH1F *h_dijetPt     = new TH1F("h_dijetPt","Dijet pT",100,0,1000);
TH1F *h_dijetMass   = new TH1F("h_dijetMass","Dijet mass",100,0,2000);
TH1F *h_jetPtCorr   = new TH1F("h_jetPtCorr","Corrected jet pT",100,0,1000);
TH1F *h_mass_VBF    = new TH1F("h_mass_VBF","VBF dimuon mass",100,0,200);
TH1F *h_mass_ggH    = new TH1F("h_mass_ggH","ggH dimuon mass",100,0,200);

// ADD THIS HERE
for (auto h : {
    h_mass,
    h_dimuonPt,
    h_dimuonEta,
    h_leadJetPt,
    h_leadJetEta,
    h_dijetPt,
    h_dijetMass,
    h_jetPtCorr,
    h_mass_VBF,
    h_mass_ggH
}) {
    h->Sumw2();
}

    // =========================================================================
    // BUILD TCHAIN
    // =========================================================================
    TChain chain("Events");
    for (auto& f : inputFiles) chain.Add(f.c_str());

    // =========================================================================
    // TREE READER
    // =========================================================================
    TTreeReader reader(&chain);

    TTreeReaderValue<bool>  HLT_IsoMu24(reader,"HLT_IsoMu24");

    TTreeReaderValue<int> nMuon(reader,"nMuon");
    TTreeReaderArray<float> Muon_pt(reader,"Muon_pt");
    TTreeReaderArray<float> Muon_eta(reader,"Muon_eta");
    TTreeReaderArray<float> Muon_phi(reader,"Muon_phi");
    TTreeReaderArray<float> Muon_iso(reader,"Muon_pfRelIso04_all");
    TTreeReaderArray<int>   Muon_charge(reader,"Muon_charge");
    TTreeReaderArray<bool>  Muon_mediumID(reader,"Muon_mediumId");
    
// =====================================================================
// TRIGGER OBJECTS
// =====================================================================

TTreeReaderValue<int> nTrigObj(reader, "nTrigObj");
TTreeReaderArray<float> TrigObj_pt(reader, "TrigObj_pt");
TTreeReaderArray<float> TrigObj_eta(reader, "TrigObj_eta");
TTreeReaderArray<float> TrigObj_phi(reader, "TrigObj_phi");

TTreeReaderArray<UShort_t> TrigObj_id(reader, "TrigObj_id");

// Remove this if unused
// TTreeReaderArray<ULong64_t> TrigObj_filterBits(reader, "TrigObj_filterBits");
//TTreeReaderArray<unsigned int> TrigObj_filterBits(reader, "TrigObj_filterBits");

    // ===================== BEAM-SPOT ADDITION =====================
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
    else if (chain.GetBranch("Jet_btagDeepB"))
        Jet_btag = new TTreeReaderArray<float>(reader,"Jet_btagDeepB");

    // =========================================================================
    // EVENT LOOP — DATA
    // =========================================================================
    while (reader.Next()) {

        if (!(*HLT_IsoMu24)) continue;
        if (*nMuon != 2) continue;

        double w = 1.0;

        // ===================== BEAM-SPOT PT =====================
        std::vector<float> pt_bs(*nMuon);
        for (int i = 0; i < *nMuon; i++) {
            if (Muon_bsConstrainedChi2[i] < 30)
                pt_bs[i] = Muon_bsConstrainedPt[i];
            else
                pt_bs[i] = Muon_pt[i];
        }
        // ========================================================

// =====================================================================
// MUON BASELINE SELECTION
// BOTH MUONS: pT >= 20, |eta| < 2.4, Medium ID, Loose Iso
// =====================================================================

double ptCorr[2];
bool passMuonSelection = true;

for (int i = 0; i < 2; i++) {

    ptCorr[i] = getCorrectedMuonPt(
        pt_bs[i],
        Muon_eta[i],
        muonScale
    );

    if (ptCorr[i] < 20.0)
        passMuonSelection = false;

    if (fabs(Muon_eta[i]) > 2.4)
        passMuonSelection = false;

    if (!Muon_mediumID[i])
        passMuonSelection = false;

    if (Muon_iso[i] > 0.25)
        passMuonSelection = false;
}

if (!passMuonSelection)
    continue;


// =====================================================================
// TAG MUON: pT >= 26 GeV
// =====================================================================

int tagIdx = -1;
double bestPt = -1.0;

for (int i = 0; i < 2; i++) {

    if (ptCorr[i] < 26.0)
        continue;

    if (ptCorr[i] > bestPt) {
        bestPt = ptCorr[i];
        tagIdx = i;
    }
}

if (tagIdx < 0)
    continue;

int probeIdx = 1 - tagIdx;


// =====================================================================
// TRIGGER OBJECT MATCHING
// Require tag muon matched to a trigger muon within DeltaR < 0.1
// =====================================================================

bool triggerMatched = false;

for (int itrig = 0; itrig < *nTrigObj; itrig++) {

    // Muon trigger object
    if (TrigObj_id[itrig] != 13)
        continue;

    double dR = deltaR(
        Muon_eta[tagIdx],
        Muon_phi[tagIdx],
        TrigObj_eta[itrig],
        TrigObj_phi[itrig]
    );

    if (dR < 0.1) {
        triggerMatched = true;
        break;
    }
}

if (!triggerMatched)
    continue;

// =====================================================================
// OPPOSITE SIGN
// =====================================================================

if (Muon_charge[tagIdx] * Muon_charge[probeIdx] >= 0)
    continue;

double pt1 = ptCorr[tagIdx];
double pt2 = ptCorr[probeIdx];

        TLorentzVector m1, m2;
        m1.SetPtEtaPhiM(pt1, Muon_eta[tagIdx], Muon_phi[tagIdx], 0.105);
        m2.SetPtEtaPhiM(pt2, Muon_eta[probeIdx], Muon_phi[probeIdx], 0.105);

        TLorentzVector dimuon = m1 + m2;

        h_mass->Fill(dimuon.M(), w);
        h_dimuonPt->Fill(dimuon.Pt(), w);
        h_dimuonEta->Fill(dimuon.Eta(), w);

        // ========================= JETS (UNCHANGED) =========================
        std::vector<TLorentzVector> jets;
        int nLooseB = 0, nMedB = 0;

        for (int j = 0; j < *nJet; j++) {

            double corrPt = jetCorr.getCorrectedJetPt(
                Jet_pt[j], Jet_eta[j], Jet_phi[j], Jet_raw[j], false
            );

            if (corrPt < 20.0) continue;
	    if (fabs(Jet_eta[j]) > 4.7) continue;

TLorentzVector jet;
jet.SetPtEtaPhiM(
    corrPt,
    Jet_eta[j],
    Jet_phi[j],
    Jet_mass[j]
);

// Jet-muon overlap veto
if (jet.DeltaR(m1) <= 0.4)
    continue;

if (jet.DeltaR(m2) <= 0.4)
    continue;

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
            TLorentzVector j1 = jets[0];
            TLorentzVector j2 = jets[1];

            double mjj  = (j1 + j2).M();
            double deta = fabs(j1.Eta() - j2.Eta());

            bool passVBF =
                (mjj > 400) &&
                (deta > 2.5) &&
                (j1.Pt() > 35) &&
                (j2.Pt() > 25);

            bool passBtagVeto = (nLooseB < 2 && nMedB < 1);

            if (passVBF && passBtagVeto) isVBF = true;
            else if (passBtagVeto)       isGGH = true;

            h_dijetPt->Fill((j1+j2).Pt(), w);
            h_dijetMass->Fill(mjj, w);

        } else if (jets.size() == 1) {
            if (nLooseB < 2) isGGH = true;
        }

        if (isVBF) h_mass_VBF->Fill(dimuon.M(), w);
        if (isGGH) h_mass_ggH->Fill(dimuon.M(), w);
    }

    // =========================================================================
    // OUTPUT
    // =========================================================================
    TFile out("/eos/user/n/nbostan/2025_Samples/output_histos_DATA_2025_F_KIT_bs_test_upd.root","RECREATE");

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

    std::cout << " Finished: DATA-only macro with beam-spot muons added" << std::endl;
}
