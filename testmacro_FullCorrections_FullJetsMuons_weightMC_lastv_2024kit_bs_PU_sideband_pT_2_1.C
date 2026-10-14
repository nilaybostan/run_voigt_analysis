// ============================================================================
// Nilay Bostan - CERN - October 2026
//
// FULL ANALYSIS MACRO — MC ONLY (HISTOGRAM VERSION)
// CMS Hmumu Run2024 MCs
//
// SYNCED WITH makeTree_MC_2024.C
//   * Corrections: BSC/Regular muon ScaRe (m_mc / a_mc), muon Medium Prompt ID SF
//     (muon_Z_2024.json, raw NanoAOD pT), PU weight (puWeights_H.json,
//     Collisions24_H_goldenJSON), UParT b-tag SF, JERC + gen-matched JER hook
//     (2.5<|eta|<3.0), JetID, jet veto map, year = 2024.
//   * Event selection = sel_full of the tree macro:
//         MET/noise filters, HLT_IsoMu24, nMuon == 2, opposite sign,
//         no mass window, trigger match (raw pT),
//         muon eta/ID/iso, corrected+FSR pT thresholds (26/20 GeV),
//         event-level jet veto map.
//   * Categories: b-veto on ANALYSIS jets, VBF / ggH (incl. 0/1-jet events),
//     pair definition = kVBFPairMode (same constant as in the DATA macro).
//   * Weight = w_total of the tree macro:
//         genWeight * xsec*lumi/sum(genWeight) * w_PU * w_muID[0]*w_muID[1] * w_btag
//   * Histograms are unchanged w.r.t. the old macro.
//
// NORMALIZATION: sum(genWeight) is taken over ALL input files of ONE call.
// Run each MC sample in a single call.
//
// KEEP IN SYNC: blocks marked [SHARED] are identical in makeTree_MC_2024.C
// and in the DATA histogram macro.
// ============================================================================

#pragma cling add_include_path("/cvmfs/cms.cern.ch/el9_amd64_gcc12/external/py3-correctionlib/2.2.2-120738cfaaf3f7c1056fe67d97e25dac/lib/python3.9/site-packages/correctionlib/include")
#pragma cling add_library_path("/cvmfs/cms.cern.ch/el9_amd64_gcc12/external/py3-correctionlib/2.2.2-120738cfaaf3f7c1056fe67d97e25dac/lib/python3.9/site-packages/correctionlib/lib")
#pragma cling load("correctionlib")

#include "JetCorrections.h"
#include "correction.h"

#include <TChain.h>
#include <TError.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLorentzVector.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

// ============================================================================
// [SHARED] GLOBAL CORRECTIONS
// ============================================================================

JetCorrections jetCorr;

std::shared_ptr<const correction::CorrectionSet> btagCorr;

// RAW correctionlib access to jet_jerc.json for the year/region-specific hooks
std::shared_ptr<const correction::CorrectionSet> jetJercCS;
std::shared_ptr<const correction::Correction> jerPtResolution;
std::shared_ptr<const correction::Correction> jerScaleFactor;
std::shared_ptr<const correction::Correction> l2l3ResidualData;
std::shared_ptr<const correction::Correction> l1FastJetMC;
std::shared_ptr<const correction::Correction> l2RelativeMC;
std::shared_ptr<const correction::Correction> l3AbsoluteMC;

std::shared_ptr<const correction::Correction>
findCorrectionBySuffix(
    const std::shared_ptr<const correction::CorrectionSet> &cs,
    const std::string &suffix) {

  for (auto const &kv : *cs) {
    const std::string &name = kv.first;

    if (name.size() >= suffix.size() &&
        name.compare(name.size() - suffix.size(), suffix.size(), suffix) == 0)
      return cs->at(name);
  }

  std::cerr << "ERROR: no correction ending in \"" << suffix
            << "\" found. Available corrections:" << std::endl;

  for (auto const &kv : *cs)
    std::cerr << "   - " << kv.first << std::endl;

  return nullptr;
}

bool loadJetJercRegionCorrections(const std::string &jetJercJsonFile) {
  try {
    jetJercCS = correction::CorrectionSet::from_file(jetJercJsonFile);

    jerPtResolution =
        findCorrectionBySuffix(jetJercCS, "_JRV2_MC_PtResolution_AK4PFPuppi");
    jerScaleFactor =
        findCorrectionBySuffix(jetJercCS, "_JRV2_MC_ScaleFactor_AK4PFPuppi");
    l2l3ResidualData =
        findCorrectionBySuffix(jetJercCS, "_DATA_L2L3Residual_AK4PFPuppi");
    l1FastJetMC =
        findCorrectionBySuffix(jetJercCS, "_MC_L1FastJet_AK4PFPuppi");
    l2RelativeMC =
        findCorrectionBySuffix(jetJercCS, "_MC_L2Relative_AK4PFPuppi");
    l3AbsoluteMC =
        findCorrectionBySuffix(jetJercCS, "_MC_L3Absolute_AK4PFPuppi");

    if (!jerPtResolution || !jerScaleFactor || !l2l3ResidualData ||
        !l1FastJetMC || !l2RelativeMC || !l3AbsoluteMC)
      return false;

  } catch (const std::exception &e) {
    std::cerr << "ERROR: cannot load jet_jerc JSON \"" << jetJercJsonFile
              << "\": " << e.what() << std::endl;
    return false;
  }

  std::cout << "Year/region-aware jet corrections loaded from: "
            << jetJercJsonFile << std::endl;
  return true;
}

// ============================================================================
// HISTOGRAMS (unchanged)
// ============================================================================

TH1F *h_mass = new TH1F("h_mass","Dimuon mass",100,0,200);
TH1F *h_dimuonPt = new TH1F("h_dimuonPt","Dimuon pT",100,0,200);
TH1F *h_dimuonEta = new TH1F("h_dimuonEta","Dimuon eta",60,-3,3);
TH1F *h_jetPtCorr = new TH1F("h_jetPtCorr","Corrected jet pT",100,0,500);
TH1F *h_leadJetPt = new TH1F("h_leadJetPt","Leading jet pT",100,0,500);
TH1F *h_leadJetEta = new TH1F("h_leadJetEta","Leading jet eta",60,-5,5);
TH1F *h_dijetMass = new TH1F("h_dijetMass","Dijet mass",100,0,2000);
TH1F *h_dijetPt = new TH1F("h_dijetPt","Dijet pT",100,0,1000);
TH1F *h_mass_VBF = new TH1F("h_mass_VBF","VBF dimuon mass",100,0,200);
TH1F *h_mass_ggH = new TH1F("h_mass_ggH","ggH dimuon mass",100,0,200);

TH1F *h_dimuonPt_115_135 = new TH1F("h_dimuonPt_115_135","Dimuon pT (115 < m_{#mu#mu} < 135 GeV)",100,0,200);
TH1F *h_dimuonEta_115_135 = new TH1F("h_dimuonEta_115_135","Dimuon eta (115 < m_{#mu#mu} < 135 GeV)",60,-3,3);
TH1F *h_dimuonPt_115_135_VBF = new TH1F("h_dimuonPt_115_135_VBF","VBF Dimuon pT (115 < m_{#mu#mu} < 135 GeV)",100,0,200);
TH1F *h_dimuonEta_115_135_VBF = new TH1F("h_dimuonEta_115_135_VBF","VBF Dimuon eta (115 < m_{#mu#mu} < 135 GeV)",60,-3,3);
TH1F *h_dimuonPt_115_135_ggH = new TH1F("h_dimuonPt_115_135_ggH","ggH Dimuon pT (115 < m_{#mu#mu} < 135 GeV)",100,0,200);
TH1F *h_dimuonEta_115_135_ggH = new TH1F("h_dimuonEta_115_135_ggH","ggH Dimuon eta (115 < m_{#mu#mu} < 135 GeV)",60,-3,3);

TH1F *h_dimuonPt_VBF = new TH1F("h_dimuonPt_VBF","VBF Dimuon pT",100,0,200);
TH1F *h_dimuonPt_ggH = new TH1F("h_dimuonPt_ggH","ggH Dimuon pT",100,0,200);
TH1F *h_dimuonEta_VBF = new TH1F("h_dimuonEta_VBF","VBF Dimuon eta",60,-3,3);
TH1F *h_dimuonEta_ggH = new TH1F("h_dimuonEta_ggH","ggH Dimuon eta",60,-3,3);
TH1F *h_dijetMass_VBF = new TH1F("h_dijetMass_VBF","VBF Dijet mass",100,0,2000);
TH1F *h_dijetMass_ggH = new TH1F("h_dijetMass_ggH","ggH Dijet mass",100,0,2000);
TH1F *h_dijetPt_VBF = new TH1F("h_dijetPt_VBF","VBF Dijet pT",100,0,1000);
TH1F *h_dijetPt_ggH = new TH1F("h_dijetPt_ggH","ggH Dijet pT",100,0,1000);
TH1F *h_leadJetPt_VBF = new TH1F("h_leadJetPt_VBF","VBF Leading jet pT",100,0,500);
TH1F *h_leadJetPt_ggH = new TH1F("h_leadJetPt_ggH","ggH Leading jet pT",100,0,500);
TH1F *h_leadJetEta_VBF = new TH1F("h_leadJetEta_VBF","VBF Leading jet eta",60,-5,5);
TH1F *h_leadJetEta_ggH = new TH1F("h_leadJetEta_ggH","ggH Leading jet eta",60,-5,5);

// ============================================================================
// [SHARED] ANALYSIS CONSTANTS
// ============================================================================

// Trigger matching on RAW Muon_pt (before ScaRe), IsoMu24 only.
static const double kTrigMuRawPtMin = 26.;
static const double kTrigObjPtMin = 24.;
static const double kTrigMatchDR = 0.1;

// TrigObj_filterBits mask for HLT_IsoMu24. 0 = not applied.
static const ULong64_t kIsoMu24FilterMask = 0;

// Loose baseline (same as the tree macro)
static const bool kRequireMETFilters = true;
static const bool kRequireHLT = true;
static const bool kRequireOS = true;

// Jets with corrected pT above this value are considered (jet-veto-map threshold)
static const double kJetPtStore = 15.;

// VBF jet-pair definition: 0 = two leading analysis jets,
//                          1 = analysis-jet pair with the highest m_jj
// Tree macro stores both; the histograms use this one. SAME value for DATA/MC.
static const int kVBFPairMode = 0;

int gYear = 2024;

double gBtagWP_L = -1.;
double gBtagWP_M = -1.;

// ============================================================================
// [SHARED] SMALL STRUCTS
// ============================================================================

struct MuInfo {
  int nanoIdx = -1; // original NanoAOD Muon index, for FSR association
  float ptNano = 0, ptBS = -1, chi2BS = -1;
  float ptUsed = 0; // pT before the muon correction (BS-constrained if used)
  float ptCorr = 0; // after the muon correction (ScaRe)
  float eta = 0, phi = 0, iso = 0;
  float trigDR = -1;
  float idSF = 1; // MC only
  int charge = 0;
  bool usedBS = false, mediumId = false, trigMatched = false;
};

struct JetInfo {
  TLorentzVector p4; // JERC-corrected
  float ptNano = 0, ptRaw = 0, btag = 0, chEmEF = 0, neEmEF = 0, dRmu = 0;
  float btagSF = 1; // MC only
  int hadFlav = 0;  // MC only
  bool idOK = false, vetoCand = false, inVetoMap = false;
  bool passRegion = false, passDRmu = false, passAnalysis = false;
};

// ============================================================================
// [SHARED] FORWARD-REGION JET RULES
// ============================================================================

bool passRegionJetCut(double corrPt, double eta) {
  double ae = fabs(eta);

  if (gYear <= 2023 && ae >= 3.0 && ae < 5.0 && corrPt < 50.)
    return false;

  if (gYear <= 2024 && ae >= 2.5 && ae < 3.0 && corrPt < 50.)
    return false;

  return true;
}

// ============================================================================
// [SHARED] JET pT CORRECTION (region aware)
//   * MC   : standard JERC + gen-matched JER in 2.5<|eta|<3.0 (<=2024)
//   * DATA : standard JERC + L2L3Residual hook in 2.0<|eta|<2.5 (2024)
// ============================================================================

static long gResidualFail = 0;

double getCorrectedJetPtRegionAware(double pt, double eta, double phi,
                                    double rawFactor, double area, double rho,
                                    unsigned int run, bool isMC,
                                    bool hasGenMatch, double genPt) {

  double corrPt = jetCorr.getCorrectedJetPt(pt, eta, phi, rawFactor, isMC);

  double ae = fabs(eta);

  if (isMC && gYear <= 2024 && ae >= 2.5 && ae < 3.0 && hasGenMatch &&
      genPt > 0.) {

    double reso = jerPtResolution->evaluate({eta, corrPt, rho});
    double sf = jerScaleFactor->evaluate({eta, corrPt});

    if (fabs(corrPt - genPt) < 3.0 * reso * corrPt)
      corrPt = std::max(0.0, genPt + sf * (corrPt - genPt));
  }

  if (!isMC && gYear == 2024 && ae >= 2.0 && ae < 2.5) {

    try {

      double rawPt = pt * (1.0 - rawFactor);

      double ptL1 = rawPt * l1FastJetMC->evaluate({area, eta, rawPt, rho});

      double ptL1L2 = ptL1 * l2RelativeMC->evaluate({eta, phi, ptL1});

      double ptL1L2L3 = ptL1L2 * l3AbsoluteMC->evaluate({eta, ptL1L2});

      if (ptL1L2L3 < 30.0) {

        double residualAt30 =
            l2l3ResidualData->evaluate({(double)run, eta, 30.0});

        corrPt = ptL1L2L3 * residualAt30;
      }

    } catch (const std::exception &e) {

      if (gResidualFail < 3)
        std::cerr << "WARNING: L2L3Residual hook failed (run " << run
                  << "): " << e.what() << " -> keeping standard corrPt"
                  << std::endl;

      gResidualFail++;
    }
  }

  return corrPt;
}

// ============================================================================
// [SHARED] B-TAG WORKING POINTS FROM CORRECTIONLIB
// ============================================================================

bool loadBtagWorkingPoints(
    const std::shared_ptr<const correction::CorrectionSet> &cs) {
  try {
    auto c = cs->at("UParTAK4_wp_values");

    gBtagWP_L = c->evaluate({std::string("L")});

    gBtagWP_M = c->evaluate({std::string("M")});

  } catch (const std::exception &e) {
    std::cerr << "ERROR: cannot read UParTAK4 working points: " << e.what()
              << "\nAvailable corrections in btagging.json:" << std::endl;

    for (auto const &corr : *cs)
      std::cerr << "   - " << corr.first << std::endl;

    return false;
  }

  std::cout << "UParTAK4 WPs from correctionlib: L = " << gBtagWP_L
            << ", M = " << gBtagWP_M << std::endl;

  return true;
}

// ============================================================================
// [SHARED] DIJET PAIR SELECTION
//   mode 0 : two leading-pT analysis jets
//   mode 1 : analysis-jet pair with the highest m_jj
// `jets` sorted by pT (descending), `ana` = indices (into `jets`) of the
// analysis jets in ascending order.
// ============================================================================

bool selectDijet(const std::vector<JetInfo> &jets, const std::vector<int> &ana,
                 int mode, int &i1, int &i2) {
  i1 = -1;
  i2 = -1;

  if (ana.size() < 2)
    return false;

  if (mode == 0) {
    i1 = ana[0];
    i2 = ana[1];
    return true;
  }

  double best = -1.;

  for (size_t a = 0; a < ana.size(); a++) {
    for (size_t b = a + 1; b < ana.size(); b++) {
      double m = (jets[ana[a]].p4 + jets[ana[b]].p4).M();

      if (m > best) {
        best = m;
        i1 = ana[a];
        i2 = ana[b];
      }
    }
  }

  return i1 >= 0;
}

// ============================================================================
// [SHARED] DELTA PHI / DELTA R
// ============================================================================

double deltaPhi(double phi1, double phi2) {

  double dphi = phi1 - phi2;

  while (dphi > M_PI)
    dphi -= 2 * M_PI;

  while (dphi <= -M_PI)
    dphi += 2 * M_PI;

  return dphi;
}

double deltaR(double eta1, double phi1, double eta2, double phi2) {

  double deta = eta1 - eta2;

  double dphi = deltaPhi(phi1, phi2);

  return sqrt(deta * deta + dphi * dphi);
}

// ============================================================================
// [SHARED] FILTER OUT UNREACHABLE / BROKEN INPUT FILES
// ============================================================================

std::vector<std::string>
FilterReadableFiles(const std::vector<std::string> &inputFiles) {

  std::vector<std::string> goodFiles;

  std::vector<std::string> badFiles;

  for (auto const &file : inputFiles) {

    TFile *tf = nullptr;

    const int maxAttempts = 2;

    for (int attempt = 1; attempt <= maxAttempts && (!tf || tf->IsZombie());
         attempt++) {

      if (tf) {
        delete tf;
        tf = nullptr;
      }

      tf = TFile::Open(file.c_str(), "READ");
    }

    if (!tf || tf->IsZombie()) {

      std::cerr << "SKIP (dosya acilamadi / unreachable): " << file
                << std::endl;

      badFiles.push_back(file);

      if (tf)
        delete tf;

      continue;
    }

    tf->Close();

    delete tf;

    goodFiles.push_back(file);
  }

  std::cout << "\nDosya kontrolu tamamlandi: " << goodFiles.size() << " / "
            << inputFiles.size() << " dosya kullanilabilir." << std::endl;

  if (!badFiles.empty()) {

    std::cout << "Atlanan (skip edilen) dosya sayisi: " << badFiles.size()
              << std::endl;

    for (auto const &b : badFiles)
      std::cout << "   - " << b << std::endl;
  }

  return goodFiles;
}

// ============================================================================
// MC-ONLY: MUON CORRECTION, ID SF, PILEUP, B-TAG SF
// ============================================================================

std::shared_ptr<const correction::CorrectionSet> muonCorrBSC;
std::shared_ptr<const correction::Correction> muon_m_mc_BSC;
std::shared_ptr<const correction::Correction> muon_a_mc_BSC;

std::shared_ptr<const correction::CorrectionSet> muonCorrRegular;
std::shared_ptr<const correction::Correction> muon_m_mc_Regular;
std::shared_ptr<const correction::Correction> muon_a_mc_Regular;

std::shared_ptr<const correction::Correction> btagSF;
std::shared_ptr<const correction::Correction> btagLightSF;

std::shared_ptr<const correction::CorrectionSet> puCorr;
std::shared_ptr<const correction::Correction> puWeightCorr;

std::shared_ptr<const correction::CorrectionSet> muonSFCorr;
std::shared_ptr<const correction::Correction> muonIDSF;

double getCorrectedMuonPt_MC(double pt, double eta, double phi, bool isBSC) {

  const std::shared_ptr<const correction::Correction> &m_corr =
      isBSC ? muon_m_mc_BSC : muon_m_mc_Regular;

  const std::shared_ptr<const correction::Correction> &a_corr =
      isBSC ? muon_a_mc_BSC : muon_a_mc_Regular;

  if (!m_corr || !a_corr)
    return pt;

  double m = m_corr->evaluate({eta, phi, "nom"});
  double a = a_corr->evaluate({eta, phi, "nom"});

  return pt * m + a;
}

double getMuonIDSF(double pt, double eta) {

  if (!muonIDSF)
    return 1.0;

  try {
    return muonIDSF->evaluate({eta, pt, "nominal"});
  } catch (const std::exception &e) {
    std::cerr << "WARNING: muon ID SF evaluation failed: " << e.what()
              << ". Using SF=1." << std::endl;
    return 1.0;
  }
}

double getPUWeight(double nTruePU) {

  if (!puWeightCorr)
    return 1.0;

  try {
    return puWeightCorr->evaluate({nTruePU, "nominal"});
  } catch (const std::exception &e) {
    std::cerr << "WARNING: PU correction evaluation failed for nTruePU = "
              << nTruePU << " : " << e.what() << ". Using PU weight = 1."
              << std::endl;
    return 1.0;
  }
}

// ============================================================================
// MAIN FUNCTION
// ============================================================================

void testmacro_FullCorrections_FullJetsMuons_weightMC_lastv_2024kit_bs_PU_sideband_pT_2_1(

    std::vector<std::string> inputFiles,

    double xsec_pb,

    double lumi_fb,

    const std::string &KITDir = "/afs/cern.ch/user/n/nbostan/new_CMS/"
                                "CMSSW_14_0_18/src/RoccoR/post2022E-update/",

    // [ADDED, optional] same defaults as makeTree_MC_2024; the output default
    // is the path that was hard-coded in the old macro.
    const std::string &outFile =
        "/eos/user/n/nbostan/2025_Samples/"
        "FullCorrections_2025KIT_ttbar_bs_PU_MuonSF_SR_2_sideband_pT_2_1_valeriaupdatedcode_less5.root",

    const std::string &BSCDir =
        "/afs/cern.ch/user/f/ferrico/cmsonly/JSON_VXBS/2024/",

    const std::string &RegularDir =
        "/afs/cern.ch/user/f/ferrico/cmsonly/JSON/2024/",

    int year = 2024,

    // empty -> KITDir + "jet_jerc.json"
    const std::string &jetJercJsonFile = ""

) {

  gErrorIgnoreLevel = kError;

  gYear = year;

  // ========================================================================
  // LOAD CORRECTIONS (same as makeTree_MC_2024)
  // ========================================================================

  std::cout << "Loading BSC muon corrections..." << std::endl;

  muonCorrBSC = correction::CorrectionSet::from_file(BSCDir + "schemaV2.json");
  muon_m_mc_BSC = muonCorrBSC->at("m_mc");
  muon_a_mc_BSC = muonCorrBSC->at("a_mc");

  std::cout << "Loading Regular muon corrections..." << std::endl;

  muonCorrRegular =
      correction::CorrectionSet::from_file(RegularDir + "schemaV2.json");
  muon_m_mc_Regular = muonCorrRegular->at("m_mc");
  muon_a_mc_Regular = muonCorrRegular->at("a_mc");

  muonSFCorr = correction::CorrectionSet::from_file(KITDir + "muon_Z_2024.json");
  muonIDSF = muonSFCorr->at("NUM_MediumPromptID_DEN_TrackerMuons");
  std::cout << "Muon ID SF loaded" << std::endl;

  btagCorr = correction::CorrectionSet::from_file(KITDir + "btagging.json");
  btagSF = btagCorr->at("UParTAK4_comb");
  btagLightSF = btagCorr->at("UParTAK4_light");

  if (!loadBtagWorkingPoints(btagCorr))
    return;

  std::cout << "Btag correction loaded" << std::endl;

  jetCorr = JetCorrections(KITDir + "jet_jerc.txt", KITDir + "jetid.txt",
                           KITDir + "jetvetomaps.txt");

  std::cout << "Jet corrections loaded (year rules: " << gYear << ")"
            << std::endl;

  const std::string jercJson =
      jetJercJsonFile.empty() ? (KITDir + "jet_jerc.json") : jetJercJsonFile;

  if (!loadJetJercRegionCorrections(jercJson))
    return;

  puCorr = correction::CorrectionSet::from_file(KITDir + "puWeights_H.json");
  puWeightCorr = puCorr->at("Collisions24_H_goldenJSON");

  std::cout << "PU correction loaded: Collisions24_H_goldenJSON" << std::endl;

  for (auto h :
       {h_mass, h_dimuonPt, h_dimuonEta, h_jetPtCorr, h_leadJetPt,
        h_leadJetEta, h_dijetMass, h_dijetPt, h_mass_VBF, h_mass_ggH,
        h_dimuonPt_115_135, h_dimuonEta_115_135, h_dimuonPt_115_135_VBF,
        h_dimuonEta_115_135_VBF, h_dimuonPt_115_135_ggH,
        h_dimuonEta_115_135_ggH, h_dimuonPt_VBF, h_dimuonPt_ggH,
        h_dimuonEta_VBF, h_dimuonEta_ggH, h_dijetMass_VBF, h_dijetMass_ggH,
        h_dijetPt_VBF, h_dijetPt_ggH, h_leadJetPt_VBF, h_leadJetPt_ggH,
        h_leadJetEta_VBF, h_leadJetEta_ggH})
    h->Sumw2();

  // ========================================================================
  // CHAIN + NORMALIZATION
  // ========================================================================

  std::vector<std::string> readableFiles = FilterReadableFiles(inputFiles);

  if (readableFiles.empty()) {

    std::cerr << "ERROR: Hicbir dosya acilamadi, cikiliyor." << std::endl;

    return;
  }

  TChain chain("Events");

  for (auto const &file : readableFiles)
    chain.Add(file.c_str());

  std::cout << "Total entries = " << chain.GetEntries() << std::endl;

  if (!chain.GetBranch("Jet_btagUParTAK4B")) {

    std::cerr << "ERROR: branch Jet_btagUParTAK4B not found." << std::endl;

    return;
  }

  double lumi_pb = lumi_fb * 1000.0;
  double sumGenWeight = 0.0;

  {
    TTreeReader sumReader(&chain);
    TTreeReaderValue<float> sumGenWeightReader(sumReader, "genWeight");

    while (sumReader.Next())
      sumGenWeight += *sumGenWeightReader;
  }

  std::cout << "Sum genWeight = " << sumGenWeight << std::endl;

  if (fabs(sumGenWeight) < 1e-6) {

    std::cerr << "ERROR: Invalid sumGenWeight" << std::endl;

    return;
  }

  double normFactor = (xsec_pb * lumi_pb) / sumGenWeight;

  std::cout << "MC normalization factor = " << normFactor << std::endl;

  // ========================================================================
  // TREE READER
  // ========================================================================

  TTreeReader reader(&chain);

  TTreeReaderValue<UInt_t> run(reader, "run");

  TTreeReaderValue<float> genWeight(reader, "genWeight");

  TTreeReaderValue<float> *nTruePU = nullptr;

  if (chain.GetBranch("Pileup_nTrueInt")) {
    nTruePU = new TTreeReaderValue<float>(reader, "Pileup_nTrueInt");
    std::cout << "Pileup_nTrueInt branch found" << std::endl;
  } else {
    std::cerr << "WARNING: Pileup_nTrueInt branch not found!" << std::endl;
  }

  TTreeReaderValue<bool> Flag_goodVertices(reader, "Flag_goodVertices");
  TTreeReaderValue<bool> Flag_globalSuperTightHalo2016Filter(
      reader, "Flag_globalSuperTightHalo2016Filter");
  TTreeReaderValue<bool> Flag_EcalDeadCellTriggerPrimitiveFilter(
      reader, "Flag_EcalDeadCellTriggerPrimitiveFilter");
  TTreeReaderValue<bool> Flag_BadPFMuonFilter(reader, "Flag_BadPFMuonFilter");
  TTreeReaderValue<bool> Flag_BadPFMuonDzFilter(reader,
                                                "Flag_BadPFMuonDzFilter");
  TTreeReaderValue<bool> Flag_hfNoisyHitsFilter(reader,
                                                "Flag_hfNoisyHitsFilter");

  TTreeReaderValue<bool> HLT_IsoMu24(reader, "HLT_IsoMu24");

  TTreeReaderValue<int> nMuon(reader, "nMuon");
  TTreeReaderArray<float> Muon_pt(reader, "Muon_pt");
  TTreeReaderArray<float> Muon_eta(reader, "Muon_eta");
  TTreeReaderArray<float> Muon_phi(reader, "Muon_phi");
  TTreeReaderArray<float> Muon_iso(reader, "Muon_pfRelIso04_all");
  TTreeReaderArray<int> Muon_charge(reader, "Muon_charge");
  TTreeReaderArray<bool> Muon_mediumId(reader, "Muon_mediumId");

  // FSR photon collection (optional in the input NanoAOD).
  bool hasFSR = chain.GetBranch("FsrPhoton_pt") &&
                chain.GetBranch("FsrPhoton_eta") &&
                chain.GetBranch("FsrPhoton_phi") &&
                chain.GetBranch("FsrPhoton_dROverEt2") &&
                chain.GetBranch("FsrPhoton_relIso03") &&
                chain.GetBranch("FsrPhoton_muonIdx");

  TTreeReaderArray<float> *FsrPhoton_pt = nullptr;
  TTreeReaderArray<float> *FsrPhoton_eta = nullptr;
  TTreeReaderArray<float> *FsrPhoton_phi = nullptr;
  TTreeReaderArray<float> *FsrPhoton_dROverEt2 = nullptr;
  TTreeReaderArray<float> *FsrPhoton_relIso03 = nullptr;
  TTreeReaderArray<int> *FsrPhoton_muonIdx = nullptr;

  if (hasFSR) {
    FsrPhoton_pt = new TTreeReaderArray<float>(reader, "FsrPhoton_pt");
    FsrPhoton_eta = new TTreeReaderArray<float>(reader, "FsrPhoton_eta");
    FsrPhoton_phi = new TTreeReaderArray<float>(reader, "FsrPhoton_phi");
    FsrPhoton_dROverEt2 =
        new TTreeReaderArray<float>(reader, "FsrPhoton_dROverEt2");
    FsrPhoton_relIso03 =
        new TTreeReaderArray<float>(reader, "FsrPhoton_relIso03");
    FsrPhoton_muonIdx = new TTreeReaderArray<int>(reader, "FsrPhoton_muonIdx");
  }

  std::cout << "FSR photon branches = " << hasFSR << std::endl;

  bool hasBS = false;

  if (chain.GetBranch("Muon_bsConstrainedPt"))
    hasBS = true;

  TTreeReaderArray<float> *Muon_bsConstrainedPt = nullptr;
  TTreeReaderArray<float> *Muon_bsConstrainedChi2 = nullptr;

  if (hasBS) {
    Muon_bsConstrainedPt =
        new TTreeReaderArray<float>(reader, "Muon_bsConstrainedPt");
    Muon_bsConstrainedChi2 =
        new TTreeReaderArray<float>(reader, "Muon_bsConstrainedChi2");
  }

  std::cout << "BeamSpot branch = " << hasBS << std::endl;

  TTreeReaderValue<int> nTrigObj(reader, "nTrigObj");
  TTreeReaderArray<float> TrigObj_pt(reader, "TrigObj_pt");
  TTreeReaderArray<float> TrigObj_eta(reader, "TrigObj_eta");
  TTreeReaderArray<float> TrigObj_phi(reader, "TrigObj_phi");
  TTreeReaderArray<UShort_t> TrigObj_id(reader, "TrigObj_id");

  TTreeReaderArray<ULong64_t> *TrigObj_filterBits = nullptr;

  if (chain.GetBranch("TrigObj_filterBits"))
    TrigObj_filterBits =
        new TTreeReaderArray<ULong64_t>(reader, "TrigObj_filterBits");

  TTreeReaderValue<int> nJet(reader, "nJet");
  TTreeReaderArray<float> Jet_pt(reader, "Jet_pt");
  TTreeReaderArray<float> Jet_eta(reader, "Jet_eta");
  TTreeReaderArray<float> Jet_phi(reader, "Jet_phi");
  TTreeReaderArray<float> Jet_mass(reader, "Jet_mass");
  TTreeReaderArray<float> Jet_rawFactor(reader, "Jet_rawFactor");
  TTreeReaderArray<float> Jet_chEmEF(reader, "Jet_chEmEF");
  TTreeReaderArray<float> Jet_neEmEF(reader, "Jet_neEmEF");
  TTreeReaderArray<float> Jet_area(reader, "Jet_area");
  TTreeReaderValue<float> Rho(reader, "Rho_fixedGridRhoFastjetAll");
  TTreeReaderArray<float> Jet_btagUParTAK4B(reader, "Jet_btagUParTAK4B");
  TTreeReaderArray<unsigned char> Jet_hadronFlavour(reader,
                                                    "Jet_hadronFlavour");

  // gen-jet matching for the 2.5<|eta|<3.0 JER hook (MC only).
  TTreeReaderArray<int> Jet_genJetIdx(reader, "Jet_genJetIdx");
  TTreeReaderArray<float> GenJet_pt(reader, "GenJet_pt");

  // ========================================================================
  // CUTFLOW COUNTERS
  // ========================================================================

  long long nTotal = 0, nMETFilter = 0, nTrigger = 0, nTwoMuon = 0, nOS = 0,
            nTrigMatch = 0, nMuonSel = 0, nJetVeto = 0,
            nFinal = 0;

  // ========================================================================
  // EVENT LOOP
  // ========================================================================

  while (reader.Next()) {

    nTotal++;

    // ---- MET / NOISE FILTERS

    bool passMET = (*Flag_goodVertices) &&
                   (*Flag_globalSuperTightHalo2016Filter) &&
                   (*Flag_EcalDeadCellTriggerPrimitiveFilter) &&
                   (*Flag_BadPFMuonFilter) && (*Flag_BadPFMuonDzFilter) &&
                   (*Flag_hfNoisyHitsFilter);

    if (kRequireMETFilters && !passMET)
      continue;

    nMETFilter++;

    // ---- TRIGGER

    if (kRequireHLT && !(*HLT_IsoMu24))
      continue;

    nTrigger++;

    // ---- EXACTLY TWO MUONS + OPPOSITE SIGN

    if (*nMuon != 2)
      continue;

    nTwoMuon++;

    if (kRequireOS && Muon_charge[0] * Muon_charge[1] >= 0)
      continue;

    nOS++;

    // ====================================================================
    // MUONS: BS-constrained pT choice, ScaRe (MC), trigger matching, ID SF
    // ====================================================================

    MuInfo mu[2];

    for (int i = 0; i < 2; i++) {

      mu[i].nanoIdx = i;
      mu[i].ptNano = Muon_pt[i];
      mu[i].eta = Muon_eta[i];
      mu[i].phi = Muon_phi[i];
      mu[i].iso = Muon_iso[i];
      mu[i].charge = Muon_charge[i];
      mu[i].mediumId = Muon_mediumId[i];

      if (hasBS) {
        mu[i].ptBS = (*Muon_bsConstrainedPt)[i];
        mu[i].chi2BS = (*Muon_bsConstrainedChi2)[i];
      }

      mu[i].usedBS = (hasBS && mu[i].chi2BS < 30.);

      mu[i].ptUsed = mu[i].usedBS ? mu[i].ptBS : mu[i].ptNano;

      mu[i].ptCorr = getCorrectedMuonPt_MC(mu[i].ptUsed, mu[i].eta, mu[i].phi,
                                           mu[i].usedBS);

      // Medium Prompt ID SF, evaluated with the raw NanoAOD muon pT
      mu[i].idSF = getMuonIDSF(mu[i].ptNano, mu[i].eta);

      // trigger matching on RAW Muon_pt (before ScaRe)
      double minDR = -1.;

      for (int t = 0; t < *nTrigObj; t++) {

        if (TrigObj_id[t] != 13)
          continue;

        if (TrigObj_pt[t] < kTrigObjPtMin)
          continue;

        if (kIsoMu24FilterMask != 0 &&
            !(TrigObj_filterBits &&
              ((*TrigObj_filterBits)[t] & kIsoMu24FilterMask)))
          continue;

        double dr = deltaR(mu[i].eta, mu[i].phi, TrigObj_eta[t], TrigObj_phi[t]);

        if (minDR < 0. || dr < minDR)
          minDR = dr;
      }

      mu[i].trigDR = minDR;

      mu[i].trigMatched = (mu[i].ptNano >= kTrigMuRawPtMin && minDR >= 0. &&
                           minDR < kTrigMatchDR);
    }

    // order by corrected pT: index 0 = leading
    if (mu[1].ptCorr > mu[0].ptCorr)
      std::swap(mu[0], mu[1]);

    bool trigMatched = mu[0].trigMatched || mu[1].trigMatched;

    // ---- muon flags

    bool passEta = true, passID = true, passIso = true;

    for (int i = 0; i < 2; i++) {

      if (fabs(mu[i].eta) > 2.4)
        passEta = false;

      if (!mu[i].mediumId)
        passID = false;

      if (mu[i].iso > 0.25)
        passIso = false;
    }

    // ---- DIMUON + FSR PHOTON RECOVERY (FSR used for the pT thresholds only,
    //      exactly as in the tree macro; the mass is the corrected one)

    TLorentzVector mu1, mu2;

    mu1.SetPtEtaPhiM(mu[0].ptCorr, mu[0].eta, mu[0].phi, 0.105);
    mu2.SetPtEtaPhiM(mu[1].ptCorr, mu[1].eta, mu[1].phi, 0.105);

    TLorentzVector mu1FSR = mu1, mu2FSR = mu2;
    TLorentzVector dimuon = mu1 + mu2;

    // FSR photon requirements: pT > 2 GeV, |eta| < 2.4, relIso03 < 1.8,
    // dROverEt2 < 0.012, DeltaR(mu,gamma) < 0.5; best = smallest dROverEt2.
    double bestDROverEt2[2] = {1.e9, 1.e9};
    TLorentzVector bestPhoton[2];
    bool muHasFSR[2] = {false, false};

    if (hasFSR) {
      for (int ig = 0; ig < (int)FsrPhoton_pt->GetSize(); ++ig) {
        const double gpt = (*FsrPhoton_pt)[ig];
        const double geta = (*FsrPhoton_eta)[ig];
        const double gphi = (*FsrPhoton_phi)[ig];
        const double drOverEt2 = (*FsrPhoton_dROverEt2)[ig];
        const double relIso = (*FsrPhoton_relIso03)[ig];
        const int muonIdx = (*FsrPhoton_muonIdx)[ig];

        if (gpt <= 2.0 || fabs(geta) >= 2.4 || relIso >= 1.8 ||
            drOverEt2 >= 0.012)
          continue;

        int im = -1;
        for (int imSorted = 0; imSorted < 2; ++imSorted) {
          if (mu[imSorted].nanoIdx == muonIdx) {
            im = imSorted;
            break;
          }
        }
        if (im < 0)
          continue;

        if (deltaR(mu[im].eta, mu[im].phi, geta, gphi) >= 0.5)
          continue;

        if (drOverEt2 < bestDROverEt2[im]) {
          bestDROverEt2[im] = drOverEt2;
          bestPhoton[im].SetPtEtaPhiM(gpt, geta, gphi, 0.0);
          muHasFSR[im] = true;
        }
      }
    }

    if (muHasFSR[0])
      mu1FSR += bestPhoton[0];

    if (muHasFSR[1])
      mu2FSR += bestPhoton[1];

    double mll = dimuon.M();

    // corrected+FSR pT thresholds
    bool passPt = (mu1FSR.Pt() >= 26.) && (mu2FSR.Pt() >= 20.);
    bool passMuonSel = passEta && passID && passIso && passPt;

    // ---- sel_full part 1: trigger match, part 2: muon selection

    if (!trigMatched)
      continue;

    nTrigMatch++;

    if (!passMuonSel)
      continue;

    nMuonSel++;

    // ====================================================================
    // JETS: JERC (+ gen-matched JER hook) + event-level jet veto + analysis jets
    // ====================================================================

    std::vector<JetInfo> jets;

    bool vetoEvent = false;

    for (int j = 0; j < *nJet; j++) {

      // gen-jet match (for the 2.5<|eta|<3.0 JER hook)
      int gi = Jet_genJetIdx[j];
      bool hasGen = (gi >= 0 && gi < (int)GenJet_pt.GetSize());
      double genPt = hasGen ? GenJet_pt[gi] : -1.;

      double corrJetPt = getCorrectedJetPtRegionAware(
          Jet_pt[j], Jet_eta[j], Jet_phi[j], Jet_rawFactor[j], Jet_area[j],
          *Rho, *run,
          true, // MC
          hasGen, genPt);

      if (corrJetPt < kJetPtStore)
        continue;

      JetInfo ji;

      ji.p4.SetPtEtaPhiM(corrJetPt, Jet_eta[j], Jet_phi[j], Jet_mass[j]);
      ji.ptNano = Jet_pt[j];
      ji.ptRaw = Jet_pt[j] * (1.0 - Jet_rawFactor[j]);
      ji.btag = Jet_btagUParTAK4B[j];
      ji.chEmEF = Jet_chEmEF[j];
      ji.neEmEF = Jet_neEmEF[j];
      ji.hadFlav = Jet_hadronFlavour[j];

      ji.idOK = jetCorr.passJetID(ji.p4);

      ji.vetoCand = ji.idOK && (Jet_chEmEF[j] + Jet_neEmEF[j]) < 0.9;

      ji.inVetoMap = false;

      if (ji.vetoCand && !jetCorr.passJetVetoMap(Jet_eta[j], Jet_phi[j])) {
        ji.inVetoMap = true;
        vetoEvent = true;
      }

      ji.dRmu = std::min(ji.p4.DeltaR(mu1), ji.p4.DeltaR(mu2));

      ji.passDRmu = (ji.dRmu >= 0.4);

      ji.passRegion = passRegionJetCut(corrJetPt, Jet_eta[j]);

      ji.passAnalysis = (corrJetPt >= 20.) && (fabs(Jet_eta[j]) <= 4.7) &&
                        ji.idOK && ji.passRegion && ji.passDRmu;

      jets.push_back(ji);
    }

    // ---- sel_full part 3: event-level jet veto map

    if (vetoEvent)
      continue;

    nJetVeto++;

    std::sort(jets.begin(), jets.end(),
              [](const JetInfo &a, const JetInfo &b) {
                return a.p4.Pt() > b.p4.Pt();
              });

    // ---- analysis jets, b-tag counting, b-tag SF (analysis jets, |eta| < 2.5)

    std::vector<int> ana;

    double btagWeight = 1.0;
    int nMediumB = 0, nLooseB = 0;

    for (size_t k = 0; k < jets.size(); k++) {

      JetInfo &ji = jets[k];

      if (!ji.passAnalysis)
        continue;

      ana.push_back((int)k);

      double absEta = fabs(ji.p4.Eta());

      if (absEta >= 2.5)
        continue;

      if (ji.btag > gBtagWP_M)
        nMediumB++;

      if (ji.btag > gBtagWP_L)
        nLooseB++;

      double sf = 1.0;

      if (ji.hadFlav == 5 || ji.hadFlav == 4) {
        sf = btagSF->evaluate({"central", "M", ji.hadFlav, absEta, ji.p4.Pt()});
      } else if (ji.hadFlav == 0) {
        sf = btagLightSF->evaluate(
            {"central", "M", ji.hadFlav, absEta, ji.p4.Pt()});
      }

      ji.btagSF = sf;
      btagWeight *= sf;
    }

    bool passBVeto = (nMediumB == 0 && nLooseB <= 1);

    // ====================================================================
    // MC WEIGHT  (= w_total of the tree macro)
    //   w_gen * w_norm * w_pu * w_muID[0] * w_muID[1] * w_btag
    // ====================================================================

    double wGen = *genWeight;
    double wPU = 1.0;

    if (nTruePU)
      wPU = getPUWeight(**nTruePU);

    double wMuTot = (double)mu[0].idSF * (double)mu[1].idSF;

    double weight = wGen * normFactor * wPU * wMuTot * btagWeight;

    // ====================================================================
    // DIMUON HISTOGRAMS (after sel_full, full MC weight)
    // ====================================================================

    h_mass->Fill(mll, weight);
    h_dimuonPt->Fill(dimuon.Pt(), weight);
    h_dimuonEta->Fill(dimuon.Eta(), weight);

    bool inSideband = (mll > 115. && mll < 135.);

    if (inSideband) {
      h_dimuonPt_115_135->Fill(dimuon.Pt(), weight);
      h_dimuonEta_115_135->Fill(dimuon.Eta(), weight);
    }

    for (int k : ana)
      h_jetPtCorr->Fill(jets[k].p4.Pt(), weight);

    // ---- LEADING (ANALYSIS) JET

    if (!ana.empty()) {
      h_leadJetPt->Fill(jets[ana[0]].p4.Pt(), weight);
      h_leadJetEta->Fill(jets[ana[0]].p4.Eta(), weight);
    }

    // ====================================================================
    // VBF / ggH CATEGORY
    // ggH = every non-VBF event passing the b-veto, incl. 0/1-jet events
    // ====================================================================

    bool isVBF = false;
    bool isGGH = false;

    int i1 = -1, i2 = -1;

    bool hasPair = selectDijet(jets, ana, kVBFPairMode, i1, i2);

    double mjj = 0., ptjj = 0.;

    if (hasPair) {

      const TLorentzVector &j1 = jets[i1].p4;
      const TLorentzVector &j2 = jets[i2].p4;

      mjj = (j1 + j2).M();
      ptjj = (j1 + j2).Pt();

      double deta = fabs(j1.Eta() - j2.Eta());

      h_dijetMass->Fill(mjj, weight);
      h_dijetPt->Fill(ptjj, weight);

      if (passBVeto && mjj > 400. && deta > 2.5 && j1.Pt() > 35. &&
          j2.Pt() > 25.)
        isVBF = true;
    }

    if (passBVeto && !isVBF)
      isGGH = true;

    // ---- category histograms

    if (isVBF) {

      nFinal++;

      h_mass_VBF->Fill(mll, weight);
      h_dimuonPt_VBF->Fill(dimuon.Pt(), weight);
      h_dimuonEta_VBF->Fill(dimuon.Eta(), weight);

      h_dijetMass_VBF->Fill(mjj, weight);
      h_dijetPt_VBF->Fill(ptjj, weight);

      if (!ana.empty()) {
        h_leadJetPt_VBF->Fill(jets[ana[0]].p4.Pt(), weight);
        h_leadJetEta_VBF->Fill(jets[ana[0]].p4.Eta(), weight);
      }

      if (inSideband) {
        h_dimuonPt_115_135_VBF->Fill(dimuon.Pt(), weight);
        h_dimuonEta_115_135_VBF->Fill(dimuon.Eta(), weight);
      }
    }

    if (isGGH) {

      nFinal++;

      h_mass_ggH->Fill(mll, weight);
      h_dimuonPt_ggH->Fill(dimuon.Pt(), weight);
      h_dimuonEta_ggH->Fill(dimuon.Eta(), weight);

      if (hasPair) {
        h_dijetMass_ggH->Fill(mjj, weight);
        h_dijetPt_ggH->Fill(ptjj, weight);
      }

      if (!ana.empty()) {
        h_leadJetPt_ggH->Fill(jets[ana[0]].p4.Pt(), weight);
        h_leadJetEta_ggH->Fill(jets[ana[0]].p4.Eta(), weight);
      }

      if (inSideband) {
        h_dimuonPt_115_135_ggH->Fill(dimuon.Pt(), weight);
        h_dimuonEta_115_135_ggH->Fill(dimuon.Eta(), weight);
      }
    }

  } // while(reader.Next())

  std::cout << "\nCutflow:"
            << "\n  total                  = " << nTotal
            << "\n  MET filters            = " << nMETFilter
            << "\n  HLT_IsoMu24            = " << nTrigger
            << "\n  nMuon == 2             = " << nTwoMuon
            << "\n  opposite sign          = " << nOS
            << "\n  trigger match          = " << nTrigMatch
            << "\n  muon selection         = " << nMuonSel
            << "\n  jet veto map (sel_full)= " << nJetVeto
            << "\n  VBF + ggH              = " << nFinal << std::endl;

  // ========================================================================
  // OUTPUT ROOT FILE
  // ========================================================================

  TFile *out = new TFile(outFile.c_str(), "RECREATE");

  if (out->IsZombie()) {
    std::cerr << "ERROR: Cannot create output ROOT file" << std::endl;
    return;
  }

  h_mass->Write();
  h_dimuonPt->Write();
  h_dimuonEta->Write();

  h_jetPtCorr->Write();
  h_leadJetPt->Write();
  h_leadJetEta->Write();

  h_dijetMass->Write();
  h_dijetPt->Write();

  h_mass_VBF->Write();
  h_mass_ggH->Write();

  h_dimuonPt_115_135->Write();
  h_dimuonEta_115_135->Write();
  h_dimuonPt_115_135_VBF->Write();
  h_dimuonEta_115_135_VBF->Write();
  h_dimuonPt_115_135_ggH->Write();
  h_dimuonEta_115_135_ggH->Write();
  h_dimuonPt_VBF->Write();
  h_dimuonPt_ggH->Write();
  h_dimuonEta_VBF->Write();
  h_dimuonEta_ggH->Write();
  h_dijetMass_VBF->Write();
  h_dijetMass_ggH->Write();
  h_dijetPt_VBF->Write();
  h_dijetPt_ggH->Write();
  h_leadJetPt_VBF->Write();
  h_leadJetPt_ggH->Write();
  h_leadJetEta_VBF->Write();
  h_leadJetEta_ggH->Write();

  out->Write();
  out->Close();

  if (TrigObj_filterBits)
    delete TrigObj_filterBits;

  if (Muon_bsConstrainedPt)
    delete Muon_bsConstrainedPt;

  if (Muon_bsConstrainedChi2)
    delete Muon_bsConstrainedChi2;

  if (nTruePU)
    delete nTruePU;

  if (FsrPhoton_pt) {
    delete FsrPhoton_pt;
    delete FsrPhoton_eta;
    delete FsrPhoton_phi;
    delete FsrPhoton_dROverEt2;
    delete FsrPhoton_relIso03;
    delete FsrPhoton_muonIdx;
  }

  std::cout << "Output written" << std::endl;

  std::cout << "\n====================================\n"
            << " MC ANALYSIS FINISHED (synced with makeTree_MC_2024)\n"
            << " Output: " << outFile << "\n"
            << "====================================\n";
}
