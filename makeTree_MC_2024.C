// ============================================================================
// Nilay Bostan 
// TREE-MAKER MACRO — MC ONLY
// CMS Hmumu Run2024
//
// Based on testmacro_FullCorrections_FullJetsMuons_weightMC_lastv_2025kit_bs_PU_sideband_pT_2_1.C
//
// KEEP IN SYNC with makeTree_DATA_2025.C: every block marked [SHARED] below is
// identical in the DATA and MC macros. Edit one, copy to the other.
//
// WHAT CHANGED w.r.t. the histogram macro (following the Pisa/Filippo sync):
//   * Output is a flat TTree ("HmmTree"), NOT a set of histograms.
//   * LOOSE BASELINE only:
//       MET/noise filters, HLT_IsoMu24, nMuon == 2, opposite sign,
//       50 < m_mumu(corrected) < 200 GeV.
//     (switchable through the kRequire* constants)
//   * Everything that used to be a hard cut (trigger matching, muon
//     eta/ID/iso/pT, jet veto map, b-veto, VBF/ggH, 115-135 window) is now
//     stored as a FLAG. The old selection is reproduced with
//         sel_full && cat_isVBF[mode]   /   sel_full && cat_isGGH[mode]
//   * Both VBF jet-pair definitions are stored (index 0 = two leading jets,
//     index 1 = pair with highest m_jj). Use the SAME index for DATA and MC.
//   * ALL MC weights are stored separately:
//       w_gen, w_norm, w_pu, w_muID[2] / w_muID_tot, w_btag, w_total
//     w_total = w_gen * w_norm * w_pu * w_muID_tot * w_btag
//
// JET CORRECTIONS: the MC macro now calls the same getCorrectedJetPtRegionAware
// as the DATA macro (isMC = true), so the 2.5<|eta|<3.0 gen-matched JER hook is
// applied. The L2L3Residual hook (2.0<|eta|<2.5) is DATA-only by construction.
//
// NORMALIZATION: w_norm = xsec*lumi / sum(genWeight) over ALL input files of
// ONE call. Run each MC sample in a single call (do not split it in pieces).
//
// Muon branches:  mu_*[2], index 0 = leading, 1 = subleading (corrected pT)
// Jet branches:   jet_* vectors, all jets with corrected pT > 15 GeV,
//                 sorted by corrected pT (descending)
//
// NOTE: no FSR correction exists in this code, so "uncorrected" = before
// the muon (ScaRe) correction only (BS-constrained pT if the BS fit was used).
//
// PU input: puWeights_H.json, correction Collisions24_H_goldenJSON (nominal)
// ============================================================================

#pragma cling add_include_path("/cvmfs/cms.cern.ch/el9_amd64_gcc12/external/py3-correctionlib/2.2.2-120738cfaaf3f7c1056fe67d97e25dac/lib/python3.9/site-packages/correctionlib/include")
#pragma cling add_library_path("/cvmfs/cms.cern.ch/el9_amd64_gcc12/external/py3-correctionlib/2.2.2-120738cfaaf3f7c1056fe67d97e25dac/lib/python3.9/site-packages/correctionlib/lib")
#pragma cling load("correctionlib")

#include "JetCorrections.h"
#include "correction.h"

#include <TChain.h>
#include <TError.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TTree.h>
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
// [SHARED] ANALYSIS CONSTANTS
// ============================================================================

// Trigger matching on RAW Muon_pt (before ScaRe), IsoMu24 only.
static const double kTrigMuRawPtMin = 26.;
static const double kTrigObjPtMin = 24.;
static const double kTrigMatchDR = 0.1;

// TrigObj_filterBits mask for HLT_IsoMu24. 0 = not applied
// (Run2025 NanoAOD filter bit needs validation).
static const ULong64_t kIsoMu24FilterMask = 0;

// LOOSE BASELINE (everything else is stored as a flag)
static const double kMassMin = 50.;
static const double kMassMax = 200.;
static const bool kRequireMETFilters = true;
static const bool kRequireHLT = true;
static const bool kRequireOS = true;

// Jets with corrected pT above this value are stored in the tree
// (15 GeV = threshold of the jet-veto-map definition).
static const double kJetPtStore = 15.;

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
// The two hooks are physically data- or MC-specific, but they live in ONE
// function so both macros always call the same code path.
// ============================================================================

// Counts failed 2.0<|eta|<2.5 L2L3Residual evaluations (run outside the
// run range of jet_jerc.json). Printed at the end of the job.
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

      // STOPGAP: run outside the JSON range -> keep the standard corrPt.
      // The right fix is a jet_jerc.json that covers the processed runs.
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
//
//   mode 0 : two leading-pT analysis jets
//   mode 1 : analysis-jet pair with the highest m_jj
//
// `jets` is sorted by pT (descending), `ana` holds the indices (into `jets`)
// of the analysis jets in ascending order. i1 < i2 are returned as indices
// into `jets`, so jets[i1] is the higher-pT jet of the pair.
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

void makeTree_MC_2024(

    std::vector<std::string> inputFiles,

    double xsec_pb,

    double lumi_fb,

    const std::string &KITDir = "/afs/cern.ch/user/n/nbostan/new_CMS/"
                                "CMSSW_14_0_18/src/RoccoR/post2022E-update/",

    const std::string &outFile =
        "/eos/user/n/nbostan/2025_Samples/tree_MC_2024_loose.root",

    const std::string &BSCDir =
        "/afs/cern.ch/user/f/ferrico/cmsonly/JSON_VXBS/2024/",

    const std::string &RegularDir =
        "/afs/cern.ch/user/f/ferrico/cmsonly/JSON/2024/",

    int year = 2024,

    // empty -> KITDir + "jet_jerc.json" (same default as the DATA macro)
    const std::string &jetJercJsonFile = ""

) {

  gErrorIgnoreLevel = kError;

  gYear = year;

  // ========================================================================
  // LOAD CORRECTIONS
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
  TTreeReaderValue<UInt_t> luminosityBlock(reader, "luminosityBlock");
  TTreeReaderValue<ULong64_t> eventNumber(reader, "event");

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
  // NanoAOD stores Jet_genJetIdx as Int_t; if ROOT complains about the type,
  // check chain.GetBranch("Jet_genJetIdx")->GetTitle() and adapt.
  TTreeReaderArray<int> Jet_genJetIdx(reader, "Jet_genJetIdx");
  TTreeReaderArray<float> GenJet_pt(reader, "GenJet_pt");

  // ========================================================================
  // OUTPUT FILE + TREE
  // ========================================================================

  TFile out(outFile.c_str(), "RECREATE");

  if (out.IsZombie()) {

    std::cerr << "ERROR: Cannot create output ROOT file: " << outFile
              << std::endl;

    return;
  }

  out.cd();

  TTree *tree = new TTree("HmmTree", "H->mumu loose-baseline tree (MC)");

  // ---- event variables
  UInt_t b_run = 0, b_lumi = 0;
  ULong64_t b_event = 0;

  Bool_t b_passMETFilters = false, b_HLT = false, b_trigMatched = false;

  Bool_t b_passMuonEta = false, b_passMuonID = false, b_passMuonIso = false,
         b_passMuonPt = false, b_passMuonSel = false;

  Bool_t b_passJetVetoMap = false, b_passBVeto = false, b_selFull = false;

  Bool_t b_inWindow_115_135 = false;
  Int_t b_massRegion = -1; // 0: m<115, 1: 115<=m<=135, 2: m>135

  Bool_t b_cat_isVBF[2] = {false, false}; // [0] two leading jets, [1] max mjj
  Bool_t b_cat_isGGH[2] = {false, false};

  Int_t b_nJetsStored = 0, b_nJetsAna = 0, b_nMediumB = 0, b_nLooseB = 0;

  // ---- MC weights
  Float_t b_w_gen = 1, b_w_norm = 1, b_w_pu = 1, b_w_muID[2] = {1, 1},
          b_w_muID_tot = 1, b_w_btag = 1, b_w_total = 1, b_nTruePU = -1;

  // ---- muons (index 0 = leading, 1 = subleading, by corrected pT)
  Float_t b_mu_pt_nano[2], b_mu_pt_bs[2], b_mu_chi2_bs[2], b_mu_pt_uncorr[2],
      b_mu_pt_corr[2], b_mu_eta[2], b_mu_phi[2], b_mu_iso[2], b_mu_trigDR[2];
  Int_t b_mu_charge[2];
  Bool_t b_mu_usedBS[2], b_mu_mediumId[2], b_mu_trigMatched[2];

  // ---- dimuon
  Float_t b_mm_mass_corr = 0, b_mm_mass_uncorr = 0, b_mm_pt_corr = 0,
          b_mm_pt_uncorr = 0, b_mm_eta_corr = 0, b_mm_phi_corr = 0;
  Float_t b_mm_mass_fsr = 0, b_mm_pt_fsr = 0;
  Float_t b_mu_fsr_pt[2] = {0, 0};
  Bool_t b_mu_hasFSR[2] = {false, false};
  Int_t b_nFSR = 0;

  // ---- dijet (index = pair mode: 0 two leading, 1 highest mjj)
  Bool_t b_jj_has[2];
  Float_t b_jj_mass[2], b_jj_pt[2], b_jj_eta[2], b_jj_deta[2];
  Int_t b_jj_idx1[2], b_jj_idx2[2];

  // ---- jets
  std::vector<float> v_jet_pt, v_jet_pt_nano, v_jet_pt_raw, v_jet_eta,
      v_jet_phi, v_jet_mass, v_jet_btag, v_jet_chEmEF, v_jet_neEmEF,
      v_jet_dRmu, v_jet_btagSF;
  std::vector<int> v_jet_hadFlav, v_jet_idOK, v_jet_vetoCand, v_jet_inVetoMap,
      v_jet_passRegion, v_jet_passDRmu, v_jet_passAnalysis;

  tree->Branch("run", &b_run, "run/i");
  tree->Branch("lumi", &b_lumi, "lumi/i");
  tree->Branch("event", &b_event, "event/l");

  tree->Branch("passMETFilters", &b_passMETFilters, "passMETFilters/O");
  tree->Branch("HLT_IsoMu24", &b_HLT, "HLT_IsoMu24/O");
  tree->Branch("trigMatched", &b_trigMatched, "trigMatched/O");

  tree->Branch("passMuonEta", &b_passMuonEta, "passMuonEta/O");
  tree->Branch("passMuonID", &b_passMuonID, "passMuonID/O");
  tree->Branch("passMuonIso", &b_passMuonIso, "passMuonIso/O");
  tree->Branch("passMuonPt", &b_passMuonPt, "passMuonPt/O");
  tree->Branch("passMuonSel", &b_passMuonSel, "passMuonSel/O");

  tree->Branch("passJetVetoMap", &b_passJetVetoMap, "passJetVetoMap/O");
  tree->Branch("passBVeto", &b_passBVeto, "passBVeto/O");
  tree->Branch("sel_full", &b_selFull, "sel_full/O");

  tree->Branch("inWindow_115_135", &b_inWindow_115_135, "inWindow_115_135/O");
  tree->Branch("massRegion", &b_massRegion, "massRegion/I");

  tree->Branch("cat_isVBF", b_cat_isVBF, "cat_isVBF[2]/O");
  tree->Branch("cat_isGGH", b_cat_isGGH, "cat_isGGH[2]/O");

  tree->Branch("nJetsStored", &b_nJetsStored, "nJetsStored/I");
  tree->Branch("nJetsAna", &b_nJetsAna, "nJetsAna/I");
  tree->Branch("nMediumB", &b_nMediumB, "nMediumB/I");
  tree->Branch("nLooseB", &b_nLooseB, "nLooseB/I");

  tree->Branch("w_gen", &b_w_gen, "w_gen/F");
  tree->Branch("w_norm", &b_w_norm, "w_norm/F");
  tree->Branch("w_pu", &b_w_pu, "w_pu/F");
  tree->Branch("w_muID", b_w_muID, "w_muID[2]/F");
  tree->Branch("w_muID_tot", &b_w_muID_tot, "w_muID_tot/F");
  tree->Branch("w_btag", &b_w_btag, "w_btag/F");
  tree->Branch("w_total", &b_w_total, "w_total/F");
  tree->Branch("nTruePU", &b_nTruePU, "nTruePU/F");

  tree->Branch("mu_pt_nano", b_mu_pt_nano, "mu_pt_nano[2]/F");
  tree->Branch("mu_pt_bs", b_mu_pt_bs, "mu_pt_bs[2]/F");
  tree->Branch("mu_chi2_bs", b_mu_chi2_bs, "mu_chi2_bs[2]/F");
  tree->Branch("mu_usedBS", b_mu_usedBS, "mu_usedBS[2]/O");
  tree->Branch("mu_pt_uncorr", b_mu_pt_uncorr, "mu_pt_uncorr[2]/F");
  tree->Branch("mu_pt_corr", b_mu_pt_corr, "mu_pt_corr[2]/F");
  tree->Branch("mu_eta", b_mu_eta, "mu_eta[2]/F");
  tree->Branch("mu_phi", b_mu_phi, "mu_phi[2]/F");
  tree->Branch("mu_charge", b_mu_charge, "mu_charge[2]/I");
  tree->Branch("mu_mediumId", b_mu_mediumId, "mu_mediumId[2]/O");
  tree->Branch("mu_iso", b_mu_iso, "mu_iso[2]/F");
  tree->Branch("mu_trigDR", b_mu_trigDR, "mu_trigDR[2]/F");
  tree->Branch("mu_trigMatched", b_mu_trigMatched, "mu_trigMatched[2]/O");

  tree->Branch("mm_mass_corr", &b_mm_mass_corr, "mm_mass_corr/F");
  tree->Branch("mm_mass_uncorr", &b_mm_mass_uncorr, "mm_mass_uncorr/F");
  tree->Branch("mm_pt_corr", &b_mm_pt_corr, "mm_pt_corr/F");
  tree->Branch("mm_pt_uncorr", &b_mm_pt_uncorr, "mm_pt_uncorr/F");
  tree->Branch("mm_eta_corr", &b_mm_eta_corr, "mm_eta_corr/F");
  tree->Branch("mm_phi_corr", &b_mm_phi_corr, "mm_phi_corr/F");
  tree->Branch("mm_mass_fsr", &b_mm_mass_fsr, "mm_mass_fsr/F");
  tree->Branch("mm_pt_fsr", &b_mm_pt_fsr, "mm_pt_fsr/F");
  tree->Branch("mu_fsr_pt", b_mu_fsr_pt, "mu_fsr_pt[2]/F");
  tree->Branch("mu_hasFSR", b_mu_hasFSR, "mu_hasFSR[2]/O");
  tree->Branch("nFSR", &b_nFSR, "nFSR/I");

  tree->Branch("jj_has", b_jj_has, "jj_has[2]/O");
  tree->Branch("jj_mass", b_jj_mass, "jj_mass[2]/F");
  tree->Branch("jj_pt", b_jj_pt, "jj_pt[2]/F");
  tree->Branch("jj_eta", b_jj_eta, "jj_eta[2]/F");
  tree->Branch("jj_deta", b_jj_deta, "jj_deta[2]/F");
  tree->Branch("jj_idx1", b_jj_idx1, "jj_idx1[2]/I");
  tree->Branch("jj_idx2", b_jj_idx2, "jj_idx2[2]/I");

  tree->Branch("jet_pt", &v_jet_pt);           // JERC-corrected
  tree->Branch("jet_pt_nano", &v_jet_pt_nano); // NanoAOD Jet_pt
  tree->Branch("jet_pt_raw", &v_jet_pt_raw);   // uncorrected: pt*(1-rawFactor)
  tree->Branch("jet_eta", &v_jet_eta);
  tree->Branch("jet_phi", &v_jet_phi);
  tree->Branch("jet_mass", &v_jet_mass);
  tree->Branch("jet_btag", &v_jet_btag);
  tree->Branch("jet_btagSF", &v_jet_btagSF); // 1 unless analysis jet |eta|<2.5
  tree->Branch("jet_hadFlav", &v_jet_hadFlav);
  tree->Branch("jet_chEmEF", &v_jet_chEmEF);
  tree->Branch("jet_neEmEF", &v_jet_neEmEF);
  tree->Branch("jet_dRmu", &v_jet_dRmu); // min dR to the two muons
  tree->Branch("jet_idOK", &v_jet_idOK);
  tree->Branch("jet_vetoCand", &v_jet_vetoCand);   // idOK && chEm+neEm < 0.9
  tree->Branch("jet_inVetoMap", &v_jet_inVetoMap); // only set for vetoCand jets
  tree->Branch("jet_passRegion", &v_jet_passRegion);
  tree->Branch("jet_passDRmu", &v_jet_passDRmu); // dR(mu) >= 0.4 for both
  tree->Branch("jet_passAnalysis", &v_jet_passAnalysis);

  // ========================================================================
  // CUTFLOW COUNTERS
  // ========================================================================

  Long64_t nTotal = 0, nMETFilter = 0, nTrigger = 0, nTwoMuon = 0, nOS = 0,
           nMassWin = 0, nSelFull = 0;

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

    bool hlt = *HLT_IsoMu24;

    if (kRequireHLT && !hlt)
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
    // MUONS: BS-constrained pT choice, ScaRe, trigger matching, ID SF
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

    // ---- muon flags (stored, NOT cut on)

    bool passEta = true, passID = true, passIso = true;

    for (int i = 0; i < 2; i++) {

      if (fabs(mu[i].eta) > 2.4)
        passEta = false;

      if (!mu[i].mediumId)
        passID = false;

      if (mu[i].iso > 0.25)
        passIso = false;
    }

    // ---- DIMUON + FSR PHOTON RECOVERY

    TLorentzVector mu1, mu2, mu1u, mu2u;

    mu1.SetPtEtaPhiM(mu[0].ptCorr, mu[0].eta, mu[0].phi, 0.105);
    mu2.SetPtEtaPhiM(mu[1].ptCorr, mu[1].eta, mu[1].phi, 0.105);

    mu1u.SetPtEtaPhiM(mu[0].ptUsed, mu[0].eta, mu[0].phi, 0.105);
    mu2u.SetPtEtaPhiM(mu[1].ptUsed, mu[1].eta, mu[1].phi, 0.105);

    TLorentzVector mu1FSR = mu1, mu2FSR = mu2;
    TLorentzVector dimuon = mu1 + mu2;
    TLorentzVector dimuonU = mu1u + mu2u;
    int nFSRUsed = 0;

    // FSR photon requirements: pT > 2 GeV, |eta| < 2.4,
    // relIso03 < 1.8, dROverEt2 < 0.012, and DeltaR(mu,gamma) < 0.5.
    // If several photons are associated with one muon, use the smallest dROverEt2.
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

    if (muHasFSR[0]) {
      mu1FSR += bestPhoton[0];
      ++nFSRUsed;
    }
    if (muHasFSR[1]) {
      mu2FSR += bestPhoton[1];
      ++nFSRUsed;
    }

    TLorentzVector dimuonFSR = mu1FSR + mu2FSR;
    double mll = dimuon.M();
    double mllFSR = dimuonFSR.M();

    // Apply the documented corrected+FSR pT thresholds for the 2024 selection.
    bool passPt = (mu1FSR.Pt() >= 26.) && (mu2FSR.Pt() >= 20.);
    bool passMuonSel = passEta && passID && passIso && passPt;

    // ---- LOOSE MASS REGION

    if (mll < kMassMin || mll > kMassMax)
      continue;

    nMassWin++;

    // ====================================================================
    // JETS
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

      // veto map is only evaluated for veto candidates (as in the old macro)
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

    std::sort(jets.begin(), jets.end(),
              [](const JetInfo &a, const JetInfo &b) {
                return a.p4.Pt() > b.p4.Pt();
              });

    // ---- analysis jets, b-tag counting, b-tag SF

    std::vector<int> ana;

    double btagWeight = 1.0;
    int nMediumB = 0, nLooseB = 0;

    for (size_t k = 0; k < jets.size(); k++) {

      JetInfo &ji = jets[k];

      if (!ji.passAnalysis)
        continue;

      ana.push_back((int)k);

      double absEta = fabs(ji.p4.Eta());

      // b-tagging only defined in the tracker
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
    // DIJET + CATEGORY FLAGS (both pair definitions)
    // ====================================================================

    for (int m = 0; m < 2; m++) {

      int i1 = -1, i2 = -1;

      bool has = selectDijet(jets, ana, m, i1, i2);

      b_jj_has[m] = has;
      b_jj_idx1[m] = i1;
      b_jj_idx2[m] = i2;
      b_jj_mass[m] = -1;
      b_jj_pt[m] = -1;
      b_jj_eta[m] = -99;
      b_jj_deta[m] = -1;
      b_cat_isVBF[m] = false;
      b_cat_isGGH[m] = false;

      if (has) {

        const TLorentzVector &j1 = jets[i1].p4;
        const TLorentzVector &j2 = jets[i2].p4;

        TLorentzVector jj = j1 + j2;

        b_jj_mass[m] = jj.M();
        b_jj_pt[m] = jj.Pt();
        b_jj_eta[m] = jj.Eta();
        b_jj_deta[m] = fabs(j1.Eta() - j2.Eta());

        if (passBVeto && jj.M() > 400. && b_jj_deta[m] > 2.5 &&
            j1.Pt() > 35. && j2.Pt() > 25.)
          b_cat_isVBF[m] = true;
      }

      // ggH = every non-VBF event, including events with 0 or 1 jet
      b_cat_isGGH[m] = passBVeto && !b_cat_isVBF[m];
    }

    // ====================================================================
    // MC WEIGHTS
    // ====================================================================

    double wGen = *genWeight;
    double wPU = 1.0;

    if (nTruePU) {
      b_nTruePU = **nTruePU;
      wPU = getPUWeight(**nTruePU);
    } else {
      b_nTruePU = -1;
    }

    double wMuTot = (double)mu[0].idSF * (double)mu[1].idSF;

    b_w_gen = wGen;
    b_w_norm = normFactor;
    b_w_pu = wPU;
    b_w_muID[0] = mu[0].idSF;
    b_w_muID[1] = mu[1].idSF;
    b_w_muID_tot = wMuTot;
    b_w_btag = btagWeight;
    b_w_total = wGen * normFactor * wPU * wMuTot * btagWeight;

    // ====================================================================
    // FILL TREE
    // ====================================================================

    b_run = *run;
    b_lumi = *luminosityBlock;
    b_event = *eventNumber;

    b_passMETFilters = passMET;
    b_HLT = hlt;
    b_trigMatched = trigMatched;

    b_passMuonEta = passEta;
    b_passMuonID = passID;
    b_passMuonIso = passIso;
    b_passMuonPt = passPt;
    b_passMuonSel = passMuonSel;

    b_passJetVetoMap = !vetoEvent;
    b_passBVeto = passBVeto;

    b_selFull = trigMatched && passMuonSel && !vetoEvent;

    if (b_selFull)
      nSelFull++;

    b_inWindow_115_135 = (mll > 115.0 && mll < 135.0);
    b_massRegion = (mll < 115.0) ? 0 : ((mll <= 135.0) ? 1 : 2);

    b_nJetsStored = (int)jets.size();
    b_nJetsAna = (int)ana.size();
    b_nMediumB = nMediumB;
    b_nLooseB = nLooseB;

    for (int i = 0; i < 2; i++) {

      b_mu_pt_nano[i] = mu[i].ptNano;
      b_mu_pt_bs[i] = mu[i].ptBS;
      b_mu_chi2_bs[i] = mu[i].chi2BS;
      b_mu_usedBS[i] = mu[i].usedBS;
      b_mu_pt_uncorr[i] = mu[i].ptUsed;
      b_mu_pt_corr[i] = mu[i].ptCorr;
      b_mu_eta[i] = mu[i].eta;
      b_mu_phi[i] = mu[i].phi;
      b_mu_charge[i] = mu[i].charge;
      b_mu_mediumId[i] = mu[i].mediumId;
      b_mu_iso[i] = mu[i].iso;
      b_mu_trigDR[i] = mu[i].trigDR;
      b_mu_trigMatched[i] = mu[i].trigMatched;
    }

    b_mm_mass_corr = mll;
    b_mm_mass_uncorr = dimuonU.M();
    b_mm_pt_corr = dimuon.Pt();
    b_mm_pt_uncorr = dimuonU.Pt();
    b_mm_eta_corr = dimuon.Eta();
    b_mm_phi_corr = dimuon.Phi();
    b_mm_mass_fsr = mllFSR;
    b_mm_pt_fsr = dimuonFSR.Pt();
    b_nFSR = nFSRUsed;
    b_mu_fsr_pt[0] = mu1FSR.Pt();
    b_mu_fsr_pt[1] = mu2FSR.Pt();
    b_mu_hasFSR[0] = muHasFSR[0];
    b_mu_hasFSR[1] = muHasFSR[1];

    v_jet_pt.clear();
    v_jet_pt_nano.clear();
    v_jet_pt_raw.clear();
    v_jet_eta.clear();
    v_jet_phi.clear();
    v_jet_mass.clear();
    v_jet_btag.clear();
    v_jet_btagSF.clear();
    v_jet_hadFlav.clear();
    v_jet_chEmEF.clear();
    v_jet_neEmEF.clear();
    v_jet_dRmu.clear();
    v_jet_idOK.clear();
    v_jet_vetoCand.clear();
    v_jet_inVetoMap.clear();
    v_jet_passRegion.clear();
    v_jet_passDRmu.clear();
    v_jet_passAnalysis.clear();

    for (auto const &ji : jets) {

      v_jet_pt.push_back(ji.p4.Pt());
      v_jet_pt_nano.push_back(ji.ptNano);
      v_jet_pt_raw.push_back(ji.ptRaw);
      v_jet_eta.push_back(ji.p4.Eta());
      v_jet_phi.push_back(ji.p4.Phi());
      v_jet_mass.push_back(ji.p4.M());
      v_jet_btag.push_back(ji.btag);
      v_jet_btagSF.push_back(ji.btagSF);
      v_jet_hadFlav.push_back(ji.hadFlav);
      v_jet_chEmEF.push_back(ji.chEmEF);
      v_jet_neEmEF.push_back(ji.neEmEF);
      v_jet_dRmu.push_back(ji.dRmu);
      v_jet_idOK.push_back(ji.idOK);
      v_jet_vetoCand.push_back(ji.vetoCand);
      v_jet_inVetoMap.push_back(ji.inVetoMap);
      v_jet_passRegion.push_back(ji.passRegion);
      v_jet_passDRmu.push_back(ji.passDRmu);
      v_jet_passAnalysis.push_back(ji.passAnalysis);
    }

    tree->Fill();

  } // while(reader.Next())

  // ========================================================================
  // CUTFLOW
  // ========================================================================

  std::cout << "\nCutflow (loose baseline -> tree):"
            << "\n  total                  = " << nTotal
            << "\n  MET filters            = " << nMETFilter
            << "\n  HLT_IsoMu24            = " << nTrigger
            << "\n  nMuon == 2             = " << nTwoMuon
            << "\n  opposite sign          = " << nOS
            << "\n  " << kMassMin << " < mmumu < " << kMassMax
            << "   = " << nMassWin << "  (= tree entries)"
            << "\n  of which sel_full      = " << nSelFull << std::endl;

  // ========================================================================
  // WRITE
  // ========================================================================

  out.cd();

  tree->Write();

  out.Close();

  if (Muon_bsConstrainedPt)
    delete Muon_bsConstrainedPt;

  if (Muon_bsConstrainedChi2)
    delete Muon_bsConstrainedChi2;

  if (TrigObj_filterBits)
    delete TrigObj_filterBits;

  if (nTruePU)
    delete nTruePU;

  std::cout << "\n====================================\n"
            << " MC TREE-MAKER FINISHED\n"
            << " Output: " << outFile << "\n"
            << "====================================\n";
}
