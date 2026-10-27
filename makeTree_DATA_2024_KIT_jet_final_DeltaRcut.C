// ============================================================================
// Nilay Bostan - October 2026
// TREE-MAKER MACRO -- DATA ONLY (AGREED SYNC SELECTION), CMS Run2024 H->mumu
//
// KEEP IN SYNC with makeTree_MC_2025.C: the [SHARED] block is identical in both,
// and the output tree has the same branch structure (DATA has no weights / MC truth,
// and no dataset_id: that branch exists only in the MC tree).
//
// Baseline written to the tree:
//   golden JSON, MET filters, HLT_IsoMu24, nMuon >= 2, opposite sign of the two
//   leading (raw pT) muons, raw pT > 20 for both, 50 < m(mumu, RAW) < 200,
//   and for DATA: events with RAW m(mumu) in [115,135] are NOT written (blinding).
// Everything else (trigger matching, muon eta/ID/iso/pT after ScaRe+FSR, third
// muon veto, electron veto, jet veto map, b-veto, VBF/ggH) is a FLAG:
//     sel_full && cat_isVBF[mode]  /  sel_full && cat_isGGH[mode]
// mode 0 = two leading analysis jets, mode 1 = pair with max m_jj.
// Analysis jets: pT > 30 (after JEC), |eta| < 4.7, tight jetID (jetid.json),
// region rules, dR(jet, mu) >= 0.4.
//
// "RAW" = NanoAOD Muon_pt (no BS constraint, no ScaRe) of Muon[0] and Muon[1].
// mm_mass_raw = m(mumu) from these RAW pTs (no BSC, no ScaRe, no FSR).
// mm_mass_uncorr is the mass with the BS-constrained pT (before ScaRe) as before.
//
// MUON / DIMUON VARIABLES (same as MC)
//   Official m(mumu) = mm_mass_fsr : pT = BSC pT if chi2_BSC < 30 else NanoAOD pT,
//   then ScaRe, then FSR. mm_mass_corr is ScaRe ONLY (no FSR) and is NOT m(mumu).
//   massRegion / inWindow_115_135 are computed from mm_mass_fsr.
//   (The blinding itself is applied on RAW m(mumu), as requested.)
//   Single-muon pT, i = 1 (leading) / 2 (subleading after ScaRe+FSR re-sort):
//     mu{i}_pt, mu{i}_pt_nano   NanoAOD Muon_pt (no BSC, no ScaRe, no FSR)
//     mu{i}_pt_nano_corr        NanoAOD pT + ScaRe
//     mu{i}_pt_nano_corr_FSR    NanoAOD pT + ScaRe + FSR
//     mu{i}_pt_raw_noCorr       BSC pT (chi2<30) else Nano pT, no ScaRe, no FSR
//     mu{i}_pt_raw_scale        same, only Scale applied (no Res, no FSR)
//     mu{i}_pt_raw_corr         same, ScaRe (no FSR)
//     mu{i}_pt_FSR_noCorr       same, FSR (no ScaRe)
//     mu{i}_pt_FSR_corr         same, ScaRe + FSR
//   Jets: jet{1,2}_{pt,pt_nano,pt_raw,eta,phi,mass,btag} = two leading analysis
//   jets; jj_{dphi,dR,y,zstar}[mode]; cat_id[mode] = -1 (b-veto fail), 0 ggH, 1 VBF.
//
// CHANGES in this version (merge of the JEC-fixed and the runext macros):
//   [A] DATA jet JEC applied explicitly from jet_jerc.json:
//       raw pT -> L1FastJet -> L2Relative -> L3Absolute -> L2L3Residual(run), all DATA
//       (l1FastJetData, l2RelativeData, l3AbsoluteData). If the chain throws, the
//       NanoAOD jet pT is used and gResidualFail is incremented.
//   [B] Jet mass is scaled with the same factor as the pT.
//   [C] JEC self-test printed after loading (expect ~74 for raw 40, eta 2.75).
//   [D] Muon NanoAOD pT upper cut in passMuonPt (kMuPtMaxNano = 1000 GeV).
//   [E] dR(jet, mu) >= 0.4 (passDRmu) is now REQUIRED in passAnalysis.
//   [F] Jet ID comes from jetid.json (JetIdJson.h): jet_idTight, jet_idTightLV,
//       jet_jetId_nano are stored per jet.
//   [G] dataset_id branch removed (MC only).
//   Kept from runext: duplicate-file removal, VBF debug counters, jetIdJsonFile arg.
//   NOT changed: blinding is on RAW m(mumu) only (as in runext).
//   The [SHARED] block changed ([A], [E] and the data JEC globals): copy it to the MC macro.
// ============================================================================

#pragma cling add_include_path("/cvmfs/cms.cern.ch/el9_amd64_gcc12/external/py3-correctionlib/2.2.2-120738cfaaf3f7c1056fe67d97e25dac/lib/python3.9/site-packages/correctionlib/include")
#pragma cling add_library_path("/cvmfs/cms.cern.ch/el9_amd64_gcc12/external/py3-correctionlib/2.2.2-120738cfaaf3f7c1056fe67d97e25dac/lib/python3.9/site-packages/correctionlib/lib")
#pragma cling load("correctionlib")

#include <nlohmann/json.hpp>
// ============================================================================
// [SHARED] BLOCK -- identical in makeTree_DATA_2025.C and makeTree_MC_2025.C
// Nilay Bostan, October 2026
//
// Everything from here down to the end of the [SHARED] block must be kept
// IDENTICAL in both macros (edit one, copy to the other).
// ============================================================================

#include <TNamed.h>
#include "JetCorrections.h"
#include "correction.h"
// ----------------------------------------------------------------------------
// Official MuonPOG ScaRe code (muonscarekit MuonScaRe.cc), inlined.
// ----------------------------------------------------------------------------
// ============================================================================
// MuonScaReKit.h
// Adapted from the official CMS MuonPOG MuonScaRe.cc (muonscarekit).
// The physics/maths is UNCHANGED; the only differences are:
//   * everything lives in namespace ScaRe,
//   * the CorrectionSet is passed as an argument instead of a global `cset`
//     (needed because the analysis uses TWO sets: BS-constrained and regular),
//   * std:: prefixes (no `using namespace std` needed).
//
// Official recipe:
//   DATA : pt_corr = ptScale(isData = true , ...)
//   MC   : pt_corr = ptScale(isData = false, ...)  then  ptResol(...)
// ============================================================================

#include <boost/math/special_functions/erf.hpp>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>


namespace ScaRe {

// Double-sided Crystal Ball (pdf / cdf / inverse cdf), used to draw the random
// number that smears MC. Identical to the official implementation (including
// pi = 3.14159, kept on purpose so numbers match the official kit).
struct CrystalBall {
  double pi = 3.14159;
  double sqrtPiOver2 = std::sqrt(pi / 2.0);
  double sqrt2 = std::sqrt(2.0);
  double m, s, a, n;
  double B, C, D, N, NA, Ns, NC, F, G, k;
  double cdfMa, cdfPa;

  CrystalBall() : m(0), s(1), a(10), n(10) { init(); }
  CrystalBall(double mean, double sigma, double alpha, double nn)
      : m(mean), s(sigma), a(alpha), n(nn) { init(); }

  void init() {
    double fa = std::fabs(a);
    double ex = std::exp(-fa * fa / 2);
    double A = std::pow(n / fa, n) * ex;
    double C1 = n / fa / (n - 1) * ex;
    double D1 = 2 * sqrtPiOver2 * std::erf(fa / sqrt2);
    B = n / fa - fa;
    C = (D1 + 2 * C1) / C1;
    D = (D1 + 2 * C1) / 2;
    N = 1.0 / s / (D1 + 2 * C1);
    k = 1.0 / (n - 1);
    NA = N * A;
    Ns = N * s;
    NC = Ns * C1;
    F = 1 - fa * fa / n;
    G = s * n / fa;
    cdfMa = cdf(m - a * s);
    cdfPa = cdf(m + a * s);
  }

  double pdf(double x) const {
    double d = (x - m) / s;
    if (d < -a) return NA * std::pow(B - d, -n);
    if (d > a) return NA * std::pow(B + d, -n);
    return N * std::exp(-d * d / 2);
  }
  double pdf(double x, double ks, double dm) const {
    double d = (x - m - dm) / (s * ks);
    if (d < -a) return NA / ks * std::pow(B - d, -n);
    if (d > a) return NA / ks * std::pow(B + d, -n);
    return N / ks * std::exp(-d * d / 2);
  }
  double cdf(double x) const {
    double d = (x - m) / s;
    if (d < -a) return NC / std::pow(F - s * d / G, n - 1);
    if (d > a) return NC * (C - std::pow(F + s * d / G, 1 - n));
    return Ns * (D - sqrtPiOver2 * std::erf(-d / sqrt2));
  }
  double invcdf(double u) const {
    if (u < cdfMa) return m + G * (F - std::pow(NC / u, k));
    if (u > cdfPa) return m - G * (F - std::pow(C - u / NC, -k));
    return m - sqrt2 * s * boost::math::erf_inv((D - u / Ns) / sqrtPiOver2);
  }
};

// ---------------------------------------------------------------------------
// MC resolution smearing helpers
// ---------------------------------------------------------------------------

// Random number following the Crystal Ball of this (|eta|, nTrackerLayers) bin.
// Deterministic: seeded by (event, lumi, phi) through the "RandomSmearing" entry.
inline double getRndm(const correction::CorrectionSet &cset, double eta, double phi,
                      float nL, int evtNumber, int lumiNumber) {
  const double aeta = std::fabs(eta);
  const double dnL = (double)nL;
  auto cb_params = cset.at("cb_params");
  double mean = cb_params->evaluate({aeta, dnL, 0});
  double sigma = cb_params->evaluate({aeta, dnL, 1});
  double n = cb_params->evaluate({aeta, dnL, 2});
  double alpha = cb_params->evaluate({aeta, dnL, 3});

  CrystalBall cb(mean, sigma, alpha, n);
  double rndm = cset.at("RandomSmearing")->evaluate({(int)evtNumber, (int)lumiNumber, phi});
  return cb.invcdf(rndm);
}

// pT-dependent width of the smearing: polynomial in pT, floored at 0.
inline double getStd(const correction::CorrectionSet &cset, double pt, double eta, float nL) {
  const double aeta = std::fabs(eta);
  const double dnL = (double)nL;
  auto poly = cset.at("poly_params");
  double p0 = poly->evaluate({aeta, dnL, 0});
  double p1 = poly->evaluate({aeta, dnL, 1});
  double p2 = poly->evaluate({aeta, dnL, 2});
  double sigma = p0 + p1 * pt + p2 * pt * pt;
  if (sigma < 0) sigma = 0;
  return sigma;
}

// Residual smearing factor k = sqrt(k_data^2 - k_mc^2); 0 if MC is already broader.
inline double getK(const correction::CorrectionSet &cset, double eta, const std::string &var) {
  const double aeta = std::fabs(eta);
  double k_data = cset.at("k_data")->evaluate({aeta, var});
  double k_mc = cset.at("k_mc")->evaluate({aeta, var});
  double k = 0;
  if (k_mc < k_data) k = std::sqrt(k_data * k_data - k_mc * k_mc);
  return k;
}

// ---------------------------------------------------------------------------
// MC ONLY: pT resolution smearing.  Input pt = pT AFTER the scale correction.
// Returns the input unchanged if pt < low_pt_threshold, pt > 200, NaN, or the
// smearing factor is outside (0.1, 2).
// ---------------------------------------------------------------------------
inline double ptResol(const correction::CorrectionSet &cset, double pt, double eta, double phi,
                      float nL, int evtNumber, int lumiNumber, double low_pt_threshold = 26) {
  double rndm = getRndm(cset, eta, phi, nL, evtNumber, lumiNumber);
  double sd = getStd(cset, pt, eta, nL);
  double k = getK(cset, eta, "nom");

  double ptc = pt * (1 + k * sd * rndm);
  if (std::isnan(ptc)) ptc = pt;
  if (ptc / pt > 2 || ptc / pt < 0.1 || ptc < 0 || pt < low_pt_threshold || pt > 200) ptc = pt;
  return ptc;
}

// ---------------------------------------------------------------------------
// DATA and MC: pT scale correction  pt' = 1 / (m/pt + charge * a)
//   isData = true  -> a_data, m_data
//   isData = false -> a_mc,   m_mc
// Returns pt unchanged if pt < low_pt_threshold.
// ---------------------------------------------------------------------------
inline double ptScale(const correction::CorrectionSet &cset, bool isData, double pt, double eta,
                      double phi, int charge, double low_pt_threshold = 26) {
  if (pt < low_pt_threshold) return pt;
  const std::string dtmc = isData ? "data" : "mc";
  const std::string nom = "nom";
  double a = cset.at("a_" + dtmc)->evaluate({eta, phi, nom});
  double m = cset.at("m_" + dtmc)->evaluate({eta, phi, nom});
  return 1. / (m / pt + charge * a);
}

}  // namespace ScaRe

#include <TChain.h>
#include <TError.h>
#include <TFile.h>
#include <TLeaf.h>
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
#include <utility>
#include <vector>

// ============================================================================
// GLOBAL CORRECTIONS
// ============================================================================
JetCorrections jetCorr;
std::shared_ptr<const correction::CorrectionSet> btagCorr;

std::shared_ptr<const correction::CorrectionSet> jetJercCS;
std::shared_ptr<const correction::Correction> jerPtResolution;
std::shared_ptr<const correction::Correction> jerScaleFactor;
std::shared_ptr<const correction::Correction> l2l3ResidualData;
std::shared_ptr<const correction::Correction> l1FastJetMC;
std::shared_ptr<const correction::Correction> l2RelativeMC;
std::shared_ptr<const correction::Correction> l3AbsoluteMC;
// [A] DATA JEC chain (loaded from jet_jerc.json)
std::shared_ptr<const correction::Correction> l1FastJetData;
std::shared_ptr<const correction::Correction> l2RelativeData;
std::shared_ptr<const correction::Correction> l3AbsoluteData;

// ============================================================================
// JET ID from jetid.json.gz (AK4PUPPI)                         [Oct 25, 2026]
//   idTight   = AK4PUPPI_Tight           -> idOK (analysis jets, "tight jetID")
//   idTightLV = AK4PUPPI_TightLeptonVeto -> veto-map candidates
//   Jet_jetId (NanoAOD) is stored as a cross-check: bit1 = tight, bit2 = tightLepVeto
//   (NanoAODv15 has no Jet_jetId -> jet_jetId_nano = -1)
// The ID inputs (fractions, multiplicities) do not depend on JEC; eta is not changed by JEC.
// ============================================================================
#include "JetIdJson.h"  // plain C++ jet ID generated from the jetid json (gen_jetid_cpp.py)
std::shared_ptr<const correction::CorrectionSet> jetIdCS;
std::shared_ptr<const correction::Correction> jetIdTight, jetIdTightLV;

bool loadJetId(const std::string &file) {
  try {
    jetIdCS = correction::CorrectionSet::from_file(file);  // .json.gz is fine
    jetIdTight = jetIdCS->at("AK4PUPPI_Tight");
    jetIdTightLV = jetIdCS->at("AK4PUPPI_TightLeptonVeto");
  } catch (const std::exception &e) {
    std::cerr << "ERROR: cannot load jetid json " << file << ": " << e.what() << std::endl;
    return false;
  }
  std::cout << "Jet ID (AK4PUPPI) loaded from: " << file << std::endl;
  return true;
}

// bit0 = Tight, bit1 = TightLeptonVeto
int evalJetId(double eta, double chHEF, double neHEF, double chEmEF, double neEmEF,
              double muEF, int chMult, int neMult, int mult) {
  const double ae = fabs(eta);
  int bits = 0;
  // cuts come from JetIdJson.h (generated from the json); correctionlib is not used here
  if (jetIdPass_AK4PUPPI_Tight(ae, chHEF, neHEF, chEmEF, neEmEF, muEF, chMult, neMult, mult))
    bits |= 1;
  if (jetIdPass_AK4PUPPI_TightLeptonVeto(ae, chHEF, neHEF, chEmEF, neEmEF, muEF, chMult, neMult, mult))
    bits |= 2;
  return bits;
}

// integer array branch that is UChar_t or Int_t depending on the NanoAOD version
struct IntArrayReader {
  TTreeReaderArray<unsigned char> *u = nullptr;
  TTreeReaderArray<int> *i = nullptr;
  bool ok = false;
  IntArrayReader(TTreeReader &reader, TChain &chain, const char *name) {
    TLeaf *l = chain.GetLeaf(name);
    const std::string t = l ? l->GetTypeName() : "";
    if (t == "UChar_t") u = new TTreeReaderArray<unsigned char>(reader, name);
    else if (t == "Int_t") i = new TTreeReaderArray<int>(reader, name);
    if (t.empty()) std::cout << name << " : not in this NanoAOD" << std::endl;
    else           std::cout << name << " type = " << t << std::endl;
    ok = (u != nullptr || i != nullptr);
  }
  ~IntArrayReader() { delete u; delete i; }
  int operator[](size_t k) const { return u ? (int)(*u)[k] : (i ? (*i)[k] : -1); }
};

std::shared_ptr<const correction::Correction>
findCorrectionBySuffix(const std::shared_ptr<const correction::CorrectionSet> &cs,
                       const std::string &suffix) {
  for (auto const &kv : *cs) {
    const std::string &name = kv.first;
    if (name.size() >= suffix.size() &&
        name.compare(name.size() - suffix.size(), suffix.size(), suffix) == 0)
      return cs->at(name);
  }
  std::cerr << "ERROR: no correction ending in " << suffix
            << " found. Available corrections:" << std::endl;
  for (auto const &kv : *cs)
    std::cerr << "   - " << kv.first << std::endl;
  return nullptr;
}

bool loadJetJercRegionCorrections(const std::string &jetJercJsonFile) {
  try {
    jetJercCS = correction::CorrectionSet::from_file(jetJercJsonFile);
    jerPtResolution = findCorrectionBySuffix(jetJercCS, "_JRV2_MC_PtResolution_AK4PFPuppi");
    jerScaleFactor = findCorrectionBySuffix(jetJercCS, "_JRV2_MC_ScaleFactor_AK4PFPuppi");
    l2l3ResidualData = findCorrectionBySuffix(jetJercCS, "_DATA_L2L3Residual_AK4PFPuppi");
    l1FastJetMC = findCorrectionBySuffix(jetJercCS, "_MC_L1FastJet_AK4PFPuppi");
    l2RelativeMC = findCorrectionBySuffix(jetJercCS, "_MC_L2Relative_AK4PFPuppi");
    l3AbsoluteMC = findCorrectionBySuffix(jetJercCS, "_MC_L3Absolute_AK4PFPuppi");
    // [A]
    l1FastJetData = findCorrectionBySuffix(jetJercCS, "_DATA_L1FastJet_AK4PFPuppi");
    l2RelativeData = findCorrectionBySuffix(jetJercCS, "_DATA_L2Relative_AK4PFPuppi");
    l3AbsoluteData = findCorrectionBySuffix(jetJercCS, "_DATA_L3Absolute_AK4PFPuppi");
    if (!jerPtResolution || !jerScaleFactor || !l2l3ResidualData ||
        !l1FastJetMC || !l2RelativeMC || !l3AbsoluteMC ||
        !l1FastJetData || !l2RelativeData || !l3AbsoluteData)
      return false;
  } catch (const std::exception &e) {
    std::cerr << "ERROR: cannot load jet_jerc JSON " << jetJercJsonFile << ": "
              << e.what() << std::endl;
    return false;
  }
  std::cout << "Year/region-aware jet corrections loaded from: " << jetJercJsonFile
            << std::endl;
  return true;
}

// ============================================================================
// ANALYSIS CONSTANTS  (agreed sync selection)
// ============================================================================
static const double kMuMass = 0.105658;

// trigger matching on RAW Muon_pt (before ScaRe), IsoMu24 only
static const double kTrigMuRawPtMin = 26.;
static const double kTrigObjPtMin = 24.;
static const double kTrigMatchDR = 0.1;
static const ULong64_t kIsoMu24FilterMask = 0; // 0 = not applied (bit needs validation)

// loose baseline
static const bool kRequireMETFilters = true;
static const bool kRequireHLT = true;
static const bool kRequireOS = true;

// muons: raw pT > 20, 50 < m(mumu, RAW) < 200, raw = NanoAOD Muon_pt (before BS/ScaRe)
static const double kRawMuPtMin = 20.;
static const double kMassMin = 50.;
static const double kMassMax = 200.;
static const double kMu1PtMin = 26.; // after ScaRe + FSR, after re-sorting
static const double kMu2PtMin = 20.;
static const double kMuEtaMax = 2.4;
static const double kMuIsoMax = 0.25;
// [D] sanity upper cut on NanoAOD muon pT (broken tracks reach 1e4-1e5 GeV)
static const double kMuPtMaxNano = 1000.;

// DATA blinding on RAW m(mumu)
static const double kBlindLo = 115.;
static const double kBlindHi = 135.;

// jets
static const double kJetPtStore = 15.;   // stored in the tree (veto-map threshold)
static const double kJetPtMin = 30.;     // analysis jet pT after JEC
static const double kJetDRMuMin = 0.4;   // [E] dR(jet, mu) >= 0.4
// Analysis jets: |eta| < 4.7 (agreed selection; needed for VBF: dEta > 2.5 and
// mjj > 400 require forward jets). There is no |eta| < 2.4 jet selection.
static const double kJetEtaMax = 4.7;
static const double kBtagEtaMax = 2.5;

// electron veto (slide: |eta| < 2.5, pT > 20, MVA WP90). Gap exclusion was in
// the old "Purdue" version but is NOT on the slide -> off by default.
static const bool kEleExcludeGap = false;

int gYear = 2024;
double gBtagWP_L = -1.;
double gBtagWP_M = -1.;

// ============================================================================
// STRUCTS
// ============================================================================
struct MuInfo {
  int nanoIdx = -1;
  int nL = 0; // Muon_nTrackerLayers (needed by the ScaRe resolution smearing)
  float ptNano = 0, ptBS = -1, chi2BS = -1;
  float ptUsed = 0; // before ScaRe (BS-constrained if chi2 < 30, else NanoAOD pT)
  float ptCorr = 0; // after ScaRe (on ptUsed), no FSR
  // --- single-muon pT variants (filled in the main loop, before the final re-sort)
  float ptNanoCorr = 0;    // NanoAOD pT + ScaRe (no BSC, no FSR)
  float ptNanoCorrFSR = 0; // NanoAOD pT + ScaRe + FSR (no BSC)
  float ptRawScale = 0;    // ptUsed + ScaRe SCALE correction only (ptScale), no smearing, no FSR
  float ptFSRnoCorr = 0;   // ptUsed + FSR (no ScaRe)
  float ptFSRcorr = 0;     // ptUsed + ScaRe + FSR   (the pT used for m_mumu)
  float eta = 0, phi = 0, iso = 0;
  float trigDR = -1;
  float idSF = 1; // MC only
  int charge = 0;
  bool usedBS = false, mediumId = false, trigMatched = false;
};

struct JetInfo {
  TLorentzVector p4; // JEC-corrected
  float ptNano = 0, ptRaw = 0, btag = 0, chEmEF = 0, neEmEF = 0, dRmu = 0;
  float btagSF = 1; // MC only
  int hadFlav = 0;  // MC only
  bool idOK = false, vetoCand = false, inVetoMap = false;
  bool idTight = false, idTightLV = false;
  int jetIdNano = 0; // NanoAOD Jet_jetId (cross-check)
  bool passRegion = false, passDRmu = false, passAnalysis = false;
};

struct DijetInfo {
  bool has = false, vbf = false;
  int i1 = -1, i2 = -1;
  float mass = -1, pt = -1, eta = -99, deta = -1, dphi = -1, dR = -1, y = -99;
};

// ============================================================================
// MUON ScaRe  (official MuonPOG recipe, see MuonScaReKit.h)
//   DATA : pT = ptScale(data)
//   MC   : pT = ptScale(mc), then ptResol (Crystal-Ball smearing, needs nTrackerLayers)
// Two correction sets: BS-constrained pT (chi2 < 30) and regular (NanoAOD) pT.
// Muons with pT < 26 GeV (and, for the smearing, pT > 200 GeV) are left unchanged
// by the official code.
// ============================================================================
std::shared_ptr<const correction::CorrectionSet> muonCorrBSC, muonCorrRegular;

// scale correction only (DATA: this IS the full correction; MC: before smearing)
double scaReScaleOnly(bool isData, double pt, double eta, double phi, int charge, bool isBSC) {
  const correction::CorrectionSet &cs = isBSC ? *muonCorrBSC : *muonCorrRegular;
  return ScaRe::ptScale(cs, isData, pt, eta, phi, charge);
}

// full correction: scale, plus smearing in MC
double scaReFull(bool isData, double pt, double eta, double phi, int charge, int nL,
                 ULong64_t evt, UInt_t lumi, bool isBSC) {
  const correction::CorrectionSet &cs = isBSC ? *muonCorrBSC : *muonCorrRegular;
  double ptS = ScaRe::ptScale(cs, isData, pt, eta, phi, charge);
  if (isData) return ptS;
  return ScaRe::ptResol(cs, ptS, eta, phi, (float)nL, (int)evt, (int)lumi);
}

// ============================================================================
// DELTA PHI / DELTA R   (must come before the templates below)
// ============================================================================
double deltaPhi(double phi1, double phi2) {
  double dphi = phi1 - phi2;
  while (dphi > M_PI) dphi -= 2 * M_PI;
  while (dphi <= -M_PI) dphi += 2 * M_PI;
  return dphi;
}

double deltaR(double eta1, double phi1, double eta2, double phi2) {
  double deta = eta1 - eta2;
  double dphi = deltaPhi(phi1, phi2);
  return sqrt(deta * deta + dphi * dphi);
}

// pT of (muon with given pT, eta, phi) + FSR photon; unchanged if no photon.
double ptWithFSR(double pt, double eta, double phi, const TLorentzVector &gamma, bool has) {
  if (!has) return pt;
  TLorentzVector m;
  m.SetPtEtaPhiM(pt, eta, phi, kMuMass);
  return (m + gamma).Pt();
}

// ============================================================================
// FORWARD-REGION JET RULES (2024: HF fixed, HE needs pT > 50)
// ============================================================================
bool passRegionJetCut(double corrPt, double eta) {
  double ae = fabs(eta);
  if (gYear <= 2023 && ae >= 3.0 && ae < 5.0 && corrPt < 50.) return false;
  if (gYear <= 2024 && ae >= 2.5 && ae < 3.0 && corrPt < 50.) return false;
  return true;
}

// ============================================================================
// JET pT CORRECTION (region aware)                                   [A]
//   DATA: raw pT -> L1FastJet -> L2Relative -> L3Absolute -> L2L3Residual(run)
//         (applied explicitly from jet_jerc.json; JetCorrections is NOT used)
//         In 2.0<|eta|<2.5 with corrected pT < 30 the residual is evaluated at
//         30 GeV. If the chain throws: NanoAOD jet pT, gResidualFail++.
//   MC  : JERC (JetCorrections) + gen-matched JER in 2.5<|eta|<3.0
// ============================================================================
static long gResidualFail = 0;

double getCorrectedJetPtRegionAware(double pt, double eta, double phi,
                                    double rawFactor, double area, double rho,
                                    unsigned int run, bool isMC,
                                    bool hasGenMatch, double genPt) {
  double corrPt;
  const double ae = fabs(eta);

  if (!isMC) {
    try {
      const double rawPt = pt * (1.0 - rawFactor);
      const double p1 = rawPt * l1FastJetData->evaluate({area, eta, rawPt, rho});
      const double p2 = p1 * l2RelativeData->evaluate({eta, phi, p1});
      const double p3 = p2 * l3AbsoluteData->evaluate({eta, p2});
      double ptForRes = p3;
      if (gYear == 2024 && ae >= 2.0 && ae < 2.5 && p3 < 30.0) ptForRes = 30.0;
      corrPt = p3 * l2l3ResidualData->evaluate({(double)run, eta, ptForRes});
    } catch (const std::exception &e) {
      if (gResidualFail < 3)
        std::cerr << "WARNING: data JEC chain failed (run " << run << ", eta " << eta
                  << "): " << e.what() << " -> using NanoAOD jet pT" << std::endl;
      gResidualFail++;
      corrPt = pt;
    }
  } else {
    corrPt = jetCorr.getCorrectedJetPt(pt, eta, phi, rawFactor, isMC);
    if (gYear <= 2024 && ae >= 2.5 && ae < 3.0 && hasGenMatch && genPt > 0.) {
      double reso = jerPtResolution->evaluate({eta, corrPt, rho});
      double sf = jerScaleFactor->evaluate({eta, corrPt});
      if (fabs(corrPt - genPt) < 3.0 * reso * corrPt)
        corrPt = std::max(0.0, genPt + sf * (corrPt - genPt));
    }
  }
  return corrPt;
}

// ============================================================================
// B-TAG WORKING POINTS
// ============================================================================
bool loadBtagWorkingPoints(const std::shared_ptr<const correction::CorrectionSet> &cs) {
  try {
    auto c = cs->at("UParTAK4_wp_values");
    gBtagWP_L = c->evaluate({std::string("L")});
    gBtagWP_M = c->evaluate({std::string("M")});
  } catch (const std::exception &e) {
    std::cerr << "ERROR: cannot read UParTAK4 working points: " << e.what()
              << "\nAvailable corrections in btagging.json:" << std::endl;
    for (auto const &corr : *cs) std::cerr << "   - " << corr.first << std::endl;
    return false;
  }
  std::cout << "UParTAK4 WPs from correctionlib: L = " << gBtagWP_L
            << ", M = " << gBtagWP_M << std::endl;
  return true;
}

// ============================================================================
// DIJET SELECTION  (mode 0: two leading analysis jets, mode 1: max m_jj pair)
// `jets` is sorted by pT, `ana` = ascending indices into `jets`.
// ============================================================================
bool selectDijet(const std::vector<JetInfo> &jets, const std::vector<int> &ana,
                 int mode, int &i1, int &i2) {
  i1 = -1; i2 = -1;
  if (ana.size() < 2) return false;
  if (mode == 0) { i1 = ana[0]; i2 = ana[1]; return true; }
  double best = -1.;
  for (size_t a = 0; a < ana.size(); a++)
    for (size_t b = a + 1; b < ana.size(); b++) {
      double m = (jets[ana[a]].p4 + jets[ana[b]].p4).M();
      if (m > best) { best = m; i1 = ana[a]; i2 = ana[b]; }
    }
  return i1 >= 0;
}

// VBF: mjj > 400, |dEta| > 2.5, pT(lead) > 35, pT(sublead) > 30, b-veto passed
DijetInfo evalDijet(const std::vector<JetInfo> &jets, const std::vector<int> &ana,
                    int mode, bool passBVeto) {
  DijetInfo d;
  int i1, i2;
  if (!selectDijet(jets, ana, mode, i1, i2)) return d;
  const TLorentzVector &j1 = jets[i1].p4;
  const TLorentzVector &j2 = jets[i2].p4;
  TLorentzVector jj = j1 + j2;
  d.has = true; d.i1 = i1; d.i2 = i2;
  d.mass = jj.M(); d.pt = jj.Pt(); d.eta = jj.Eta();
  d.deta = fabs(j1.Eta() - j2.Eta());
  d.dphi = fabs(deltaPhi(j1.Phi(), j2.Phi()));
  d.dR = j1.DeltaR(j2);
  d.y = jj.Rapidity();
  d.vbf = passBVeto && d.mass > 400. && d.deta > 2.5 && j1.Pt() > 35. && j2.Pt() > 30.;
  return d;
}

// ============================================================================
// TRIGGER-OBJECT MATCHING (returns min dR, -1 if no candidate object)
// ============================================================================
template <class AF, class AU, class AB>
double minTrigObjDR(int nObj, AF &pt, AF &eta, AF &phi, AU &id, AB *filterBits,
                    double mEta, double mPhi) {
  double minDR = -1.;
  for (int t = 0; t < nObj; t++) {
    if (id[t] != 13) continue;
    if (pt[t] < kTrigObjPtMin) continue;
    if (kIsoMu24FilterMask != 0 &&
        !(filterBits && (((*filterBits)[t]) & kIsoMu24FilterMask)))
      continue;
    double dr = deltaR(mEta, mPhi, eta[t], phi[t]);
    if (minDR < 0. || dr < minDR) minDR = dr;
  }
  return minDR;
}

// ============================================================================
// FSR PHOTONS. `mu` must be the UNSORTED muons (mu[i].nanoIdx == i).
// pT > 2, |eta| < 2.4, relIso03 < 1.8, dROverEt2 < 0.012, dR(mu,gamma) < 0.5,
// best photon per muon = smallest dROverEt2.
// ============================================================================
template <class AF, class AI>
void findFSRPhotons(AF &pt, AF &eta, AF &phi, AF &drOverEt2, AF &relIso,
                    AI &muonIdx, const MuInfo *mu, TLorentzVector *best, bool *has) {
  double bestD[2] = {1.e9, 1.e9};
  has[0] = false; has[1] = false;
  for (size_t ig = 0; ig < (size_t)pt.GetSize(); ++ig) {
    const double gpt = pt[ig], geta = eta[ig], gphi = phi[ig];
    const double d2 = drOverEt2[ig], iso = relIso[ig];
    const int im = muonIdx[ig];
    if (gpt <= 2.0 || fabs(geta) >= 2.4 || iso >= 1.8 || d2 >= 0.012) continue;
    if (im < 0 || im > 1) continue;
    if (deltaR(mu[im].eta, mu[im].phi, geta, gphi) >= 0.5) continue;
    if (d2 < bestD[im]) {
      bestD[im] = d2;
      best[im].SetPtEtaPhiM(gpt, geta, gphi, 0.0);
      has[im] = true;
    }
  }
}

// ============================================================================
// ELECTRON VETO COUNT (pT > 20, |eta| < 2.5, MVA-iso WP90)
// ============================================================================
template <class AF, class AB>
int countVetoElectrons(AF &pt, AF &eta, AB &wp90) {
  int n = 0;
  for (size_t ie = 0; ie < (size_t)pt.GetSize(); ie++) {
    double ae = fabs(eta[ie]);
    if (pt[ie] <= 20. || ae >= 2.5) continue;
    if (kEleExcludeGap && ae > 1.44 && ae < 1.57) continue;
    if (wp90[ie]) n++;
  }
  return n;
}

// ============================================================================
// FILTER OUT UNREACHABLE / BROKEN INPUT FILES
// ============================================================================
std::vector<std::string> FilterReadableFiles(const std::vector<std::string> &inputFiles) {
  std::vector<std::string> goodFiles, badFiles;
  for (auto const &file : inputFiles) {
    TFile *tf = nullptr;
    for (int attempt = 1; attempt <= 2 && (!tf || tf->IsZombie()); attempt++) {
      if (tf) { delete tf; tf = nullptr; }
      tf = TFile::Open(file.c_str(), "READ");
    }
    if (!tf || tf->IsZombie()) {
      std::cerr << "SKIP (dosya acilamadi / unreachable): " << file << std::endl;
      badFiles.push_back(file);
      if (tf) delete tf;
      continue;
    }
    tf->Close();
    delete tf;
    goodFiles.push_back(file);
  }
  std::cout << "\nDosya kontrolu tamamlandi: " << goodFiles.size() << " / "
            << inputFiles.size() << " dosya kullanilabilir." << std::endl;
  if (!badFiles.empty()) {
    std::cout << "Atlanan (skip edilen) dosya sayisi: " << badFiles.size() << std::endl;
    for (auto const &b : badFiles) std::cout << "   - " << b << std::endl;
  }
  return goodFiles;
}

// ============================================================================
// [END SHARED]
// ============================================================================


#include <set>

using json = nlohmann::json;

static const bool kBlindData = true;

// ============================================================================
// DATA-ONLY: GOLDEN JSON
// ============================================================================

std::map<unsigned int, std::vector<std::pair<unsigned int, unsigned int>>> goodLumis;

void LoadGoldenJSON(const std::string &filename) {
  goodLumis.clear();
  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cerr << "ERROR: Cannot open Golden JSON: " << filename << std::endl;
    exit(1);
  }
  json j;
  file >> j;
  for (auto &run : j.items()) {
    unsigned int runNumber = std::stoul(run.key());
    for (auto &lumi : run.value())
      goodLumis[runNumber].push_back({(unsigned int)lumi[0], (unsigned int)lumi[1]});
  }
  std::cout << "Golden JSON loaded. Runs = " << goodLumis.size() << std::endl;
}

bool PassGoldenJSON(unsigned int run, unsigned int lumi) {
  auto it = goodLumis.find(run);
  if (it == goodLumis.end()) return false;
  for (auto &range : it->second)
    if (lumi >= range.first && lumi <= range.second) return true;
  return false;
}

// ============================================================================
// MAIN
// ============================================================================
void makeTree_DATA_2024_KIT_jet_final_DeltaRcut(
    std::vector<std::string> inputFiles,
    bool isData = true,
    const std::string &KITDir = "/afs/cern.ch/user/n/nbostan/new_CMS/CMSSW_14_0_18/src/RoccoR/post2022E-update/",
    const std::string &BSCDir = "/afs/cern.ch/user/f/ferrico/cmsonly/JSON_VXBS/2024/",
    const std::string &RegularDir = "/afs/cern.ch/user/f/ferrico/cmsonly/JSON/2024/",
    int year = 2024,
    const std::string &jetJercJsonFile = "", // empty -> KITDir + "jet_jerc.json"
    // NOTE: full 2024 (no era splitting) -> use the FULL-YEAR golden JSON here,
    // not the 2024H-only file.
    const std::string &goldenJSON = "/cvmfs/cms-griddata.cern.ch/cat/metadata/DC/Collisions24/latest/2024H_Muon.json",
    const std::string &outFile = "/eos/user/n/nbostan/2025_Samples/tree_DATA_2024_loose_round2_KIT_JECfixed_DRmu_latest.root",
    const std::string &sampleName = "",
    const std::string &jetIdJsonFile = "/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/2026-07-16/jetid.json.gz") {

  gErrorIgnoreLevel = kError;
  gYear = year;

  // ---- corrections
  std::cout << "Loading BSC Muon ScaRe corrections..." << std::endl;
  muonCorrBSC = correction::CorrectionSet::from_file(BSCDir + "schemaV2.json");

  std::cout << "Loading Regular Muon ScaRe corrections..." << std::endl;
  muonCorrRegular = correction::CorrectionSet::from_file(RegularDir + "schemaV2.json");

  btagCorr = correction::CorrectionSet::from_file(KITDir + "btagging.json");
  if (!loadBtagWorkingPoints(btagCorr)) return;

  jetCorr = JetCorrections(KITDir + "jet_jerc.txt", KITDir + "jetid.txt",
                           KITDir + "jetvetomaps.txt");
  std::cout << "Jet corrections loaded (year rules: " << gYear << ")" << std::endl;

  const std::string jercJson =
      jetJercJsonFile.empty() ? (KITDir + "jet_jerc.json") : jetJercJsonFile;
  if (!loadJetJercRegionCorrections(jercJson)) return;

  if (!loadJetId(jetIdJsonFile)) return;

  // golden JSON must be loaded BEFORE the self-test (the test run is taken from it)
  if (isData) LoadGoldenJSON(goldenJSON);

  // [C] JEC self-test: rawFactor = 0 -> the input pT is already the raw pT.
  // The run number comes from the golden JSON (first run), so the L2L3Residual IOV
  // is one that really exists in the processed data.
  {
    const unsigned int testRun = goodLumis.empty() ? 385000u : goodLumis.begin()->first;
    const double testRaw = 40., testEta = 2.75;
    double c = getCorrectedJetPtRegionAware(testRaw, testEta, 0.5, 0., 0.5, 30., testRun,
                                            false, false, -1.);
    std::cout << "JEC self-test (DATA, run " << testRun << "): raw " << testRaw
              << ", eta " << testEta << " -> " << c << "  (expect ~74)" << std::endl;
  }

  // ---- input files: remove files listed more than once (a repeated file would
  //      double-count events)
  {
    std::set<std::string> seen;
    std::vector<std::string> uniq;
    for (auto const &f : inputFiles)
      if (seen.insert(f).second) uniq.push_back(f);
    if (uniq.size() != inputFiles.size())
      std::cout << "WARNING: " << (inputFiles.size() - uniq.size())
                << " duplicated file name(s) removed from the input list" << std::endl;
    inputFiles = uniq;
  }

  // ---- chain
  std::vector<std::string> readableFiles = FilterReadableFiles(inputFiles);
  if (readableFiles.empty()) {
    std::cerr << "ERROR: Hicbir dosya acilamadi, cikiliyor." << std::endl;
    return;
  }
  TChain chain("Events");
  for (auto const &file : readableFiles) chain.Add(file.c_str());
  std::cout << "Total entries = " << chain.GetEntries() << std::endl;

  if (!chain.GetBranch("Jet_btagUParTAK4B")) {
    std::cerr << "ERROR: branch Jet_btagUParTAK4B not found." << std::endl;
    return;
  }

  // ---- tree reader
  TTreeReader reader(&chain);

  TTreeReaderValue<UInt_t> run(reader, "run");
  TTreeReaderValue<UInt_t> luminosityBlock(reader, "luminosityBlock");
  TTreeReaderValue<ULong64_t> eventNumber(reader, "event");

  TTreeReaderValue<bool> Flag_goodVertices(reader, "Flag_goodVertices");
  TTreeReaderValue<bool> Flag_globalSuperTightHalo2016Filter(reader, "Flag_globalSuperTightHalo2016Filter");
  TTreeReaderValue<bool> Flag_EcalDeadCellTriggerPrimitiveFilter(reader, "Flag_EcalDeadCellTriggerPrimitiveFilter");
  TTreeReaderValue<bool> Flag_BadPFMuonFilter(reader, "Flag_BadPFMuonFilter");
  TTreeReaderValue<bool> Flag_BadPFMuonDzFilter(reader, "Flag_BadPFMuonDzFilter");
  TTreeReaderValue<bool> Flag_hfNoisyHitsFilter(reader, "Flag_hfNoisyHitsFilter");

  TTreeReaderValue<bool> HLT_IsoMu24(reader, "HLT_IsoMu24");

  bool hasEle = chain.GetBranch("Electron_pt") && chain.GetBranch("Electron_eta") &&
                chain.GetBranch("Electron_mvaIso_WP90");
  TTreeReaderArray<float> *Electron_pt = nullptr, *Electron_eta = nullptr;
  TTreeReaderArray<bool> *Electron_mvaIso_WP90 = nullptr;
  if (hasEle) {
    Electron_pt = new TTreeReaderArray<float>(reader, "Electron_pt");
    Electron_eta = new TTreeReaderArray<float>(reader, "Electron_eta");
    Electron_mvaIso_WP90 = new TTreeReaderArray<bool>(reader, "Electron_mvaIso_WP90");
  } else {
    std::cerr << "WARNING: Electron branches not found, passEleVeto = true" << std::endl;
  }

  bool hasPuppiMET = chain.GetBranch("PuppiMET_pt") && chain.GetBranch("PuppiMET_phi");
  bool hasPFMET = chain.GetBranch("MET_pt") && chain.GetBranch("MET_phi");
  TTreeReaderValue<float> *PuppiMET_pt = nullptr, *PuppiMET_phi = nullptr;
  TTreeReaderValue<float> *PFMET_pt = nullptr, *PFMET_phi = nullptr;
  if (hasPuppiMET) {
    PuppiMET_pt = new TTreeReaderValue<float>(reader, "PuppiMET_pt");
    PuppiMET_phi = new TTreeReaderValue<float>(reader, "PuppiMET_phi");
  } else {
    std::cerr << "WARNING: PuppiMET branches not found" << std::endl;
  }
  if (hasPFMET) {
    PFMET_pt = new TTreeReaderValue<float>(reader, "MET_pt");
    PFMET_phi = new TTreeReaderValue<float>(reader, "MET_phi");
  }

  TTreeReaderValue<int> nMuon(reader, "nMuon");
  TTreeReaderArray<float> Muon_pt(reader, "Muon_pt");
  TTreeReaderArray<float> Muon_eta(reader, "Muon_eta");
  TTreeReaderArray<float> Muon_phi(reader, "Muon_phi");
  TTreeReaderArray<float> Muon_iso(reader, "Muon_pfRelIso04_all");
  TTreeReaderArray<int> Muon_charge(reader, "Muon_charge");
  TTreeReaderArray<bool> Muon_mediumId(reader, "Muon_mediumId");

  // Muon_nTrackerLayers: needed for the ScaRe resolution (cb_params / poly_params bins).
  // Stored as Int_t in NanoAOD; fall back to UChar_t if the file says so.
  if (!chain.GetBranch("Muon_nTrackerLayers")) {
    std::cerr << "ERROR: branch Muon_nTrackerLayers not found (needed by ScaRe)." << std::endl;
    return;
  }
  TLeaf *leafNL = chain.GetLeaf("Muon_nTrackerLayers");
  const bool nLisUChar = leafNL && std::string(leafNL->GetTypeName()) == "UChar_t";
  std::cout << "Muon_nTrackerLayers type = " << (leafNL ? leafNL->GetTypeName() : "?") << std::endl;
  TTreeReaderArray<int> *Muon_nTrackerLayersI = nullptr;
  TTreeReaderArray<unsigned char> *Muon_nTrackerLayersU = nullptr;
  if (nLisUChar) Muon_nTrackerLayersU = new TTreeReaderArray<unsigned char>(reader, "Muon_nTrackerLayers");
  else           Muon_nTrackerLayersI = new TTreeReaderArray<int>(reader, "Muon_nTrackerLayers");
  auto getNL = [&](int i) -> int {
    return nLisUChar ? (int)(*Muon_nTrackerLayersU)[i] : (*Muon_nTrackerLayersI)[i];
  };

  bool hasFSR = chain.GetBranch("FsrPhoton_pt") && chain.GetBranch("FsrPhoton_eta") &&
                chain.GetBranch("FsrPhoton_phi") && chain.GetBranch("FsrPhoton_dROverEt2") &&
                chain.GetBranch("FsrPhoton_relIso03") && chain.GetBranch("FsrPhoton_muonIdx");
  TTreeReaderArray<float> *FsrPhoton_pt = nullptr, *FsrPhoton_eta = nullptr,
                          *FsrPhoton_phi = nullptr, *FsrPhoton_dROverEt2 = nullptr,
                          *FsrPhoton_relIso03 = nullptr;
  TTreeReaderArray<Short_t> *FsrPhoton_muonIdx = nullptr;
  if (hasFSR) {
    FsrPhoton_pt = new TTreeReaderArray<float>(reader, "FsrPhoton_pt");
    FsrPhoton_eta = new TTreeReaderArray<float>(reader, "FsrPhoton_eta");
    FsrPhoton_phi = new TTreeReaderArray<float>(reader, "FsrPhoton_phi");
    FsrPhoton_dROverEt2 = new TTreeReaderArray<float>(reader, "FsrPhoton_dROverEt2");
    FsrPhoton_relIso03 = new TTreeReaderArray<float>(reader, "FsrPhoton_relIso03");
    FsrPhoton_muonIdx = new TTreeReaderArray<Short_t>(reader, "FsrPhoton_muonIdx");
  }
  std::cout << "FSR photon branches = " << hasFSR << std::endl;

  bool hasBS = chain.GetBranch("Muon_bsConstrainedPt");
  TTreeReaderArray<float> *Muon_bsConstrainedPt = nullptr, *Muon_bsConstrainedChi2 = nullptr;
  if (hasBS) {
    Muon_bsConstrainedPt = new TTreeReaderArray<float>(reader, "Muon_bsConstrainedPt");
    Muon_bsConstrainedChi2 = new TTreeReaderArray<float>(reader, "Muon_bsConstrainedChi2");
  }
  std::cout << "BeamSpot branch = " << hasBS << std::endl;

  TTreeReaderValue<int> nTrigObj(reader, "nTrigObj");
  TTreeReaderArray<float> TrigObj_pt(reader, "TrigObj_pt");
  TTreeReaderArray<float> TrigObj_eta(reader, "TrigObj_eta");
  TTreeReaderArray<float> TrigObj_phi(reader, "TrigObj_phi");
  TTreeReaderArray<UShort_t> TrigObj_id(reader, "TrigObj_id");
  TTreeReaderArray<ULong64_t> *TrigObj_filterBits = nullptr;
  if (chain.GetBranch("TrigObj_filterBits"))
    TrigObj_filterBits = new TTreeReaderArray<ULong64_t>(reader, "TrigObj_filterBits");

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

  // ---- jet ID inputs (jetid.json). Jet_jetId is only stored as a cross-check.
  for (const char *bn : {"Jet_chHEF", "Jet_neHEF", "Jet_muEF"})
    if (!chain.GetBranch(bn)) {
      std::cerr << "ERROR: branch " << bn << " not found (needed by the jet ID)." << std::endl;
      return;
    }
  TTreeReaderArray<float> Jet_chHEF(reader, "Jet_chHEF");
  TTreeReaderArray<float> Jet_neHEF(reader, "Jet_neHEF");
  TTreeReaderArray<float> Jet_muEF(reader, "Jet_muEF");
  IntArrayReader Jet_chMultiplicity(reader, chain, "Jet_chMultiplicity");
  IntArrayReader Jet_neMultiplicity(reader, chain, "Jet_neMultiplicity");
  IntArrayReader Jet_nConstituents(reader, chain, "Jet_nConstituents");
  IntArrayReader Jet_jetId(reader, chain, "Jet_jetId");
  if (!Jet_jetId.ok)
    std::cout << "INFO: Jet_jetId is not stored in NanoAODv15 -> jet_jetId_nano = -1 "
                 "(jet ID is computed from jetid.json)" << std::endl;
  if (!Jet_chMultiplicity.ok || !Jet_neMultiplicity.ok || !Jet_nConstituents.ok) {
    std::cerr << "ERROR: unsupported/missing type of a jet ID branch (see types printed above)." << std::endl;
    return;
  }

  // ---- output
  TFile out(outFile.c_str(), "RECREATE");
  if (out.IsZombie()) {
    std::cerr << "ERROR: Cannot create output ROOT file: " << outFile << std::endl;
    return;
  }
  out.cd();
  TTree *tree = new TTree("HmmTree", "H->mumu loose-baseline tree (DATA)");
  tree->SetAutoSave(50000000); // ~50 MB'ta bir kaydet (crash olursa veri kurtarilabilir)

  UInt_t b_run = 0, b_lumi = 0;
  ULong64_t b_event = 0;
  Bool_t b_passGolden = false, b_passMETFilters = false, b_HLT = false, b_trigMatched = false;
  Bool_t b_passMuonEta = false, b_passMuonID = false, b_passMuonIso = false,
         b_passMuonPt = false, b_passMuonSel = false;
  Bool_t b_passMu3rdVeto = true;
  Int_t b_nExtraGoodMu = 0;
  Bool_t b_passJetVetoMap = false, b_passBVeto = false, b_selFull = false, b_selNoTrig = false;
  Bool_t b_inWindow_115_135 = false;

  std::string b_sample = sampleName;
  Bool_t b_passEleVeto = true;
  Int_t b_nEleVeto = 0;
  Float_t b_puppimet_pt = -1, b_puppimet_phi = -99, b_pfmet_pt = -1, b_pfmet_phi = -99;
  Int_t b_massRegion = -1;

  Bool_t b_cat_isVBF[2] = {false, false}, b_cat_isGGH[2] = {false, false};
  Int_t b_nJetsStored = 0, b_nJetsAna = 0, b_nMediumB = 0, b_nLooseB = 0;

  Float_t b_mu_pt_nano[2], b_mu_pt_bs[2], b_mu_chi2_bs[2], b_mu_pt_uncorr[2], b_mu_pt_corr[2],
      b_mu_eta[2], b_mu_phi[2], b_mu_iso[2], b_mu_trigDR[2];
  Int_t b_mu_charge[2];
  Bool_t b_mu_usedBS[2], b_mu_mediumId[2], b_mu_trigMatched[2];

  Float_t b_mm_mass_raw = 0, b_mm_mass_corr = 0, b_mm_mass_uncorr = 0, b_mm_pt_corr = 0,
          b_mm_pt_uncorr = 0, b_mm_eta_corr = 0, b_mm_phi_corr = 0;
  Float_t b_mm_mass_fsr = 0, b_mm_pt_fsr = 0;
  Float_t b_mu_fsr_pt[2] = {0, 0};
  Bool_t b_mu_hasFSR[2] = {false, false};
  Int_t b_nFSR = 0;
  Float_t b_mm_eta_fsr = -99, b_mm_phi_fsr = -99, b_mm_y_fsr = -99, b_mm_dR = -1;

  // single-muon pT variants: [variant][0 = leading, 1 = subleading]
  const int kNMuVar = 9;
  const char *kMuVarNames[9] = {"pt", "pt_nano", "pt_nano_corr", "pt_nano_corr_FSR",
                                "pt_raw_noCorr", "pt_raw_scale", "pt_raw_corr",
                                "pt_FSR_noCorr", "pt_FSR_corr"};
  Float_t b_muv[9][2] = {};
  // two leading analysis jets: [variable][0 = leading, 1 = subleading]
  const int kNJetVar = 7;
  const char *kJetVarNames[7] = {"pt", "pt_nano", "pt_raw", "eta", "phi", "mass", "btag"};
  Float_t b_jv[7][2];
  Float_t b_jj_dphi[2], b_jj_dR[2], b_jj_y[2], b_jj_zstar[2];
  Int_t b_cat_id[2] = {-1, -1};

  Bool_t b_jj_has[2];
  Float_t b_jj_mass[2], b_jj_pt[2], b_jj_eta[2], b_jj_deta[2];
  Int_t b_jj_idx1[2], b_jj_idx2[2];
  std::vector<int> v_jet_idTight, v_jet_idTightLV, v_jet_jetId_nano;

  std::vector<float> v_jet_pt, v_jet_pt_nano, v_jet_pt_raw, v_jet_eta, v_jet_phi, v_jet_mass,
      v_jet_btag, v_jet_chEmEF, v_jet_neEmEF, v_jet_dRmu;
  std::vector<int> v_jet_idOK, v_jet_vetoCand, v_jet_inVetoMap, v_jet_passRegion,
      v_jet_passDRmu, v_jet_passAnalysis;

  tree->Branch("run", &b_run, "run/i");
  tree->Branch("lumi", &b_lumi, "lumi/i");
  tree->Branch("event", &b_event, "event/l");
  tree->Branch("sample", &b_sample);

  tree->Branch("passGolden", &b_passGolden, "passGolden/O");
  tree->Branch("passMETFilters", &b_passMETFilters, "passMETFilters/O");
  tree->Branch("HLT_IsoMu24", &b_HLT, "HLT_IsoMu24/O");
  tree->Branch("trigMatched", &b_trigMatched, "trigMatched/O");

  tree->Branch("passMuonEta", &b_passMuonEta, "passMuonEta/O");
  tree->Branch("passMuonID", &b_passMuonID, "passMuonID/O");
  tree->Branch("passMuonIso", &b_passMuonIso, "passMuonIso/O");
  tree->Branch("passMuonPt", &b_passMuonPt, "passMuonPt/O");
  tree->Branch("passMuonSel", &b_passMuonSel, "passMuonSel/O");
  tree->Branch("nExtraGoodMu", &b_nExtraGoodMu, "nExtraGoodMu/I");
  tree->Branch("passMu3rdVeto", &b_passMu3rdVeto, "passMu3rdVeto/O");
  tree->Branch("passEleVeto", &b_passEleVeto, "passEleVeto/O");
  tree->Branch("nEleVeto", &b_nEleVeto, "nEleVeto/I");

  tree->Branch("passJetVetoMap", &b_passJetVetoMap, "passJetVetoMap/O");
  tree->Branch("passBVeto", &b_passBVeto, "passBVeto/O");
  tree->Branch("sel_full", &b_selFull, "sel_full/O");       // incl. trigger matching
  tree->Branch("sel_noTrigMatch", &b_selNoTrig, "sel_noTrigMatch/O");

  tree->Branch("puppimet_pt", &b_puppimet_pt, "puppimet_pt/F");
  tree->Branch("puppimet_phi", &b_puppimet_phi, "puppimet_phi/F");
  tree->Branch("pfmet_pt", &b_pfmet_pt, "pfmet_pt/F");
  tree->Branch("pfmet_phi", &b_pfmet_phi, "pfmet_phi/F");
  tree->Branch("inWindow_115_135", &b_inWindow_115_135, "inWindow_115_135/O");
  tree->Branch("massRegion", &b_massRegion, "massRegion/I");

  tree->Branch("cat_isVBF", b_cat_isVBF, "cat_isVBF[2]/O");
  tree->Branch("cat_isGGH", b_cat_isGGH, "cat_isGGH[2]/O");

  tree->Branch("nJetsStored", &b_nJetsStored, "nJetsStored/I");
  tree->Branch("nJetsAna", &b_nJetsAna, "nJetsAna/I");
  tree->Branch("nMediumB", &b_nMediumB, "nMediumB/I");
  tree->Branch("nLooseB", &b_nLooseB, "nLooseB/I");

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

  tree->Branch("mm_mass_raw", &b_mm_mass_raw, "mm_mass_raw/F");
  tree->Branch("mm_mass_corr", &b_mm_mass_corr, "mm_mass_corr/F"); // ScaRe only, NO FSR (use mm_mass_fsr)
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
  tree->Branch("mm_eta_fsr", &b_mm_eta_fsr, "mm_eta_fsr/F");
  tree->Branch("mm_phi_fsr", &b_mm_phi_fsr, "mm_phi_fsr/F");
  tree->Branch("mm_y_fsr", &b_mm_y_fsr, "mm_y_fsr/F");
  tree->Branch("mm_dR", &b_mm_dR, "mm_dR/F");
  for (int v = 0; v < kNMuVar; v++)
    for (int i = 0; i < 2; i++) {
      std::string n = "mu" + std::to_string(i + 1) + "_" + kMuVarNames[v];
      tree->Branch(n.c_str(), &b_muv[v][i], (n + "/F").c_str());
    }
  for (int v = 0; v < kNJetVar; v++)
    for (int k = 0; k < 2; k++) {
      std::string n = "jet" + std::to_string(k + 1) + "_" + kJetVarNames[v];
      tree->Branch(n.c_str(), &b_jv[v][k], (n + "/F").c_str());
    }
  tree->Branch("jj_dphi", b_jj_dphi, "jj_dphi[2]/F");
  tree->Branch("jj_dR", b_jj_dR, "jj_dR[2]/F");
  tree->Branch("jj_y", b_jj_y, "jj_y[2]/F");
  tree->Branch("jj_zstar", b_jj_zstar, "jj_zstar[2]/F");
  tree->Branch("cat_id", b_cat_id, "cat_id[2]/I");

  tree->Branch("jj_has", b_jj_has, "jj_has[2]/O");
  tree->Branch("jj_mass", b_jj_mass, "jj_mass[2]/F");
  tree->Branch("jj_pt", b_jj_pt, "jj_pt[2]/F");
  tree->Branch("jj_eta", b_jj_eta, "jj_eta[2]/F");
  tree->Branch("jj_deta", b_jj_deta, "jj_deta[2]/F");
  tree->Branch("jj_idx1", b_jj_idx1, "jj_idx1[2]/I");
  tree->Branch("jj_idx2", b_jj_idx2, "jj_idx2[2]/I");

  tree->Branch("jet_pt", &v_jet_pt);
  tree->Branch("jet_pt_nano", &v_jet_pt_nano);
  tree->Branch("jet_pt_raw", &v_jet_pt_raw);
  tree->Branch("jet_eta", &v_jet_eta);
  tree->Branch("jet_phi", &v_jet_phi);
  tree->Branch("jet_mass", &v_jet_mass);
  tree->Branch("jet_btag", &v_jet_btag);
  tree->Branch("jet_chEmEF", &v_jet_chEmEF);
  tree->Branch("jet_neEmEF", &v_jet_neEmEF);
  tree->Branch("jet_dRmu", &v_jet_dRmu);
  tree->Branch("jet_idOK", &v_jet_idOK);
  tree->Branch("jet_idTight", &v_jet_idTight);
  tree->Branch("jet_idTightLV", &v_jet_idTightLV);
  tree->Branch("jet_jetId_nano", &v_jet_jetId_nano);
  tree->Branch("jet_vetoCand", &v_jet_vetoCand);
  tree->Branch("jet_inVetoMap", &v_jet_inVetoMap);
  tree->Branch("jet_passRegion", &v_jet_passRegion);
  tree->Branch("jet_passDRmu", &v_jet_passDRmu);
  tree->Branch("jet_passAnalysis", &v_jet_passAnalysis);

  // ---- reader sanity check (branch tip uyusmazligini erken yakala)
  if (reader.SetEntry(0) != TTreeReader::kEntryValid) {
    std::cerr << "ERROR: reader setup failed (branch type mismatch?)" << std::endl;
    return;
  }
  reader.Restart();

  Long64_t nTotal = 0, nGolden = 0, nMETFilter = 0, nTrigger = 0, nTwoMuon = 0, nOS = 0,
           nPre = 0, nBlind = 0, nWritten = 0, nSelFull = 0, nVBF0 = 0, nVBF1 = 0;
  // VBF debug counters among sel_full events, per dijet mode
  Long64_t nJJhas[2] = {0, 0}, nJJmass[2] = {0, 0}, nJJdeta[2] = {0, 0}, nJJbfail[2] = {0, 0};

  // ========================================================================
  // EVENT LOOP
  // ========================================================================
  while (reader.Next()) {
    nTotal++;

    bool passGolden = true;
    if (isData) {
      passGolden = PassGoldenJSON(*run, *luminosityBlock);
      if (!passGolden) continue;
      nGolden++;
    }

    bool passMET = (*Flag_goodVertices) && (*Flag_globalSuperTightHalo2016Filter) &&
                   (*Flag_EcalDeadCellTriggerPrimitiveFilter) && (*Flag_BadPFMuonFilter) &&
                   (*Flag_BadPFMuonDzFilter) && (*Flag_hfNoisyHitsFilter);
    if (kRequireMETFilters && !passMET) continue;
    nMETFilter++;

    bool hlt = *HLT_IsoMu24;
    if (kRequireHLT && !hlt) continue;
    nTrigger++;

    // ---- two leading (raw pT) muons, opposite sign
    if (*nMuon < 2) continue;
    nTwoMuon++;
    if (kRequireOS && Muon_charge[0] * Muon_charge[1] >= 0) continue;
    nOS++;

    // ---- RAW preselection: raw pT > 20, 50 < m(mumu, raw) < 200
    if (Muon_pt[0] <= kRawMuPtMin || Muon_pt[1] <= kRawMuPtMin) continue;
    TLorentzVector r0, r1;
    r0.SetPtEtaPhiM(Muon_pt[0], Muon_eta[0], Muon_phi[0], kMuMass);
    r1.SetPtEtaPhiM(Muon_pt[1], Muon_eta[1], Muon_phi[1], kMuMass);
    const double mRaw = (r0 + r1).M();
    if (mRaw <= kMassMin || mRaw >= kMassMax) continue;
    nPre++;

    // ---- DATA blinding on RAW m(mumu)
    if (kBlindData && mRaw >= kBlindLo && mRaw <= kBlindHi) {
      nBlind++;
      continue;
    }

    // ---- third muon veto: exactly 2 good muons
    int nExtraGoodMu = 0;
    for (int k = 2; k < *nMuon; k++)
      if (Muon_pt[k] > kRawMuPtMin && fabs(Muon_eta[k]) < kMuEtaMax &&
          Muon_mediumId[k] && Muon_iso[k] < kMuIsoMax)
        nExtraGoodMu++;

    // ====================================================================
    // MUONS (index == NanoAOD index until the final re-sorting)
    // ====================================================================
    MuInfo mu[2];
    for (int i = 0; i < 2; i++) {
      mu[i].nanoIdx = i;
      mu[i].nL = getNL(i);
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
      mu[i].ptCorr = scaReFull(true, mu[i].ptUsed, mu[i].eta, mu[i].phi, mu[i].charge, mu[i].nL,
                               *eventNumber, *luminosityBlock, mu[i].usedBS);

      double minDR = minTrigObjDR(*nTrigObj, TrigObj_pt, TrigObj_eta, TrigObj_phi,
                                  TrigObj_id, TrigObj_filterBits, mu[i].eta, mu[i].phi);
      mu[i].trigDR = minDR;
      mu[i].trigMatched = (mu[i].ptNano >= kTrigMuRawPtMin && minDR >= 0. && minDR < kTrigMatchDR);
    }

    // ---- ScaRe-corrected and FSR-corrected 4-vectors
    TLorentzVector mC[2], mF[2], bestPhoton[2];
    bool muHasFSR[2] = {false, false};
    for (int i = 0; i < 2; i++)
      mC[i].SetPtEtaPhiM(mu[i].ptCorr, mu[i].eta, mu[i].phi, kMuMass);

    if (hasFSR)
      findFSRPhotons(*FsrPhoton_pt, *FsrPhoton_eta, *FsrPhoton_phi, *FsrPhoton_dROverEt2,
                     *FsrPhoton_relIso03, *FsrPhoton_muonIdx, mu, bestPhoton, muHasFSR);

    for (int i = 0; i < 2; i++) {
      mF[i] = mC[i];
      if (muHasFSR[i]) mF[i] += bestPhoton[i];
    }
    const int nFSRUsed = (int)muHasFSR[0] + (int)muHasFSR[1];

    // ---- single-muon pT variants (before the re-sort, so they travel with mu[i])
    for (int i = 0; i < 2; i++) {
      mu[i].ptNanoCorr = scaReFull(true, mu[i].ptNano, mu[i].eta, mu[i].phi, mu[i].charge, mu[i].nL,
                                    *eventNumber, *luminosityBlock, false);
      mu[i].ptNanoCorrFSR = ptWithFSR(mu[i].ptNanoCorr, mu[i].eta, mu[i].phi, bestPhoton[i], muHasFSR[i]);
      mu[i].ptRawScale = scaReScaleOnly(true, mu[i].ptUsed, mu[i].eta, mu[i].phi, mu[i].charge, mu[i].usedBS);
      mu[i].ptFSRnoCorr = ptWithFSR(mu[i].ptUsed, mu[i].eta, mu[i].phi, bestPhoton[i], muHasFSR[i]);
      mu[i].ptFSRcorr = mF[i].Pt();
    }

    // ---- re-sort by pT AFTER ScaRe + FSR: index 0 = leading
    if (mF[1].Pt() > mF[0].Pt()) {
      std::swap(mu[0], mu[1]);
      std::swap(mC[0], mC[1]);
      std::swap(mF[0], mF[1]);
      std::swap(muHasFSR[0], muHasFSR[1]);
    }

    bool trigMatched = mu[0].trigMatched || mu[1].trigMatched;

    bool passEta = true, passID = true, passIso = true;
    for (int i = 0; i < 2; i++) {
      if (fabs(mu[i].eta) > kMuEtaMax) passEta = false;
      if (!mu[i].mediumId) passID = false;
      if (mu[i].iso > kMuIsoMax) passIso = false;
    }
    // [D] + NanoAOD pT upper sanity cut (flag only, event is kept)
    bool passPt = (mF[0].Pt() > kMu1PtMin) && (mF[1].Pt() > kMu2PtMin) &&
                  (mu[0].ptNano < kMuPtMaxNano) && (mu[1].ptNano < kMuPtMaxNano);
    bool passMuonSel = passEta && passID && passIso && passPt;

    TLorentzVector mU0, mU1;
    mU0.SetPtEtaPhiM(mu[0].ptUsed, mu[0].eta, mu[0].phi, kMuMass);
    mU1.SetPtEtaPhiM(mu[1].ptUsed, mu[1].eta, mu[1].phi, kMuMass);

    TLorentzVector dimuon = mC[0] + mC[1];
    TLorentzVector dimuonU = mU0 + mU1;
    TLorentzVector dimuonFSR = mF[0] + mF[1];
    const double mll = dimuon.M();          // ScaRe only, no FSR
    const double mllFSR = dimuonFSR.M();    // official m(mumu): ScaRe + FSR

    nWritten++;

    // ====================================================================
    // JETS
    // ====================================================================
    std::vector<JetInfo> jets;
    bool vetoEvent = false;

    for (int j = 0; j < *nJet; j++) {
      double corrJetPt = getCorrectedJetPtRegionAware(
          Jet_pt[j], Jet_eta[j], Jet_phi[j], Jet_rawFactor[j], Jet_area[j], *Rho, *run,
          false, false, -1.);
      if (corrJetPt < kJetPtStore) continue;

      // [B] scale the jet mass with the same factor as the pT
      // (Jet_mass is consistent with the NanoAOD-corrected Jet_pt)
      const double jecScale = (Jet_pt[j] > 0.) ? corrJetPt / Jet_pt[j] : 1.;

      JetInfo ji;
      ji.p4.SetPtEtaPhiM(corrJetPt, Jet_eta[j], Jet_phi[j], Jet_mass[j] * jecScale);
      ji.ptNano = Jet_pt[j];
      ji.ptRaw = Jet_pt[j] * (1.0 - Jet_rawFactor[j]);
      ji.btag = Jet_btagUParTAK4B[j];
      ji.chEmEF = Jet_chEmEF[j];
      ji.neEmEF = Jet_neEmEF[j];

      // [F] jet ID from jetid.json (AK4PUPPI)
      const int idBits = evalJetId(Jet_eta[j], Jet_chHEF[j], Jet_neHEF[j], Jet_chEmEF[j],
                                   Jet_neEmEF[j], Jet_muEF[j], Jet_chMultiplicity[j],
                                   Jet_neMultiplicity[j], Jet_chMultiplicity[j] + Jet_neMultiplicity[j]);
      ji.idTight = (idBits & 1) != 0;
      ji.idTightLV = (idBits & 2) != 0;
      ji.jetIdNano = Jet_jetId[j];
      ji.idOK = ji.idTight;  // "tight jetID" of the agreed selection
      ji.vetoCand = ji.idTightLV && (Jet_chEmEF[j] + Jet_neEmEF[j]) < 0.9;
      ji.inVetoMap = false;
      if (ji.vetoCand && !jetCorr.passJetVetoMap(Jet_eta[j], Jet_phi[j])) {
        ji.inVetoMap = true;
        vetoEvent = true;
      }

      // [E] dR(jet, mu) >= 0.4 is part of the analysis jet selection
      ji.dRmu = std::min(ji.p4.DeltaR(mC[0]), ji.p4.DeltaR(mC[1]));
      ji.passDRmu = (ji.dRmu >= kJetDRMuMin);
      ji.passRegion = passRegionJetCut(corrJetPt, Jet_eta[j]);
      ji.passAnalysis = (corrJetPt > kJetPtMin) && (fabs(Jet_eta[j]) < kJetEtaMax) &&
                        ji.idOK && ji.passRegion && ji.passDRmu;
      jets.push_back(ji);
    }

    // re-sort after JEC
    std::sort(jets.begin(), jets.end(),
              [](const JetInfo &a, const JetInfo &b) { return a.p4.Pt() > b.p4.Pt(); });

    std::vector<int> ana;
    int nMediumB = 0, nLooseB = 0;
    for (size_t k = 0; k < jets.size(); k++) {
      if (!jets[k].passAnalysis) continue;
      ana.push_back((int)k);
      if (fabs(jets[k].p4.Eta()) >= kBtagEtaMax) continue; // b-veto only in the tracker
      if (jets[k].btag > gBtagWP_M) nMediumB++;
      if (jets[k].btag > gBtagWP_L) nLooseB++;
    }
    bool passBVeto = (nMediumB == 0 && nLooseB < 2);

    // ---- dijet + categories (both pair definitions)
    for (int m = 0; m < 2; m++) {
      DijetInfo d = evalDijet(jets, ana, m, passBVeto);
      b_jj_has[m] = d.has; b_jj_idx1[m] = d.i1; b_jj_idx2[m] = d.i2;
      b_jj_mass[m] = d.mass; b_jj_pt[m] = d.pt; b_jj_eta[m] = d.eta; b_jj_deta[m] = d.deta;
      b_jj_dphi[m] = d.dphi; b_jj_dR[m] = d.dR; b_jj_y[m] = d.y;
      b_jj_zstar[m] = -99;
      if (d.has) {
        const double y1 = jets[d.i1].p4.Rapidity(), y2 = jets[d.i2].p4.Rapidity();
        const double dy = fabs(y1 - y2);
        if (dy > 0.) b_jj_zstar[m] = (dimuonFSR.Rapidity() - 0.5 * (y1 + y2)) / dy;
      }
      b_cat_id[m] = !passBVeto ? -1 : (d.vbf ? 1 : 0); // -1 b-veto fail, 0 ggH, 1 VBF
      b_cat_isVBF[m] = d.vbf;
      b_cat_isGGH[m] = passBVeto && !d.vbf; // ggH = every non-VBF event (0/1 jets included)
    }

    // ---- electron veto
    int nEleVeto = hasEle ? countVetoElectrons(*Electron_pt, *Electron_eta, *Electron_mvaIso_WP90) : 0;

    // ====================================================================
    // FILL
    // ====================================================================
    b_run = *run;
    b_lumi = *luminosityBlock;
    b_event = *eventNumber;

    b_passGolden = passGolden;
    b_passMETFilters = passMET;
    b_HLT = hlt;
    b_trigMatched = trigMatched;

    b_passMuonEta = passEta;
    b_passMuonID = passID;
    b_passMuonIso = passIso;
    b_passMuonPt = passPt;
    b_passMuonSel = passMuonSel;
    b_nExtraGoodMu = nExtraGoodMu;
    b_passMu3rdVeto = (nExtraGoodMu == 0);
    b_nEleVeto = nEleVeto;
    b_passEleVeto = (nEleVeto == 0);

    b_passJetVetoMap = !vetoEvent;
    b_passBVeto = passBVeto;

    b_selNoTrig = passMuonSel && b_passMu3rdVeto && b_passEleVeto && !vetoEvent;
    b_selFull = trigMatched && b_selNoTrig;
    if (b_selFull) {
      nSelFull++;
      if (b_cat_isVBF[0]) nVBF0++;
      if (b_cat_isVBF[1]) nVBF1++;
      // VBF debug counters (per dijet mode)
      for (int m = 0; m < 2; m++) {
        if (!b_jj_has[m]) continue;
        nJJhas[m]++;
        if (b_jj_mass[m] > 400.) {
          nJJmass[m]++;
          if (b_jj_deta[m] > 2.5) {
            nJJdeta[m]++;
            if (!b_passBVeto) nJJbfail[m]++;
          }
        }
      }
    }

    // regions from the OFFICIAL m(mumu) = ScaRe + FSR (BSC/Nano pT)
    b_inWindow_115_135 = (mllFSR > 115.0 && mllFSR < 135.0);
    b_massRegion = (mllFSR < 115.0) ? 0 : ((mllFSR <= 135.0) ? 1 : 2);

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
      b_mu_fsr_pt[i] = mF[i].Pt();
      b_mu_hasFSR[i] = muHasFSR[i];
    }

    for (int i = 0; i < 2; i++) {
      const float vals[9] = {mu[i].ptNano, mu[i].ptNano, mu[i].ptNanoCorr, mu[i].ptNanoCorrFSR,
                             mu[i].ptUsed, mu[i].ptRawScale, mu[i].ptCorr,
                             mu[i].ptFSRnoCorr, mu[i].ptFSRcorr};
      for (int v = 0; v < kNMuVar; v++) b_muv[v][i] = vals[v];
    }

    b_mm_mass_raw = mRaw;
    b_mm_mass_corr = mll;
    b_mm_mass_uncorr = dimuonU.M();
    b_mm_pt_corr = dimuon.Pt();
    b_mm_pt_uncorr = dimuonU.Pt();
    b_mm_eta_corr = dimuon.Eta();
    b_mm_phi_corr = dimuon.Phi();
    b_mm_mass_fsr = mllFSR;
    b_mm_pt_fsr = dimuonFSR.Pt();
    b_nFSR = nFSRUsed;
    b_mm_eta_fsr = dimuonFSR.Eta();
    b_mm_phi_fsr = dimuonFSR.Phi();
    b_mm_y_fsr = dimuonFSR.Rapidity();
    b_mm_dR = mC[0].DeltaR(mC[1]);

    // two leading analysis jets (JEC-corrected pT; Nano and raw pT for comparison)
    for (int k = 0; k < 2; k++) {
      if ((int)ana.size() > k) {
        const JetInfo &ji = jets[ana[k]];
        const float vals[7] = {(float)ji.p4.Pt(), ji.ptNano, ji.ptRaw, (float)ji.p4.Eta(),
                               (float)ji.p4.Phi(), (float)ji.p4.M(), ji.btag};
        for (int v = 0; v < kNJetVar; v++) b_jv[v][k] = vals[v];
      } else {
        for (int v = 0; v < kNJetVar; v++) b_jv[v][k] = -99.f;
      }
    }

    v_jet_pt.clear(); v_jet_pt_nano.clear(); v_jet_pt_raw.clear(); v_jet_eta.clear();
    v_jet_phi.clear(); v_jet_mass.clear(); v_jet_btag.clear(); v_jet_chEmEF.clear();
    v_jet_neEmEF.clear(); v_jet_dRmu.clear(); v_jet_idOK.clear(); v_jet_vetoCand.clear();
    v_jet_inVetoMap.clear(); v_jet_passRegion.clear(); v_jet_passDRmu.clear();
    v_jet_passAnalysis.clear();

    v_jet_idTight.clear(); v_jet_idTightLV.clear(); v_jet_jetId_nano.clear();

    for (auto const &ji : jets) {
      v_jet_pt.push_back(ji.p4.Pt());
      v_jet_pt_nano.push_back(ji.ptNano);
      v_jet_pt_raw.push_back(ji.ptRaw);
      v_jet_eta.push_back(ji.p4.Eta());
      v_jet_phi.push_back(ji.p4.Phi());
      v_jet_mass.push_back(ji.p4.M());
      v_jet_btag.push_back(ji.btag);
      v_jet_chEmEF.push_back(ji.chEmEF);
      v_jet_neEmEF.push_back(ji.neEmEF);
      v_jet_dRmu.push_back(ji.dRmu);
      v_jet_idOK.push_back(ji.idOK);
      v_jet_idTight.push_back(ji.idTight);
      v_jet_idTightLV.push_back(ji.idTightLV);
      v_jet_jetId_nano.push_back(ji.jetIdNano);
      v_jet_vetoCand.push_back(ji.vetoCand);
      v_jet_inVetoMap.push_back(ji.inVetoMap);
      v_jet_passRegion.push_back(ji.passRegion);
      v_jet_passDRmu.push_back(ji.passDRmu);
      v_jet_passAnalysis.push_back(ji.passAnalysis);
    }

    b_puppimet_pt = hasPuppiMET ? **PuppiMET_pt : -1.f;
    b_puppimet_phi = hasPuppiMET ? **PuppiMET_phi : -99.f;
    b_pfmet_pt = hasPFMET ? **PFMET_pt : -1.f;
    b_pfmet_phi = hasPFMET ? **PFMET_phi : -99.f;

    tree->Fill();
  }

  // ========================================================================
  std::cout << "\nCutflow (DATA):"
            << "\n  total                    = " << nTotal
            << "\n  golden JSON              = " << nGolden
            << "\n  MET filters              = " << nMETFilter
            << "\n  HLT_IsoMu24              = " << nTrigger
            << "\n  nMuon >= 2               = " << nTwoMuon
            << "\n  opposite sign            = " << nOS
            << "\n  raw pT>20, 50<mraw<200   = " << nPre
            << "\n  blinded (raw 115-135)    = " << nBlind << "  (removed)"
            << "\n  tree entries             = " << nWritten
            << "\n  of which sel_full        = " << nSelFull
            << "\n    sel_full && VBF[0] / [1]            = " << nVBF0 << " / " << nVBF1
            << "\n  VBF steps among sel_full (has dijet / mjj>400 / +deta>2.5 / +bveto-fail):"
            << "\n    mode 0 (two leading jets) = " << nJJhas[0] << " / " << nJJmass[0] << " / "
            << nJJdeta[0] << " / " << nJJbfail[0]
            << "\n    mode 1 (max mjj pair)     = " << nJJhas[1] << " / " << nJJmass[1] << " / "
            << nJJdeta[1] << " / " << nJJbfail[1]
            << "\n  JEC-chain failures       = " << gResidualFail << std::endl;

  out.cd();
  tree->Write();
  TNamed("METCollection",
         "PuppiMET_pt/phi stored as puppimet_*; PF MET stored as pfmet_*; no re-propagation of JEC")
      .Write();
  TNamed("sampleName", sampleName.c_str()).Write();
  TNamed("jetEtaMax", std::to_string(kJetEtaMax).c_str()).Write();
  TNamed("blinding", "events with RAW (Muon_pt) m(mumu) in [115,135] GeV not written").Write();
  TNamed("jetJEC", "DATA jets: raw pT -> L1FastJet -> L2Relative -> L3Absolute -> L2L3Residual (jet_jerc.json); jet mass scaled with pT")
      .Write();
  TNamed("jetSelection", "pT>30 (JEC), |eta|<4.7, tight jetID (jetid.json), region rules, dR(jet,mu)>=0.4")
      .Write();
  out.Close();

  delete Muon_nTrackerLayersI;
  delete Muon_nTrackerLayersU;
  delete Muon_bsConstrainedPt;
  delete Muon_bsConstrainedChi2;
  delete TrigObj_filterBits;

  std::cout << "\n====================================\n"
            << " DATA TREE-MAKER FINISHED\n Output: " << outFile
            << "\n====================================\n";
}
