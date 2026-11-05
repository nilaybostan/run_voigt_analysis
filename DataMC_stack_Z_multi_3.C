#include "TFile.h"
#include "TTree.h"
#include "TNamed.h"
#include "TH1D.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLatex.h"
#include "TObjArray.h"
#include "TObjString.h"
#include <cmath>
#include <cstdio>
#include <vector>
#include <string>
#include <map>
#include <algorithm>

// ----------------------------- AYARLAR --------------------------------------
const std::string kDir      = "/eos/user/n/nbostan/2025_Samples/";
const std::string kTreeName = "HmmTree";
const std::string kOut      = "plots_stack/";

const std::string kDataPat  = kDir + "tree_DATA_2024H_final_Muon01_chunk*.root";

struct Group {
  std::string name;
  std::vector<std::string> pats;
  int color;
};
std::vector<Group> kGroups = {
  // stack sirasi: ilk eleman en altta, son eleman en ustte
  {"VBF H", { kDir + "tree_VBFH_2024H_final_nominal.root" },              kRed - 4},
  {"ggH",               { kDir + "tree_GGH_2024H_final_nominal.root" },               kOrange + 1},
  {"Diboson",           { kDir + "tree_ZZto2L2Q_2024H_final_chunk*.root",
                          kDir + "tree_ZZto2L2Nu_2024H_final_chunk*.root",
                          kDir + "tree_WZto3LNu_2024H_final_chunk1.root",
                          kDir + "tree_WZto2L2Q_2024H_final_chunk1.root",
                          kDir + "tree_ZZto4L_2024H_final_chunk*.root" },               kViolet + 1},
  {"EWK",               { kDir + "tree_EWK_2024H_final_chunk*.root" },                kGreen + 2},
  {"t#bar{t}",          { kDir + "tree_TTto2L2Nu_2024H_final_chunk*.root" },          kYellow + 1},
  {"DY",                { kDir + "tree_DY_2024H_final_chunk1_*.root" },               kAzure + 1}
};

const char* kMCWeight   = "w_total";   // sync agirligi (b-tag SF icermez)

// Ortak secim (data ve MC ayni): baseline + sel_full
const std::string kBaseSel  = "passMETFilters&&HLT_IsoMu24&&sel_full";
const std::string kSelMC    = kBaseSel;
const std::string kSelData  = kBaseSel + "&&passGolden";

// Z bolgesi kutle penceresi (resmi m(mumu) = mm_mass_fsr). Kapatmak icin false yap.
const bool   kApplyZWindow = true;
const double kZMassLo      = 70.;
const double kZMassHi      = 110.;
// Sinyal penceresi (sadece MC)
const double kSRMassLo     = 115.;
const double kSRMassHi     = 135.;
// Blind kenarlarinda (RAW vs FSR kutle farki) maskelenen ek genislik [GeV]; 0 yaparsan kapanir
const double kBlindMargin  = 1.0;

const bool   kNormToData = true;
const bool   kUseInclusiveK = false; // true: ggH/VBF plotlarinda k = ayni degiskenin inclusive k'si
const bool   kRatioLog   = false;    // false: lineer ratio (aralik kategoriye gore, bkz. MakeCats)
const double kYMin       = 1e-4;     // tum plotlarda ana pad y minimum
const int    kJJIdx      = 0;        // jj_* dizilerinde index: 0 = iki leading jet, 1 = max mjj cifti
const char*  kLumiLabel  = "2024 Data H, 5.44 fb^{-1} (13.6 TeV)";

// ------------------------- KATEGORI TANIMLARI -------------------------------
struct Cat {
  std::string tag;     // dosya/histogram adi oneki ("" = inclusive)
  std::string cut;     // ek kesim
  std::string label;   // plotta gosterilen etiket ("" = yok)
  int    rebin;        // dusuk istatistik icin binleri birlestir (bolen bulunamazsa degismez)
  double rLo, rHi;     // ratio pad araligi
};

std::vector<Cat> MakeCats() {
  const std::string j = std::to_string(kJJIdx);
  std::vector<Cat> c;
  c.push_back({"",    "",                         "",                         1, 0.5, 1.5});  // inclusive
  c.push_back({"ggH", "cat_isGGH[" + j + "]",     "ggH category (b-veto)",    1, 0.5, 1.5});
  c.push_back({"VBF", "cat_isVBF[" + j + "]",     "VBF category (b-veto)",    2, 0.0, 2.0});
  return c;
}

// ------------------------- DEGISKEN TANIMLARI -------------------------------
struct Var {
  std::string name;       // dosya adi icin
  std::string expr;       // cizilecek ifade
  std::string extraCut;   // ek kesim ("" olabilir)
  int    nb; double lo, hi;
  std::string xTitle;
  std::string unit;       // "GeV" veya ""
};

std::string JJ(const char* b) { return Form("%s[%d]", b, kJJIdx); }

std::vector<Var> MakeVars() {
  std::string jjCut = Form("jj_has[%d]", kJJIdx);
  std::vector<Var> v;
  // ---- Dimuon (resmi: ScaRe + FSR) ----
  v.push_back({"mmass", "mm_mass_fsr", "", 160, 70, 110, "m_{#mu#mu} [GeV]", "GeV"});
  v.push_back({"mmpt",  "mm_pt_fsr",   "", 100, 0, 200,  "p_{T}^{#mu#mu} [GeV]", "GeV"});
  v.push_back({"mmeta", "mm_eta_fsr",  "", 30, -3, 3,    "#eta^{#mu#mu}", ""});
  // ---- Muonlar: pT1 > 26, pT2 > 20, |eta| < 2.4 ----
  v.push_back({"mu1pt", "mu1_pt_FSR_corr", "", 124, 26, 150, "Leading muon p_{T} [GeV]", "GeV"});
  v.push_back({"mu2pt", "mu2_pt_FSR_corr", "", 80, 20, 100,  "Subleading muon p_{T} [GeV]", "GeV"});
  v.push_back({"mu1eta", "mu_eta[0]", "", 48, -2.4, 2.4, "Leading muon #eta", ""});
  v.push_back({"mu2eta", "mu_eta[1]", "", 48, -2.4, 2.4, "Subleading muon #eta", ""});
  // ---- Jetler: pT > 30, |eta| < 4.7 ----
  v.push_back({"jet1pt",  "jet1_pt",  "jet1_pt>0", 54, 30, 300, "Leading jet p_{T} [GeV]", "GeV"});
  v.push_back({"jet2pt",  "jet2_pt",  "jet2_pt>0", 54, 30, 300, "Subleading jet p_{T} [GeV]", "GeV"});
  v.push_back({"jet1eta", "jet1_eta", "jet1_pt>0", 47, -4.7, 4.7, "Leading jet #eta", ""});
  v.push_back({"jet2eta", "jet2_eta", "jet2_pt>0", 47, -4.7, 4.7, "Subleading jet #eta", ""});
  // ---- Dijet ----
  v.push_back({"jjmass", JJ("jj_mass"), jjCut, 75, 0, 1500, "m_{jj} [GeV]", "GeV"});
  v.push_back({"jjpt",   JJ("jj_pt"),   jjCut, 50, 0, 500,  "p_{T}^{jj} [GeV]", "GeV"});
  v.push_back({"jjeta",  JJ("jj_eta"),  jjCut, 40, -5, 5,   "#eta^{jj}", ""});
  v.push_back({"jjdeta", JJ("jj_deta"), jjCut, 18, 0, 9,    "|#Delta#eta_{jj}|", ""});   // YENI (VBF icin onemli)
  return v;
}

// bin sayisi r'ye (veya r'den kucuk bir bolene) tam bolunuyorsa binleri birlestir
Var Rebinned(Var v, int r) {
  for (int d = r; d > 1; --d)
    if (v.nb % d == 0) { v.nb /= d; break; }
  return v;
}
// ----------------------------------------------------------------------------

static int gCnt = 0;
static std::map<std::string, double> gKincl;   // inclusive k'ler (degisken adina gore)

std::vector<std::string> Glob(const std::string& pat) {
  std::vector<std::string> v;
  TString out = gSystem->GetFromPipe(Form("ls -1 %s 2>/dev/null", pat.c_str()));
  TObjArray* arr = out.Tokenize("\n");
  for (int i = 0; i < arr->GetEntries(); ++i) {
    std::string s = ((TObjString*)arr->At(i))->GetString().Data();
    if (!s.empty()) v.push_back(s);
  }
  delete arr;
  return v;
}

double ReadSumGen(const std::string& fn) {
  TFile* f = TFile::Open(fn.c_str());
  if (!f || f->IsZombie()) return 0;
  double S = 0;
  TNamed* nm = (TNamed*)f->Get("normalization");
  if (!nm || sscanf(nm->GetTitle(), "sumGenWeight = %lf", &S) != 1) S = 0;
  f->Close();
  return S;
}

void FillPattern(TH1D* h, const std::string& pat, const std::string& expr,
                 const std::string& cut, const char* w, bool isMC) {
  std::vector<std::string> files = Glob(pat);
  if (files.empty()) { printf("UYARI: dosya bulunamadi -> %s\n", pat.c_str()); return; }

  std::vector<double> frac(files.size(), 1.0);
  if (isMC && files.size() > 1) {
    double Stot = 0; bool ok = true;
    std::vector<double> S(files.size(), 0.);
    for (size_t i = 0; i < files.size(); ++i) {
      S[i] = ReadSumGen(files[i]);
      if (S[i] == 0) { ok = false; break; }
      Stot += S[i];
    }
    if (ok) for (size_t i = 0; i < files.size(); ++i) frac[i] = S[i] / Stot;
    else printf("UYARI: %s icin normalization okunamadi -> chunk duzeltmesi YOK\n", pat.c_str());
  }

  std::string c = (isMC && w) ? Form("(%s)*(%s)", cut.c_str(), w) : Form("(%s)", cut.c_str());
  for (size_t i = 0; i < files.size(); ++i) {
    TFile* f = TFile::Open(files[i].c_str());
    if (!f || f->IsZombie()) { printf("UYARI: acilamadi %s\n", files[i].c_str()); continue; }
    TTree* t = (TTree*)f->Get(kTreeName.c_str());
    if (!t) { printf("UYARI: %s icinde %s yok\n", files[i].c_str(), kTreeName.c_str()); f->Close(); continue; }
    TH1D* tmp = new TH1D(Form("tmp_%d", gCnt++), "", h->GetNbinsX(), h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
    tmp->Sumw2();
    t->Draw(Form("%s>>%s", expr.c_str(), tmp->GetName()), c.c_str(), "goff");
    tmp->Scale(frac[i]);
    h->Add(tmp);
    delete tmp;
    f->Close();
  }
  printf("   %-60s : %zu dosya, integral = %.6g\n", pat.c_str(), files.size(), h->Integral());
}

TH1D* MakeHist(const std::vector<std::string>& pats, const std::string& expr, const std::string& cut,
               const char* w, bool isMC, const char* name, int nb, double lo, double hi) {
  TH1D* h = new TH1D(name, "", nb, lo, hi);
  h->Sumw2();
  h->SetDirectory(nullptr);
  for (auto& p : pats) FillPattern(h, p, expr, cut, w, isMC);
  return h;
}

// s = yazi olcegi: TLatex NDC boyutu pad yuksekligine gore olceklenir. Ratio'lu plotta
// ana pad canvas'in %70'i, MC-only plotta %100'u -> MC-only icin s = 0.7 kullan.
void CMSLabel(double s = 1.0, const char* extra = "Work in progress") {
  TLatex l; l.SetNDC();
  l.SetTextFont(62); l.SetTextSize(0.05 * s); l.DrawLatex(0.14, 0.935, "CMS");
  l.SetTextFont(52); l.SetTextSize(0.04 * s); l.DrawLatex(0.23, 0.935, extra);
  l.SetTextFont(42); l.SetTextSize(0.04 * s); l.SetTextAlign(31); l.DrawLatex(0.95, 0.935, kLumiLabel);
}

// ---------------------------------------------------------------------------
// mcOnly = false : Data/MC stack + ratio (MC data'ya normalize)
// mcOnly = true  : sadece MC stack (blind bolge), data ve ratio yok
// massLo/massHi  : mm_mass_fsr penceresi (massLo<massHi ise uygulanir)
// cat            : kategori (cat.cut secime eklenir; cat.tag dosya adina girer)
// ---------------------------------------------------------------------------
// blindLo/blindHi : (blindLo<blindHi ise) bu aralikta data bin'leri silinir, MC normalizasyonu
//                   ve chi2 sadece blind disindaki (sideband) bin'lerden hesaplanir
void StackPlot(const Var& var, const std::string& id0, const char* lab,
               double massLo, double massHi, bool mcOnly, const Cat& cat,
               double blindLo = 0., double blindHi = 0.) {
  const std::string id = cat.tag.empty() ? id0 : (cat.tag + "_" + id0);  // dosya + ROOT nesne adlari
  const bool blind = (!mcOnly && blindLo < blindHi);
  printf("\n===== %s (%s) %s%s =====\n", id.c_str(), var.expr.c_str(), mcOnly ? "[MC only]" : "",
         cat.tag.empty() ? "" : (" [" + cat.tag + "]").c_str());

  std::string cutMC   = kSelMC;
  std::string cutData = kSelData;
  if (massLo < massHi) {
    std::string mc = Form("mm_mass_fsr>%g&&mm_mass_fsr<%g", massLo, massHi);
    cutMC += "&&" + mc; cutData += "&&" + mc;
  }
  if (!var.extraCut.empty()) { cutMC += "&&" + var.extraCut; cutData += "&&" + var.extraCut; }
  if (!cat.cut.empty())      { cutMC += "&&" + cat.cut;      cutData += "&&" + cat.cut; }

  TH1D* hD = nullptr;
  if (!mcOnly) {
    printf(" Data:\n");
    hD = MakeHist({kDataPat}, var.expr, cutData, nullptr, false, Form("hD_%s", id.c_str()),
                  var.nb, var.lo, var.hi);
    // blind bolge: tree'de RAW mass'e gore blind var, fsr mass'e gore kenarda sizan olay olabilir -> bin'leri temizle
    // RAW-FSR kutle farki yuzunden blind kenarlarinda (+-kBlindMargin) data eksik/bozuk:
    // bu bin'ler de maskelenir (cizilmez, normalizasyona ve chi2'ye girmez)
    if (blind)
      for (int i = 1; i <= hD->GetNbinsX(); ++i)
        if (hD->GetBinCenter(i) > blindLo - kBlindMargin && hD->GetBinCenter(i) < blindHi + kBlindMargin) {
          hD->SetBinContent(i, 0); hD->SetBinError(i, 0);
        }
  }

  std::vector<TH1D*> hG;
  std::vector<Group*> used;
  for (auto& g : kGroups) {
    printf(" %s:\n", g.name.c_str());
    TH1D* h = MakeHist(g.pats, var.expr, cutMC, kMCWeight, true, Form("h_%s_%d", id.c_str(), gCnt++),
                       var.nb, var.lo, var.hi);
    if (h->Integral() <= 0) { printf("   -> bos, grup atlandi: %s\n", g.name.c_str()); delete h; continue; }
    h->SetFillColor(g.color); h->SetLineColor(kBlack); h->SetLineWidth(1);
    hG.push_back(h); used.push_back(&g);
  }
  if (hG.empty() || (!mcOnly && hD->Integral() <= 0)) { printf("HATA: data veya MC bos (%s)\n", id.c_str()); return; }

  TH1D* hMC = (TH1D*)hG[0]->Clone(Form("hMC_%s", id.c_str()));
  hMC->SetDirectory(nullptr);
  for (size_t i = 1; i < hG.size(); ++i) hMC->Add(hG[i]);

  // (Hibrit RAW-kesimli stack kaldirildi: 115/135'te basamak olusturuyordu. Stack her yerde
  //  tam MC; blind kenarlarindaki bozuk data bin'leri yukarida maskelendi.)
  TH1D* hCmp = hMC;

  double k = 1.0;
  if (!mcOnly) {
    // normalizasyon: blind bolge haric (sadece sideband bin'leri)
    double nD = 0, nM = 0;
    for (int i = 1; i <= hD->GetNbinsX(); ++i) {
      double x = hD->GetBinCenter(i);
      if (blind && x > blindLo - kBlindMargin && x < blindHi + kBlindMargin) continue;
      nD += hD->GetBinContent(i); nM += hCmp->GetBinContent(i);
    }
    k = kNormToData ? nD / nM : 1.0;
    printf(" Data integral = %.6g, MC integral = %.6g, Data/MC = %.4f%s\n", nD, nM, nD / nM,
           kNormToData ? "  (MC data'ya normalize edildi)" : "");
    if (kNormToData && kUseInclusiveK) {
      if (cat.tag.empty()) gKincl[id0] = k;                       // inclusive: k'yi sakla
      else if (gKincl.count(id0)) {                               // kategori: inclusive k'yi kullan
        printf(" k(kategori) = %.4f -> inclusive k = %.4f kullaniliyor\n", k, gKincl[id0]);
        k = gKincl[id0];
      }
    }
  } else {
    printf(" MC integral = %.6g (lumi'ye normalize, data yok)\n", hMC->Integral());
  }
  for (auto h : hG) h->Scale(k);
  hMC->Scale(k);
  if (hCmp != hMC) hCmp->Scale(k);

  THStack* hs = new THStack(Form("hs_%s", id.c_str()), "");
  for (auto h : hG) hs->Add(h);

  const double s = mcOnly ? 0.7 : 1.0;   // yazi olcegi (bkz. CMSLabel)
  TCanvas* c = new TCanvas(Form("c_%s", id.c_str()), "", 800, 800);
  TPad* p1; TPad* p2 = nullptr;
  if (mcOnly) {
    p1 = new TPad(Form("p1_%s", id.c_str()), "", 0, 0, 1, 1);
    p1->SetBottomMargin(0.13); p1->SetLeftMargin(0.13); p1->SetTopMargin(0.07); p1->Draw();
  } else {
    p1 = new TPad(Form("p1_%s", id.c_str()), "", 0, 0.30, 1, 1);
    p1->SetBottomMargin(0.02); p1->SetLeftMargin(0.13); p1->SetTopMargin(0.07); p1->Draw();
    p2 = new TPad(Form("p2_%s", id.c_str()), "", 0, 0.00, 1, 0.30);
    p2->SetTopMargin(0.03); p2->SetBottomMargin(0.35); p2->SetLeftMargin(0.13); p2->Draw();
  }

  p1->cd();
  p1->SetLogy();
  hs->Draw("hist");
  double mx = hs->GetMaximum();
  if (!mcOnly) mx = std::max(mx, hD->GetMaximum());
  hs->SetMaximum(200 * mx);
  hs->SetMinimum(kYMin);
  double bw = (var.hi - var.lo) / var.nb;
  std::string yt = var.unit.empty() ? Form("Events / %.3g", bw)
                                    : Form("Events / %.3g %s", bw, var.unit.c_str());
  hs->GetYaxis()->SetTitle(yt.c_str());
  hs->GetYaxis()->SetTitleOffset(1.5);
  if (mcOnly) {
    hs->GetXaxis()->SetTitle(var.xTitle.c_str());
    hs->GetXaxis()->SetTitleSize(0.045); hs->GetXaxis()->SetLabelSize(0.040);
  } else {
    hs->GetXaxis()->SetLabelSize(0);
  }

  if (!mcOnly) {
    hD->SetMarkerStyle(20); hD->SetMarkerSize(0.7); hD->SetLineColor(kBlack); hD->SetMarkerColor(kBlack);
    hD->Draw("E same");
  }

  // iki sutunlu legend, ust kosede (stack ustundeki bos alanda kalir)
  TLegend* lg = new TLegend(0.50, 0.74, 0.93, 0.91);
  lg->SetNColumns(2);
  lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.034 * s);
  if (!mcOnly) lg->AddEntry(hD, "Data 2024H", "lep");
  for (int i = (int)hG.size() - 1; i >= 0; --i) lg->AddEntry(hG[i], used[i]->name.c_str(), "f");
  lg->Draw();

  // sol ust yazilar: bolge etiketi, (kategori etiketi), (blind notu) -- alt alta
  double yTxt = 0.85;
  TLatex l; l.SetNDC(); l.SetTextSize(0.036 * s); l.DrawLatex(0.16, yTxt, lab);
  yTxt -= 0.05;
  if (!cat.label.empty()) {
    TLatex lc; lc.SetNDC(); lc.SetTextFont(62); lc.SetTextSize(0.034 * s);
    lc.DrawLatex(0.16, yTxt, cat.label.c_str());
    yTxt -= 0.05;
  }
  CMSLabel(s);

  if (mcOnly) {
    TLatex lb; lb.SetNDC(); lb.SetTextSize(0.032 * s);
    lb.DrawLatex(0.16, yTxt, "Data blinded (MC only)");
  } else {
    // ratio
    p2->cd();
    if (kRatioLog) p2->SetLogy();
    TH1D* r = (TH1D*)hD->Clone(Form("hR_%s", id.c_str()));
    r->Divide(hCmp);
    r->SetTitle(""); r->SetStats(0);
    r->GetYaxis()->SetTitle("Data / MC");
    if (kRatioLog) { r->GetYaxis()->SetRangeUser(0.5, 2.0); r->GetYaxis()->SetMoreLogLabels(); }
    else           { r->GetYaxis()->SetRangeUser(cat.rLo, cat.rHi); r->GetYaxis()->SetNdivisions(505); }
    r->GetYaxis()->SetTitleSize(0.11); r->GetYaxis()->SetTitleOffset(0.42); r->GetYaxis()->SetLabelSize(0.09);
    r->GetXaxis()->SetTitle(var.xTitle.c_str());
    r->GetXaxis()->SetTitleSize(0.13); r->GetXaxis()->SetLabelSize(0.10);
    r->Draw("E");

    TLine* one = new TLine(var.lo, 1, var.hi, 1); one->SetLineStyle(2); one->Draw();

    // chi2/ndf sadece terminale yazdirilir (plota konmaz)
    double chi2 = 0; int ndf = 0;
    for (int i = 1; i <= r->GetNbinsX(); ++i)
      if (r->GetBinError(i) > 0 && hCmp->GetBinContent(i) > 0) {
        double d = (r->GetBinContent(i) - 1) / r->GetBinError(i); chi2 += d * d; ++ndf;
      }
    printf(" chi2/ndf = %.1f/%d = %.2f\n", chi2, ndf, ndf ? chi2 / ndf : 0.);
    if (blind) {
      p1->cd();
      TLatex l2; l2.SetNDC(); l2.SetTextSize(0.032);
      l2.DrawLatex(0.16, yTxt, Form("Data blinded: %g-%g GeV", blindLo, blindHi));
    }
  }

  std::string prefix = mcOnly ? "MCstack_" : "DataMC_stack_Z_";
  c->SaveAs((kOut + prefix + id + ".png").c_str());
  c->SaveAs((kOut + prefix + id + ".pdf").c_str());
}

// ============================================================================
void DataMC_stack_Z_multi_3() {
  gSystem->mkdir(kOut.c_str(), true);
  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.13);
  TH1::SetDefaultSumw2(true);

  std::vector<Var> vars = MakeVars();
  std::vector<Cat> cats = MakeCats();      // inclusive ILK olmali (kUseInclusiveK icin)

  for (const auto& cat : cats) {
    printf("\n########## KATEGORI: %s ##########\n", cat.tag.empty() ? "inclusive" : cat.tag.c_str());

    // 1) Z bolgesi: Data/MC stack + ratio
    for (const auto& v0 : vars) {
      Var v = Rebinned(v0, cat.rebin);
      double lo = kApplyZWindow ? kZMassLo : 0., hi = kApplyZWindow ? kZMassHi : 0.;
      // etiket uygulanan kesimden uretilir (pencere kapaliysa kesim yazilmaz)
      std::string zlab = kApplyZWindow
          ? std::string(Form("Z region: %g < m_{#mu#mu} < %g GeV", kZMassLo, kZMassHi))
          : std::string("Z region (no mass window)");
      StackPlot(v, v0.name, zlab.c_str(), lo, hi, false, cat);
    }

    // 2) 115-135 GeV: dimuon pT, dijet pT (+ VBF/ggH icin dijet mass ve dEta) -- sadece MC, data blind
    for (const auto& v0 : vars) {
      bool keep = (v0.name == "mmpt" || v0.name == "jjpt");
      if (!cat.tag.empty() && (v0.name == "jjmass" || v0.name == "jjdeta")) keep = true;
      if (!keep) continue;
      StackPlot(Rebinned(v0, cat.rebin), v0.name + "_SR115135", "115 < m_{#mu#mu} < 135 GeV",
                kSRMassLo, kSRMassHi, true, cat);
    }

    // 3) 110-150 GeV dimuon mass: data 115-135 blind (MC normalizasyonu sideband'lerden)
    Var vw = {"mmass_110_150", "mm_mass_fsr", "", 80, 110, 150, "m_{#mu#mu} [GeV]", "GeV"};
    StackPlot(Rebinned(vw, cat.rebin), vw.name, "110 < m_{#mu#mu} < 150 GeV", 110., 150., false, cat,
              kSRMassLo, kSRMassHi);
  }

  printf("\nBitti. Plotlar: %s\n", kOut.c_str());
}
