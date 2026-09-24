// ============================================================================
// CMS H->mumu Data/MC plots
//
// Produces:
//   1) Dimuon mass  : h_mass, h_mass_VBF, h_mass_ggH
//   2) Dimuon pT    : h_dimuonPt, h_dimuonPt_VBF, h_dimuonPt_ggH
//   3) Dimuon eta   : h_dimuonEta, h_dimuonEta_VBF, h_dimuonEta_ggH
//
// Main changes:
//   - CMS-style canvas/pad margins so ALL axes are visible
//   - CMS-style legend and labels
//   - Data points + stacked backgrounds
//   - VBF/ggH as line overlays
//   - MC statistical uncertainty band in ratio
//   - Data/MC ratio
//   - Mass blinding: 115 < m_mumu < 135 GeV
//   - Saves both PDF and PNG
//
// ============================================================================

#include <iostream>
#include <string>
#include <vector>
#include <cmath>

#include "TFile.h"
#include "TH1F.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TGraphAsymmErrors.h"
#include "TStyle.h"
#include "TAxis.h"

using namespace std;


// ============================================================================
// Global CMS-style settings
// ============================================================================

void SetCMSStyle()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    // Font
    gStyle->SetTextFont(42);
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetTitleFont(42, "XYZ");

    // General axis sizes
    gStyle->SetLabelSize(0.045, "XYZ");
    gStyle->SetTitleSize(0.050, "XYZ");
    gStyle->SetTitleOffset(1.10, "X");
    gStyle->SetTitleOffset(1.35, "Y");

    // Ticks
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    // Margins are set explicitly for each pad below.
}


// ============================================================================
// Histogram reader
// ============================================================================

TH1F* GetHist(TFile *f, const string &name)
{
    if (!f) {
        cerr << "ERROR: null TFile while looking for " << name << endl;
        return nullptr;
    }

    TH1F *h = dynamic_cast<TH1F*>(f->Get(name.c_str()));

    if (!h) {
        cerr << "ERROR: histogram " << name
             << " not found in " << f->GetName() << endl;
        return nullptr;
    }

    TH1F *hc = dynamic_cast<TH1F*>(h->Clone());
    hc->SetDirectory(nullptr);

    return hc;
}


// ============================================================================
// Bin-width normalization
// ============================================================================

void NormalizeByBinWidth(TH1F *h)
{
    if (!h) return;

    for (int i = 1; i <= h->GetNbinsX(); ++i) {

        const double width = h->GetBinWidth(i);

        if (width <= 0.0) continue;

        h->SetBinContent(
            i,
            h->GetBinContent(i) / width
        );

        h->SetBinError(
            i,
            h->GetBinError(i) / width
        );
    }
}


// ============================================================================
// Common histogram styling
// ============================================================================

void StyleHistograms(
    TH1F *hDY,
    TH1F *hEWK,
    TH1F *hTOP,
    TH1F *hVV,
    TH1F *hVBF,
    TH1F *hGGH,
    TH1F *hData
)
{
    // ------------------------------------------------------------------------
    // Background colors
    // ------------------------------------------------------------------------

    // Reference-like CMS colors:
    // DY  = gray
    // EWK = orange
    // TOP = dark red
    // VV  = blue

    hDY->SetFillColor(kGray + 1);
    hEWK->SetFillColor(kOrange + 1);
    hTOP->SetFillColor(kRed + 1);
    hVV->SetFillColor(kBlue + 1);

    for (TH1F *h : {hDY, hEWK, hTOP, hVV}) {

        h->SetFillStyle(1001);

        h->SetLineColor(kBlack);
        h->SetLineWidth(1);

        h->SetMarkerStyle(0);
    }

    // ------------------------------------------------------------------------
    // Signal lines
    // ------------------------------------------------------------------------

    hGGH->SetLineColor(kMagenta + 1);
    hGGH->SetLineWidth(3);
    hGGH->SetLineStyle(1);
    hGGH->SetFillStyle(0);

    hVBF->SetLineColor(kBrown + 2);
    hVBF->SetLineWidth(3);
    hVBF->SetLineStyle(1);
    hVBF->SetFillStyle(0);

    // ------------------------------------------------------------------------
    // Data
    // ------------------------------------------------------------------------

    hData->SetMarkerStyle(20);
    hData->SetMarkerSize(0.85);
    hData->SetMarkerColor(kBlack);

    hData->SetLineColor(kBlack);
    hData->SetLineWidth(1);
}


// ============================================================================
// Build total MC
// ============================================================================

TH1F* BuildTotalMC(
    TH1F *hDY,
    TH1F *hEWK,
    TH1F *hTOP,
    TH1F *hVV,
    TH1F *hVBF,
    TH1F *hGGH,
    const string &name
)
{
    TH1F *hMC = dynamic_cast<TH1F*>(hDY->Clone(name.c_str()));
    hMC->SetDirectory(nullptr);

    hMC->Reset();
    hMC->Sumw2();

    hMC->Add(hDY);
    hMC->Add(hEWK);
    hMC->Add(hTOP);
    hMC->Add(hVV);
    hMC->Add(hVBF);
    hMC->Add(hGGH);

    return hMC;
}


// ============================================================================
// MC uncertainty band for ratio
// ============================================================================

TGraphAsymmErrors* MakeMCUncertaintyBand(TH1F *hMC)
{
    if (!hMC) return nullptr;

    const int n = hMC->GetNbinsX();

    TGraphAsymmErrors *g = new TGraphAsymmErrors(n);

    for (int i = 1; i <= n; ++i) {

        const double x  = hMC->GetBinCenter(i);
        const double dx = hMC->GetBinWidth(i) / 2.0;

        const double mc  = hMC->GetBinContent(i);
        const double err = hMC->GetBinError(i);

        double relErr = 0.0;

        if (mc > 0.0)
            relErr = err / mc;

        g->SetPoint(i - 1, x, 1.0);

        g->SetPointError(
            i - 1,
            dx,
            dx,
            relErr,
            relErr
        );
    }

    g->SetFillColor(kGray + 1);
    g->SetFillStyle(1001);

    g->SetLineColor(kGray + 1);
    g->SetLineWidth(0);

    return g;
}


// ============================================================================
// CMS label
// ============================================================================

void DrawCMSLabel(
    const string &lumiText,
    const string &yearText
)
{
    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(42);

    // CMS
    latex.SetTextSize(0.055);
    latex.SetTextFont(62);
    latex.DrawLatex(
        0.17,
        0.96,
        "CMS Preliminary"
    );

    // Preliminary
    ///latex.SetTextFont(52);
    //latex.SetTextSize(0.047);
    //latex.DrawLatex(
      //  0.235,
        //0.96,
        //"Preliminary"
    //);

    // Energy + luminosity
    latex.SetTextFont(42);
    latex.SetTextSize(0.040);
    latex.DrawLatex(
        0.69,
        0.96,
        lumiText.c_str()
    );

    // Optional year
    // Uncomment if you want "2025" explicitly.
    //
    // latex.SetTextSize(0.038);
    // latex.DrawLatex(0.17, 0.855, yearText.c_str());
}


// ============================================================================
// Draw ratio panel
// ============================================================================

void DrawRatioPanel(
    TPad *pad,
    TH1F *hData,
    TH1F *hMC,
    double xMin,
    double xMax,
    const string &xTitle,
    bool massPlot = false
)
{
    pad->cd();

    // ------------------------------------------------------------------------
    // Ratio
    // ------------------------------------------------------------------------

    TH1F *hRatio =
        dynamic_cast<TH1F*>(hData->Clone("hRatio_tmp"));

    hRatio->SetDirectory(nullptr);

    hRatio->Divide(hMC);

    // Mass blinding
    if (massPlot) {

        for (int i = 1; i <= hRatio->GetNbinsX(); ++i) {

            const double x = hRatio->GetBinCenter(i);

            if (x > 115.0 && x < 135.0) {

                hRatio->SetBinContent(i, 0.0);
                hRatio->SetBinError(i, 0.0);
            }
        }
    }

    // ------------------------------------------------------------------------
    // Axis setup
    // ------------------------------------------------------------------------

    hRatio->SetTitle("");

    hRatio->GetXaxis()->SetRangeUser(xMin, xMax);

    hRatio->GetYaxis()->SetRangeUser(0.50, 1.50);

    hRatio->GetYaxis()->SetTitle("Data / MC");

    hRatio->GetYaxis()->SetTitleSize(0.095);
    hRatio->GetYaxis()->SetLabelSize(0.085);
    hRatio->GetYaxis()->SetTitleOffset(0.62);

    hRatio->GetXaxis()->SetTitle(xTitle.c_str());

    hRatio->GetXaxis()->SetTitleSize(0.105);
    hRatio->GetXaxis()->SetLabelSize(0.085);
    hRatio->GetXaxis()->SetTitleOffset(1.05);

    hRatio->GetXaxis()->SetNdivisions(505);
    hRatio->GetYaxis()->SetNdivisions(505);

    // Draw axes first
    hRatio->Draw("AXIS");

    // ------------------------------------------------------------------------
    // MC uncertainty band
    // ------------------------------------------------------------------------

    TGraphAsymmErrors *gMC =
        MakeMCUncertaintyBand(hMC);

    if (gMC)
        gMC->Draw("2 SAME");

    // ------------------------------------------------------------------------
    // Data ratio
    // ------------------------------------------------------------------------

    hRatio->SetMarkerStyle(20);
    hRatio->SetMarkerSize(0.75);
    hRatio->SetMarkerColor(kBlack);
    hRatio->SetLineColor(kBlack);
    hRatio->SetLineWidth(1);

    hRatio->Draw("E SAME");

    // ------------------------------------------------------------------------
    // Horizontal reference lines
    // ------------------------------------------------------------------------

    TLine *line1 =
        new TLine(xMin, 1.0, xMax, 1.0);

    line1->SetLineStyle(2);
    line1->SetLineColor(kGray + 2);
    line1->SetLineWidth(1);
    line1->Draw("SAME");

    // Reference lines at 0.6 and 1.4
    TLine *lineLow =
        new TLine(xMin, 0.6, xMax, 0.6);

    TLine *lineHigh =
        new TLine(xMin, 1.4, xMax, 1.4);

    lineLow->SetLineStyle(2);
    lineHigh->SetLineStyle(2);

    lineLow->SetLineColor(kGray + 1);
    lineHigh->SetLineColor(kGray + 1);

    lineLow->SetLineWidth(1);
    lineHigh->SetLineWidth(1);

    lineLow->Draw("SAME");
    lineHigh->Draw("SAME");

    // ------------------------------------------------------------------------
    // Mass blinding markers
    // ------------------------------------------------------------------------

    if (massPlot) {

        TLine *blindLow =
            new TLine(115.0, 0.5, 115.0, 1.5);

        TLine *blindHigh =
            new TLine(135.0, 0.5, 135.0, 1.5);

        blindLow->SetLineStyle(2);
        blindHigh->SetLineStyle(2);

        blindLow->SetLineColor(kGray + 2);
        blindHigh->SetLineColor(kGray + 2);

        blindLow->SetLineWidth(1);
        blindHigh->SetLineWidth(1);

        blindLow->Draw("SAME");
        blindHigh->Draw("SAME");
    }

    pad->RedrawAxis();
}


// ============================================================================
// Draw one variable plot
// ============================================================================

void DrawVariablePlot(
    TFile *fDY,
    TFile *fEWK,
    TFile *fTOP,
    TFile *fVV,
    TFile *fVBF,
    TFile *fGGH,
    TFile *fData,

    const string &varSuffix,
    const string &xTitle,

    double xMin,
    double xMax,

    double yMin,
    double yMax,

    const string &yearText,
    const string &lumiText,

    const string &pdfName,
    const string &pngName
)
{
    const string hInc  = "h_" + varSuffix;
    const string hVBFn = "h_" + varSuffix + "_VBF";
    const string hGGHn = "h_" + varSuffix + "_ggH";

    // ------------------------------------------------------------------------
    // Read histograms
    // ------------------------------------------------------------------------

    TH1F *hDY   = GetHist(fDY,   hInc);
    TH1F *hEWK  = GetHist(fEWK,  hInc);
    TH1F *hTOP  = GetHist(fTOP,  hInc);
    TH1F *hVV   = GetHist(fVV,   hInc);
    TH1F *hVBF  = GetHist(fVBF,  hVBFn);
    TH1F *hGGH  = GetHist(fGGH,  hGGHn);
    TH1F *hData = GetHist(fData, hInc);

    if (!hDY || !hEWK || !hTOP || !hVV ||
        !hVBF || !hGGH || !hData) {

        cerr << "ERROR: missing histogram(s) for "
             << varSuffix << endl;

        return;
    }

    // ------------------------------------------------------------------------
    // Bin width
    // ------------------------------------------------------------------------

    for (TH1F *h :
         {hDY, hEWK, hTOP, hVV, hVBF, hGGH, hData}) {

        NormalizeByBinWidth(h);
    }

    // ------------------------------------------------------------------------
    // Normalize total MC to Data
    // ------------------------------------------------------------------------

    TH1F *hMCbefore =
        BuildTotalMC(
            hDY,
            hEWK,
            hTOP,
            hVV,
            hVBF,
            hGGH,
            "hMCbefore"
        );

    const double dataIntegral =
        hData->Integral();

    const double mcIntegral =
        hMCbefore->Integral();

    if (mcIntegral <= 0.0) {

        cerr << "ERROR: MC integral <= 0 for "
             << varSuffix << endl;

        delete hMCbefore;
        return;
    }

    const double scale =
        dataIntegral / mcIntegral;

    for (TH1F *h :
         {hDY, hEWK, hTOP, hVV, hVBF, hGGH}) {

        h->Scale(scale);
    }

    delete hMCbefore;

    // ------------------------------------------------------------------------
    // Total MC after scaling
    // ------------------------------------------------------------------------

    TH1F *hMC =
        BuildTotalMC(
            hDY,
            hEWK,
            hTOP,
            hVV,
            hVBF,
            hGGH,
            "hMC_total"
        );

    // ------------------------------------------------------------------------
    // Style
    // ------------------------------------------------------------------------

    StyleHistograms(
        hDY,
        hEWK,
        hTOP,
        hVV,
        hVBF,
        hGGH,
        hData
    );

    // ------------------------------------------------------------------------
    // Canvas
    //
    // The larger left margin is the important fix for the missing
    // Y-axis labels/titles.
    // ------------------------------------------------------------------------

    const string cname = "c_" + varSuffix;

    TCanvas *c =
        new TCanvas(
            cname.c_str(),
            cname.c_str(),
            900,
            900
        );

    c->SetFillColor(kWhite);

    // ------------------------------------------------------------------------
    // Pads
    // ------------------------------------------------------------------------

    TPad *pad1 =
        new TPad(
            "pad1",
            "upper",
            0.0,
            0.30,
            1.0,
            1.0
        );

    TPad *pad2 =
        new TPad(
            "pad2",
            "ratio",
            0.0,
            0.0,
            1.0,
            0.30
        );

    // IMPORTANT:
    // Do NOT use left margin = 0.
    // This was the reason the Y-axis disappeared in the previous plot.

    pad1->SetLeftMargin(0.16);
    pad1->SetRightMargin(0.04);
    pad1->SetTopMargin(0.055);
    pad1->SetBottomMargin(0.015);

    pad2->SetLeftMargin(0.16);
    pad2->SetRightMargin(0.04);
    pad2->SetTopMargin(0.015);
    pad2->SetBottomMargin(0.30);

    pad1->SetTicks(1,1);
    pad2->SetTicks(1,1);

    pad1->Draw();
    pad2->Draw();

    // ========================================================================
    // UPPER PAD
    // ========================================================================

    pad1->cd();

    pad1->SetLogy();

    // ------------------------------------------------------------------------
    // Stack
    // ------------------------------------------------------------------------

    THStack *hs =
        new THStack(
            ("hs_" + varSuffix).c_str(),
            ""
        );

    // Bottom -> top
    hs->Add(hVV);
    hs->Add(hTOP);
    hs->Add(hEWK);
    hs->Add(hDY);

    hs->Draw("HIST");

    hs->GetXaxis()->SetRangeUser(xMin, xMax);

    hs->GetYaxis()->SetTitle(
        "Events / 1 GeV"
    );

    hs->GetYaxis()->SetTitleSize(0.050);
    hs->GetYaxis()->SetLabelSize(0.045);
    hs->GetYaxis()->SetTitleOffset(1.35);

    hs->GetYaxis()->SetNdivisions(510);

    hs->SetMinimum(yMin);
    hs->SetMaximum(yMax);

    // No X labels on upper pad
    hs->GetXaxis()->SetLabelSize(0.0);
    hs->GetXaxis()->SetTitleSize(0.0);

    // ------------------------------------------------------------------------
    // Signal
    // ------------------------------------------------------------------------

    hGGH->Draw("HIST SAME");
    hVBF->Draw("HIST SAME");

    // ------------------------------------------------------------------------
    // Data
    // ------------------------------------------------------------------------

    hData->Draw("E SAME");

    // ------------------------------------------------------------------------
    // Legend
    //
    // Two columns, similar to the reference CMS plot.
    // ------------------------------------------------------------------------

    TLegend *leg =
        new TLegend(
            0.55,
            0.58,
            0.94,
            0.91
        );

    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.040);
    leg->SetNColumns(2);
    leg->SetColumnSeparation(0.12);

    leg->AddEntry(hDY,   "DY",     "f");
    leg->AddEntry(hEWK,  "EWK",    "f");
    leg->AddEntry(hTOP,  "TOP",    "f");
    leg->AddEntry(hVV,   "VV",     "f");
    leg->AddEntry(hGGH,  "ggH",    "l");
    leg->AddEntry(hVBF,  "VBF",    "l");
    leg->AddEntry(hData, "Data",   "lep");

    leg->Draw();

    // ------------------------------------------------------------------------
    // CMS labels
    // ------------------------------------------------------------------------

    DrawCMSLabel(
        lumiText,
        yearText
    );

    pad1->RedrawAxis();

    // ========================================================================
    // RATIO PAD
    // ========================================================================

    DrawRatioPanel(
        pad2,
        hData,
        hMC,
        xMin,
        xMax,
        xTitle,
        false
    );

    // ------------------------------------------------------------------------
    // Save
    // ------------------------------------------------------------------------

    c->cd();

    c->SaveAs(pdfName.c_str());
    c->SaveAs(pngName.c_str());

    cout << "Saved: " << pdfName << endl;
    cout << "Saved: " << pngName << endl;

    // ------------------------------------------------------------------------
    // Cleanup
    // ------------------------------------------------------------------------

    delete hDY;
    delete hEWK;
    delete hTOP;
    delete hVV;
    delete hVBF;
    delete hGGH;
    delete hData;
    delete hMC;

    delete leg;
    delete hs;
    delete pad1;
    delete pad2;
    delete c;
}


// ============================================================================
// MASS PLOT
// ============================================================================

void DrawMassPlot(
    TFile *fDY,
    TFile *fEWK,
    TFile *fTOP,
    TFile *fVV,
    TFile *fVBF,
    TFile *fGGH,
    TFile *fData,

    const string &yearText,
    const string &lumiText,

    const string &pdfName,
    const string &pngName
)
{
    // ------------------------------------------------------------------------
    // Read
    // ------------------------------------------------------------------------

    TH1F *hDY   = GetHist(fDY,   "h_mass");
    TH1F *hEWK  = GetHist(fEWK,  "h_mass");
    TH1F *hTOP  = GetHist(fTOP,  "h_mass");
    TH1F *hVV   = GetHist(fVV,   "h_mass");
    TH1F *hVBF  = GetHist(fVBF,  "h_mass_VBF");
    TH1F *hGGH  = GetHist(fGGH,  "h_mass_ggH");
    TH1F *hData = GetHist(fData, "h_mass");

    if (!hDY || !hEWK || !hTOP || !hVV ||
        !hVBF || !hGGH || !hData) {

        cerr << "ERROR: missing mass histogram(s)." << endl;
        return;
    }

    // ------------------------------------------------------------------------
    // Bin-width normalization
    // ------------------------------------------------------------------------

    for (TH1F *h :
         {hDY, hEWK, hTOP, hVV, hVBF, hGGH, hData}) {

        NormalizeByBinWidth(h);
    }

    // ------------------------------------------------------------------------
    // MC normalization
    // ------------------------------------------------------------------------

    TH1F *hMCbefore =
        BuildTotalMC(
            hDY,
            hEWK,
            hTOP,
            hVV,
            hVBF,
            hGGH,
            "hMassMCbefore"
        );

    const double dataIntegral =
        hData->Integral();

    const double mcIntegral =
        hMCbefore->Integral();

    if (mcIntegral <= 0.0) {

        cerr << "ERROR: mass MC integral <= 0." << endl;

        delete hMCbefore;
        return;
    }

    const double scale =
        dataIntegral / mcIntegral;

    for (TH1F *h :
         {hDY, hEWK, hTOP, hVV, hVBF, hGGH}) {

        h->Scale(scale);
    }

    delete hMCbefore;

    TH1F *hMC =
        BuildTotalMC(
            hDY,
            hEWK,
            hTOP,
            hVV,
            hVBF,
            hGGH,
            "hMassMC"
        );

    // ------------------------------------------------------------------------
    // Blinded Data
    // ------------------------------------------------------------------------

    TH1F *hDataBlind =
        dynamic_cast<TH1F*>(hData->Clone("hDataBlind"));

    hDataBlind->SetDirectory(nullptr);

    for (int i = 1; i <= hDataBlind->GetNbinsX(); ++i) {

        const double x =
            hDataBlind->GetBinCenter(i);

        if (x > 115.0 && x < 135.0) {

            hDataBlind->SetBinContent(i, 0.0);
            hDataBlind->SetBinError(i, 0.0);
        }
    }

    // ------------------------------------------------------------------------
    // Style
    // ------------------------------------------------------------------------

    StyleHistograms(
        hDY,
        hEWK,
        hTOP,
        hVV,
        hVBF,
        hGGH,
        hDataBlind
    );

    // ------------------------------------------------------------------------
    // Canvas
    // ------------------------------------------------------------------------

    TCanvas *c =
        new TCanvas(
            "c_mass",
            "Dimuon mass",
            900,
            900
        );

    c->SetFillColor(kWhite);

    TPad *pad1 =
        new TPad(
            "pad1_mass",
            "upper",
            0.0,
            0.30,
            1.0,
            1.0
        );

    TPad *pad2 =
        new TPad(
            "pad2_mass",
            "ratio",
            0.0,
            0.0,
            1.0,
            0.30
        );

    // IMPORTANT FIX:
    // Large enough left margin for Y-axis title and labels.

    pad1->SetLeftMargin(0.16);
    pad1->SetRightMargin(0.04);
    pad1->SetTopMargin(0.055);
    pad1->SetBottomMargin(0.015);

    pad2->SetLeftMargin(0.16);
    pad2->SetRightMargin(0.04);
    pad2->SetTopMargin(0.015);
    pad2->SetBottomMargin(0.30);

    pad1->SetTicks(1,1);
    pad2->SetTicks(1,1);

    pad1->Draw();
    pad2->Draw();

    // ========================================================================
    // Upper
    // ========================================================================

    pad1->cd();

    pad1->SetLogy();

    THStack *hs =
        new THStack(
            "hs_mass",
            ""
        );

    hs->Add(hVV);
    hs->Add(hTOP);
    hs->Add(hEWK);
    hs->Add(hDY);

    hs->Draw("HIST");

    hs->GetXaxis()->SetRangeUser(70.0, 150.0);

    hs->GetYaxis()->SetTitle(
        "Events / 1 GeV"
    );

    hs->GetYaxis()->SetTitleSize(0.050);
    hs->GetYaxis()->SetLabelSize(0.045);
    hs->GetYaxis()->SetTitleOffset(1.35);

    hs->GetYaxis()->SetNdivisions(510);

    hs->SetMinimum(1e-3);
    hs->SetMaximum(1e7);

    // Hide upper X-axis labels
    hs->GetXaxis()->SetLabelSize(0.0);
    hs->GetXaxis()->SetTitleSize(0.0);

    // Signal
    hGGH->Draw("HIST SAME");
    hVBF->Draw("HIST SAME");

    // Blinded data
    hDataBlind->Draw("E SAME");

    // ------------------------------------------------------------------------
    // Legend
    // ------------------------------------------------------------------------

    TLegend *leg =
        new TLegend(
            0.55,
            0.58,
            0.94,
            0.91
        );

    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.040);
    leg->SetNColumns(2);
    leg->SetColumnSeparation(0.12);

    leg->AddEntry(hDY,       "DY",   "f");
    leg->AddEntry(hEWK,      "EWK",  "f");
    leg->AddEntry(hTOP,      "TOP",  "f");
    leg->AddEntry(hVV,       "VV",   "f");
    leg->AddEntry(hGGH,      "ggH",  "l");
    leg->AddEntry(hVBF,      "VBF",  "l");
    leg->AddEntry(hDataBlind,"Data", "lep");

    leg->Draw();

    // ------------------------------------------------------------------------
    // CMS label
    // ------------------------------------------------------------------------

    DrawCMSLabel(
        lumiText,
        yearText
    );

    pad1->RedrawAxis();

    // ========================================================================
    // Ratio
    // ========================================================================

    DrawRatioPanel(
        pad2,
        hData,
        hMC,
        70.0,
        150.0,
        "m_{#mu#mu} [GeV]",
        true
    );

    // ------------------------------------------------------------------------
    // Save
    // ------------------------------------------------------------------------

    c->cd();

    c->SaveAs(pdfName.c_str());
    c->SaveAs(pngName.c_str());

    cout << "Saved: " << pdfName << endl;
    cout << "Saved: " << pngName << endl;

    // ------------------------------------------------------------------------
    // Cleanup
    // ------------------------------------------------------------------------

    delete hDY;
    delete hEWK;
    delete hTOP;
    delete hVV;
    delete hVBF;
    delete hGGH;
    delete hData;
    delete hDataBlind;
    delete hMC;

    delete leg;
    delete hs;
    delete pad1;
    delete pad2;
    delete c;
}


// ============================================================================
// MAIN
// ============================================================================

void DimuonMass_pt_eta_CMSstyle_reference()
{
    SetCMSStyle();

    // =========================================================================
    // INPUT FILES
    // =========================================================================

    TFile *fDY =
        TFile::Open(
            "FullCorrections_2025KIT_DY_bs_PU_MuonSF_SR_2_2024.root"
        );

    TFile *fEWK =
        TFile::Open(
            "FullCorrections_2025KIT_ewk_bs_PU_MuonSF_SR_2_2024.root"
        );

    TFile *fTOP =
        TFile::Open(
            "FullCorrections_2025KIT_ttbar_bs_PU_MuonSF_SR_2_2024.root"
        );

    TFile *fVV =
        TFile::Open(
            "FullCorrections_2025KIT_vv_bs_PU_MuonSF_SR_2_2024.root"
        );

    TFile *fVBF =
        TFile::Open(
            "FullCorrections_2025KIT_vbf_bs_PU_MuonSF_SR_2_2024.root"
        );

    TFile *fGGH =
        TFile::Open(
            "FullCorrections_2025KIT_ggH_bs_PU_MuonSF_SR_2_2024.root"
        );

    TFile *fData =
        TFile::Open(
            "output_histos_DATA_2025_KIT_bs_golden_2_SR_2_2024_large.root"
        );

    // =========================================================================
    // Check files
    // =========================================================================

    if (!fDY || fDY->IsZombie() ||
        !fEWK || fEWK->IsZombie() ||
        !fTOP || fTOP->IsZombie() ||
        !fVV || fVV->IsZombie() ||
        !fVBF || fVBF->IsZombie() ||
        !fGGH || fGGH->IsZombie() ||
        !fData || fData->IsZombie()) {

        cerr << endl;
        cerr << "ERROR: one or more input ROOT files could not be opened."
             << endl;
        cerr << endl;

        return;
    }

    // =========================================================================
    // Output directory
    //
    // Make sure this directory exists:
    //
    //   mkdir -p plots
    //
    // =========================================================================

    gSystem->Exec("mkdir -p plots");

    // =========================================================================
    // Labels
    //
    // These are visual labels only.
    // =========================================================================

    const string yearText = "2025";

    // Change this if the actual luminosity for your sample is different.
    const string lumiText =
        "110.73 fb^{-1} (13.6 TeV)";

    // =========================================================================
    // 1. DIMUON MASS
    // =========================================================================

    DrawMassPlot(
        fDY,
        fEWK,
        fTOP,
        fVV,
        fVBF,
        fGGH,
        fData,

        yearText,
        lumiText,

        "plots/DimuonMass_CMSstyle_final.pdf",
        "plots/DimuonMass_CMSstyle_final.png"
    );

    // =========================================================================
    // 2. DIMUON pT
    //
    // Reference-like range: 0 - 300 GeV
    // =========================================================================

    DrawVariablePlot(
        fDY,
        fEWK,
        fTOP,
        fVV,
        fVBF,
        fGGH,
        fData,

        "dimuonPt",
        "p_{T}^{#mu#mu} [GeV]",

        0.0,
        200.0,

        1e-3,
        1e9,

        yearText,
        lumiText,

        "plots/DimuonPt_CMSstyle_final.pdf",
        "plots/DimuonPt_CMSstyle_final.png"
    );

    // =========================================================================
    // 3. DIMUON ETA
    // =========================================================================

    DrawVariablePlot(
        fDY,
        fEWK,
        fTOP,
        fVV,
        fVBF,
        fGGH,
        fData,

        "dimuonEta",
        "#eta_{#mu#mu}",

        -2.4,
        2.4,

        1e-3,
        1e6,

        yearText,
        lumiText,

        "plots/DimuonEta_CMSstyle_final.pdf",
        "plots/DimuonEta_CMSstyle_final.png"
    );

    // =========================================================================
    // Close files
    // =========================================================================

    fDY->Close();
    fEWK->Close();
    fTOP->Close();
    fVV->Close();
    fVBF->Close();
    fGGH->Close();
    fData->Close();

    cout << endl;
    cout << "==============================================" << endl;
    cout << " All plots have been produced." << endl;
    cout << " Output directory: plots/" << endl;
    cout << "==============================================" << endl;
}
