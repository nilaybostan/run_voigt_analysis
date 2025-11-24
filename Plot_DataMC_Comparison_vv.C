#include <TFile.h>
#include <TH1.h>
#include <TH1F.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPad.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TString.h>
#include <TLine.h>
#include <iostream>
#include <vector>
#include <algorithm>

// --- Helper function to get histogram ---
TH1F* tryGetHist(TFile* f, const std::vector<TString>& names) {
    for(const auto &n : names) {
        if(!f) continue;
        TH1* h = (TH1*)f->Get(n);
        if(h && h->InheritsFrom("TH1")) {
            TH1F *hf = dynamic_cast<TH1F*>(h);
            if(hf) return hf;
            TH1F *conv = (TH1F*)h->Clone(Form("%s_conv", n.Data()));
            conv->SetDirectory(nullptr);
            return conv;
        }
    }
    return nullptr;
}

void Plot_DataMC_Comparison_vv() {

    // --- Input files ---
    TFile *fMCDY   = new TFile("output_FullCorrections_DY_MC_weighted.root");
    TFile *fMCEW   = new TFile("output_FullCorrections_ewk_MC_weighted.root");
    TFile *fMCTT   = new TFile("output_FullCorrections_tt_MC_weighted.root");
    TFile *fMCVBF  = new TFile("output_FullCorrections_vbf_MC_weighted.root");
    TFile *fMCGGH  = new TFile("output_FullCorrections_ggH_MC_weighted.root");
    TFile *fMCVV   = new TFile("output_FullCorrections_VV_MC_weighted.root");
    TFile *fData   = new TFile("output_FullCorrections_Data.root");

    if(!fMCDY || !fMCEW || !fMCTT || !fMCVBF || !fMCGGH || !fMCVV || !fData) {
        std::cerr << "❌ One or more input files could not be opened!" << std::endl;
        return;
    }

    // --- Histogram list ---
    std::vector<std::pair<std::string,std::string>> histList = {
        {"h_mass",        "Dimuon invariant mass;M_{#mu#mu} [GeV];Events"},
        {"h_dimuonPt",    "Dimuon p_{T};p_{T}^{#mu#mu} [GeV];Events"},
        {"h_dimuonEta",   "Dimuon #eta;#eta^{#mu#mu};Events"},
        {"h_leadJetPt",   "Leading jet p_{T};p_{T}^{jet} [GeV];Events"},
        {"h_leadJetEta",  "Leading jet #eta;#eta^{jet};Events"},
        {"h_dijetPt",     "Dijet p_{T};p_{T}^{jj} [GeV];Events"},
        {"h_dijetMass",   "Dijet invariant mass;M_{jj} [GeV];Events"},
        {"h_jetPtCorr",   "Corrected jet p_{T};p_{T}^{corr} [GeV];Events"}
    };

    // --- CMS style ---
    gStyle->SetOptStat(0);
    gStyle->SetTitleFont(42,"XYZ");
    gStyle->SetLabelFont(42,"XYZ");
    gStyle->SetLegendFont(42);
    gStyle->SetPadTopMargin(0.05);
    gStyle->SetPadBottomMargin(0.12);
    gStyle->SetPadLeftMargin(0.16);
    gStyle->SetPadRightMargin(0.04);
    gStyle->SetTitleSize(0.05,"XYZ");
    gStyle->SetLabelSize(0.045,"XYZ");
    gStyle->SetTitleOffset(1.6,"Y");

    auto fixZeros = [](TH1F* h){
        if(!h) return;
        for(int i=1;i<=h->GetNbinsX();++i){
            if(h->GetBinContent(i) <= 0) h->SetBinContent(i, 1e-6);
        }
    };

    auto getFromFileWithFallback = [&](TFile* f, const std::string &base, const std::vector<std::string>& suffixes) -> TH1F* {
        std::vector<TString> candidates;
        candidates.push_back(base.c_str());
        for(const auto &s : suffixes) candidates.push_back(Form("%s%s", base.c_str(), s.c_str()));
        return tryGetHist(f, candidates);
    };

    for(auto &hp : histList) {
        const std::string name = hp.first;
        const std::string titleWithAxes = hp.second;

        TH1F *hDY   = getFromFileWithFallback(fMCDY,  name, {"","_DY"});
        TH1F *hEW   = getFromFileWithFallback(fMCEW,  name, {"","_EW","_ewk"});
        TH1F *hTT   = getFromFileWithFallback(fMCTT,  name, {"","_TT","_tt"});
        TH1F *hVBF  = getFromFileWithFallback(fMCVBF, name, { "", "_VBF"});
        TH1F *hGGH  = getFromFileWithFallback(fMCGGH, name, { "", "_ggH", "_GGH"});
        TH1F *hVV   = getFromFileWithFallback(fMCVV,  name, { "", "_VV"});
        TH1F *hData = getFromFileWithFallback(fData,  name, {"", "_DATA"});

        if(!hDY || !hEW || !hTT || !hVBF || !hGGH || !hVV || !hData) continue;

        TH1F *hDYc   = (TH1F*)hDY->Clone(Form("cDY_%s", name.c_str()));
        TH1F *hEWc   = (TH1F*)hEW->Clone(Form("cEW_%s", name.c_str()));
        TH1F *hTTc   = (TH1F*)hTT->Clone(Form("cTT_%s", name.c_str()));
        TH1F *hVBFc  = (TH1F*)hVBF->Clone(Form("cVBF_%s", name.c_str()));
        TH1F *hGGHc  = (TH1F*)hGGH->Clone(Form("cGGH_%s", name.c_str()));
        TH1F *hVVc   = (TH1F*)hVV->Clone(Form("cVV_%s", name.c_str()));
        TH1F *hDatac = (TH1F*)hData->Clone(Form("cDATA_%s", name.c_str()));

        hDYc->SetDirectory(0); hEWc->SetDirectory(0); hTTc->SetDirectory(0);
        hVBFc->SetDirectory(0); hGGHc->SetDirectory(0); hVVc->SetDirectory(0); hDatac->SetDirectory(0);

        fixZeros(hDYc); fixZeros(hEWc); fixZeros(hTTc);
        fixZeros(hVBFc); fixZeros(hGGHc); fixZeros(hVVc); fixZeros(hDatac);

        // --- CMS official colors and thick lines ---
        hDYc->SetLineColor(kBlue+1);    hDYc->SetFillColorAlpha(kBlue+1,0.35); hDYc->SetLineWidth(3);
        hEWc->SetLineColor(kOrange+1);  hEWc->SetFillColorAlpha(kOrange+1,0.35); hEWc->SetLineWidth(3);
        hTTc->SetLineColor(kRed+1);     hTTc->SetFillColorAlpha(kRed+1,0.35); hTTc->SetLineWidth(3);
        hVVc->SetLineColor(kCyan+2);    hVVc->SetFillColorAlpha(kCyan+2,0.35); hVVc->SetLineWidth(3);
        hVBFc->SetLineColor(kGreen+2);  hVBFc->SetLineWidth(3); hVBFc->SetFillStyle(0);
        hGGHc->SetLineColor(kMagenta+2); hGGHc->SetLineWidth(3); hGGHc->SetFillStyle(0);
        hDatac->SetMarkerStyle(20); hDatac->SetMarkerSize(1.2); hDatac->SetLineWidth(3);

        bool isDijet = (name=="h_dijetPt" || name=="h_dijetMass");

        // --- Main linear/log plots with ratio ---
        for(int ilog=0; ilog<2; ++ilog){
            bool doLog = (ilog==1);
            TString canvName = Form("c_%s_%s", name.c_str(), doLog ? "log" : "lin");
            TCanvas *c = new TCanvas(canvName, canvName, 900, 750);

            // Top pad
            TPad *pad1 = new TPad("pad1","",0,0.33,1,1);
            pad1->SetBottomMargin(0.02); pad1->SetTopMargin(0.06);
            pad1->SetLeftMargin(0.16); pad1->SetRightMargin(0.04); pad1->SetTicks(1,1);
            pad1->Draw(); pad1->cd();

            THStack *hs = new THStack("hs","");
            if(isDijet) { hs->Add(hVBFc); hs->Add(hGGHc); }
            else { hs->Add(hDYc); hs->Add(hEWc); hs->Add(hTTc); hs->Add(hVVc); }
            hs->Draw("HIST");

            hs->GetXaxis()->SetLabelSize(0);
            hs->GetXaxis()->SetTitleSize(0);
            hs->GetYaxis()->SetTitle("Events");
            hs->GetYaxis()->SetTitleSize(0.05);
            hs->GetYaxis()->SetLabelSize(0.045);

            if(!isDijet){ hVBFc->Draw("HIST SAME"); hGGHc->Draw("HIST SAME"); }
            hDatac->Draw("E SAME");

            if(doLog) pad1->SetLogy();

            TH1F *hTotalMC = (TH1F*)hDYc->Clone(); hTotalMC->SetDirectory(0);
            hTotalMC->Add(hEWc); hTotalMC->Add(hTTc); hTotalMC->Add(hVVc);
            if(!isDijet){ hTotalMC->Add(hVBFc); hTotalMC->Add(hGGHc); }

            double ymax = std::max(hDatac->GetMaximum(), hTotalMC->GetMaximum())*1.4;
            hs->SetMaximum(ymax); hs->SetMinimum(doLog ? 1e-3 : 0.0);

            TLegend *leg = new TLegend(0.66,0.62,0.95,0.88);
            leg->AddEntry(hDatac,"Data","lep");
            if(!isDijet){
                leg->AddEntry(hDYc,"DY","f"); leg->AddEntry(hEWc,"EWK","f");
                leg->AddEntry(hTTc,"ttbar","f"); leg->AddEntry(hVVc,"VV","f");
                leg->AddEntry(hVBFc,"VBF","l"); leg->AddEntry(hGGHc,"ggH","l");
            } else { leg->AddEntry(hVBFc,"VBF","l"); leg->AddEntry(hGGHc,"ggH","l"); }
            leg->SetBorderSize(0); leg->SetFillStyle(0); leg->Draw();

            TLatex latex; latex.SetNDC(); latex.SetTextSize(0.045);
            latex.DrawLatex(0.12,0.94,"#bf{CMS} #it{Preliminary} Run 3  #sqrt{s}=13.6 TeV");

            // Ratio pad
            c->cd();
            TPad *pad2 = new TPad("pad2","",0,0,1,0.33);
            pad2->SetTopMargin(0.06); pad2->SetBottomMargin(0.3);
            pad2->SetLeftMargin(0.16); pad2->SetRightMargin(0.04); pad2->SetTicks(1,1);
            pad2->Draw(); pad2->cd();

            TH1F *hRatio = (TH1F*)hDatac->Clone(); hRatio->Divide(hTotalMC);
            hRatio->SetMarkerStyle(20); hRatio->SetMarkerSize(1.0); hRatio->SetLineWidth(2);
            hRatio->SetTitle("");
            hRatio->GetYaxis()->SetTitle("Data / MC");
            hRatio->GetYaxis()->SetTitleSize(0.11);
            hRatio->GetYaxis()->SetLabelSize(0.10);
            hRatio->GetYaxis()->SetTitleOffset(0.5);
            hRatio->GetYaxis()->SetRangeUser(0.,3.5);
            hRatio->GetXaxis()->SetTitle(hDatac->GetXaxis()->GetTitle());
            hRatio->GetXaxis()->SetTitleSize(0.12);
            hRatio->GetXaxis()->SetLabelSize(0.10);
            hRatio->Draw("E");

            TLine *l1 = new TLine(hRatio->GetXaxis()->GetXmin(),1,hRatio->GetXaxis()->GetXmax(),1);
            l1->SetLineColor(kBlack); l1->SetLineStyle(2); l1->SetLineWidth(2); l1->Draw();

            c->SaveAs(Form("Plot_%s_%s.pdf", name.c_str(), doLog?"log":"lin"));

            delete leg; delete hs; delete pad1; delete pad2; delete hTotalMC; delete hRatio; delete l1; delete c;
        }

        // --- Cut plots section ---
        std::vector<std::pair<std::string,double>> cutInfo;
        if(name=="h_mass") cutInfo.push_back({"Mass70to115", 70.0});
        if(name=="h_leadJetPt") cutInfo.push_back({"leadJetPt25", 25.0});

        for(auto &cut : cutInfo) {
            TString tag = cut.first;
            double xmin = cut.second;
            double xmax = (name=="h_mass")?115.0:hDatac->GetXaxis()->GetXmax();

            TH1F *hDY_cut   = (TH1F*)hDYc->Clone(); hDY_cut->SetDirectory(0);
            TH1F *hEW_cut   = (TH1F*)hEWc->Clone(); hEW_cut->SetDirectory(0);
            TH1F *hTT_cut   = (TH1F*)hTTc->Clone(); hTT_cut->SetDirectory(0);
            TH1F *hVBF_cut  = (TH1F*)hVBFc->Clone(); hVBF_cut->SetDirectory(0);
            TH1F *hGGH_cut  = (TH1F*)hGGHc->Clone(); hGGH_cut->SetDirectory(0);
            TH1F *hVV_cut   = (TH1F*)hVVc->Clone(); hVV_cut->SetDirectory(0);
            TH1F *hData_cut = (TH1F*)hDatac->Clone(); hData_cut->SetDirectory(0);

            int nbins=hDatac->GetNbinsX();
            for(int ib=1; ib<=nbins; ++ib){
                double x = hDatac->GetBinCenter(ib);
                if(x<xmin || x>xmax){
                    hDY_cut->SetBinContent(ib,0); hEW_cut->SetBinContent(ib,0); hTT_cut->SetBinContent(ib,0);
                    hVBF_cut->SetBinContent(ib,0); hGGH_cut->SetBinContent(ib,0); hVV_cut->SetBinContent(ib,0);
                    hData_cut->SetBinContent(ib,0);
                }
            }

            fixZeros(hDY_cut); fixZeros(hEW_cut); fixZeros(hTT_cut);
            fixZeros(hVBF_cut); fixZeros(hGGH_cut); fixZeros(hVV_cut); fixZeros(hData_cut);

            for(int ilog=0; ilog<2; ++ilog){
                bool doLog = (ilog==1);
                TString cName = Form("c_cut_%s_%s", name.c_str(), doLog?"log":"lin");
                TCanvas *c = new TCanvas(cName, cName, 900,750);

                // Decide if Mass70to115 should have ratio pad or not
                bool addRatioPad = !(name=="h_mass" && tag=="Mass70to115");

                if(addRatioPad){
                    // --- same ratio pad as above for leadJetPt25 ---
                    TPad *pad1 = new TPad("pad1","",0,0.33,1,1);
                    pad1->SetBottomMargin(0.02); pad1->SetTopMargin(0.06);
                    pad1->SetLeftMargin(0.16); pad1->SetRightMargin(0.04); pad1->SetTicks(1,1);
                    pad1->Draw(); pad1->cd();

                    THStack *hs_cut = new THStack("hs_cut","");
                    hs_cut->Add(hDY_cut); hs_cut->Add(hEW_cut); hs_cut->Add(hTT_cut); hs_cut->Add(hVV_cut);
                    hs_cut->Draw("HIST");
                    hVBF_cut->Draw("HIST SAME"); hGGH_cut->Draw("HIST SAME");
                    hData_cut->Draw("E SAME");

                    hs_cut->GetXaxis()->SetLabelSize(0);
                    hs_cut->GetXaxis()->SetTitleSize(0);
                    hs_cut->GetYaxis()->SetTitle("Events");
                    hs_cut->GetYaxis()->SetTitleSize(0.05);
                    hs_cut->GetYaxis()->SetLabelSize(0.045);

                    if(doLog) pad1->SetLogy();

                    TH1F *hTotalMC_cut = (TH1F*)hDY_cut->Clone(); hTotalMC_cut->SetDirectory(0);
                    hTotalMC_cut->Add(hEW_cut); hTotalMC_cut->Add(hTT_cut); hTotalMC_cut->Add(hVV_cut);
                    hTotalMC_cut->Add(hVBF_cut); hTotalMC_cut->Add(hGGH_cut);
                    double ymax = std::max(hData_cut->GetMaximum(), hTotalMC_cut->GetMaximum())*1.4;
                    hs_cut->SetMaximum(ymax); hs_cut->SetMinimum(doLog ? 1e-3 : 0.0);

                    TLegend *leg = new TLegend(0.66,0.62,0.95,0.88);
                    leg->AddEntry(hData_cut,"Data","lep");
                    leg->AddEntry(hDY_cut,"DY","f"); leg->AddEntry(hEW_cut,"EWK","f");
                    leg->AddEntry(hTT_cut,"ttbar","f"); leg->AddEntry(hVV_cut,"VV","f");
                    leg->AddEntry(hVBF_cut,"VBF","l"); leg->AddEntry(hGGH_cut,"ggH","l");
                    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->Draw();

                    TLatex latex; latex.SetNDC(); latex.SetTextSize(0.035);
                    latex.DrawLatex(0.12,0.94,"#bf{CMS} #it{Preliminary} Run 3  #sqrt{s}=13.6 TeV");

                    c->cd();
                    TPad *pad2 = new TPad("pad2","",0,0,1,0.33);
                    pad2->SetTopMargin(0.06); pad2->SetBottomMargin(0.3);
                    pad2->SetLeftMargin(0.16); pad2->SetRightMargin(0.04); pad2->SetTicks(1,1);
                    pad2->Draw(); pad2->cd();

                    TH1F *hRatio = (TH1F*)hData_cut->Clone(); hRatio->Divide(hTotalMC_cut);
                    hRatio->SetMarkerStyle(20); hRatio->SetMarkerSize(1.0); hRatio->SetLineWidth(2);
                    hRatio->SetTitle("");
                    hRatio->GetYaxis()->SetTitle("Data / MC");
                    hRatio->GetYaxis()->SetTitleSize(0.11);
                    hRatio->GetYaxis()->SetLabelSize(0.10);
                    hRatio->GetYaxis()->SetTitleOffset(0.5);
                    hRatio->GetYaxis()->SetRangeUser(0.,3.5);
                    hRatio->GetXaxis()->SetTitle(hData_cut->GetXaxis()->GetTitle());
                    hRatio->GetXaxis()->SetTitleSize(0.12);
                    hRatio->GetXaxis()->SetLabelSize(0.10);
                    hRatio->Draw("E");

                    TLine *l1 = new TLine(hRatio->GetXaxis()->GetXmin(),1,hRatio->GetXaxis()->GetXmax(),1);
                    l1->SetLineColor(kBlack); l1->SetLineStyle(2); l1->SetLineWidth(2); l1->Draw();

                    c->SaveAs(Form("Plot_%s_%s_%s.pdf", name.c_str(), tag.Data(), doLog?"log":"lin"));

                    delete pad1; delete pad2; delete hs_cut; delete hTotalMC_cut; delete hRatio; delete l1; delete leg;
                } else {
                    // --- Single pad without ratio ---
                    TPad *pad = new TPad("pad","",0,0,1,1);
                    pad->SetBottomMargin(0.12); pad->SetTopMargin(0.06);
                    pad->SetLeftMargin(0.16); pad->SetRightMargin(0.04); pad->SetTicks(1,1);
                    pad->Draw(); pad->cd();

                    THStack *hs_cut = new THStack("hs_cut","");
                    hs_cut->Add(hDY_cut); hs_cut->Add(hEW_cut); hs_cut->Add(hTT_cut); hs_cut->Add(hVV_cut);
                    hs_cut->Draw("HIST");
                    hVBF_cut->Draw("HIST SAME"); hGGH_cut->Draw("HIST SAME");
                    hData_cut->Draw("E SAME");

                    hs_cut->GetXaxis()->SetTitle(hData_cut->GetXaxis()->GetTitle());
                    hs_cut->GetXaxis()->SetTitleSize(0.05);
                    hs_cut->GetXaxis()->SetLabelSize(0.045);
                    hs_cut->GetYaxis()->SetTitle("Events");
                    hs_cut->GetYaxis()->SetTitleSize(0.05);
                    hs_cut->GetYaxis()->SetLabelSize(0.045);

                    if(doLog) pad->SetLogy();

                    TH1F *hTotalMC_cut = (TH1F*)hDY_cut->Clone(); hTotalMC_cut->SetDirectory(0);
                    hTotalMC_cut->Add(hEW_cut); hTotalMC_cut->Add(hTT_cut); hTotalMC_cut->Add(hVV_cut);
                    hTotalMC_cut->Add(hVBF_cut); hTotalMC_cut->Add(hGGH_cut);
                    double ymax = std::max(hData_cut->GetMaximum(), hTotalMC_cut->GetMaximum())*1.4;
                    hs_cut->SetMaximum(ymax); hs_cut->SetMinimum(doLog ? 1e-3 : 0.0);

                    TLegend *leg = new TLegend(0.66,0.62,0.95,0.88);
                    leg->AddEntry(hData_cut,"Data","lep");
                    leg->AddEntry(hDY_cut,"DY","f"); leg->AddEntry(hEW_cut,"EWK","f");
                    leg->AddEntry(hTT_cut,"ttbar","f"); leg->AddEntry(hVV_cut,"VV","f");
                    leg->AddEntry(hVBF_cut,"VBF","l"); leg->AddEntry(hGGH_cut,"ggH","l");
                    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->Draw();

                    TLatex latex; latex.SetNDC(); latex.SetTextSize(0.035);
                    latex.DrawLatex(1.02,0.92,"#bf{CMS} #it{Preliminary} Run 3  #sqrt{s}=13.6 TeV");

                    c->SaveAs(Form("Plot_%s_%s_%s.pdf", name.c_str(), tag.Data(), doLog?"log":"lin"));

                    delete pad; delete hs_cut; delete hTotalMC_cut; delete leg;
                }

                delete c;
            }
        }
    }
}
