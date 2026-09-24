void pT_withRatio_binWidth_SR(){
   
    
   
    gStyle->SetOptStat(0);
    gStyle->SetPadLeftMargin(0.14);
    gStyle->SetPadRightMargin(0.05);
    
    //double scaleFactor = 50 / 923.6;
    //double scaleFactor1 = 50 / 923.6;

    // --- Input files ---
  TFile *fDY = TFile::Open("FullCorrections_2025KIT_DY_bs_PU_MuonSF_SR_2_2024.root");
  TFile *fEWK = TFile::Open("FullCorrections_2025KIT_ewk_bs_PU_MuonSF_SR_2_2024.root");
  TFile *fTT = TFile::Open("FullCorrections_2025KIT_ttbar_bs_PU_MuonSF_SR_2_2024.root");
  TFile *fVV = TFile::Open("FullCorrections_2025KIT_vv_bs_PU_MuonSF_SR_2_2024.root");
  TFile *fVBF = TFile::Open("FullCorrections_2025KIT_vbf_bs_PU_MuonSF_SR_2_2024.root");
  TFile *fGGH = TFile::Open("FullCorrections_2025KIT_ggH_bs_PU_MuonSF_SR_2_2024.root");
  TFile *fData = TFile::Open("output_histos_DATA_2025_KIT_bs_golden_2_SR_2_2024_large.root");


    if(!fDY || !fEWK || !fTT || !fVV || !fVBF || !fGGH || !fData){
        std::cerr << "ERROR: missing input files\n";
        return;
    }

    auto getHist = [&](TFile* f, const char* histName = "h_dimuonPt_115_135"){
        TH1F* h = (TH1F*)f->Get(histName);
        return h ? (TH1F*)h->Clone() : nullptr;
    };

    TH1F *hDY   = getHist(fDY);
    TH1F *hEWK  = getHist(fEWK);
    TH1F *hTT   = getHist(fTT);
    TH1F *hVV   = getHist(fVV);
    // VBF and ggH pT now come from the sideband (115 < m_mumu < 135) histograms
    TH1F *hVBF  = getHist(fVBF, "h_dimuonPt_115_135_VBF");
    TH1F *hGGH  = getHist(fGGH, "h_dimuonPt_115_135_ggH");
    TH1F *hData = getHist(fData);

    if(!hDY || !hEWK || !hTT || !hVV || !hVBF || !hGGH || !hData){
        std::cerr << "ERROR: missing input histograms\n";
        return;
    }
    
    //hTT->Scale(scaleFactor);
    

    // --- Bin-width normalization ---
    auto normBW = [](TH1F* h){
        for(int i=1;i<=h->GetNbinsX();++i){
            double w = h->GetBinWidth(i);
            h->SetBinContent(i, h->GetBinContent(i)/w);
            h->SetBinError(i,   h->GetBinError(i)/w);
        }
    };

    for(auto h : {hDY,hEWK,hTT,hVV,hVBF,hGGH,hData}) normBW(h);
    
    

    // --- Normalize MC to Data ---
    TH1F *hMC = (TH1F*)hDY->Clone("hMC");
    hMC->Add(hEWK); hMC->Add(hTT); hMC->Add(hVV);
    hMC->Add(hVBF); hMC->Add(hGGH);

    double mcScale = hData->Integral() / hMC->Integral();
    for(auto h : {hDY,hEWK,hTT,hVV,hVBF,hGGH}) h->Scale(mcScale);

    // ===================== Styles =====================

    // --- Filled backgrounds ---
    hDY ->SetFillColor(kBlue-7);
    hEWK->SetFillColor(kTeal-7);
    hTT ->SetFillColor(kRed-7);
    hVV ->SetFillColor(kOrange-3);

    // Solid fill + black borders
    for(auto h : {hDY,hEWK,hTT,hVV}){
        h->SetFillStyle(1001);
        h->SetLineColor(kBlack);
        h->SetLineWidth(1);
    }

    // --- Signal lines ---
    hVBF->SetLineColor(kGreen+2);
    hVBF->SetLineWidth(3);
    hVBF->SetFillStyle(0);

    hGGH->SetLineColor(kMagenta+2);
    hGGH->SetLineWidth(3);
    hGGH->SetFillStyle(0);

    // --- Data ---
    hData->SetMarkerStyle(20);
    hData->SetMarkerSize(1.1);
    hData->SetLineWidth(2);

    // ===================== Canvas & pads =====================
    TCanvas *c = new TCanvas("c","Dimuon mass",900,800);

    TPad *pad1 = new TPad("pad1","",0,0.30,1,1);
    TPad *pad2 = new TPad("pad2","",0,0,1,0.30);

    pad1->SetBottomMargin(0.02);
    pad1->SetTopMargin(0.08);
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.35);

    pad1->Draw();
    pad2->Draw();

    // ===================== Upper pad =====================
    pad1->cd();
    pad1->SetLogy();

    THStack *hs = new THStack("hs","");

    // sensible order: small → large
    hs->Add(hVV);
    hs->Add(hTT);
    hs->Add(hEWK);
    hs->Add(hDY);

    hs->Draw("HIST");
    hs->GetXaxis()->SetRangeUser(0,200);
    hs->GetYaxis()->SetTitle("Events / GeV");
    hs->SetMinimum(1e-4);

    hVBF->Draw("HIST SAME");
    hGGH->Draw("HIST SAME");
    hData->Draw("P SAME");

    // --- Legend ---
    TLegend *leg = new TLegend(0.62,0.63,0.90,0.87);
    leg->AddEntry(hData,"Data","lep");
    leg->AddEntry(hDY,"DY","f");
    leg->AddEntry(hEWK,"EWK","f");
    leg->AddEntry(hTT,"t#bar{t}","f");
    leg->AddEntry(hVV,"VV","f");
    leg->AddEntry(hVBF,"VBF","l");
    leg->AddEntry(hGGH,"ggH","l");
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->Draw();

    // --- CMS label ---
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.05);
    latex.DrawLatex(0.15,0.92,"#bf{CMS} #it{Preliminary} 13.6 TeV, 2025, Signal Region");

    // ===================== Ratio pad =====================
    pad2->cd();

    TH1F *hRatio = (TH1F*)hData->Clone("hRatio");
    hRatio->Divide(hMC);

    hRatio->SetTitle("");
    hRatio->GetYaxis()->SetTitle("Data / MC");
    hRatio->GetYaxis()->SetRangeUser(0.5,1.5);
    hRatio->GetYaxis()->SetTitleSize(0.10);
    hRatio->GetYaxis()->SetLabelSize(0.10);
    hRatio->GetYaxis()->SetTitleOffset(0.45);

    hRatio->GetXaxis()->SetTitle("p_{T} [GeV]");
    hRatio->GetXaxis()->SetTitleSize(0.12);
    hRatio->GetXaxis()->SetLabelSize(0.10);
    hRatio->GetXaxis()->SetRangeUser(0,200);

    hRatio->Draw("P");

    TLine *l = new TLine(0,1,200,1);
    l->SetLineStyle(2);
    l->SetLineWidth(2);
    l->Draw();

    c->SaveAs("plots/pT_withRatio_binWidth_2025.pdf");
}
