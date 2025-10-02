void ratio_histos() {
    TFile *f1 = TFile::Open("output_FullCorrections_Data_2.root");
    TFile *f2 = TFile::Open("output_FullCorrections_MC_2.root");

    TH1D *h1 = (TH1D*)f1->Get("hMassweighted"); // Data
    TH1D *h2 = (TH1D*)f2->Get("hMassweighted"); // MC
    TH1D *hZ_data = (TH1D*)f1->Get("hZWindow");   // Data Z window
    TH1D *hZ_mc   = (TH1D*)f2->Get("hZWindow");   // MC Z window

    // --- Full normalization ---
    double nDataFull = h1->Integral();
    double nMCFull   = h2->Integral();
    TH1D *h2_fullNorm = (TH1D*)h2->Clone("h2_fullNorm");
    if (nMCFull > 0) h2_fullNorm->Scale(nDataFull / nMCFull);
    TH1D *h_ratio_full = (TH1D*)h1->Clone("h_ratio_full");
    h_ratio_full->Divide(h2_fullNorm);

    // --- Z-window normalization ---
    double nDataZ = hZ_data->Integral();
    double nMCZ   = hZ_mc->Integral();
    TH1D *h2_ZNorm = (TH1D*)h2->Clone("h2_ZNorm");
    if (nMCZ > 0) h2_ZNorm->Scale(nDataZ / nMCZ);
    TH1D *h_ratio_Z = (TH1D*)h1->Clone("h_ratio_Z");
    h_ratio_Z->Divide(h2_ZNorm);

    // --- Plotting function ---
    auto makeCanvas = [&](TH1D *hData, TH1D *hMC, TH1D *hRatio,
                          const char *legMC, const char *outName) {
        hData->SetStats(0);
        hMC->SetStats(0);
        hRatio->SetStats(0);

        hData->SetLineColor(kRed);
        hData->SetLineWidth(2);
        hMC->SetLineColor(kBlue);
        hMC->SetLineWidth(2);

        TCanvas *c = new TCanvas(outName, outName, 800, 800);
        // Top pad bigger, bottom pad smaller
        TPad *pad1 = new TPad("pad1","Top Pad",0,0.25,1,1.0); // now 75% of canvas
        TPad *pad2 = new TPad("pad2","Bottom Pad",0,0.0,1,0.25); // 25% of canvas
        pad1->SetBottomMargin(0.02);
        pad2->SetTopMargin(0.02);
        pad2->SetBottomMargin(0.4);
        pad1->Draw();
        pad2->Draw();

        // --- Top pad ---
        pad1->cd();
        hData->SetTitle("Full Muon Corrections Dimuon Mass");
        hData->GetYaxis()->SetTitle("Events (normalized)");
        hData->GetXaxis()->SetLabelSize(0); // remove x labels from top pad
        hData->Draw("HIST");
        hMC->Draw("HIST SAME");

        auto legend = new TLegend(0.65,0.75,0.88,0.88);
        legend->AddEntry(hData,"Data","l");
        legend->AddEntry(hMC,legMC,"l");
        legend->Draw();

        // Add CMS Preliminary
        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.04);
        latex.DrawLatex(0.15,0.85,"CMS Preliminary");

        // --- Bottom pad ---
        pad2->cd();
        hRatio->SetTitle("");  // no title
        hRatio->GetYaxis()->SetRangeUser(0.7,1.3);
        hRatio->GetYaxis()->SetTitle("Data / MC");
        hRatio->GetXaxis()->SetTitle("Dimuon Mass [GeV]");
        hRatio->GetYaxis()->SetNdivisions(505);
        hRatio->GetYaxis()->SetTitleSize(0.09);
        hRatio->GetYaxis()->SetTitleOffset(0.5);
        hRatio->GetYaxis()->SetLabelSize(0.08);
        hRatio->GetXaxis()->SetTitleSize(0.12);
        hRatio->GetXaxis()->SetLabelSize(0.10);
        hRatio->Draw("E");

        // optional: draw line at ratio=1
        TLine *line = new TLine(hRatio->GetXaxis()->GetXmin(),1.0,
                                hRatio->GetXaxis()->GetXmax(),1.0);
        line->SetLineColor(kBlack);
        line->SetLineStyle(2);
        line->Draw("SAME");

        c->SaveAs(Form("%s.png", outName));
    };

    // --- Make both canvases ---
    makeCanvas(h1, h2_fullNorm, h_ratio_full,
               "DY MC (full normalization)", "comparison_ratio_fullNorm");

    makeCanvas(h1, h2_ZNorm, h_ratio_Z,
               "DY MC (Z-peak normalization)", "comparison_ratio_ZNorm");
}
