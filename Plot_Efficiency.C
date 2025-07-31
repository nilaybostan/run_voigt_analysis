void Plot_Efficiency() {
    
    TFile *f = new TFile("eff_output.root");

    auto *g_eff_pt = (TGraphAsymmErrors*)f->Get("g_eff_pt");
    auto *g_eff_eta = (TGraphAsymmErrors*)f->Get("g_eff_eta");
    auto *g_eff_eta_barrel = (TGraphAsymmErrors*)f->Get("g_eff_eta_barrel");
    auto *g_eff_eta_endcap = (TGraphAsymmErrors*)f->Get("g_eff_eta_endcap");
    auto *h2_eff_pt_eta = (TH2F*)f->Get("h2_eff_pt_eta");

    h2_eff_pt_eta->SetStats(kFALSE);

    TCanvas *c_pt = new TCanvas("c_pt", "Muon pT  Efficiency", 800, 600);
    g_eff_pt->SetMarkerStyle(20);
    g_eff_pt->SetMarkerColor(kBlue);
    g_eff_pt->SetMinimum(0.6);  // Set y-axis min
    g_eff_pt->SetMaximum(1.05);  // Set y-axis max
    g_eff_pt->Draw("AP");
    c_pt->SaveAs("muon_pt_efficiency.pdf");

    TCanvas *c_eta = new TCanvas("c_eta", "Muon Eta Efficiency", 800, 600);
    g_eff_eta->SetMarkerStyle(20);
    g_eff_eta->SetMarkerColor(kRed);
    g_eff_eta->SetMinimum(0.6);
    g_eff_eta->SetMaximum(1.05);
    g_eff_eta->Draw("AP");
    c_eta->SaveAs("muon_eta_efficiency.pdf");

    TCanvas *c_eta_barrel = new TCanvas("c_eta_barrel", "Muon Barrel Eta Efficiency", 800, 600);
    g_eff_eta_barrel->SetMarkerStyle(20);
    g_eff_eta_barrel->SetMarkerColor(kGreen+2);
    g_eff_eta_barrel->SetMinimum(0.6);
    g_eff_eta_barrel->SetMaximum(1.05);
    g_eff_eta_barrel->Draw("AP");
    c_eta_barrel->SaveAs("muon_eta_efficiency_barrel.pdf");

    TCanvas *c_eta_endcap = new TCanvas("c_eta_endcap", "Muon Endcap Eta Efficiency", 800, 600);
    g_eff_eta_endcap->SetMarkerStyle(20);
    g_eff_eta_endcap->SetMarkerColor(kMagenta+2);
    g_eff_eta_endcap->SetMinimum(0.6);
    g_eff_eta_endcap->SetMaximum(1.05);
    g_eff_eta_endcap->Draw("AP");
    c_eta_endcap->SaveAs("muon_eta_efficiency_endcap.pdf");

    TCanvas *c2d = new TCanvas("c2d", "Muon Efficiency vs pT and Eta", 800, 600);
    gStyle->SetPalette(1);
    h2_eff_pt_eta->SetMinimum(0.6);
    h2_eff_pt_eta->SetMaximum(1.05);
    h2_eff_pt_eta->Draw("COLZ");
    c2d->SaveAs("muon_efficiency_pt_eta_2D.pdf");

    std::cout << "Efficiency plots saved as muon_pt_efficiency.pdf, muon_eta_efficiency.pdf,"
              << " muon_eta_efficiency_barrel.pdf, muon_eta_efficiency_endcap.pdf,"
              << " muon_efficiency_pt_eta_2D.pdf" << std::endl;
              
}
