void Plot_Voigt_Fits() {
    TFile *f = new TFile("eff_output.root");
    auto *h_mass_pass = (TH1F*)f->Get("h_mass_pass");
    auto *h_mass_fail = (TH1F*)f->Get("h_mass_fail");

    // Disable stat boxes on these histograms
    h_mass_pass->SetStats(kFALSE);
    h_mass_fail->SetStats(kFALSE);

    // --- RooFit Voigtian fits (existing code) ---
    RooRealVar mass("mass", "Invariant Mass [GeV/c^{2}]", 70, 115);
    RooDataHist dataPass("dataPass", "Pass", mass, h_mass_pass);
    RooDataHist dataFail("dataFail", "Fail", mass, h_mass_fail);

    RooRealVar mean("mean", "Mean", 91.2, 88, 94);
    RooRealVar width("width", "Width", 2.5, 1.0, 5.0);
    RooRealVar sigma("sigma", "Sigma", 1.5, 0.5, 5.0);

    RooVoigtian voigtPass("voigtPass", "Voigt Pass", mass, mean, width, sigma);
    RooVoigtian voigtFail("voigtFail", "Voigt Fail", mass, mean, width, sigma);

    voigtPass.fitTo(dataPass, RooFit::Save(), RooFit::PrintLevel(-1));
    voigtFail.fitTo(dataFail, RooFit::Save(), RooFit::PrintLevel(-1));

    // Canvas for Voigt Pass fit
    TCanvas *c1 = new TCanvas("c1", "Voigt Fit Pass", 800, 600);
    RooPlot *frame1 = mass.frame();
    frame1->SetTitle("Invariant Mass Fit with Voigt Function (Pass)");
    dataPass.plotOn(frame1);
    voigtPass.plotOn(frame1);
    frame1->Draw();
    c1->SaveAs("voigt_fit_pass_effcut.pdf");

    // Canvas for Voigt Fail fit
    TCanvas *c2 = new TCanvas("c2", "Voigt Fit Fail", 800, 600);
    RooPlot *frame2 = mass.frame();
    frame2->SetTitle("Invariant Mass Fit with Voigt Function (Fail)");
    dataFail.plotOn(frame2);
    voigtFail.plotOn(frame2);
    frame2->Draw();
    c2->SaveAs("voigt_fit_fail_effcut.pdf");

    // --- Add Gaussian fits ---
    TF1 *gausPass = new TF1("gausPass", "gaus", 85, 96);
    h_mass_pass->Fit(gausPass, "RQ"); // Quiet fit, reduced printout
    TCanvas *c3 = new TCanvas("c3", "Gaussian Fit Pass", 800, 600);
    h_mass_pass->SetStats(kFALSE);  // Also disable here before drawing
    h_mass_pass->Draw();
    gausPass->SetLineColor(kRed);
    gausPass->Draw("SAME");
    c3->SetTitle("Invariant Mass Gaussian Fit (Pass)");
    c3->SaveAs("gaussian_fit_pass_effcut.pdf");

    TF1 *gausFail = new TF1("gausFail", "gaus", 85, 96);
    h_mass_fail->Fit(gausFail, "RQ");
    TCanvas *c4 = new TCanvas("c4", "Gaussian Fit Fail", 800, 600);
    h_mass_fail->SetStats(kFALSE); 
    h_mass_fail->Draw();
    gausFail->SetLineColor(kRed);
    gausFail->Draw("SAME");
    c4->SetTitle("Invariant Mass Gaussian Fit (Fail)");
    c4->SaveAs("gaussian_fit_fail_effcut.pdf");

    std::cout << "Fits saved as voigt_fit_pass_effcut.pdf, voigt_fit_fail_effcut.pdf,"
              << " gaussian_fit_pass_effcut.pdf, gaussian_fit_fail_effcut.pdf" << std::endl;
}
