//NilayBostan // July31/2025

#include <TChain.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphAsymmErrors.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <iostream>
#include <TF1.h> 
#include <RooRealVar.h>
#include <RooVoigtian.h>
#include <RooDataHist.h>
#include <RooPlot.h>
#include <RooFitResult.h>
#include <TStyle.h>

int main() {
    TChain *chain = new TChain("Events");
    chain->Add("root://xrootd-cms.infn.it//store/data/Run2022F/Muon/NANOAOD/*.root");
	
    // --- Tag-and-Probe Selection Code Start ---
    TTreeReader reader(chain);
    TTreeReaderValue<int> nMuon(reader, "nMuon");
    TTreeReaderArray<float> Muon_pt(reader, "Muon_pt");
    TTreeReaderArray<float> Muon_eta(reader, "Muon_eta");
    TTreeReaderArray<float> Muon_phi(reader, "Muon_phi");
    TTreeReaderArray<int> Muon_charge(reader, "Muon_charge");
    TTreeReaderArray<bool> Muon_isTracker(reader, "Muon_isTracker");
    TTreeReaderArray<bool> Muon_isGlobal(reader, "Muon_isGlobal");
    TTreeReaderArray<unsigned char> Muon_nStations(reader, "Muon_nStations");
    TTreeReaderArray<bool> Muon_tightId(reader, "Muon_tightId");
    TTreeReaderArray<float> Muon_pfRelIso04_all(reader, "Muon_pfRelIso04_all");
    TTreeReaderArray<bool> Muon_isStandalone(reader, "Muon_isStandalone");
    TTreeReaderValue<Bool_t> HLT_IsoMu24(reader, "HLT_IsoMu24");
    TTreeReaderArray<float> TrigObj_pt(reader, "TrigObj_pt");
    TTreeReaderArray<float> TrigObj_eta(reader, "TrigObj_eta");
    TTreeReaderArray<float> TrigObj_phi(reader, "TrigObj_phi");
    TTreeReaderArray<unsigned short> TrigObj_id(reader, "TrigObj_id");

    // Histograms for mass, pt, eta, phi pass/fail/total
    TH1F* h_mass_pass = new TH1F("h_mass_pass", ";Invariant Mass [GeV]", 70, 70, 115);
    TH1F* h_mass_fail = new TH1F("h_mass_fail", ";Invariant Mass [GeV]", 70, 70, 115);
    TH1F* h_pt_pass = new TH1F("h_pt_pass", ";Probe  p_{T} [GeV/c]", 120, 0, 120);
    TH1F* h_pt_fail = new TH1F("h_pt_fail", ";Probe  p_{T} [GeV/c]", 120, 0, 120);
    TH1F* h_pt_total = new TH1F("h_pt_total", ";Probe  p_{T} [GeV/c]", 120, 0, 120);
    TH1F* h_eta_pass = new TH1F("h_eta_pass", ";Probe  #eta", 96, -2.4, 2.4);
    TH1F* h_eta_fail = new TH1F("h_eta_fail", ";Probe  #eta", 96, -2.4, 2.4);
    TH1F* h_eta_total = new TH1F("h_eta_total", ";Probe  #eta", 96, -2.4, 2.4);
    TH1F* h_phi_pass = new TH1F("h_phi_pass", ";Probe  #phi", 128, -3.2, 3.2);
    TH1F* h_phi_fail = new TH1F("h_phi_fail", ";Probe  #phi", 128, -3.2, 3.2);
    TH1F* h_phi_total = new TH1F("h_phi_total", ";Probe  #phi", 128, -3.2, 3.2);

    // pT: 100 bins (1 GeV/bin), eta: 96 bins (~0.05/bin)
    TH2F* h2_pt_eta_pass = new TH2F("h2_pt_eta_pass", ";p_{T} [GeV/c];#eta", 100, 20, 120, 96, -2.4, 2.4);
    TH2F* h2_pt_eta_total = new TH2F("h2_pt_eta_total", ";p_{T} [GeV/c];#eta", 100, 20, 120, 96, -2.4, 2.4);

    // Barrel and Endcap eta histograms for pass/fail/total (binning same as h_eta histograms)
    TH1F* h_eta_pass_barrel = new TH1F("h_eta_pass_barrel", "Barrel Muon #eta Efficiency;Probe #eta;Events", 48, -1.2, 1.2);
    TH1F* h_eta_fail_barrel = new TH1F("h_eta_fail_barrel", "Barrel Muon #eta Failures;Probe #eta;Events", 48, -1.2, 1.2);
    TH1F* h_eta_total_barrel = new TH1F("h_eta_total_barrel", "Barrel Muon #eta Total;Probe #eta;Events", 28, -1.2, 1.2);

    TH1F* h_eta_pass_endcap = new TH1F("h_eta_pass_endcap", "Endcap Muon #eta Efficiency;Probe #eta;Events", 48, 1.2, 2.4);
    TH1F* h_eta_fail_endcap = new TH1F("h_eta_fail_endcap", "Endcap Muon #eta Failures;Probe #eta;Events", 48, 1.2, 2.4);
    TH1F* h_eta_total_endcap = new TH1F("h_eta_total_endcap", "Endcap Muon #eta Total;Probe #eta;Events", 48, 1.2, 2.4);

    while (reader.Next()) {
        if (!*HLT_IsoMu24 || *nMuon < 2) continue;

        int tagIdx = -1;
        float maxTagPt = -1;
        for (size_t i = 0; i < Muon_pt.GetSize(); ++i) {
            if (!Muon_tightId[i]) continue;
            if (Muon_pt[i] < 29 || std::abs(Muon_eta[i]) >= 2.4 || Muon_pfRelIso04_all[i] >= 0.15) continue;
            if (!Muon_isGlobal[i] || !Muon_isTracker[i]) continue;

            bool matched = false;
            for (size_t j = 0; j < TrigObj_pt.GetSize(); ++j) {
                if (TrigObj_id[j] != 13 || TrigObj_pt[j] < 24) continue;
                float deta = Muon_eta[i] - TrigObj_eta[j];
                float dphi = std::abs(Muon_phi[i] - TrigObj_phi[j]);
                if (dphi > M_PI) dphi = 2 * M_PI - dphi;
                float dr2 = deta*deta + dphi*dphi;
                if (dr2 < 0.01) matched = true;
            }

            if (matched && Muon_pt[i] > maxTagPt) {
                maxTagPt = Muon_pt[i];
                tagIdx = i;
            }
        }

        if (tagIdx < 0) continue;

        for (size_t j = 0; j < Muon_pt.GetSize(); ++j) {
            if ((int)j == tagIdx) continue;
            if (Muon_charge[tagIdx] * Muon_charge[j] >= 0) continue;
            if (!Muon_isTracker[j] || !Muon_isStandalone[j] ||  Muon_pt[j] < 20 || std::abs(Muon_eta[j]) >= 2.4 || Muon_nStations[j] <= 1) continue;

            TLorentzVector tag, probe;
            tag.SetPtEtaPhiM(Muon_pt[tagIdx], Muon_eta[tagIdx], Muon_phi[tagIdx], 0.105);
            probe.SetPtEtaPhiM(Muon_pt[j], Muon_eta[j], Muon_phi[j], 0.105);
            float mass = (tag + probe).M();
            if (mass < 70 || mass > 115) continue;

            h_pt_total->Fill(Muon_pt[j]);
            h_eta_total->Fill(Muon_eta[j]);
            h_phi_total->Fill(Muon_phi[j]);
            h2_pt_eta_total->Fill(Muon_pt[j], Muon_eta[j]);

            // Barrel / Endcap filling for total
            float absEta = std::abs(Muon_eta[j]);
            if (absEta < 1.2) {
                h_eta_total_barrel->Fill(Muon_eta[j]);
            } else if (absEta >= 1.2 && absEta < 2.4) {
                // Fill positive and negative separately for endcap (optional: mirror to positive eta)
                if (Muon_eta[j] > 0) h_eta_total_endcap->Fill(Muon_eta[j]);
                else h_eta_total_endcap->Fill(-Muon_eta[j]); // mirror negative eta into positive side
            }

            if (Muon_tightId[j]) {
                h_mass_pass->Fill(mass);
                h_pt_pass->Fill(Muon_pt[j]);
                h_eta_pass->Fill(Muon_eta[j]);
                h_phi_pass->Fill(Muon_phi[j]);
                h2_pt_eta_pass->Fill(Muon_pt[j], Muon_eta[j]);

                // Barrel / Endcap filling for pass
                if (absEta < 1.2) {
                    h_eta_pass_barrel->Fill(Muon_eta[j]);
                } else if (absEta >= 1.2 && absEta < 2.4) {
                    if (Muon_eta[j] > 0) h_eta_pass_endcap->Fill(Muon_eta[j]);
                    else h_eta_pass_endcap->Fill(-Muon_eta[j]);
                }

            } else {
                h_mass_fail->Fill(mass);
                h_pt_fail->Fill(Muon_pt[j]);
                h_eta_fail->Fill(Muon_eta[j]);
                h_phi_fail->Fill(Muon_phi[j]);

                // Barrel / Endcap filling for fail
                if (absEta < 1.2) {
                    h_eta_fail_barrel->Fill(Muon_eta[j]);
                } else if (absEta >= 1.2 && absEta < 2.4) {
                    if (Muon_eta[j] > 0) h_eta_fail_endcap->Fill(Muon_eta[j]);
                    else h_eta_fail_endcap->Fill(-Muon_eta[j]);
                }
            }
        }
    }

    // Efficiency histograms: passing / total with binomial errors
    TGraphAsymmErrors* g_eff_pt = new TGraphAsymmErrors(h_pt_pass, h_pt_total, "cl=0.683 b(1,1)");
    g_eff_pt->SetName("g_eff_pt");
    g_eff_pt->SetTitle("Muon  p_{T}  Efficiency;Probe  p_{T} [GeV/c];Efficiency");

    TGraphAsymmErrors* g_eff_eta = new TGraphAsymmErrors(h_eta_pass, h_eta_total, "cl=0.683 b(1,1)");
    g_eff_eta->SetName("g_eff_eta");
    g_eff_eta->SetTitle("Muon  #eta  Efficiency;Probe  #eta;Efficiency");

    // Barrel eta efficiency
    TGraphAsymmErrors* g_eff_eta_barrel = new TGraphAsymmErrors(h_eta_pass_barrel, h_eta_total_barrel, "cl=0.683 b(1,1)");
    g_eff_eta_barrel->SetName("g_eff_eta_barrel");
    g_eff_eta_barrel->SetTitle("Muon  #eta  Efficiency (Barrel);Probe #eta;Efficiency");

    // Endcap eta efficiency
    TGraphAsymmErrors* g_eff_eta_endcap = new TGraphAsymmErrors(h_eta_pass_endcap, h_eta_total_endcap, "cl=0.683 b(1,1)");
    g_eff_eta_endcap->SetName("g_eff_eta_endcap");
    g_eff_eta_endcap->SetTitle("Muon  #eta  Efficiency (Endcap);Probe  #eta;Efficiency");

    TH2F* h2_eff_pt_eta = (TH2F*)h2_pt_eta_pass->Clone("h2_eff_pt_eta");
    h2_eff_pt_eta->SetTitle("Muon Efficiency vs  p_{T}  and  #eta;Probe   p_{T} [GeV/c];Probe  #eta;Efficiency");
    h2_eff_pt_eta->Divide(h2_pt_eta_total);

    TFile fout("eff_output.root", "RECREATE");
    h_mass_pass->Write(); h_mass_fail->Write();
    h_pt_pass->Write(); h_pt_fail->Write(); h_pt_total->Write();
    h_eta_pass->Write(); h_eta_fail->Write(); h_eta_total->Write();
    h_phi_pass->Write(); h_phi_fail->Write(); h_phi_total->Write();
    h2_pt_eta_pass->Write(); h2_pt_eta_total->Write();

    h_eta_pass_barrel->Write(); h_eta_fail_barrel->Write(); h_eta_total_barrel->Write();
    h_eta_pass_endcap->Write(); h_eta_fail_endcap->Write(); h_eta_total_endcap->Write();

    g_eff_pt->Write();
    g_eff_eta->Write();
    g_eff_eta_barrel->Write();
    g_eff_eta_endcap->Write();

    h2_eff_pt_eta->Write();

    fout.Close();

    std::cout << "Histograms and efficiencies written to eff_output.root" << std::endl;
}

void Plot_Voigt_Fits() {
    TFile *f = new TFile("eff_output.root");
    auto *h_mass_pass = (TH1F*)f->Get("h_mass_pass");
    auto *h_mass_fail = (TH1F*)f->Get("h_mass_fail");

    // Disable stat boxes on these histograms
    h_mass_pass->SetStats(kFALSE);
    h_mass_fail->SetStats(kFALSE);

    // --- RooFit Voigtian fits ---
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
    h_mass_fail->SetStats(kFALSE);  // Disable stat box here too
    h_mass_fail->Draw();
    gausFail->SetLineColor(kRed);
    gausFail->Draw("SAME");
    c4->SetTitle("Invariant Mass Gaussian Fit (Fail)");
    c4->SaveAs("gaussian_fit_fail_effcut.pdf");

    std::cout << "Fits saved as voigt_fit_pass_effcut.pdf, voigt_fit_fail_effcut.pdf,"
              << " gaussian_fit_pass_effcut.pdf, gaussian_fit_fail_effcut.pdf" << std::endl;
}

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
    g_eff_pt->SetMinimum(0.5);
    g_eff_pt->SetMaximum(1.0);
    g_eff_pt->Draw("AP");
    c_pt->SaveAs("muon_pt_efficiency.pdf");

    TCanvas *c_eta = new TCanvas("c_eta", "Muon Eta Efficiency", 800, 600);
    g_eff_eta->SetMarkerStyle(20);
    g_eff_eta->SetMarkerColor(kRed);
    g_eff_eta->SetMinimum(0.5);
    g_eff_eta->SetMaximum(1.0);
    g_eff_eta->Draw("AP");
    c_eta->SaveAs("muon_eta_efficiency.pdf");

    TCanvas *c_eta_barrel = new TCanvas("c_eta_barrel", "Muon Barrel Eta Efficiency", 800, 600);
    g_eff_eta_barrel->SetMarkerStyle(20);
    g_eff_eta_barrel->SetMarkerColor(kGreen+2);
    g_eff_eta_barrel->SetMinimum(0.5);
    g_eff_eta_barrel->SetMaximum(1.0);
    g_eff_eta_barrel->Draw("AP");
    c_eta_barrel->SaveAs("muon_eta_efficiency_barrel.pdf");

    TCanvas *c_eta_endcap = new TCanvas("c_eta_endcap", "Muon Endcap Eta Efficiency", 800, 600);
    g_eff_eta_endcap->SetMarkerStyle(20);
    g_eff_eta_endcap->SetMarkerColor(kMagenta+2);
    g_eff_eta_endcap->SetMinimum(0.5);
    g_eff_eta_endcap->SetMaximum(1.0);
    g_eff_eta_endcap->Draw("AP");
    c_eta_endcap->SaveAs("muon_eta_efficiency_endcap.pdf");

    TCanvas *c2d = new TCanvas("c2d", "Muon Efficiency vs pT and Eta", 800, 600);
    gStyle->SetPalette(1);
    h2_eff_pt_eta->SetMinimum(0.5);
    h2_eff_pt_eta->SetMaximum(1.0);
    h2_eff_pt_eta->Draw("COLZ");
    c2d->SaveAs("muon_efficiency_pt_eta_2D.pdf");

    std::cout << "Efficiency plots saved as muon_pt_efficiency.pdf, muon_eta_efficiency.pdf,"
              << " muon_eta_efficiency_barrel.pdf, muon_eta_efficiency_endcap.pdf,"
              << " muon_efficiency_pt_eta_2D.pdf" << std::endl;
}
