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

void run3_TagAndProbe() {
    // --- Chain all Run 3 Muon NanoAOD files ---
    TChain *chain = new TChain("Events");
    chain->Add("root://xrootd-cms.infn.it//store/data/Run2022F/Muon/NANOAOD/22Sep2023-v2/2550000/0584d50d-f062-401b-98bf-dd9458915cd3.root");
	chain->Add("root://xrootd-cms.infn.it//store/data/Run2022F/Muon/NANOAOD/22Sep2023-v2/2550000/073d93ef-2644-4f55-8e09-597964cf9b4d.root");
	chain->Add("root://xrootd-cms.infn.it//store/data/Run2022F/Muon/NANOAOD/22Sep2023-v2/2550000/1026994c-892e-4f29-8a42-4d8452370772.root");
	chain->Add("root://xrootd-cms.infn.it//store/data/Run2022F/Muon/NANOAOD/22Sep2023-v2/2550000/187620f5-a4ae-467a-a56f-4569e7b32801.root");
	

    // --- Reader setup ---
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
    TTreeReaderArray<bool> Muon_mediumId(reader, "Muon_mediumId");  // <-- added
    TTreeReaderArray<float> Muon_pfRelIso04_all(reader, "Muon_pfRelIso04_all");
    TTreeReaderArray<bool> Muon_isStandalone(reader, "Muon_isStandalone");
    TTreeReaderValue<Bool_t> HLT_IsoMu24(reader, "HLT_IsoMu24");
    TTreeReaderArray<float> TrigObj_pt(reader, "TrigObj_pt");
    TTreeReaderArray<float> TrigObj_eta(reader, "TrigObj_eta");
    TTreeReaderArray<float> TrigObj_phi(reader, "TrigObj_phi");
    TTreeReaderArray<unsigned short> TrigObj_id(reader, "TrigObj_id");

    // Histograms (same as before)
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

    TH2F* h2_pt_eta_pass = new TH2F("h2_pt_eta_pass", ";p_{T} [GeV/c];#eta", 100, 20, 120, 96, -2.4, 2.4);
    TH2F* h2_pt_eta_total = new TH2F("h2_pt_eta_total", ";p_{T} [GeV/c];#eta", 100, 20, 120, 96, -2.4, 2.4);

    TH1F* h_eta_pass_barrel  = new TH1F("h_eta_pass_barrel",  "Barrel Muon #eta Efficiency;Probe #eta;Events", 48, -1.2, 1.2);
    TH1F* h_eta_fail_barrel  = new TH1F("h_eta_fail_barrel",  "Barrel Muon #eta Failures;Probe #eta;Events", 48, -1.2, 1.2);
    TH1F* h_eta_total_barrel = new TH1F("h_eta_total_barrel", "Barrel Muon #eta Total;Probe #eta;Events", 48, -1.2, 1.2);

    TH1F* h_eta_pass_endcap  = new TH1F("h_eta_pass_endcap",  "Endcap Muon #eta Efficiency;Probe #eta;Events", 24, 1.2, 2.4);
    TH1F* h_eta_fail_endcap  = new TH1F("h_eta_fail_endcap",  "Endcap Muon #eta Failures;Probe #eta;Events", 24, 1.2, 2.4);
    TH1F* h_eta_total_endcap = new TH1F("h_eta_total_endcap", "Endcap Muon #eta Total;Probe #eta;Events", 24, 1.2, 2.4);

    while (reader.Next()) {
        if (!*HLT_IsoMu24 || *nMuon < 2) continue;

        int tagIdx = -1;
        float maxTagPt = -1;
        for (size_t i = 0; i < Muon_pt.GetSize(); ++i) {

            // --- UPDATED SELECTION ---
            // Medium ID instead of Tight ID
            if (!Muon_mediumId[i]) continue;

            // Tag muon pT threshold: 26 GeV instead of 29 GeV
            if (Muon_pt[i] < 26 || std::abs(Muon_eta[i]) >= 2.4 || Muon_pfRelIso04_all[i] >= 0.15) continue;
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
            if (!Muon_isTracker[j] || !Muon_isStandalone[j] || Muon_pt[j] < 20 || std::abs(Muon_eta[j]) >= 2.4 || Muon_nStations[j] <= 1) continue;

            TLorentzVector tag, probe;
            tag.SetPtEtaPhiM(Muon_pt[tagIdx], Muon_eta[tagIdx], Muon_phi[tagIdx], 0.105);
            probe.SetPtEtaPhiM(Muon_pt[j], Muon_eta[j], Muon_phi[j], 0.105);
            float mass = (tag + probe).M();
            if (mass < 70 || mass > 115) continue;

            h_pt_total->Fill(Muon_pt[j]);
            h_eta_total->Fill(Muon_eta[j]);
            h_phi_total->Fill(Muon_phi[j]);
            h2_pt_eta_total->Fill(Muon_pt[j], Muon_eta[j]);

            float absEta = std::abs(Muon_eta[j]);
            if (absEta < 1.2) {
                h_eta_total_barrel->Fill(Muon_eta[j]);
            } else if (absEta >= 1.2 && absEta < 2.4) {
                if (Muon_eta[j] > 0) h_eta_total_endcap->Fill(Muon_eta[j]);
                else h_eta_total_endcap->Fill(-Muon_eta[j]);
            }

            // --- UPDATED PROBE CRITERIA ---
            // Medium ID used for passing probe
            if (Muon_mediumId[j]) {
                h_mass_pass->Fill(mass);
                h_pt_pass->Fill(Muon_pt[j]);
                h_eta_pass->Fill(Muon_eta[j]);
                h_phi_pass->Fill(Muon_phi[j]);
                h2_pt_eta_pass->Fill(Muon_pt[j], Muon_eta[j]);

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

                if (absEta < 1.2) {
                    h_eta_fail_barrel->Fill(Muon_eta[j]);
                } else if (absEta >= 1.2 && absEta < 2.4) {
                    if (Muon_eta[j] > 0) h_eta_fail_endcap->Fill(Muon_eta[j]);
                    else h_eta_fail_endcap->Fill(-Muon_eta[j]);
                }
            }
        }
    }

    // Efficiency plots (same)
    TGraphAsymmErrors* g_eff_pt = new TGraphAsymmErrors(h_pt_pass, h_pt_total, "cl=0.683 b(1,1)");
    g_eff_pt->SetName("g_eff_pt");
    g_eff_pt->SetTitle("Muon  p_{T}  Efficiency;Probe  p_{T} [GeV/c];Efficiency");

    TGraphAsymmErrors* g_eff_eta = new TGraphAsymmErrors(h_eta_pass, h_eta_total, "cl=0.683 b(1,1)");
    g_eff_eta->SetName("g_eff_eta");
    g_eff_eta->SetTitle("Muon  #eta  Efficiency;Probe  #eta;Efficiency");

    TGraphAsymmErrors* g_eff_eta_barrel = new TGraphAsymmErrors(h_eta_pass_barrel, h_eta_total_barrel, "cl=0.683 b(1,1)");
    g_eff_eta_barrel->SetName("g_eff_eta_barrel");
    g_eff_eta_barrel->SetTitle("Muon  #eta  Efficiency (Barrel);Probe #eta;Efficiency");

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
