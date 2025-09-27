#include "RoccoR.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TChain.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"

// ---------------- Structs for SF tables ----------------
struct MuonSF {
    double ptMin, ptMax;
    double etaMin, etaMax;
    double sf;
};

struct PUWeight {
    int nMin, nMax;
    double weight;
};

// ---------------- Helper functions --------------------
double getMuonSF(double pt, double eta, const std::vector<MuonSF>& sfTable){
    for(const auto &bin : sfTable){
        if(pt>=bin.ptMin && pt<bin.ptMax && eta>=bin.etaMin && eta<bin.etaMax)
            return bin.sf;
    }
    return 1.0;
}

double getPUWeight(int nTruePU, const std::vector<PUWeight>& puTable){
    for(const auto &bin : puTable){
        if(nTruePU>=bin.nMin && nTruePU<bin.nMax)
            return bin.weight;
    }
    return 1.0;
}

// ---------------- File readers ------------------------
std::vector<MuonSF> readMuonSF(const std::string& filename){
    std::vector<MuonSF> table;
    std::ifstream infile(filename);
    if(!infile.is_open()){ std::cerr<<"Cannot open "<<filename<<std::endl; return table; }
    std::string line;
    while(std::getline(infile,line)){
        if(line.empty() || line[0]=='#') continue;
        std::stringstream ss(line);
        MuonSF bin;
        ss >> bin.ptMin >> bin.ptMax >> bin.etaMin >> bin.etaMax >> bin.sf;
        table.push_back(bin);
    }
    infile.close();
    return table;
}

std::vector<PUWeight> readPUWeights(const std::string& filename){
    std::vector<PUWeight> table;
    std::ifstream infile(filename);
    if(!infile.is_open()){ std::cerr<<"Cannot open "<<filename<<std::endl; return table; }
    std::string line;
    while(std::getline(infile,line)){
        if(line.empty() || line[0]=='#') continue;
        std::stringstream ss(line);
        PUWeight bin;
        ss >> bin.nMin >> bin.nMax >> bin.weight;
        table.push_back(bin);
    }
    infile.close();
    return table;
}

// ---------------- CMS-style comparison plotting -----------------
void DrawComparisonWithRatio(TH1D* hUnweighted, TH1D* hWeighted, TH1D* hRatio, const std::string& name, const std::string& xTitle) {
    TCanvas* c = new TCanvas(("c_"+name).c_str(), name.c_str(), 800, 800);

    // Upper pad
    TPad* pad1 = new TPad("pad1","pad1",0,0.3,1,1);
    pad1->SetBottomMargin(0.02);
    pad1->Draw();
    pad1->cd();

    hUnweighted->SetLineColor(kRed);
    hUnweighted->SetMarkerColor(kRed);
    hUnweighted->SetMarkerStyle(20);
    hWeighted->SetLineColor(kBlue);
    hWeighted->SetMarkerColor(kBlue);
    hWeighted->SetMarkerStyle(21);

    hUnweighted->SetTitle(("Comparison of "+name).c_str());
    hUnweighted->GetXaxis()->SetTitle("");  // X-axis hidden in upper pad
    hUnweighted->GetYaxis()->SetTitle("Events");
    hUnweighted->Draw("E");
    hWeighted->Draw("E SAME");

    TLegend* legend = new TLegend(0.6,0.7,0.88,0.88);
    legend->AddEntry(hUnweighted,"Uncorr","lep");
    legend->AddEntry(hWeighted,"Corr","lep");
    legend->Draw();

    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(42);
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.15,0.92,"CMS Preliminary");

    // Lower pad
    c->cd();
    TPad* pad2 = new TPad("pad2","pad2",0,0.05,1,0.3);
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.35);
    pad2->Draw();
    pad2->cd();

    hRatio->SetLineColor(kBlack);
    hRatio->SetMarkerStyle(20);
    hRatio->SetTitle("");
    hRatio->GetXaxis()->SetTitle(xTitle.c_str());
    hRatio->GetXaxis()->SetTitleSize(0.12);
    hRatio->GetXaxis()->SetTitleOffset(1.0);
    hRatio->GetXaxis()->SetLabelSize(0.10);
    hRatio->GetYaxis()->SetTitle("Uncorr / Corr");
    hRatio->GetYaxis()->SetNdivisions(505);
    hRatio->GetYaxis()->SetTitleSize(0.12);
    hRatio->GetYaxis()->SetTitleOffset(0.5);
    hRatio->GetYaxis()->SetLabelSize(0.10);
    hRatio->Draw("E");

    c->SaveAs((name+".pdf").c_str());
    delete c;
}

// ---------------- Main Macro --------------------------
void testmacro_FullCorrections_NoJets_2(
    std::vector<std::string> inputFiles, 
    bool isData = true,
    const std::string &puFileName = "/afs/cern.ch/user/n/nbostan/new_CMS/CMSSW_14_0_18/src/RoccoR/data/puWeights.txt",
    const std::string &muPOGDir = "/afs/cern.ch/user/n/nbostan/new_CMS/CMSSW_14_0_18/src/RoccoR/data/POG/MUO/"
)
{
    bool isMC = !isData; // internally consistent

    // --------------------- Rochester ---------------------
    RoccoR rc;
    rc.init("/afs/cern.ch/user/n/nbostan/new_CMS/CMSSW_14_0_18/src/RoccoR/data/RoccoR2022.txt");

    // --------------------- MuonSF tables -----------------
    std::vector<MuonSF> muonJPsi       = readMuonSF(muPOGDir+"muon_JPsi.txt");
    std::vector<MuonSF> muonZ          = readMuonSF(muPOGDir+"muon_Z.txt");
    std::vector<MuonSF> muonHighPt     = readMuonSF(muPOGDir+"muon_HighPt.txt");
    std::vector<MuonSF> muonScaleSmear = readMuonSF(muPOGDir+"muon_scalesmearing.txt");

    // --------------------- PU -----------------------------
    std::vector<PUWeight> puWeights;
    if(isMC) puWeights = readPUWeights(puFileName);

    // --------------------- Histograms --------------------
    TH1D* hMassCorrected = new TH1D("hMassCorrected","Dimuon Mass FullCorrected; M_{#mu#mu} [GeV]; Events",80,70,110);
    TH1D* hJetPt = new TH1D("hJetPt","Jet Pt (Raw, no corrections); p_{T} [GeV]; Jets",100,0,1000);
    TH2D* hEffPtEta = new TH2D("hEffPtEta","Probe Efficiency Map; p_{T} [GeV]; #eta",50,0,200,48,-2.4,2.4);
    TH1D* hZWindow = new TH1D("hZWindow","Z mass window (81–101 GeV); Events",1,0,1);

    TH1D* hMassUnweighted = new TH1D("hMassUnweighted","Dimuon Mass (Unweighted); M_{#mu#mu} [GeV]; Events",80,70,110);
    TH1D* hMassWeighted   = new TH1D("hMassWeighted","Dimuon Mass (Weighted); M_{#mu#mu} [GeV]; Events",80,70,110);

    TH1D* hProbePtUnweighted = new TH1D("hProbePtUnweighted","Probe p_{T} (Unweighted); p_{T} [GeV]; Events",50,0,200);
    TH1D* hProbePtWeighted   = new TH1D("hProbePtWeighted","Probe p_{T} (Weighted); p_{T} [GeV]; Events",50,0,200);

    TH1D* hProbeEtaUnweighted = new TH1D("hProbeEtaUnweighted","Probe #eta (Unweighted); #eta; Events",48,-2.4,2.4);
    TH1D* hProbeEtaWeighted   = new TH1D("hProbeEtaWeighted","Probe #eta (Weighted); #eta; Events",48,-2.4,2.4);

    TH1D* hMassRatio = new TH1D("hMassRatio","Mass ratio (Uncorr / Corr); M_{#mu#mu} [GeV]; Ratio",80,70,110);
    TH1D* hPtRatio   = new TH1D("hPtRatio","pT ratio (Uncorr / Corr); p_{T} [GeV]; Ratio",50,0,200);
    TH1D* hEtaRatio  = new TH1D("hEtaRatio","#eta ratio (Uncorr / Corr); #eta; Ratio",48,-2.4,2.4);

    // --------------------- TChain setup ------------------
    TChain chain("Events");
    for(const auto &f : inputFiles) chain.Add(f.c_str());

    // --------------------- Branches ----------------------
    Int_t nMuon, nJet;
    Float_t Muon_pt[100], Muon_eta[100], Muon_phi[100], Muon_pfRelIso04_all[100];
    Int_t Muon_charge[100];
    Bool_t Muon_tightId[100], Muon_isGlobal[100], Muon_isTracker[100], Muon_isStandalone[100];
    UChar_t Muon_nStations[100];
    Float_t Muon_genPt[100];
    Bool_t HLT_IsoMu24;
    Float_t Pileup_nTrueInt=0; // <-- fixed type

    chain.SetBranchAddress("nMuon",&nMuon);
    chain.SetBranchAddress("Muon_pt",Muon_pt);
    chain.SetBranchAddress("Muon_eta",Muon_eta);
    chain.SetBranchAddress("Muon_phi",Muon_phi);
    chain.SetBranchAddress("Muon_charge",Muon_charge);
    chain.SetBranchAddress("Muon_tightId",Muon_tightId);
    chain.SetBranchAddress("Muon_isGlobal",Muon_isGlobal);
    chain.SetBranchAddress("Muon_isTracker",Muon_isTracker);
    chain.SetBranchAddress("Muon_isStandalone",Muon_isStandalone);
    chain.SetBranchAddress("Muon_nStations",Muon_nStations);
    chain.SetBranchAddress("Muon_pfRelIso04_all",Muon_pfRelIso04_all);
    if(isMC && chain.GetBranch("Muon_genPt")) chain.SetBranchAddress("Muon_genPt",Muon_genPt);
    chain.SetBranchAddress("HLT_IsoMu24",&HLT_IsoMu24);
    if(isMC && chain.GetBranch("Pileup_nTrueInt")) chain.SetBranchAddress("Pileup_nTrueInt",&Pileup_nTrueInt);

    // Jets
    Float_t Jet_pt[200], Jet_eta[200], Jet_phi[200], Jet_mass[200], Jet_rawFactor[200];
    chain.SetBranchAddress("nJet",&nJet);
    chain.SetBranchAddress("Jet_pt",Jet_pt);
    chain.SetBranchAddress("Jet_eta",Jet_eta);
    chain.SetBranchAddress("Jet_phi",Jet_phi);
    chain.SetBranchAddress("Jet_mass",Jet_mass);
    chain.SetBranchAddress("Jet_rawFactor",Jet_rawFactor);

    // --------------------- Event loop -------------------
    Long64_t nEntries = chain.GetEntries();
    for(Long64_t i=0;i<nEntries;i++){
        chain.GetEntry(i);
        if(!HLT_IsoMu24 || nMuon<2) continue;

        double puWeight=1.0;
        if(isMC) puWeight=getPUWeight(Pileup_nTrueInt,puWeights);

        // ---------- Tag selection -----------------
        int tagIdx=-1; float maxPt=-1;
        for(int m=0;m<nMuon;m++){
            if(!Muon_tightId[m]) continue;
            if(Muon_pt[m]<29 || fabs(Muon_eta[m])>=2.4 || Muon_pfRelIso04_all[m]>=0.15) continue;
            if(!Muon_isGlobal[m] || !Muon_isTracker[m]) continue;
            if(Muon_pt[m]>maxPt){ maxPt=Muon_pt[m]; tagIdx=m; }
        }
        if(tagIdx<0) continue;

        // ---------- Probe loop -----------------
        for(int p=0;p<nMuon;p++){
            if(p==tagIdx) continue;
            if(Muon_charge[tagIdx]*Muon_charge[p]>=0) continue;
            if(!Muon_isTracker[p] || !Muon_isStandalone[p] || Muon_pt[p]<20 || fabs(Muon_eta[p])>=2.4 || Muon_nStations[p]<=1) continue;

            TLorentzVector tag, probe;
            tag.SetPtEtaPhiM(Muon_pt[tagIdx], Muon_eta[tagIdx], Muon_phi[tagIdx], 0.105);
            probe.SetPtEtaPhiM(Muon_pt[p], Muon_eta[p], Muon_phi[p], 0.105);

            // Rochester SF
            double sf_tag=1.0, sf_probe=1.0;
            if(isMC){
                if(Muon_genPt[tagIdx]>0)
                    sf_tag = rc.kSpreadMC(Muon_charge[tagIdx], Muon_pt[tagIdx], Muon_eta[tagIdx], Muon_phi[tagIdx], Muon_genPt[tagIdx],0,0);
                else sf_tag = rc.kScaleMC(Muon_charge[tagIdx], Muon_pt[tagIdx], Muon_eta[tagIdx], Muon_phi[tagIdx],0,0);

                if(Muon_genPt[p]>0)
                    sf_probe = rc.kSpreadMC(Muon_charge[p], Muon_pt[p], Muon_eta[p], Muon_phi[p], Muon_genPt[p],0,0);
                else sf_probe = rc.kScaleMC(Muon_charge[p], Muon_pt[p], Muon_eta[p], Muon_phi[p],0,0);
            } else {
                sf_tag = rc.kScaleDT(Muon_charge[tagIdx], Muon_pt[tagIdx], Muon_eta[tagIdx], Muon_phi[tagIdx],0,0);
                sf_probe = rc.kScaleDT(Muon_charge[p], Muon_pt[p], Muon_eta[p], Muon_phi[p],0,0);
            }

            // MuonPOG SFs
            double sf_tag_pog= getMuonSF(Muon_pt[tagIdx], Muon_eta[tagIdx], muonZ);  
            double sf_probe_pog= getMuonSF(Muon_pt[p], Muon_eta[p], muonZ);

            // ScaleSmearing SF
            double sf_tag_smear  = getMuonSF(Muon_pt[tagIdx], Muon_eta[tagIdx], muonScaleSmear);
            double sf_probe_smear= getMuonSF(Muon_pt[p], Muon_eta[p], muonScaleSmear);

            // --- Fill unweighted histograms ---
            float massUnweighted = (tag+probe).M();
            hMassUnweighted->Fill(massUnweighted);
            hProbePtUnweighted->Fill(probe.Pt());
            hProbeEtaUnweighted->Fill(probe.Eta());

            // Apply corrections
            tag.SetPtEtaPhiM(Muon_pt[tagIdx]*sf_tag*sf_tag_pog*sf_tag_smear, Muon_eta[tagIdx], Muon_phi[tagIdx], 0.105);
            probe.SetPtEtaPhiM(Muon_pt[p]*sf_probe*sf_probe_pog*sf_probe_smear, Muon_eta[p], Muon_phi[p], 0.105);

            double totalWeight = puWeight * sf_tag * sf_probe * sf_tag_pog * sf_probe_pog * sf_tag_smear * sf_probe_smear;

            // --- Fill weighted histograms ---
            float massCorr=(tag+probe).M();
            hMassCorrected->Fill(massCorr, totalWeight);
            hMassWeighted->Fill(massCorr, totalWeight);
            hProbePtWeighted->Fill(probe.Pt(), totalWeight);
            hProbeEtaWeighted->Fill(probe.Eta(), totalWeight);

            if(massCorr>81 && massCorr<101) hZWindow->Fill(0.5, totalWeight);
            hEffPtEta->Fill(probe.Pt(), probe.Eta(), totalWeight);
        }

        // ---------------- Jets -----------------
        for(int j=0;j<nJet;j++){
            double rawPt = Jet_pt[j]*(1.0/Jet_rawFactor[j]);
            hJetPt->Fill(rawPt, puWeight);
        }
    }

    // ---------------- Compute ratios -----------------
    hMassRatio->Divide(hMassUnweighted, hMassWeighted, 1.0, 1.0, "B");
    hPtRatio->Divide(hProbePtUnweighted, hProbePtWeighted, 1.0, 1.0, "B");
    hEtaRatio->Divide(hProbeEtaUnweighted, hProbeEtaWeighted, 1.0, 1.0, "B");

    // ---------------- Save histograms -----------------
    std::string outFileName = isMC ? "output_FullCorrections_MC.root" : "output_FullCorrections_Data.root";
    TFile outFile(outFileName.c_str(),"RECREATE");
    hMassCorrected->Write();
    hEffPtEta->Write();
    hJetPt->Write();
    hZWindow->Write();

    hMassUnweighted->Write();
    hMassWeighted->Write();
    hProbePtUnweighted->Write();
    hProbePtWeighted->Write();
    hProbeEtaUnweighted->Write();
    hProbeEtaWeighted->Write();
    hMassRatio->Write();
    hPtRatio->Write();
    hEtaRatio->Write();
    outFile.Close();
    std::cout<<"Saved histograms to "<<outFileName<<std::endl;

    // ---------------- Draw CMS-style comparison plots -----------------
    DrawComparisonWithRatio(hMassUnweighted, hMassWeighted, hMassRatio, "DimuonMass", "M_{#mu#mu} [GeV]");
    DrawComparisonWithRatio(hProbePtUnweighted, hProbePtWeighted, hPtRatio, "ProbePt", "p_{T} [GeV]");
    DrawComparisonWithRatio(hProbeEtaUnweighted, hProbeEtaWeighted, hEtaRatio, "ProbeEta", "#eta");
}

