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

// ---------------- Main Macro --------------------------
void testmacro_FullCorrections_NoJets(
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
    std::vector<MuonSF> muonScaleSmear = readMuonSF(muPOGDir+"muon_scalesmearing.txt"); // new

    // --------------------- PU -----------------------------
    std::vector<PUWeight> puWeights;
    if(isMC) puWeights = readPUWeights(puFileName);

    // --------------------- Histograms --------------------
    TH1D* hMassCorrected = new TH1D("hMassCorrected","Dimuon Mass FullCorrected; M_{#mu#mu} [GeV]; Events",80,70,110);
    TH1D* hJetPt = new TH1D("hJetPt","Jet Pt (Raw, no corrections); p_{T} [GeV]; Jets",100,0,1000);
    TH2D* hEffPtEta = new TH2D("hEffPtEta","Probe Efficiency Map; p_{T} [GeV]; #eta",50,0,200,48,-2.4,2.4);
    TH1D* hZWindow = new TH1D("hZWindow","Z mass window (81–101 GeV); Events",1,0,1);

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
    Int_t Pileup_nTrueInt=0;

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
            double sf_tag_pog= getMuonSF(Muon_pt[tagIdx], Muon_eta[tagIdx], muonZ);  // example: using Z table
            double sf_probe_pog= getMuonSF(Muon_pt[p], Muon_eta[p], muonZ);

            // ScaleSmearing SF
            double sf_tag_smear  = getMuonSF(Muon_pt[tagIdx], Muon_eta[tagIdx], muonScaleSmear);
            double sf_probe_smear= getMuonSF(Muon_pt[p], Muon_eta[p], muonScaleSmear);

            // Apply corrections
            tag.SetPtEtaPhiM(Muon_pt[tagIdx]*sf_tag*sf_tag_pog*sf_tag_smear, Muon_eta[tagIdx], Muon_phi[tagIdx], 0.105);
            probe.SetPtEtaPhiM(Muon_pt[p]*sf_probe*sf_probe_pog*sf_probe_smear, Muon_eta[p], Muon_phi[p], 0.105);

            double totalWeight = puWeight * sf_tag * sf_probe * sf_tag_pog * sf_probe_pog * sf_tag_smear * sf_probe_smear;
            float massCorr=(tag+probe).M();
            hMassCorrected->Fill(massCorr, totalWeight);
            if(massCorr>81 && massCorr<101) hZWindow->Fill(0.5, totalWeight);
            hEffPtEta->Fill(probe.Pt(), probe.Eta(), totalWeight);
        }

        // ---------------- Jets -----------------
        for(int j=0;j<nJet;j++){
            double rawPt = Jet_pt[j]*(1.0/Jet_rawFactor[j]);
            hJetPt->Fill(rawPt, puWeight);
        }
    }

    // ---------------- Save -----------------
    std::string outFileName = isMC ? "output_FullCorrections_MC.root" : "output_FullCorrections_Data.root";
    TFile outFile(outFileName.c_str(),"RECREATE");
    hMassCorrected->Write();
    hEffPtEta->Write();
    hJetPt->Write();
    hZWindow->Write();
    outFile.Close();
    std::cout<<"Saved histograms to "<<outFileName<<std::endl;
}
