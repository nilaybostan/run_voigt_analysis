// ============================================================================
// NilayBostan/Aug 6 - 2026 - CERN
// FULL ANALYSIS MACRO — MC ONLY
// CMS Run2025 KIT
//
// Muon correctionlib
// BeamSpot muon pT
// Muon ID/ISO/Trigger SF
// PU weight
// JERC
// Jet ID
// b-tag SF
// VBF / ggH categorization
//
// COMPILE VERSION
// ============================================================================


#pragma cling add_include_path("/cvmfs/cms.cern.ch/el9_amd64_gcc12/external/py3-correctionlib/2.2.2-120738cfaaf3f7c1056fe67d97e25dac/lib/python3.9/site-packages/correctionlib/include")

#include "JetCorrections.h"
#include "correction.h"

#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TLorentzVector.h>
#include <TError.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <memory>



// ============================================================================
// GLOBAL CORRECTIONS
// ============================================================================


std::shared_ptr<correction::CorrectionSet> muonCorr;
std::shared_ptr<correction::CorrectionSet> muonSF;
std::shared_ptr<correction::CorrectionSet> btagCorr;


std::shared_ptr<const correction::Correction>
muon_m_mc;


std::shared_ptr<const correction::Correction>
muon_a_mc;


std::shared_ptr<const correction::Correction>
btagSF;


std::shared_ptr<const correction::Correction>
btagLightSF;

JetCorrections jetCorr;

// ============================================================================
// DELTA FUNCTIONS
// ============================================================================


double deltaPhi(
double phi1,
double phi2
)
{
    double dphi = phi1 - phi2;

    while(dphi > M_PI)
        dphi -= 2*M_PI;

    while(dphi <= -M_PI)
        dphi += 2*M_PI;

    return dphi;
}



double deltaR(
double eta1,
double phi1,
double eta2,
double phi2
)
{

    double deta =
    eta1-eta2;


    double dphi =
    deltaPhi(phi1,phi2);


    return sqrt(
        deta*deta +
        dphi*dphi
    );

}




// ============================================================================
// MUON PT CORRECTION
// ============================================================================


double getCorrectedMuonPt_MC(
double pt,
double eta,
double phi
)
{

    if(!muon_m_mc || !muon_a_mc)
        return pt;


    double m =
    muon_m_mc->evaluate(
        {
            eta,
            phi,
            "nominal"
        }
    );


    double a =
    muon_a_mc->evaluate(
        {
            eta,
            phi,
            "nominal"
        }
    );


    return pt*m+a;

}






// ============================================================================
// MUON ID SF
// ============================================================================


double getMuonIDSF(
double pt,
double eta
)
{

    auto corr =
    muonSF->at(
    "NUM_MediumID_DEN_TrackerMuons"
    );


   return corr->evaluate(
{
    eta,
    pt,
    "nominal"
});

}



// ============================================================================
// MUON ISO SF
// ============================================================================


double getMuonIsoSF(
double pt,
double eta
)
{

    auto corr =
    muonSF->at(
    "NUM_MediumPFIso_DEN_MediumID"
    );


   return corr->evaluate(
{
    eta,
    pt,
    "nominal"
});

}




// ============================================================================
// TRIGGER SF
// ============================================================================


double getTriggerSF(
double pt,
double eta
)
{

    auto corr =
    muonSF->at(
    "NUM_IsoMu24_DEN_CutBasedIdMedium_and_PFIsoMedium"
    );


   return corr->evaluate(
{
    eta,
    pt,
    "nominal"
});

}



// ============================================================================
// PU TABLE
// ============================================================================


struct PUWeight
{
    int nMin;
    int nMax;
    double weight;
};




std::vector<PUWeight>
readPUWeights(
const std::string& filename
)
{

    std::vector<PUWeight> table;


    std::ifstream f(filename);


    if(!f.is_open())
    {
        std::cerr
        <<"Cannot open PU file "
        <<filename
        <<std::endl;

        return table;
    }


    PUWeight x;


    while(
        f >>
        x.nMin >>
        x.nMax >>
        x.weight
    )
    {
        table.push_back(x);
    }


    return table;

}






double getPUWeight(
float nTruePU,
const std::vector<PUWeight>& table
)
{

    for(auto const& x : table)
    {

        if(
        nTruePU>=x.nMin &&
        nTruePU<x.nMax
        )
            return x.weight;

    }


    return 1.0;

}

// ============================================================================
// MAIN FUNCTION
// ============================================================================


void testmacro_FullCorrections_FullJetsMuons_weightMC_lastv_2025kit_bs
(
std::vector<std::string> inputFiles,

double xsec_pb,

double lumi_fb = 26.50,

const std::string& KITDir =
"/afs/cern.ch/user/n/nbostan/new_CMS/CMSSW_14_0_18/"
"src/RoccoR/post2022E-update/"
)

{


gErrorIgnoreLevel = kError;




// ============================================================================
// LOAD CORRECTIONS
// ============================================================================


std::cout<<"Loading muon correction"<<std::endl;

muonCorr =
correction::CorrectionSet::from_file
(
KITDir+"schemaV2_2025.json"
);

std::cout<<"Muon correction OK"<<std::endl;


muonSF =
correction::CorrectionSet::from_file
(
KITDir+"muon_Z.json"
);

std::cout<<"Muon SF OK"<<std::endl;


btagCorr =
correction::CorrectionSet::from_file
(
KITDir+"btagging.json"
);

std::cout<<"Btag OK"<<std::endl;


btagSF =
btagCorr->at(
"UParTAK4_comb"
);


btagLightSF =
btagCorr->at(
"UParTAK4_light"
);




std::cout
<<"Corrections loaded"
<<std::endl;







// ============================================================================
// PU
// ============================================================================


auto puWeights =
readPUWeights
(
KITDir+
"puWeights_2025pp_Golden_Summer24_25ns_69200ub.txt"
);






// ============================================================================
// BUILD CHAIN
// ============================================================================


TChain chain("Events");



for(auto const& f : inputFiles)
{
    chain.Add(f.c_str());
}



std::cout
<<"Entries = "
<<chain.GetEntries()
<<std::endl;






// ============================================================================
// NORMALIZATION
// ============================================================================


double lumi_pb =
lumi_fb*1000.0;



double sumGenWeight=0.0;



{

TTreeReader r(&chain);


TTreeReaderValue<float>
gw(r,"genWeight");



while(r.Next())
{

    sumGenWeight += *gw;

}

}




double normFactor =
(xsec_pb*lumi_pb)
/
sumGenWeight;




std::cout
<<"sumGenWeight = "
<<sumGenWeight
<<std::endl;



std::cout
<<"normFactor = "
<<normFactor
<<std::endl;







// ============================================================================
// TREE READER
// ============================================================================


TTreeReader reader(&chain);







// ============================================================================
// EVENT BRANCHES
// ============================================================================


// -------------------------
// Generator
// -------------------------


TTreeReaderValue<float>
genWeight(
reader,
"genWeight"
);




// -------------------------
// PU
// -------------------------


TTreeReaderValue<float>* nTruePU=nullptr;



if(chain.GetBranch("Pileup_nTrueInt"))
{

nTruePU =
new TTreeReaderValue<float>
(
reader,
"Pileup_nTrueInt"
);

}






// ============================================================================
// MUONS
// ============================================================================


TTreeReaderValue<int>
nMuon(
reader,
"nMuon"
);



TTreeReaderArray<float>
Muon_pt(
reader,
"Muon_pt"
);



TTreeReaderArray<float>
Muon_eta(
reader,
"Muon_eta"
);



TTreeReaderArray<float>
Muon_phi(
reader,
"Muon_phi"
);



TTreeReaderArray<float>
Muon_iso(
reader,
"Muon_pfRelIso04_all"
);



TTreeReaderArray<int>
Muon_charge(
reader,
"Muon_charge"
);



TTreeReaderArray<bool>
Muon_mediumID(
reader,
"Muon_mediumId"
);







// ============================================================================
// BEAM SPOT MUON
// ============================================================================


bool hasBS=false;



if(chain.GetBranch("Muon_bsConstrainedPt"))
    hasBS=true;




TTreeReaderArray<float>* Muon_bsConstrainedPt=nullptr;

TTreeReaderArray<float>* Muon_bsConstrainedChi2=nullptr;



if(hasBS)
{

Muon_bsConstrainedPt =
new TTreeReaderArray<float>
(
reader,
"Muon_bsConstrainedPt"
);



Muon_bsConstrainedChi2 =
new TTreeReaderArray<float>
(
reader,
"Muon_bsConstrainedChi2"
);

}







// ============================================================================
// TRIGGER OBJECT
// ============================================================================


TTreeReaderValue<int>
nTrigObj(
reader,
"nTrigObj"
);



TTreeReaderArray<float>
TrigObj_pt(
reader,
"TrigObj_pt"
);



TTreeReaderArray<float>
TrigObj_eta(
reader,
"TrigObj_eta"
);



TTreeReaderArray<float>
TrigObj_phi(
reader,
"TrigObj_phi"
);



TTreeReaderArray<UShort_t>
TrigObj_id(
reader,
"TrigObj_id"
);



TTreeReaderArray<ULong64_t>* TrigObj_filterBits=nullptr;



if(chain.GetBranch("TrigObj_filterBits"))
{

TrigObj_filterBits =
new TTreeReaderArray<ULong64_t>
(
reader,
"TrigObj_filterBits"
);

}





// ============================================================================
// HLT
// ============================================================================


TTreeReaderValue<bool>
HLT_IsoMu24(
reader,
"HLT_IsoMu24"
);







// ============================================================================
// JETS
// ============================================================================


TTreeReaderValue<int>
nJet(
reader,
"nJet"
);



TTreeReaderArray<float>
Jet_pt(
reader,
"Jet_pt"
);



TTreeReaderArray<float>
Jet_eta(
reader,
"Jet_eta"
);



TTreeReaderArray<float>
Jet_phi(
reader,
"Jet_phi"
);



TTreeReaderArray<float>
Jet_mass(
reader,
"Jet_mass"
);



TTreeReaderArray<float>
Jet_rawFactor(
reader,
"Jet_rawFactor"
);





TTreeReaderArray<float>* Jet_btag=nullptr;



if(chain.GetBranch("Jet_btagDeepFlavB"))
{

Jet_btag =
new TTreeReaderArray<float>
(
reader,
"Jet_btagDeepFlavB"
);

}



TTreeReaderArray<unsigned char> Jet_hadronFlavour(
    reader,
    "Jet_hadronFlavour"
);



// ============================================================================
// HISTOGRAMS
// ============================================================================


TH1F *h_mass =
new TH1F(
"h_mass",
"Dimuon mass",
100,
0,
200
);


TH1F *h_dimuonPt =
new TH1F(
"h_dimuonPt",
"Dimuon pT",
100,
0,
200
);


TH1F *h_dimuonEta =
new TH1F(
"h_dimuonEta",
"Dimuon eta",
50,
-2.5,
2.5
);


TH1F *h_leadJetPt =
new TH1F(
"h_leadJetPt",
"Leading jet pT",
100,
0,
500
);


TH1F *h_leadJetEta =
new TH1F(
"h_leadJetEta",
"Leading jet eta",
50,
-5,
5
);


TH1F *h_dijetMass =
new TH1F(
"h_dijetMass",
"mjj",
100,
0,
2000
);


TH1F *h_dijetPt =
new TH1F(
"h_dijetPt",
"dijet pT",
100,
0,
1000
);


TH1F *h_mass_VBF =
new TH1F(
"h_mass_VBF",
"VBF dimuon mass",
100,
0,
200
);


TH1F *h_mass_ggH =
new TH1F(
"h_mass_ggH",
"ggH dimuon mass",
100,
0,
200
);



for(auto h :
{
h_mass,
h_dimuonPt,
h_dimuonEta,
h_leadJetPt,
h_leadJetEta,
h_dijetMass,
h_dijetPt,
h_mass_VBF,
h_mass_ggH
})
{
    h->Sumw2();
}





// ============================================================================
// CUT FLOW
// ============================================================================


long long nTotal=0;
long long nTrigger=0;
long long nTwoMuon=0;
long long nMuonSel=0;
long long nTrigMatch=0;
long long nFinal=0;






// ============================================================================
// EVENT LOOP
// ============================================================================


while(reader.Next())
{


nTotal++;




// ==========================================================================
// TRIGGER
// ==========================================================================


if(!(*HLT_IsoMu24))
    continue;


nTrigger++;




// ==========================================================================
// EXACTLY TWO MUONS
// ==========================================================================


if(*nMuon != 2)
    continue;


nTwoMuon++;






// ==========================================================================
// EVENT WEIGHT
// ==========================================================================


double weight = 1.0;



weight *=
(*genWeight);



weight *=
normFactor;



if(nTruePU)
{

weight *=
getPUWeight(
**nTruePU,
puWeights
);

}





// ==========================================================================
// BEAM SPOT PT
// ==========================================================================


double muPt[2];



for(int i=0;i<2;i++)
{


    if(
    hasBS &&
    (*Muon_bsConstrainedChi2)[i] < 30.
    )
    {

        muPt[i] =
        (*Muon_bsConstrainedPt)[i];

    }

    else
    {

        muPt[i] =
        Muon_pt[i];

    }


}







// ==========================================================================
// MUON CORRECTION + ID + ISO
// ==========================================================================


double corrPt[2];


bool passMuon=true;



for(int i=0;i<2;i++)
{


    corrPt[i] =
    getCorrectedMuonPt_MC(
        muPt[i],
        Muon_eta[i],
        Muon_phi[i]
    );



    if(corrPt[i] < 20.)
        passMuon=false;



    if(fabs(Muon_eta[i]) > 2.4)
        passMuon=false;



    if(!Muon_mediumID[i])
        passMuon=false;



    if(Muon_iso[i] > 0.25)
        passMuon=false;


}



if(!passMuon)
    continue;


nMuonSel++;






// ==========================================================================
// TAG MUON
// ==========================================================================


int tag=-1;

double lead=-1;



for(int i=0;i<2;i++)
{

    if(
    corrPt[i] > 26 &&
    corrPt[i] > lead
    )
    {

        lead=corrPt[i];
        tag=i;

    }

}



if(tag<0)
    continue;



int probe =
1-tag;







// ==========================================================================
// TRIGGER MATCH
// ==========================================================================


bool matched=false;



for(int i=0;i<*nTrigObj;i++)
{


    if(TrigObj_id[i]!=13)
        continue;



    if(TrigObj_pt[i]<24)
        continue;




    double dr =
    deltaR(
        Muon_eta[tag],
        Muon_phi[tag],
        TrigObj_eta[i],
        TrigObj_phi[i]
    );



    if(dr<0.1)
    {

        matched=true;
        break;

    }


}



if(!matched)
    continue;


nTrigMatch++;







// ==========================================================================
// OPPOSITE SIGN
// ==========================================================================


if(
Muon_charge[tag]*
Muon_charge[probe]
>=0
)
continue;







// ==========================================================================
// MUON SCALE FACTORS
// ==========================================================================


weight *=
getMuonIDSF(
corrPt[tag],
Muon_eta[tag]
);


weight *=
getMuonIDSF(
corrPt[probe],
Muon_eta[probe]
);



weight *=
getMuonIsoSF(
corrPt[tag],
Muon_eta[tag]
);



weight *=
getMuonIsoSF(
corrPt[probe],
Muon_eta[probe]
);



weight *=
getTriggerSF(
corrPt[tag],
Muon_eta[tag]
);







// ==========================================================================
// DIMUON
// ==========================================================================


TLorentzVector mu1;
TLorentzVector mu2;



mu1.SetPtEtaPhiM(
corrPt[tag],
Muon_eta[tag],
Muon_phi[tag],
0.105
);



mu2.SetPtEtaPhiM(
corrPt[probe],
Muon_eta[probe],
Muon_phi[probe],
0.105
);



TLorentzVector dimuon =
mu1+mu2;



nFinal++;




h_mass->Fill(
dimuon.M(),
weight
);



h_dimuonPt->Fill(
dimuon.Pt(),
weight
);



h_dimuonEta->Fill(
dimuon.Eta(),
weight
);

// ==========================================================================
// JETS
// ==========================================================================


std::vector<TLorentzVector> jets;


int nMediumB = 0;



for(int j=0;j<*nJet;j++)
{


    // ----------------------------------------------------------
    // JERC
    // ----------------------------------------------------------


    double correctedPt =
    jetCorr.getCorrectedJetPt(
        Jet_pt[j],
        Jet_eta[j],
        Jet_phi[j],
        Jet_rawFactor[j],
        true
    );



    if(correctedPt < 20.)
        continue;



    if(fabs(Jet_eta[j]) > 4.7)
        continue;





    TLorentzVector jet;



   jet.SetPtEtaPhiM(
correctedPt,
Jet_eta[j],
Jet_phi[j],
Jet_mass[j]
);




    // ----------------------------------------------------------
    // Jet ID
    // ----------------------------------------------------------


    if(!jetCorr.passJetID(jet))
        continue;






    // ----------------------------------------------------------
    // Muon overlap removal
    // ----------------------------------------------------------


    if(jet.DeltaR(mu1)<0.4)
        continue;


    if(jet.DeltaR(mu2)<0.4)
        continue;






    jets.push_back(jet);






    // ----------------------------------------------------------
    // B TAGGING
    // DeepFlavour
    // Medium WP 2025
    // ----------------------------------------------------------


    if(
    Jet_btag &&
    Jet_btag->GetSize() > j &&
    Jet_hadronFlavour.GetSize() > j
)
{


        double discr =
        (*Jet_btag)[j];



        if(discr > 0.1272)
        {
            nMediumB++;
        }



    }



}








// ==========================================================================
// B-TAG SCALE FACTOR
// ==========================================================================


double btagWeight = 1.0;


if(btagSF)
{


for(int j=0;j<*nJet;j++)
{


    double correctedPt =
    jetCorr.getCorrectedJetPt(
        Jet_pt[j],
        Jet_eta[j],
        Jet_phi[j],
        Jet_rawFactor[j],
        true
    );


    if(correctedPt < 20.)
        continue;


    if(fabs(Jet_eta[j])>4.7)
        continue;



    TLorentzVector jet;


    jet.SetPtEtaPhiM(
        correctedPt,
        Jet_eta[j],
        Jet_phi[j],
        Jet_mass[j]
    );



    if(!jetCorr.passJetID(jet))
        continue;


    if(jet.DeltaR(mu1)<0.4)
        continue;


    if(jet.DeltaR(mu2)<0.4)
    continue;


// BURAYA GELİYOR


int hadronFlavour = 0;

if(Jet_hadronFlavour.GetSize() > j)
{
    hadronFlavour =
    static_cast<int>(Jet_hadronFlavour[j]);
}




int flavor = hadronFlavour;

double absEta = fabs(Jet_eta[j]);
double pt = correctedPt;


double sf = 1.0;

if(flavor==5 || flavor==4)
{

    if(absEta < 2.5)
    {

        sf =
        btagSF->evaluate({
            "central",
            "M",
            flavor,
            absEta,
            pt
        });

    }

}
else if(flavor==0)
{

    if(absEta < 2.5)
    {

        sf =
        btagLightSF->evaluate({
            "central",
            "M",
            flavor,
            absEta,
            pt
        });

    }
    else
    {
        sf = 1.0;
    }

}


btagWeight *= sf;

}

}




weight *= btagWeight;







// ==========================================================================
// SORT JETS
// ==========================================================================


std::sort(
jets.begin(),
jets.end(),
[](const TLorentzVector&a,
   const TLorentzVector&b)
{
    return a.Pt()>b.Pt();
}
);






// ==========================================================================
// LEADING JET
// ==========================================================================


if(!jets.empty())
{

h_leadJetPt->Fill(
jets[0].Pt(),
weight
);



h_leadJetEta->Fill(
jets[0].Eta(),
weight
);


}







// ==========================================================================
// VBF / ggH CATEGORY
// ==========================================================================


bool isVBF=false;
bool isGGH=false;





if(jets.size()>=2)
{


    TLorentzVector j1=jets[0];
    TLorentzVector j2=jets[1];



    double mjj =
    (j1+j2).M();



    double deta =
    fabs(
    j1.Eta()-j2.Eta()
    );




    h_dijetMass->Fill(
        mjj,
        weight
    );


    h_dijetPt->Fill(
        (j1+j2).Pt(),
        weight
    );







    if(
        mjj>400 &&
        deta>2.5 &&
        j1.Pt()>35 &&
        j2.Pt()>25 &&
        nMediumB==0
    )
    {

        isVBF=true;

    }

    else if(nMediumB==0)
    {

        isGGH=true;

    }



}

else if(jets.size()==1)
{


    if(nMediumB==0)
        isGGH=true;


}



// ==========================================================================
// CATEGORY MASS
// ==========================================================================


if(isVBF)
{

h_mass_VBF->Fill(
dimuon.M(),
weight
);

}



if(isGGH)
{

h_mass_ggH->Fill(
dimuon.M(),
weight
);

}



} // end event loop








// ==========================================================================
// OUTPUT
// ==========================================================================


std::cout
<<"=============================="
<<std::endl;


std::cout
<<"Total events : "
<<nTotal
<<std::endl;


std::cout
<<"Trigger      : "
<<nTrigger
<<std::endl;


std::cout
<<"Two muons    : "
<<nTwoMuon
<<std::endl;


std::cout
<<"Muon sel     : "
<<nMuonSel
<<std::endl;


std::cout
<<"Trig match   : "
<<nTrigMatch
<<std::endl;


std::cout
<<"Final        : "
<<nFinal
<<std::endl;


std::cout
<<"=============================="
<<std::endl;




TFile *out =
new TFile(
"FullCorrections_2025KIT_MC.root",
"RECREATE"
);



h_mass->Write();
h_dimuonPt->Write();
h_dimuonEta->Write();

h_leadJetPt->Write();
h_leadJetEta->Write();

h_dijetMass->Write();
h_dijetPt->Write();

h_mass_VBF->Write();
h_mass_ggH->Write();



out->Close();




std::cout
<<"Output written"
<<std::endl;

}
