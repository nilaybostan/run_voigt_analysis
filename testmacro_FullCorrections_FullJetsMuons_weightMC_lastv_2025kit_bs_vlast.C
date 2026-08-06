// ============================================================================
// Nilay Bostan
// CERN - August 2026
//
// FULL ANALYSIS MACRO — MC ONLY
// CMS Hmumu Run2025 MCs
//
// KIT correctionlib muon correction
// BeamSpot constrained muons
// PU weight
// Muon SF
// JERC MC
// Jet ID
// Jet veto map
// VBF / ggH categorization
//
// CLEAN MC VERSION
// ============================================================================

#pragma cling add_include_path("/cvmfs/cms.cern.ch/el9_amd64_gcc12/external/py3-correctionlib/2.2.2-120738cfaaf3f7c1056fe67d97e25dac/lib/python3.9/site-packages/correctionlib/include")

#pragma cling add_library_path("/cvmfs/cms.cern.ch/el9_amd64_gcc12/external/py3-correctionlib/2.2.2-120738cfaaf3f7c1056fe67d97e25dac/lib/python3.9/site-packages/correctionlib/lib")

#pragma cling load("correctionlib")

#include "correction.h"
#include "JetCorrections.h"


#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>

#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

#include <TLorentzVector.h>

#include <TRandom.h>
#include <TError.h>


#include <iostream>
#include <fstream>

#include <vector>
#include <string>

#include <cmath>
#include <algorithm>

#include <memory>




// ============================================================================
// GLOBAL CORRECTION OBJECTS
// ============================================================================


std::shared_ptr<correction::CorrectionSet>
muonCorr;



std::shared_ptr<const correction::Correction>
muon_m_mc;



std::shared_ptr<const correction::Correction>
muon_a_mc;



JetCorrections jetCorr;

// ============================================================================
// B TAG CORRECTIONS
// ============================================================================

std::shared_ptr<correction::CorrectionSet>
btagCorr;


std::shared_ptr<const correction::Correction>
btagSF;


std::shared_ptr<const correction::Correction>
btagLightSF;

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
60,
-3,
3
);



TH1F *h_jetPtCorr =
new TH1F(
"h_jetPtCorr",
"Corrected jet pT",
100,
0,
500
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
60,
-5,
5
);



TH1F *h_dijetMass =
new TH1F(
"h_dijetMass",
"Dijet mass",
100,
0,
2000
);



TH1F *h_dijetPt =
new TH1F(
"h_dijetPt",
"Dijet pT",
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





// ============================================================================
// DELTA PHI
// ============================================================================


double deltaPhi(
double phi1,
double phi2
)
{

    double dphi =
    phi1 - phi2;


    while(dphi > M_PI)
        dphi -= 2*M_PI;


    while(dphi <= -M_PI)
        dphi += 2*M_PI;


    return dphi;

}




// ============================================================================
// DELTA R
// ============================================================================


double deltaR(
double eta1,
double phi1,
double eta2,
double phi2
)
{

    double deta =
    eta1 - eta2;


    double dphi =
    deltaPhi(
        phi1,
        phi2
    );


    return sqrt(
        deta*deta +
        dphi*dphi
    );

}




// ============================================================================
// KIT MUON CORRECTION MC
//
// pT corrected = pT * m_mc + a_mc
//
// ============================================================================


double getCorrectedMuonPt_MC(

double pt,

double eta,

double phi

)
{


    if(
        !muon_m_mc ||
        !muon_a_mc
    )
    {
        return pt;
    }



    double m =
    muon_m_mc->evaluate(
    {
        eta,
        phi,
        "nom"
    });



    double a =
    muon_a_mc->evaluate(
    {
        eta,
        phi,
        "nom"
    });



    return pt*m+a;


}




// ============================================================================
// MC EVENT WEIGHT
// ============================================================================


double getMCWeight(
double genWeight
)
{

    if(genWeight>=0)
        return 1.0;

    else
        return -1.0;

}

// ============================================================================
// MAIN FUNCTION
// ============================================================================


void testmacro_FullCorrections_FullJetsMuons_weightMC_lastv_2025kit_bs_vlast
(
std::vector<std::string> inputFiles,

double xsec_pb,

double lumi_fb,

const std::string& KITDir
)

{


gErrorIgnoreLevel = kError;



// ============================================================================
// LOAD KIT MUON CORRECTIONS
// ============================================================================


std::cout
<<"Loading KIT muon corrections..."
<<std::endl;



muonCorr =

correction::CorrectionSet::from_file

(

KITDir+

"schemaV2_2025.json"

);



muon_m_mc =

muonCorr->at(

"m_mc"

);



muon_a_mc =

muonCorr->at(

"a_mc"

);



std::cout
<<"KIT MC muon correction loaded"
<<std::endl;


// ============================================================================
// LOAD BTAG CORRECTIONS
// ============================================================================


btagCorr =

correction::CorrectionSet::from_file

(
    KITDir+
    "btagging.json"
);



btagSF =

btagCorr->at(
    "UParTAK4_comb"
);



btagLightSF =

btagCorr->at(
    "UParTAK4_light"
);



std::cout
<<"Btag correction loaded"
<<std::endl;


// ============================================================================
// LOAD JET CORRECTIONS
// ============================================================================


jetCorr =

JetCorrections

(

KITDir+"jet_jerc.txt",

KITDir+"jetid.txt",

KITDir+"jetvetomaps.txt"

);



std::cout
<<"Jet corrections loaded"
<<std::endl;






// ============================================================================
// HISTOGRAM SUMW2
// ============================================================================


for(auto h:
{

h_mass,

h_dimuonPt,

h_dimuonEta,

h_jetPtCorr,

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
// BUILD TCHAIN
// ============================================================================


TChain chain(

"Events"

);



for(auto const& file:
inputFiles)

{

    chain.Add(

        file.c_str()

    );

}



std::cout
<<"Total entries = "
<<chain.GetEntries()
<<std::endl;



// ============================================================================
// MC NORMALIZATION
// ============================================================================

double lumi_pb = lumi_fb * 1000.0;


double sumGenWeight = 0.0;


TTreeReader sumReader(&chain);


TTreeReaderValue<float> sumGenWeightReader
(
    sumReader,
    "genWeight"
);



while(sumReader.Next())
{

    sumGenWeight += *sumGenWeightReader;

}



std::cout
<<"Sum genWeight = "
<<sumGenWeight
<<std::endl;

if(fabs(sumGenWeight)<1e-6)
{
    std::cerr
    <<"ERROR: Invalid sumGenWeight"
    <<std::endl;
    return;
}

double normFactor =
(
    xsec_pb *
    lumi_pb
)
/
sumGenWeight;



std::cout
<<"MC normalization factor = "
<<normFactor
<<std::endl;

// ============================================================================
// TREE READER
// ============================================================================


TTreeReader reader(

&chain

);





// ============================================================================
// EVENT INFORMATION
// ============================================================================


TTreeReaderValue<UInt_t>

run

(

reader,

"run"

);



TTreeReaderValue<UInt_t>

luminosityBlock

(

reader,

"luminosityBlock"

);





// ============================================================================
// MC WEIGHT
// ============================================================================


TTreeReaderValue<float>

genWeight

(

reader,

"genWeight"

);





// ============================================================================
// TRIGGER
// ============================================================================


TTreeReaderValue<bool>
HLT_IsoMu24(
    reader,
    "HLT_IsoMu24"
);




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
Muon_pfRelIso04_all(
    reader,
    "Muon_pfRelIso04_all"
);



TTreeReaderArray<int>
Muon_charge(
    reader,
    "Muon_charge"
);



TTreeReaderArray<bool>
Muon_mediumId(
    reader,
    "Muon_mediumId"
);





// ============================================================================
// BEAMSPOT CONSTRAINED MUONS
// ============================================================================


bool hasBS=false;



if(
    chain.GetBranch(
        "Muon_bsConstrainedPt"
    )
)
{
    hasBS=true;
}




TTreeReaderArray<float>*
Muon_bsConstrainedPt=nullptr;



TTreeReaderArray<float>*
Muon_bsConstrainedChi2=nullptr;




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



std::cout
<<"BeamSpot branch = "
<<hasBS
<<std::endl;







// ============================================================================
// TRIGGER OBJECTS
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



TTreeReaderArray<ULong64_t>* TrigObj_filterBits = nullptr;



if(
    chain.GetBranch(
        "TrigObj_filterBits"
    )
)
{

TrigObj_filterBits =
    new TTreeReaderArray<ULong64_t>(
        reader,
        "TrigObj_filterBits"
    );

}






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

// ============================================================================
// B TAG BRANCHES
// ============================================================================


TTreeReaderArray<float>* Jet_btagDeepFlavB=nullptr;

TTreeReaderArray<unsigned char>* Jet_hadronFlavour=nullptr;



if(chain.GetBranch("Jet_btagDeepFlavB"))
{

    Jet_btagDeepFlavB =
    new TTreeReaderArray<float>(
        reader,
        "Jet_btagDeepFlavB"
    );

}




if(chain.GetBranch("Jet_hadronFlavour"))
{

    Jet_hadronFlavour =
    new TTreeReaderArray<unsigned char>
    (
        reader,
        "Jet_hadronFlavour"
    );

}

// ============================================================================
// EVENT LOOP
// ============================================================================

long long nTotal   = 0;
long long nTwoMuon = 0;
long long nMuonSel = 0;
long long nFinal   = 0;

while(reader.Next())
{


    nTotal++;


    // ========================================================================
    // MC EVENT WEIGHT
    // ========================================================================


    double weight = 1.0;



    if(genWeight.GetSetupStatus()==0)
    {

        weight =
        (*genWeight)
        *
        normFactor;

    }





    // ========================================================================
    // EXACTLY TWO MUONS
    // ========================================================================


    if(*nMuon != 2)
        continue;



    nTwoMuon++;

// ============================================================================
// BEAMSPOT MUON PT
// ============================================================================


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






// ============================================================================
// KIT MUON CORRECTION
// MC
//
// pT_corrected = pT_BS*m_mc + a_mc
// ============================================================================


double corrPt[2];

bool passMuon=true;


for(int i=0;i<2;i++)
{

    corrPt[i] =
    getCorrectedMuonPt_MC
    (
        muPt[i],
        Muon_eta[i],
        Muon_phi[i]
    );


    // eta

    if(
        fabs(Muon_eta[i]) > 2.4
    )
    {
        passMuon=false;
    }


    // Medium ID

    if(
        !Muon_mediumId[i]
    )
    {
        passMuon=false;
    }


    // Isolation

    if(
        Muon_pfRelIso04_all[i] > 0.25
    )
    {
        passMuon=false;
    }

}

double leadPt =
std::max(corrPt[0], corrPt[1]);


double subleadPt =
std::min(corrPt[0], corrPt[1]);


if(leadPt < 26.)
{
    passMuon=false;
}


if(subleadPt < 20.)
{
    passMuon=false;
}
if(!passMuon)
    continue;



nMuonSel++;







// ============================================================================
// OPPOSITE SIGN
// ============================================================================


if(
    Muon_charge[0] *
    Muon_charge[1]
    >=0
)
{
    continue;
}






// ============================================================================
// TRIGGER MATCHING
// ============================================================================


bool triggerMatched=false;




for(int i=0;i<*nTrigObj;i++)
{


    if(
        TrigObj_id[i] != 13
    )
        continue;



    if(
        TrigObj_pt[i] < 24.
    )
        continue;




    if(TrigObj_filterBits)
    {

        if(
        (((*TrigObj_filterBits)[i]
        &
        (1ULL<<3))
        ==0)
        )
        {
            continue;
        }

    }




    double dr0 =
    deltaR(

        Muon_eta[0],
        Muon_phi[0],

        TrigObj_eta[i],
        TrigObj_phi[i]

    );



    double dr1 =
    deltaR(

        Muon_eta[1],
        Muon_phi[1],

        TrigObj_eta[i],
        TrigObj_phi[i]

    );




    if(
        dr0 < 0.1 ||
        dr1 < 0.1
    )
    {

        triggerMatched=true;

        break;

    }


}




if(!triggerMatched)
    continue;







// ============================================================================
// BUILD MUON FOUR VECTORS
// ============================================================================


TLorentzVector mu1;



TLorentzVector mu2;




mu1.SetPtEtaPhiM(

    corrPt[0],

    Muon_eta[0],

    Muon_phi[0],

    0.105

);



mu2.SetPtEtaPhiM(

    corrPt[1],

    Muon_eta[1],

    Muon_phi[1],

    0.105

);






TLorentzVector dimuon =
mu1 + mu2;







// ============================================================================
// DIMUON HISTOGRAMS
// ============================================================================


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

// ============================================================================
// JETS
// ============================================================================



std::vector<TLorentzVector> jets;


int nMediumB = 0;

for(int j=0; j<*nJet; j++)
{


    // ========================================================================
    // JERC CORRECTION
    //
    // MC MODE
    // ========================================================================


    double correctedPt =
    jetCorr.getCorrectedJetPt(

        Jet_pt[j],

        Jet_eta[j],

        Jet_phi[j],

        Jet_rawFactor[j],

        true       // MC

    );




    if(correctedPt < 20.)
        continue;



    if(
        fabs(Jet_eta[j]) > 4.7
    )
        continue;





    TLorentzVector jet;



    jet.SetPtEtaPhiM(

        correctedPt,

        Jet_eta[j],

        Jet_phi[j],

        Jet_mass[j]

    );






    // ========================================================================
    // JET ID
    // ========================================================================


    if(
        !jetCorr.passJetID(jet)
    )
        continue;





    // ========================================================================
    // MC:
    //
    // NO Jet veto map
    //
    // Jet veto map only DATA
    //
    // ========================================================================






    // ========================================================================
    // MUON-JET OVERLAP REMOVAL
    // ========================================================================


    if(
        jet.DeltaR(mu1)<0.4
    )
        continue;




    if(
        jet.DeltaR(mu2)<0.4
    )
        continue;





    // ========================================================================
    // STORE JET
    // ========================================================================


    jets.push_back(jet);


// ========================================================================
// B TAGGING
// DeepFlavour Medium WP 2025
// ========================================================================



double discr = 0.0;


if(
    Jet_btagDeepFlavB &&
    Jet_btagDeepFlavB->GetSize() > j
)
{

    discr =
    (*Jet_btagDeepFlavB)[j];

}



if(discr > 0.1272)
{
    nMediumB++;
}



h_jetPtCorr->Fill(

    correctedPt,

    weight

);


}






// ============================================================================
// B TAG SCALE FACTOR
// ============================================================================


double btagWeight = 1.0;



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


    if(fabs(Jet_eta[j]) > 4.7)
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




    int hadronFlavour = 0;


if(
    Jet_hadronFlavour &&
    Jet_hadronFlavour->GetSize() > j
)
{

    hadronFlavour =
    static_cast<int>(
        (*Jet_hadronFlavour)[j]
    );

}


    double absEta =
    fabs(Jet_eta[j]);



    double pt =
    correctedPt;



    double sf=1.0;


if(hadronFlavour==5 || hadronFlavour==4)
{

    if(absEta<2.5)
    {

        sf =
        btagSF->evaluate({

            "central",
            "M",
            hadronFlavour,
            absEta,
            pt

        });

    }

}


else if(hadronFlavour==0)
{

    if(absEta<2.5)
    {

        sf =
        btagLightSF->evaluate({

            "central",
            "M",
            hadronFlavour,
            absEta,
            pt

        });

    }

}


btagWeight *= sf;



}



weight *= btagWeight;

// ============================================================================
// SORT JETS BY PT
// ============================================================================


sort(

    jets.begin(),

    jets.end(),

    [](const TLorentzVector& a,
       const TLorentzVector& b)

    {

        return a.Pt() > b.Pt();

    }

);









// ============================================================================
// LEADING JET
// ============================================================================


if(
    !jets.empty()
)
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









// ============================================================================
// VBF / ggH CATEGORY
// ============================================================================


bool isVBF=false;

bool isGGH=false;





if(
    jets.size() >= 2
)
{


    TLorentzVector j1 =
    jets[0];



    TLorentzVector j2 =
    jets[1];




    double mjj =
    (j1+j2).M();




    double deta =
    fabs(
        j1.Eta()
        -
        j2.Eta()
    );





    h_dijetMass->Fill(

        mjj,

        weight

    );




    h_dijetPt->Fill(

        (j1+j2).Pt(),

        weight

    );







    // ========================================================================
    // VBF SELECTION
    // ========================================================================


    if(

mjj > 400.

&&

deta > 2.5

&&

j1.Pt() > 35.

&&

j2.Pt() > 25.

&&

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




else if(
    jets.size()==1
)
{

if(nMediumB==0)
    isGGH=true;

}








// ============================================================================
// CATEGORY MASS
// ============================================================================


if(isVBF)
{

    nFinal++;



    h_mass_VBF->Fill(

        dimuon.M(),

        weight

    );


}





if(isGGH)
{

    nFinal++;



    h_mass_ggH->Fill(

        dimuon.M(),

        weight

    );


}

// ============================================================================
// END OF EVENT LOOP
// ============================================================================


} // while(reader.Next())




// ============================================================================
// OUTPUT ROOT FILE
// ============================================================================


TFile *out =
new TFile(

"/eos/user/n/nbostan/2025_Samples/"
"FullCorrections_2025KIT_MC_bs.root",

"RECREATE"

);





if(
out->IsZombie()
)
{

    cerr
    <<"ERROR: Cannot create output ROOT file"
    <<endl;


    return;

}









// ============================================================================
// WRITE HISTOGRAMS
// ============================================================================


h_mass->Write();


h_dimuonPt->Write();


h_dimuonEta->Write();



h_jetPtCorr->Write();


h_leadJetPt->Write();


h_leadJetEta->Write();



h_dijetMass->Write();


h_dijetPt->Write();



h_mass_VBF->Write();


h_mass_ggH->Write();





out->Write();

out->Close();


// ==========================================================================
// CLEAN POINTERS
// ==========================================================================


if(Jet_btagDeepFlavB)
    delete Jet_btagDeepFlavB;


if(Jet_hadronFlavour)
    delete Jet_hadronFlavour;


if(TrigObj_filterBits)
    delete TrigObj_filterBits;


if(Muon_bsConstrainedPt)
    delete Muon_bsConstrainedPt;


if(Muon_bsConstrainedChi2)
    delete Muon_bsConstrainedChi2;



std::cout
<<"Output written"
<<std::endl;



std::cout
<<"\n====================================\n"
<<" MC ANALYSIS FINISHED\n"
<<" KIT muon correction applied\n"
<<" BeamSpot muons used\n"
<<" JERC MC applied\n"
<<" Jet ID applied\n"
<<" VBF/ggH categorization done\n"
<<" Output:\n"
<<" FullCorrections_2025KIT_MC_bs.root\n"
<<"====================================\n";

}

