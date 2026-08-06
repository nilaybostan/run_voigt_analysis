// ============================================================================
// Nilay Bostan
// CERN - August 2026
//
// FULL ANALYSIS MACRO — DATA ONLY
// CMS Run2025
//
// KIT correctionlib muon correction
// BeamSpot constrained muons
// Golden JSON
// JERC
// Jet ID
// Jet veto maps
// Trigger matching
// VBF / ggH categorization
//
// CLEAN DATA VERSION
// ============================================================================


#pragma cling add_include_path("/cvmfs/cms.cern.ch/el9_amd64_gcc12/external/py3-correctionlib/2.2.2-120738cfaaf3f7c1056fe67d97e25dac/lib/python3.9/site-packages/correctionlib/include")



#include "JetCorrections.h"
#include "correction.h"


#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>
#include <TLorentzVector.h>
#include <TError.h>


#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <memory>
#include <map>


#include <nlohmann/json.hpp>


using json = nlohmann::json;



// ============================================================================
// GLOBAL CORRECTIONS
// ============================================================================


std::shared_ptr<correction::CorrectionSet>
muonCorr;


std::shared_ptr<const correction::Correction>
muon_m_data;


std::shared_ptr<const correction::Correction>
muon_a_data;



JetCorrections jetCorr;



// ============================================================================
// HISTOGRAMS
// ============================================================================


TH1F *h_mass =
new TH1F(
"h_mass",
"Dimuon Mass",
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
"Corrected Jet pT",
100,
0,
300
);



TH1F *h_leadJetPt =
new TH1F(
"h_leadJetPt",
"Leading Jet pT",
100,
0,
300
);



TH1F *h_leadJetEta =
new TH1F(
"h_leadJetEta",
"Leading Jet eta",
60,
-5,
5
);



TH1F *h_dijetMass =
new TH1F(
"h_dijetMass",
"Dijet invariant mass",
100,
0,
1000
);



TH1F *h_dijetPt =
new TH1F(
"h_dijetPt",
"Dijet pT",
100,
0,
500
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

    double dphi = phi1 - phi2;


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
    deltaPhi(phi1,phi2);


    return sqrt(
        deta*deta +
        dphi*dphi
    );

}



// ============================================================================
// KIT MUON CORRECTION — DATA
//
// pT_corrected = pT_BS * m_data + a_data
//
// ============================================================================

double getCorrectedMuonPt_DATA(
    double pt,
    double eta,
    double phi
)
{

   double m =
muon_m_data->evaluate({
    eta,
    phi,
    "nom"
});

double a =
muon_a_data->evaluate({
    eta,
    phi,
    "nom"
});


    return pt*m+a;

}



// ============================================================================
// GOLDEN JSON STORAGE
// ============================================================================


std::map<
unsigned int,
std::vector<
std::pair<unsigned int,unsigned int>
>
> goodLumis;



// ============================================================================
// LOAD GOLDEN JSON
// ============================================================================


void LoadGoldenJSON(
    const std::string& filename
)
{

    goodLumis.clear();



    std::ifstream file(filename);



    if(!file.is_open())
    {

        std::cerr
        <<"ERROR: Cannot open Golden JSON: "
        <<filename
        <<std::endl;


        exit(1);

    }



    json j;

    file >> j;



    for(auto& run : j.items())
    {

        unsigned int runNumber =
        std::stoul(run.key());



        for(auto& lumi : run.value())
        {

            unsigned int start =
            lumi[0];


            unsigned int end =
            lumi[1];



            goodLumis[runNumber]
            .push_back(
                {
                    start,
                    end
                }
            );

        }

    }



    std::cout
    <<"Golden JSON loaded. Runs = "
    <<goodLumis.size()
    <<std::endl;


}




// ============================================================================
// CHECK GOLDEN JSON
// ============================================================================


bool PassGoldenJSON(
    unsigned int run,
    unsigned int lumi
)
{


    auto it =
    goodLumis.find(run);



    if(it == goodLumis.end())
        return false;



    for(auto& range : it->second)
    {

        if(
            lumi >= range.first &&
            lumi <= range.second
        )
        {
            return true;
        }

    }



    return false;

}
// ============================================================================
// MAIN FUNCTION
// ============================================================================


void testmacro_FullCorrections_FullJetsMuons_2025kit_bs_golden_vlast(

    std::vector<std::string> inputFiles,


    bool isData=true,


    const std::string& KITDir =
    "/afs/cern.ch/user/n/nbostan/new_CMS/"
    "CMSSW_14_0_18/src/RoccoR/post2022E-update/"

)
{


    gErrorIgnoreLevel = kError;



    // ========================================================================
    // LOAD KIT MUON CORRECTIONS
    // ========================================================================


    std::cout
    <<"Loading KIT muon corrections..."
    <<std::endl;



    muonCorr =
    correction::CorrectionSet::from_file(

        KITDir +
        "schemaV2_2025.json"

    );



    muon_m_data =
    muonCorr->at(
        "m_data"
    );


    muon_a_data =
    muonCorr->at(
        "a_data"
    );



    std::cout
    <<"KIT muon correction loaded"
    <<std::endl;




    // ========================================================================
    // LOAD JET CORRECTIONS
    // ========================================================================


    jetCorr =
    JetCorrections(

        KITDir+"jet_jerc.txt",

        KITDir+"jetid.txt",

        KITDir+"jetvetomaps.txt"

    );



    std::cout
    <<"Jet corrections loaded"
    <<std::endl;





    // ========================================================================
    // LOAD GOLDEN JSON
    // ========================================================================


    if(isData)
    {

        LoadGoldenJSON(

        "/cvmfs/cms-griddata.cern.ch/cat/"
        "metadata/DC/Collisions25/latest/"
        "Cert_Collisions2025_391658_398903_Muon.json"

        );

    }





    // ========================================================================
    // HISTOGRAM ERRORS
    // ========================================================================


    for(auto h :
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





    // ========================================================================
    // BUILD TCHAIN
    // ========================================================================


    TChain chain(
        "Events"
    );



    for(auto const& file : inputFiles)
    {

        chain.Add(
            file.c_str()
        );

    }



    std::cout
    <<"Total entries = "
    <<chain.GetEntries()
    <<std::endl;





    // ========================================================================
    // TREE READER
    // ========================================================================


    TTreeReader reader(
        &chain
    );
    
    // ============================================================================
// EVENT INFORMATION
// ============================================================================


TTreeReaderValue<UInt_t>
run(
    reader,
    "run"
);



TTreeReaderValue<UInt_t>
luminosityBlock(
    reader,
    "luminosityBlock"
);




// ============================================================================
// DATA QUALITY FLAGS
// ============================================================================


TTreeReaderValue<bool>
Flag_goodVertices(
    reader,
    "Flag_goodVertices"
);



TTreeReaderValue<bool>
Flag_globalSuperTightHalo2016Filter(
    reader,
    "Flag_globalSuperTightHalo2016Filter"
);



TTreeReaderValue<bool>
Flag_EcalDeadCellTriggerPrimitiveFilter(
    reader,
    "Flag_EcalDeadCellTriggerPrimitiveFilter"
);



TTreeReaderValue<bool>
Flag_BadPFMuonFilter(
    reader,
    "Flag_BadPFMuonFilter"
);



TTreeReaderValue<bool>
Flag_BadPFMuonDzFilter(
    reader,
    "Flag_BadPFMuonDzFilter"
);



TTreeReaderValue<bool>
Flag_hfNoisyHitsFilter(
    reader,
    "Flag_hfNoisyHitsFilter"
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
Muon_mediumId(
    reader,
    "Muon_mediumId"
);





// ============================================================================
// BEAMSPOT CONSTRAINED MUON
// ============================================================================


bool hasBS = false;



if(chain.GetBranch("Muon_bsConstrainedPt"))
{
    hasBS = true;
}



TTreeReaderArray<float>*
Muon_bsConstrainedPt = nullptr;



TTreeReaderArray<float>*
Muon_bsConstrainedChi2 = nullptr;




if(hasBS)
{

    Muon_bsConstrainedPt =
    new TTreeReaderArray<float>(
        reader,
        "Muon_bsConstrainedPt"
    );



    Muon_bsConstrainedChi2 =
    new TTreeReaderArray<float>(
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



TTreeReaderArray<ULong64_t>*
TrigObj_filterBits = nullptr;



if(chain.GetBranch("TrigObj_filterBits"))
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
// CUTFLOW COUNTERS
// ============================================================================


Long64_t nTotal      = 0;

Long64_t nGolden     = 0;

Long64_t nMETFilter  = 0;

Long64_t nTrigger    = 0;

Long64_t nTwoMuon    = 0;

Long64_t nMuonSel    = 0;

Long64_t nTrigMatch  = 0;

Long64_t nFinal      = 0;
// ============================================================================
// EVENT LOOP
// ============================================================================


while(reader.Next())
{


    nTotal++;




    // ========================================================================
    // GOLDEN JSON
    // ========================================================================


    if(isData)
    {

        if(
            !PassGoldenJSON(
                *run,
                *luminosityBlock
            )
        )
            continue;


        nGolden++;

    }





    // ========================================================================
    // DATA QUALITY FILTERS
    // ========================================================================


    if(
        !(*Flag_goodVertices) ||
        !(*Flag_globalSuperTightHalo2016Filter) ||
        !(*Flag_EcalDeadCellTriggerPrimitiveFilter) ||
        !(*Flag_BadPFMuonFilter) ||
        !(*Flag_BadPFMuonDzFilter) ||
        !(*Flag_hfNoisyHitsFilter)
    )
        continue;



    nMETFilter++;






    // ========================================================================
    // TRIGGER
    // ========================================================================


    if(!(*HLT_IsoMu24))
        continue;


    nTrigger++;






    // ========================================================================
    // EXACTLY TWO MUONS
    // ========================================================================


    if(*nMuon != 2)
        continue;


    nTwoMuon++;




    double weight = 1.0;





    // ========================================================================
    // BEAMSPOT CONSTRAINED MUON PT
    // ========================================================================


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





    // ========================================================================
    // KIT MUON CORRECTION
    //
    // pT = pT_BS*m_data+a_data
    //
    // ========================================================================


    double corrPt[2];

    bool passMuon=true;



    for(int i=0;i<2;i++)
    {


        corrPt[i] =
        getCorrectedMuonPt_DATA(

            muPt[i],

            Muon_eta[i],

            Muon_phi[i]

        );


        // pT

        if(
            corrPt[i] < 20.
        )
            passMuon=false;



        // eta

        if(
            fabs(Muon_eta[i]) > 2.4
        )
            passMuon=false;



        // ID

        if(
            !Muon_mediumId[i]
        )
            passMuon=false;



        // isolation

        if(
            Muon_iso[i] > 0.25
        )
            passMuon=false;


    }




    if(!passMuon)
        continue;



    nMuonSel++;







    // ========================================================================
    // TAG MUON
    // ========================================================================


    int tag=-1;

    double leadPt=-1;



    for(int i=0;i<2;i++)
    {


        if(
            corrPt[i] > 26. &&
            corrPt[i] > leadPt
        )
        {

            leadPt =
            corrPt[i];

            tag=i;

        }

    }




    if(tag<0)
        continue;



    int probe = 1-tag;







    // ========================================================================
    // TRIGGER MATCHING
    // ========================================================================


    bool triggerMatched=false;



    for(int i=0;i<*nTrigObj;i++)
    {


        if(
            TrigObj_id[i]!=13
        )
            continue;



        if(
            TrigObj_pt[i]<24.
        )
            continue;




        if(TrigObj_filterBits)
        {

            if(
            (((*TrigObj_filterBits)[i]
            &
            (1ULL<<3))==0)
            )
                continue;

        }




        double dr =
        deltaR(

            Muon_eta[tag],

            Muon_phi[tag],

            TrigObj_eta[i],

            TrigObj_phi[i]

        );



        if(dr<0.1)
        {

            triggerMatched=true;

            break;

        }


    }





    if(!triggerMatched)
        continue;



    nTrigMatch++;







    // ========================================================================
    // OPPOSITE SIGN
    // ========================================================================


    if(
        Muon_charge[tag] *
        Muon_charge[probe]
        >=0
    )
        continue;







    // ========================================================================
    // DIMUON FOUR VECTOR
    // ========================================================================


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
    mu1 + mu2;





    // ========================================================================
    // DIMUON HISTOGRAMS
    // ========================================================================


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



for(int j=0; j<*nJet; j++)
{


    // ========================================================================
    // JEC CORRECTION — DATA
    // ========================================================================


    double corrJetPt =
    jetCorr.getCorrectedJetPt(

        Jet_pt[j],

        Jet_eta[j],

        Jet_phi[j],

        Jet_rawFactor[j],

        false   // DATA

    );



    if(corrJetPt < 20.)
        continue;



    if(
        fabs(Jet_eta[j]) > 4.7
    )
        continue;





    TLorentzVector jet;



    jet.SetPtEtaPhiM(

        corrJetPt,

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
    // JET VETO MAP
    // ========================================================================


    if(
        !jetCorr.passJetVetoMap(
            Jet_eta[j],
            Jet_phi[j]
        )
    )
    {
        continue;
    }





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





    // keep jet

    jets.push_back(jet);




    h_jetPtCorr->Fill(

        corrJetPt,

        weight

    );


}







// ============================================================================
// SORT JETS BY PT
// ============================================================================


std::sort(

    jets.begin(),

    jets.end(),

    [](const TLorentzVector &a,
       const TLorentzVector &b)

    {

        return a.Pt() > b.Pt();

    }

);








// ============================================================================
// LEADING JET
// ============================================================================


if(jets.size()>0)
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






if(jets.size()>=2)
{


    TLorentzVector j1 =
    jets[0];


    TLorentzVector j2 =
    jets[1];



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





    // VBF selection

    if(

        mjj > 400.

        &&

        deta > 2.5

        &&

        j1.Pt() > 35.

        &&

        j2.Pt() > 25.

    )
    {

        isVBF=true;

    }

    else
    {

        isGGH=true;

    }


}



else if(jets.size()==1)
{

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
// CLOSE EVENT LOOP
// ============================================================================


} // while(reader.Next())






// ============================================================================
// OUTPUT ROOT FILE
// ============================================================================


TFile out(

    "/eos/user/n/nbostan/2025_Samples/"
    "output_histos_DATA_2025_KIT_bs_golden.root",

    "RECREATE"

);





if(out.IsZombie())
{

    std::cerr
    <<"ERROR: Cannot create output ROOT file"
    <<std::endl;


    return;

}





// ============================================================================
// WRITE HISTOGRAMS
// ============================================================================


std::cout
<<"\nWriting histograms..."
<<std::endl;



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







std::cout
<<"\n====================================\n"
<<" DATA ANALYSIS FINISHED\n"
<<" KIT MUON CORRECTION USED\n"
<<" BeamSpot muons used\n"
<<" Golden JSON applied\n"
<<" Output:\n"
<<" output_histos_DATA_2025_KIT_bs_golden.root\n"
<<"====================================\n";

}
