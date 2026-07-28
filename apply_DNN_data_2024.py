import numpy as np
import awkward as ak
import uproot

import torch
import torch.nn as nn

import joblib
import json
import glob

import matplotlib.pyplot as plt

import datetime

# ============================================================
# CONFIGURATION
# ============================================================

MODEL_FILE = "/afs/cern.ch/user/n/nbostan/new_CMS/CMSSW_14_0_18/src/RoccoR/condor_BDT_DNN/large_BDT/best_model_pre2022.pth"
SCALER_FILE = "/afs/cern.ch/user/n/nbostan/new_CMS/CMSSW_14_0_18/src/RoccoR/condor_BDT_DNN/large_BDT/scaler_pre2022.pkl"
THRESHOLD_FILE = "/afs/cern.ch/user/n/nbostan/new_CMS/CMSSW_14_0_18/src/RoccoR/condor_BDT_DNN/large_BDT/threshold_pre2022.json"

OUTPUT_FILE = "/eos/user/n/nbostan/plots/Data_DNNScore_pre2022.root"


PLOT_DIR = "/eos/user/n/nbostan/plots/"
SUMMARY_FILE = PLOT_DIR + "DNN_summary_pre2022.txt"


# ============================================================
# INPUT DATA FILES
# ============================================================


DATA_FILES = ["root://cms-xrd-global.cern.ch//store/data/Run2022C/Muon/NANOAOD/16Dec2023-v1/2550000/99577992-2ea3-4ad3-8e87-69351341b12d.root","root://cms-xrd-global.cern.ch//store/data/Run2022C/Muon/NANOAOD/16Dec2023-v1/2550000/e664be04-0621-4ca2-889f-869dde5b763b.root","root://cms-xrd-global.cern.ch//store/data/Run2022C/Muon/NANOAOD/16Dec2023-v1/2550000/f1bfe9dc-259c-4a11-b436-749e6a18ddea.root","root://cms-xrd-global.cern.ch//store/data/Run2022C/Muon/NANOAOD/16Dec2023-v1/2550000/f0eb2bf9-d610-4c81-9fdd-07b7fd247cf4.root","root://cms-xrd-global.cern.ch//store/data/Run2022C/Muon/NANOAOD/16Dec2023-v1/2550000/3220e70c-61a7-4dbb-a479-ce1d0948976c.root","root://cms-xrd-global.cern.ch//store/data/Run2022C/Muon/NANOAOD/16Dec2023-v1/2550000/ab62ccfc-efcb-49ed-8334-ce8d716be234.root","root://cms-xrd-global.cern.ch//store/data/Run2022C/Muon/NANOAOD/16Dec2023-v1/2550000/30824fd2-bbc1-4c06-bc29-4214d534be48.root","root://cms-xrd-global.cern.ch//store/data/Run2022C/Muon/NANOAOD/16Dec2023-v1/2550000/c4c82c99-25e6-4bb7-9976-b1c814811665.root","root://cms-xrd-global.cern.ch//store/data/Run2022C/Muon/NANOAOD/16Dec2023-v1/2550000/257de6b5-ad85-43cb-b40e-24e209040b41.root","root://cms-xrd-global.cern.ch//store/data/Run2022C/Muon/NANOAOD/16Dec2023-v1/2550000/69f9d1c8-e490-49f7-b59a-99ab2f643ac4.root","root://cms-xrd-global.cern.ch//store/data/Run2022C/Muon/NANOAOD/16Dec2023-v1/2550000/f4488c49-5873-49d1-96d8-6ce84b4e8aeb.root","root://cms-xrd-global.cern.ch//store/data/Run2022C/Muon/NANOAOD/16Dec2023-v1/2550000/5c373a95-ef61-4784-ac8d-be899595d681.root","root://cms-xrd-global.cern.ch//store/data/Run2022C/Muon/NANOAOD/16Dec2023-v1/2550000/a235b9ac-0342-472b-97ec-15fceac0fada.root","root://cms-xrd-global.cern.ch//store/data/Run2022C/Muon/NANOAOD/16Dec2023-v1/2550000/e100ae34-ca16-41bc-8004-59f0e82ee94b.root","root://cms-xrd-global.cern.ch//store/data/Run2022C/Muon/NANOAOD/16Dec2023-v1/2550000/dd091748-99e0-4d48-8ab3-0fbe4544f07b.root","root://cms-xrd-global.cern.ch//store/data/Run2022C/Muon/NANOAOD/16Dec2023-v1/2550000/164fbe32-be2c-49b2-b58f-9fc0c8fec268.root","root://cms-xrd-global.cern.ch//store/data/Run2022C/Muon/NANOAOD/16Dec2023-v1/2550000/4cea3500-9607-4e91-9fb4-cc71ef2eb3df.root","root://cms-xrd-global.cern.ch//store/data/Run2022C/Muon/NANOAOD/16Dec2023-v1/2550000/41b17077-d6f4-4655-83b3-b7222ab9e06d.root","root://cms-xrd-global.cern.ch//store/data/Run2022C/Muon/NANOAOD/16Dec2023-v1/2550000/8cae1ae4-639a-49ec-b2cc-26cbb5ef14e2.root","root://cms-xrd-global.cern.ch//store/data/Run2022C/Muon/NANOAOD/16Dec2023-v1/2550000/8bcece20-cd24-4d35-85eb-8bd8d7d8f3c3.root","root://cms-xrd-global.cern.ch//store/data/Run2022C/Muon/NANOAOD/16Dec2023-v1/2550000/765717a9-da8a-4eaf-afe3-17a8dad9621b.root","root://cms-xrd-global.cern.ch//store/data/Run2022C/Muon/NANOAOD/16Dec2023-v1/2550000/9a7df399-4e8d-41c6-8b65-98dafad7a871.root","root://cms-xrd-global.cern.ch//store/data/Run2022C/Muon/NANOAOD/16Dec2023-v1/2550000/0f7e88f4-95ce-431e-8f01-40480a2ac878.root","root://cms-xrd-global.cern.ch//store/data/Run2022C/Muon/NANOAOD/16Dec2023-v1/2550000/019b0054-1684-410d-9fb2-837f92da4956.root","root://cms-xrd-global.cern.ch//store/data/Run2022C/Muon/NANOAOD/16Dec2023-v1/2550000/ad6b4e1d-acf1-4a58-85f0-47c521787684.root","root://cms-xrd-global.cern.ch//store/data/Run2022C/SingleMuon/NANOAOD/16Dec2023-v1/50000/57749259-0917-494b-b140-0ea02999b16f.root","root://cms-xrd-global.cern.ch//store/data/Run2022C/SingleMuon/NANOAOD/16Dec2023-v1/50000/f6fd2850-c736-4720-b274-96f75408fa35.root","root://cms-xrd-global.cern.ch//store/data/Run2022C/SingleMuon/NANOAOD/16Dec2023-v1/50000/3b6ba063-4f81-44e5-934c-4f301724e5bc.root","root://cms-xrd-global.cern.ch//store/data/Run2022C/SingleMuon/NANOAOD/16Dec2023-v1/50000/2e7fd442-fded-42be-abc0-30f01d30e56e.root","root://cms-xrd-global.cern.ch//store/data/Run2022C/SingleMuon/NANOAOD/16Dec2023-v1/50000/e040babe-5669-471e-af11-f0a8c0c8575c.root","root://cms-xrd-global.cern.ch//store/data/Run2022C/SingleMuon/NANOAOD/16Dec2023-v1/50000/ba50071d-3700-4e09-be47-b96fc36ea58c.root","root://cms-xrd-global.cern.ch//store/data/Run2022C/SingleMuon/NANOAOD/16Dec2023-v1/50000/a85eb637-9ea7-4f6b-b14a-cacec9ae6e4a.root","root://cms-xrd-global.cern.ch//store/data/Run2022C/SingleMuon/NANOAOD/16Dec2023-v1/50000/0a0b3638-91b0-481f-80e8-5db38c1f9961.root","root://cms-xrd-global.cern.ch//store/data/Run2022C/SingleMuon/NANOAOD/16Dec2023-v1/50000/9255782a-4e41-43df-9f4a-cf3c46b3527a.root","root://cms-xrd-global.cern.ch//store/data/Run2022C/SingleMuon/NANOAOD/16Dec2023-v1/50000/46ff4e59-985f-4f63-9266-92f1de0b4cd7.root","root://cms-xrd-global.cern.ch//store/data/Run2022C/SingleMuon/NANOAOD/16Dec2023-v1/50000/20ad69d6-a37e-420a-a396-0932b49e2f0b.root","root://cms-xrd-global.cern.ch//store/data/Run2022C/SingleMuon/NANOAOD/16Dec2023-v1/50000/31f18b8a-334f-41a9-83a9-77f9723b896f.root","root://cms-xrd-global.cern.ch//store/data/Run2022C/SingleMuon/NANOAOD/16Dec2023-v1/50000/4d0327f8-4bb5-478a-a300-83e8686eab5b.root","root://cms-xrd-global.cern.ch//store/data/Run2022C/SingleMuon/NANOAOD/16Dec2023-v1/50000/578b7213-edd2-4821-9e19-318603197cf8.root","root://cms-xrd-global.cern.ch//store/data/Run2022C/SingleMuon/NANOAOD/16Dec2023-v1/50000/2dd3054f-bf93-48d8-a6a5-8e5d89806ca8.root","root://cms-xrd-global.cern.ch//store/data/Run2022D/Muon/NANOAOD/16Dec2023-v1/50000/c9f48dfb-f13b-43b9-8705-b319d1aaf0bf.root","root://cms-xrd-global.cern.ch//store/data/Run2022D/Muon/NANOAOD/16Dec2023-v1/50000/6a06a004-a6c8-4875-bd6f-0286fb249bde.root","root://cms-xrd-global.cern.ch//store/data/Run2022D/Muon/NANOAOD/16Dec2023-v1/50000/1ab7fa69-c84a-47e4-93e6-6d3f78936c74.root"
]


# Example:
#
# DATA_FILES = glob.glob(
#     "/your/path/*.root"
# )


print("Number of files:",len(DATA_FILES))



# ============================================================
# MODEL
# ============================================================


class DNN(nn.Module):
    def __init__(self):
        super().__init__()

        self.add_module(
            "0",
            nn.Linear(13,225)
        )
        self.add_module(
            "1",
            nn.ReLU()
        )
        self.add_module(
            "2",
            nn.Dropout(0.1)
        )
        self.add_module(
            "3",
            nn.Linear(225,100)
        )
        self.add_module(
            "4",
            nn.ReLU()
        )
        self.add_module(
            "5",
            nn.Dropout(0.1)
        )
        self.add_module(
            "6",
            nn.Linear(100,64)
        )
        self.add_module(
            "7",
            nn.ReLU()
        )
        self.add_module(
            "8",
            nn.Dropout(0.1)
        )
        self.add_module(
            "9",
            nn.Linear(64,1)
        )

    def forward(self,x):
        x = self._modules["0"](x)
        x = self._modules["1"](x)
        x = self._modules["2"](x)
        x = self._modules["3"](x)
        x = self._modules["4"](x)
        x = self._modules["5"](x)
        x = self._modules["6"](x)
        x = self._modules["7"](x)
        x = self._modules["8"](x)
        x = self._modules["9"](x)
        return x


# ============================================================
# LOAD MODEL
# ============================================================


model = DNN()

model.load_state_dict(
    torch.load(
        MODEL_FILE,
        map_location="cpu"
    )
)

model.eval()


scaler = joblib.load(
    SCALER_FILE
)


with open(THRESHOLD_FILE) as f:

    threshold = json.load(f)["threshold"]



print("Threshold =",threshold)



# ============================================================
# HELPER
# ============================================================


def alt(arr,i,fill):

    return ak.fill_none(
        ak.pad_none(
            arr,
            i+1,
            clip=True
        )[:,i],
        fill
    )



# ============================================================
# READ ROOT FILES
# ============================================================


branches=[

"nMuon",
"nElectron",

"Muon_charge",
"Muon_pt",
"Muon_eta",
"Muon_phi",

"MET_pt",
"MET_phi",

"Jet_pt",
"Jet_eta",
"Jet_phi",

"nJet"

]



arrays = []

for i, f in enumerate(DATA_FILES):
    print(f"Reading {i+1}/{len(DATA_FILES)}")

    try:
        arr = uproot.open(f)["Events"].arrays(
            branches,
            library="ak"
        )
        arrays.append(arr)

    except Exception as e:
        print(f"FAILED: {f}")
        print(e)
        continue

data = ak.concatenate(arrays)


print(
"Total events:",
len(data)
)




# ============================================================
# FEATURES
# ============================================================


mu1pt = alt(data.Muon_pt,0,0.)
mu2pt = alt(data.Muon_pt,1,0.)

mu1eta = alt(data.Muon_eta,0,0.)
mu2eta = alt(data.Muon_eta,1,0.)

mu1phi = alt(data.Muon_phi,0,0.)
mu2phi = alt(data.Muon_phi,1,0.)

mu1q = alt(data.Muon_charge,0,0)
mu2q = alt(data.Muon_charge,1,0)



mll = np.sqrt(

np.maximum(

0,

2*mu1pt*mu2pt*

(
np.cosh(mu1eta-mu2eta)
-
np.cos(mu1phi-mu2phi)

)

)

)



dphi = np.arccos(

np.clip(
np.cos(mu1phi-mu2phi),
-1,
1)

)


dR = np.sqrt(

(mu1eta-mu2eta)**2+dphi**2

)



j1pt = alt(data.Jet_pt,0,-999.)
j2pt = alt(data.Jet_pt,1,-999.)

j1eta = alt(data.Jet_eta,0,-999.)
j2eta = alt(data.Jet_eta,1,-999.)

j1phi = alt(data.Jet_phi,0,-999.)
j2phi = alt(data.Jet_phi,1,-999.)



mjj=np.full(len(data),-999.)
dEta_jj=np.full(len(data),-999.)


valid=(j1pt>0)&(j2pt>0)


deta=j1eta-j2eta
dphi_jj=j1phi-j2phi


mjj[valid]=np.sqrt(

2*j1pt[valid]*j2pt[valid]*

(
np.cosh(deta[valid])
-
np.cos(dphi_jj[valid])

)

)


dEta_jj[valid]=np.abs(deta[valid])



# ============================================================
# SELECTION
# ============================================================


sel=(

(data.nMuon>=2)

&

(data.nElectron==0)

&

(mu1q!=mu2q)

&

(mu1pt>26)

&

(mu2pt>20)

&

(abs(mu1eta)<2.4)

&

(abs(mu2eta)<2.4)

&

(mll>110)

&

(mll<150)

)


sel=ak.to_numpy(sel)


print(
"Selected:",
np.sum(sel)
)




# ============================================================
# CREATE DNN INPUT
# ============================================================


X=np.column_stack([

ak.to_numpy(mu1pt),

ak.to_numpy(mu2pt),

ak.to_numpy(mu1eta),

ak.to_numpy(mu2eta),

ak.to_numpy(dR),

ak.to_numpy(data.MET_pt),

np.sin(
ak.to_numpy(data.MET_phi)
),

np.cos(
ak.to_numpy(data.MET_phi)
),

ak.to_numpy(j1pt),

ak.to_numpy(j2pt),

dEta_jj,

mjj,

ak.to_numpy(data.nJet)

])



X=X[sel]



# ============================================================
# SCALE
# ============================================================


X_scaled=scaler.transform(X)



# ============================================================
# DNN PREDICTION
# ============================================================


with torch.no_grad():

    score=torch.sigmoid(

        model(
        torch.tensor(
        X_scaled,
        dtype=torch.float32
        )
        )

    ).numpy().ravel()



DNN_pass=(score>threshold)



print(
"Events passing DNN:",
np.sum(DNN_pass)
)



# ============================================================
# SAVE ROOT
# ============================================================


with uproot.recreate(
    OUTPUT_FILE
) as fout:


    fout["Events"]={

        "DNN_score":score,

        "DNN_pass":
        DNN_pass.astype(np.int32),

        "mll":
        ak.to_numpy(mll[sel]),

        "mu1_pt":
        ak.to_numpy(mu1pt[sel]),

        "mu2_pt":
        ak.to_numpy(mu2pt[sel]),

        "mjj":
        mjj[sel],

        "dEta_jj":
        dEta_jj[sel],

        "nJet":
        ak.to_numpy(data.nJet[sel])

    }



print("Saved:",OUTPUT_FILE)


# ============================================================
# SUMMARY VARIABLES
# ============================================================

total_events = len(data)

selected_events = np.sum(sel)

passing_events = np.sum(DNN_pass)

selection_efficiency = (
    100.0 * selected_events / total_events
    if total_events > 0 else 0.0
)

passing_efficiency = (
    100.0 * passing_events / selected_events
    if selected_events > 0 else 0.0
)

overall_efficiency = (
    100.0 * passing_events / total_events
    if total_events > 0 else 0.0
)

# ============================================================
# PLOTS
# ============================================================


import os

os.makedirs(
    PLOT_DIR,
    exist_ok=True
)



# ------------------------------------------------------------
# DNN SCORE
# ------------------------------------------------------------


plt.figure(figsize=(7,6))

plt.hist(
    score,
    bins=50,
    histtype="step",
    linewidth=2
)

#plt.yscale("log") 
plt.xlabel("DNN score")
plt.ylabel("Events")

plt.grid(alpha=0.3)

plt.savefig(
PLOT_DIR+"DNNScore_Data_pre2022.png",
dpi=300
)

plt.close()



# ------------------------------------------------------------
# DNN SCORE + THRESHOLD
# ------------------------------------------------------------


plt.figure(figsize=(7,6))


plt.hist(
    score,
    bins=50,
    histtype="step",
    linewidth=2
)


plt.axvline(
threshold,
linestyle="--",
label=f"Threshold={threshold:.3f}"
)

#plt.yscale("log") 
plt.xlabel("DNN score")
plt.ylabel("Events")

plt.legend()

plt.grid(alpha=0.3)


plt.savefig(
PLOT_DIR+"DNNScore_Data_Threshold_pre2022.png",
dpi=300
)

plt.close()



# ------------------------------------------------------------
# PASS FAIL
# ------------------------------------------------------------


plt.figure(figsize=(5,5))


plt.bar(
["Fail","Pass"],
[
len(score)-np.sum(DNN_pass),
np.sum(DNN_pass)
]
)


#plt.yscale("log") 
plt.xlabel("DNN Selection")
plt.ylabel("Events")

plt.savefig(
PLOT_DIR+"DNNPass_pre2022.png",
dpi=300
)

plt.close()



# ------------------------------------------------------------
# KINEMATICS
# ------------------------------------------------------------


mll_np = ak.to_numpy(mll)

blind_sel = sel & ~((mll_np > 115) & (mll_np < 135))

plots = [

    (
        ak.to_numpy(mll[blind_sel]),
        r"$m_{\mu\mu}$ (GeV)",
        "mll_Data_pre2022.png",
        (110, 150),
        40
    ),

    (
        ak.to_numpy(mll[sel]),
        r"$m_{\mu\mu}$ (GeV)",
        "mll_Data_pre2022_full.png",
        (0, 200),
        55
    ),

    (
        mjj[sel],
        r"$m_{jj}$ (GeV)",
        "mjj_Data_pre2022.png",
        (0, 2000),
        50
    ),

    (
        dEta_jj[sel],
        r"$\Delta\eta_{jj}$",
        "DeltaEtaJJ_Data_pre2022.png",
        (0, 10),
        100
    ),

    (
        ak.to_numpy(mu1pt[sel]),
        "Leading muon $p_T$ (GeV)",
        "MuonPt1_Data_pre2022.png",
        (0, 200),
        40
    ),

    (
        ak.to_numpy(mu2pt[sel]),
        "Subleading muon $p_T$ (GeV)",
        "MuonPt2_Data_pre2022.png",
        (0, 200),
        40
    )

]

for values, xlabel, name, xrange, bins in plots:

    plt.figure(figsize=(7,6))

    plt.hist(
        values,
        bins=bins,
        range=xrange,
        histtype="step",
        linewidth=2
    )

    #plt.yscale("log")
    plt.ylim(0.8, None)          # <-- buraya
    plt.minorticks_on()          # <-- buraya

    plt.xlabel(xlabel)
    plt.ylabel("Events")

    plt.grid(which="both", alpha=0.3)   # <-- buraya

    plt.tight_layout()

    plt.savefig(
        PLOT_DIR + name,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()



print("\nAll plots produced successfully!")

# ============================================================
# SAVE SUMMARY REPORT
# ============================================================

try:

    with open(SUMMARY_FILE, "w") as f:

        f.write("=====================================================\n")
        f.write("         PyTorch DNN Application Summary\n")
        f.write("=====================================================\n\n")

        f.write(f"Date : {datetime.datetime.now()}\n\n")

        f.write("INPUT INFORMATION\n")
        f.write("-----------------------------\n")

        f.write(f"Input ROOT files          : {len(DATA_FILES)}\n")
        f.write(f"Total events read         : {total_events:,}\n")
        f.write(f"Selected events           : {selected_events:,}\n")
        f.write(f"Selection efficiency (%)  : {selection_efficiency:.3f}\n\n")


        f.write("DNN INFORMATION\n")
        f.write("-----------------------------\n")

        f.write(f"DNN threshold             : {threshold:.6f}\n")
        f.write(f"Events passing DNN        : {passing_events:,}\n")
        f.write(f"DNN efficiency (%)        : {passing_efficiency:.3f}\n")


    print(f"Summary written to {SUMMARY_FILE}")


except Exception as e:

    print("ERROR writing summary:")
    print(e)
