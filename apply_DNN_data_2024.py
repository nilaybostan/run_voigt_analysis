import numpy as np
import awkward as ak
import uproot

import torch
import torch.nn as nn

import joblib
import json
import glob

import matplotlib.pyplot as plt



# ============================================================
# CONFIGURATION
# ============================================================

MODEL_FILE = "/afs/cern.ch/user/n/nbostan/new_CMS/CMSSW_14_0_18/src/RoccoR/condor_BDT_DNN/large_BDT/best_model_2024.pth"
SCALER_FILE = "/afs/cern.ch/user/n/nbostan/new_CMS/CMSSW_14_0_18/src/RoccoR/condor_BDT_DNN/large_BDT/scaler_2024.pkl"
THRESHOLD_FILE = "/afs/cern.ch/user/n/nbostan/new_CMS/CMSSW_14_0_18/src/RoccoR/condor_BDT_DNN/large_BDT/threshold_2024.json"

OUTPUT_FILE = "/afs/cern.ch/user/n/nbostan/new_CMS/CMSSW_14_0_18/src/RoccoR/condor_BDT_DNN/large_BDT/Data_DNNScore_2024.root"


PLOT_DIR = "/afs/cern.ch/user/n/nbostan/new_CMS/CMSSW_14_0_18/src/RoccoR/condor_BDT_DNN/large_BDT/plots/"



# ============================================================
# INPUT DATA FILES
# ============================================================


DATA_FILES = ["root://cms-xrd-global.cern.ch//store/data/Run2024C/Muon0/NANOAOD/MINIv6NANOv15-v1/2530000/677e3bb0-8199-4ffb-83af-165410a7b7a6.root","root://cms-xrd-global.cern.ch//store/data/Run2024C/Muon0/NANOAOD/MINIv6NANOv15-v1/2530000/a4a34841-e26f-455c-8984-06703b2e23d7.root","root://cms-xrd-global.cern.ch//store/data/Run2024C/Muon0/NANOAOD/MINIv6NANOv15-v1/2530000/94b5a3dc-54ea-418a-8fa3-41a9e38e2cf2.root","root://cms-xrd-global.cern.ch//store/data/Run2024C/Muon0/NANOAOD/MINIv6NANOv15-v1/2530000/9be46da1-631e-409a-9218-075ea4537cad.root","root://cms-xrd-global.cern.ch//store/data/Run2024C/Muon0/NANOAOD/MINIv6NANOv15-v1/2530000/13bdae97-f65a-418b-aed8-2ec264d42be0.root","root://cms-xrd-global.cern.ch//store/data/Run2024C/Muon0/NANOAOD/MINIv6NANOv15-v1/2530000/390b1d89-09d1-406f-8f6c-417d59419d65.root","root://cms-xrd-global.cern.ch//store/data/Run2024C/Muon0/NANOAOD/MINIv6NANOv15-v1/2530000/d26e3d05-22e5-4e89-ab02-98181acd50f2.root","root://cms-xrd-global.cern.ch//store/data/Run2024C/Muon0/NANOAOD/MINIv6NANOv15-v1/2530000/fb8a66c0-e9e8-4ddf-8810-5f5b71ebe1fe.root","root://cms-xrd-global.cern.ch//store/data/Run2024C/Muon0/NANOAOD/MINIv6NANOv15-v1/2530000/8c8ea958-4a49-41cc-89cd-9fd583500f9c.root","root://cms-xrd-global.cern.ch//store/data/Run2024C/Muon0/NANOAOD/MINIv6NANOv15-v1/2530000/fc9f1c1a-6699-4d70-a04b-40d501ad4373.root","root://cms-xrd-global.cern.ch//store/data/Run2024C/Muon0/NANOAOD/MINIv6NANOv15-v1/2530000/3622c672-ed80-448e-8e8f-58a65506166f.root","root://cms-xrd-global.cern.ch//store/data/Run2024C/Muon0/NANOAOD/MINIv6NANOv15-v1/2530000/3343ce06-8bfa-49a9-a2f4-276504e26307.root","root://cms-xrd-global.cern.ch//store/data/Run2024C/Muon0/NANOAOD/MINIv6NANOv15-v1/2530000/ab7c36e4-66c9-481b-8d88-114ab2c818ac.root","root://cms-xrd-global.cern.ch//store/data/Run2024C/Muon0/NANOAOD/MINIv6NANOv15-v1/2530000/202efc0a-927a-475a-b4b6-7ed9aabad668.root","root://cms-xrd-global.cern.ch//store/data/Run2024C/Muon0/NANOAOD/MINIv6NANOv15-v1/2530000/e2546003-cca7-4229-a540-d56b807399c9.root","root://cms-xrd-global.cern.ch//store/data/Run2024D/Muon0/NANOAOD/MINIv6NANOv15-v1/2540000/59956cf4-d632-4fa4-b930-4483e7a7f81c.root","root://cms-xrd-global.cern.ch//store/data/Run2024D/Muon0/NANOAOD/MINIv6NANOv15-v1/110000/9ab649c7-bdb6-46f9-9b76-74f61e40fa47.root","root://cms-xrd-global.cern.ch//store/data/Run2024D/Muon0/NANOAOD/MINIv6NANOv15-v1/110000/d190f0d2-d9a1-4ddf-9431-5f3943fb49a8.root","root://cms-xrd-global.cern.ch//store/data/Run2024D/Muon0/NANOAOD/MINIv6NANOv15-v1/110000/cc76ab1f-6a2f-4ca6-a393-f2fe16a4bd8d.root","root://cms-xrd-global.cern.ch//store/data/Run2024D/Muon0/NANOAOD/MINIv6NANOv15-v1/110000/713cc4ef-72a3-4cce-8dd6-8d24c080ddd3.root","root://cms-xrd-global.cern.ch//store/data/Run2024D/Muon0/NANOAOD/MINIv6NANOv15-v1/2540000/030d0310-0476-45ff-87b6-4752bdfac00d.root","root://cms-xrd-global.cern.ch//store/data/Run2024D/Muon0/NANOAOD/MINIv6NANOv15-v1/110000/a7b201ca-7c2e-4974-8712-d32d7a04f396.root","root://cms-xrd-global.cern.ch//store/data/Run2024D/Muon0/NANOAOD/MINIv6NANOv15-v1/110000/bc389e26-d603-4978-b3b0-17cb627ed052.root","root://cms-xrd-global.cern.ch//store/data/Run2024D/Muon0/NANOAOD/MINIv6NANOv15-v1/110000/b605eb6f-4b88-4ee0-bd9a-c0a7e4576eb5.root","root://cms-xrd-global.cern.ch//store/data/Run2024D/Muon0/NANOAOD/MINIv6NANOv15-v1/110000/6c28d2ba-a32f-4661-a1e8-8a743585604f.root","root://cms-xrd-global.cern.ch//store/data/Run2024E/Muon0/NANOAOD/MINIv6NANOv15-v1/2520000/7531adac-d24d-4784-b028-e8b988bc578c.root","root://cms-xrd-global.cern.ch//store/data/Run2024E/Muon0/NANOAOD/MINIv6NANOv15-v1/2520000/6423ef94-41bd-426e-ae57-1c9eeee26313.root","root://cms-xrd-global.cern.ch//store/data/Run2024E/Muon0/NANOAOD/MINIv6NANOv15-v1/2520000/612c7579-9520-4132-8a9f-051637fe6bae.root","root://cms-xrd-global.cern.ch//store/data/Run2024E/Muon0/NANOAOD/MINIv6NANOv15-v1/2520000/4a99092d-5b9e-4c98-965b-1f81eda5fb8e.root","root://cms-xrd-global.cern.ch//store/data/Run2024E/Muon0/NANOAOD/MINIv6NANOv15-v1/2810000/b9fc5dcd-10af-48c1-bc75-7880d0fcef9e.root","root://cms-xrd-global.cern.ch//store/data/Run2024E/Muon0/NANOAOD/MINIv6NANOv15-v1/2810000/df135a40-af5a-42f8-a270-3cb7c0bb6d1d.root","root://cms-xrd-global.cern.ch//store/data/Run2024E/Muon0/NANOAOD/MINIv6NANOv15-v1/2810000/76078b63-6241-467d-80cc-05d316fcefe5.root","root://cms-xrd-global.cern.ch//store/data/Run2024E/Muon0/NANOAOD/MINIv6NANOv15-v1/2520000/c8ecc5bf-d3c8-4135-9496-5f72ced2c15d.root","root://cms-xrd-global.cern.ch//store/data/Run2024E/Muon0/NANOAOD/MINIv6NANOv15-v1/2520000/30f91c56-9c39-4349-a37d-cd292c3502f6.root","root://cms-xrd-global.cern.ch//store/data/Run2024E/Muon0/NANOAOD/MINIv6NANOv15-v1/2520000/c679b348-9bc0-444a-b49e-58f28274e085.root","root://cms-xrd-global.cern.ch//store/data/Run2024F/Muon0/NANOAOD/MINIv6NANOv15-v1/2540000/7f4a3771-0ef4-4a84-a530-a0d4574056d4.root","root://cms-xrd-global.cern.ch//store/data/Run2024F/Muon0/NANOAOD/MINIv6NANOv15-v1/2540000/4eda15c6-786a-4381-bead-7b5d18c41a2c.root","root://cms-xrd-global.cern.ch//store/data/Run2024F/Muon0/NANOAOD/MINIv6NANOv15-v1/2540000/e1fd7044-4993-4469-8258-55bd6ddeb1b2.root","root://cms-xrd-global.cern.ch//store/data/Run2024F/Muon0/NANOAOD/MINIv6NANOv15-v1/2530000/9af32129-aef3-424e-b2ef-36689fd0e506.root","root://cms-xrd-global.cern.ch//store/data/Run2024F/Muon0/NANOAOD/MINIv6NANOv15-v1/2530000/a4a69c55-2f56-4787-9695-93d59744e431.root","root://cms-xrd-global.cern.ch//store/data/Run2024F/Muon0/NANOAOD/MINIv6NANOv15-v1/2530000/6ba191d0-8e16-4f87-9483-4b9e11d1ccc9.root","root://cms-xrd-global.cern.ch//store/data/Run2024F/Muon0/NANOAOD/MINIv6NANOv15-v1/2530000/6ac104c1-e2ef-409d-8c1d-ace5167a32f3.root","root://cms-xrd-global.cern.ch//store/data/Run2024F/Muon0/NANOAOD/MINIv6NANOv15-v1/2530000/d4e5ef15-36a1-4e54-92e8-35f446a41aa1.root","root://cms-xrd-global.cern.ch//store/data/Run2024F/Muon0/NANOAOD/MINIv6NANOv15-v1/2530000/fbef915c-0469-47cf-830e-70938dcdd3d7.root","root://cms-xrd-global.cern.ch//store/data/Run2024F/Muon0/NANOAOD/MINIv6NANOv15-v1/2530000/71250f00-226c-4b83-87fe-1ae5ca0db899.root","root://cms-xrd-global.cern.ch//store/data/Run2024G/Muon0/NANOAOD/MINIv6NANOv15-v1/2540000/03f2d829-00be-4c48-8d87-49ffc6825a92.root","root://cms-xrd-global.cern.ch//store/data/Run2024G/Muon0/NANOAOD/MINIv6NANOv15-v1/2540000/4a74a403-1653-4ac5-8896-ab6f961ad32b.root","root://cms-xrd-global.cern.ch//store/data/Run2024G/Muon0/NANOAOD/MINIv6NANOv15-v1/2540000/f8bd4af4-689e-48bd-a32f-49d119937bfb.root","root://cms-xrd-global.cern.ch//store/data/Run2024G/Muon0/NANOAOD/MINIv6NANOv15-v1/2540000/ef1e01fd-e68e-4388-a91e-b720e9ea8046.root","root://cms-xrd-global.cern.ch//store/data/Run2024G/Muon0/NANOAOD/MINIv6NANOv15-v1/2540000/fd29aa87-bc34-468a-963f-452e30a621bb.root","root://cms-xrd-global.cern.ch//store/data/Run2024G/Muon0/NANOAOD/MINIv6NANOv15-v1/2540000/ac50650d-bf19-483d-8bfe-1bd8ab35877e.root","root://cms-xrd-global.cern.ch//store/data/Run2024G/Muon0/NANOAOD/MINIv6NANOv15-v1/2540000/c6c64c05-c964-4f58-8aab-08e17749342a.root","root://cms-xrd-global.cern.ch//store/data/Run2024G/Muon0/NANOAOD/MINIv6NANOv15-v1/2540000/a39769ef-fffa-44e7-b8e5-b088662dc28c.root","root://cms-xrd-global.cern.ch//store/data/Run2024G/Muon0/NANOAOD/MINIv6NANOv15-v1/2540000/9fb277a0-8828-4b7e-93a8-89dfcc1486fb.root","root://cms-xrd-global.cern.ch//store/data/Run2024G/Muon0/NANOAOD/MINIv6NANOv15-v1/2540000/5bac465c-7c4a-4342-ba02-56248f756e0c.root","root://cms-xrd-global.cern.ch//store/data/Run2024H/Muon0/NANOAOD/MINIv6NANOv15-v1/90000/244fa588-5f93-4a47-b711-32b5ded4852d.root","root://cms-xrd-global.cern.ch//store/data/Run2024H/Muon0/NANOAOD/MINIv6NANOv15-v1/90000/0aab5444-026c-4ac7-a629-59cf11182aee.root","root://cms-xrd-global.cern.ch//store/data/Run2024H/Muon0/NANOAOD/MINIv6NANOv15-v1/90000/22d6e4c3-1aa2-4162-ada3-05869758bc62.root","root://cms-xrd-global.cern.ch//store/data/Run2024H/Muon0/NANOAOD/MINIv6NANOv15-v1/90000/c11fb5b7-a95f-4c89-9cde-2324a10e086e.root","root://cms-xrd-global.cern.ch//store/data/Run2024H/Muon0/NANOAOD/MINIv6NANOv15-v1/90000/7235c962-57ea-434a-812f-ab2dae5d51c8.root","root://cms-xrd-global.cern.ch//store/data/Run2024H/Muon0/NANOAOD/MINIv6NANOv15-v1/90000/420235c8-4adc-47ae-affa-0cfc85beb125.root","root://cms-xrd-global.cern.ch//store/data/Run2024H/Muon0/NANOAOD/MINIv6NANOv15-v1/90000/3799ec42-aabf-4a85-8fec-f6628d8e1dd2.root","root://cms-xrd-global.cern.ch//store/data/Run2024H/Muon0/NANOAOD/MINIv6NANOv15-v1/90000/eb6a523d-6ff4-4c39-96b5-b3e026673849.root","root://cms-xrd-global.cern.ch//store/data/Run2024H/Muon0/NANOAOD/MINIv6NANOv15-v1/90000/406d9ec2-1a07-4bb8-a73b-87c779282ed5.root","root://cms-xrd-global.cern.ch//store/data/Run2024H/Muon0/NANOAOD/MINIv6NANOv15-v1/90000/6b2ba2fa-224e-4918-bfd6-3a7d21b305f4.root","root://cms-xrd-global.cern.ch//store/data/Run2024I/Muon0/NANOAOD/MINIv6NANOv15_v2-v1/2530000/6c392d88-b78b-42bb-a6af-44746fdcff11.root","root://cms-xrd-global.cern.ch//store/data/Run2024I/Muon0/NANOAOD/MINIv6NANOv15_v2-v1/2530000/21d16008-96bc-4096-82b8-2247844dfc8e.root","root://cms-xrd-global.cern.ch//store/data/Run2024I/Muon0/NANOAOD/MINIv6NANOv15_v2-v1/2530000/47d3b12b-ab93-41c0-a284-e5dc0c848047.root","root://cms-xrd-global.cern.ch//store/data/Run2024I/Muon0/NANOAOD/MINIv6NANOv15_v2-v1/2530000/a86c55a7-cb18-4ec3-a321-ec5aded80dee.root","root://cms-xrd-global.cern.ch//store/data/Run2024I/Muon0/NANOAOD/MINIv6NANOv15_v2-v1/2530000/a0d2f532-14d2-4f53-b94a-6ef65b2d8faa.root","root://cms-xrd-global.cern.ch//store/data/Run2024I/Muon0/NANOAOD/MINIv6NANOv15_v2-v1/2530000/d1d401cf-0635-488c-8d4d-9c3b17543a12.root","root://cms-xrd-global.cern.ch//store/data/Run2024I/Muon0/NANOAOD/MINIv6NANOv15_v2-v1/2530000/ad5a3efe-9b4e-4228-a33f-bc7eeb8c4cd1.root","root://cms-xrd-global.cern.ch//store/data/Run2024I/Muon0/NANOAOD/MINIv6NANOv15_v2-v1/2530000/83560a0d-1b88-4640-808f-ef7a0a6dd34f.root","root://cms-xrd-global.cern.ch//store/data/Run2024I/Muon0/NANOAOD/MINIv6NANOv15_v2-v1/2530000/448fbb8f-4534-4636-863b-58ede98b5bbf.root","root://cms-xrd-global.cern.ch//store/data/Run2024I/Muon0/NANOAOD/MINIv6NANOv15_v2-v1/2530000/77ff9f7e-9f5f-4484-b98c-9e4dbeaa17e3.root"
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

"PuppiMET_pt",
"PuppiMET_phi",

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

ak.to_numpy(data.PuppiMET_pt),

np.sin(
ak.to_numpy(data.PuppiMET_phi)
),

np.cos(
ak.to_numpy(data.PuppiMET_phi)
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

plt.yscale("log") 
plt.xlabel("DNN score")
plt.ylabel("Events")

plt.grid(alpha=0.3)

plt.savefig(
PLOT_DIR+"DNNScore_Data_2024.png",
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

plt.yscale("log") 
plt.xlabel("DNN score")
plt.ylabel("Events")

plt.legend()

plt.grid(alpha=0.3)


plt.savefig(
PLOT_DIR+"DNNScore_Data_Threshold_2024.png",
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


plt.yscale("log") 
plt.xlabel("DNN Selection")
plt.ylabel("Events")

plt.savefig(
PLOT_DIR+"DNNPass_2024.png",
dpi=300
)

plt.close()



# ------------------------------------------------------------
# KINEMATICS
# ------------------------------------------------------------


plots=[

    (
        mll,
        "$m_{\\mu\\mu}$ (GeV)",
        "mll_Data_2024.png"
    ),

    (
        mjj,
        "$m_{jj}$ (GeV)",
        "mjj_Data_2024.png"
    ),

    (
        dEta_jj,
        "$\\Delta\\eta_{jj}$",
        "DeltaEtaJJ_Data_2024.png"
    ),

    (
        ak.to_numpy(mu1pt[sel]),
        "Leading muon pT",
        "MuonPt1_Data_2024.png"
    ),

    (
        ak.to_numpy(mu2pt[sel]),
        "Subleading muon pT",
        "MuonPt2_Data_2024.png"
    )

]



for values,xlabel,name in plots:


    plt.figure(figsize=(7,6))


    plt.hist(
        values,
        bins=50,
        histtype="step",
        linewidth=2
    )


    plt.xlabel(xlabel)
    plt.ylabel("Events")

    plt.grid(alpha=0.3)

    plt.tight_layout()

    plt.savefig(
        PLOT_DIR+name,
        dpi=300
    )

    plt.close()



print("\nAll plots produced successfully!")
