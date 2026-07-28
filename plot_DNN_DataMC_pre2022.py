#!/usr/bin/env python3

# ============================================================
# plot_DNN_DataMC_pre2022.py
#
# Data vs MC DNN score plot
#
# Based on:
#   train_pytorch_updated_last_version_split_pre2022.py
#   apply_DNN_data_pre2022.py
#
# Produces:
#   - stacked MC histogram
#   - Data points
#   - ggH signal line
#   - VBF signal line
#   - Data / MC ratio
#   - MC statistical uncertainty band
#
# ============================================================

import os
import json
import time
import numpy as np
import awkward as ak
import uproot

import torch
import torch.nn as nn

import joblib

import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator


# ============================================================
# CONFIGURATION
# ============================================================

ERA = "pre2022"

# ------------------------------------------------------------
# Luminosity
# ------------------------------------------------------------

LUMI_FB = 7.9804
LUMI_PB = LUMI_FB * 1000.0


# ------------------------------------------------------------
# Model files
# ------------------------------------------------------------
#
# These are the same files used in
# apply_DNN_data_pre2022.py
#
# Change this path if necessary.
# ------------------------------------------------------------

MODEL_FILE = (
    "/afs/cern.ch/user/n/nbostan/"
    "new_CMS/CMSSW_14_0_18/src/RoccoR/"
    "condor_BDT_DNN/large_BDT/best_model_pre2022.pth"
)

SCALER_FILE = (
    "/afs/cern.ch/user/n/nbostan/"
    "new_CMS/CMSSW_14_0_18/src/RoccoR/"
    "condor_BDT_DNN/large_BDT/scaler_pre2022.pkl"
)

THRESHOLD_FILE = (
    "/afs/cern.ch/user/n/nbostan/"
    "new_CMS/CMSSW_14_0_18/src/RoccoR/"
    "condor_BDT_DNN/large_BDT/threshold_pre2022.json"
)


# ------------------------------------------------------------
# Output directory
# ------------------------------------------------------------

PLOT_DIR = (
    "/afs/cern.ch/user/n/nbostan/"
    "new_CMS/CMSSW_14_0_18/src/RoccoR/"
    "condor_BDT_DNN/large_BDT/plot/"
)

os.makedirs(PLOT_DIR, exist_ok=True)


OUTPUT_PNG = os.path.join(
    PLOT_DIR,
    "DNN_DataMC_pre2022.png"
)

OUTPUT_PDF = os.path.join(
    PLOT_DIR,
    "DNN_DataMC_pre2022.pdf"
)


# ============================================================
# CROSS SECTIONS
# ============================================================
#
# Exactly as in the training script
# ============================================================

XSECS_PB = {

    "ggH": 0.0135096,

    "VBF": 0.00105742,

    "DY": 5991.00,

    "EWK": 1.421,

    "TT": 923.6,

    "WZ": 7.568,

    "ZZ": 6.788,

}


# ============================================================
# HISTOGRAM CONFIGURATION
# ============================================================

# DNN output is sigmoid probability:
# 0 <= score <= 1

N_BINS = 20

XMIN = 0.0
XMAX = 1.0

BINS = np.linspace(
    XMIN,
    XMAX,
    N_BINS + 1
)

BIN_CENTERS = (
    0.5 *
    (BINS[:-1] + BINS[1:])
)

BIN_WIDTH = (
    BINS[1] - BINS[0]
)


# ============================================================
# DATA FILE
# ============================================================
#
# Instead of re-reading all Data NanoAOD files,
# we use the ROOT file produced by
# apply_DNN_data_pre2022.py
#
# This ROOT file contains:
#
#   DNN_score
#   DNN_pass
#   mll
#   mu1_pt
#   mu2_pt
#   mjj
#   dEta_jj
#   nJet
#
# This is exactly the output structure of the
# Data application script.
# ============================================================

DATA_DNN_FILE = (
    "/eos/user/n/nbostan/"
    "plots/Data_DNNScore_pre2022.root"
)


# ============================================================
# MC FILE LISTS
# ============================================================
#
# IMPORTANT:
#
# These are taken from the same process separation logic
# as the training code.
#
# The easiest and safest way to keep the exact same file list
# is to import the training file as a module.
#
# HOWEVER:
#
# The training script executes training immediately when imported.
#
# Therefore we do NOT import it.
#
# Instead, you should copy the following four lists from your
# training script if your file lists change:
#
# SIGNAL_FILES
# BACKGROUND_FILES
#
# The code below assumes that this plotting script is placed
# in the same directory as the training script and extracts
# the literal list assignments from that file.
#
# This avoids executing the training.
# ============================================================

import ast


TRAINING_SCRIPT = (
    "train_pytorch_updated_last_version_split_pre2022.py"
)


def extract_literal_variable(
    filename,
    variable_name
):
    """
    Extract a literal Python variable assignment
    from the training script without executing it.

    Used for:
        SIGNAL_FILES
        BACKGROUND_FILES
    """

    with open(
        filename,
        "r",
        encoding="utf-8"
    ) as f:

        source = f.read()

    tree = ast.parse(source)

    for node in tree.body:

        if isinstance(
            node,
            ast.Assign
        ):

            for target in node.targets:

                if (
                    isinstance(
                        target,
                        ast.Name
                    )
                    and
                    target.id == variable_name
                ):

                    return ast.literal_eval(
                        node.value
                    )

    raise RuntimeError(
        f"Could not find literal assignment "
        f"'{variable_name}' in {filename}"
    )


print(
    "\n========================================"
)

print(
    "READING MC FILE LISTS"
)

print(
    "========================================"
)


SIGNAL_FILES = extract_literal_variable(
    TRAINING_SCRIPT,
    "SIGNAL_FILES"
)

BACKGROUND_FILES = extract_literal_variable(
    TRAINING_SCRIPT,
    "BACKGROUND_FILES"
)


# Remove duplicates

SIGNAL_FILES = list(
    dict.fromkeys(
        SIGNAL_FILES
    )
)

BACKGROUND_FILES = list(
    dict.fromkeys(
        BACKGROUND_FILES
    )
)


# ============================================================
# PROCESS CLASSIFICATION
# ============================================================

GGH_FILES = [

    f for f in SIGNAL_FILES

    if "GluGluHto2Mu" in f

]


VBF_FILES = [

    f for f in SIGNAL_FILES

    if "VBFHto2Mu" in f

]


DY_FILES = [

    f for f in BACKGROUND_FILES

    if "DYto2L-2Jets" in f

]


EWK_FILES = [

    f for f in BACKGROUND_FILES

    if "EWK_2L2J" in f

]


TT_FILES = [

    f for f in BACKGROUND_FILES

    if "TT_" in f

]


WZ_FILES = [

    f for f in BACKGROUND_FILES

    if "WZto2L2Q" in f

]


ZZ_FILES = [

    f for f in BACKGROUND_FILES

    if "ZZto2L2Q" in f

]


SAMPLES = {

    "ggH": GGH_FILES,

    "VBF": VBF_FILES,

    "DY": DY_FILES,

    "EWK": EWK_FILES,

    "TT": TT_FILES,

    "WZ": WZ_FILES,

    "ZZ": ZZ_FILES,

}


print(
    "\n========================================"
)

print(
    "FILES PER PROCESS"
)

print(
    "========================================"
)


for process, files in SAMPLES.items():

    print(
        f"{process:8s}: "
        f"{len(files)} files"
    )


# ============================================================
# XROOTD CONFIGURATION
# ============================================================

XRD_RETRIES = 4

XRD_TIMEOUT = 120

XRD_RETRY_SLEEP = 8


XRD_REDIRECTORS = [

    "root://cms-xrd-global.cern.ch//",

    "root://xrootd-cms.infn.it//",

    "root://cmsxrootd.fnal.gov//",

]


def candidate_urls(
    url
):

    """
    Return original URL plus alternative CMS redirectors.
    """

    out = [
        url
    ]

    if "/store/" in url:

        store_path = (
            "/store/"
            +
            url.split(
                "/store/",
                1
            )[1]
        )

        for redirector in XRD_REDIRECTORS:

            alt_url = (

                redirector

                +

                store_path.lstrip("/")

            )

            if alt_url not in out:

                out.append(
                    alt_url
                )

    return out


# ============================================================
# READ ROOT TREE WITH RETRIES
# ============================================================

def read_tree_arrays(
    url,
    tree_name,
    branches,
    library
):

    last_error = None

    for candidate in candidate_urls(
        url
    ):

        for trial in range(
            1,
            XRD_RETRIES + 1
        ):

            try:

                with uproot.open(
                    candidate,
                    timeout=XRD_TIMEOUT
                ) as root_file:

                    if tree_name not in root_file:

                        raise KeyError(
                            f"Tree '{tree_name}' not found"
                        )

                    tree = root_file[
                        tree_name
                    ]

                    return (
                        tree.arrays(
                            branches,
                            library=library
                        ),
                        candidate
                    )

            except KeyboardInterrupt:

                raise

            except Exception as exc:

                last_error = exc

                print(
                    f"[{tree_name}] "
                    f"attempt {trial}/{XRD_RETRIES} failed: "
                    f"{candidate}"
                )

                print(
                    type(exc).__name__,
                    exc
                )

                if (
                    trial
                    <
                    XRD_RETRIES
                ):

                    time.sleep(
                        XRD_RETRY_SLEEP
                        *
                        trial
                    )


    raise OSError(

        f"Could not read tree "
        f"'{tree_name}' after retries. "
        f"Original file: {url}. "
        f"Last error: {last_error}"

    )


# ============================================================
# SUM GEN WEIGHT
# ============================================================

def get_sum_genweight(
    files
):

    total_sumw = 0.0

    usable_files = []


    for i, f in enumerate(
        files,
        1
    ):

        print(
            f"Checking Runs "
            f"{i}/{len(files)}"
        )


        try:

            runs_arrays, working_url = (

                read_tree_arrays(

                    f,

                    "Runs",

                    [
                        "genEventSumw"
                    ],

                    library="np"

                )

            )


            file_sumw = float(

                np.sum(
                    runs_arrays[
                        "genEventSumw"
                    ]
                )

            )


            with uproot.open(
                working_url,
                timeout=XRD_TIMEOUT
            ) as root_file:

                if (
                    "Events"
                    not in root_file
                ):

                    raise KeyError(
                        "Events tree missing"
                    )


            total_sumw += file_sumw

            usable_files.append(
                working_url
            )


        except Exception as exc:

            print(
                "WARNING: skipping file"
            )

            print(
                f,
                exc
            )


    if len(
        usable_files
    ) == 0:

        raise RuntimeError(
            "No usable ROOT files"
        )


    if (
        not np.isfinite(
            total_sumw
        )
        or
        total_sumw == 0
    ):

        raise RuntimeError(
            f"Invalid sumGenWeight: "
            f"{total_sumw}"
        )


    print(
        "Usable files:",
        len(usable_files)
    )

    print(
        "sumGenWeight:",
        total_sumw
    )


    return (
        total_sumw,
        usable_files
    )


# ============================================================
# DNN MODEL
# ============================================================

class DNN(
    nn.Module
):

    def __init__(
        self
    ):

        super().__init__()


        self.add_module(
            "0",
            nn.Linear(
                13,
                225
            )
        )


        self.add_module(
            "1",
            nn.ReLU()
        )


        self.add_module(
            "2",
            nn.Dropout(
                0.1
            )
        )


        self.add_module(
            "3",
            nn.Linear(
                225,
                100
            )
        )


        self.add_module(
            "4",
            nn.ReLU()
        )


        self.add_module(
            "5",
            nn.Dropout(
                0.1
            )
        )


        self.add_module(
            "6",
            nn.Linear(
                100,
                64
            )
        )


        self.add_module(
            "7",
            nn.ReLU()
        )


        self.add_module(
            "8",
            nn.Dropout(
                0.1
            )
        )


        self.add_module(
            "9",
            nn.Linear(
                64,
                1
            )
        )


    def forward(
        self,
        x
    ):

        x = self._modules[
            "0"
        ](x)

        x = self._modules[
            "1"
        ](x)

        x = self._modules[
            "2"
        ](x)

        x = self._modules[
            "3"
        ](x)

        x = self._modules[
            "4"
        ](x)

        x = self._modules[
            "5"
        ](x)

        x = self._modules[
            "6"
        ](x)

        x = self._modules[
            "7"
        ](x)

        x = self._modules[
            "8"
        ](x)

        x = self._modules[
            "9"
        ](x)

        return x


# ============================================================
# LOAD MODEL
# ============================================================

print(
    "\nLoading DNN model..."
)


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


with open(
    THRESHOLD_FILE,
    "r"
) as f:

    threshold = json.load(
        f
    )[
        "threshold"
    ]


print(
    "DNN threshold =",
    threshold
)


# ============================================================
# HELPER
# ============================================================

def alt(
    arr,
    i,
    fill
):

    return ak.fill_none(

        ak.pad_none(

            arr,

            i + 1,

            clip=True

        )[:, i],

        fill

    )


# ============================================================
# CREATE FEATURES
# ============================================================

def make_features(
    data
):

    # --------------------------------------------------------
    # Muons
    # --------------------------------------------------------

    mu1pt = alt(
        data.Muon_pt,
        0,
        0.0
    )

    mu2pt = alt(
        data.Muon_pt,
        1,
        0.0
    )


    mu1eta = alt(
        data.Muon_eta,
        0,
        0.0
    )

    mu2eta = alt(
        data.Muon_eta,
        1,
        0.0
    )


    mu1phi = alt(
        data.Muon_phi,
        0,
        0.0
    )

    mu2phi = alt(
        data.Muon_phi,
        1,
        0.0
    )


    mu1q = alt(
        data.Muon_charge,
        0,
        0
    )

    mu2q = alt(
        data.Muon_charge,
        1,
        0
    )


    # --------------------------------------------------------
    # Dilepton mass
    # --------------------------------------------------------

    mll = np.sqrt(

        np.maximum(

            0,

            2
            *
            mu1pt
            *
            mu2pt
            *
            (

                np.cosh(
                    mu1eta
                    -
                    mu2eta
                )

                -

                np.cos(
                    mu1phi
                    -
                    mu2phi
                )

            )

        )

    )


    # --------------------------------------------------------
    # Delta phi
    # --------------------------------------------------------

    dphi = np.arccos(

        np.clip(

            np.cos(
                mu1phi
                -
                mu2phi
            ),

            -1,

            1

        )

    )


    # --------------------------------------------------------
    # Delta R
    # --------------------------------------------------------

    dR = np.sqrt(

        (
            mu1eta
            -
            mu2eta
        ) ** 2

        +

        dphi ** 2

    )


    # --------------------------------------------------------
    # Jets
    # --------------------------------------------------------

    j1pt = alt(
        data.Jet_pt,
        0,
        -999.0
    )

    j2pt = alt(
        data.Jet_pt,
        1,
        -999.0
    )


    j1eta = alt(
        data.Jet_eta,
        0,
        -999.0
    )

    j2eta = alt(
        data.Jet_eta,
        1,
        -999.0
    )


    j1phi = alt(
        data.Jet_phi,
        0,
        -999.0
    )

    j2phi = alt(
        data.Jet_phi,
        1,
        -999.0
    )


    # --------------------------------------------------------
    # mjj and delta eta jj
    # --------------------------------------------------------

    n_events = len(
        data
    )


    mjj = np.full(

        n_events,

        -999.0,

        dtype=np.float64

    )


    dEta_jj = np.full(

        n_events,

        -999.0,

        dtype=np.float64

    )


    valid = (

        (j1pt > 0)

        &

        (j2pt > 0)

    )


    valid_np = ak.to_numpy(
        valid
    )


    if np.any(
        valid_np
    ):

        j1pt_np = ak.to_numpy(
            j1pt
        )

        j2pt_np = ak.to_numpy(
            j2pt
        )

        j1eta_np = ak.to_numpy(
            j1eta
        )

        j2eta_np = ak.to_numpy(
            j2eta
        )

        j1phi_np = ak.to_numpy(
            j1phi
        )

        j2phi_np = ak.to_numpy(
            j2phi
        )


        deta = (

            j1eta_np
            -
            j2eta_np

        )


        dphi_jj = (

            j1phi_np
            -
            j2phi_np

        )


        mjj[valid_np] = np.sqrt(

            2
            *
            j1pt_np[valid_np]
            *
            j2pt_np[valid_np]
            *
            (

                np.cosh(
                    deta[valid_np]
                )

                -

                np.cos(
                    dphi_jj[valid_np]
                )

            )

        )


        dEta_jj[valid_np] = np.abs(

            deta[
                valid_np
            ]

        )


    # --------------------------------------------------------
    # Selection
    # --------------------------------------------------------

    sel = (

        (data.nMuon >= 2)

        &

        (data.nElectron == 0)

        &

        (mu1q != mu2q)

        &

        (mu1pt > 26)

        &

        (mu2pt > 20)

        &

        (abs(mu1eta) < 2.4)

        &

        (abs(mu2eta) < 2.4)

        &

        (mll > 110)

        &

        (mll < 150)

    )


    sel = ak.to_numpy(
        sel
    )


    # --------------------------------------------------------
    # DNN input
    # --------------------------------------------------------

    X = np.column_stack(

        [

            ak.to_numpy(
                mu1pt
            ),

            ak.to_numpy(
                mu2pt
            ),

            ak.to_numpy(
                mu1eta
            ),

            ak.to_numpy(
                mu2eta
            ),

            ak.to_numpy(
                dR
            ),

            ak.to_numpy(
                data.MET_pt
            ),

            np.sin(
                ak.to_numpy(
                    data.MET_phi
                )
            ),

            np.cos(
                ak.to_numpy(
                    data.MET_phi
                )
            ),

            ak.to_numpy(
                j1pt
            ),

            ak.to_numpy(
                j2pt
            ),

            dEta_jj,

            mjj,

            ak.to_numpy(
                data.nJet
            )

        ]

    )


    return (

        X[sel],

        mll[sel],

        sel

    )


# ============================================================
# DNN PREDICTION
# ============================================================

def predict_score(
    X
):

    if len(
        X
    ) == 0:

        return np.array(
            [],
            dtype=np.float64
        )


    X_scaled = scaler.transform(
        X
    )


    with torch.no_grad():

        logits = model(

            torch.tensor(

                X_scaled,

                dtype=torch.float32

            )

        )


        score = torch.sigmoid(
            logits
        ).numpy().ravel()


    return score


# ============================================================
# LOAD DATA
# ============================================================

print(
    "\n========================================"
)

print(
    "READING DATA"
)

print(
    "========================================"
)


data_branches = [

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

    "nJet",

]


data_arrays = []


with uproot.open(
    DATA_DNN_FILE
) as f:

    tree = f[
        "Events"
    ]

    data_score = tree[
        "DNN_score"
    ].array(
        library="np"
    )


print(
    "Data selected events:",
    len(data_score)
)


data_score = np.asarray(
    data_score,
    dtype=np.float64
)


data_score = data_score[
    np.isfinite(
        data_score
    )
]


# ============================================================
# LOAD MC SAMPLE
# ============================================================

MC_RESULTS = {}


def load_mc_sample(
    process,
    files
):

    print(
        "\n----------------------------------------"
    )

    print(
        f"Loading MC process: {process}"
    )

    print(
        "----------------------------------------"
    )


    if len(
        files
    ) == 0:

        print(
            f"WARNING: no files for {process}"
        )

        return {

            "score": np.array(
                []
            ),

            "weight": np.array(
                []
            )

        }


    # --------------------------------------------------------
    # Get sum gen weight
    # --------------------------------------------------------

    sumw, usable_files = (

        get_sum_genweight(
            files
        )

    )


    normalization = (

        XSECS_PB[
            process
        ]

        *

        LUMI_PB

        /

        sumw

    )


    print(
        "Normalization =",
        normalization
    )


    # --------------------------------------------------------
    # Read Events
    # --------------------------------------------------------

    arrays = []

    loaded_sumw = 0.0


    branches = [

        "genWeight",

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

        "nJet",

    ]


    for i, f in enumerate(

        usable_files,

        1

    ):

        print(

            f"Reading {process}: "
            f"{i}/{len(usable_files)}"

        )


        try:

            arr, working_url = (

                read_tree_arrays(

                    f,

                    "Events",

                    branches,

                    library="ak"

                )

            )


            # Read sumw for exactly
            # the successfully loaded file

            runs_arrays, _ = (

                read_tree_arrays(

                    working_url,

                    "Runs",

                    [
                        "genEventSumw"
                    ],

                    library="np"

                )

            )


            loaded_sumw += float(

                np.sum(

                    runs_arrays[
                        "genEventSumw"
                    ]

                )

            )


            arrays.append(
                arr
            )


        except Exception as exc:

            print(
                "WARNING: skipping Events file"
            )

            print(
                f,
                exc
            )


    if len(
        arrays
    ) == 0:

        return {

            "score": np.array(
                []
            ),

            "weight": np.array(
                []
            )

        }


    # --------------------------------------------------------
    # Recalculate normalization using
    # successfully loaded files
    # --------------------------------------------------------

    normalization = (

        XSECS_PB[
            process
        ]

        *

        LUMI_PB

        /

        loaded_sumw

    )


    data = ak.concatenate(
        arrays
    )


    # --------------------------------------------------------
    # genWeight
    # --------------------------------------------------------

    gen_weight = ak.to_numpy(

        data.genWeight

    ).astype(
        np.float64
    )


    physics_weight = (

        gen_weight

        *

        normalization

    )


    # --------------------------------------------------------
    # Features and selection
    # --------------------------------------------------------

    X, mll, sel = make_features(
        data
    )


    # Apply same selection to genWeight

    physics_weight = (

        physics_weight[
            sel
        ]

    )


    # --------------------------------------------------------
    # Remove non-finite events
    # --------------------------------------------------------

    finite_mask = (

        np.isfinite(
            X
        ).all(
            axis=1
        )

        &

        np.isfinite(
            physics_weight
        )

        &

        np.isfinite(
            mll
        )

    )


    X = X[
        finite_mask
    ]


    physics_weight = physics_weight[
        finite_mask
    ]


    mll = mll[
        finite_mask
    ]


    # --------------------------------------------------------
    # DNN score
    # --------------------------------------------------------

    score = predict_score(
        X
    )


    print(
        f"{process}: "
        f"{len(score)} selected events"
    )


    print(
        f"{process}: "
        f"expected yield = "
        f"{np.sum(physics_weight):.4f}"
    )


    return {

        "score": score,

        "weight": physics_weight,

        "mll": mll

    }


# ============================================================
# RUN MC
# ============================================================

for process in SAMPLES:

    MC_RESULTS[
        process
    ] = load_mc_sample(

        process,

        SAMPLES[
            process
        ]

    )


# ============================================================
# BUILD HISTOGRAMS
# ============================================================

print(
    "\n========================================"
)

print(
    "BUILDING HISTOGRAMS"
)

print(
    "========================================"
)


def hist_weighted(
    score,
    weight
):

    h, _ = np.histogram(

        score,

        bins=BINS,

        weights=weight

    )

    h2, _ = np.histogram(

        score,

        bins=BINS,

        weights=weight ** 2

    )

    return (

        h,

        np.sqrt(
            h2
        )

    )


def hist_data(
    score
):

    h, _ = np.histogram(

        score,

        bins=BINS

    )

    return h


# ------------------------------------------------------------
# Data
# ------------------------------------------------------------

data_hist = hist_data(
    data_score
)


data_error = np.sqrt(
    data_hist
)


# ------------------------------------------------------------
# MC
# ------------------------------------------------------------

mc_hists = {}

mc_errors = {}


for process in MC_RESULTS:

    score = MC_RESULTS[
        process
    ][
        "score"
    ]

    weight = MC_RESULTS[
        process
    ][
        "weight"
    ]


    h, e = hist_weighted(

        score,

        weight

    )


    mc_hists[
        process
    ] = h


    mc_errors[
        process
    ] = e


# ============================================================
# BACKGROUND SUM
# ============================================================

background_processes = [

    "DY",

    "EWK",

    "TT",

    "WZ",

    "ZZ",

]


total_background = np.zeros(
    N_BINS
)


total_background_variance = np.zeros(
    N_BINS
)


for process in background_processes:

    total_background += (

        mc_hists[
            process
        ]

    )


    total_background_variance += (

        mc_errors[
            process
        ] ** 2

    )


total_background_error = np.sqrt(

    total_background_variance

)


# ============================================================
# SIGNAL
# ============================================================

ggh_hist = mc_hists[
    "ggH"
]

vbf_hist = mc_hists[
    "VBF"
]


# ============================================================
# PLOT
# ============================================================

print(
    "\n========================================"
)

print(
    "MAKING PLOT"
)

print(
    "========================================"
)


fig, (
    ax,
    ax_ratio
) = plt.subplots(

    2,

    1,

    figsize=(
        9,
        9
    ),

    gridspec_kw={

        "height_ratios": [
            3.5,
            1
        ],

        "hspace": 0.05

    },

    sharex=True

)


# ============================================================
# TOP PANEL
# ============================================================


# ------------------------------------------------------------
# Stacked background
# ------------------------------------------------------------

bottom = np.zeros(
    N_BINS
)


# Use CMS-like ordering:
# Other bkg -> ZZ/WZ -> EWK -> TT -> DY

stack_order = [

    (
        "ZZ",
        "Other bkg.",
    ),

    (
        "WZ",
        "Diboson",
    ),

    (
        "EWK",
        "EWK",
    ),

    (
        "TT",
        "Top quark",
    ),

    (
        "DY",
        "DY",
    ),

]


stack_colors = [

    "lightgray",

    "deepskyblue",

    "magenta",

    "limegreen",

    "orange",

]


for (
    process,
    label
), color in zip(

    stack_order,

    stack_colors

):

    ax.bar(

        BIN_CENTERS,

        mc_hists[
            process
        ],

        width=BIN_WIDTH,

        bottom=bottom,

        align="center",

        color=color,

        edgecolor="black",

        linewidth=0.4,

        label=label

    )


    bottom += mc_hists[
        process
    ]


# ------------------------------------------------------------
# MC uncertainty band
# ------------------------------------------------------------

total_mc = bottom


mc_unc_low = (

    total_mc

    -

    total_background_error

)


mc_unc_high = (

    total_mc

    +

    total_background_error

)


ax.fill_between(

    BIN_CENTERS,

    mc_unc_low,

    mc_unc_high,

    step="mid",

    alpha=0.30,

    color="gray",

    linewidth=0,

    label="MC stat. unc."

)


# ------------------------------------------------------------
# Signal lines
# ------------------------------------------------------------

ax.step(

    BINS,

    np.r_[

        ggh_hist,

        ggh_hist[-1]

    ],

    where="post",

    color="black",

    linewidth=2.0,

    label="ggH"

)


ax.step(

    BINS,

    np.r_[

        vbf_hist,

        vbf_hist[-1]

    ],

    where="post",

    color="red",

    linewidth=2.0,

    label="VBF"

)


# ------------------------------------------------------------
# Data
# ------------------------------------------------------------

ax.errorbar(

    BIN_CENTERS,

    data_hist,

    yerr=data_error,

    fmt="o",

    color="black",

    markersize=5,

    markeredgewidth=1,

    linewidth=1,

    capsize=0,

    label="Data",

    zorder=10

)


# ============================================================
# AXIS
# ============================================================

ax.set_yscale(
    "log"
)


ax.set_ylabel(
    "Events / "
    f"{BIN_WIDTH:.2f}",
    fontsize=14
)


ax.set_xlim(
    XMIN,
    XMAX
)


ax.set_ylim(
    0.1,
    None
)


ax.tick_params(

    axis="both",

    which="both",

    direction="in",

    top=True,

    right=True,

    labelsize=12

)


ax.xaxis.set_minor_locator(
    AutoMinorLocator()
)


ax.yaxis.set_minor_locator(
    AutoMinorLocator()
)


# ============================================================
# CMS LABEL
# ============================================================

ax.text(

    0.05,

    0.93,

    "CMS",

    transform=ax.transAxes,

    fontsize=24,

    fontweight="bold",

    va="top"

)


ax.text(

    0.95,

    1.04,

    r"$7.98\ \mathrm{fb}^{-1}\ (13.6\ \mathrm{TeV})$",

    transform=ax.transAxes,

    fontsize=14,

    ha="right",

    va="bottom"

)


# ============================================================
# LEGEND
# ============================================================

handles, labels = ax.get_legend_handles_labels()


ax.legend(

    handles,

    labels,

    loc="upper center",

    bbox_to_anchor=(

        0.60,

        0.99

    ),

    ncol=2,

    fontsize=10,

    frameon=False

)


# ============================================================
# RATIO PANEL
# ============================================================

ratio = np.divide(

    data_hist,

    total_mc,

    out=np.zeros_like(

        data_hist,

        dtype=float

    ),

    where=total_mc != 0

)


ratio_error = np.divide(

    data_error,

    total_mc,

    out=np.zeros_like(

        data_error,

        dtype=float

    ),

    where=total_mc != 0

)


# MC relative uncertainty

mc_rel_error = np.divide(

    total_background_error,

    total_mc,

    out=np.zeros_like(

        total_background_error,

        dtype=float

    ),

    where=total_mc != 0

)


ax_ratio.axhline(

    1.0,

    color="black",

    linewidth=1,

    linestyle="--"

)


ax_ratio.fill_between(

    BIN_CENTERS,

    1.0 - mc_rel_error,

    1.0 + mc_rel_error,

    step="mid",

    alpha=0.30,

    color="gray"

)


ax_ratio.errorbar(

    BIN_CENTERS,

    ratio,

    yerr=ratio_error,

    fmt="o",

    color="black",

    markersize=5,

    capsize=0,

    zorder=10

)


# ============================================================
# RATIO AXIS
# ============================================================

ax_ratio.set_ylabel(

    "Data / MC",

    fontsize=13

)


ax_ratio.set_xlabel(

    "DNN output",

    fontsize=15

)


ax_ratio.set_ylim(

    0.5,

    1.5

)


ax_ratio.set_xlim(

    XMIN,

    XMAX

)


ax_ratio.set_yticks(

    [
        0.6,
        0.8,
        1.0,
        1.2,
        1.4
    ]

)


ax_ratio.tick_params(

    axis="both",

    which="both",

    direction="in",

    top=True,

    right=True,

    labelsize=12

)


ax_ratio.xaxis.set_minor_locator(

    AutoMinorLocator()

)


# ============================================================
# GRID
# ============================================================

ax.grid(

    axis="y",

    which="major",

    linestyle=":",

    alpha=0.4

)


ax_ratio.grid(

    axis="y",

    which="major",

    linestyle=":",

    alpha=0.4

)


# ============================================================
# FINAL LAYOUT
# ============================================================

plt.tight_layout()


# ============================================================
# SAVE
# ============================================================

print(
    "\nSaving:"
)


print(
    OUTPUT_PNG
)


plt.savefig(

    OUTPUT_PNG,

    dpi=300,

    bbox_inches="tight"

)


print(
    OUTPUT_PDF
)


plt.savefig(

    OUTPUT_PDF,

    bbox_inches="tight"

)


plt.close()


print(
    "\n========================================"
)

print(
    "DONE"
)

print(
    "========================================"
)
