#!/usr/bin/env python3
import os
import math
import numpy as np
import uproot
import awkward as ak
import ROOT
import json

ROOT.gErrorIgnoreLevel = ROOT.kFatal  # quiet ROOT
print('uproot:', uproot.__version__, '| awkward:', ak.__version__)
print("User Tet 013:")
"""
ALTERNATIVE RUN IN VEV
----------------------
python3 -m venv $HOME/venvs/uproot_env
source $HOME/venvs/uproot_env/bin/activate
pip install --upgrade pip
pip install uproot awkward numpy matplotlib
pip install --force-reinstall --ignore-installed numpy==2.0.2
python -c "import uproot; print(uproot.__version__)"
python -c "import uproot, awkward as ak; print('uproot:', uproot.__version__, '| awkward:', ak.__version__)"
deactivate
"""

# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------
def make_log_bins(xmin, xmax, nbins):
    """Log-spaced bin edges (like C++ make_log_bins)."""
    xmin = max(xmin, 1e-6)
    return np.logspace(np.log10(xmin), np.log10(xmax), nbins + 1)


def basename(path: str) -> str:
    return os.path.basename(path)


def extract_meta_topkey(path: str) -> str:
    """
    BKG: DSID f.ex. "601190".
    SIG: generator tag f.ex. "MGPy8EG_LRSM_WR_HNL_5000_500_2L2J_MadSpin".
    """
    tag = "TeV."
    pos = path.find(tag)
    if pos == -1:
        return ""
    pos += len(tag)
    pos2 = path.find(".", pos)
    if pos2 == -1:
        return ""
    return path[pos:pos2]


def extract_meta_filekey(path: str, top_block: dict) -> str:
    base = basename(path)
    if isinstance(top_block, dict) and base in top_block:
        return base
    if isinstance(top_block, dict) and len(top_block) == 1:
        return next(iter(top_block.keys()))
    return base


def get_total_sumw_from_metadata(sample, meta) -> float:
    totalW = 0.0
    for fullPath in sample["files"]:
        topkey = extract_meta_topkey(fullPath)
        if not topkey:
            print(f"WARNING: could not extract metadata topkey from path: {fullPath}")
            continue
        if topkey not in meta:
            print(f"WARNING: topkey {topkey} not found in metadata")
            continue
        top_block = meta[topkey]
        fileKey = extract_meta_filekey(fullPath, top_block)
        if not isinstance(top_block, dict) or fileKey not in top_block:
            print(f"WARNING: file key {fileKey} not found in metadata for topkey {topkey}")
            continue
        fileBlock = top_block[fileKey]
        if "AllExecutedEvents" not in fileBlock:
            print(f"WARNING: AllExecutedEvents not found for {fileKey} (topkey {topkey})")
            continue
        allExec = fileBlock["AllExecutedEvents"]
        if "W" not in allExec:
            print(f"WARNING: AllExecutedEvents has no W for {fileKey} (topkey {topkey})")
            continue
        W = float(allExec["W"])
        totalW += W
    return totalW


# ------------------------------------------------------------
# Global selections (mirroring C++)
# ------------------------------------------------------------
lep_pt_min  = 10.0  # GeV
lep_eta_max = 2.5
jet_pt_min  = 10.0  # GeV
jet_eta_max = 2.5
jvt_min = 0.5
jWidth_cut = 0.05
jNumTrkPt500_cut = 10.0

# isolation
lep_rel_ptvarcone20_max  = 0.15
lep_rel_topoetcone20_max = 0.15

# ------------------------------------------------------------
# Metadata
# ------------------------------------------------------------
use_metadata_totalW = True
metadata_json_path = "/home/alexanum/WORKSPACE/testing/MadGraph5/analysis/SIGvsBKG_RECO/metadataLLP1_backgrounds.json"
metadata = {}
if use_metadata_totalW:
    try:
        with open(metadata_json_path, "r") as jf:
            metadata = json.load(jf)
    except OSError as e:
        print(f"ERROR: cannot open metadata JSON: {metadata_json_path} ({e})")
        use_metadata_totalW = False

# ------------------------------------------------------------
# Sample definition (mirrors C++ list)
# ------------------------------------------------------------
#Define samples
'''
MCprofile: CrossSection[pb], Int.Lum.[pb^-], FilterEff, k-factor
Signal: CS(50GeV, 5TeV, 0.1)=4.072e-4pb 
        Int. Lum: 29.0493 fb^- = 29049.3 pb^-
        Filter Eff: MG (#tau events: 6%) 
Signal: CS(70GeV, 5TeV, 0.1)=3.968-4pb  
Signal: CS(500GeV, 5TeV, 0.1)=2.32479e-4pb  
Signal: CS(500GeV, 10TeV, 0.1)=4.672e-6pb     (e-15mb)
Signal: CS(50GeV, 10TeV, 0.1)=1.303-5pb  
'''
samples = [
    dict(
        name="WR5000_N500_2L2J",  # match your C++ signal 500 GeV, 5 TeV example
        type="SIG",
        sigma_pb=1.303e-5,
        lumi_pb=500000.0,  # pb^-1
        filter_eff=0.94,
        k_factor=1.0,
        files=[
            "/storage/alexanum/DAOD_LLP_WORKFLOW/ROOT_files/DAOD_LLP1.DAOD_LLP1.mc23_13p6TeV.MGPy8EG_LRSM_WR_HNL_10000_50_2L2J_MadSpin.DAOD_LLP1.pool.root"
        ],
    ),
    dict(
        name="Zmumu",
        type="BKG_Z",
        sigma_pb=1998.8,
        lumi_pb=500000.0,
        filter_eff=1.0,
        k_factor=1.0,
        files=[
            "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.601190.PhPy8EG_AZNLO_Zmumu.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762661._000102.pool.root.1",
            "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.601190.PhPy8EG_AZNLO_Zmumu.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762661._000278.pool.root.1",
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.601190.PhPy8EG_AZNLO_Zmumu.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762668._000497.pool.root.1",
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.601190.PhPy8EG_AZNLO_Zmumu.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762671._000561.pool.root.1",
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.601190.PhPy8EG_AZNLO_Zmumu.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762661._000142.pool.root.1",
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.601190.PhPy8EG_AZNLO_Zmumu.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762666._000115.pool.root.1",
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.601190.PhPy8EG_AZNLO_Zmumu.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762668._000525.pool.root.1",
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.601190.PhPy8EG_AZNLO_Zmumu.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762661._000157.pool.root.1",
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.601190.PhPy8EG_AZNLO_Zmumu.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762666._000254.pool.root.1",
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.601190.PhPy8EG_AZNLO_Zmumu.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762671._000421.pool.root.1",
        ],
    ),
    dict(
        name="Zmumu_CVBV",
        type="BKG_Z",
        sigma_pb=2335.9,
        lumi_pb=500000.0,
        filter_eff=8.444979e-1,
        k_factor=0.9332,
        files=[
            "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700791.Sh_2214_Zmumu_maxHTpTV2_CVetoBVeto.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762549._000584.pool.root.1",
            "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700791.Sh_2214_Zmumu_maxHTpTV2_CVetoBVeto.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762554._000186.pool.root.1",  
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700791.Sh_2214_Zmumu_maxHTpTV2_CVetoBVeto.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762557._000423.pool.root.1",  
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700791.Sh_2214_Zmumu_maxHTpTV2_CVetoBVeto.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762562._000039.pool.root.1",
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700791.Sh_2214_Zmumu_maxHTpTV2_CVetoBVeto.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762552._000255.pool.root.1",  
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700791.Sh_2214_Zmumu_maxHTpTV2_CVetoBVeto.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762554._000310.pool.root.1",  
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700791.Sh_2214_Zmumu_maxHTpTV2_CVetoBVeto.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762557._000664.pool.root.1",
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700791.Sh_2214_Zmumu_maxHTpTV2_CVetoBVeto.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762552._000434.pool.root.1",  
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700791.Sh_2214_Zmumu_maxHTpTV2_CVetoBVeto.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762554._000458.pool.root.1",  
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700791.Sh_2214_Zmumu_maxHTpTV2_CVetoBVeto.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762559._000104.pool.root.1",
        ],
    ),
    dict(
        name="Zmumu_BF",
        type="BKG_Z",
        sigma_pb=2336.1,
        lumi_pb=500000.0,
        filter_eff=2.502305e-2,
        k_factor=0.9332,
        files=[
            "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700789.Sh_2214_Zmumu_maxHTpTV2_BFilter.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762363._000036.pool.root.1",  
            "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700789.Sh_2214_Zmumu_maxHTpTV2_BFilter.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762365._000150.pool.root.1",  
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700789.Sh_2214_Zmumu_maxHTpTV2_BFilter.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762365._000303.pool.root.1",  
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700789.Sh_2214_Zmumu_maxHTpTV2_BFilter.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762365._000351.pool.root.1",
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700789.Sh_2214_Zmumu_maxHTpTV2_BFilter.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762365._000068.pool.root.1",  
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700789.Sh_2214_Zmumu_maxHTpTV2_BFilter.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762365._000158.pool.root.1",  
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700789.Sh_2214_Zmumu_maxHTpTV2_BFilter.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762365._000304.pool.root.1",
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700789.Sh_2214_Zmumu_maxHTpTV2_BFilter.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762365._000071.pool.root.1",  
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700789.Sh_2214_Zmumu_maxHTpTV2_BFilter.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762365._000221.pool.root.1",  
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700789.Sh_2214_Zmumu_maxHTpTV2_BFilter.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762365._000327.pool.root.1",
        ],
    ),
    dict(
        name="Zmumu_CFBV",
        type="BKG_Z",
        sigma_pb=2336.1,
        lumi_pb=500000.0,
        filter_eff=1.304099e-1,
        k_factor=0.9332,
        files=[
            "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700790.Sh_2214_Zmumu_maxHTpTV2_CFilterBVeto.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762751._000027.pool.root.1",  
            "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700790.Sh_2214_Zmumu_maxHTpTV2_CFilterBVeto.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762751._000256.pool.root.1",  
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700790.Sh_2214_Zmumu_maxHTpTV2_CFilterBVeto.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762751._000315.pool.root.1", 
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700790.Sh_2214_Zmumu_maxHTpTV2_CFilterBVeto.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762751._000541.pool.root.1",
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700790.Sh_2214_Zmumu_maxHTpTV2_CFilterBVeto.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762751._000167.pool.root.1",  
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700790.Sh_2214_Zmumu_maxHTpTV2_CFilterBVeto.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762751._000263.pool.root.1",  
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700790.Sh_2214_Zmumu_maxHTpTV2_CFilterBVeto.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762751._000418.pool.root.1",
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700790.Sh_2214_Zmumu_maxHTpTV2_CFilterBVeto.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762751._000236.pool.root.1",  
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700790.Sh_2214_Zmumu_maxHTpTV2_CFilterBVeto.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762751._000306.pool.root.1",  
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700790.Sh_2214_Zmumu_maxHTpTV2_CFilterBVeto.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762751._000489.pool.root.1",
        ],
    ),
    dict(
        name="Zee_CVBV",
        type="BKG_Z",
        sigma_pb=2336.1,
        lumi_pb=500000.0,
        filter_eff=8.442598e-1,
        k_factor=0.9332,
        files=[
            "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700788.Sh_2214_Zee_maxHTpTV2_CVetoBVeto.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762860._000141.pool.root.1",
            "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700788.Sh_2214_Zee_maxHTpTV2_CVetoBVeto.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762860._000953.pool.root.1",
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700788.Sh_2214_Zee_maxHTpTV2_CVetoBVeto.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762867._000327.pool.root.1", 
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700788.Sh_2214_Zee_maxHTpTV2_CVetoBVeto.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762872._000793.pool.root.1",
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700788.Sh_2214_Zee_maxHTpTV2_CVetoBVeto.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762860._000145.pool.root.1",  
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700788.Sh_2214_Zee_maxHTpTV2_CVetoBVeto.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762860._001124.pool.root.1",  
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700788.Sh_2214_Zee_maxHTpTV2_CVetoBVeto.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762870._000286.pool.root.1",
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700788.Sh_2214_Zee_maxHTpTV2_CVetoBVeto.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762860._000186.pool.root.1",  
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700788.Sh_2214_Zee_maxHTpTV2_CVetoBVeto.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762865._000087.pool.root.1",  
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700788.Sh_2214_Zee_maxHTpTV2_CVetoBVeto.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762872._000196.pool.root.1",
        ],
    ),
    dict(
        name="Zee_BF",
        type="BKG_Z",
        sigma_pb=2336.1,
        lumi_pb=500000.0,
        filter_eff=2.558475e-2,
        k_factor=0.9332,
        files=[
            "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700786.Sh_2214_Zee_maxHTpTV2_BFilter.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47656889._000048.pool.root.1", 
            "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700786.Sh_2214_Zee_maxHTpTV2_BFilter.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47656889._000091.pool.root.1",  
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700786.Sh_2214_Zee_maxHTpTV2_BFilter.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47656889._000203.pool.root.1",  
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700786.Sh_2214_Zee_maxHTpTV2_BFilter.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47656889._000345.pool.root.1",
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700786.Sh_2214_Zee_maxHTpTV2_BFilter.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47656889._000059.pool.root.1",  
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700786.Sh_2214_Zee_maxHTpTV2_BFilter.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47656889._000143.pool.root.1",  
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700786.Sh_2214_Zee_maxHTpTV2_BFilter.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47656889._000212.pool.root.1",
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700786.Sh_2214_Zee_maxHTpTV2_BFilter.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47656889._000073.pool.root.1",  
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700786.Sh_2214_Zee_maxHTpTV2_BFilter.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47656889._000185.pool.root.1",  
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700786.Sh_2214_Zee_maxHTpTV2_BFilter.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47656889._000337.pool.root.1",
        ],
    ),
    dict(
        name="Zee_CFBV",
        type="BKG_Z",
        sigma_pb=2336.1,
        lumi_pb=500000.0,
        filter_eff=1.301249e-1,
        k_factor=0.9332,
        files=[
            "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700787.Sh_2214_Zee_maxHTpTV2_CFilterBVeto.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47658627._000207.pool.root.1",  
            "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700787.Sh_2214_Zee_maxHTpTV2_CFilterBVeto.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47658627._000294.pool.root.1",  
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700787.Sh_2214_Zee_maxHTpTV2_CFilterBVeto.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47658627._000451.pool.root.1",  
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700787.Sh_2214_Zee_maxHTpTV2_CFilterBVeto.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47658627._000601.pool.root.1",
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700787.Sh_2214_Zee_maxHTpTV2_CFilterBVeto.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47658627._000246.pool.root.1",  
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700787.Sh_2214_Zee_maxHTpTV2_CFilterBVeto.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47658627._000331.pool.root.1",  
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700787.Sh_2214_Zee_maxHTpTV2_CFilterBVeto.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47658627._000563.pool.root.1",
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700787.Sh_2214_Zee_maxHTpTV2_CFilterBVeto.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47658627._000274.pool.root.1",  
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700787.Sh_2214_Zee_maxHTpTV2_CFilterBVeto.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47658627._000392.pool.root.1",  
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700787.Sh_2214_Zee_maxHTpTV2_CFilterBVeto.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47658627._000572.pool.root.1",
        ],
    ),
    dict(
        name="ttbar_llqq",
        type="BKG_T",
        sigma_pb=811.29,
        lumi_pb=500000.0,
        filter_eff=4.384566e-1,
        k_factor=1.138433852,
        files=[
            "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.601229.PhPy8EG_A14_ttbar_hdamp258p75_SingleLep.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47659319._000078.pool.root.1",  
            "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.601229.PhPy8EG_A14_ttbar_hdamp258p75_SingleLep.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47659319._000245.pool.root.1",  
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.601229.PhPy8EG_A14_ttbar_hdamp258p75_SingleLep.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47659319._000814.pool.root.1",  
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.601229.PhPy8EG_A14_ttbar_hdamp258p75_SingleLep.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47659319._000865.pool.root.1",
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.601229.PhPy8EG_A14_ttbar_hdamp258p75_SingleLep.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47659319._000137.pool.root.1",  
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.601229.PhPy8EG_A14_ttbar_hdamp258p75_SingleLep.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47659319._000682.pool.root.1",  
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.601229.PhPy8EG_A14_ttbar_hdamp258p75_SingleLep.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47659319._000848.pool.root.1",
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.601229.PhPy8EG_A14_ttbar_hdamp258p75_SingleLep.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47659319._000187.pool.root.1",  
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.601229.PhPy8EG_A14_ttbar_hdamp258p75_SingleLep.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47659319._000723.pool.root.1",  
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.601229.PhPy8EG_A14_ttbar_hdamp258p75_SingleLep.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47659319._000850.pool.root.1",
        ],
    ),
    dict(
        name="ttbar_llll",
        type="BKG_T",
        sigma_pb=85.482,
        lumi_pb=500000.0,
        filter_eff=1.0,
        k_factor=1.138433852,
        files=[
            "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.601230.PhPy8EG_A14_ttbar_hdamp258p75_dil.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762465._000196.pool.root.1",  
            "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.601230.PhPy8EG_A14_ttbar_hdamp258p75_dil.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762465._000356.pool.root.1",  
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.601230.PhPy8EG_A14_ttbar_hdamp258p75_dil.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762465._000667.pool.root.1",  
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.601230.PhPy8EG_A14_ttbar_hdamp258p75_dil.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762465._000872.pool.root.1",
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.601230.PhPy8EG_A14_ttbar_hdamp258p75_dil.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762465._000331.pool.root.1",  
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.601230.PhPy8EG_A14_ttbar_hdamp258p75_dil.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762465._000497.pool.root.1",  
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.601230.PhPy8EG_A14_ttbar_hdamp258p75_dil.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762465._000789.pool.root.1",
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.601230.PhPy8EG_A14_ttbar_hdamp258p75_dil.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762465._000342.pool.root.1",  
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.601230.PhPy8EG_A14_ttbar_hdamp258p75_dil.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762465._000624.pool.root.1",  
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.601230.PhPy8EG_A14_ttbar_hdamp258p75_dil.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762465._000819.pool.root.1",
        ],
    ),
    dict(
        name="ttbar_qqqq",
        type="BKG_T",
        sigma_pb=811.29,
        lumi_pb=3549.0,
        filter_eff=4.561727e-1,
        k_factor=1.1384339,
        files=[
            "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.601237.PhPy8EG_A14_ttbar_hdamp258p75_allhad.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47657677._000301.pool.root.1",  
            "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.601237.PhPy8EG_A14_ttbar_hdamp258p75_allhad.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47657677._000475.pool.root.1",  
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.601237.PhPy8EG_A14_ttbar_hdamp258p75_allhad.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47657677._000570.pool.root.1",  
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.601237.PhPy8EG_A14_ttbar_hdamp258p75_allhad.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47657677._000741.pool.root.1",
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.601237.PhPy8EG_A14_ttbar_hdamp258p75_allhad.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47657677._000424.pool.root.1",  
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.601237.PhPy8EG_A14_ttbar_hdamp258p75_allhad.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47657677._000506.pool.root.1",  
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.601237.PhPy8EG_A14_ttbar_hdamp258p75_allhad.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47657677._000584.pool.root.1",
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.601237.PhPy8EG_A14_ttbar_hdamp258p75_allhad.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47657677._000427.pool.root.1",  
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.601237.PhPy8EG_A14_ttbar_hdamp258p75_allhad.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47657677._000532.pool.root.1",  
            # "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.601237.PhPy8EG_A14_ttbar_hdamp258p75_allhad.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47657677._000595.pool.root.1",
        ],
    ),
]

# Group definitions
group_definitions = {
    "SIG":   ["WR5000_N500_2L2J"],
    "Zmumu": ["Zmumu", "Zmumu_CVBV", "Zmumu_BF", "Zmumu_CFBV"],
    "Zee":   ["Zee_CVBV", "Zee_BF", "Zee_CFBV"],
    "ttbar": ["ttbar_llqq", "ttbar_llll", "ttbar_qqqq"],
}

output_dir = "/home/alexanum/WORKSPACE/testing/MadGraph5/analysis/SIGvsBKG_RECO/Figures/SIGvsBKG_recoDV_500fb"
os.makedirs(output_dir, exist_ok=True)

# ------------------------------------------------------------
# Histogram binning (match C++ as close as practical)
# ------------------------------------------------------------
lep1_bins = make_log_bins(lep_pt_min, 4350.0, 20)   #3950
lep2_bins = make_log_bins(lep_pt_min, 2950.0, 15)
jet_bins  = make_log_bins(jet_pt_min, 3350.0, 20)
disp_bins = make_log_bins(0.01, 20.0, 20)
m_lljj_bins = make_log_bins(95.0, 1950.0, 20)

# ------------------------------------------------------------
# Histogram helpers
# ------------------------------------------------------------
def make_hist_1d(name, title, nbins, edges=None, xmin=None, xmax=None):
    if edges is not None:
        arr = np.array(edges, dtype="float64")
        h = ROOT.TH1F(name, title, len(arr) - 1, arr)
    else:
        h = ROOT.TH1F(name, title, nbins, xmin, xmax)
    ROOT.TH1.SetDefaultSumw2(True)
    return h


def book_sample_histograms(sname, stype):
    """Return dict of histograms for this sample."""
    y_title = "event weight: #sigma #scale[0.7]{#int} #font[12]{L} #varepsilon k w^{i} / W_{tot}"

    histos = {}
    histos["lep1_pt"] = make_hist_1d(
        f"h_{sname}_lep1_pt",
        f" ;pT(#it{{l_{{1}}}} ) / GeV;{y_title}",
        nbins=len(lep1_bins) - 1,
        edges=lep1_bins,
    )
    histos["lep2_pt"] = make_hist_1d(
        f"h_{sname}_lep2_pt",
        f" ;pT(#it{{l_{{2}}}} ) / GeV;{y_title}",
        nbins=len(lep2_bins) - 1,
        edges=lep2_bins,
    )
    histos["jet_pt"] = make_hist_1d(
        f"h_{sname}_jet_pt",
        f" ;pT(#it{{j}} ) / GeV;{y_title}",
        nbins=len(jet_bins) - 1,
        edges=jet_bins,
    )
    histos["disp"] = make_hist_1d(
        f"h_{sname}_disp",
        f" ;max |d_{0}(#it{{l}} )| / mm;{y_title}",
        nbins=len(disp_bins) - 1,
        edges=disp_bins,
    )
    # dilepton charge
    edges_q = np.array([-0.75, -0.25, 0.25, 0.75], dtype="float64")
    histos["dilep_q"] = ROOT.TH1F(
        f"h_{sname}_dilep_q",
        f" ;;{y_title}",
        3,
        edges_q,
    )
    ROOT.TH1.SetDefaultSumw2(True)

    histos["dr_lep1_jj"] = make_hist_1d(
        f"h_{sname}_dr_lep1_jj",
        f" ;#DeltaR(#it{{l_{{1}}, jj}} );{y_title}",
        20,
        xmin=0.0,
        xmax=5.5,
    )
    histos["dr_lep2_jj"] = make_hist_1d(
        f"h_{sname}_dr_lep2_jj",
        f" ;#DeltaR(#it{{l_{{2}}, jj}} );{y_title}",
        20,
        xmin=0.0,
        xmax=5.5,
    )
    histos["m_l2jj"] = make_hist_1d(
        f"h_{sname}_m_l2jj",
        f" ;m(#it{{l_{{2}}jj}}) / GeV;{y_title}",
        25,
        xmin=20.0,
        xmax=145.0,
    )
    histos["m_lljj"] = make_hist_1d(
        f"h_{sname}_m_lljj",
        f" ;m(#it{{lljj}}) / GeV;{y_title}",
        nbins=len(m_lljj_bins) - 1,
        edges=m_lljj_bins,
    )
    histos["m_jj"] = make_hist_1d(
        f"h_{sname}_m_jj",
        f" ;m(#it{{jj}}) / GeV;{y_title}",
        20,
        xmin=15.0,
        xmax=195.0,
    )

    # if stype == "SIG":
    #     histos["disp_reco_d0"] = make_hist_1d(
    #         f"h_{sname}_disp_reco_d0",
    #         f" ;L_{{reco}}^{{d0}} / mm;{y_title}",
    #         nbins=len(disp_bins) - 1,
    #         edges=disp_bins,
    #     )
    #     histos["disp_truth_Lxyz"] = make_hist_1d(
    #         f"h_{sname}_disp_truth_Lxyz",
    #         f" ;L_{{truth}}^{{xyz}} / mm;{y_title}",
    #         nbins=len(disp_bins) - 1,
    #         edges=disp_bins,
    #     )

    return histos


def fill_reco_variables(
    sname,
    leptons,
    jets_sel,
    w_event,
    buff_weight,
    histos
):
    """
    Python translation of the C++ fill_reco_variables, with an extra
    third criterion based on plain d0 (no PV correction):

      1) log(pT) difference > 1.5  -> use lepton[1] for m_l2jj
      2) else, |ΔR(l1,jj) - ΔR(l2,jj)| >= 1.0  -> use lepton closer to jj
      3) else, |d0_1 - d0_2| >= 0.01 mm -> use lepton with larger |d0_plain|
      4) else, skip event.
    """

    # need at least 2 jets, 2 leptons
    if len(jets_sel) < 2 or len(leptons) < 2:
        return

    # build jj
    j1 = jets_sel[0]
    j2 = jets_sel[1]
    vj1 = ROOT.TLorentzVector()
    vj2 = ROOT.TLorentzVector()
    vj1.SetPtEtaPhiM(j1["pt"], j1["eta"], j1["phi"], j1["m"])
    vj2.SetPtEtaPhiM(j2["pt"], j2["eta"], j2["phi"], j2["m"])
    vjj = vj1 + vj2
    dR_jj = vj1.DeltaR(vj2)

    # categorical validation of jj
    if dR_jj > 1.0:
        return

    # leptons must be sorted by pT: leptons[0] = l1, leptons[1] = l2
    l1 = leptons[0]
    l2 = leptons[1]

    # same sign signature
    qprod = l1["q"] * l2["q"]
    if qprod < 0:
        return

    pt1 = l1["pt"]
    pt2 = l2["pt"]
    if pt1 <= 0 or pt2 <= 0:
        return

    # decision logic
    use_event = False
    idx_lep_for_m = -1

    # 1) pT-based criterion
    logdiff = math.log(pt1) - math.log(pt2)
    if logdiff > 1.5:
        # leptons discriminated by pT, keep ordering, use [1] for m_N
        idx_lep_for_m = 1
        use_event = True
    else:
        # 2) ΔR-based criterion
        vl1 = ROOT.TLorentzVector()
        vl2 = ROOT.TLorentzVector()
        vl1.SetPtEtaPhiM(l1["pt"], l1["eta"], l1["phi"], 0.0)
        vl2.SetPtEtaPhiM(l2["pt"], l2["eta"], l2["phi"], 0.0)
        dR1 = vl1.DeltaR(vjj)
        dR2 = vl2.DeltaR(vjj)

        if abs(dR1 - dR2) >= 1.0:
            # assumption "the closer lepton is l2" from C++ logic:
            # but we explicitly choose whichever is closer
            idx_lep_for_m = 0 if dR1 < dR2 else 1
            use_event = True
        else:
            # 3) Δd0-based criterion (plain d0, no PV correction)
            d0_1 = float(l1.get("d0_plain", 0.0))
            d0_2 = float(l2.get("d0_plain", 0.0))
            if abs(d0_1 - d0_2) >= 0.03:    #abs(d0_1 - d0_2) >= 0.03:
                # choose lepton with larger |d0_plain|
                idx_lep_for_m = 0 if abs(d0_1) >= abs(d0_2) else 1
                use_event = True
            else:
                # 4) neither pT nor ΔR nor Δd0 criteria satisfied -> skip
                return

    # histogram map for this sample
    if sname not in histos:
        return
    varMap = histos[sname]

    # ---- fill m_lljj (reco-invariant) ----
    vl1 = ROOT.TLorentzVector()
    vl2 = ROOT.TLorentzVector()
    vl1.SetPtEtaPhiM(leptons[0]["pt"], leptons[0]["eta"], leptons[0]["phi"], 0.0)
    vl2.SetPtEtaPhiM(leptons[1]["pt"], leptons[1]["eta"], leptons[1]["phi"], 0.0)
    v_lljj = vl1 + vl2 + vj1 + vj2
    m_lljj = v_lljj.M()

    if "m_lljj" in varMap:
        if buff_weight:
            varMap["m_lljj"].Fill(m_lljj, w_event*1e6)
        else:
            varMap["m_lljj"].Fill(m_lljj, w_event)

    if not use_event or idx_lep_for_m < 0:
        return

    # ---- fill m_l2jj (reco-logic dependent) ----
    lep_sel = leptons[idx_lep_for_m]
    v_l2 = ROOT.TLorentzVector()
    v_l2.SetPtEtaPhiM(lep_sel["pt"], lep_sel["eta"], lep_sel["phi"], 0.0)
    v_l2jj = vjj + v_l2
    m_l2jj = v_l2jj.M()

    if "m_l2jj" in varMap:
        if buff_weight:
            varMap["m_l2jj"].Fill(m_l2jj, w_event*1e6)
        else:
            varMap["m_l2jj"].Fill(m_l2jj, w_event)


# ------------------------------------------------------------
# Process one sample: full reco (leptons + jets + displacement)
# ------------------------------------------------------------
def process_sample(sample, histos, metadata=None, use_metadata_totalW=True):
    name = sample["name"]
    stype = sample["type"]
    files = sample["files"]

    print(f"Processing sample: {name} ({stype})")

    # Branch list (union of leptons, jets, PV, tracks, truth)
    branches_common = [
        "EventInfoAuxDyn.mcEventWeights",
        "EventInfoAuxDyn.DFCommonJets_eventClean_LooseBad",
        "ElectronsAuxDyn.DFCommonElectronsLHLoose",
        "ElectronsAuxDyn.pt",
        "ElectronsAuxDyn.eta",
        "ElectronsAuxDyn.phi",
        "ElectronsAuxDyn.charge",
        "ElectronsAuxDyn.ptvarcone20",
        "ElectronsAuxDyn.topoetcone20",
        "ElectronsAuxDyn.trackParticleLinks",

        "MuonsAuxDyn.pt",
        "MuonsAuxDyn.eta",
        "MuonsAuxDyn.phi",
        "MuonsAuxDyn.charge",
        "MuonsAuxDyn.ptvarcone20",
        "MuonsAuxDyn.topoetcone20",
        "MuonsAuxDyn.inDetTrackParticleLink.m_persIndex",

        "AntiKt4EMPFlowJetsAuxDyn.JetConstitScaleMomentum_pt",
        "AntiKt4EMPFlowJetsAuxDyn.JetConstitScaleMomentum_eta",
        "AntiKt4EMPFlowJetsAuxDyn.JetConstitScaleMomentum_phi",
        "AntiKt4EMPFlowJetsAuxDyn.JetConstitScaleMomentum_m",
        "AntiKt4EMPFlowJetsAuxDyn.Jvt",
        "AntiKt4EMPFlowJetsAuxDyn.Width",
        "AntiKt4EMPFlowJetsAuxDyn.NumTrkPt500",

        "PrimaryVerticesAuxDyn.vertexType",
        "PrimaryVerticesAuxDyn.x",
        "PrimaryVerticesAuxDyn.y",
        "PrimaryVerticesAuxDyn.z",

        "InDetTrackParticlesAuxDyn.d0",
        "InDetTrackParticlesAuxDyn.z0",
        "InDetTrackParticlesAuxDyn.phi",
        "InDetTrackParticlesAuxDyn.TTVA_AMVFVertices",
        "InDetTrackParticlesAuxDyn.TTVA_AMVFWeights",
    ]

    branches_sig = [
        "TruthBSMWithDecayParticlesAuxDyn.pdgId",
        "TruthBSMWithDecayParticlesAuxDyn.decayVtxLink.m_persIndex",
        "TruthBSMWithDecayParticlesAuxDyn.prodVtxLink.m_persIndex",
        "TruthBSMWithDecayVerticesAuxDyn.x",
        "TruthBSMWithDecayVerticesAuxDyn.y",
        "TruthBSMWithDecayVerticesAuxDyn.z",
    ]

    branch_list = branches_common + (branches_sig if stype == "SIG" else [])

    tree = uproot.concatenate(
        {fname: "CollectionTree" for fname in files},
        filter_name=branch_list,
        library="ak",
    )

    # event weights
    w_all = tree["EventInfoAuxDyn.mcEventWeights"]
    event_clean_all = tree["EventInfoAuxDyn.DFCommonJets_eventClean_LooseBad"]
    e_clean_all = tree["ElectronsAuxDyn.DFCommonElectronsLHLoose"]
    w_nom = w_all[:, 0]
    total_sumw_loop = float(ak.sum(w_nom))

    # reco leptons
    e_pt   = tree["ElectronsAuxDyn.pt"]
    e_eta  = tree["ElectronsAuxDyn.eta"]
    e_phi  = tree["ElectronsAuxDyn.phi"]
    e_q    = tree["ElectronsAuxDyn.charge"]
    e_ptvar = tree["ElectronsAuxDyn.ptvarcone20"]
    e_topo  = tree["ElectronsAuxDyn.topoetcone20"]
    e_trk   = tree["ElectronsAuxDyn.trackParticleLinks"]

    m_pt   = tree["MuonsAuxDyn.pt"]
    m_eta  = tree["MuonsAuxDyn.eta"]
    m_phi  = tree["MuonsAuxDyn.phi"]
    m_q    = tree["MuonsAuxDyn.charge"]
    m_ptvar = tree["MuonsAuxDyn.ptvarcone20"]
    m_topo  = tree["MuonsAuxDyn.topoetcone20"]
    m_idx   = tree["MuonsAuxDyn.inDetTrackParticleLink.m_persIndex"]

    # jets
    j_pt   = tree["AntiKt4EMPFlowJetsAuxDyn.JetConstitScaleMomentum_pt"]
    j_eta  = tree["AntiKt4EMPFlowJetsAuxDyn.JetConstitScaleMomentum_eta"]
    j_phi  = tree["AntiKt4EMPFlowJetsAuxDyn.JetConstitScaleMomentum_phi"]
    j_m    = tree["AntiKt4EMPFlowJetsAuxDyn.JetConstitScaleMomentum_m"]
    j_jvt  = tree["AntiKt4EMPFlowJetsAuxDyn.Jvt"]
    j_width = tree["AntiKt4EMPFlowJetsAuxDyn.Width"]
    j_ntrk  = tree["AntiKt4EMPFlowJetsAuxDyn.NumTrkPt500"]

    # PV and tracks
    pv_type = tree["PrimaryVerticesAuxDyn.vertexType"]
    pv_x = tree["PrimaryVerticesAuxDyn.x"]
    pv_y = tree["PrimaryVerticesAuxDyn.y"]
    pv_z = tree["PrimaryVerticesAuxDyn.z"]

    trk_d0  = tree["InDetTrackParticlesAuxDyn.d0"]
    trk_z0  = tree["InDetTrackParticlesAuxDyn.z0"]
    trk_phi = tree["InDetTrackParticlesAuxDyn.phi"]
    trk_vtx = tree["InDetTrackParticlesAuxDyn.TTVA_AMVFVertices"]
    trk_wgt = tree["InDetTrackParticlesAuxDyn.TTVA_AMVFWeights"]

    # truth for SIG
    if stype == "SIG":
        pdgId       = tree["TruthBSMWithDecayParticlesAuxDyn.pdgId"]
        decay_idx   = tree["TruthBSMWithDecayParticlesAuxDyn.decayVtxLink.m_persIndex"]
        prod_idx    = tree["TruthBSMWithDecayParticlesAuxDyn.prodVtxLink.m_persIndex"]
        vx_truth    = tree["TruthBSMWithDecayVerticesAuxDyn.x"]
        vy_truth    = tree["TruthBSMWithDecayVerticesAuxDyn.y"]
        vz_truth    = tree["TruthBSMWithDecayVerticesAuxDyn.z"]

    n_events = len(w_nom)
    print(f"  Entries: {n_events}")

    h = histos[name]  # local alias

    # Event loop
    PVrejected = 0
    CleanEventrejected = 0
    tlv = ROOT.TLorentzVector()

    for ievt in range(n_events):
        # if ievt <10:
        #     print("evt", ievt, "vertexType:", pv_type[ievt])
        w = float(w_nom[ievt])
        if (event_clean_all[ievt] == 0):
            CleanEventrejected += 1 
            continue
        # --- Primary vertex selection (at least one good PV with vertexType==1) ---
        types_evt = pv_type[ievt]
        hasGoodPV = bool(ak.any(types_evt == 1))
        if not hasGoodPV:
            PVrejected += 1
            continue

        # # --- Truth displacement (SIG only) ---
        # if stype == "SIG":
        #     pdg_evt  = pdgId[ievt]
        #     dec_evt  = decay_idx[ievt]
        #     prod_evt = prod_idx[ievt]
        #     vx_evt   = vx_truth[ievt]
        #     vy_evt   = vy_truth[ievt]
        #     vz_evt   = vz_truth[ievt]

        #     mask = ak.where(abs(pdg_evt) == 9900012)[0]
        #     for p_idx in mask:
        #         d_idx = int(dec_evt[p_idx])
        #         p_idx_v = int(prod_evt[p_idx])
        #         n_vtx = len(vx_evt)
        #         if d_idx < 0 or d_idx >= n_vtx:
        #             continue
        #         if p_idx_v < 0 or p_idx_v >= n_vtx:
        #             continue
        #         xd = float(vx_evt[d_idx])
        #         yd = float(vy_evt[d_idx])
        #         zd = float(vz_evt[d_idx])
        #         xp = float(vx_evt[p_idx_v])
        #         yp = float(vy_evt[p_idx_v])
        #         zp = float(vz_evt[p_idx_v])
        #         dx = xd - xp
        #         dy = yd - yp
        #         dz = zd - zp
        #         Lxyz = math.sqrt(dx*dx + dy*dy + dz*dz)
        #         if Lxyz > 1e-10:
        #             h["disp_truth_Lxyz"].Fill(Lxyz, w)

        # --- Build leptons with displacement ---
        leptons = []  # each: dict(pt,eta,phi,q,type,disp3D)

        # PV coordinates: pick PV with type==1, highest sum of weights, etc. Here just the first type==1.
        pvx_evt = pv_x[ievt]
        pvy_evt = pv_y[ievt]
        pvz_evt = pv_z[ievt]

        good_pv_indices = ak.where(pv_type[ievt] == 1)[0]
        if len(good_pv_indices) > 0:
            pv_idx0 = int(good_pv_indices[0])
        else:
            pv_idx0 = 0 if len(pvx_evt) > 0 else -1

        xPV = float(pvx_evt[pv_idx0]) if pv_idx0 >= 0 and pv_idx0 < len(pvx_evt) else 0.0
        yPV = float(pvy_evt[pv_idx0]) if pv_idx0 >= 0 and pv_idx0 < len(pvy_evt) else 0.0
        zPV = float(pvz_evt[pv_idx0]) if pv_idx0 >= 0 and pv_idx0 < len(pvz_evt) else 0.0

        trk_d0_evt  = trk_d0[ievt]
        trk_z0_evt  = trk_z0[ievt]
        trk_phi_evt = trk_phi[ievt]
        trk_vtx_evt = trk_vtx[ievt]
        trk_wgt_evt = trk_wgt[ievt]

        # electrons
        e_quality_evt = e_clean_all[ievt]
        e_pt_evt   = e_pt[ievt]
        e_eta_evt  = e_eta[ievt]
        e_phi_evt  = e_phi[ievt]
        e_q_evt    = e_q[ievt]
        e_ptvar_evt = e_ptvar[ievt]
        e_topo_evt  = e_topo[ievt]
        e_trk_evt   = e_trk[ievt]

        for i in range(len(e_pt_evt)):
            pt_mev = float(e_pt_evt[i])
            if e_quality_evt[i] != 1:
                continue
            if pt_mev <= 0:
                continue
            pt = pt_mev * 0.001  # MeV -> GeV
            eta = float(e_eta_evt[i])
            phi = float(e_phi_evt[i])
            q   = int(e_q_evt[i])

            if pt <= lep_pt_min or abs(eta) > lep_eta_max:
                continue

            rel_ptvar = (float(e_ptvar_evt[i]) / pt_mev) if pt_mev > 0 else 999.0
            rel_topo  = (float(e_topo_evt[i]) / pt_mev) if pt_mev > 0 else 999.0
            if rel_ptvar > lep_rel_ptvarcone20_max:
                continue
            if rel_topo > lep_rel_topoetcone20_max:
                continue

            link_rec = e_trk_evt[i]
            if len(link_rec) == 0:
                continue
            trk_idx = int(link_rec[0]["m_persIndex"])
            if trk_idx < 0 or trk_idx >= len(trk_d0_evt):
                continue

            d0_beam = float(trk_d0_evt[trk_idx])
            z0_beam = float(trk_z0_evt[trk_idx])

            # # Find best associated vertex via TTVA weights
            # if len(trk_vtx_evt[trk_idx]) > 0:
            #     vert_idx_list = trk_vtx_evt[trk_idx]
            #     weight_list   = trk_wgt_evt[trk_idx]
            #     best = ak.argmax(weight_list)
            #     best_vertex_link  = vert_idx_list[best]
            #     best_vertex_index = int(best_vertex_link["m_persIndex"])
            #     if 0 <= best_vertex_index < len(pvx_evt):
            #         xPV = float(pvx_evt[best_vertex_index])
            #         yPV = float(pvy_evt[best_vertex_index])
            #         zPV = float(pvz_evt[best_vertex_index])

            # # shift from beam to PV
            # shift = xPV * math.sin(phi) - yPV * math.cos(phi)
            # d0_PV = d0_beam - shift
            # z0_PV = z0_beam - zPV
            # x_rel = -d0_PV * math.sin(phi)
            # y_rel =  d0_PV * math.cos(phi)
            # z_rel =  z0_PV
            # disp3D = math.sqrt(x_rel*x_rel + y_rel*y_rel + z_rel*z_rel)

            leptons.append(dict(
                pt=pt,
                eta=eta,
                phi=phi,
                q=q,
                type="e",
                d0_plain=d0_beam
            ))

            # if stype == "SIG":
            #     # fill reco displacement histogram for signal
            #     h["disp_reco_d0"].Fill(d0_beam, w)

        # muons
        m_pt_evt   = m_pt[ievt]
        m_eta_evt  = m_eta[ievt]
        m_phi_evt  = m_phi[ievt]
        m_q_evt    = m_q[ievt]
        m_ptvar_evt = m_ptvar[ievt]
        m_topo_evt  = m_topo[ievt]
        m_idx_evt   = m_idx[ievt]

        for i in range(len(m_pt_evt)):
            pt_mev = float(m_pt_evt[i])
            if pt_mev <= 0:
                continue
            pt = pt_mev * 0.001
            eta = float(m_eta_evt[i])
            phi = float(m_phi_evt[i])
            q   = int(m_q_evt[i])

            if pt <= lep_pt_min or abs(eta) > lep_eta_max:
                continue

            rel_ptvar = (float(m_ptvar_evt[i]) / pt_mev) if pt_mev > 0 else 999.0
            rel_topo  = (float(m_topo_evt[i]) / pt_mev) if pt_mev > 0 else 999.0
            if rel_ptvar > lep_rel_ptvarcone20_max:
                continue
            if rel_topo > lep_rel_topoetcone20_max:
                continue

            trk_idx = int(m_idx_evt[i])
            if trk_idx < 0 or trk_idx >= len(trk_d0_evt):
                continue

            d0_beam = float(trk_d0_evt[trk_idx])
            z0_beam = float(trk_z0_evt[trk_idx])

            # if len(trk_vtx_evt[trk_idx]) > 0:
            #     vert_idx_list = trk_vtx_evt[trk_idx]
            #     weight_list   = trk_wgt_evt[trk_idx]
            #     best = ak.argmax(weight_list)
            #     best_vertex_link  = vert_idx_list[best]
            #     best_vertex_index = int(best_vertex_link["m_persIndex"])
            #     if 0 <= best_vertex_index < len(pvx_evt):
            #         xPV = float(pvx_evt[best_vertex_index])
            #         yPV = float(pvy_evt[best_vertex_index])
            #         zPV = float(pvz_evt[best_vertex_index])

            # shift = xPV * math.sin(phi) - yPV * math.cos(phi)
            # d0_PV = d0_beam - shift
            # z0_PV = z0_beam - zPV
            # x_rel = -d0_PV * math.sin(phi)
            # y_rel =  d0_PV * math.cos(phi)
            # z_rel =  z0_PV
            # disp3D = math.sqrt(x_rel*x_rel + y_rel*y_rel + z_rel*z_rel)

            leptons.append(dict(
                pt=pt,
                eta=eta,
                phi=phi,
                q=q,
                type="m",
                d0_plain=d0_beam,
            ))

            # if stype == "SIG":
            #     h["disp_reco_d0"].Fill(d0_beam, w)

        # sort leptons by pt (for lep1/lep2 and charges)
        leptons.sort(key=lambda L: L["pt"], reverse=True)

        # --- Build and select jets ---
        jets = []
        j_pt_evt   = j_pt[ievt]
        j_eta_evt  = j_eta[ievt]
        j_phi_evt  = j_phi[ievt]
        j_m_evt    = j_m[ievt]
        j_jvt_evt  = j_jvt[ievt]
        j_width_evt = j_width[ievt]
        j_ntrk_evt  = j_ntrk[ievt]

        for i in range(len(j_pt_evt)):
            pt = float(j_pt_evt[i]) * 0.001
            eta = float(j_eta_evt[i])
            phi = float(j_phi_evt[i])
            m_j = float(j_m_evt[i]) * 0.001
            jvt_val = float(j_jvt_evt[i])
            width_raw = j_width_evt[i]
            ntrk_raw  = j_ntrk_evt[i]

            def _as_scalar(x):
            # unwrap width_raw and ntrk_raw if they are 0‑ or 1‑D arrays/lists
                if isinstance(x, (ak.highlevel.Array, np.ndarray, list)):
                    if len(x) == 0:
                        return 0.0
                    return float(x[0])
                return float(x)

            width_val = _as_scalar(width_raw)
            ntrk_val  = _as_scalar(ntrk_raw)

            if pt <= jet_pt_min:
                continue
            if abs(eta) > jet_eta_max:
                continue
            if jvt_val < jvt_min:
                continue
            if width_val > jWidth_cut:
                continue
            if ntrk_val > jNumTrkPt500_cut:
                continue

            jets.append(dict(pt=pt, eta=eta, phi=phi, m=m_j))

        jets.sort(key=lambda J: J["pt"], reverse=True)

        # --- Fill basic histos: pts, dilep_q, DR(l,jj) ---
        if len(jets) >= 2:
            j1 = jets[0]
            j2 = jets[1]
            vj1 = ROOT.TLorentzVector()
            vj2 = ROOT.TLorentzVector()
            vj1.SetPtEtaPhiM(j1["pt"], j1["eta"], j1["phi"], j1["m"])
            vj2.SetPtEtaPhiM(j2["pt"], j2["eta"], j2["phi"], j2["m"])
            dR_jj = vj1.DeltaR(vj2)
            if dR_jj <= 2.0:
                h["jet_pt"].Fill(j1["pt"], w)
                h["jet_pt"].Fill(j2["pt"], w)
                vjj = vj1+vj2
                m_jj = vjj.M()
                h["m_jj"].Fill(m_jj, w)

        if len(leptons) > 0:
            h["lep1_pt"].Fill(leptons[0]["pt"], w)
        if len(leptons) > 1:
            h["lep2_pt"].Fill(leptons[1]["pt"], w)
            # dilepton charge
            hq = h["dilep_q"]
            hq.GetXaxis().SetBinLabel(1, "OS")
            hq.GetXaxis().SetBinLabel(2, "")
            hq.GetXaxis().SetBinLabel(3, "SS")
            qprod = leptons[0]["q"] * leptons[1]["q"]
            if qprod < 0:
                hq.Fill(-0.5, w)
            elif qprod > 0:
                hq.Fill(+0.5, w)

            d0_1 = abs(float(leptons[0].get("d0_plain", 0.0)))
            d0_2 = abs(float(leptons[1].get("d0_plain", 0.0)))
            dd0 = max(d0_1, d0_2)  #abs
            if (stype == "SIG"):
                h["disp"].Fill(dd0, w*1e7)
            else:
                h["disp"].Fill(dd0, w)

        # DR(l, jj) with selected jets
        if len(jets) >= 2 and len(leptons) > 0:
            vj1 = ROOT.TLorentzVector()
            vj2 = ROOT.TLorentzVector()
            vj1.SetPtEtaPhiM(jets[0]["pt"], jets[0]["eta"], jets[0]["phi"], jets[0]["m"])
            vj2.SetPtEtaPhiM(jets[1]["pt"], jets[1]["eta"], jets[1]["phi"], jets[1]["m"])
            dR_jj = vj1.DeltaR(vj2)
            if dR_jj <= 2.0:
                vjj = vj1 + vj2
                # leading lepton
                l1 = leptons[0]
                vl1 = ROOT.TLorentzVector()
                vl1.SetPtEtaPhiM(l1["pt"], l1["eta"], l1["phi"], 0.0)
                dR1 = vl1.DeltaR(vjj)
                h["dr_lep1_jj"].Fill(dR1, w)
                if len(leptons) > 1:
                    l2 = leptons[1]
                    vl2 = ROOT.TLorentzVector()
                    vl2.SetPtEtaPhiM(l2["pt"], l2["eta"], l2["phi"], 0.0)
                    dR2 = vl2.DeltaR(vjj)
                    h["dr_lep2_jj"].Fill(dR2, w)

        # --- Reconstruction of m_lljj and m_l2jj ---
        jets_sel = jets
        buff_weight = (stype == "SIG")
        if stype != "SIG":
            fill_reco_variables(
                sname=name,
                leptons=leptons,
                jets_sel=jets_sel,
                w_event=w,
                buff_weight=False,
                histos=histos,
            )
        else:
            fill_reco_variables(
                sname=name,
                leptons=leptons,
                jets_sel=jets_sel,
                w_event=w,
                buff_weight=True,
                histos=histos,
            )

    print(f"  total_sumw_loop = {total_sumw_loop}")
    print(f"  #PV-rejected events: {PVrejected}")
    print(f"  #False-EventClean events: {CleanEventrejected}")

    # Metadata totalW
    total_sumw_used = total_sumw_loop
    if use_metadata_totalW and metadata is not None and stype != "DATA":
        total_sumw_meta = get_total_sumw_from_metadata(sample, metadata)
        if total_sumw_meta > 0.0:
            total_sumw_used = total_sumw_meta
            print(f"  Using metadata total_sumw (sum AllExecutedEvents.W) = {total_sumw_meta}")
        else:
            print(f"  WARNING: metadata total_sumw <= 0 for {name}, falling back to loop total_sumw.")

    return total_sumw_loop, total_sumw_used


# ------------------------------------------------------------
# Book histos for all samples, process them, and scale
# ------------------------------------------------------------
histos = {}  # histos[sname][var]
super_factors = {}

for s in samples:
    sname = s["name"]
    stype = s["type"]
    histos[sname] = book_sample_histograms(sname, stype)

for s in samples:
    sname = s["name"]
    total_loop, total_used = process_sample(
        s,
        histos,
        metadata=metadata if use_metadata_totalW else None,
        use_metadata_totalW=use_metadata_totalW,
    )
    sigma = s["sigma_pb"]
    lumi  = s["lumi_pb"]
    eff   = s["filter_eff"]
    kfac  = s["k_factor"]

    if s["type"] == "DATA":
        sf = 1.0
    elif total_used != 0.0:
        sf = sigma * lumi * eff * kfac / total_used
    else:
        sf = 0.0
        print(f"  WARNING: total_sumw_used=0 for {sname}, super_factor=0")

    print(f"  super_factor({sname}) = {sf}")
    super_factors[sname] = sf

    for var, h in histos[sname].items():
        if s["type"] != "DATA":
            h.Scale(sf)

# ------------------------------------------------------------
# Build pooled S,B,Z histograms for all variables
# ------------------------------------------------------------
variables = ["lep1_pt", "lep2_pt", "disp", "jet_pt",
             "dilep_q", "dr_lep1_jj", "dr_lep2_jj",
             "m_l2jj", "m_lljj", "m_jj"]

pooled_S = {}
pooled_B = {}
pooled_Z = {}

for var in variables:
    sig_template = None
    bkg_template = None
    for s in samples:
        if var not in histos[s["name"]]:
            continue
        h = histos[s["name"]][var]
        if s["type"].startswith("SIG") and sig_template is None:
            sig_template = h
        if s["type"].startswith("BKG") and bkg_template is None:
            bkg_template = h
    if not sig_template or not bkg_template:
        continue

    hS = sig_template.Clone(f"h_pooled_S_{var}")
    hS.Reset("ICES")
    hB = bkg_template.Clone(f"h_pooled_B_{var}")
    hB.Reset("ICES")

    for s in samples:
        if var not in histos[s["name"]]:
            continue
        h = histos[s["name"]][var]
        if s["type"] == "SIG":
            hS.Add(h)
        elif s["type"].startswith("BKG"):
            hB.Add(h)

    pooled_S[var] = hS
    pooled_B[var] = hB

    hZ = hS.Clone(f"h_significance_{var}")
    hZ.Reset("ICES")
    hZ.SetTitle(f";{hS.GetXaxis().GetTitle()};Z [#sigma]")

    b_min = 1e-3
    nb = hS.GetNbinsX()
    for ib in range(1, nb + 1):
        s_val = hS.GetBinContent(ib)

        if (var == "m_lljj"):
            s_val /= 1e6 
        elif (var == "m_l2jj"):
            s_val /= 1e6
        elif (var == "disp"):
            s_val /= 1e7

        b_val = hB.GetBinContent(ib)
        Z = 0.0
        if b_val >= b_min and s_val > 0:
            arg = 1.0 + s_val / b_val
            if arg > 0:
                val = 2.0 * ((s_val + b_val) * math.log(arg) - s_val)
                if val > 0:
                    Z = math.sqrt(val)
        hZ.SetBinContent(ib, Z)

    pooled_Z[var] = hZ

# ------------------------------------------------------------
# Build group histos (SIG, Zmumu, Zee, ttbar) for each variable
# ------------------------------------------------------------
group_histos = {var: {} for var in variables}

for var in variables:
    for gname, members in group_definitions.items():
        template = None
        for m in members:
            if m in histos and var in histos[m]:
                template = histos[m][var]
                break
        if template is None:
            continue
        h_group = template.Clone(f"h_{gname}_{var}")
        h_group.Reset("ICES")
        for m in members:
            if m in histos and var in histos[m]:
                h_group.Add(histos[m][var])
        group_histos[var][gname] = h_group

# Smooth backgrounds a bit (like in C++)
# for var in variables:
#     for gname, h in group_histos[var].items():
#         if gname != "DATA":  
#             if h.GetNbinsX() > 10:
#                 h.Smooth(1)

# ------------------------------------------------------------
# Plotting: generic canvas (like C++ make_sig_canvas)
# ------------------------------------------------------------
def make_sig_canvas(var, title_top, outname, logy_top, draw_text_bottom,
                    group_histo_var, hZ):
    ROOT.gStyle.SetOptStat(0)
    ROOT.TGaxis.SetMaxDigits(3)
    c = ROOT.TCanvas(f"c_{var}", title_top, 1000, 800)
    c.Divide(1, 2)

    logx = (var not in ["dilep_q", "dr_lep1_jj", "dr_lep2_jj", "m_l2jj", "m_jj"])

    # top pad
    c.cd(1)
    ROOT.gPad.SetPad(0.0, 0.30, 1.0, 1.0)
    ROOT.gPad.SetBottomMargin(0.02)
    if logy_top:
        ROOT.gPad.SetLogy()
    if logx:
        ROOT.gPad.SetLogx()
    ROOT.gPad.SetTicky(1)

    group_colors = {
        "Zmumu": ROOT.kRed + 2,
        "Zee":   ROOT.kBlue + 2,
        "ttbar": ROOT.kGreen + 3,
        "SIG":   ROOT.kOrange + 1,
    }

    # legend
    if var != "disp":
        x1, y1, x2, y2 = 0.76, 0.65, 0.94, 0.90
    else:
        x1, y1, x2, y2 = 0.76, 0.65, 0.94, 0.90   #0.37, 0.68, 0.52, 0.90
    leg = ROOT.TLegend(x1, y1, x2, y2)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.04)

    bkg_order = ["ttbar", "Zmumu", "Zee"]
    sig_name = "SIG"

    xmin_vis = 0.0
    xmax_vis = 0.0
    if var in ["lep1_pt", "lep2_pt", "jet_pt"]:
        xmin_vis = lep_pt_min
        xmax_vis = list(group_histo_var.values())[0].GetXaxis().GetXmax()
    elif var == "disp":
        xmin_vis = disp_bins[0]
        xmax_vis = disp_bins[-1]
    elif var == "dilep_q":
        xmin_vis = -1.0
        xmax_vis = 1.0
    else:
        hx = list(group_histo_var.values())[0]
        xmin_vis = hx.GetXaxis().GetXmin()
        xmax_vis = hx.GetXaxis().GetXmax()

    # build stack
    hs_bkg = ROOT.THStack(f"hs_{var}", "")
    h_sig = None

    for gname in bkg_order:
        if gname not in group_histo_var:
            continue
        h = group_histo_var[gname]
        col = group_colors[gname]
        h.SetLineColor(ROOT.kBlack)
        h.SetLineWidth(1)
        h.SetFillColor(col)
        h.SetFillStyle(1001)
        h.SetMarkerStyle(0)
        hs_bkg.Add(h, "HIST")
        leg.AddEntry(h, gname, "f")

    if sig_name in group_histo_var:
        h_sig = group_histo_var[sig_name]
        col = group_colors[sig_name]
        h_sig.SetLineColor(col)
        h_sig.SetLineWidth(4)
        h_sig.SetLineStyle(1)
        h_sig.SetFillColor(col)
        h_sig.SetFillStyle(3004)
        h_sig.SetMarkerStyle(20)
        h_sig.SetMarkerSize(0.8)
        h_sig.SetMarkerColor(col)
        leg.AddEntry(h_sig, "model#it{A*}", "l")
        if (var == "m_lljj" or var == "disp" ):
            leg.AddEntry("", "#times 10^{6}", " ")
            

    h_axis = None
    for gname in bkg_order:
        if gname in group_histo_var:
            h_axis = group_histo_var[gname]
            break
    if h_axis is None and h_sig is not None:
        h_axis = h_sig

    hs_bkg.Draw("HIST")
    if h_axis is not None:
        hs_bkg.GetXaxis().SetTitle(h_axis.GetXaxis().GetTitle())
        hs_bkg.GetYaxis().SetTitle(h_axis.GetYaxis().GetTitle())
        hs_bkg.GetXaxis().SetRangeUser(xmin_vis, xmax_vis)
        if var == "disp":
            hs_bkg.SetMaximum(2.24e8)
            hs_bkg.SetMinimum(1.5e-4)
        elif var == "m_l2jj":
            hs_bkg.SetMaximum(14.9e3)
        elif var == "m_jj":
            hs_bkg.SetMaximum(1.5e5)
            hs_bkg.SetMinimum(1.5e-4)
        elif var == "m_lljj":
            hs_bkg.SetMaximum(14.9e3)
        elif var == "dilep_q":
            hs_bkg.SetMaximum(1.5e10)
            hs_bkg.SetMinimum(1.5e-4)
        elif var == "jet_pt":
            hs_bkg.SetMaximum(1.5e8)
            hs_bkg.SetMinimum(1.5e-4)
        elif var == "lep1_pt":
            hs_bkg.SetMaximum(1.5e9)
            hs_bkg.SetMinimum(1.5e-4)
        else:
            hs_bkg.SetMaximum(1.5e8)
            hs_bkg.SetMinimum(1.5e-4)

        hs_bkg.GetXaxis().SetLabelSize(0.0)
        hs_bkg.GetXaxis().SetTitleSize(0.0)
        hs_bkg.GetYaxis().SetTitleOffset(0.95)
        hs_bkg.GetYaxis().SetTitleSize(0.042)
        hs_bkg.GetYaxis().SetLabelSize(0.04)
        hs_bkg.GetYaxis().SetTitleFont(32)
        hs_bkg.GetXaxis().SetTitleFont(32)

    if h_sig is not None:
        h_sig.Draw("HIST SAME")

    leg.Draw()

    # ATLAS label
    lat = ROOT.TLatex()
    lat.SetNDC()
    lat.SetTextSize(0.04)
    lat.SetTextFont(72)
    
    lat.DrawLatex(0.55, 0.85, "ATLAS")
    lat.SetTextFont(42)
    lat.DrawLatex(0.63, 0.85, "Internal")
    lat.SetTextFont(32)
    lat.DrawLatex(0.52, 0.8, "pp #rightarrow W_{R} #rightarrow l (N #rightarrow l j j )")
    lat.DrawLatex(0.51, 0.73, "#sqrt{s} = 13.6 TeV, #scale[0.7]{#int} L = 500 fb^{-1}")
    # if var != "disp": else:
    #     lat.DrawLatex(0.16, 0.85, "ATLAS")
    #     lat.SetTextFont(42)
    #     lat.DrawLatex(0.24, 0.85, "Internal")
    #     lat.SetTextFont(32)
    #     lat.DrawLatex(0.13, 0.8, "pp #rightarrow W_{R} #rightarrow l (N #rightarrow l j j )")
    #     lat.DrawLatex(0.12, 0.73, "#sqrt{s} = 13.6 TeV, #scale[0.7]{#int} L = 500 fb^{-1}")   #3.55
    ROOT.gPad.RedrawAxis()

    # bottom pad
    c.cd(2)
    ROOT.gPad.SetPad(0.0, 0.0, 1.0, 0.30)
    ROOT.gPad.SetTopMargin(0.05)
    ROOT.gPad.SetBottomMargin(0.35)
    ROOT.gPad.SetGridy(True)
    ROOT.gPad.SetLogy(False)
    if logx:
        ROOT.gPad.SetLogx()

    hZ.SetStats(0)
    hZ.SetMinimum(0.0)
    hZ.SetMaximum(1.2 * hZ.GetMaximum())
    hZ.GetXaxis().SetRangeUser(xmin_vis, xmax_vis)
    hZ.SetLineColor(ROOT.kBlack)
    hZ.SetLineWidth(2)
    hZ.SetFillColor(0)
    hZ.SetMarkerStyle(20)
    hZ.SetMarkerSize(0.8)

    xaxis = hZ.GetXaxis()
    yaxis = hZ.GetYaxis()
    xaxis.SetTitle(h_axis.GetXaxis().GetTitle() if h_axis else "")
    xaxis.SetTitleSize(0.11)
    xaxis.SetTitleOffset(1.2)
    xaxis.SetLabelSize(0.10)
    xaxis.SetNdivisions(510)
    xaxis.SetMoreLogLabels(True)
    xaxis.SetNoExponent(True)
    ROOT.TGaxis.SetExponentOffset(-0.045, -0.03, "y")
    yaxis.SetTitle("asimov sign. /#sigma")
    yaxis.SetTitleFont(32)
    xaxis.SetTitleFont(32)
    yaxis.SetTitleSize(0.10)
    yaxis.SetMaxDigits(3)
    yaxis.SetTitleOffset(0.4)
    yaxis.SetLabelSize(0.08)
    yaxis.SetNdivisions(505)

    hZ.Draw("PE1")
    ROOT.gPad.Modified()
    ROOT.gPad.Update()
    ROOT.gPad.RedrawAxis()
    ROOT.gPad.Update()

    c.SaveAs(outname)
    print(f"Saved {outname}")


# ------------------------------------------------------------
# Make all canvases (like C++)
# ------------------------------------------------------------
print(f"Plotting results to {output_dir}")

make_sig_canvas(
    "lep1_pt",
    "Leading lepton p_{T}",
    os.path.join(output_dir, "leading_lepton_pt_with_significance.png"),
    True,
    False,
    group_histos["lep1_pt"],
    pooled_Z["lep1_pt"],
)
make_sig_canvas(
    "lep2_pt",
    "Subleading lepton p_{T}",
    os.path.join(output_dir, "subleading_lepton_pt_with_significance.png"),
    True,
    False,
    group_histos["lep2_pt"],
    pooled_Z["lep2_pt"],
)
make_sig_canvas(
    "jet_pt",
    "Di-jet p_{T}",
    os.path.join(output_dir, "dijet_pt_with_significance.png"),
    True,
    False,
    group_histos["jet_pt"],
    pooled_Z["jet_pt"],
)
make_sig_canvas(
    "dilep_q",
    "Dilepton charge (OS / SS)",
    os.path.join(output_dir, "dilepton_charge_with_significance.png"),
    True,
    True,
    group_histos["dilep_q"],
    pooled_Z["dilep_q"],
)
make_sig_canvas(
    "dr_lep1_jj",
    "#DeltaR(l_{1}, jj)",
    os.path.join(output_dir, "dr_lep1_jj_with_significance.png"),
    True,
    False,
    group_histos["dr_lep1_jj"],
    pooled_Z["dr_lep1_jj"],
)
make_sig_canvas(
    "dr_lep2_jj",
    "#DeltaR(l_{2}, jj)",
    os.path.join(output_dir, "dr_lep2_jj_with_significance.png"),
    True,
    False,
    group_histos["dr_lep2_jj"],
    pooled_Z["dr_lep2_jj"],
)
make_sig_canvas(
    "m_l2jj",
    "Invariant mass m(l_{2}jj) (L2 = most displaced lepton)",
    os.path.join(output_dir, "m_l2jj_with_significance.png"),
    False,
    False,
    group_histos["m_l2jj"],
    pooled_Z["m_l2jj"],
)
make_sig_canvas(
    "m_lljj",
    "Invariant mass m(l_{1}l_{2}jj)",
    os.path.join(output_dir, "m_lljj_with_significance.png"),
    False,
    False,
    group_histos["m_lljj"],
    pooled_Z["m_lljj"],
)
make_sig_canvas(
    "m_jj",
    "Invariant mass m(jj)",
    os.path.join(output_dir, "m_jj_with_significance.png"),
    True,
    False,
    group_histos["m_jj"],
    pooled_Z["m_jj"],
)

make_sig_canvas(
    "disp",
    rf"max |d_{0}(l_{1}, l_{2})|",
    os.path.join(output_dir, "dd0_lep_with_significance.png"),
    False,  
    False,
    group_histos["disp"],
    pooled_Z["disp"],
)

print(f"\nAnalysis successfully completed. Results saved to {output_dir}")
