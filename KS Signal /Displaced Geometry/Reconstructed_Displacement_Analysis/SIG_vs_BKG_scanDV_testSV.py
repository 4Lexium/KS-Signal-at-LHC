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
# Define helper functions: 
# ------------------------------------------------------------
def make_log_bins(xmin, xmax, nbins):
    """make log-spaced bin edges"""
    xmin = max(xmin, 1e-6)
    return np.logspace(np.log10(xmin), np.log10(xmax), nbins + 1)

def basename(path: str) -> str:
    """basename helper"""
    return os.path.basename(path)

def extract_meta_topkey(path: str) -> str:
    """
    Extract the top-level JSON key from the file path.
    BKG: DSID f.ex. "601190".
    Signal: generator tag f.ex. "MGPy8EG_LRSM_WR_HNL_5000_30_2L2J_MadSpin".
    """
    tag = "TeV."
    pos = path.find(tag)
    if pos == -1:
        return ""
    pos += len(tag)
    pos2 = path.find(".", pos)
    if pos2 == -1:
        return ""
    token = path[pos:pos2]
    return token  

def extract_meta_filekey(path: str, top_block: dict) -> str:
    """
    Decide which file key to use inside meta[topkey].
    """
    base = basename(path)
    if isinstance(top_block, dict) and base in top_block:
        return base
    if isinstance(top_block, dict) and len(top_block) == 1:
        return next(iter(top_block.keys()))
    return base

def get_total_sumw_from_metadata(sample, meta) -> float:
    """
    Sum AllExecutedEvents.W over all files in a sample
    Inputs:
        topkey = token after 'TeV.' in the filename (BKG: DSID, SIG: generator tag)
        filekey = basename or (if only one key) the single inner key
    """
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
# Extract Total weight from Metadata, otherwise fallback to loop-sum
# ------------------------------------------------------------
use_metadata_totalW = True
metadata_json_path = "/home/alexanum/WORKSPACE/RUN_CONTROLS/ATLAS_jobdir/ROOT_files/metadataLLP1_backgrounds.json"
metadata = {}
if use_metadata_totalW:
    try:
        with open(metadata_json_path, "r") as jf:
            metadata = json.load(jf)
    except OSError as e:
        print(f"ERROR: cannot open metadata JSON: {metadata_json_path} ({e})")
        use_metadata_totalW = False

# ------------------------------------------------------------
# Sample definition (mirrors C++ Sample vector)
# ------------------------------------------------------------
"""
MCprofile: CrossSection[pb], Int.Lum.[pb^-], FilterEff, k-factor
Signal: CS(50GeV, 5TeV, 0.1)=4.072e-4pb 
        Int. Lum: 29.0493 fb^- = 29049.3 pb^-
        Filter Eff: MG (#tau events: 6%) 
Signal: CS(70GeV, 5TeV, 0.1)=3.968e-4pb  
Signal: CS(500GeV, 5TeV, 0.1)=2.32479e-4pb  
Signal: CS(500GeV, 10TeV, 0.1)=4.672e-6pb     (e-15mb)
Signal: CS(50GeV, 10TeV, 0.1)=1.303e-5pb  
"""
samples = [
    dict(
        name="WR5000_N50_2L2J",
        type="SIG",
        sigma_pb=1.303e-5,
        lumi_pb=29049.3,
        filter_eff=0.94,
        k_factor=1.0,
        files=[
            "/storage/alexanum/DAOD_LLP_WORKFLOW/ROOT_files/DAOD_LLP1.DAOD_LLP1.mc23_13p6TeV.MGPy8EG_LRSM_WR_HNL_5000_70_2L2J_MadSpin.DAOD_LLP1.pool.root"
        ],
    ),
    dict(
        name="Zmumu",
        type="BKG_Z",
        sigma_pb=1998.8,
        lumi_pb=29049.3,
        filter_eff=1.0,
        k_factor=1.0,
        files=[
            "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.601190.PhPy8EG_AZNLO_Zmumu.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762661._000102.pool.root.1"
        ],
    ),
    # dict(
    #     name="Zmumu_CVBV",
    #     type="BKG_Z",
    #     sigma_pb=2335.9,
    #     lumi_pb=29049.3,
    #     filter_eff=8.444979e-1,
    #     k_factor=0.9332,
    #     files=[
    #         "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700791.Sh_2214_Zmumu_maxHTpTV2_CVetoBVeto.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762549._000584.pool.root.1"
    #     ],
    # ),
    # dict(
    #     name="Zmumu_BF",
    #     type="BKG_Z",
    #     sigma_pb=2336.1,
    #     lumi_pb=29049.3,
    #     filter_eff=2.502305e-2,
    #     k_factor=0.9332,
    #     files=[
    #         "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700789.Sh_2214_Zmumu_maxHTpTV2_BFilter.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762363._000036.pool.root.1"
    #     ],
    # ),
    # dict(
    #     name="Zmumu_CFBV",
    #     type="BKG_Z",
    #     sigma_pb=2336.1,
    #     lumi_pb=29049.3,
    #     filter_eff=1.304099e-1,
    #     k_factor=0.9332,
    #     files=[
    #         "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700790.Sh_2214_Zmumu_maxHTpTV2_CFilterBVeto.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762751._000027.pool.root.1"
    #     ],
    # ),
    # dict(
    #     name="Zee_CVBV",
    #     type="BKG_Z",
    #     sigma_pb=2336.1,
    #     lumi_pb=29049.3,
    #     filter_eff=8.442598e-1,
    #     k_factor=0.9332,
    #     files=[
    #         "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700788.Sh_2214_Zee_maxHTpTV2_CVetoBVeto.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762860._000141.pool.root.1"
    #     ],
    # ),
    # dict(
    #     name="Zee_BF",
    #     type="BKG_Z",
    #     sigma_pb=2336.1,
    #     lumi_pb=29049.3,
    #     filter_eff=2.558475e-2,
    #     k_factor=0.9332,
    #     files=[
    #         "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700786.Sh_2214_Zee_maxHTpTV2_BFilter.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47656889._000048.pool.root.1"
    #     ],
    # ),
    # dict(
    #     name="Zee_CFBV",
    #     type="BKG_Z",
    #     sigma_pb=2336.1,
    #     lumi_pb=29049.3,
    #     filter_eff=1.301249e-1,
    #     k_factor=0.9332,
    #     files=[
    #         "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.700787.Sh_2214_Zee_maxHTpTV2_CFilterBVeto.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47658627._000207.pool.root.1"
    #     ],
    # ),
    # dict(
    #     name="ttbar_llqq",
    #     type="BKG_T",
    #     sigma_pb=811.29,
    #     lumi_pb=29049.3,
    #     filter_eff=4.384566e-1,
    #     k_factor=1.138433852,
    #     files=[
    #         "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.601229.PhPy8EG_A14_ttbar_hdamp258p75_SingleLep.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47659319._000078.pool.root.1"
    #     ],
    # ),
    # dict(
    #     name="ttbar_llll",
    #     type="BKG_T",
    #     sigma_pb=85.482,
    #     lumi_pb=29049.3,
    #     filter_eff=1.0,
    #     k_factor=1.138433852,
    #     files=[
    #         "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.601230.PhPy8EG_A14_ttbar_hdamp258p75_dil.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762465._000196.pool.root.1"
    #     ],
    # ),
    # dict(
    #     name="ttbar_qqqq",
    #     type="BKG_T",
    #     sigma_pb=811.29,
    #     lumi_pb=29049.3,
    #     filter_eff=4.561727e-1,
    #     k_factor=1.1384339,
    #     files=[
    #         "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.601237.PhPy8EG_A14_ttbar_hdamp258p75_allhad.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47657677._000301.pool.root.1"
    #     ],
    # ),
]

# group definitions (same as C++)
group_definitions = {
    "SIG":   ["WR5000_N50_2L2J"],
    "Zmumu": ["Zmumu", "Zmumu_CVBV", "Zmumu_BF", "Zmumu_CFBV"],
    "Zee":   ["Zee_CVBV", "Zee_BF", "Zee_CFBV"],
    "ttbar": ["ttbar_llqq", "ttbar_llll", "ttbar_qqqq"],
}

output_dir = "/home/alexanum/WORKSPACE/RUN_CONTROLS/ATLAS_jobdir/EVENT_analysis/Figures/SIGvsBKG_DV_testnewDV_70_5"
os.makedirs(output_dir, exist_ok=True)

# ------------------------------------------------------------
# Histogram binning for displacement
# ------------------------------------------------------------
L_min = 0.01
L_max = 110.0
nbins = 100
edges = make_log_bins(L_min, L_max, nbins)

def make_disp_hist(name, title):
    h = ROOT.TH1F(name, title, nbins, edges)
    ROOT.TH1.SetDefaultSumw2(True)
    return h

# -----------------------------------------------------------
# Global Selection Criteria for SIG
# -----------------------------------------------------------
lep_pt_min  = 10.0  # GeV
lep_eta_max = 2.5
jet_pt_min  = 10.0   # GeV

# isolation cuts
lep_rel_ptvarcone20_max  = 0.15
lep_rel_topoetcone20_max = 0.15

# -----------------------------------------------------------
# Process one sample: build multiple reco histos + truth
# ------------------------------------------------------------
def process_sample(sample, metadata=None, use_metadata_totalW=True):
    name = sample["name"]
    stype = sample["type"]
    filenames = sample["files"]

    print(f"Processing sample: {name} ({stype})")

    branches_common = [
        "EventInfoAuxDyn.mcEventWeights",
        "ElectronsAuxDyn.pt",
        "ElectronsAuxDyn.trackParticleLinks",
        "MuonsAuxDyn.pt",
        "MuonsAuxDyn.eta",
        "ElectronsAuxDyn.eta",
        "ElectronsAuxDyn.ptvarcone20",
        "ElectronsAuxDyn.topoetcone20",
        "MuonsAuxDyn.ptvarcone20",
        "MuonsAuxDyn.topoetcone20",
        "MuonsAuxDyn.inDetTrackParticleLink.m_persIndex",
        "InDetTrackParticlesAuxDyn.d0",
        "InDetTrackParticlesAuxDyn.z0",
        "InDetTrackParticlesAuxDyn.phi",
        "InDetTrackParticlesAuxDyn.TTVA_AMVFVertices",
        "InDetTrackParticlesAuxDyn.TTVA_AMVFWeights",
        "PrimaryVerticesAuxDyn.x", 
        "PrimaryVerticesAuxDyn.y", 
        "PrimaryVerticesAuxDyn.z",
        # SV-related branches
        "InDetTrackParticlesAuxDyn.eta_wrtSV_IDAndGSF_LepTrack",
        "InDetTrackParticlesAuxDyn.phi_wrtSV_IDAndGSF_LepTrack",
        "InDetTrackParticlesAuxDyn.d0_wrtSV_IDAndGSF_LepTrack",
        "InDetTrackParticlesAuxDyn.z0_wrtSV_IDAndGSF_LepTrack",
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
        {fname: "CollectionTree" for fname in filenames},
        filter_name=branch_list,
        library="ak",
    )

    if name == "WR5000_N50_2L2J":
        print("Fields in concatenated tree:")
        for f in tree.fields:
            if "Muon" in f or "inDetTrackParticleLink" in f:
                print("  ", f)

    # nominal weights
    w_all = tree["EventInfoAuxDyn.mcEventWeights"]
    w_nom = w_all[:, 0]
    total_sumw_loop = float(ak.sum(w_nom))

    # reco branches
    el_pt   = tree["ElectronsAuxDyn.pt"]
    el_eta  = tree["ElectronsAuxDyn.eta"]
    el_ptvarcone= tree["ElectronsAuxDyn.ptvarcone20"]
    el_topoetacone = tree["ElectronsAuxDyn.topoetcone20"]
    el_trk  = tree["ElectronsAuxDyn.trackParticleLinks"]
    mu_pt   = tree["MuonsAuxDyn.pt"]
    mu_eta  = tree["MuonsAuxDyn.eta"]
    mu_ptvarcone= tree["MuonsAuxDyn.ptvarcone20"]
    mu_topoetacone = tree["MuonsAuxDyn.topoetcone20"]
    mu_idx_all = tree["MuonsAuxDyn.inDetTrackParticleLink.m_persIndex"]

    trk_TTVA_vertices = tree["InDetTrackParticlesAuxDyn.TTVA_AMVFVertices"]
    trk_TTVA_weights  = tree["InDetTrackParticlesAuxDyn.TTVA_AMVFWeights"]

    trk_d0   = tree["InDetTrackParticlesAuxDyn.d0"]
    trk_z0   = tree["InDetTrackParticlesAuxDyn.z0"]
    trk_phi  = tree["InDetTrackParticlesAuxDyn.phi"]

    # Ltrk_d0    = tree["InDetLargeD0TrackParticlesAuxDyn.d0"].array(library="ak")
    # Ltrk_z0    = tree["InDetLargeD0TrackParticlesAuxDyn.z0"].array(library="ak")
    # Ltrk_phi   = tree["InDetLargeD0TrackParticlesAuxDyn.phi"].array(library="ak")
    # Ltrk_pt    = tree['InDetLargeD0TrackParticlesAuxDyn.pt_wrtSV'].array(library="ak")

    # SV-frame branches
    trk_eta_sv = tree["InDetTrackParticlesAuxDyn.eta_wrtSV_IDAndGSF_LepTrack"]
    trk_phi_sv = tree["InDetTrackParticlesAuxDyn.phi_wrtSV_IDAndGSF_LepTrack"]
    trk_d0_sv  = tree["InDetTrackParticlesAuxDyn.d0_wrtSV_IDAndGSF_LepTrack"]
    trk_z0_sv  = tree["InDetTrackParticlesAuxDyn.z0_wrtSV_IDAndGSF_LepTrack"]

    PV_x = tree["PrimaryVerticesAuxDyn.x"]
    PV_y = tree["PrimaryVerticesAuxDyn.y"]
    PV_z = tree["PrimaryVerticesAuxDyn.z"]

    n_events = len(w_nom)
    print(f"  Entries: {n_events}")

    # histograms per sample: three reco variants
    y_title = "event weight: #sigma #int #font[12]{L} #epsilon k W^{i} /#sum W^{i}"
    h_disp_reco_plain = make_disp_hist(
        f"h_{name}_disp_reco_plain",
        f" ;L_{{reco}}^{{plain}} / mm;{y_title}",
    )
    h_disp_reco_pv = make_disp_hist(
        f"h_{name}_disp_reco_pv",
        f" ;L_{{reco}}^{{d0,PV}} / mm;{y_title}",
    )
    h_disp_reco_sv = make_disp_hist(
        f"h_{name}_disp_reco_sv",
        f" ;L_{{reco}}^{{SV}} / mm;{y_title}",
    )
    h_Lxyz_reco_L = make_disp_hist(
        f"h_{name}_disp_reco_LD0",
        f" ;L_{{reco}}^{{LD0,PV}} / mm;{y_title}",
    )
    h_disp_truth = None
    if stype == "SIG":
        h_disp_truth = make_disp_hist(
            f"h_{name}_disp_truth",
            f" ;L_{{truth}}^{{xyz}} / mm;{y_title}",
        )

    # truth branches
    if stype == "SIG":
        pdgId       = tree["TruthBSMWithDecayParticlesAuxDyn.pdgId"]
        decay_index = tree["TruthBSMWithDecayParticlesAuxDyn.decayVtxLink.m_persIndex"]
        prod_index  = tree["TruthBSMWithDecayParticlesAuxDyn.prodVtxLink.m_persIndex"]
        vx_all      = tree["TruthBSMWithDecayVerticesAuxDyn.x"]
        vy_all      = tree["TruthBSMWithDecayVerticesAuxDyn.y"]
        vz_all      = tree["TruthBSMWithDecayVerticesAuxDyn.z"]

    # event loop
    # ------------------------------------------------------------
    # xAOD invalid ElementLink sentinel
    # ------------------------------------------------------------
    INVALID_LINK = 4294967295

    # ------------------------------------------------------------
    # event loop
    # ------------------------------------------------------------
    for evt in range(n_events):

        w = float(w_nom[evt])

        # ========================================================
        # truth HNL displacement
        # ========================================================
        if stype == "SIG":

            pdg_evt = pdgId[evt]
            dec_evt = decay_index[evt]
            prod_evt = prod_index[evt]

            vx_evt = vx_all[evt]
            vy_evt = vy_all[evt]
            vz_evt = vz_all[evt]

            mask = ak.where(abs(pdg_evt) == 9900012)[0]

            for p_idx in mask:

                # ------------------------------------------------
                # RAW xAOD ElementLinks
                # ------------------------------------------------
                d_idx_raw = dec_evt[p_idx]
                p_idx_raw = prod_evt[p_idx]

                # ------------------------------------------------
                # protect against invalid xAOD links
                # ------------------------------------------------
                if d_idx_raw == INVALID_LINK:
                    continue

                if p_idx_raw == INVALID_LINK:
                    continue

                d_idx   = int(d_idx_raw)
                p_idx_v = int(p_idx_raw)

                n_vtx = len(vx_evt)

                if d_idx >= n_vtx:
                    continue

                if p_idx_v >= n_vtx:
                    continue

                # ------------------------------------------------
                # decay vertex
                # ------------------------------------------------
                xd = float(vx_evt[d_idx])
                yd = float(vy_evt[d_idx])
                zd = float(vz_evt[d_idx])

                # ------------------------------------------------
                # production vertex
                # ------------------------------------------------
                xp = float(vx_evt[p_idx_v])
                yp = float(vy_evt[p_idx_v])
                zp = float(vz_evt[p_idx_v])

                # ------------------------------------------------
                # truth Lxyz
                # ------------------------------------------------
                dx = xd - xp
                dy = yd - yp
                dz = zd - zp

                Lxyz = math.sqrt(dx*dx + dy*dy + dz*dz)

                if Lxyz > 1e-10:
                    h_disp_truth.Fill(Lxyz, w)

        # ========================================================
        # reco tracks
        # ========================================================
        leptons_std = []

        trk_d0_evt   = trk_d0[evt]
        trk_z0_evt   = trk_z0[evt]
        trk_phi_evt  = trk_phi[evt]

        # --------------------------------------------------------
        # NEW SV-frame branches
        # --------------------------------------------------------
        trk_eta_sv_evt = trk_eta_sv[evt]
        trk_phi_sv_evt = trk_phi_sv[evt]
        trk_d0_sv_evt  = trk_d0_sv[evt]
        trk_z0_sv_evt  = trk_z0_sv[evt]

        pv_x_evt = PV_x[evt]
        pv_y_evt = PV_y[evt]
        pv_z_evt = PV_z[evt]

        # ========================================================
        # ELECTRONS
        # ========================================================
        el_pt_evt  = el_pt[evt]
        el_trk_evt = el_trk[evt]

        for i in range(len(el_pt_evt)):

            link = el_trk_evt[i]

            if len(link) == 0:
                continue

            trk_idx = link[0]["m_persIndex"]
            trk_idx_int = int(trk_idx)

            # ----------------------------------------------------
            # guard track index
            # ----------------------------------------------------
            if trk_idx_int == INVALID_LINK:
                continue

            if trk_idx_int < 0 or trk_idx_int >= len(trk_d0_evt):
                continue

            vert_idx_list = trk_TTVA_vertices[evt][trk_idx_int]
            weight_list   = trk_TTVA_weights[evt][trk_idx_int]

            if len(vert_idx_list) == 0:
                continue

            d0_beam = float(trk_d0_evt[trk_idx_int])
            z0_beam = float(trk_z0_evt[trk_idx_int])
            phi     = float(trk_phi_evt[trk_idx_int])

            # Ld0 = Ltrk_d0[evt][trk_idx_int]
            # Lz0 = Ltrk_z0[evt][trk_idx_int]
            # Lphi = Ltrk_phi[evt][trk_idx_int]
            # Lx = -Ld0*np.sin(Lphi)
            # Ly = Ld0*np.cos(Lphi)
            # leptons_LD0.append((Ltrk_pt[evt][i], Lx, Ly, Lz0, "e")) 

            pt = el_pt[evt][i]*0.001 # MeV->GeV
            eta = el_eta[evt][i]
            if pt <= lep_pt_min or abs(eta) > lep_eta_max:
                continue

            rel_ptvar = (el_ptvarcone[evt][i] / el_pt[evt][i]) if el_pt[evt][i] > 0 else 999.0
            rel_topo  = (el_topoetacone[evt][i] / el_pt[evt][i]) if el_pt[evt][i] > 0 else 999.0

            if rel_ptvar > lep_rel_ptvarcone20_max:
                continue
            if rel_topo > lep_rel_topoetcone20_max:
                continue

            best = ak.argmax(weight_list)

            best_vertex_link = vert_idx_list[best]
            best_vertex_index = int(best_vertex_link["m_persIndex"])

            # ----------------------------------------------------
            # guard PV link
            # ----------------------------------------------------
            if best_vertex_index == INVALID_LINK:
                continue

            if best_vertex_index < 0 or best_vertex_index >= len(pv_x_evt):
                continue

            xPV = float(pv_x_evt[best_vertex_index])
            yPV = float(pv_y_evt[best_vertex_index])
            zPV = float(pv_z_evt[best_vertex_index])

            # ====================================================
            # plain coordinates
            # ====================================================
            x_plain = -d0_beam * math.sin(phi)
            y_plain =  d0_beam * math.cos(phi)
            z_plain =  z0_beam

            # ====================================================
            # PV-corrected coordinates
            # ====================================================
            shift = xPV * math.sin(phi) - yPV * math.cos(phi)

            d0_PV = d0_beam + shift
            z0_PV = z0_beam - zPV

            x_rel = -d0_PV * math.sin(phi)
            y_rel =  d0_PV * math.cos(phi)
            z_rel =  z0_PV

            # ====================================================
            # NEW SV-frame coordinates
            # ====================================================
            x_sv = y_sv = z_sv = None

            if (
                len(trk_d0_sv_evt)  > trk_idx_int and
                len(trk_phi_sv_evt) > trk_idx_int and
                len(trk_z0_sv_evt)  > trk_idx_int
            ):

                d0_sv_val  = float(trk_d0_sv_evt[trk_idx_int])
                z0_sv_val  = float(trk_z0_sv_evt[trk_idx_int])
                phi_sv_val = float(trk_phi_sv_evt[trk_idx_int])

                x_sv = -d0_sv_val * math.sin(phi_sv_val)
                y_sv =  d0_sv_val * math.cos(phi_sv_val)
                z_sv =  z0_sv_val

            leptons_std.append(
                (
                    float(el_pt_evt[i]),
                    (x_plain, y_plain, z_plain),
                    (x_rel,   y_rel,   z_rel),
                    (x_sv,    y_sv,    z_sv),
                    best_vertex_index,
                    "e",
                )
            )

        # ========================================================
        # MUONS
        # ========================================================
        mu_pt_evt  = mu_pt[evt]
        mu_idx_evt = mu_idx_all[evt]

        for i in range(len(mu_idx_evt)):

            trk_idx_int = int(mu_idx_evt[i])

            # ----------------------------------------------------
            # guard invalid xAOD links
            # ----------------------------------------------------
            if trk_idx_int == INVALID_LINK:
                continue

            if trk_idx_int < 0 or trk_idx_int >= len(trk_d0_evt):
                continue

            vert_idx_list = trk_TTVA_vertices[evt][trk_idx_int]
            weight_list   = trk_TTVA_weights[evt][trk_idx_int]

            if len(vert_idx_list) == 0:
                continue

            d0_beam = float(trk_d0_evt[trk_idx_int])
            z0_beam = float(trk_z0_evt[trk_idx_int])
            phi     = float(trk_phi_evt[trk_idx_int])

            pt = mu_pt[evt][i]*0.001 # MeV->GeV
            eta = mu_eta[evt][i]
            if pt <= lep_pt_min or abs(eta) > lep_eta_max:
                continue

            rel_ptvar = (mu_ptvarcone[evt][i] / mu_pt[evt][i]) if mu_pt[evt][i] > 0 else 999.0
            rel_topo  = (mu_topoetacone[evt][i] / mu_pt[evt][i]) if mu_pt[evt][i] > 0 else 999.0

            if rel_ptvar > lep_rel_ptvarcone20_max:
                continue
            if rel_topo > lep_rel_topoetcone20_max:
                continue

            best = ak.argmax(weight_list)

            best_vertex_link = vert_idx_list[best]
            best_vertex_index = int(best_vertex_link["m_persIndex"])

            # ----------------------------------------------------
            # guard invalid PV link
            # ----------------------------------------------------
            if best_vertex_index == INVALID_LINK:
                continue

            if best_vertex_index < 0 or best_vertex_index >= len(pv_x_evt):
                continue

            xPV = float(pv_x_evt[best_vertex_index])
            yPV = float(pv_y_evt[best_vertex_index])
            zPV = float(pv_z_evt[best_vertex_index])

            # ====================================================
            # plain coordinates
            # ====================================================
            x_plain = -d0_beam * math.sin(phi)
            y_plain =  d0_beam * math.cos(phi)
            z_plain =  z0_beam

            # ====================================================
            # PV-corrected coordinates
            # ====================================================
            shift = xPV * math.sin(phi) - yPV * math.cos(phi)

            d0_PV = d0_beam + shift
            z0_PV = z0_beam - zPV

            x_rel = -d0_PV * math.sin(phi)
            y_rel =  d0_PV * math.cos(phi)
            z_rel =  z0_PV

            # ====================================================
            # NEW SV-frame coordinates
            # ====================================================
            x_sv = y_sv = z_sv = None

            if (
                len(trk_d0_sv_evt)  > trk_idx_int and
                len(trk_phi_sv_evt) > trk_idx_int and
                len(trk_z0_sv_evt)  > trk_idx_int
            ):

                d0_sv_val  = float(trk_d0_sv_evt[trk_idx_int])
                z0_sv_val  = float(trk_z0_sv_evt[trk_idx_int])
                phi_sv_val = float(trk_phi_sv_evt[trk_idx_int])

                x_sv = -d0_sv_val * math.sin(phi_sv_val)
                y_sv =  d0_sv_val * math.cos(phi_sv_val)
                z_sv =  z0_sv_val

            leptons_std.append(
                (
                    float(mu_pt_evt[i]),
                    (x_plain, y_plain, z_plain),
                    (x_rel,   y_rel,   z_rel),
                    (x_sv,    y_sv,    z_sv),
                    best_vertex_index,
                    "m",
                )
            )

        # ========================================================
        # two-lepton DV reconstruction
        # ========================================================
        if len(leptons_std) >= 2:

            leptons_sorted = sorted(
                leptons_std,
                key=lambda x: x[0],
                reverse=True
            )

            l1 = leptons_sorted[0]
            l2 = leptons_sorted[1]

            pv1 = l1[4]
            pv2 = l2[4]

            # require same PV
            if pv1 != pv2:
                continue

            # ----------------------------------------------------
            # unpack coordinates
            # ----------------------------------------------------
            (x1_plain, y1_plain, z1_plain) = l1[1]
            (x2_plain, y2_plain, z2_plain) = l2[1]

            (x1_pv, y1_pv, z1_pv) = l1[2]
            (x2_pv, y2_pv, z2_pv) = l2[2]

            (x1_sv, y1_sv, z1_sv) = l1[3]
            (x2_sv, y2_sv, z2_sv) = l2[3]

            # ====================================================
            # plain separation
            # ====================================================
            dx_plain = x1_plain - x2_plain
            dy_plain = y1_plain - y2_plain
            dz_plain = z1_plain - z2_plain

            dL_plain = math.sqrt(
                dx_plain*dx_plain +
                dy_plain*dy_plain +
                dz_plain*dz_plain
            )

            # ====================================================
            # PV-corrected separation
            # ====================================================
            dx_pv = x1_pv - x2_pv
            dy_pv = y1_pv - y2_pv
            dz_pv = z1_pv - z2_pv

            dL_pv = math.sqrt(
                dx_pv*dx_pv +
                dy_pv*dy_pv +
                dz_pv*dz_pv
            )

            # ====================================================
            # NEW ATLAS-LIKE DV RECONSTRUCTION
            # ====================================================
            dL_sv = None

            if (x1_sv is not None) and (x2_sv is not None):

                # reconstructed SV position
                xSV = 0.5 * (x1_sv + x2_sv)
                ySV = 0.5 * (y1_sv + y2_sv)
                zSV = 0.5 * (z1_sv + z2_sv)

                # associated PV
                xPV = float(pv_x_evt[pv1])
                yPV = float(pv_y_evt[pv1])
                zPV = float(pv_z_evt[pv1])

                dx = xSV - xPV
                dy = ySV - yPV
                dz = zSV - zPV

                dL_sv = math.sqrt(dx*dx + dy*dy + dz*dz)

            # ====================================================
            # fill histograms
            # ====================================================
            if dL_plain > 0:
                h_disp_reco_plain.Fill(dL_plain, w)

            if dL_pv > 0:
                h_disp_reco_pv.Fill(dL_pv, w)

            if dL_sv is not None and dL_sv > 0:
                h_disp_reco_sv.Fill(dL_sv, w)

    print(f"  total_sumw_loop = {total_sumw_loop}")
    total_sumw_used = total_sumw_loop
    if use_metadata_totalW and metadata is not None:
        total_sumw_meta = get_total_sumw_from_metadata(sample, metadata)
        if total_sumw_meta > 0.0:
            total_sumw_used = total_sumw_meta
            print(f"  Using metadata total_sumw (sum AllExecutedEvents.W) = {total_sumw_meta}")
        else:
            print(f"  WARNING: metadata total_sumw <= 0 for {name}, falling back to loop total_sumw.")

    return (
        h_disp_reco_plain,
        h_disp_reco_pv,
        h_disp_reco_sv,
        h_disp_truth,
        total_sumw_loop,
        total_sumw_used,
    )

# ------------------------------------------------------------
# Process all samples, compute super_factors, scale histos
# ------------------------------------------------------------
histos_disp_reco_plain = {}
histos_disp_reco_pv    = {}
histos_disp_reco_sv    = {}
histos_disp_truth      = {}
super_factors          = {}

for s in samples:
    (
        h_reco_plain,
        h_reco_pv,
        h_reco_sv,
        h_truth,
        sumw_loop,
        sumw_used,
    ) = process_sample(
        s,
        metadata=metadata if use_metadata_totalW else None,
        use_metadata_totalW=use_metadata_totalW,
    )

    sigma = s["sigma_pb"]
    lumi  = s["lumi_pb"]
    eff   = s["filter_eff"]
    kfac  = s["k_factor"]

    if sumw_used != 0.0:
        sf = sigma * lumi * eff * kfac / sumw_used
    else:
        sf = 0.0
        print(f"  WARNING: total_sumw_used=0 for {s['name']}, super_factor=0")

    print(f"  super_factor = {sf}")
    super_factors[s["name"]] = sf

    h_reco_plain.Scale(sf)
    h_reco_pv.Scale(sf)
    h_reco_sv.Scale(sf)
    histos_disp_reco_plain[s["name"]] = h_reco_plain
    histos_disp_reco_pv[s["name"]]    = h_reco_pv
    histos_disp_reco_sv[s["name"]]    = h_reco_sv

    if h_truth is not None:
        h_truth.Scale(sf)
        histos_disp_truth[s["name"]] = h_truth

# ------------------------------------------------------------
# Build group histos for displacement (using PV-corrected histos)
# ------------------------------------------------------------
group_histos_reco = {}
for gname, members in group_definitions.items():
    template = None
    for s in members:
        if s in histos_disp_reco_pv:
            template = histos_disp_reco_pv[s]
            break
    if template is None:
        continue
    h_group = template.Clone(f"h_{gname}_disp_reco")
    h_group.Reset("ICES")
    for s in members:
        if s in histos_disp_reco_pv:
            h_group.Add(histos_disp_reco_pv[s])
    group_histos_reco[gname] = h_group

# ------------------------------------------------------------
# Build pooled S and B histos for each reco definition
# ------------------------------------------------------------
# PV-corrected (baseline)
h_pooled_S      = histos_disp_reco_pv["WR5000_N50_2L2J"].Clone("h_pooled_S_disp")
h_pooled_S.Reset("ICES")
h_pooled_B      = histos_disp_reco_pv["WR5000_N50_2L2J"].Clone("h_pooled_B_disp")
h_pooled_B.Reset("ICES")

# plain
h_pooled_S_plain = histos_disp_reco_plain["WR5000_N50_2L2J"].Clone("h_pooled_S_disp_plain")
h_pooled_S_plain.Reset("ICES")
h_pooled_B_plain = histos_disp_reco_plain["WR5000_N50_2L2J"].Clone("h_pooled_B_disp_plain")
h_pooled_B_plain.Reset("ICES")

# SV-frame
h_pooled_S_sv = histos_disp_reco_sv["WR5000_N50_2L2J"].Clone("h_pooled_S_disp_sv")
h_pooled_S_sv.Reset("ICES")
h_pooled_B_sv = histos_disp_reco_sv["WR5000_N50_2L2J"].Clone("h_pooled_B_disp_sv")
h_pooled_B_sv.Reset("ICES")

for s in samples:
    name = s["name"]
    h_plain = histos_disp_reco_plain[name]
    h_pv    = histos_disp_reco_pv[name]
    h_sv    = histos_disp_reco_sv[name]

    if s["type"] == "SIG":
        h_pooled_S.Add(h_pv)
        h_pooled_S_plain.Add(h_plain)
        h_pooled_S_sv.Add(h_sv)
    elif s["type"].startswith("BKG"):
        h_pooled_B.Add(h_pv)
        h_pooled_B_plain.Add(h_plain)
        h_pooled_B_sv.Add(h_sv)

# ------------------------------------------------------------
# significance histogram (use PV-corrected as baseline)
# ------------------------------------------------------------
hZ = h_pooled_S.Clone("h_significance_disp")
hZ.Reset("ICES")
hZ.SetTitle(";L_{reco}^{d0} / mm;Z [#sigma]")

b_min = 1e-3
nbins = h_pooled_S.GetNbinsX()
for ib in range(1, nbins + 1):
    s_val = h_pooled_S.GetBinContent(ib)
    b_val = h_pooled_B.GetBinContent(ib)
    Z = 0.0
    if b_val >= b_min and s_val > 0:
        arg = 1.0 + s_val / b_val
        if arg > 0:
            val = 2.0 * ((s_val + b_val) * math.log(arg) - s_val)
            if val > 0:
                Z = math.sqrt(val)
    hZ.SetBinContent(ib, Z)

# ------------------------------------------------------------
# Plotting: displacement with significance (PV-corrected histos)
# ------------------------------------------------------------
def plot_disp_with_significance():
    ROOT.gStyle.SetOptStat(0)
    c = ROOT.TCanvas("c_disp", "SV / HNL displacement", 1000, 800)
    c.Divide(1, 2)

    # top pad
    c.cd(1)
    ROOT.gPad.SetPad(0.0, 0.30, 1.0, 1.0)
    ROOT.gPad.SetBottomMargin(0.02)
    ROOT.gPad.SetLogy(True)
    ROOT.gPad.SetLogx(True)
    ROOT.gPad.SetTicky(1)

    bkg_order = ["ttbar", "Zmumu", "Zee"]
    sig_name  = "SIG"
    group_colors = {
        "Zmumu": ROOT.kRed + 2,
        "Zee":   ROOT.kBlue + 2,
        "ttbar": ROOT.kGreen + 3,
        "SIG":   ROOT.kOrange + 1,
    }
    xmin_vis = L_min
    xmax_vis = L_max

    x1, y1, x2, y2 = 0.37, 0.68, 0.52, 0.90
    leg = ROOT.TLegend(x1, y1, x2, y2)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.04)

    hs_bkg = ROOT.THStack("hs_disp", "")
    max_bkg = 0.0
    h_sig = None

    # backgrounds
    for gname in bkg_order:
        if gname not in group_histos_reco:
            continue
        h = group_histos_reco[gname]
        col = group_colors[gname]
        h.SetLineColor(ROOT.kBlack)
        h.SetLineWidth(1)
        h.SetFillColor(col)
        h.SetFillStyle(1001)
        h.SetMarkerStyle(0)
        max_bkg = max(max_bkg, h.GetMaximum())
        hs_bkg.Add(h, "HIST")
        leg.AddEntry(h, gname, "f")

    # signal on top
    if "SIG" in group_histos_reco:
        h_sig = group_histos_reco["SIG"]
        col = group_colors["SIG"]
        h_sig.SetLineColor(col)
        h_sig.SetLineWidth(4)
        h_sig.SetLineStyle(1)
        h_sig.SetFillColor(col)
        h_sig.SetFillStyle(3004)
        h_sig.SetMarkerStyle(20)
        h_sig.SetMarkerSize(0.8)
        h_sig.SetMarkerColor(col)
        leg.AddEntry(h_sig, "model#it{A}", "l")

    # choose axis histogram
    h_axis = None
    for gname in bkg_order:
        if gname in group_histos_reco:
            h_axis = group_histos_reco[gname]
            break
    if h_axis is None and h_sig is not None:
        h_axis = h_sig

    hs_bkg.Draw("HIST")
    if h_axis is not None:
        hs_bkg.GetXaxis().SetTitle(h_axis.GetXaxis().GetTitle())
        hs_bkg.GetYaxis().SetTitle(h_axis.GetYaxis().GetTitle())
        hs_bkg.GetXaxis().SetRangeUser(xmin_vis, xmax_vis)
        hs_bkg.SetMaximum(1.5e9)
        hs_bkg.SetMinimum(1.5e-10)
        hs_bkg.GetXaxis().SetLabelSize(0.0)
        hs_bkg.GetXaxis().SetTitleSize(0.0)
        hs_bkg.GetYaxis().SetTitleOffset(0.8)
        hs_bkg.GetYaxis().SetTitleSize(0.045)
        hs_bkg.GetYaxis().SetLabelSize(0.045)
        hs_bkg.GetYaxis().SetTitleFont(32)
        hs_bkg.GetXaxis().SetTitleFont(32)
        hs_bkg.GetXaxis().SetNdivisions(510)
        hs_bkg.GetXaxis().SetMoreLogLabels(True)

    if h_sig is not None:
        h_sig.Draw("HIST SAME")

    leg.Draw()

    lat = ROOT.TLatex()
    lat.SetNDC()
    lat.SetTextSize(0.04)
    lat.SetTextFont(72)
    lat.DrawLatex(0.16, 0.85, "ATLAS")
    lat.SetTextFont(42)
    lat.DrawLatex(0.24, 0.85, "Internal")
    lat.SetTextFont(32)
    lat.DrawLatex(0.13, 0.8, "pp #rightarrow W_{R} #rightarrow l (N #rightarrow l j j )")
    lat.DrawLatex(0.12, 0.73, "#sqrt{s} = 13.6 TeV, #int L = 29.0 fb^{-1}")
    ROOT.gPad.RedrawAxis()

    # bottom pad: significance
    c.cd(2)
    ROOT.gPad.SetPad(0.0, 0.0, 1.0, 0.30)
    ROOT.gPad.SetTopMargin(0.05)
    ROOT.gPad.SetBottomMargin(0.35)
    ROOT.gPad.SetGridy(True)
    ROOT.gPad.SetLogy(False)
    ROOT.gPad.SetLogx(True)

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
    xaxis.SetTitle("L_{reco}^{d0} / mm")
    xaxis.SetTitleSize(0.11)
    xaxis.SetTitleOffset(1.2)
    xaxis.SetLabelSize(0.10)
    ROOT.TGaxis.SetExponentOffset(-0.045, -0.03, "y")
    yaxis.SetTitle("significance /#sigma")
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

    outname = os.path.join(output_dir, "sv_displacement_with_significance_uproot.png")
    c.SaveAs(outname)
    print(f"Saved {outname}")

# ------------------------------------------------------------
# Debug: truth vs multiple reco definitions for signal only
# ------------------------------------------------------------
def plot_signal_truth_vs_reco():
    sig_name = "WR5000_N50_2L2J"
    if sig_name not in histos_disp_truth:
        print("Signal truth hist missing; skip truth-vs-reco plot.")
        return

    hTruth = histos_disp_truth[sig_name]
    global h_pooled_S, h_pooled_S_plain, h_pooled_S_sv

    hRecoS_plain = h_pooled_S_plain.Clone("hRecoS_plain")
    hRecoS_pv    = h_pooled_S.Clone("hRecoS_pv")
    hRecoS_sv    = h_pooled_S_sv.Clone("hRecoS_sv")

    if hRecoS_plain.Integral() > 0:
        hRecoS_plain.Scale(1.0 / hRecoS_plain.Integral())
    if hRecoS_pv.Integral() > 0:
        hRecoS_pv.Scale(1.0 / hRecoS_pv.Integral())
    if hRecoS_sv.Integral() > 0:
        hRecoS_sv.Scale(1.0 / hRecoS_sv.Integral())
    if hTruth.Integral() > 0:
        hTruth.Scale(1.0 / hTruth.Integral())

    L_mean = 0.135662 #0.135662   9.32
    Ntot   = 1
    hExp   = hRecoS_pv.Clone("h_disp_exponential")
    hExp.Reset("ICES")

    nb = hExp.GetNbinsX()
    for ib in range(1, nb + 1):
        x     = hExp.GetBinCenter(ib)
        width = hExp.GetBinWidth(ib)
        if L_mean > 0:
            fL = (1.0 / L_mean) * math.exp(-x / L_mean)
            val = Ntot * fL * width
            hExp.SetBinContent(ib, val)

    hRecoS_plain.SetLineColor(ROOT.kRed + 1)
    hRecoS_plain.SetLineWidth(2)
    hRecoS_plain.SetFillStyle(0)
    hRecoS_plain.SetMarkerStyle(20)
    hRecoS_plain.SetMarkerColor(ROOT.kRed + 1)
    hRecoS_plain.SetMarkerSize(0.8)

    hRecoS_pv.SetLineColor(ROOT.kGreen + 2)
    hRecoS_pv.SetLineWidth(2)
    hRecoS_pv.SetFillStyle(0)
    hRecoS_pv.SetMarkerStyle(21)
    hRecoS_pv.SetMarkerColor(ROOT.kGreen + 2)
    hRecoS_pv.SetMarkerSize(0.8)

    hRecoS_sv.SetLineColor(ROOT.kMagenta + 1)
    hRecoS_sv.SetLineWidth(2)
    hRecoS_sv.SetFillStyle(0)
    hRecoS_sv.SetMarkerStyle(22)
    hRecoS_sv.SetMarkerColor(ROOT.kMagenta + 1)
    hRecoS_sv.SetMarkerSize(0.8)

    hTruth.SetLineColor(ROOT.kBlue + 1)
    hTruth.SetLineWidth(2)
    hTruth.SetFillStyle(0)
    hTruth.SetMarkerStyle(24)
    hTruth.SetMarkerColor(ROOT.kBlue + 1)
    hTruth.SetMarkerSize(0.8)

    hExp.SetLineColor(ROOT.kBlack)
    hExp.SetLineWidth(2)
    hExp.SetFillStyle(0)
    hExp.SetMarkerStyle(1)
    hExp.SetMarkerSize(0.0)

    cL = ROOT.TCanvas("c_disp_reco_vs_truth", "Debug: reco vs truth displacement (signal)", 1000, 600)
    cL.cd()
    ROOT.gPad.SetLogy(False)
    ROOT.gPad.SetLogx(True)

    maxData = max(
        hRecoS_plain.GetMaximum(),
        hRecoS_pv.GetMaximum(),
        hRecoS_sv.GetMaximum(),
        hTruth.GetMaximum(),
    )
    maxExp  = hExp.GetMaximum()
    maxY    = max(maxData, maxExp)

    hRecoS_plain.SetMinimum(0.0)
    hRecoS_plain.SetMaximum(1.2 * maxY)

    hRecoS_plain.GetXaxis().SetTitle("L / mm")
    hRecoS_plain.GetYaxis().SetTitle("normalized event weight")
    hRecoS_plain.GetXaxis().SetTitleSize(0.045)
    hRecoS_plain.GetXaxis().SetLabelSize(0.04)
    hRecoS_plain.GetYaxis().SetTitleSize(0.045)
    hRecoS_plain.GetYaxis().SetLabelSize(0.04)

    hRecoS_plain.Draw("HIST E1")
    hRecoS_pv.Draw("HIST E1 SAME")
    hRecoS_sv.Draw("HIST E1 SAME")
    hTruth.Draw("HIST E1 SAME")
    hExp.Draw("HIST SAME")

    leg = ROOT.TLegend(0.15, 0.65, 0.65, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.035)
    leg.AddEntry(hRecoS_plain, "Rec: perigree wrt. beamline", "l")
    leg.AddEntry(hRecoS_pv,    "Rec: perigee wrt. PV", "l")
    leg.AddEntry(hRecoS_sv,    "Rec: DCA wrt. SV", "l")
    leg.AddEntry(hTruth,       "Truth: Vprod-Vdec", "l")
    leg.AddEntry(hExp,         "P(L; #it{L_{0}} = 0.135662 mm)", "l")
    leg.Draw()

    outname = os.path.join(output_dir, "debug_disp_reco_multi_vs_truth_Lxyz_signal_uproot.png")
    cL.SaveAs(outname)
    print(f"Saved {outname}")

# ------------------------------------------------------------
# Main
# ------------------------------------------------------------
if __name__ == "__main__":
    print("User Test 050:")
    plot_disp_with_significance()
    plot_signal_truth_vs_reco()
    print(f"Analysis successfully completed. Results saved to {output_dir}")
