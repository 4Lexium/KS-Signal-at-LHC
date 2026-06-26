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
            "/storage/alexanum/DAOD_LLP_WORKFLOW/ROOT_files/DAOD_LLP1.DAOD_LLP1.mc23_13p6TeV.MGPy8EG_LRSM_WR_HNL_10000_50_2L2J_MadSpin.DAOD_LLP1.pool.root"
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
    dict(
        name="ttbar_llqq",
        type="BKG_T",
        sigma_pb=811.29,
        lumi_pb=29049.3,
        filter_eff=4.384566e-1,
        k_factor=1.138433852,
        files=[
            "/storage/alexanum/DAOD_LLP_WORKFLOW/SMbackgrounds/mc23_13p6TeV.601229.PhPy8EG_A14_ttbar_hdamp258p75_SingleLep.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47659319._000078.pool.root.1"
        ],
    ),
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

output_dir = "/home/alexanum/WORKSPACE/RUN_CONTROLS/ATLAS_jobdir/EVENT_analysis/Figures/SIGvsBKG_DV_testnewDV"
os.makedirs(output_dir, exist_ok=True)

# ------------------------------------------------------------
# Histogram binning for displacement (use your Python scheme)
# ------------------------------------------------------------
L_min = 0.01
L_max = 110.0 #5.0   used to be 10
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

# ------------------------------------------------------------
# Process one sample: fill h_disp_reco (and h_disp_truth for SIG),
# accumulate total_sumw of w_nom.
# ------------------------------------------------------------
def process_sample(sample, metadata=None, use_metadata_totalW=True):
    name = sample["name"]
    stype = sample["type"]
    filenames = sample["files"]

    print(f"Processing sample: {name} ({stype})")

    # Merge all files for this sample via uproot.concatenate
    # Note: we select only branches needed for displacement
    branches_common = [
        "EventInfoAuxDyn.mcEventWeights",
        "ElectronsAuxDyn.pt",
        "ElectronsAuxDyn.eta",
        "ElectronsAuxDyn.trackParticleLinks",
        "MuonsAuxDyn.pt",
        "MuonsAuxDyn.eta",
        "ElectronsAuxDyn.ptvarcone20",
        "ElectronsAuxDyn.topoetcone20",
        "MuonsAuxDyn.ptvarcone20",
        "MuonsAuxDyn.topoetcone20",
        "MuonsAuxDyn.inDetTrackParticleLink.m_persIndex",  # <- use leaf name
        "InDetTrackParticlesAuxDyn.d0",
        "InDetTrackParticlesAuxDyn.z0",
        "InDetTrackParticlesAuxDyn.phi",
        "InDetTrackParticlesAuxDyn.TTVA_AMVFVertices",
        "InDetTrackParticlesAuxDyn.TTVA_AMVFWeights",
        "PrimaryVerticesAuxDyn.x", 
        "PrimaryVerticesAuxDyn.y", 
        "PrimaryVerticesAuxDyn.z",
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
    # take first weight as nominal, like C++
    w_nom = w_all[:, 0]
    total_sumw_loop = float(ak.sum(w_nom))

    # reco branches (electrons + muons + tracks)
    el_pt   = tree["ElectronsAuxDyn.pt"]
    el_eta  = tree["ElectronsAuxDyn.eta"]
    el_ptvarcone= tree["ElectronsAuxDyn.ptvarcone20"]
    el_topoetacone = tree["ElectronsAuxDyn.topoetcone20"]
    el_trk  = tree["ElectronsAuxDyn.trackParticleLinks"]  # still a record, as in your working code
    mu_pt   = tree["MuonsAuxDyn.pt"]
    mu_eta  = tree["MuonsAuxDyn.eta"]
    mu_ptvarcone= tree["MuonsAuxDyn.ptvarcone20"]
    mu_topoetacone = tree["MuonsAuxDyn.topoetcone20"]
    mu_idx_all = tree["MuonsAuxDyn.inDetTrackParticleLink.m_persIndex"]  # flattened leaf

    trk_TTVA_vertices = tree["InDetTrackParticlesAuxDyn.TTVA_AMVFVertices"]
    trk_TTVA_weights  = tree["InDetTrackParticlesAuxDyn.TTVA_AMVFWeights"]

    trk_d0  = tree["InDetTrackParticlesAuxDyn.d0"]
    trk_z0  = tree["InDetTrackParticlesAuxDyn.z0"]
    trk_phi = tree["InDetTrackParticlesAuxDyn.phi"]

    PV_x = tree["PrimaryVerticesAuxDyn.x"]
    PV_y = tree["PrimaryVerticesAuxDyn.y"]
    PV_z = tree["PrimaryVerticesAuxDyn.z"]

    n_events = len(w_nom)
    print(f"  Entries: {n_events}")

    # histograms per sample
    y_title = "event weight: #sigma #int #font[12]{L} #epsilon k W^{i} /#sum W^{i}"
    h_disp_reco = make_disp_hist(
        f"h_{name}_disp_reco",
        f" ;L_{{reco}}^{{d0}} / mm;{y_title}",
    )

    h_disp_truth = None
    if stype == "SIG":
        h_disp_truth = make_disp_hist(
            f"h_{name}_disp_truth",
            f" ;L_{{truth}}^{{xyz}} / mm;{y_title}",
        )

    # for SIG: truth branches
    if stype == "SIG":
        pdgId       = tree["TruthBSMWithDecayParticlesAuxDyn.pdgId"]
        decay_index = tree["TruthBSMWithDecayParticlesAuxDyn.decayVtxLink.m_persIndex"]
        prod_index  = tree["TruthBSMWithDecayParticlesAuxDyn.prodVtxLink.m_persIndex"]
        vx_all      = tree["TruthBSMWithDecayVerticesAuxDyn.x"]
        vy_all      = tree["TruthBSMWithDecayVerticesAuxDyn.y"]
        vz_all      = tree["TruthBSMWithDecayVerticesAuxDyn.z"]

    # event loop
    for evt in range(n_events):
        w = float(w_nom[evt])

        # ------------- truth HNL Lxyz (SIG only) -------------
        if stype == "SIG":
            pdg_evt   = pdgId[evt]
            dec_evt   = decay_index[evt]
            prod_evt  = prod_index[evt]
            vx_evt    = vx_all[evt]
            vy_evt    = vy_all[evt]
            vz_evt    = vz_all[evt]

            # mask HNL
            mask = ak.where(abs(pdg_evt) == 9900012)[0]
            for p_idx in mask:
                d_idx   = int(dec_evt[p_idx])
                p_idx_v = int(prod_evt[p_idx])

                n_vtx = len(vx_evt)
                if d_idx < 0 or d_idx >= n_vtx:
                    continue
                if p_idx_v < 0 or p_idx_v >= n_vtx:
                    continue

                xd = float(vx_evt[d_idx])
                yd = float(vy_evt[d_idx])
                zd = float(vz_evt[d_idx])
                xp = float(vx_evt[p_idx_v])
                yp = float(vy_evt[p_idx_v])
                zp = float(vz_evt[p_idx_v])

                dx = xd - xp
                dy = yd - yp
                dz = zd - zp
                Lxyz = math.sqrt(dx*dx + dy*dy + dz*dz)
                if Lxyz > 1e-10:
                    h_disp_truth.Fill(Lxyz)

        # ------------- reco lepton tracks (e + mu) -------------
        leptons_std = []

        # electrons
        el_pt_evt  = el_pt[evt]
        el_eta_evt = el_eta[evt]
        el_trk_evt = el_trk[evt]
        trk_d0_evt  = trk_d0[evt]
        trk_z0_evt  = trk_z0[evt]
        trk_phi_evt = trk_phi[evt]
        pv_x_evt = PV_x[evt]
        pv_y_evt = PV_y[evt]
        pv_z_evt = PV_z[evt]

        for i in range(len(el_pt_evt)):
            link = el_trk_evt[i]
            if len(link) == 0:
                continue
            trk_idx = link[0]["m_persIndex"]
            trk_idx_int = int(trk_idx)
            if trk_idx_int < 0 or trk_idx_int >= len(trk_d0_evt):
                continue
            # electron > track > vertex > pv_coordinates
            # vert_idx_list = trk_TTVA_vertices[evt][trk_idx_int]
            # weight_list = trk_TTVA_weights[evt][trk_idx_int]
            # if len(vert_idx_list) == 0:
            #     continue
            d0_beam = float(trk_d0_evt[trk_idx_int])
            z0_beam = float(trk_z0_evt[trk_idx_int])
            phi     = float(trk_phi_evt[trk_idx_int])

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

            x = -d0_beam*np.sin(phi)
            y = d0_beam*np.cos(phi)
            leptons_std.append((el_pt[evt][i]*0.001, x, y, z0_beam, "e"))

            # best = ak.argmax(weight_list)
            # best_vertex_link  = vert_idx_list[best]  
            # best_vertex_index = int(best_vertex_link["m_persIndex"])
            # best_weight = float(weight_list[best])
            # if best_vertex_index < 0 or best_vertex_index >= len(pv_x_evt):
            #     continue
            # xPV = float(pv_x_evt[best_vertex_index])
            # yPV = float(pv_y_evt[best_vertex_index])
            # zPV = float(pv_z_evt[best_vertex_index])

            # shift = xPV * math.sin(phi) - yPV * math.cos(phi)
            # d0_PV = d0_beam - shift
            # z0_PV = z0_beam - zPV
            # x_rel = -d0_PV * math.sin(phi)
            # y_rel =  d0_PV * math.cos(phi)
            # z_rel =  z0_PV
            # leptons_std.append((float(el_pt_evt[i]), x_rel, y_rel, z_rel, "e"))

        # muons
        mu_pt_evt  = mu_pt[evt]
        mu_idx_evt = mu_idx_all[evt]  # per-event indices

        for i in range(len(mu_idx_evt)):
            trk_idx_int = int(mu_idx_evt[i])
            if trk_idx_int < 0 or trk_idx_int >= len(trk_d0_evt):
                continue
            # muon > track > vertex > pv_coordinates
            vert_idx_list = trk_TTVA_vertices[evt][trk_idx_int]
            weight_list = trk_TTVA_weights[evt][trk_idx_int]
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

            x = -d0_beam*np.sin(phi)
            y = d0_beam*np.cos(phi)
            leptons_std.append((mu_pt[evt][i]*0.001, x, y, z0_beam, "m"))
            # best = ak.argmax(weight_list)
            # best_vertex_link  = vert_idx_list[best]  
            # best_vertex_index = int(best_vertex_link["m_persIndex"])
            # best_weight = float(weight_list[best])
            # if best_vertex_index < 0 or best_vertex_index >= len(pv_x_evt):
            #     continue
            # xPV = float(pv_x_evt[best_vertex_index])
            # yPV = float(pv_y_evt[best_vertex_index])
            # zPV = float(pv_z_evt[best_vertex_index])

            # shift = xPV * math.sin(phi) - yPV * math.cos(phi)
            # d0_PV = d0_beam - shift
            # z0_PV = z0_beam - zPV
            # x_rel = -d0_PV * math.sin(phi)
            # y_rel =  d0_PV * math.cos(phi)
            # z_rel =  z0_PV
            # leptons_std.append((float(mu_pt_evt[i]), x_rel, y_rel, z_rel, "m"))

        # require at least 2 leptons; use 2 leading pT
        if len(leptons_std) >= 2:
            leptons_sorted = sorted(leptons_std, key=lambda x: x[0], reverse=True)
            l1 = leptons_sorted[0]
            l2 = leptons_sorted[1]
            dx = l1[1] - l2[1]
            dy = l1[2] - l2[2]
            dz = l1[3] - l2[3]
            dL = math.sqrt(dx*dx + dy*dy + dz*dz)
            if dL > 0:
                h_disp_reco.Fill(dL)

    print(f"  total_sumw_loop = {total_sumw_loop}")
    total_sumw_used = total_sumw_loop
    if use_metadata_totalW and metadata is not None:
        total_sumw_meta = get_total_sumw_from_metadata(sample, metadata)
        if total_sumw_meta > 0.0:
            total_sumw_used = total_sumw_meta
            print(f"  Using metadata total_sumw (sum AllExecutedEvents.W) = {total_sumw_meta}")
        else:
            print(f"  WARNING: metadata total_sumw <= 0 for {name}, falling back to loop total_sumw.")
    return h_disp_reco, h_disp_truth, total_sumw_loop, total_sumw_used


# ------------------------------------------------------------
# Process all samples, compute super_factors, scale histos
# ------------------------------------------------------------
histos_disp_reco = {}
histos_disp_truth = {}
super_factors = {}

for s in samples:
    h_reco, h_truth, sumw_loop, sumw_used = process_sample(s, metadata=metadata if use_metadata_totalW else None, use_metadata_totalW=use_metadata_totalW)

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

    # scale histos
    h_reco.Scale(sf)
    histos_disp_reco[s["name"]] = h_reco

    if h_truth is not None:
        h_truth.Scale(sf)
        histos_disp_truth[s["name"]] = h_truth


# ------------------------------------------------------------
# Build group histos for displacement
# ------------------------------------------------------------
group_histos_reco = {}
for gname, members in group_definitions.items():
    # find a template
    template = None
    for s in members:
        if s in histos_disp_reco:
            template = histos_disp_reco[s]
            break
    if template is None:
        continue
    h_group = template.Clone(f"h_{gname}_disp_reco")
    h_group.Reset("ICES")
    for s in members:
        if s in histos_disp_reco:
            h_group.Add(histos_disp_reco[s])
    group_histos_reco[gname] = h_group

# pooled S (SIG) and B (all BKG_*)
h_pooled_S = histos_disp_reco["WR5000_N50_2L2J"].Clone("h_pooled_S_disp")
h_pooled_S.Reset("ICES")
h_pooled_B = histos_disp_reco["WR5000_N50_2L2J"].Clone("h_pooled_B_disp")
h_pooled_B.Reset("ICES")

for s in samples:
    h = histos_disp_reco[s["name"]]
    if s["type"] == "SIG":
        h_pooled_S.Add(h)
    elif s["type"].startswith("BKG"):
        h_pooled_B.Add(h)

# significance histogram
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
# Plotting: displacement with significance (copy C++ style)
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

    # define drawing order like C++: backgrounds then signal
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

    # legend position 
    x1, y1, x2, y2 = 0.37, 0.68, 0.52, 0.90   #0.37, 0.68, 0.52, 0.90
    leg = ROOT.TLegend(x1, y1, x2, y2)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.04)

    # build THStack for backgrounds
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

    # choose axis histogram (for fonts/titles)
    h_axis = None
    for gname in bkg_order:
        if gname in group_histos_reco:
            h_axis = group_histos_reco[gname]
            break
    if h_axis is None and h_sig is not None:
        h_axis = h_sig

    # draw backgrounds stack
    hs_bkg.Draw("HIST")
    if h_axis is not None:
        # reuse axis settings from C++ (hard-coded for aesthetics)
        hs_bkg.GetXaxis().SetTitle(h_axis.GetXaxis().GetTitle())
        hs_bkg.GetYaxis().SetTitle(h_axis.GetYaxis().GetTitle())
        hs_bkg.GetXaxis().SetRangeUser(xmin_vis, xmax_vis)
        hs_bkg.SetMaximum(1.5e9)  # as in C++, tune if needed
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

    # signal overlay
    if h_sig is not None:
        h_sig.Draw("HIST SAME")

    leg.Draw()

    # ATLAS label
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

    # bottom pad: significance (unchanged except that we now match C++ style more closely)
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
# Debug: truth vs reco for signal only
# ------------------------------------------------------------
def plot_signal_truth_vs_reco():
    sig_name = "WR5000_N50_2L2J"
    if sig_name not in histos_disp_truth:
        print("Signal reco or truth hist missing; skip truth-vs-reco plot.")
        return

    # hReco  = histos_disp_reco[sig_name]
    hTruth = histos_disp_truth[sig_name]
    global h_pooled_S, h_pooled_B
    hRecoS = h_pooled_S
    hRecoB = h_pooled_B
    
    if hRecoS.Integral() > 0:
        hRecoS.Scale(1.0 / hRecoS.Integral())

    if hRecoB.Integral() > 0:
        hRecoB.Scale(1.0 / hRecoB.Integral())

    if hTruth.Integral() > 0:
        hTruth.Scale(1.0 / hTruth.Integral())

    # exponential model using mean from truth (or fixed)
    L_mean = 0.135662  
    Ntot   = 1
    hExp   = hRecoS.Clone("h_disp_exponential")
    hExp.Reset("ICES")

    nb = hExp.GetNbinsX()
    for ib in range(1, nb + 1):
        x     = hExp.GetBinCenter(ib)
        width = hExp.GetBinWidth(ib)
        if L_mean > 0:
            fL = (1.0 / L_mean) * math.exp(-x / L_mean)
            val = Ntot * fL * width
            hExp.SetBinContent(ib, val)

    # style
    hRecoS.SetLineColor(ROOT.kRed + 1)
    hRecoS.SetLineWidth(2)
    hRecoS.SetFillStyle(0)
    hRecoS.SetMarkerStyle(20)
    hRecoS.SetMarkerColor(ROOT.kRed + 1)
    hRecoS.SetMarkerSize(0.8)

    hRecoB.SetLineColor(ROOT.kGreen + 2)
    hRecoB.SetLineWidth(2)
    hRecoB.SetFillStyle(0)
    hRecoB.SetMarkerStyle(22)
    hRecoB.SetMarkerColor(ROOT.kGreen + 2)
    hRecoB.SetMarkerSize(0.8)

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

    maxData = max(hRecoS.GetMaximum(), hTruth.GetMaximum(), hRecoB.GetMaximum())
    maxExp  = hExp.GetMaximum()
    maxY    = max(maxData, maxExp)
    hRecoS.SetMinimum(0.0)
    hRecoS.SetMaximum(1.2 * maxY)

    hRecoS.GetXaxis().SetTitle("L / mm")
    hRecoS.GetYaxis().SetTitle("event-rate (normalized)")
    hRecoS.GetXaxis().SetTitleSize(0.045)
    hRecoS.GetXaxis().SetLabelSize(0.04)
    hRecoS.GetYaxis().SetTitleSize(0.045)
    hRecoS.GetYaxis().SetLabelSize(0.04)

    hRecoS.Draw("HIST E1")
    hRecoB.Draw("HIST E1 SAME")
    hTruth.Draw("HIST E1 SAME")
    hExp.Draw("HIST SAME")

    leg = ROOT.TLegend(0.4, 0.70, 0.85, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.04)
    leg.AddEntry(hRecoS, "L_{reco}^{d0} SIG: mN=50GeV, mWR=10TeV", "l")
    leg.AddEntry(hRecoB, "L_{reco}^{d0} BKG: Zmumu/ttbar", "l")
    leg.AddEntry(hTruth, "L_{truth}^{xyz} BSMtruth: prodV-decV", "l")
    leg.AddEntry(hExp,   "P(L; #it{L_{0}} = 0.135662 mm)", "l")   #=#frac{1}{#mu}exp(-L/#it{L_{0}})
    leg.Draw()

    outname = os.path.join(output_dir, "debug_disp_reco_d0_vs_truth_Lxyz_signal_uproot.png")
    cL.SaveAs(outname)
    print(f"Saved {outname}")


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------
if __name__ == "__main__":
    print("User Test 049:")
    plot_disp_with_significance()
    plot_signal_truth_vs_reco()
    print(f"Analysis successfully completed. Results saved to {output_dir}")
