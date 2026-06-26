import uproot 
import ROOT
import awkward as ak
import numpy as np
from collections import defaultdict
import math
import os

print('uproot:', uproot.__version__, '| awkward:', ak.__version__)
print("User Test 027")

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

# Scan Configuration
# ---------------------------------------------------------------------
max_events = None      # set 100_000 to set max limit < N_events
step_size  = 50_000    # events per uproot.iterate basket

# FILE HANDLING
#----------------------------------------------------------------------------
output_dir = "/home/alexanum/WORKSPACE/testing/MadGraph5/analysis/SIGvsBKG_RECO/Figures/Multiplicities67"
os.makedirs(output_dir, exist_ok=True)

# file_path = "/storage/alexanum/DAOD_LLP_WORKFLOW/ROOT_files/DAOD_LLP1.DAOD_LLP1.mc23_13p6TeV.MGPy8EG_LRSM_WR_HNL_5000_30_2L2J_MadSpin.DAOD_LLP1.pool.root"
# file_path = "/home/alexanum/WORKSPACE/RUN_CONTROLS/ATLAS_jobdir/ROOT_files/DAOD_LLP1.mc23_13p6TeV.150001.MGPy8EG_LRSM_WR_HNL_5000_50_2L2J_MadSpin.DAOD_LLP1.pool.root"
# file_path = "/storage/alexanum/DAOD_LLP_WORKFLOW/ROOT_files/DAOD_LLP1.DAOD_LLP1.mc23_13p6TeV.MGPy8EG_LRSM_WR_HNL_5000_70_2L2J_MadSpin.DAOD_LLP1.pool.root"
file_path = "/storage/alexanum/DAOD_LLP_WORKFLOW/ROOT_files/DAOD_LLP1.DAOD_LLP1.mc23_13p6TeV.MGPy8EG_LRSM_WR_HNL_5000_500_2L2J_MadSpin.DAOD_LLP1.pool.root"
# file_path = "/storage/alexanum/DAOD_LLP_WORKFLOW/ROOT_files/DAOD_LLP1.DAOD_LLP1.mc23_13p6TeV.MGPy8EG_LRSM_WR_HNL_5000_500_2L2J_MadSpin.DAOD_LLP1.pool.root"
# file_path = "/storage/alexanum/DAOD_LLP_WORKFLOW/ROOT_files/DAOD_LLP1.DAOD_LLP1.mc23_13p6TeV.MGPy8EG_LRSM_WR_HNL_10000_50_2L2J_MadSpin.DAOD_LLP1.pool.root"

# DEFINE GLOBAL SELECTION CUTS
#----------------------------------------------------------------------------
lep_pt_min   = 10.0 #GeV  # choose < mN/3
lep_eta_max  = 2.5
jet_pt_min   = 10.0 #GeV
jet_eta_max  = 2.5

# isolation cuts (relative)
lep_rel_ptvarcone20_max  = 0.15
lep_rel_topoetcone20_max = 0.15

# jet-quality / pileup rejection 
jvt_min = 0.5 
jet_width_cut = 0.05   #0.05
jet_ntrk_cut = 10      #10
dr_jj_cut = 1
# ---------------------------------------------------------------------
# Functions: Kinematics & Awk-Type Conversion
# ---------------------------------------------------------------------
def calculate_pt(px, py):
    """Calculate transverse momentum in GeV from MeV components"""
    return np.sqrt(px**2 + py**2) / 1000.0

def deltaR(eta1, phi1, eta2, phi2):
    dphi = abs(phi1 - phi2)
    if dphi > math.pi:
        dphi = 2*math.pi - dphi
    deta = eta1 - eta2
    return math.sqrt(deta*deta + dphi*dphi)
    # return abs(deta)

def calculate_mass(px, py, pz, e):
    """Calculate invariant mass GeV from MeV components"""
    return np.sqrt(np.maximum(0, e**2 - px**2 - py**2 - pz**2)) / 1000.0

def sort_by_pt(objects):
    """Sort list of dicts by pT (descending)"""
    return sorted(objects, key=lambda x: x['pt'], reverse=True)

def ak_scalar(x):
    """Convert Awkward/NumPy 0-d or length-1 array-like into a Python float."""
    arr = ak.to_numpy(x)
    arr = np.asarray(arr)
    if arr.shape == ():        # 0-d array (scalar)
        return float(arr)
    return float(arr.reshape(-1)[0])

def calculate_eta_phi(px, py, pz):
    """Return (eta, phi) from MeV 3-momentum components."""
    px = np.asarray(px, dtype=float)
    py = np.asarray(py, dtype=float)
    pz = np.asarray(pz, dtype=float)
    phi = np.arctan2(py, px)
    pT = np.sqrt(px**2 + py**2)
    eta = np.zeros_like(pT)
    mask = pT > 0
    eta[mask] = np.arcsinh(pz[mask] / pT[mask])
    return eta, phi

# ---------------------------------------------------------------------
# HISTOGRAMS (GLOBAL)
# ---------------------------------------------------------------------
# multiplicities (reconstructed before vs after basic-cuts)
h_lep_before = ROOT.TH1F("h_lep_before", " ; ; ", 20, 0, 10)
h_lep_after  = ROOT.TH1F("h_lep_after",  " ; ; ", 20, 0, 10)
h_jet_before = ROOT.TH1F("h_jet_before", " ; ; ", 20, 0, 10)
h_jet_after  = ROOT.TH1F("h_jet_after",  " ; ; ", 20, 0, 10)

# multiplicity-counters
lep_before_counts = defaultdict(int)
lep_after_counts  = defaultdict(int)
jet_before_counts = defaultdict(int)
jet_after_counts  = defaultdict(int)
Tlep_counts       = defaultdict(int)
Tquark_counts     = defaultdict(int)

# pT distributions: reco vs truth 
h_lep_lead_reco  = ROOT.TH1F("h_lep_lead_reco",  "Dilepton (leading+subleading leptons), post- Kin/Iso cuts;p_{T} [GeV];Events", 100, 0, 4000)
h_lep_lead_truth = ROOT.TH1F("h_lep_lead_truth", "Leading leptons (truth);p_{T} [GeV];Events", 100, 0, 4000)
h_lep_other      = ROOT.TH1F("h_lep_other",      "Other leptons (reco);p_{T} [GeV];Events", 100, 0, 1000)

h_jet_lead_reco  = ROOT.TH1F("h_jet_lead_reco",  "Dijet (leading+subleading small-R jets), post- Kin/JVT cuts;p_{T} [GeV];Events", 100, 0, 3000)
h_jet_lead_reco2  = ROOT.TH1F("h_jet_lead_reco2",  "Dijet (leading+subleading small-R jets), post- Kin/JVT cuts;p_{T} [GeV];Events", 100, 0, 3000)
h_jet_lead_reco3  = ROOT.TH1F("h_jet_lead_reco3",  "Dijet (leading+subleading small-R jets), post- Kin/JVT cuts;p_{T} [GeV];Events", 100, 0, 3000)
h_jet_lead_reco4  = ROOT.TH1F("h_jet_lead_reco4",  "Dijet (leading+subleading small-R jets), post- Kin/JVT cuts;p_{T} [GeV];Events", 100, 0, 3000)
h_jet_lead_truth = ROOT.TH1F("h_jet_lead_truth", "Reconstrcuted Leading Dilepton; p_{T} [GeV];Events", 100, 0, 3000)
h_jet_other      = ROOT.TH1F("h_jet_other",      "Other jets (reco);p_{T} [GeV];Events", 100, 0, 1200)
h_mjj_truth = ROOT.TH1F(
    "h_mjj_truth",
    "Inv. mass of leading jet-pair;m_{inv.} [GeV];Events",
    75, 0, 5000
)
h_mjj_reco = ROOT.TH1F(
    "h_mjj_reco",
    "Inv. mass of leading jet-pair with #Delta R<2;m_{jj} [GeV];Events",
    75, 0, 5000    #USE: 75, 0, 5000 to zoom in on small 50,10000 region
)
h_mjj_reco2 = ROOT.TH1F(
    "h_mjj_reco2",
    "Inv. mass of leading jet-pair with #Delta R<2;m_{jj} [GeV];Events",
    75, 0, 5000   #75, 0, 5000
)
h_mjj_reco3 = ROOT.TH1F(
    "h_mjj_reco3",
    "Inv. mass of leading jet-pair with #Delta R<2;m_{jj} [GeV];Events",
    75, 0, 5000
)
h_mjj_truthjet = ROOT.TH1F(
    "h_mjj_truthjet",
    "Truth jet dijet invariant mass;m_{jj}^{truth\ jet} [GeV];Events",
    75, 0, 5000
)

h_dR_jet_reco  = ROOT.TH1F(
    "h_dR_jet_reco",
    "Reco: #DeltaR(j_{1}, j_{2}) (leading/subleading small-R jets);#DeltaR(j_{1}, j_{2});Events",
    50, 0, 5
)
h_dR_jet_truth = ROOT.TH1F(
    "h_dR_jet_truth",
    "Truth: #DeltaR(j_{1}, j_{2}) (leading/subleading AK4 truth jets);#DeltaR(j_{1}, j_{2});Events",
    50, 0, 5
)

# isolation distributions
h_lep_rel_iso_ptvar = ROOT.TH1F("h_lep_rel_iso_ptvar", "Lepton rel ptvarcone20;ptvarcone20/p_{T};Leptons", 50, 0, 0.5)
h_lep_rel_iso_topo  = ROOT.TH1F("h_lep_rel_iso_topo",  "Lepton rel topoetcone20;topoetcone20/p_{T};Leptons", 50, 0, 0.5)

# jet shape distributions for leading jet
h_jet_lead_mass  = ROOT.TH1F("h_jet_lead_mass",  "Leading small-R jet mass;m_{jet} [GeV];Events", 50, 0, 500)
h_jet_lead_width = ROOT.TH1F("h_jet_lead_width", "Leading small-R jet width;Width;Events", 50, 0, 0.5)
h_jet_lead_ntrk  = ROOT.TH1F("h_jet_lead_ntrk",  "Leading small-R jet NumTrkPt500;NumTrkPt500;Events", 50, 0, 50)

# large-R jet multiplicity
h_largeR_jet_mult = ROOT.TH1F("h_largeR_jet_mult", " ; ; ", 10, 0, 5)

# HNL gamma vs jet multiplicity
h2_hnlgamma_vs_nlep = ROOT.TH2F(
    "h2_hnlgamma_vs_nlep",
    "Correlation: HNL-boost to Lepton Isolation (mN=500GeV);#gamma_{N} ; ",
    109, 0, 109, 2, 0.5, 2.5
)

h_lep_truth_rank = [
    ROOT.TH1F("h_lep_truth_lead",   "Truth leptons; p_{T} [GeV];Events", 100, 0, 4000),
    ROOT.TH1F("h_lep_truth_sub",    "Truth leptons; p_{T} [GeV];Events", 100, 0, 4000),
    ROOT.TH1F("h_lep_truth_subsub", "Truth leptons; p_{T} [GeV];Events", 100, 0, 4000),
]

h_lep_reco_rank = [
    ROOT.TH1F("h_lep_reco_lead",   "Reco leptons; p_{T} [GeV];Events", 100, 0, 4000),
    ROOT.TH1F("h_lep_reco_sub",    "Reco leptons; p_{T} [GeV];Events", 100, 0, 4000),
    ROOT.TH1F("h_lep_reco_subsub", "Reco leptons; p_{T} [GeV];Events", 100, 0, 4000),
]

# Jets: truth (from quarks) and reco
h_jet_truth_rank = [
    ROOT.TH1F("h_jet_truth_lead",   "Truth jets (quark-based); p_{T} [GeV];Events", 100, 0, 3000),
    ROOT.TH1F("h_jet_truth_sub",    "Truth jets (quark-based); p_{T} [GeV];Events", 100, 0, 3000),
    ROOT.TH1F("h_jet_truth_subsub", "Truth jets (quark-based); p_{T} [GeV];Events", 100, 0, 3000),
]

h_jet_reco_rank = [
    ROOT.TH1F("h_jet_reco_lead",   "Reco jets; p_{T} [GeV];Events", 100, 0, 3000),
    ROOT.TH1F("h_jet_reco_sub",    "Reco jets; p_{T} [GeV];Events", 100, 0, 3000),
    ROOT.TH1F("h_jet_reco_subsub", "Reco jets; p_{T} [GeV];Events", 100, 0, 3000),
]

#119
# ---------------------------------------------------------------------
# EVENT LOOP WITH uproot.iterate
# ---------------------------------------------------------------------
print("Starting event loop...")

branches = [
    # electrons
    "ElectronsAuxDyn.pt",
    "ElectronsAuxDyn.eta",
    "ElectronsAuxDyn.m",
    "ElectronsAuxDyn.ptvarcone20",
    "ElectronsAuxDyn.topoetcone20",
    # muons
    "MuonsAuxDyn.pt",
    "MuonsAuxDyn.eta",
    "MuonsAuxDyn.ptvarcone20",
    "MuonsAuxDyn.topoetcone20",
    # small-R jets
    "AntiKt4EMPFlowJetsAuxDyn.JetConstitScaleMomentum_pt",
    "AntiKt4EMPFlowJetsAuxDyn.JetConstitScaleMomentum_eta",
    "AntiKt4EMPFlowJetsAuxDyn.JetConstitScaleMomentum_phi",
    "AntiKt4EMPFlowJetsAuxDyn.JetConstitScaleMomentum_m",
    "AntiKt4EMPFlowJetsAuxDyn.Jvt",
    "AntiKt4EMPFlowJetsAuxDyn.NumTrkPt500",
    "AntiKt4EMPFlowJetsAuxDyn.Width",
    # truth jets
    "AntiKt4TruthJetsAuxDyn.JetConstitScaleMomentum_pt", 
    "AntiKt4TruthJetsAuxDyn.JetConstitScaleMomentum_eta", 
    "AntiKt4TruthJetsAuxDyn.JetConstitScaleMomentum_phi", 
    "AntiKt4TruthJetsAuxDyn.JetConstitScaleMomentum_m",
    # large-R jets
    "AntiKt10EMTopoRCJetsAuxDyn.JetConstitScaleMomentum_pt",
    "AntiKt10EMTopoRCJetsAuxDyn.JetConstitScaleMomentum_eta",
    "AntiKt10EMTopoRCJetsAuxDyn.JetConstitScaleMomentum_phi",
    "AntiKt10EMTopoRCJetsAuxDyn.JetConstitScaleMomentum_m",
    # truth particles
    "TruthParticlesAuxDyn.pdgId",
    "TruthParticlesAuxDyn.px",
    "TruthParticlesAuxDyn.py",
    "TruthParticlesAuxDyn.pz",
    "TruthParticlesAuxDyn.e",
    # BSM / HNL
    "TruthBSMWithDecayParticlesAuxDyn.pdgId",
    "TruthBSMWithDecayParticlesAuxDyn.px",
    "TruthBSMWithDecayParticlesAuxDyn.py",
    "TruthBSMWithDecayParticlesAuxDyn.pz",
    "TruthBSMWithDecayParticlesAuxDyn.e",
]

total_events = 0
chunk_idx    = 0

for arrays in uproot.iterate(
    file_path + ":CollectionTree",
    branches,
    step_size=step_size,
    library="ak"
):
    chunk_idx += 1

    # Limit total events if requested
    n_events_chunk = len(arrays["ElectronsAuxDyn.pt"])
    if max_events is not None:
        if total_events >= max_events:
            break
        if total_events + n_events_chunk > max_events:
            n_events_chunk = max_events - total_events
            arrays = arrays[:n_events_chunk]

    total_events += n_events_chunk
    print(f"Chunk {chunk_idx}, events in chunk: {n_events_chunk}, total: {total_events}")

    # ---------------------- convenience aliases ----------------------
    e_pt           = arrays["ElectronsAuxDyn.pt"]
    e_eta          = arrays["ElectronsAuxDyn.eta"]
    e_m            = arrays["ElectronsAuxDyn.m"]
    e_ptvarcone20  = arrays["ElectronsAuxDyn.ptvarcone20"]
    e_topoetcone20 = arrays["ElectronsAuxDyn.topoetcone20"]

    m_pt           = arrays["MuonsAuxDyn.pt"]
    m_eta          = arrays["MuonsAuxDyn.eta"]
    m_ptvarcone20  = arrays["MuonsAuxDyn.ptvarcone20"]
    m_topoetcone20 = arrays["MuonsAuxDyn.topoetcone20"]

    j_pt          = arrays["AntiKt4EMPFlowJetsAuxDyn.JetConstitScaleMomentum_pt"]
    j_eta         = arrays["AntiKt4EMPFlowJetsAuxDyn.JetConstitScaleMomentum_eta"]
    j_phi         = arrays["AntiKt4EMPFlowJetsAuxDyn.JetConstitScaleMomentum_phi"]
    j_m           = arrays["AntiKt4EMPFlowJetsAuxDyn.JetConstitScaleMomentum_m"]
    j_Jvt         = arrays["AntiKt4EMPFlowJetsAuxDyn.Jvt"]
    j_NumTrkPt500 = arrays["AntiKt4EMPFlowJetsAuxDyn.NumTrkPt500"]
    j_Width       = arrays["AntiKt4EMPFlowJetsAuxDyn.Width"]

    tj_pt  = arrays["AntiKt4TruthJetsAuxDyn.JetConstitScaleMomentum_pt"]
    tj_eta = arrays["AntiKt4TruthJetsAuxDyn.JetConstitScaleMomentum_eta"]
    tj_phi = arrays["AntiKt4TruthJetsAuxDyn.JetConstitScaleMomentum_phi"]
    tj_m   = arrays["AntiKt4TruthJetsAuxDyn.JetConstitScaleMomentum_m"]

    J10_pt  = arrays["AntiKt10EMTopoRCJetsAuxDyn.JetConstitScaleMomentum_pt"]
    J10_eta = arrays["AntiKt10EMTopoRCJetsAuxDyn.JetConstitScaleMomentum_eta"]
    J10_m   = arrays["AntiKt10EMTopoRCJetsAuxDyn.JetConstitScaleMomentum_m"]

    truth_pdg = arrays["TruthParticlesAuxDyn.pdgId"]
    truth_px  = arrays["TruthParticlesAuxDyn.px"]
    truth_py  = arrays["TruthParticlesAuxDyn.py"]
    truth_pz  = arrays["TruthParticlesAuxDyn.pz"]
    truth_e   = arrays["TruthParticlesAuxDyn.e"]

    bsm_pdg = arrays["TruthBSMWithDecayParticlesAuxDyn.pdgId"]
    bsm_px  = arrays["TruthBSMWithDecayParticlesAuxDyn.px"]
    bsm_py  = arrays["TruthBSMWithDecayParticlesAuxDyn.py"]
    bsm_pz  = arrays["TruthBSMWithDecayParticlesAuxDyn.pz"]
    bsm_e   = arrays["TruthBSMWithDecayParticlesAuxDyn.e"]

    # -----------------------------------------------------------------
    # EVENT LOOP over this chunk
    # -----------------------------------------------------------------
    for evt in range(n_events_chunk):
        # print less often
        if evt % 5000 == 0:
            print(f"  processing event {evt}/{n_events_chunk} in chunk {chunk_idx}")

        # ======================= RECO LEPTONS ========================
        leptons_all = []

        # ELECTRONS
        e_pt_evt   = e_pt[evt]
        e_eta_evt  = e_eta[evt]
        e_m_evt    = e_m[evt]
        e_ptv_evt  = e_ptvarcone20[evt]
        e_topo_evt = e_topoetcone20[evt]

        # precompute arrays (vectorized) then only build dicts once per electron
        if len(e_pt_evt) > 0:
            pt_e   = ak.to_numpy(e_pt_evt) / 1000.0
            eta_e  = ak.to_numpy(e_eta_evt)
            m_e    = ak.to_numpy(e_m_evt) / 1000.0
            ptv_e  = ak.to_numpy(e_ptv_evt)
            topo_e = ak.to_numpy(e_topo_evt)

            valid_mask_e = pt_e > 0
            pt_e   = pt_e[valid_mask_e]
            eta_e  = eta_e[valid_mask_e]
            m_e    = m_e[valid_mask_e]
            ptv_e  = ptv_e[valid_mask_e]
            topo_e = topo_e[valid_mask_e]

            if pt_e.size > 0:
                rel_ptvar_e = np.where(pt_e > 0, ptv_e / (pt_e * 1000.0), 999)
                rel_topo_e  = np.where(pt_e > 0, topo_e / (pt_e * 1000.0), 999)

                for pt, eta, mll, rpv, rtop in zip(pt_e, eta_e, m_e, rel_ptvar_e, rel_topo_e):
                    leptons_all.append({
                        "pt": float(pt),
                        "eta": float(eta),
                        "m": float(mll),
                        "rel_ptvar": float(rpv),
                        "rel_topo": float(rtop),
                        "type": "e",
                    })
                    h_lep_rel_iso_ptvar.Fill(rpv)
                    h_lep_rel_iso_topo.Fill(rtop)

        # MUONS
        m_pt_evt   = m_pt[evt]
        m_eta_evt  = m_eta[evt]
        m_ptv_evt  = m_ptvarcone20[evt]
        m_topo_evt = m_topoetcone20[evt]

        if len(m_pt_evt) > 0:
            pt_m   = ak.to_numpy(m_pt_evt) / 1000.0
            eta_m  = ak.to_numpy(m_eta_evt)
            ptv_m  = ak.to_numpy(m_ptv_evt)
            topo_m = ak.to_numpy(m_topo_evt)

            valid_mask_m = pt_m > 0
            pt_m   = pt_m[valid_mask_m]
            eta_m  = eta_m[valid_mask_m]
            ptv_m  = ptv_m[valid_mask_m]
            topo_m = topo_m[valid_mask_m]

            if pt_m.size > 0:
                m_mu = 0.10566  # GeV
                rel_ptvar_m = np.where(pt_m > 0, ptv_m / (pt_m * 1000.0), 999)
                rel_topo_m  = np.where(pt_m > 0, topo_m / (pt_m * 1000.0), 999)

                for pt, eta, rpv, rtop in zip(pt_m, eta_m, rel_ptvar_m, rel_topo_m):
                    leptons_all.append({
                        "pt": float(pt),
                        "eta": float(eta),
                        "m": m_mu,
                        "rel_ptvar": float(rpv),
                        "rel_topo": float(rtop),
                        "type": "mu",
                    })
                    h_lep_rel_iso_ptvar.Fill(rpv)
                    h_lep_rel_iso_topo.Fill(rtop)

        # count leptons before cuts
        n_lep_before = len(leptons_all)
        lep_before_counts[n_lep_before] += 1
        h_lep_before.Fill(n_lep_before)

        # apply kinematic + isolation cuts
        leptons_sel = [
            l for l in leptons_all
            if (
                l["pt"] > lep_pt_min
                and abs(l["eta"]) < lep_eta_max
                and l["rel_ptvar"] < lep_rel_ptvarcone20_max
                and l["rel_topo"]  < lep_rel_topoetcone20_max
            )
        ]
        leptons_sel = sort_by_pt(leptons_sel)

        n_lep_after = len(leptons_sel)
        lep_after_counts[n_lep_after] += 1
        h_lep_after.Fill(n_lep_after)

        # ======================= RECO SMALL-R JETS ========================
        jets_all = []

        j_pt_evt   = j_pt[evt]
        j_eta_evt  = j_eta[evt]
        j_phi_evt  = j_phi[evt]
        j_m_evt    = j_m[evt]
        j_Jvt_evt  = j_Jvt[evt]
        j_Ntrk_evt = j_NumTrkPt500[evt]
        j_W_evt    = j_Width[evt]

        n_jets_evt = len(j_pt_evt)
        if n_jets_evt > 0:
            # pt/eta/m are simple numbers in MeV; convert to NumPy vectors for speed
            pt_j  = ak.to_numpy(j_pt_evt)  / 1000.0
            eta_j = ak.to_numpy(j_eta_evt)
            phi_j = ak.to_numpy(j_phi_evt)  
            m_j   = ak.to_numpy(j_m_evt)   / 1000.0

            # Jvt, NumTrkPt500, Width may still be scalar-like or length-1 arrays per jet
            for i in range(n_jets_evt):
                pt   = float(pt_j[i])
                eta  = float(eta_j[i])
                phi  = float(phi_j[i]) 
                mass = float(m_j[i])

                jvt_val   = ak_scalar(j_Jvt_evt[i])
                ntrk_val  = ak_scalar(j_Ntrk_evt[i])
                width_val = ak_scalar(j_W_evt[i])

                jets_all.append({
                    "pt": pt,
                    "eta": eta,
                    "phi": phi,
                    "m": mass,
                    "Jvt": jvt_val,
                    "NumTrkPt500": ntrk_val,
                    "Width": width_val,
                })


        n_jet_before = len(jets_all)
        jet_before_counts[n_jet_before] += 1
        h_jet_before.Fill(n_jet_before)

        # apply kinematic + simple JVT cut
        jets_sel = [
            j for j in jets_all
            if (
                j["pt"] > jet_pt_min
                and abs(j["eta"]) < jet_eta_max
                # and j["Jvt"] > jvt_min
                # and j["Width"] < jet_width_cut
                # and j["NumTrkPt500"] < jet_ntrk_cut
            )
        ]
        jets_sel2 = [
            j for j in jets_sel
            if (
                j["Jvt"] > jvt_min
                # and j["Width"] < jet_width_cut
                and j["NumTrkPt500"] < jet_ntrk_cut
            )
        ]
        jets_sel = sort_by_pt(jets_sel)
        jets_sel2 = sort_by_pt(jets_sel2)
        if len(jets_sel) >= 2:
            j0 = jets_sel[0]
            j1 = jets_sel[1]
            # build 4-vectors in GeV using ROOT.TLorentzVector
            vj0 = ROOT.TLorentzVector()
            vj1 = ROOT.TLorentzVector()
            vj0.SetPtEtaPhiM(j0["pt"], j0["eta"], j0["phi"], j0["m"])
            vj1.SetPtEtaPhiM(j1["pt"], j1["eta"], j1["phi"], j1["m"])
            mjj_reco = (vj0 + vj1).M()   # GeV
            if deltaR(j0["eta"], j0["phi"], j1["eta"], j1["phi"]) <= dr_jj_cut:
                h_jet_lead_reco4.Fill(j0["pt"])
                h_jet_lead_reco4.Fill(j1["pt"])
            h_mjj_reco2.Fill(mjj_reco)
        
        if len(jets_sel2) >= 2:
            j0 = jets_sel2[0]
            j1 = jets_sel2[1]
            # build 4-vectors in GeV using ROOT.TLorentzVector
            vj0 = ROOT.TLorentzVector()
            vj1 = ROOT.TLorentzVector()
            vj0.SetPtEtaPhiM(j0["pt"], j0["eta"], j0["phi"], j0["m"])
            vj1.SetPtEtaPhiM(j1["pt"], j1["eta"], j1["phi"], j1["m"])
            mjj_reco = (vj0 + vj1).M()   # GeV
            if deltaR(j0["eta"], j0["phi"], j1["eta"], j1["phi"]) <= dr_jj_cut:
                h_mjj_reco.Fill(mjj_reco)
                h_jet_lead_reco2.Fill(j0["pt"])
                h_jet_lead_reco2.Fill(j1["pt"])
            h_jet_lead_reco3.Fill(j0["pt"])
            h_jet_lead_reco3.Fill(j1["pt"])
            h_mjj_reco3.Fill(mjj_reco)


        n_jet_after = len(jets_sel)
        jet_after_counts[n_jet_after] += 1
        h_jet_after.Fill(n_jet_after)

        # ======================= TRUTH JETS (AK4) m_jj =====================
        truth_jets = []
        tj_pt_evt  = tj_pt[evt]
        tj_eta_evt = tj_eta[evt]
        tj_phi_evt = tj_phi[evt]
        tj_m_evt   = tj_m[evt]

        if len(tj_pt_evt) > 0:
            pt_tj  = ak.to_numpy(tj_pt_evt)  / 1000.0  # GeV
            eta_tj = ak.to_numpy(tj_eta_evt)
            phi_tj = ak.to_numpy(tj_phi_evt)
            m_tj   = ak.to_numpy(tj_m_evt)   / 1000.0  # GeV

            for pt, eta, phi, mass in zip(pt_tj, eta_tj, phi_tj, m_tj):
                truth_jets.append({
                    "pt": float(pt),
                    "eta": float(eta),
                    "phi": float(phi),
                    "m": float(mass),
                })

            # sort by pt and take two leading truth jets
            truth_jets = sorted(truth_jets, key=lambda j: j["pt"], reverse=True)
            if len(truth_jets) >= 2:
                tj0 = truth_jets[0]
                tj1 = truth_jets[1]
                v0 = ROOT.TLorentzVector()
                v1 = ROOT.TLorentzVector()
                v0.SetPtEtaPhiM(tj0["pt"], tj0["eta"], tj0["phi"], tj0["m"])
                v1.SetPtEtaPhiM(tj1["pt"], tj1["eta"], tj1["phi"], tj1["m"])
                mjj_truthjet = (v0 + v1).M()  # GeV
                h_mjj_truthjet.Fill(mjj_truthjet)

        # ======================= LARGE-R JETS =========================
        largeR_jets = []
        J10_pt_evt  = J10_pt[evt]
        J10_eta_evt = J10_eta[evt]
        J10_m_evt   = J10_m[evt]

        if len(J10_pt_evt) > 0:
            pt_J  = ak.to_numpy(J10_pt_evt) / 1000.0
            eta_J = ak.to_numpy(J10_eta_evt)
            m_J   = ak.to_numpy(J10_m_evt) / 1000.0
            for pt, eta, mass in zip(pt_J, eta_J, m_J):
                    if (pt > jet_pt_min and abs(eta) < jet_eta_max):
                        largeR_jets.append({
                            "pt": float(pt),
                            "eta": float(eta),
                            "m": float(mass),
                        })

        h_largeR_jet_mult.Fill(len(largeR_jets))

        # ======================= TRUTH SELECTION ======================
        truth_leptons = []
        truth_quarks  = []

        truth_pdg_evt = truth_pdg[evt]
        truth_px_evt  = truth_px[evt]
        truth_py_evt  = truth_py[evt]
        truth_pz_evt = truth_pz[evt]
        truth_e_evt  = truth_e[evt]

        if len(truth_pdg_evt) > 0:
            pdg_evt = ak.to_numpy(truth_pdg_evt)
            px_evt  = ak.to_numpy(truth_px_evt)
            py_evt  = ak.to_numpy(truth_py_evt)
            pz_evt  = ak.to_numpy(truth_pz_evt)
            e_evt   = ak.to_numpy(truth_e_evt)

            # leptons (e, mu)
            mask_lep = np.isin(np.abs(pdg_evt), [11, 13])
            if np.any(mask_lep):
                pt_lep = calculate_pt(px_evt[mask_lep], py_evt[mask_lep])
                for pt in pt_lep:
                    truth_leptons.append({"pt": float(pt)})

            # quarks
            mask_q = np.isin(np.abs(pdg_evt), [1, 2, 3, 4, 5])
            if np.any(mask_q):
                px_q = px_evt[mask_q]
                py_q = py_evt[mask_q]
                pz_q = pz_evt[mask_q]
                e_q  = e_evt[mask_q]
                pt_q = calculate_pt(px_q, py_q)
                eta_q, phi_q = calculate_eta_phi(px_q, py_q, pz_q)
                order_q = np.argsort(-pt_q)
                if order_q.size >= 2:
                    i1, i2 = order_q[0], order_q[1]
                    px1, py1, pz1, e1 = px_q[i1], py_q[i1], pz_q[i1], e_q[i1]
                    px2, py2, pz2, e2 = px_q[i2], py_q[i2], pz_q[i2], e_q[i2]
                    px_tot = px1 + px2
                    py_tot = py1 + py2
                    pz_tot = pz1 + pz2
                    e_tot  = e1  + e2

                    m2 = e_tot**2 - px_tot**2 - py_tot**2 - pz_tot**2
                    m2 = np.maximum(0.0, m2)
                    m_qql = math.sqrt(m2) / 1000.0  # GeV
                    h_mjj_truth.Fill(m_qql)

                for pt, eta, phi in zip(pt_q, eta_q, phi_q):
                    truth_quarks.append({
                        "pt":  float(pt),
                        "eta": float(eta),
                        "phi": float(phi),
                    })

        truth_leptons = sort_by_pt(truth_leptons)
        n_Tlep_after = len(truth_leptons)
        Tlep_counts[n_Tlep_after] += 1

        truth_quarks  = sort_by_pt(truth_quarks)
        n_Tquark_after = len(truth_quarks)
        Tquark_counts[n_Tquark_after] += 1

        # ======================= HNL TRUTH ============================
        HNL_PDG = 9900012
        bsm_pdg_evt = bsm_pdg[evt]
        if len(bsm_pdg_evt) > 0:
            pdg_bsm_evt = ak.to_numpy(bsm_pdg_evt)
            idxs = np.where(pdg_bsm_evt == HNL_PDG)[0]
            if idxs.size > 0:
                idx = idxs[0]
                px_hnl = bsm_px[evt][idx]
                py_hnl = bsm_py[evt][idx]
                pz_hnl = bsm_pz[evt][idx]
                e_hnl  = bsm_e[evt][idx]

                pt_hnl = calculate_pt(px_hnl, py_hnl)
                m_hnl  = calculate_mass(px_hnl, py_hnl, pz_hnl, e_hnl)
                gamma_hnl = (e_hnl / 1000.0) / m_hnl if m_hnl > 0 else 0.0

                # h2_hnlgamma_vs_njet.Fill(gamma_hnl, n_jet_after)
                h2_hnlgamma_vs_nlep.Fill(gamma_hnl, n_lep_after)

        # ======================= FILL LEP / JET HISTOS =================
        # leptons (reco)
        for i, l in enumerate(leptons_sel):
            if i < 2:
                h_lep_lead_reco.Fill(l["pt"])
            else:
                h_lep_other.Fill(l["pt"])

        # leptons (truth)
        for i, l in enumerate(truth_leptons):
            if i < 2:
                h_lep_lead_truth.Fill(l["pt"])

        # jets (reco)
        for i, j in enumerate(jets_sel):
            if i < 2:
                h_jet_lead_reco.Fill(j["pt"])
            else:
                h_jet_other.Fill(j["pt"])

        # leading-jet shape info (if at least one jet)
        if n_jet_after > 0:
            j0 = jets_sel[0]
            h_jet_lead_mass.Fill(j0["m"])
            h_jet_lead_width.Fill(j0["Width"])
            h_jet_lead_ntrk.Fill(j0["NumTrkPt500"])

        # jets (truth from quarks)
        for i, j in enumerate(truth_quarks):
            if i < 2:
                h_jet_lead_truth.Fill(j["pt"])

        # ---------------- AUX: fill rank-based pT histos ----------------
        # leptons: reco (0:leading,1:subleading,2:sub-subleading)
        for i in range(min(3, len(leptons_sel))):
            h_lep_reco_rank[i].Fill(leptons_sel[i]["pt"])

        # leptons: truth (same convention)
        for i in range(min(3, len(truth_leptons))):
            h_lep_truth_rank[i].Fill(truth_leptons[i]["pt"])

        # jets: reco
        for i in range(min(3, len(jets_sel))):
            h_jet_reco_rank[i].Fill(jets_sel[i]["pt"])

        # jets: truth (from quarks)
        for i in range(min(3, len(truth_quarks))):
            h_jet_truth_rank[i].Fill(truth_quarks[i]["pt"])   

        # --- ΔR between leading and subleading jets: reco ---
        if len(jets_sel) >= 2:
            j0 = jets_sel[0]
            j1 = jets_sel[1]
            dR_reco = deltaR(j0["eta"], j0["phi"], j1["eta"], j1["phi"])
            h_dR_jet_reco.Fill(dR_reco)

        # --- ΔR between leading and subleading jets: truth AK4 ---
        if len(truth_jets) >= 2:
            # tj0 = truth_jets[0]
            # tj1 = truth_jets[1]
            tj0 = truth_quarks[0]
            tj1 = truth_quarks[1]
            dR_truth = deltaR(tj0["eta"], tj0["phi"], tj1["eta"], tj1["phi"])
            h_dR_jet_truth.Fill(dR_truth)  

print("\n=== Lepton Multiplicity (Before Iso Cuts) ===")
for k in sorted(lep_before_counts):
    print(f"N={k}: {lep_before_counts[k]}")

print("\n=== Lepton Multiplicity (After Iso Cuts) ===")
for k in sorted(lep_after_counts):
    print(f"N={k}: {lep_after_counts[k]}")

print("\n=== Jet Multiplicity (Before Cuts) ===")
for k in sorted(jet_before_counts):
    print(f"N={k}: {jet_before_counts[k]}")

print("\n=== Jet Multiplicity (After Cuts) ===")
for k in sorted(jet_after_counts):
    print(f"N={k}: {jet_after_counts[k]}")

print("\n=== False Truth Lepton Multiplicity Signatures ===")
for k in sorted(Tlep_counts):
    if k == 0 or k == 1:
        print(f"N={k}: {Tlep_counts[k]}") 

print("\n=== False Truth Quark Multiplicity Signatures ===")
for k in sorted(Tquark_counts):
    if k == 0 or k == 1:
        print(f"N={k}: {Tquark_counts[k]}")    

# ---------------------------------------------------------------------
# SAVE HISTOGRAMS
# ---------------------------------------------------------------------
ROOT.gStyle.SetOptStat(0)

# Multiplicities
c1 = ROOT.TCanvas("c1","",610,270)
c1.Divide(2,2)
navy = ROOT.TColor.GetColor("#000080")

for ipad in range(1, 5):
    pad = c1.cd(ipad)
    if ipad == 1 or  ipad == 2:
        pad.SetTopMargin(0.13)
    else:
        pad.SetTopMargin(0.05)   # check this!
    pad.SetBottomMargin(0.085)
    pad.SetLeftMargin(0.088)
    pad.SetRightMargin(0.05)
    pad.Update()

for h in [h_lep_before, h_lep_after, h_jet_before, h_jet_after]:
    h.SetLineColor(navy)
    h.SetFillColor(navy)
    h.SetFillStyle(1001)
    h.GetXaxis().SetTitleSize(0.07)
    h.GetXaxis().SetLabelSize(0.1)
    h.GetYaxis().SetTitleSize(0.07)
    h.GetYaxis().SetLabelSize(0.1)
    h.GetYaxis().SetNdivisions(505)

c1.cd(1); 
# h_lep_before.GetYaxis().SetNdivisions(505)
h_lep_before.Draw()
mini1 = ROOT.TPaveText(0.12, 0.65, 0.4, 0.85, "NDC")  # x1,y1,x2,y2
mini1.AddText(r"Dilepton+Rad.")
mini1.SetFillColor(0)
mini1.SetFillStyle(0)
mini1.SetBorderSize(0)
mini1.SetTextFont(32)
mini1.SetTextSize(0.13)
mini1.Draw()

c1.cd(2); 
# h_lep_after.GetYaxis().SetNdivisions(505)
h_lep_after.Draw()
mini2 = ROOT.TPaveText(0.6, 0.65, 0.8, 0.85, "NDC")  # x1,y1,x2,y2
mini2.AddText(r"post- Kin/Iso cuts")
mini2.SetFillColor(0)
mini2.SetFillStyle(0)
mini2.SetBorderSize(0)
mini2.SetTextFont(32)
mini2.SetTextSize(0.13)
mini2.Draw()

c1.cd(3); 
# h_jet_before.GetYaxis().SetNdivisions(505)
h_jet_before.Draw()
mini3 = ROOT.TPaveText(0.06, 0.65, 0.4, 0.85, "NDC")  # x1,y1,x2,y2
mini3.AddText(r"Dijet+Rad.")
mini3.SetFillColor(0)
mini3.SetFillStyle(0)
mini3.SetBorderSize(0)
mini3.SetTextFont(32)
mini3.SetTextSize(0.13)
mini3.Draw()

c1.cd(4); 
# h_jet_after.GetYaxis().SetNdivisions(505)
h_jet_after.Draw()
mini4 = ROOT.TPaveText(0.6, 0.73, 0.8, 0.93, "NDC")  # x1,y1,x2,y2  //65, 85
mini4.AddText(r"post- Kin/JVT cuts")
mini4.SetFillColor(0)
mini4.SetFillStyle(0)
mini4.SetBorderSize(0)
mini4.SetTextFont(32)
mini4.SetTextSize(0.13)
mini4.Draw()

c1.cd()  # go back to canvas, not a sub-pad
title = ROOT.TPaveText(0.15, 0.95, 0.85, 0.99, "NDC")  # (x1,y1,x2,y2) in NDC
title.AddText("Lepton/Jet Multiplicities for (mN=500GeV, mWR=5TeV)")
title.SetFillColor(0)   # no fill
title.SetFillStyle(0)   # transparent
title.SetBorderSize(0)  # no border
title.SetTextFont(32)   # or 32, etc.
title.SetTextSize(0.06)
title.Draw()
c1.Update()
c1.SaveAs(f"{output_dir}/multiplicities.png")
# c1.SaveAs(f"{output_dir}/30_multiplicities.png")

# Leading leptons reco vs truth
c2 = ROOT.TCanvas("c2","",800,600)
h_lep_lead_reco.SetLineColor(ROOT.kRed)
h_lep_lead_truth.SetLineColor(ROOT.kBlue)
h_lep_lead_reco.Draw("HIST")
h_lep_lead_truth.Draw("HIST SAME")
leg2 = ROOT.TLegend(0.7,0.7,0.9,0.9)
leg2.AddEntry(h_lep_lead_reco, "Reco leading leptons", "l")
leg2.AddEntry(h_lep_lead_truth, "Truth leading leptons", "l")
leg2.Draw()
c2.SaveAs(f"{output_dir}/leading_leptons.png")
# c2.SaveAs(f"{output_dir}/30_leading_leptons.png")

# Leading jets reco vs truth
c3 = ROOT.TCanvas("c3","",800,600)
h_jet_lead_reco.SetLineColor(ROOT.kRed)
h_jet_lead_truth.SetLineColor(ROOT.kBlue)
h_jet_lead_reco.Draw("HIST")
h_jet_lead_truth.Draw("HIST SAME")
leg3 = ROOT.TLegend(0.7,0.7,0.9,0.9)
leg3.AddEntry(h_jet_lead_reco, "Reco leading jets", "l")
leg3.AddEntry(h_jet_lead_truth, "Truth leading jets", "l")
leg3.Draw()
# c3.SaveAs(f"{output_dir}/30_leading_jets.png")
c3.SaveAs(f"{output_dir}/leading_jets.png")

c_mjj = ROOT.TCanvas("c_mjj", "", 800, 600)
ROOT.gPad.SetLogy()
h_mjj_truth.SetLineColor(ROOT.kBlue)
h_mjj_truth.SetLineWidth(3)
h_mjj_truthjet.SetLineColor(ROOT.kGreen+2)
h_mjj_truthjet.SetLineWidth(3)
h_mjj_reco.SetLineColor(ROOT.kRed)
h_mjj_reco.SetLineWidth(3)
h_mjj_reco.SetLineStyle(3)
h_mjj_reco.SetFillColor(ROOT.kRed)
h_mjj_reco.SetFillStyle(3004) 
h_mjj_reco2.SetLineColor(ROOT.kRed)
h_mjj_reco2.SetLineWidth(3)
h_mjj_reco3.SetLineColor(ROOT.kRed)
h_mjj_reco3.SetLineWidth(3)
h_mjj_reco3.SetLineStyle(2)

# common y-range
maxY = max(
    h_mjj_truth.GetMaximum(),
    h_mjj_truthjet.GetMaximum(),
    h_mjj_reco.GetMaximum()
)
h_mjj_truth.SetMaximum(1.2 * maxY)
h_mjj_truth.SetMinimum(2)
h_mjj_truth.Draw("HIST")
h_mjj_truthjet.Draw("HIST SAME")
h_mjj_reco.Draw("HIST SAME")
h_mjj_reco2.Draw("HIST SAME")
h_mjj_reco3.Draw("HIST SAME")
ROOT.gPad.Update()
pt = ROOT.gPad.GetPrimitive("title")
if pt:
    pt.SetTextSize(0.045)   # relative size in pad
    pt.SetTextFont(32)  
leg_mjj = ROOT.TLegend(0.25, 0.7, 0.9, 0.90) 
leg_mjj.SetBorderSize(0)
leg_mjj.SetNColumns(3)
leg_mjj.SetFillStyle(0)
leg_mjj.SetTextFont(32)
leg_mjj.SetTextSize(0.04)
leg_mjj.AddEntry(h_mjj_reco2, "reco dijet", "l")
leg_mjj.AddEntry(h_mjj_truthjet, "truth dijet", "l")
leg_mjj.AddEntry(h_mjj_truth, "gen. truth diquark", "l")
leg_mjj.AddEntry(h_mjj_reco3, "jet obj. cuts*", "l")
leg_mjj.AddEntry(0, "", "")
leg_mjj.AddEntry(0, "", "")
leg_mjj.AddEntry(h_mjj_reco, f"#DeltaR#leq{dr_jj_cut}", "f")
leg_mjj.Draw()
c_mjj.SaveAs(f"{output_dir}/mjj_truth_quark_truthjet_reco.png")


# Other leptons/jets
c4 = ROOT.TCanvas("c4","",1200,600)
c4.Divide(2,1)
c4.cd(1)
h_lep_other.SetMarkerStyle(20)
h_lep_other.SetMarkerColor(ROOT.kGreen+2)
h_lep_other.SetLineColor(ROOT.kGreen+2)
h_lep_other.Draw("HIST")  
c4.cd(2); h_jet_other.Draw("HIST")
# c4.SaveAs(f"{output_dir}/30_others.png")
c4.SaveAs(f"{output_dir}/others.png")

# Overlay leptons, log-y
c5 = ROOT.TCanvas("c5","",800,600)
ROOT.gPad.SetLogy()
ROOT.gPad.SetLineWidth(3)
h_lep_lead_reco.SetLineColor(ROOT.kRed)
h_lep_lead_truth.SetLineColor(ROOT.kBlue)
h_lep_other.SetLineColor(ROOT.kGreen+2)
h_lep_lead_reco.SetLineWidth(2)
h_lep_lead_truth.SetLineWidth(2)
h_lep_other.SetLineWidth(2)

h_lep_lead_reco.Draw("HIST")
h_lep_lead_truth.Draw("HIST SAME")
h_lep_other.Draw("HIST SAME")
ROOT.gPad.Update()
pt = ROOT.gPad.GetPrimitive("title")
if pt:
    pt.SetTextSize(0.045)   # relative size in pad
    pt.SetTextFont(32)  
leg5 = ROOT.TLegend(0.25, 0.75, 0.85, 0.90) 
leg5.SetTextFont(32)     
leg5.SetTextSize(0.04)   
leg5.SetBorderSize(0)    
leg5.SetFillStyle(0)     
leg5.SetNColumns(3)
leg5.AddEntry(h_lep_lead_reco, "reconstructed", "l")
leg5.AddEntry(h_lep_lead_truth, "truth", "l")
leg5.AddEntry(h_lep_other, "ISR/FSR", "l")
leg5.Draw()
# c5.SaveAs(f"{output_dir}/30_overlay_log_lep.png")
c5.SaveAs(f"{output_dir}/overlay_log_lep.png")

# Overlay jets, log-y
c6 = ROOT.TCanvas("c6","",800,600)
ROOT.gPad.SetLogy()
ROOT.gPad.SetLineWidth(3)
h_jet_lead_reco.SetLineColor(ROOT.kRed)
h_jet_lead_truth.SetLineColor(ROOT.kBlue)
h_jet_other.SetLineColor(ROOT.kGreen+2)
h_jet_lead_reco.SetLineWidth(2)
h_jet_lead_truth.SetLineWidth(2)
h_jet_other.SetLineWidth(2)

h_jet_lead_reco2.SetLineColor(ROOT.kRed)
h_jet_lead_reco2.SetLineWidth(2)
h_jet_lead_reco2.SetLineStyle(3)
h_jet_lead_reco2.SetFillColor(ROOT.kRed)
h_jet_lead_reco2.SetFillStyle(3004) 

h_jet_lead_reco3.SetLineColor(ROOT.kRed)
h_jet_lead_reco3.SetLineWidth(3)
h_jet_lead_reco3.SetLineStyle(3)

h_jet_lead_reco4.SetLineColor(ROOT.kGray+2)
h_jet_lead_reco4.SetLineWidth(3)
h_jet_lead_reco4.SetLineStyle(2)

h_jet_lead_reco.Draw("HIST")
h_jet_lead_reco2.Draw("HIST SAME")
h_jet_lead_reco3.Draw("HIST SAME")
h_jet_lead_reco4.Draw("HIST SAME")
h_jet_lead_truth.Draw("HIST SAME")
# h_jet_other.Draw("HIST SAME")
ROOT.gPad.Update()
pt = ROOT.gPad.GetPrimitive("title")
if pt:
    pt.SetTextSize(0.045)   # relative size in pad
    pt.SetTextFont(32)  
leg6 = ROOT.TLegend(0.2, 0.75, 0.9, 0.90) 
leg6.SetTextFont(32)     
leg6.SetTextSize(0.04)   
leg6.SetBorderSize(0)    
leg6.SetFillStyle(0)     
leg6.SetNColumns(3)     
leg6.AddEntry(h_jet_lead_truth, "gen. truth (quarks)", "l") 
# leg6.AddEntry(h_jet_other, "ISR/FSR", "l")
leg6.AddEntry(h_jet_lead_reco, "reconstructed", "l")
leg6.AddEntry(h_jet_lead_reco3, "jet obj. cuts", "l") 
leg6.AddEntry(0, "", "")
leg6.AddEntry(h_jet_lead_reco4, f"#DeltaR#leq{dr_jj_cut} (alone)", "l") 
leg6.AddEntry(h_jet_lead_reco2, f"#DeltaR#leq{dr_jj_cut}", "f")
leg6.Draw()
# c6.SaveAs(f"{output_dir}/30_overlay_log_jet.png")
c6.SaveAs(f"{output_dir}/overlay_log_jet.png")

# Lepton isolation distributions
c7 = ROOT.TCanvas("c7", "", 1200, 600)
c7.Divide(2,1)
c7.cd(1); h_lep_rel_iso_ptvar.Draw("HIST")
c7.cd(2); h_lep_rel_iso_topo.Draw("HIST")
# c7.SaveAs(f"{output_dir}/30_lepton_isolation.png")
c7.SaveAs(f"{output_dir}/lepton_isolation.png")

# Jet shape distributions (leading jet)
# c8 = ROOT.TCanvas("c8", "", 1200, 600)
# c8.Divide(3,1)
# c8.cd(1); h_jet_lead_mass.Draw("HIST")
# c8.cd(2); h_jet_lead_width.Draw("HIST")
# c8.cd(3); h_jet_lead_ntrk.Draw("HIST")
# # c8.SaveAs(f"{output_dir}/30_leading_jet_shapes.png")
c8 = ROOT.TCanvas("c8", "", 720, 320)
c8.Divide(2,1)
c8.cd(1)
ROOT.gPad.SetTopMargin(0.08)
ROOT.gPad.SetBottomMargin(0.12)
ROOT.gPad.SetLeftMargin(0.14)
miniLR = ROOT.TPaveText(0.55, 0.65, 0.9, 0.85, "NDC")
miniLR.AddText(r"post-Kin/JVT cuts")
miniLR.Draw()
h_jet_lead_width.SetTitle("")
h_jet_lead_width.SetLineWidth(3)
h_jet_lead_width.GetYaxis().SetTitle("Events")
h_jet_lead_width.GetYaxis().SetTitleSize(0.05)
h_jet_lead_width.GetYaxis().SetLabelSize(0.04)
h_jet_lead_width.GetXaxis().SetLabelSize(0.05)
h_jet_lead_width.GetXaxis().SetTitleSize(0.05)
h_jet_lead_width.GetYaxis().SetNdivisions(505)
h_jet_lead_width.Draw("HIST")
c8.cd(2)
ROOT.gPad.SetTopMargin(0.08)
ROOT.gPad.SetBottomMargin(0.12)
ROOT.gPad.SetLeftMargin(0.14)
miniLR = ROOT.TPaveText(0.55, 0.65, 0.9, 0.85, "NDC")
miniLR.AddText(r"post-Kin/JVT cuts")
miniLR.Draw()
h_jet_lead_ntrk.SetTitle("")
h_jet_lead_ntrk.SetLineWidth(3)
h_jet_lead_ntrk.GetYaxis().SetTitle("Events")
h_jet_lead_ntrk.GetYaxis().SetTitleSize(0.05)
h_jet_lead_ntrk.GetYaxis().SetLabelSize(0.04)
h_jet_lead_ntrk.GetXaxis().SetLabelSize(0.05)
h_jet_lead_ntrk.GetXaxis().SetTitleSize(0.05)
h_jet_lead_ntrk.GetYaxis().SetNdivisions(505)
h_jet_lead_ntrk.Draw("HIST")
c8.cd()
title = ROOT.TPaveText(0.25, 0.95, 0.85, 0.99, "NDC")
title.AddText(r"Low-R Jet Properties for (mN=50GeV, mWR=10TeV)")
title.SetFillColor(0)
title.SetFillStyle(0)
title.SetBorderSize(0)
title.SetTextFont(32)
title.SetTextSize(0.08)
title.Draw()
c8.Update()
c8.SaveAs(f"{output_dir}/leading_jet_shapes.png")

# Large-R jet multiplicity
# c10 = ROOT.TCanvas("c10", "", 800, 600)
# h_largeR_jet_mult.Draw("HIST")
# # c10.SaveAs(f"{output_dir}/30_largeR_jet_mult.png")
# c10.SaveAs(f"{output_dir}/largeR_jet_mult.png")
c10 = ROOT.TCanvas("c10", "", 610, 320)  # SAME SIZE AS c1
c10.SetTopMargin(0.12)
c10.SetBottomMargin(0.1)
c10.SetLeftMargin(0.088)
c10.SetRightMargin(0.05)
h_largeR_jet_mult.SetLineColor(navy)
h_largeR_jet_mult.SetFillColor(navy)
h_largeR_jet_mult.SetFillStyle(1001)
h_largeR_jet_mult.GetXaxis().SetTitleSize(0.005)
h_largeR_jet_mult.GetXaxis().SetLabelSize(0.07)
h_largeR_jet_mult.GetYaxis().SetTitleSize(0.005)
h_largeR_jet_mult.GetYaxis().SetLabelSize(0.07)
h_largeR_jet_mult.GetXaxis().SetNdivisions(505)
h_largeR_jet_mult.GetYaxis().SetNdivisions(505)
h_largeR_jet_mult.Draw("HIST")
miniLR = ROOT.TPaveText(0.55, 0.65, 0.9, 0.85, "NDC")
miniLR.AddText(r"post-Kin/JVT cuts")
miniLR.SetFillColor(0)
miniLR.SetFillStyle(0)
miniLR.SetBorderSize(0)
miniLR.SetTextFont(32)
miniLR.SetTextSize(0.1)
miniLR.Draw()
title = ROOT.TPaveText(0.25, 0.95, 0.85, 0.99, "NDC")
title.AddText(r"Large-R Jet Multiplicity for (mN=50GeV, mWR=10TeV)")
title.SetFillColor(0)
title.SetFillStyle(0)
title.SetBorderSize(0)
title.SetTextFont(32)
title.SetTextSize(0.08)
title.Draw()
c10.Update()
c10.SaveAs(f"{output_dir}/largeR_jet_mult.png")
# HNL gamma vs jet multiplicity
c13 = ROOT.TCanvas("c13", "", 800, 310)
c13.cd().SetTopMargin(0.08)
# c13.SetFillColor(ROOT.kBlack)
ROOT.gStyle.SetPalette(ROOT.kBird)   #kBird , kBlueRedYellow , kGreyScale, kVirdis
# ROOT.gStyle.SetNumberContours(500)
# c13.SetLogz()  
h2_hnlgamma_vs_nlep.Draw("COLZ")
# h2_hnlgamma_vs_nlep.SetMaximum(1000)
h2_hnlgamma_vs_nlep.SetMinimum(1)
h2_hnlgamma_vs_nlep.GetXaxis().SetTitleSize(0.05)
h2_hnlgamma_vs_nlep.GetXaxis().SetLabelSize(0.05)
h2_hnlgamma_vs_nlep.GetYaxis().SetTitleSize(0.1)
h2_hnlgamma_vs_nlep.GetYaxis().SetLabelSize(0.1)
h2_hnlgamma_vs_nlep.GetXaxis().SetTitleOffset(0.8) 
h2_hnlgamma_vs_nlep.GetYaxis().SetBinLabel(1, "1 Iso lep")
h2_hnlgamma_vs_nlep.GetYaxis().SetBinLabel(2, "2 Iso lep")
ROOT.gPad.Update()
pt = ROOT.gPad.GetPrimitive("title")
if pt:
    pt.SetTextSize(0.08)   # relative size in pad
    pt.SetTextFont(32)  
# c13.SaveAs(f"{output_dir}/30_hnl_gamma_vs_nlep.png")
c13.SaveAs(f"{output_dir}/hnl_gamma_vs_nlep2.png")

# Auxiliary: overlay truth vs reco leptons by rank (leading, sub, subsub)
c_lep_rank = ROOT.TCanvas("c_lep_rank", "", 800, 600)
ROOT.gPad.SetLogy()
ROOT.gPad.SetLineWidth(3)

colors = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2]

# Style reco
for i, h in enumerate(h_lep_reco_rank):
    h.SetLineColor(colors[i])
    h.SetLineWidth(2)

# Style truth
for i, h in enumerate(h_lep_truth_rank):
    h.SetLineColor(colors[i])
    h.SetLineStyle(2)  # dashed for truth
    h.SetLineWidth(2)

# Normalize y-axis range
maxY = max(
    max(h.GetMaximum() for h in h_lep_reco_rank),
    max(h.GetMaximum() for h in h_lep_truth_rank),
)
h_lep_reco_rank[0].SetMaximum(1.2 * maxY)
h_lep_reco_rank[0].SetMinimum(1e-3)

# Draw reco as solid, then truth as dashed
h_lep_reco_rank[0].Draw("HIST")
h_lep_reco_rank[1].Draw("HIST SAME")
h_lep_reco_rank[2].Draw("HIST SAME")
h_lep_truth_rank[0].Draw("HIST SAME")
h_lep_truth_rank[1].Draw("HIST SAME")
h_lep_truth_rank[2].Draw("HIST SAME")

leg_lep_rank = ROOT.TLegend(0.45, 0.65, 0.88, 0.90)
leg_lep_rank.SetBorderSize(0)
leg_lep_rank.SetFillStyle(0)
leg_lep_rank.SetTextFont(32)
leg_lep_rank.SetTextSize(0.035)
leg_lep_rank.AddEntry(h_lep_reco_rank[0],  "Reco lep: leading",          "l")
leg_lep_rank.AddEntry(h_lep_reco_rank[1],  "Reco lep: subleading",       "l")
leg_lep_rank.AddEntry(h_lep_reco_rank[2],  "Reco lep: sub-subleading",   "l")
leg_lep_rank.AddEntry(h_lep_truth_rank[0], "Truth lep: leading (dashed)", "l")
leg_lep_rank.AddEntry(h_lep_truth_rank[1], "Truth lep: subleading",      "l")
leg_lep_rank.AddEntry(h_lep_truth_rank[2], "Truth lep: sub-subleading",  "l")
leg_lep_rank.Draw()

c_lep_rank.SaveAs(f"{output_dir}/aux_lepton_pt_rank_truth_vs_reco.png")


# Auxiliary: overlay truth vs reco jets by rank
c_jet_rank = ROOT.TCanvas("c_jet_rank", "", 800, 600)
ROOT.gPad.SetLogy()
ROOT.gPad.SetLineWidth(3)

# reuse same colors
# Style reco
for i, h in enumerate(h_jet_reco_rank):
    h.SetLineColor(colors[i])
    h.SetLineWidth(2)

# Style truth
for i, h in enumerate(h_jet_truth_rank):
    h.SetLineColor(colors[i])
    h.SetLineStyle(2)  # dashed for truth
    h.SetLineWidth(2)

maxY = max(
    max(h.GetMaximum() for h in h_jet_reco_rank),
    max(h.GetMaximum() for h in h_jet_truth_rank),
)
h_jet_reco_rank[0].SetMaximum(1.2 * maxY)
h_jet_reco_rank[0].SetMinimum(1e-3)

h_jet_reco_rank[0].Draw("HIST")
h_jet_reco_rank[1].Draw("HIST SAME")
h_jet_reco_rank[2].Draw("HIST SAME")
h_jet_truth_rank[0].Draw("HIST SAME")
h_jet_truth_rank[1].Draw("HIST SAME")
h_jet_truth_rank[2].Draw("HIST SAME")

leg_jet_rank = ROOT.TLegend(0.45, 0.65, 0.88, 0.90)
leg_jet_rank.SetBorderSize(0)
leg_jet_rank.SetFillStyle(0)
leg_jet_rank.SetTextFont(32)
leg_jet_rank.SetTextSize(0.035)
leg_jet_rank.AddEntry(h_jet_reco_rank[0],  "Reco jet: leading",          "l")
leg_jet_rank.AddEntry(h_jet_reco_rank[1],  "Reco jet: subleading",       "l")
leg_jet_rank.AddEntry(h_jet_reco_rank[2],  "Reco jet: sub-subleading",   "l")
leg_jet_rank.AddEntry(h_jet_truth_rank[0], "Truth jet: leading (dashed)", "l")
leg_jet_rank.AddEntry(h_jet_truth_rank[1], "Truth jet: subleading",      "l")
leg_jet_rank.AddEntry(h_jet_truth_rank[2], "Truth jet: sub-subleading",  "l")
leg_jet_rank.Draw()

c_jet_rank.SaveAs(f"{output_dir}/aux_jet_pt_rank_truth_vs_reco.png")



# ΔR(j1,j2) debug: reco vs truth
c_dR = ROOT.TCanvas("c_dR", "", 800, 600)
ROOT.gPad.SetLineWidth(3)
h_dR_jet_reco.SetLineColor(ROOT.kRed)
h_dR_jet_truth.SetLineColor(ROOT.kBlue)
h_dR_jet_reco.SetLineWidth(2)
h_dR_jet_truth.SetLineWidth(2)
# Set common y-range
maxY_dR = max(h_dR_jet_reco.GetMaximum(), h_dR_jet_truth.GetMaximum())
h_dR_jet_reco.SetMaximum(1.2 * maxY_dR)
h_dR_jet_reco.SetMinimum(0)
h_dR_jet_reco.Draw("HIST")
h_dR_jet_truth.Draw("HIST SAME")
leg_dR = ROOT.TLegend(0.55, 0.7, 0.88, 0.88)
leg_dR.SetBorderSize(0)
leg_dR.SetFillStyle(0)
leg_dR.SetTextFont(32)
leg_dR.SetTextSize(0.04)
leg_dR.AddEntry(h_dR_jet_reco,  "Reco jets: #DeltaR(j_{1}, j_{2})",  "l")
leg_dR.AddEntry(h_dR_jet_truth, "Truth AK4 jets: #DeltaR(j_{1}, j_{2})", "l")
leg_dR.Draw()
c_dR.SaveAs(f"{output_dir}/aux_deltaR_lead_sublead_jets_truth_vs_reco.png")

print(f"\nAnalysis Successfully completed! Results saved to {output_dir}/")
