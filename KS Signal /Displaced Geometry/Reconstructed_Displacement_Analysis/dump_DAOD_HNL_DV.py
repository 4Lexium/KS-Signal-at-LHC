#!/usr/bin/env python3
import uproot
import awkward as ak
import numpy as np

# ============================
# User configuration
# ============================
signal_file = "/home/alexanum/WORKSPACE/RUN_CONTROLS/ATLAS_jobdir/ROOT_files/DAOD_LLP1.mc23_13p6TeV.150001.MGPy8EG_LRSM_WR_HNL_5000_50_2L2J_MadSpin.DAOD_LLP1.pool.root"
tree_name   = "CollectionTree"
HNL_PDG     = 9900012   # change if your HNL uses different PDG id

# ============================
# Open file and load branches
# ============================
branches = [
    # Truth BSM particles and vertices <ref: index=3732371 firstWord=1 lastWord=40/>
    "TruthBSMWithDecayParticlesAuxDyn.pdgId",
    "TruthBSMWithDecayParticlesAuxDyn.status",
    "TruthBSMWithDecayParticlesAuxDyn.px",
    "TruthBSMWithDecayParticlesAuxDyn.py",
    "TruthBSMWithDecayParticlesAuxDyn.pz",
    "TruthBSMWithDecayParticlesAuxDyn.e",

    "TruthBSMWithDecayVerticesAuxDyn.x",
    "TruthBSMWithDecayVerticesAuxDyn.y",
    "TruthBSMWithDecayVerticesAuxDyn.z",
    "TruthBSMWithDecayVerticesAuxDyn.outgoingParticleLinks",

    # Inclusive SV collections (reco) and their dCloseVrt / positions
    # 1) "default" inclusive <ref: index=3732513 firstWord=1 lastWord=40/>
    "VrtSecInclusive_SecondaryVerticesAuxDyn.vtx_px",
    "VrtSecInclusive_SecondaryVerticesAuxDyn.vtx_py",
    "VrtSecInclusive_SecondaryVerticesAuxDyn.vtx_pz",
    "VrtSecInclusive_SecondaryVerticesAuxDyn.dCloseVrt",

    # 2) Lepton‑based <ref: index=3732533 firstWord=40 lastWord=80/>
    "VrtSecInclusive_SecondaryVertices_LeptonsAuxDyn.vtx_px",
    "VrtSecInclusive_SecondaryVertices_LeptonsAuxDyn.vtx_py",
    "VrtSecInclusive_SecondaryVertices_LeptonsAuxDyn.vtx_pz",
    "VrtSecInclusive_SecondaryVertices_LeptonsAuxDyn.dCloseVrt",

    # 3) LepTrack‑based <ref: index=3732512 firstWord=1 lastWord=40/>
    "VrtSecInclusive_SecondaryVertices_LepTrackAuxDyn.vtx_px",
    "VrtSecInclusive_SecondaryVertices_LepTrackAuxDyn.vtx_py",
    "VrtSecInclusive_SecondaryVertices_LepTrackAuxDyn.vtx_pz",
    "VrtSecInclusive_SecondaryVertices_LepTrackAuxDyn.dCloseVrt",

    # 4) shortLifetime <ref: index=3732512 firstWord=40 lastWord=80/>, <ref: index=3732400 firstWord=1 lastWord=40/>
    "VrtSecInclusive_SecondaryVertices_shortLifetimeAuxDyn.vtx_px",
    "VrtSecInclusive_SecondaryVertices_shortLifetimeAuxDyn.vtx_py",
    "VrtSecInclusive_SecondaryVertices_shortLifetimeAuxDyn.vtx_pz",
    "VrtSecInclusive_SecondaryVertices_shortLifetimeAuxDyn.dCloseVrt",

    # 5) disappearing <ref: index=3727075 firstWord=1 lastWord=40/>
    "VrtSecInclusive_SecondaryVertices_disappearingAuxDyn.vtx_px",
    "VrtSecInclusive_SecondaryVertices_disappearingAuxDyn.vtx_py",
    "VrtSecInclusive_SecondaryVertices_disappearingAuxDyn.vtx_pz",
    "VrtSecInclusive_SecondaryVertices_disappearingAuxDyn.dCloseVrt",

    # 6) BoostedElectrons <ref: index=3727082 firstWord=1 lastWord=40/>
    "VrtSecInclusive_SecondaryVertices_BoostedElectronsAuxDyn.vtx_px",
    "VrtSecInclusive_SecondaryVertices_BoostedElectronsAuxDyn.vtx_py",
    "VrtSecInclusive_SecondaryVertices_BoostedElectronsAuxDyn.vtx_pz",
    "VrtSecInclusive_SecondaryVertices_BoostedElectronsAuxDyn.dCloseVrt",

    # 7) BoostedMuons <ref: index=3727120 firstWord=1 lastWord=40/>
    "VrtSecInclusive_SecondaryVertices_BoostedMuonsAuxDyn.vtx_px",
    "VrtSecInclusive_SecondaryVertices_BoostedMuonsAuxDyn.vtx_py",
    "VrtSecInclusive_SecondaryVertices_BoostedMuonsAuxDyn.vtx_pz",
    "VrtSecInclusive_SecondaryVertices_BoostedMuonsAuxDyn.dCloseVrt",
]

print(f"Opening signal file: {signal_file}")
with uproot.open(signal_file) as f:
    tree = f[tree_name]
    arrays = tree.arrays(branches, library="ak")

# shorthand variables
bsmp_pdg  = arrays["TruthBSMWithDecayParticlesAuxDyn.pdgId"]
bsmp_stat = arrays["TruthBSMWithDecayParticlesAuxDyn.status"]
bsmp_px   = arrays["TruthBSMWithDecayParticlesAuxDyn.px"]
bsmp_py   = arrays["TruthBSMWithDecayParticlesAuxDyn.py"]
bsmp_pz   = arrays["TruthBSMWithDecayParticlesAuxDyn.pz"]
bsmp_e    = arrays["TruthBSMWithDecayParticlesAuxDyn.e"]

bsmv_x    = arrays["TruthBSMWithDecayVerticesAuxDyn.x"]
bsmv_y    = arrays["TruthBSMWithDecayVerticesAuxDyn.y"]
bsmv_z    = arrays["TruthBSMWithDecayVerticesAuxDyn.z"]
bsmv_out  = arrays["TruthBSMWithDecayVerticesAuxDyn.outgoingParticleLinks"]

# SV collections
sv_defs = {
    "SV_default": {
        "px": arrays["VrtSecInclusive_SecondaryVerticesAuxDyn.vtx_px"],
        "py": arrays["VrtSecInclusive_SecondaryVerticesAuxDyn.vtx_py"],
        "pz": arrays["VrtSecInclusive_SecondaryVerticesAuxDyn.vtx_pz"],
        "dR": arrays["VrtSecInclusive_SecondaryVerticesAuxDyn.dCloseVrt"],
    },
    "SV_Leptons": {
        "px": arrays["VrtSecInclusive_SecondaryVertices_LeptonsAuxDyn.vtx_px"],
        "py": arrays["VrtSecInclusive_SecondaryVertices_LeptonsAuxDyn.vtx_py"],
        "pz": arrays["VrtSecInclusive_SecondaryVertices_LeptonsAuxDyn.vtx_pz"],
        "dR": arrays["VrtSecInclusive_SecondaryVertices_LeptonsAuxDyn.dCloseVrt"],
    },
    "SV_LepTrack": {
        "px": arrays["VrtSecInclusive_SecondaryVertices_LepTrackAuxDyn.vtx_px"],
        "py": arrays["VrtSecInclusive_SecondaryVertices_LepTrackAuxDyn.vtx_py"],
        "pz": arrays["VrtSecInclusive_SecondaryVertices_LepTrackAuxDyn.vtx_pz"],
        "dR": arrays["VrtSecInclusive_SecondaryVertices_LepTrackAuxDyn.dCloseVrt"],
    },
    "SV_shortLifetime": {
        "px": arrays["VrtSecInclusive_SecondaryVertices_shortLifetimeAuxDyn.vtx_px"],
        "py": arrays["VrtSecInclusive_SecondaryVertices_shortLifetimeAuxDyn.vtx_py"],
        "pz": arrays["VrtSecInclusive_SecondaryVertices_shortLifetimeAuxDyn.vtx_pz"],
        "dR": arrays["VrtSecInclusive_SecondaryVertices_shortLifetimeAuxDyn.dCloseVrt"],
    },
    "SV_disappearing": {
        "px": arrays["VrtSecInclusive_SecondaryVertices_disappearingAuxDyn.vtx_px"],
        "py": arrays["VrtSecInclusive_SecondaryVertices_disappearingAuxDyn.vtx_py"],
        "pz": arrays["VrtSecInclusive_SecondaryVertices_disappearingAuxDyn.vtx_pz"],
        "dR": arrays["VrtSecInclusive_SecondaryVertices_disappearingAuxDyn.dCloseVrt"],
    },
    "SV_BoostedElectrons": {
        "px": arrays["VrtSecInclusive_SecondaryVertices_BoostedElectronsAuxDyn.vtx_px"],
        "py": arrays["VrtSecInclusive_SecondaryVertices_BoostedElectronsAuxDyn.vtx_py"],
        "pz": arrays["VrtSecInclusive_SecondaryVertices_BoostedElectronsAuxDyn.vtx_pz"],
        "dR": arrays["VrtSecInclusive_SecondaryVertices_BoostedElectronsAuxDyn.dCloseVrt"],
    },
    "SV_BoostedMuons": {
        "px": arrays["VrtSecInclusive_SecondaryVertices_BoostedMuonsAuxDyn.vtx_px"],
        "py": arrays["VrtSecInclusive_SecondaryVertices_BoostedMuonsAuxDyn.vtx_py"],
        "pz": arrays["VrtSecInclusive_SecondaryVertices_BoostedMuonsAuxDyn.vtx_pz"],
        "dR": arrays["VrtSecInclusive_SecondaryVertices_BoostedMuonsAuxDyn.dCloseVrt"],
    },
}

n_events = len(bsmp_pdg)
print(f"\nNumber of events in file: {n_events}")

# ============================
# Per‑event debugging for first N events
# ============================
N_DEBUG = min(10, n_events)

for ievt in range(N_DEBUG):
    print("\n" + "="*70)
    print(f"Event {ievt}")

    # ---- Truth BSM particles (HNLs etc.) ----
    pdgs_evt  = bsmp_pdg[ievt]
    stat_evt  = bsmp_stat[ievt]
    px_evt    = bsmp_px[ievt]
    py_evt    = bsmp_py[ievt]
    pz_evt    = bsmp_pz[ievt]
    e_evt     = bsmp_e[ievt]

    n_bsm = len(pdgs_evt)
    print(f"  #Truth BSM particles: {n_bsm}")

    # Find HNLs
    hnl_indices = [i for i, pdg in enumerate(pdgs_evt) if pdg == HNL_PDG]
    print(f"  HNL indices (pdgId == {HNL_PDG}): {hnl_indices}")

    # Print a small table of BSM particles for this event
    for i in range(n_bsm):
        print(f"    BSM[{i}]: pdgId={pdgs_evt[i]}, status={stat_evt[i]}, "
              f"px={px_evt[i]}, py={py_evt[i]}, pz={pz_evt[i]}, e={e_evt[i]}")

    # ---- Truth BSM vertices ----
    xv = bsmv_x[ievt]
    yv = bsmv_y[ievt]
    zv = bsmv_z[ievt]
    out_links = bsmv_out[ievt]

    n_vtx = len(xv)
    print(f"  #Truth BSM vertices: {n_vtx}")
    for iv in range(n_vtx):
        print(f"    Vtx[{iv}]: (x,y,z)=({xv[iv]}, {yv[iv]}, {zv[iv]})")

    # ---- SV collections ----
    for name, coll in sv_defs.items():
        dR_evt = coll["dR"][ievt]
        n_sv_evt = len(dR_evt)
        print(f"  {name}: N(SV) = {n_sv_evt}")
        if n_sv_evt > 0:
            # Print first few dCloseVrt values
            dR_np = np.asarray(ak.to_numpy(dR_evt), dtype=float)
            print(f"    dCloseVrt (first up to 5): {dR_np[:5]}")

# ============================
# Global summary: how many events have any SVs,
# and any non‑zero dCloseVrt per collection
# ============================
print("\n" + "="*70)
print("Global SV summary over all events:")

for name, coll in sv_defs.items():
    dR_all = coll["dR"]
    n_events_total = len(dR_all)

    # count events with at least one vertex
    n_evt_with_sv = sum(len(dR_all[ievt]) > 0 for ievt in range(n_events_total))

    # flatten all dCloseVrt values
    flat = ak.flatten(dR_all, axis=None)
    flat_np = np.asarray(ak.to_numpy(flat), dtype=float)
    n_total_sv = flat_np.size
    n_nonzero  = np.count_nonzero(flat_np)
    max_dR     = flat_np.max() if flat_np.size > 0 else None

    print(f"\nCollection: {name}")
    print(f"  Events with >=1 vertex: {n_evt_with_sv} / {n_events_total}")
    print(f"  Total #SV entries: {n_total_sv}")
    print(f"  #SV with non‑zero dCloseVrt: {n_nonzero}")
    print(f"  max(dCloseVrt) = {max_dR}")

print("\nDone.")
