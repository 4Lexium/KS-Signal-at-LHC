# IMPORT LIBRARIES
#------------------------------------------------------------------------------
import uproot
import ROOT
import math
import awkward as ak
import numpy as np
import os
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

print("User Test 026:")

# FILE HANDLING (LLP1 FILE)
#----------------------------------------------------------------------------
# file = uproot.open("/storage/alexanum/DAOD_LLP_WORKFLOW/ROOT_files/DAOD_LLP1.DAOD_LLP1.mc23_13p6TeV.MGPy8EG_LRSM_WR_HNL_5000_30_2L2J_MadSpin.DAOD_LLP1.pool.root")
file = uproot.open("/storage/alexanum/DAOD_LLP_WORKFLOW/ROOT_files/DAOD_LLP1.DAOD_LLP1.mc23_13p6TeV.MGPy8EG_LRSM_WR_HNL_5000_70_2L2J_MadSpin.DAOD_LLP1.pool.root")
# file = uproot.open("/storage/alexanum/DAOD_LLP_WORKFLOW/ROOT_files/DAOD_LLP1.DAOD_LLP1.mc23_13p6TeV.MGPy8EG_LRSM_WR_HNL_10000_50_2L2J_MadSpin.DAOD_LLP1.pool.root")
# file = uproot.open("/home/alexanum/WORKSPACE/RUN_CONTROLS/ATLAS_jobdir/ROOT_files/DAOD_LLP1.mc23_13p6TeV.150001.MGPy8EG_LRSM_WR_HNL_5000_50_2L2J_MadSpin.DAOD_LLP1.pool.root")
tree = file["CollectionTree"]
output_dir = "Figures/SIGvsBKG_DVtest4"
os.makedirs(output_dir, exist_ok=True)

# READ BSM- DECAY/VERTEX BRANCHES
#----------------------------------------------------------------------------
pdgId       = tree["TruthBSMWithDecayParticlesAuxDyn.pdgId"].array(library="ak")
decay_index = tree["TruthBSMWithDecayParticlesAuxDyn.decayVtxLink.m_persIndex"].array(library="ak")
prod_index  = tree["TruthBSMWithDecayParticlesAuxDyn.prodVtxLink.m_persIndex"].array(library="ak")

vx_all = tree["TruthBSMWithDecayVerticesAuxDyn.x"].array(library="ak")
vy_all = tree["TruthBSMWithDecayVerticesAuxDyn.y"].array(library="ak")
vz_all = tree["TruthBSMWithDecayVerticesAuxDyn.z"].array(library="ak")

# RECO LEPTONS 
#----------------------------------------------------------------------------
el_pt   = tree["ElectronsAuxDyn.pt"].array(library="ak")
el_eta  = tree["ElectronsAuxDyn.eta"].array(library="ak")
el_charge = tree["ElectronsAuxDyn.charge"].array(library="ak")
el_trk  = tree["ElectronsAuxDyn.trackParticleLinks"].array(library="ak")
el_ptvarcone= tree["ElectronsAuxDyn.ptvarcone20"].array(library="ak")
el_topoetacone = tree["ElectronsAuxDyn.topoetcone20"].array(library="ak")
mu_pt   = tree["MuonsAuxDyn.pt"].array(library="ak")
mu_charge = tree["MuonsAuxDyn.charge"].array(library="ak")
mu_eta  = tree["MuonsAuxDyn.eta"].array(library="ak")
mu_trk  = tree["MuonsAuxDyn.inDetTrackParticleLink"].array(library="ak")
mu_ptvarcone= tree["MuonsAuxDyn.ptvarcone20"].array(library="ak")
mu_topoetacone = tree["MuonsAuxDyn.topoetcone20"].array(library="ak")
# 'MuonsAuxDyn.inDetTrackParticleLink', 
# 'MuonsAuxDyn.inDetTrackParticleLink/MuonsAuxDyn.inDetTrackParticleLink.m_persKey', 
# 'MuonsAuxDyn.inDetTrackParticleLink/MuonsAuxDyn.inDetTrackParticleLink.m_persIndex',



# TRACK PARTICLES
trk_d0    = tree["InDetTrackParticlesAuxDyn.d0"].array(library="ak")
trk_z0    = tree["InDetTrackParticlesAuxDyn.z0"].array(library="ak")
trk_phi   = tree["InDetTrackParticlesAuxDyn.phi"].array(library="ak")

Ltrk_d0    = tree["InDetLargeD0TrackParticlesAuxDyn.d0"].array(library="ak")
Ltrk_z0    = tree["InDetLargeD0TrackParticlesAuxDyn.z0"].array(library="ak")
Ltrk_phi   = tree["InDetLargeD0TrackParticlesAuxDyn.phi"].array(library="ak")
Ltrk_pt    = tree['InDetLargeD0TrackParticlesAuxDyn.pt_wrtSV'].array(library="ak")

# -----------------------------------------------------------
# Global Selection Criteria for SIG
# -----------------------------------------------------------
lep_pt_min  = 10.0  # GeV
lep_eta_max = 2.5
jet_pt_min  = 10.0   # GeV

# isolation cuts
lep_rel_ptvarcone20_max  = 0.15
lep_rel_topoetcone20_max = 0.15

# -----------------------------
# Histograms (truth vs reco)
# -----------------------------
L_min = 0.01
L_max = 110.0 #5.0  #50.0
nbins_log = 100
edges = np.logspace(np.log10(L_min), np.log10(L_max), nbins_log + 1)
edges2 = np.logspace(np.log10(1e-6), np.log10(9), nbins_log + 1)     #9
edgespT = np.logspace(np.log10(lep_pt_min), np.log10(5000), nbins_log + 1)    #5000 and 2500

h_Lxyz_truth = ROOT.TH1F(
    "h_Lxyz_truth",
    "Displacement: Reconstructed Lepton-tracks, Truth-BSM; L_{xyz} / mm; Normalized Events",
    nbins_log, edges
)

h_Lxyz_reco = ROOT.TH1F(
    "h_Lxyz_reco",
    "L_{xyz} Truth vs Reco; L_{xyz} [mm]; Normalized Events",
    nbins_log, edges
)
h_Lxyz_reco_e = ROOT.TH1F(
    "h_Lxyz_reco_e",
    "L_{xyz} Truth vs Reco; L_{xyz} [mm]; Normalized Events",
    nbins_log, edges
)
h_Lxyz_reco_m = ROOT.TH1F(
    "h_Lxyz_reco_m",
    "L_{xyz} Truth vs Reco; L_{xyz} [mm]; Normalized Events",
    nbins_log, edges
)

h_Lxyz_reco_L = ROOT.TH1F(
    "h_Lxyz_reco_L",
    "L_{xyz} Truth vs Reco; L_{xyz} [mm]; Normalized Events",
    nbins_log, edges
)

h_d0_l1 = ROOT.TH1F(
    "h_d0_l1",
    "Leading and sub-leading leptons;|d_{0}(l_{1})| [mm];Events",
    nbins_log, edges2
)

h_d0_l2 = ROOT.TH1F(
    "h_d0_l2",
    "Leading and sub-leading leptons;|d_{0}(l_{2})| [mm];Events",
    nbins_log, edges2
)

h2_pt_d0_l1 = ROOT.TH2F(
    "h2_pt_d0_l1",
    "#it{mN=70GeV, mWR=5TeV};p_{T} / GeV;d_{0} / mm",
    nbins_log, edgespT,
    nbins_log, edges2
)

g_l1 = ROOT.TGraph()
g_l1.SetName("g_l1_d0_vs_z0")
g_l1.SetTitle("#it{m_{N}=50GeV, m_{W_{R}}=10TeV};d_{0} / mm;z_{0} / mm")

g_l2 = ROOT.TGraph()
g_l2.SetName("g_l2_d0_vs_z0")
g_l2.SetTitle("#it{mN=50GeV, mWR=10TeV};d_{0} / mm;z_{0} / mm")

# Styling (initial)
# h_Lxyz_truth.SetLineColor(ROOT.kRed + 1)
# h_Lxyz_truth.SetLineWidth(3)
# h_Lxyz_reco.SetLineColor(ROOT.kBlue + 1)
# h_Lxyz_reco.SetLineWidth(3)
# h_Lxyz_reco_L.SetLineColor(ROOT.kGreen + 2)
# h_Lxyz_reco_L.SetLineWidth(3)


h_Lxyz_reco_SS = ROOT.TH1F(
    "h_Lxyz_reco_SS",
    "Same-sign reco",
    nbins_log,
    edges
)
h2_pt_d0_l1_SS = ROOT.TH2F(
    "h2_pt_d0_l1_SS",
    "same-sign;p_{T};d_{0}",
    nbins_log,
    edgespT,
    nbins_log,
    edges2
)
g_l1_SS = ROOT.TGraph()
g_l2_SS = ROOT.TGraph()
g_l1_SS.SetName("g_l1_SS")
g_l2_SS.SetName("g_l2_SS")
# -----------------------------
# Event loop
# -----------------------------
print("starting event loop...")
for evt in range(len(pdgId)):
    # --- TRUTH HNL DV ---
    mask = np.abs(pdgId[evt]) == 9900012  # HNL
    for p_idx in np.where(mask)[0]:
        d_idx   = decay_index[evt][p_idx]
        p_idx_v = prod_index[evt][p_idx]

        # convert to plain Python ints (may be uint32)
        d_idx_int   = int(d_idx)
        p_idxv_int  = int(p_idx_v)

        # skip invalid (null) links, e.g. -1 -> 4294967295
        n_vtx = len(vx_all[evt])
        if d_idx_int < 0 or d_idx_int >= n_vtx:
            continue
        if p_idxv_int < 0 or p_idxv_int >= n_vtx:
            continue

        xd = vx_all[evt][d_idx_int]
        yd = vy_all[evt][d_idx_int]
        zd = vz_all[evt][d_idx_int]

        xp = vx_all[evt][p_idxv_int]
        yp = vy_all[evt][p_idxv_int]
        zp = vz_all[evt][p_idxv_int]

        dx = xd - xp
        dy = yd - yp
        dz = zd - zp

        Lxyz = np.sqrt(dx*dx + dy*dy + dz*dz)
        if Lxyz > 1e-10:
            h_Lxyz_truth.Fill(float(Lxyz))

    #------------------------------------------------------------
    # RECONSTRUCTED LEPTON TRACK ORIGINS (combined e + mu)
    # ------------------------------------------------------------
    leptons_std = []
    leptons_LD0 = []

    # ELECTRONS (standard tracks)
    for i in range(len(el_pt[evt])):
        link = el_trk[evt][i]
        if len(link) == 0:
            continue
        trk_idx = link[0]["m_persIndex"]
        trk_idx_int = int(trk_idx)
        if trk_idx_int < 0 or trk_idx_int >= len(trk_d0[evt]):
            continue
        d0 = trk_d0[evt][trk_idx_int]
        z0 = trk_z0[evt][trk_idx_int]
        phi = trk_phi[evt][trk_idx_int]

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

        x = -d0*np.sin(phi)
        y = d0*np.cos(phi)
        leptons_std.append((el_pt[evt][i]*0.001, x, y, z0, "e", d0, int(el_charge[evt][i])))

    # MUONS (Different Branch-link Structure)
    # print("DEBUG)")
    # print("mu_trk type:", type(mu_trk[evt]))
    # print(mu_trk[evt])
    # print(mu_trk.fields)
    mu_evt = mu_trk[evt]
    mu_idx = mu_evt["MuonsAuxDyn.inDetTrackParticleLink.m_persIndex"]
    for i in range(len(mu_idx)):
        trk_idx_int = int(mu_idx[i])
        if trk_idx_int < 0 or trk_idx_int >= len(trk_d0[evt]):
            continue
        d0 = trk_d0[evt][trk_idx_int]
        z0 = trk_z0[evt][trk_idx_int]
        phi = trk_phi[evt][trk_idx_int]

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

        x = -d0*np.sin(phi)
        y = d0*np.cos(phi)
        leptons_std.append((mu_pt[evt][i]*0.001, x, y, z0, "m", d0, int(mu_charge[evt][i])))

    # LD0 Tracks
    # for i in range(len(Ltrk_pt[evt])):
    #     d0 = Ltrk_d0[evt][i]
    #     z0 = Ltrk_z0[evt][i]
    #     phi = Ltrk_phi[evt][i]
    #     x = -d0*np.sin(phi)
    #     y = d0*np.cos(phi)
    #     leptons_LD0.append((Ltrk_pt[evt][i]*0.001, x, y, z0, "LD0"))

        # Ld0 = Ltrk_d0[evt][trk_idx_int]
        # Lz0 = Ltrk_z0[evt][trk_idx_int]
        # Lphi = Ltrk_phi[evt][trk_idx_int]
        # Lx = -Ld0*np.sin(Lphi)
        # Ly = Ld0*np.cos(Lphi)
        # leptons_LD0.append((Ltrk_pt[evt][i], Lx, Ly, Lz0, "e")) 


    # Require 2 Leptons per event
    if len(leptons_std) >= 2:
        leptons_sorted = sorted(leptons_std, key=lambda x: abs(x[0]), reverse=True)
        l1 = leptons_sorted[0]
        l2 = leptons_sorted[1]

        q1 = l1[6]
        q2 = l2[6]
        isSS = (q1*q2 > 0)

        dx = l1[1] - l2[1]
        dy = l1[2] - l2[2]
        dz = l1[3] - l2[3]
        dL = np.sqrt(dx*dx + dy*dy + dz*dz)
        h_Lxyz_reco.Fill(float(dL))
        if isSS:
            h_Lxyz_reco_SS.Fill(float(dL))
        # unpack d0 and z0
        d0_1 = abs(float(l1[5]))
        z0_1 = abs(float(l1[3]))
        d0_2 = abs(float(l2[5]))
        z0_2 = abs(float(l2[3]))

        # d0_1 = float(l1[5])
        # z0_1 = float(l1[3])
        # d0_2 = float(l2[5])
        # z0_2 = float(l2[3])

        h_d0_l1.Fill(abs(d0_1))
        h_d0_l2.Fill(abs(d0_2))

        n1 = g_l1.GetN()
        g_l1.SetPoint(n1, d0_1, z0_1)
        n2 = g_l2.GetN()
        g_l2.SetPoint(n2, d0_2, z0_2)
        if isSS:
            nSS1 = g_l1_SS.GetN()
            g_l1_SS.SetPoint(nSS1, d0_1, z0_1)
            nSS2 = g_l2_SS.GetN()
            g_l2_SS.SetPoint(nSS2, d0_2, z0_2)
        # if (l1[4]=="e" and l2[4]=="e"):
        #     h_Lxyz_reco_e.Fill(float(dL))
        # elif (l1[4]=="m" and l2[4]=="m"):
        #     h_Lxyz_reco_m.Fill(float(dL))
        pt1 = float(l1[0])  
        pt2 = float(l2[0])
        h2_pt_d0_l1.Fill(pt1, d0_1)
        h2_pt_d0_l1.Fill(pt2, d0_2)
        if isSS:
            h2_pt_d0_l1_SS.Fill(pt1,d0_1)
            h2_pt_d0_l1_SS.Fill(pt2,d0_2)
    # if len(leptons_LD0) >= 2:
    #     leptons_LD0_sorted = sorted(leptons_LD0, key=lambda x: x[0], reverse=True)
    #     l1 = leptons_LD0_sorted[0]
    #     l2 = leptons_LD0_sorted[1]
    #     dx = l1[1] - l2[1]
    #     dy = l1[2] - l2[2]
    #     dz = l1[3] - l2[3]
    #     dL = np.sqrt(dx*dx + dy*dy + dz*dz)
    #     h_Lxyz_reco_L.Fill(float(dL))

# -----------------------------
# Normalize histograms
# -----------------------------
if h_Lxyz_truth.Integral() > 0:
    h_Lxyz_truth.Scale(1.0 / h_Lxyz_truth.Integral())

if h_Lxyz_reco.Integral() > 0:
    h_Lxyz_reco.Scale(1.0 / h_Lxyz_reco.Integral())

# if h_Lxyz_reco_e.Integral() > 0:
#     h_Lxyz_reco_e.Scale(1.0 / h_Lxyz_reco_e.Integral())

# if h_Lxyz_reco_m.Integral() > 0:
#     h_Lxyz_reco_m.Scale(1.0 / h_Lxyz_reco_m.Integral())

# if h_Lxyz_reco_L.Integral() > 0:
#     h_Lxyz_reco_L.Scale(1.0 / h_Lxyz_reco_L.Integral())

if h_d0_l1.Integral() > 0:
    h_d0_l1.Scale(1.0 / h_d0_l1.Integral())

if h_d0_l2.Integral() > 0:
    h_d0_l2.Scale(1.0 / h_d0_l2.Integral())

if h_Lxyz_reco_SS.Integral() > 0:
    h_Lxyz_reco_SS.Scale(
        1.0 / h_Lxyz_reco_SS.Integral())
# -----------------------------
# Exponential model: P(L; mu) = 1/mu * exp(-L/mu)
# -----------------------------
L_mean = 9.32  # 0.135662 (70G, 5T)  9.32 (50G, 10T)
Ntot = 1

hExp = h_Lxyz_reco.Clone("h_disp_exponential")
hExp.Reset("ICES")
hExp.SetTitle("")

nbins = hExp.GetNbinsX()
for ibin in range(1, nbins + 1):
    x = hExp.GetBinCenter(ibin)   # in mm
    width = hExp.GetBinWidth(ibin)
    fL = (1.0/L_mean)*np.exp(-x/L_mean)
    val = Ntot* fL*width
    hExp.SetBinContent(ibin, val)

# -----------------------------
# Styling (match C++ macro)
# -----------------------------
# Style reco
h_Lxyz_reco.SetLineColor(ROOT.kRed + 1)
h_Lxyz_reco.SetLineWidth(2)
h_Lxyz_reco.SetFillStyle(0)
h_Lxyz_reco.SetMarkerStyle(20)
h_Lxyz_reco.SetMarkerColor(ROOT.kRed + 1)
h_Lxyz_reco.SetMarkerSize(0.8)

h_Lxyz_reco_SS.SetLineColor(ROOT.kRed + 1)
h_Lxyz_reco_SS.SetLineWidth(2)
h_Lxyz_reco_SS.SetFillColor(ROOT.kRed + 1)
h_Lxyz_reco_SS.SetFillStyle(3354); 

# h_Lxyz_reco_L.SetLineColor(ROOT.kGreen + 1)
# h_Lxyz_reco_L.SetLineWidth(2)
# h_Lxyz_reco_L.SetFillStyle(0)
# h_Lxyz_reco_L.SetMarkerStyle(22)
# h_Lxyz_reco_L.SetMarkerColor(ROOT.kGreen + 1)
# h_Lxyz_reco_L.SetMarkerSize(0.8)

# h_Lxyz_reco_e.SetLineColor(ROOT.kGreen + 2)
# h_Lxyz_reco_e.SetLineWidth(2)
# h_Lxyz_reco_m.SetLineColor(ROOT.kYellow + 2)
# h_Lxyz_reco_m.SetLineWidth(2)

# Style truth
h_Lxyz_truth.SetLineColor(ROOT.kBlue + 1)
h_Lxyz_truth.SetLineWidth(2)
h_Lxyz_truth.SetFillStyle(0)
h_Lxyz_truth.SetMarkerStyle(24)
h_Lxyz_truth.SetMarkerColor(ROOT.kBlue + 1)
h_Lxyz_truth.SetMarkerSize(0.8)

# Style exponential model
hExp.SetLineColor(ROOT.kBlack)
hExp.SetLineWidth(2)
hExp.SetFillStyle(0)
hExp.SetMarkerStyle(1)
hExp.SetMarkerSize(0.0)

# -----------------------------
# Plot (reco vs truth vs exponential model)
# -----------------------------
cL = ROOT.TCanvas("c_disp_reco_vs_truth","Debug: reco vs truth displacement (signal)",1000, 600)
cL.cd()
ROOT.gPad.SetLogy(False)
ROOT.gPad.SetLogx(True)
ROOT.gStyle.SetOptStat(0)
maxData = max(h_Lxyz_reco.GetMaximum(), h_Lxyz_truth.GetMaximum())
# maxData = max(maxData, h_Lxyz_reco_L.GetMaximum())
maxExp = hExp.GetMaximum()
maxY = max(maxData, maxExp)

h_Lxyz_truth.SetMinimum(0.0)
h_Lxyz_truth.SetMaximum(1.1 * maxY)
h_Lxyz_truth.SetTitle("Displacement: Reconstructed Lepton-tracks, Truth-BSM;L_{xyz} / mm; event-rate (normalized)")
h_Lxyz_truth.GetXaxis().SetTitleSize(0.045)
h_Lxyz_truth.GetXaxis().SetLabelSize(0.04)
h_Lxyz_truth.GetYaxis().SetTitleSize(0.045)
h_Lxyz_truth.GetYaxis().SetLabelSize(0.04)
h_Lxyz_truth.GetYaxis().SetTitleOffset(1.15)   
h_Lxyz_truth.GetYaxis().SetTitleFont(42)
h_Lxyz_truth.GetXaxis().SetTitleFont(42)
h_Lxyz_truth.Draw("HIST E")
h_Lxyz_reco.Draw("HIST E SAME")
h_Lxyz_reco_L.Draw("HIST E SAME")
h_Lxyz_reco_SS.Draw("HIST SAME"); 
# h_Lxyz_reco_e.Draw("HIST SAME")
# h_Lxyz_reco_m.Draw("HIST SAME")
hExp.Draw("HIST SAME")

line_mean = ROOT.TLine(L_mean, 0, L_mean, maxY)
line_mean.SetLineColor(ROOT.kBlack)
line_mean.SetLineStyle(2)   
line_mean.SetLineWidth(2)
line_mean.Draw("SAME")

# legL = ROOT.TLegend(0.52, 0.62, 0.95, 0.88)
legL = ROOT.TLegend(0.12, 0.68, 0.55, 0.88)
# legL = ROOT.TLegend(0.32, 0.62, 0.55, 0.88)
legL.SetBorderSize(0)
legL.SetFillStyle(0)
legL.SetTextSize(0.04)
legL.AddEntry(h_Lxyz_reco, "L_{reco}^{dz#phi} between (#it{l_{1}}, #it{l_{2}}) tracks", "pe")
# legL.AddEntry(h_Lxyz_reco_L, "L_{reco}^{dz#phi} between LongD0-tracks", "pe")
legL.AddEntry(h_Lxyz_truth, "L_{truth}^{xyz} BSMtruth: prodV-decV", "pe")
legL.AddEntry(hExp, "P(L; #it{L_{0}}=9.32 mm)", "l")
legL.AddEntry(h_Lxyz_reco_SS,"(same-sign lep)","f")
# legL.AddEntry(hExp, "P(L; #it{L_{0}})=#frac{1}{#it{L_{0}}}exp(-L/#it{L_{0}})", "l")
# legL.AddEntry(line_mean, "#it{L_{0}} = 0.135662 mm", "l")  #9.32  #0.135662
# legL.AddEntry(h_Lxyz_reco_e,  "electron tracks", "l")
# legL.AddEntry(h_Lxyz_reco_m,  "muon tracks", "l")
legL.Draw()

ROOT.gPad.RedrawAxis()
pt = ROOT.gPad.GetPrimitive("title")
if pt:
    pt.SetTextSize(0.045)
    pt.SetTextFont(32)
cL.Modified()
cL.Update()
out_name = os.path.join(output_dir, "LLP_RECOvsBSMTRUTH_with_exp_model.png")
cL.SaveAs(out_name)


c_scatter = ROOT.TCanvas("c_d0_z0_l1_l2", "d0 vs z0 for leading and sub-leading leptons", 1000, 600)
c_scatter.cd()
ROOT.gPad.SetLogx(True)
ROOT.gPad.SetLogy(True)
# ROOT.gPad.SetGrid(0)

# style leading lepton: red
g_l1.SetMarkerStyle(20)
g_l1.SetMarkerSize(0.7)
g_l1.SetMarkerColor(ROOT.kRed + 1)

# style sub-leading: blue
g_l2.SetMarkerStyle(24)
g_l2.SetMarkerSize(0.7)
g_l2.SetMarkerColor(ROOT.kBlue + 1)

g_l1.Draw("AP")
frame = g_l1.GetHistogram()
g_l2.Draw("P SAME")

g_l1.GetXaxis().SetTitle("d_{0} / mm")
g_l1.GetYaxis().SetTitle("z_{0} / mm")
g_l1.GetXaxis().SetTitleSize(0.045)
g_l1.GetYaxis().SetTitleSize(0.045)

# legend
leg_sc = ROOT.TLegend(0.15, 0.15, 0.45, 0.45)
leg_sc.SetBorderSize(0)
leg_sc.SetFillStyle(0)
leg_sc.SetTextSize(0.04)
leg_sc.AddEntry(g_l1, "leading pT: #it{l_{1}}", "p")
leg_sc.AddEntry(g_l2, "sub-leading pT: #it{l_{2}}", "p")
leg_sc.Draw()

out_scatter = os.path.join(output_dir, "scatter_d0_vs_z0_l1_l2.png")
c_scatter.SaveAs(out_scatter)

c_d0 = ROOT.TCanvas("c_d0_l1_l2","|d0| for leading and sub-leading leptons", 1000, 600)
c_d0.cd()
ROOT.gPad.SetLogx(True)  
ROOT.gPad.SetLogy(False)  
ROOT.gStyle.SetOptStat(0)

# style
h_d0_l1.SetLineColor(ROOT.kRed + 1)
h_d0_l1.SetLineWidth(2)
h_d0_l1.SetMarkerStyle(20)
h_d0_l1.SetMarkerColor(ROOT.kRed + 1)
h_d0_l1.SetMarkerSize(0.8)

h_d0_l2.SetLineColor(ROOT.kBlue + 1)
h_d0_l2.SetLineWidth(2)
h_d0_l2.SetMarkerStyle(24)
h_d0_l2.SetMarkerColor(ROOT.kBlue + 1)
h_d0_l2.SetMarkerSize(0.8)

# axis labels
h_d0_l1.GetXaxis().SetTitle("|d_{0}| / mm")
h_d0_l1.GetYaxis().SetTitle("event rate (normalized)")
h_d0_l1.GetXaxis().SetTitleSize(0.045)
h_d0_l1.GetXaxis().SetLabelSize(0.04)
h_d0_l1.GetYaxis().SetTitleSize(0.045)
h_d0_l1.GetYaxis().SetLabelSize(0.04)

# set y-range
max_d0 = max(h_d0_l1.GetMaximum(), h_d0_l2.GetMaximum())
h_d0_l1.SetMinimum(0.0)
h_d0_l1.SetMaximum(1.2 * max_d0)

# draw
h_d0_l1.Draw("HIST E")
h_d0_l2.Draw("HIST E SAME")

leg_d0 = ROOT.TLegend(0.55, 0.70, 0.90, 0.88)
leg_d0.SetBorderSize(0)
leg_d0.SetFillStyle(0)
leg_d0.SetTextSize(0.04)
leg_d0.AddEntry(h_d0_l1, "leading pT: #it{l_{1}}", "l")
leg_d0.AddEntry(h_d0_l2, "sub-leading pT: #it{l_{2}}", "l")
leg_d0.Draw()

out_d0 = os.path.join(output_dir, "d0_l1_l2_logbins.png")
c_d0.SaveAs(out_d0)

# -----------------------------
# 2D plots: pT vs |d0|
# -----------------------------
c_pt_d0_l1 = ROOT.TCanvas("c_pt_d0_l1", "pT vs |d0| for leading lepton", 600, 600)
c_pt_d0_l1.cd()
ROOT.gPad.SetLogx(True)   
ROOT.gPad.SetLogy(True)   
ROOT.gPad.SetLogz(False)  
ROOT.gStyle.SetOptStat(0)

h2_pt_d0_l1.GetXaxis().SetTitle("#it{p_{T}} / GeV")
h2_pt_d0_l1.GetYaxis().SetTitle("|#it{d_{0}}| / mm")
h2_pt_d0_l1.GetXaxis().SetTitleSize(0.04)
h2_pt_d0_l1.GetXaxis().SetLabelSize(0.03)
h2_pt_d0_l1.GetYaxis().SetTitleSize(0.04)
h2_pt_d0_l1.GetYaxis().SetLabelSize(0.03)
h2_pt_d0_l1.GetYaxis().SetTitleOffset(0.8) 
ROOT.gStyle.SetTitleFontSize(0.05) 
h2_pt_d0_l1.Draw("COL")
out_pt_d0_l1 = os.path.join(output_dir, "pt_vs_d0.png")
c_pt_d0_l1.SaveAs(out_pt_d0_l1)


c_pt_d0_l1_SS = ROOT.TCanvas("c_pt_d0_l1_SS", "pT vs |d0| for leading lepton", 600, 600)
c_pt_d0_l1_SS.cd()
ROOT.gPad.SetLogx(True)   
ROOT.gPad.SetLogy(True)   
ROOT.gPad.SetLogz(False)  
ROOT.gStyle.SetOptStat(0)

h2_pt_d0_l1_SS.GetXaxis().SetTitle("#it{p_{T}} / GeV")
h2_pt_d0_l1_SS.GetYaxis().SetTitle("|#it{d_{0}}| / mm")
h2_pt_d0_l1_SS.GetXaxis().SetTitleSize(0.04)
h2_pt_d0_l1_SS.GetXaxis().SetLabelSize(0.03)
h2_pt_d0_l1_SS.GetYaxis().SetTitleSize(0.04)
h2_pt_d0_l1_SS.GetYaxis().SetLabelSize(0.03)
h2_pt_d0_l1_SS.GetYaxis().SetTitleOffset(0.8) 
ROOT.gStyle.SetTitleFontSize(0.05) 
h2_pt_d0_l1_SS.Draw("COL")
out_pt_d0_l1_SS = os.path.join(output_dir, "pt_vs_d0_SS.png")
c_pt_d0_l1_SS.SaveAs(out_pt_d0_l1_SS)

print(f"Output successfully written to '{out_name}'")
